#include "bcifit.hpp"
#include "eleclj.hpp"
#include "dmat.hpp"
#include "linalg.hpp"
#include "fns.h"
#include "cmplx.hpp"

MolSpec::MolSpec() : number_of_symmetry_copies(1)
{
  dipole.x = UNDEF_VAL;
  quadrupole.XX = UNDEF_VAL;
  octopole.XXX = UNDEF_VAL;
  hexadecapole.XXXX = UNDEF_VAL;
}

BCIFit::BCIFit() : 
  verbose(1),
  kspace_cutoff(5),
  elec_1_4_scale_factor(1.0),
  atomic_charge_weight(1.0),
  dipole_moment_weight(1.0),
  quadrupole_moment_weight(1.0),
  octopole_moment_weight(1.0),
  hexadecapole_moment_weight(1.0),
  structure_factor_weight(1.0),
  svd_tolerance(1e-3),
  use_traceless_multipoles(true),
  tolerance(1e-5),
  bci_minimum(0.0),
  max_iterations(10000)
{ }

static double help(const gsl_vector *x, void *p)
{
  return ((BCIFit *) p)->help(x);
}

static double distance(const Tensor &a, const Tensor &b)
{
  const Tensor tmp = a - b;
  return tmp.col(0).sq() + tmp.col(1).sq() + tmp.col(2).sq();
}

double BCIFit::help(const gsl_vector *x)
{
  int i;
  HTabIterator<BondChargeIncrementParam> b;
  for (i = 0, b.init(bci_to_fit); b.ok(); b.next(), i++)
    if ((bcitab[b.key()].k = b.val().k = gsl_vector_get(x,i)) < bci_minimum)
      return 1e20;
  LstItrMod<MolSpec> ms;
  double cost = 0;
  for (ms.init(molecule_specification), i = 0; ms.ok(); ms.next(), i++) {
    ms().elec.bond_charge_increment = bcitab;
    ms().elec.types_changed();
    cost += ::distance(ms().elec.polarizability(),ms().polarizability);
  }
  return cost;
}

void BCIFit::fit_bci()
{
  Out() << "Starting BCIFit::fit_bci...\n";
  HTabIterator<BondChargeIncrementParam> b;
  const int nparam = bci_to_fit.size();
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *fmin = gsl_multimin_fminimizer_alloc(T, nparam);
  gsl_multimin_function f;
  gsl_vector *x = gsl_vector_alloc(nparam);
  gsl_vector *step = gsl_vector_alloc(nparam);
  f.n = nparam;
  f.f = &::help;
  f.params = this;
  int i;
  bcitab = bci_to_fix;
  for (i = 0, b.init(bci_to_fit); b.ok(); b.next(), i++) {
    gsl_vector_set(x, i, bcitab[b.key()].k = b.val().k);
    gsl_vector_set(step, i, b.val().k*0.1);
  }
  LstItrMod<MolSpec> ms;
  for (ms.init(molecule_specification), i = 0; ms.ok(); ms.next(), i++) {
    ms().elec.verbose = verbose;
    ms().elec.kspace_cutoff = kspace_cutoff;
    ms().elec.msite = msite;
    ms().elec.sp3 = sp3;
    ms().elec.sp2 = sp2;
    ms().elec.elec_1_4_scale_factor = elec_1_4_scale_factor;
    ms().elec.fixed_charge = fixed_charge;
    ms().elec.screening_radius = screening_radius;
    ms().elec.bond_charge_increment = bcitab;
    ms().elec.one_three_interaction = one_three_interaction;
    ms().elec.init(&ms().msys);
    ms().elec.coordinates_changed();
  }
  gsl_multimin_fminimizer_set(fmin, &f, x, step);
  for (i = 1; i <= max_iterations; i++)
    if (gsl_multimin_fminimizer_iterate(fmin) ||
	gsl_multimin_test_size(gsl_multimin_fminimizer_size(fmin), tolerance) ==
	GSL_SUCCESS)
      break;
  for (i = 0, b.init(bci_to_fit); b.ok(); b.next(), i++)
    b.val().k = gsl_vector_get(fmin->x,i);
  gsl_multimin_fminimizer_free(fmin);
  gsl_vector_free(x);
  gsl_vector_free(step);
  if (strlen(write_parameters_to) > 0) {
    ostream *fp = FileStream(write_parameters_to);
    *fp << bcitab;
    delete fp;
  }
  Out() << "BCI parameters: " << bcitab << "\n\n" << flush;
  double cost = 0;
  for (ms.init(molecule_specification), i = 0; ms.ok(); ms.next(), i++) {
    ms().elec.bond_charge_increment = bcitab;
    ms().elec.types_changed();
    const Tensor a = ms().elec.polarizability();
    Out() << "Target polarizability (A^3): " << ms().polarizability << "\n"
	  << "Model polarizability (A^3): " << a << "\n";
    cost += ::distance(ms().polarizability,a);
  }
  Out() << "Cost function: " << cost << "\n";
}

void BCIFit::fit()
{
  Out() << "Starting BCIFit::fit...\n"
	<< "Using traceless multipoles: " << use_traceless_multipoles << "\n"
	<< flush;
  int i, ntot = 0;
  const Cartesian zero(0,0,0);
  LstItrMod<MolSpec> ms;
  bcitab = bci_to_fix;
  HTabIterator<BondChargeIncrementParam> b;
  for (b.init(bci_to_fit); b.ok(); b.next()) {
    b.val().q0 = 0.0;
    bcitab[b.key()] = b.val();
  }
  for (i = 0, ms.init(molecule_specification); ms.ok(); i++, ms.next()) {
    const MSys &msys = ms().msys;
    const int nat = msys.atom.size();
    if (ms().charges.size() != nat)
      die("BCIFit::fit: %d atoms, but %d charges, for system %d",
	  nat, ms().charges.size(), i+1);
    ntot += nat;
    if (msys.boundary_conditions->type == non_periodic) {
      if (dipole_moment_weight > 0)
	ntot += 3;
      if (quadrupole_moment_weight > 0)
	ntot += 6;
      if (octopole_moment_weight > 0)
	ntot += 10;
      if (hexadecapole_moment_weight > 0)
	ntot += 15;
    } else {
      ntot += 2*ms().electronic_structure_factors.size();
    }
  }
  const int nbci = bci_to_fit.size();
  DMat A(ntot,nbci,0.0);
  DVec c(ntot,0.0), q(nbci,0.0), target(ntot,0.0);
  int itot = 0;
  for (ms.init(molecule_specification), i = 0; ms.ok(); ms.next(), i++) {
    MSys &msys = ms().msys;
    ms().elec.verbose = verbose;
    ms().elec.kspace_cutoff = kspace_cutoff;
    ms().elec.msite = msite;
    ms().elec.sp3 = sp3;
    ms().elec.sp2 = sp2;
    ms().elec.elec_1_4_scale_factor = elec_1_4_scale_factor;
    ms().elec.fixed_charge = fixed_charge;
    ms().elec.screening_radius = screening_radius;
    ms().elec.bond_charge_increment = bcitab;
    ms().elec.one_three_interaction = one_three_interaction;
    ms().elec.init(&msys);
    ms().elec.coordinates_changed();
    ms().elec.solve_for_bci(zero);
    const int nsymm = ms().number_of_symmetry_copies;
    const double sum_fixed = ms().elec.net_charge();
    double sum_target = 0;
    const int nat = msys.atom.size();
    insist(nat % nsymm == 0);
    for (int j = 0; j < nat; j++, itot++) {
      const double q0 = e_to_charge_unit(ms().charges[j]);
      c[itot] = atomic_charge_weight * (q0 - ms().elec.charge_on_atom(j))/nsymm;
      sum_target += q0;
    }
    if (fabs(sum_fixed-sum_target) > 1e-8) {
      ms().elec.write();
      die("BCIFit::fit: target charges add to %f, but fixed charges add to %f for system %d",
	  charge_unit_to_e(sum_target), 
	  charge_unit_to_e(sum_fixed), i+1);
    }
    if (msys.boundary_conditions->type == non_periodic) {
      if (nsymm != 1)
        die("BCIFit::fit: number_of_symmetry_copies must be 1 for non-periodic boundary conditions");
#define TMP(alpha) c[itot] = w * M.alpha; target[itot] = w * M0.alpha; itot++
      if (dipole_moment_weight > 0) {
	if (are_approximately_equal(ms().dipole.x, UNDEF_VAL))
	  die("BCIFit::fit: must specify dipole for system %d", i+1);
	const double w = nat * dipole_moment_weight;
	const Cartesian M = ms().elec.dipole_moment();
	const Cartesian M0 = ms().dipole.map(debye_to_dipole_unit);
	TMP(x); TMP(y); TMP(z);
      }
      if (quadrupole_moment_weight > 0) {
	if (are_approximately_equal(ms().quadrupole.XX, UNDEF_VAL))
	  die("BCIFit::fit: must specify quadrupole for system %d", i+1);
	const double w = nat * quadrupole_moment_weight;
	Quadrupole M = ms().elec.quadrupole_moment();
	Quadrupole M0 = ms().quadrupole.map(debye_to_dipole_unit);
	if (use_traceless_multipoles) {
	  M.make_traceless();
	  M0.make_traceless();
	}
	TMP(XX); TMP(XY); TMP(XZ); 
	TMP(YY); TMP(YZ); TMP(ZZ);
      }
      if (octopole_moment_weight > 0) {
	if (are_approximately_equal(ms().octopole.XXX, UNDEF_VAL))
	  die("BCIFit::fit: must specify quadrupole for system %d", i+1);
	const double w = nat * octopole_moment_weight;
	Octopole M = ms().elec.octopole_moment();
	Octopole M0 = ms().octopole.map(debye_to_dipole_unit);
	if (use_traceless_multipoles) {
	  M.make_traceless();
	  M0.make_traceless();
	}
	TMP(XXX); TMP(YYY); TMP(ZZZ);
	TMP(XYY); TMP(XXY); TMP(XXZ);
	TMP(XZZ); TMP(YZZ); TMP(YYZ);
	TMP(XYZ);
      }
      if (hexadecapole_moment_weight > 0) {
	if (are_approximately_equal(ms().hexadecapole.XXXX, UNDEF_VAL))
	  die("BCIFit::fit: must specify quadrupole for system %d", i+1);
	const double w = nat * hexadecapole_moment_weight;
	const Hexadecapole M = ms().elec.hexadecapole_moment();
	const Hexadecapole M0 = ms().hexadecapole.map(debye_to_dipole_unit);
	TMP(XXXX); TMP(YYYY); TMP(ZZZZ); TMP(XXXY); 
	TMP(XXXZ); TMP(YYYX); TMP(YYYZ); TMP(ZZZX); 
	TMP(ZZZY); TMP(XXYY); TMP(XXZZ); TMP(YYZZ);
	TMP(XXYZ); TMP(YYXZ); TMP(ZZXY);
      }
#undef TMP
    } else { /* periodic boundary conditions */
      const double w = nat * structure_factor_weight/nsymm;
      LstItr<StructureFactor> sk;
      const Tensor h = msys.boundary_conditions->reciprocal_lattice_vectors(); 
      for (sk.init(ms().electronic_structure_factors); sk.ok(); sk.next()) {
	const Cartesian k =  h * sk().k;
	const Complex Sk = sk().Sk.map(e_to_charge_unit) + msys.nuclear_structure_factor(k) - 
	  ms().elec.structure_factor(k);
	c[itot++] = w * Sk.re;
	c[itot++] = w * Sk.im;
      }
    }
  }
  insist(itot == ntot);
  int ib;
  for (ib = 0, b.init(bci_to_fit); b.ok(); ib++, b.next()) {
    itot = 0;
    for (ms.init(molecule_specification); ms.ok(); ms.next()) {
      const MSys &msys = ms().msys;
      ms().elec.bond_charge_increment[b.key()].q0 = 1.0;
      ms().elec.types_changed();
      ms().elec.update_q();
      ms().elec.solve_for_bci(zero);
      const int nat = msys.atom.size();
      const int nsymm = ms().number_of_symmetry_copies;
      for (int j = 0; j < nat; j++, itot++)
	A(itot,ib) = atomic_charge_weight * ms().elec.charge_on_atom(j)/nsymm;
      if (msys.boundary_conditions->type == non_periodic) {
#define TMP(alpha) A(itot,ib) = w*M.alpha - c[itot]; itot++;
	if (dipole_moment_weight > 0) {
	  const double w = msys.atom.size() * dipole_moment_weight;
	  const Cartesian M = ms().elec.dipole_moment();
	  TMP(x); TMP(y); TMP(z);
	}
	if (quadrupole_moment_weight > 0) {
	  const double w = msys.atom.size() * quadrupole_moment_weight;
	  Quadrupole M = ms().elec.quadrupole_moment();
	  if (use_traceless_multipoles)
	    M.make_traceless();
	  TMP(XX); TMP(XY); TMP(XZ); 
	  TMP(YY); TMP(YZ); TMP(ZZ);
	}
	if (octopole_moment_weight > 0) {
	  const double w = msys.atom.size() * octopole_moment_weight;
	  Octopole M = ms().elec.octopole_moment();
	  if (use_traceless_multipoles)
	    M.make_traceless();
	  TMP(XXX); TMP(YYY); TMP(ZZZ);
	  TMP(XYY); TMP(XXY); TMP(XXZ);
	  TMP(XZZ); TMP(YZZ); TMP(YYZ);
	  TMP(XYZ);
	}
	if (hexadecapole_moment_weight > 0) {
	  const double w = msys.atom.size() * hexadecapole_moment_weight;
	  const Hexadecapole M = ms().elec.hexadecapole_moment();
	  TMP(XXXX); TMP(YYYY); TMP(ZZZZ); TMP(XXXY); 
	  TMP(XXXZ); TMP(YYYX); TMP(YYYZ); TMP(ZZZX); 
	  TMP(ZZZY); TMP(XXYY); TMP(XXZZ); TMP(YYZZ);
	  TMP(XXYZ); TMP(YYXZ); TMP(ZZXY);
	}
#undef TMP
      } else { /* periodic boundary conditions */
        const Tensor h = msys.boundary_conditions->reciprocal_lattice_vectors();
	const double w = nat * structure_factor_weight/nsymm;
	LstItr<StructureFactor> sk;
	for (sk.init(ms().electronic_structure_factors); sk.ok(); sk.next()) {
	  const Cartesian k = h * sk().k;
	  const Complex Sk = ms().elec.structure_factor(k);
	  A(itot++,ib) = w * Sk.re;
	  A(itot++,ib) = w * Sk.im;
	}
      }
      ms().elec.bond_charge_increment[b.key()].q0 = 0;
    }
  }
  insist(itot == ntot);
  Out() << "Fitting " << nbci << " parameters to " << ntot << " target values.\n";
  if (nbci > ntot)
    die("BCIFit::fit: number of parameters is greater than or equal to number of target values");
  if (nbci > 0) {
    target -= c;
    if (verbose > 2) {
      Out() << "A: " << A << "\n";
      Out() << "target: " << target << "\n";
    }
    DMat Acopy = A.copy();
    DVec tcopy = target.copy();
    SVDLeastSquares(A,q,target,svd_tolerance,true);
    DVec tmp = Acopy*q;
    tmp -= tcopy;
    Out() << "Residual: " << tmp.sq() << "\n";
  }
  for (ib = 0, b.init(bci_to_fit); b.ok(); ib++, b.next())
    b.val().q0 = bcitab[b.key()].q0 = is_almost_zero(q[ib]) ? 0 : q[ib];
  if (strlen(write_parameters_to) > 0) {
    ostream *fp = FileStream(write_parameters_to);
    *fp << bcitab;
    delete fp;
  }
  Out() << "BCI parameters: " << bcitab << "\n\n";
  double msdtot = 0;
  int nmsdtot = 0;
  for (i = 0, ms.init(molecule_specification); ms.ok(); i++, ms.next()) {
    const MSys &msys = ms().msys;
    ms().elec.fixed_charge = fixed_charge;
    ms().elec.screening_radius = screening_radius;
    ms().elec.bond_charge_increment = bcitab;
    ms().elec.types_changed();
    ms().elec.update_q();
    ms().elec.solve_for_bci(zero);
    if (verbose > 1)
      ms().elec.write();
    Out() << "System " << i+1 << ":\n";
    double msd = 0;
    const int nat = ms().msys.atom.size();
    OutPrintf("%-8s %-8s %12s %12s\n", "symbol", "type", "charge (e)", "target (e)");
    for (int j = 0; j < nat; j++) {
      const double qmodel = charge_unit_to_e(ms().elec.charge_on_atom(j));
      const double qtarget = ms().charges[j];
      OutPrintf("%-8s %-8s %12.8f %12.8f\n", 
		(const char *) ms().msys.atom[j].symbol, 
		(const char *) ms().msys.atom[j].type, 
		qmodel, qtarget);
      msd += sq(qmodel-qtarget);
    }
    Out() << "Net charge: " << ms().elec.net_charge() << " e\n"
	  << "RMSD: " << sqrt(msd/nat) << " e\n\n";
    if (msys.boundary_conditions->type == non_periodic) {
      const Cartesian mu = ms().elec.dipole_moment().map(dipole_unit_to_debye);
      OutPrintf("Model dipole moment (Debye):  %12.8f %12.8f %12.8f\n", mu.x, mu.y, mu.z);
      OutPrintf("Target dipole moment (Debye): %12.8f %12.8f %12.8f\n\n", 
		ms().dipole.x, ms().dipole.y, ms().dipole.z);
      if (quadrupole_moment_weight > 0) {
	Quadrupole M = ms().elec.quadrupole_moment().map(dipole_unit_to_debye);
	Quadrupole M0 = ms().quadrupole;
	if (use_traceless_multipoles) {
	  M.make_traceless();
	  M0.make_traceless();
	}
	Out() << "Model quadrupole moment (Debye A):\n" 
	      << M
	      << "\n"
	      << "Target quadrupole moment (Debye A):\n" 
	      << M0
	      << "\n";
      }
      if (octopole_moment_weight > 0) {
	Octopole M = ms().elec.octopole_moment().map(dipole_unit_to_debye);
	Octopole M0 = ms().octopole;
	if (use_traceless_multipoles) {
	  M.make_traceless();
	  M0.make_traceless();
	}
	Out() << "Model octopole moment (Debye A):\n" 
	      << M
	      << "\n"
	      << "Target octopole moment (Debye A):\n" 
	      << M0
	      << "\n";
      }
      if (hexadecapole_moment_weight > 0) {
	Hexadecapole M = ms().elec.hexadecapole_moment().map(dipole_unit_to_debye);
	Hexadecapole M0 = ms().hexadecapole;
	if (use_traceless_multipoles) {
	  M.make_traceless();
	  M0.make_traceless();
	}
	Out() << "Model hexadecapole moment (Debye A):\n" 
	      << M
	      << "\n"
	      << "Target hexadecapole moment (Debye A):\n" 
	      << M0
	      << "\n";
      }
    } else {
      if (ms().electronic_structure_factors.size() > 0) {
	OutPrintf("Electronic structure factors (atomic units):\n");
	OutPrintf("%12s %12s %12s   %12s %12s   %12s %12s\n",
		  "kx","ky","kz", "Model Re","Im", "Target Re","Im");
	LstItr<StructureFactor> sk;
        const Tensor h = msys.boundary_conditions->reciprocal_lattice_vectors();
	for (sk.init(ms().electronic_structure_factors); sk.ok(); sk.next()) {
	  const Cartesian k = h * sk().k;
	  const Complex selec = (ms().elec.structure_factor(k) -
				 msys.nuclear_structure_factor(k)).map(charge_unit_to_e);
	  OutPrintf("%12.8f %12.8f %12.8f   %12.8f %12.8f   %12.8f %12.8f\n",
		    sk().k.x, sk().k.y, sk().k.z, 
		    selec.re, selec.im, 
		    sk().Sk.re, sk().Sk.im);
	}
      }
    }
    const int nsymm = ms().number_of_symmetry_copies;
    msdtot += msd/nsymm;
    nmsdtot += nat/nsymm;
    Out() << "\n";
  }
  Out() << "Total RMSD: " << sqrt(msdtot/nmsdtot) << " e" << "\n\n" << flush;
}
