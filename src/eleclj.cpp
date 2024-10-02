#include "eleclj.hpp"
#include "msys.hpp"
#include "rspace.h"
#include "kspace.h"
#include "timing.h"
#include "linalg.hpp"
#include "nbrlist.hpp"
#include "pbin.h"
#include "cmplx.hpp"

#define DEFAULT_P3M_GRID_SPACING 0.5

ElecLJ::ElecLJ() : 
  elec_1_4_scale_factor(1), lj_1_4_scale_factor(1),
  rspace_cutoff(0), smoothing_width(1), neighbor_list_width(0), 
  kspace_cutoff(5), ewald_screening(-1), tolerance(1e-10),
  bond_charge_increment_frequency(1500), 
  p3m_grid_spacing(DEFAULT_P3M_GRID_SPACING), p3m_grid_density(-1),
  verbose(1), max_iterations(1000), p3m_alias_sum_limit(2),
  p3m_assignment_order(3), reorder_sites(false),
  geometric_combining_for_sigma(false), include_lj_correction(true),
  include_rspace(true), include_kspace(true), include_self(true), 
  use_p3m(false), p3m_ik_differentiate(false), 
  kT(K_to_kcal_mol(298.15)),
  nsite(0), nbci(0), ngroup(0),
  any_fixed_charge_params(false), any_bci_params(false), 
  any_lj_params(false), any_virtual_sites(false),
  any_site_has_q(false), any_site_has_lj(false), 
  any_bci(false), have_solved_for_bci(false),
  dyn(false), totA(0), totB(0),
  lj_correction(0), elec_energy(0), lj_energy(0), rc(0), nlw(0), eta(0),
  nlist(0), excluded(0), group_excluded(0), one_four(0),
  rspace_chebyshev(0), kspace(0), p3m(0)
{ 
  sigma["Lj"] = epsilon["Lj"] = 1.0;
  sigma["dummy"] = 1.0;
  epsilon["dummy"] = 0.0;
  screening_radius["dummy"] = 2.0;
  fixed_charge["dummy"] = 0.0;
}

ElecLJ::ElecLJ(const ElecLJ &t) : 
  msite(t.msite),
  sp3(t.sp3),
  sp2(t.sp2),
  fixed_charge(t.fixed_charge),
  screening_radius(t.screening_radius),
  bond_charge_increment(t.bond_charge_increment),
  one_three_interaction(t.one_three_interaction),
  sigma(t.sigma),
  epsilon(t.epsilon),
  general_type(t.general_type),
  elec_1_4_scale_factor(t.elec_1_4_scale_factor),
  lj_1_4_scale_factor(t.lj_1_4_scale_factor),
  rspace_cutoff(t.rspace_cutoff),
  smoothing_width(t.smoothing_width),
  neighbor_list_width(t.neighbor_list_width),
  kspace_cutoff(t.kspace_cutoff),
  ewald_screening(t.ewald_screening),
  tolerance(t.tolerance),
  bond_charge_increment_frequency(t.bond_charge_increment_frequency),
  p3m_grid_spacing(t.p3m_grid_spacing),
  p3m_grid_density(t.p3m_grid_density),
  verbose(t.verbose),
  max_iterations(t.max_iterations),
  p3m_alias_sum_limit(t.p3m_alias_sum_limit),
  p3m_assignment_order(t.p3m_assignment_order),
  reorder_sites(t.reorder_sites),
  geometric_combining_for_sigma(t.geometric_combining_for_sigma),
  include_lj_correction(t.include_lj_correction),
  include_rspace(t.include_rspace),
  include_kspace(t.include_kspace),
  include_self(t.include_self),
  use_p3m(t.use_p3m),
  p3m_ik_differentiate(t.p3m_ik_differentiate),
  kT(t.kT),
  nsite(0), ngroup(0),
  any_fixed_charge_params(false), any_bci_params(false), 
  any_lj_params(false), any_virtual_sites(false),
  any_site_has_q(false), any_site_has_lj(false), 
  any_bci(false), have_solved_for_bci(false),
  dyn(false), totA(0), totB(0),
  lj_correction(0), elec_energy(0), lj_energy(0), rc(0), nlw(0), eta(0),
  nlist(0), excluded(0), group_excluded(0), one_four(0),
  rspace_chebyshev(0), kspace(0), p3m(0)
{ }

ElecLJ::~ElecLJ()
{
  cleanup();
}

void ElecLJ::cleanup()
{
  if (excluded)
    PairBinary_delete(excluded);
  if (group_excluded)
    PairBinary_delete(group_excluded);
  if (one_four)
    PairBinary_delete(one_four);
  if (rspace_chebyshev)
    free(rspace_chebyshev);
  if (kspace)
    kspace_delete(kspace);
  if (p3m)
    p3m_delete(p3m);
  if (nlist)
    delete nlist;
  excluded = 0;
  group_excluded = 0;
  one_four = 0;
  rspace_chebyshev = 0;
  kspace = 0;
  p3m = 0;
  nlist = 0;
  any_site_has_q = false;
  any_site_has_lj = false;
  dyn = false;
  any_bci = false;
  ngroup = 0;
  for (int i = 0; i < nsite; i++) {
    delete site[i];
    site[i] = 0;
  }
  nsite = 0;
}

/* Error estimates from kolafa92 */
static double estimated_dimensionless_kspace_error(double kc, double et, double volume)
{
  const double l = cbrt(volume);
  return et*l/M_PI * sqrt(8/kc) * exp(-sq(M_PI*kc/(et*l)));
}

static double estimated_dimensionless_rspace_error(double rc, double et, double volume)
{
  return 2/sqrt(rc/cbrt(volume)) * exp(-sq(rc*et));
}

void ElecLJ::init(MSys *m)
{
  Potential::init(m);
  cleanup();
  any_fixed_charge_params = fixed_charge.size() > 0;
  any_bci_params = bond_charge_increment.size() > 0;
  any_lj_params = sigma.size() > 0;
  any_virtual_sites = msite.size() > 0 || sp3.size() > 0 || sp2.size() > 0;
  if (!(any_fixed_charge_params || any_bci_params || any_lj_params))
    return;
  Out() << "Starting ElecLJ::init...\n" << flush;
  /* Check some input parameters */
  insist(tolerance > 0);
  insist(max_iterations > 0);
  site_on_atom = Vec<Sites>(msys->atom.size());
  /* Initialize all sites and parameters */
  types_changed();
  /* Initialize cutoffs and Ewald */
  const BoundaryConditions &bc = *msys->boundary_conditions;
  if (bc.type == non_periodic) {
    if (rspace_cutoff > 0)
      Out() << "*** ElecLJ::init: cutoff will be ignored with nonperiodic boundary conditions";
  } else {
    /* Periodic boundary conditions */
    rc = rspace_cutoff;
    const double min_diameter = bc.min_diameter();
    if (is_almost_zero(rc))
      rc = 0.45 * min_diameter;
    if (rc > 0.49*min_diameter)
      die("ElecLJ::init: rspace_cutoff (%f) must be <= 0.49 * minimum diameter of box (%f)",
	  rc, min_diameter);
    insist(rc > 0);
    nlw = neighbor_list_width;
    if (is_almost_zero(nlw)) {
      nlw = 0.2*rc;
      if (rc + nlw >= 0.5*min_diameter)
	nlw = 0.499*min_diameter - rc;
    }
    nlist = new NeighborList(ngroup, rc, nlw);
    nlist->set_boundary_conditions(bc);
    if (smoothing_width >= rc)
      die("ElecLJ::init: smoothing_width (%f) must be less than real-space cutoff (%f)",
	  smoothing_width,rc);
    Out() << "Real-space cutoff: " << rc << " A\n";
    if (smoothing_width > 0)
      Out() << "Smoothing width: " << smoothing_width << " A\n";
    else
      Out() << "No smoothing.\n";
    Out() << "Neighbor list width: " << nlw << " A\n";
    eta = ewald_screening;
    if (any_site_has_q) {
      double kerr = 0;
      if (include_kspace && use_p3m) {
	icart_t ngrid;
	const Tensor a = bc.lattice_vectors();
	double spacing = p3m_grid_spacing;
	if (p3m_grid_density > 0) {
	  if (!are_approximately_equal(spacing,DEFAULT_P3M_GRID_SPACING))
	    die("ElecLJ::init: only one of 'p3m_grid_spacing','p3m_grid_density' should be specified");
	  spacing = bc.min_diameter()/cbrt(p3m_grid_density*bc.volume());
	}
	p3m = p3m_new(nsite, (const double *) &a, 
		      spacing, p3m_alias_sum_limit, p3m_ik_differentiate, 
		      p3m_assignment_order, (int *) &ngrid);
	if (eta < 0) {
	  Out() << "Optimizing ewald screening parameter...\n" << flush;
	  double chi_rspace = 0;
	  p3m_optimize(p3m, rc, &eta, &chi_rspace, &kerr);
	  Out() << "done.\n" << flush;
	} else {
	  p3m_init(p3m, eta, &kerr);
	}
	const int ngridtot = ngrid.x*ngrid.y*ngrid.z;
	Out() << "Using P3M with " 
	      << (p3m_ik_differentiate ? "reciprocal space" : "real space")
	      << " differentiation.\n"
	      << "Number of grid points: " 
	      << ngrid.x << " x " << ngrid.y << " x " << ngrid.z << " = " << ngridtot << "\n"
	      << "Grid spacing (A): "
	      << a.col(0).magnitude()/ngrid.x << " " 
	      << a.col(1).magnitude()/ngrid.y << " " 
	      << a.col(2).magnitude()/ngrid.z << "\n"
	      << "Grid density: " << ngridtot/bc.volume() << " A^-3\n"
	      << "Assignment order: " << p3m_assignment_order << "\n"
	      << "Limit of alias sum: " << p3m_alias_sum_limit << "\n"
	      << "Estimated dimensionless P3M RMS force error: " << kerr << "\n"
	      << flush;
      } else { /* ordinary Ewald */
        Out() << "Using ordinary Ewald.\n";
	if (eta < 0)
	  eta = sqrt((M_PI*kspace_cutoff)/(rc*min_diameter));
        if (include_kspace) {
	  kspace = kspace_new(kspace_cutoff, nsite);
          kerr = estimated_dimensionless_kspace_error(kspace_cutoff, eta, bc.volume());
          Out() << "k-space cutoff: " << kspace_cutoff << "\n"
                << "Number of wavevectors: " << kspace_number_of_wavevectors(kspace) << "\n"
	        << "Estimated dimensionless k-space RMS force error: " << kerr << "\n"
	        << flush;
        }
      }
      double rerr = 0;
      if (include_rspace) {
	rerr = estimated_dimensionless_rspace_error(rc, eta, bc.volume());
	Out() << "Estimated dimensionless r-space RMS force error: " << rerr << "\n";
	double maxerr = 0, xerr = 0;
	rspace_chebyshev = RSpaceChebyshev_new(rc, eta, &maxerr, &xerr);
	if (maxerr > 0)
	  Out() << "Maximum absolute error for Chebyshev approximation: " 
		<< maxerr << " at r^2 = " << xerr << "\n"
		<< flush;
      } 
      double sumq2 = 0;
      for (int i = 0; i < nsite; i++)
	sumq2 += sq(site[i]->q0);
      const double err0 = sqrt(sq(rerr)+sq(kerr));
      const double err = nsite > 0 ?
	sumq2/(sqrt(nsite)*pow(bc.volume(),2.0/3.0)) * err0 : 0;
      Out() << "Estimated dimensionless total RMS force error: " << err0 << "\n"
	    << "Estimated RMS force error (based on fixed charges): " << err << " kcal/(mol A)\n"
	    << flush;
      insist(eta > 0);
      Out() << "Ewald screening parameter: " << eta << " A^-1\n" << flush;
    }
  }
  Out() << "\n" << flush;
}

void ElecLJ::types_changed()
{
  types_changed(*msys);
  set_lambda(0);
}

void ElecLJ::set_lambda(double lambda)
{
  const double tmp = 1 - lambda;
  int i;
  HTab<int> se;
  /* Calculate current parameters for each site */
  for (i = 0; i < nsite; i++) {
    Site &s = *site[i];
    s.set_lambda(lambda);
    rscreen[i] = tmp*s.rscreen1 + lambda*s.rscreen2;
    if (!(site_mask[i] & SITE_HAS_LJ))
      continue;
    sig[i] = tmp*s.sig1 + lambda*s.sig2;
    eps[i] = tmp*s.eps1 + lambda*s.eps2;
    char key[128];
    snprintf(key, 128, "%12.6g %12.6g", sig[i], eps[i]);
    if (se.exists(key))
      se[key]++;
    else
      se[key] = 1;
  }
  /* Calculate totA and totB for LJ correction */
  totA = totB = 0;
  HTabIterator<int> ti,tj;
  for (ti.init(se); ti.ok(); ti.next()) {
    double si, ei;
    sscanf(ti.key(), "%lf%lf", &si,&ei);
    for (tj.init(se); tj.ok();  tj.next()) {
      double sj, ej;
      sscanf(tj.key(), "%lf%lf", &sj,&ej);
      const double s6 = pow(geometric_combining_for_sigma ? sqrt(si*sj) : 0.5*(si+sj), 6);
      const double b = 4*sqrt(ei*ej)*ti.val()*tj.val()*s6;
      totB -= b;
      totA += b*s6;
    }
  }
  /* Calculate current bond charge increment params */
  HTabIterator<BondChargeIncrement,Int2> bc;
  for (bc.init(bci); bc.ok(); bc.next()) {
    BondChargeIncrement &b = bc.val();
    b.k = tmp*b.k1 + lambda*b.k2;
    b.q0 = tmp*b.q01 + lambda*b.q02;
  }
  HTabIterator<OneThree,Int2> k;
  for (k.init(one_three); k.ok(); k.next())
    k.val().param = tmp*k.val().p1 + lambda*k.val().p2;
}

const char *ElecLJ::bond_type_for_sites(const Atom *a, int i, int j) const
{
  const AtomSite *ai = dynamic_cast<const AtomSite *> (site[i]);
  const AtomSite *aj = dynamic_cast<const AtomSite *> (site[j]);
  return ai && aj ? a[ai->iatom].bond_type_for_neighbor(aj->iatom) : 0;
}

void ElecLJ::types_changed(const MSys &pert)
{
  TIMESTART("Eleclj::types_changed");
  const int natom = msys->atom.size();
  insist(pert.atom.size() == natom);
  const Atom *at = msys->atom;
  const Atom *pat = pert.atom;
  bool any_new_site = false;
  /* Create new sites, if necessary */
  int i;
  for (i = 0; i < natom; i++) {
    const Atom &a = at[i];
    const Atom &p = pat[i];
    if (!site_on_atom[i].atom) {
      site_on_atom[i].atom = new AtomSite(i);
      any_new_site = true;
    }
    site_on_atom[i].atom->type1 = a.type;
    site_on_atom[i].atom->type2 = p.type;
    if (!any_virtual_sites)
      continue;
    if (a.neighbors.size() == 2) {
      LstPairItr<int> n;
      n.init(a.neighbors);
      const Str3 s1a(at[n.i()].type,a.type,at[n.j()].type);
      const Str3 s2a(at[n.j()].type,a.type,at[n.i()].type);
      const Str3 s1p(pat[n.i()].type,p.type,pat[n.j()].type);
      const Str3 s2p(pat[n.j()].type,p.type,pat[n.i()].type);
      /* M Sites */
      MSiteSpec *msa = msite.get(s1a);
      if (!msa)
	msa = msite.get(s2a);
      MSiteSpec *msp = msite.get(s1p);
      if (!msp)
	msp = msite.get(s2p);
      if (msa || msp) {
	if (!site_on_atom[i].m) {
	  site_on_atom[i].m = new MSite(n.i(),i,n.j());
	  any_new_site = true;
	}
	site_on_atom[i].m->init(msa,msp);
      }
      /* Sp3 sites */
      Sp3SiteSpec *sp3a = sp3.get(s1a);
      if (!sp3a)
	sp3a = sp3.get(s2a);
      Sp3SiteSpec *sp3p = sp3.get(s1p);
      if (!sp3p)
	sp3p = sp3.get(s2p);
      if (sp3a || sp3p) {
	if (!site_on_atom[i].sp3a) {
	  insist(!site_on_atom[i].sp3b);
	  site_on_atom[i].sp3a = new Sp3Site(n.i(),i,n.j(),1);
	  site_on_atom[i].sp3b = new Sp3Site(n.i(),i,n.j(),-1);
	  any_new_site = true;
	}
	site_on_atom[i].sp3a->init(sp3a,sp3p);
	site_on_atom[i].sp3b->init(sp3a,sp3p);
      }
    } else if (a.neighbors.size() == 3) {
      /* Sp2 sites */
      LstItr<int> n;
      int ni, nn[3];
      for (ni = 0, n.init(a.neighbors); n.ok(); ni++, n.next())
	nn[ni] = n();
      for (n.init(a.neighbors); n.ok(); n.next()) {
	Sp2SiteSpec *sp2a = sp2.get(Str2(at[n()].type,a.type));
	Sp2SiteSpec *sp2p = sp2.get(Str2(pat[n()].type,p.type));
	if (sp2a || sp2p) {
	  if (!site_on_atom[n()].sp2a) {
	    insist(!site_on_atom[n()].sp2b);
	    int b,c;
	    if (nn[0] == n()) {
	      b = nn[1];
	      c = nn[2];
	    } else if (nn[1] == n()) {
	      b = nn[0];
	      c = nn[2];
	    } else {
	      b = nn[0];
	      c = nn[1];
	    }
	    site_on_atom[n()].sp2a = new Sp2Site(n(),b,c,i,1);
	    site_on_atom[n()].sp2b = new Sp2Site(n(),b,c,i,-1);
	    any_new_site = true;
	  }
	  site_on_atom[n()].sp2a->init(sp2a,sp2p);
	  site_on_atom[n()].sp2b->init(sp2a,sp2p);
	}
      }
    }
  }
  if (any_new_site) {
    init_sites();
    if (p3m) {
      p3m_t tmp = p3m;
      p3m = p3m_copy(p3m,nsite);
      p3m_delete(tmp);
    }
    if (kspace) {
      struct kspace *tmp = kspace;
      kspace = kspace_new(kspace_cutoff,nsite);
      kspace_delete(tmp);
    }
    if (nlist) {
      NeighborList *tmp = nlist;
      nlist = new NeighborList(ngroup, rc, nlw);
      nlist->set_boundary_conditions(*msys->boundary_conditions);
      delete tmp;
    }
  }
  /* Initialize fixed charge, screening, and LJ params */
  for (i = 0; i < nsite; i++) {
    Site &s = *site[i];
    site_mask[i] = 0;
    s.rezero();
    const double *par;
    if ((par = fixed_charge.get(s.type1)) && !is_almost_zero(*par)) {
      s.q01 = e_to_charge_unit(*par);
      site_mask[i] |= SITE_HAS_Q;
      any_site_has_q = true;
    }
    if ((par = fixed_charge.get(s.type2)) && !is_almost_zero(*par)) {
      s.q02 = e_to_charge_unit(*par);
      site_mask[i] |= SITE_HAS_Q;
      any_site_has_q = true;
    }
    if ((par = screening_radius.get(s.type1)) != 0) {
      if (*par < 0)
	die("ElecLJ::types_changed: screening radius for type %s (%f) is < 0", (const char *) s.type1, *par);
      s.rscreen1 = *par;
    }
    if ((par = screening_radius.get(s.type2)) != 0) {
      if (*par < 0)
	die("ElecLJ::types_changed: screening radius for type %s (%f) is < 0", (const char *) s.type1, *par);
      s.rscreen2 = *par;
    }
    const Str *sljtmp;
    const char *slj;
    double *tsig, *teps;
    sljtmp = general_type.get(s.type1);
    slj = sljtmp ? *sljtmp : s.type1;
    tsig = sigma.get(slj);
    teps = epsilon.get(slj);
    if (tsig && teps) {
      any_site_has_lj = true;
      if (!(tsig && teps))
	die("ElecLJ: sigma or epsilon, but not both, exist for type %s", slj);
      if (is_almost_zero(*teps))
	continue;
      if (*teps < 0)
	die("ElecLJ: epsilon must be >= 0 for type %s", slj);
      if (*tsig <= 0)
	die("ElecLJ: sigma must be > 0 for type %s", slj);
      site_mask[i] |= SITE_HAS_LJ;
      s.sig1 = *tsig;
      s.eps1 = *teps;
    }
    sljtmp = general_type.get(s.type2);
    slj = sljtmp ? *sljtmp : s.type2;
    tsig = sigma.get(slj);
    teps = epsilon.get(slj);
    if (tsig && teps) {
      if (!(tsig && teps))
	die("ElecLJ: sigma or epsilon, but not both, exist for type %s", slj);
      if (is_almost_zero(*teps))
	continue;
      if (*teps < 0)
	die("ElecLJ: epsilon must be >= 0 for type %s", slj);
      if (*tsig <= 0)
	die("ElecLJ: sigma must be > 0 for type %s", slj);
      site_mask[i] |= SITE_HAS_LJ;
      s.sig2 = *tsig;
      s.eps2 = *teps;
    }
  }
  /***
   * Keys for bond charge increment parameters are of the form
   * type1_bondtype_type2
   * for bonds between atoms, with bond types defined.
   * Otherwise they are of the form
   * type1_type2
   ***/
  one_three.remove_all();
  if (any_bci_params)
    for (i = 0; i < nsite; i++) {
      LstItr<int> j;
      for (j.init(neighbors[i]); j.ok(); j.next())
	if (i < j()) {
	  const Str ti1 = site[i]->type1, tj1 = site[j()]->type1;
	  const Str ti2 = site[i]->type2, tj2 = site[j()]->type2;
	  const char *bt1 = bond_type_for_sites(at,i,j());
	  const char *bt2 = bond_type_for_sites(pat,i,j());
	  const Str tij1 = bt1 ? ti1 + "_" + bt1 + "_" + tj1 : ti1 + "_" + tj1;
	  const Str tji1 = bt1 ? tj1 + "_" + bt1 + "_" + ti1 : tj1 + "_" + ti1;
	  const Str tij2 = bt2 ? ti2 + "_" + bt2 + "_" + tj2 : ti2 + "_" + tj2;
	  const Str tji2 = bt2 ? tj2 + "_" + bt2 + "_" + ti2 : tj2 + "_" + ti2;
          const bool same1 = !strcmp(ti1,tj1);
          const bool same2 = !strcmp(ti2,tj2);
	  double s1 = 1, s2 = 1;
	  BondChargeIncrementParam *p1 = bond_charge_increment.get(tij1);
	  if (p1 && p1->k < 0)
	    die("ElecLJ::types_changed: k < 0 for bond charge increment %s", (const char *) tij1);
          if (p1 && same1 && !is_almost_zero(p1->q0))
            die("ElecLJ::types_changed: q0 nonzero for bond charge increment %s", (const char *) tij1);
	  if (!p1) {
	    p1 = bond_charge_increment.get(tji1);
	    if (p1 && p1->k < 0)
	      die("ElecLJ::types_changed: k < 0 for bond charge increment %s", (const char *) tji1);
	    s1 = -1;
	  }
	  BondChargeIncrementParam *p2 = bond_charge_increment.get(tij2);
	  if (p2 && p2->k < 0)
	    die("ElecLJ::types_changed: k < 0 for bond charge increment %s", (const char *) tij2);
          if (p2 && same2 && !is_almost_zero(p2->q0))
            die("ElecLJ::types_changed: q0 nonzero for bond charge increment %s", (const char *) tij2);
	  if (!p2) {
	    p2 = bond_charge_increment.get(tji2);
	    if (p2 && p2->k < 0)
	      die("ElecLJ::types_changed: k < 0 for bond charge increment %s", (const char *) tji2);
	    s2 = -1;
	  }
	  if (p1 || p2) {
	    BondChargeIncrement &b = bci[Int2(i,j())];
	    b.k1 = p1 ? p1->k : 1e8;
	    b.k2 = p2 ? p2->k : 1e8;
	    b.q01 = p1 ? e_to_charge_unit(s1*p1->q0) : 0;
	    b.q02 = p2 ? e_to_charge_unit(s2*p2->q0) : 0;
	    b.omega2 = sq(bond_charge_increment_frequency);
	    b.dq = 0;
	    b.velocity = 0;
	    site_mask[i] |= SITE_HAS_Q;
	    site_mask[j()] |= SITE_HAS_Q;
	    any_site_has_q = true;
	    any_bci = true;
	  }
	}
      LstPairItr<int> n2;
      for (n2.init(neighbors[i]); n2.ok(); n2.next()) {
	double *t1 = one_three_interaction.get(Str3(site[n2.i()]->type1,site[i]->type1,site[n2.j()]->type1));
	if (!t1)
	  t1 = one_three_interaction.get(Str3(site[n2.j()]->type1,site[i]->type1,site[n2.i()]->type1));
	double *t2 = one_three_interaction.get(Str3(site[n2.i()]->type2,site[i]->type2,site[n2.j()]->type2));
	if (!t2)
	  t2 = one_three_interaction.get(Str3(site[n2.j()]->type2,site[i]->type2,site[n2.i()]->type2));
	if (t1 || t2) {
	  one_three[Int2(n2.i(),n2.j())].p1 = t1 ? *t1 : 0;
	  one_three[Int2(n2.i(),n2.j())].p2 = t2 ? *t2 : 0;
	}
      }
    }
  TIMESTOP("Eleclj::types_changed");
}

int Sites::how_many() const
{
  int n = 1;
  insist(atom);
  if (m)
    n++;
  if (sp3a) {
    insist(sp3b);
    n += 2;
  }
  if (sp2a) {
    insist(sp2b);
    n += 2;
  }
  return n;
}

void Sites::add_to(Site **s, int &n) const
{
  insist(atom);
  s[n++] = atom;
  if (m)
    s[n++] = m;
  if (sp3a) {
    insist(sp3b);
    s[n++] = sp3a;
    s[n++] = sp3b;
  }
  if (sp2a) {
    insist(sp2b);
    s[n++] = sp2a;
    s[n++] = sp2b;
  }
}

void Sites::add_virtual_neighbors(const Sites *site_on_atom, Lst<int> *nei) const
{
  if (m) {
    nei[m->index].add(site_on_atom[m->a].atom->index);
    nei[m->index].add(site_on_atom[m->b].atom->index);
    nei[m->index].add(site_on_atom[m->c].atom->index);
    nei[site_on_atom[m->a].atom->index].add(m->index);
    nei[site_on_atom[m->b].atom->index].add(m->index);
    nei[site_on_atom[m->c].atom->index].add(m->index);
  }
  if (sp3a) {
    insist(sp3b);
    nei[atom->index].add(sp3a->index);
    nei[atom->index].add(sp3b->index);
    nei[sp3a->index].add(atom->index);
    nei[sp3b->index].add(atom->index);
  }
  if (sp2a) {
    insist(sp2b);
    nei[atom->index].add(sp2a->index);
    nei[atom->index].add(sp2b->index);
    nei[sp2a->index].add(atom->index);
    nei[sp2b->index].add(atom->index);
  }
}

void ElecLJ::init_sites()
{
  int i;
  const int natom = msys->atom.size();
  const Atom *at = msys->atom;
  nsite = 0;
  for (i = 0; i < natom; i++)
    nsite += site_on_atom[i].how_many();
  site = Vec<Site *>(nsite);
  int n = 0;
  Lst<Int2> gtmp;
  for (i = 0; i < natom; i++) {
    const Lst<int> &nbr = at[i].neighbors;
    if (nbr.size() == 1) { 
      const int j = nbr.first();
      if (at[j].neighbors.size() == 1) {
	/* Dimer */
	insist(at[j].neighbors.first() == i);
	if (i < j) {
	  const int n0 = n;
	  site_on_atom[i].add_to(site,n);
	  site_on_atom[j].add_to(site,n);
	  gtmp.add(Int2(n0,n));
	}
      }
    } else {
      const int n0 = n;
      site_on_atom[i].add_to(site,n);
      LstItr<int> j;
      for (j.init(nbr); j.ok(); j.next())
	if (at[j()].neighbors.size() == 1)
	  site_on_atom[j()].add_to(site,n);
      gtmp.add(Int2(n0,n));
    }
  }
  insist(nsite == n);
  for (i = 0; i < n; i++)
    site[i]->index = i;
  group.copy_from_list(gtmp);
  ngroup = group.size();
  Out() << "Number of sites: " << nsite << "\n"
	<< "Number of groups: " << ngroup << "\n"
	<< flush;
  /* Allocate memory */
  gr = CVec(ngroup);
  q = DVec(nsite, 0.0);
  sig = DVec(nsite);
  eps = DVec(nsite);
  rscreen = DVec(nsite, 0.0);
  phi = DVec(nsite, 0.0);
  evec = CVec(nsite);
  evec.zero();
  r = CVec(nsite);
  lj_force = CVec(nsite);
  lj_force.zero();
  site_mask = Vec<int>(nsite, 0);
  lj_type = Vec<int>(nsite, -1);
  any_site_has_q = false;
  any_site_has_lj = false;
  if (any_bci_params)
    bci = HTab<BondChargeIncrement,Int2>(prime_near(nsite));
  /* Initialize site covalent neighbors */
  neighbors = Vec<Lst<int> >(nsite);
  for (i = 0; i < natom; i++) {
    LstItr<int> j;
    for (j.init(at[i].neighbors); j.ok(); j.next())
      neighbors[site_on_atom[i].atom->index].add(site_on_atom[j()].atom->index);
    site_on_atom[i].add_virtual_neighbors(site_on_atom,neighbors);
  }
  /* Initialize excluded list - 1-2 or 1-3 sites */
  if (excluded)
    PairBinary_delete(excluded);
  excluded = PairBinary_new(nsite);
  LstItr<int> j;
  LstPairItr<int> nn;
  for (i = 0; i < nsite; i++) {
    for (j.init(neighbors[i]); j.ok(); j.next())
      PairBinary_add(excluded, i, j());
    for (nn.init(neighbors[i]); nn.ok(); nn.next())
      PairBinary_add(excluded, nn.i(), nn.j());
  }
  PairBinary_init(excluded);
  /* Initialize 1-4 */
  if (!are_approximately_equal(elec_1_4_scale_factor,1) || 
      !are_approximately_equal(lj_1_4_scale_factor,1)) {
    if (one_four)
      PairBinary_delete(one_four);
    one_four = PairBinary_new(nsite);
    Lst<Int2> bonds;
    LstItr<int> ni, nl;
    for (i = 0; i < nsite; i++) {
      for (ni.init(neighbors[i]); ni.ok(); ni.next())
	if (i < ni())
	  bonds.add(Int2(i,ni()));
    }
    LstItr<Int2> b;
    for (b.init(bonds); b.ok(); b.next())
      for (ni.init(neighbors[b().a]); ni.ok(); ni.next())
	if (ni() != b().b)
	  for (nl.init(neighbors[b().b]); nl.ok(); nl.next())
	    if (nl() != b().a && nl() != ni() && 
		!PairBinary_exists(excluded,ni(),nl()))
	      PairBinary_add(one_four, ni(), nl());
    PairBinary_init(one_four);
  }
  /* Initialize excluded groups (i.e. groups containing
     sites that are excluded) */
  group_for_site = Vec<int>(nsite, -1);
  for (i = 0; i < ngroup; i++) {
    int k;
    for (k = group[i].a; k < group[i].b; k++)
      group_for_site[k] = i;
  }
  for (i = 0; i < nsite; i++)
    assert(group_for_site[i] > -1);
  group_excluded = PairBinary_new(ngroup);
  for (i = 0; i < nsite; i++) {
    const int ig = group_for_site[i];
    const int *k = PairBinary_elements(excluded)[i];
    const int *kend = &k[PairBinary_number_of_elements(excluded)[i]];
    while (k < kend) {
      const int jg = group_for_site[*k];
      if (ig != jg)
	PairBinary_add(group_excluded,ig,jg);
      k++;
    }
  }
  PairBinary_init(group_excluded);
}

void ElecLJ::coordinates_changed()
{
  int i;
  for (i = 0; i < nsite; i++)
    r[i] = site[i]->position(msys->atom);
  for (i = 0; i < ngroup; i++)
    gr[i] = r[group[i].a];
  if (nlist)
    nlist->coordinates_changed(gr);
  update_q();
}

void ElecLJ::add_to_energy_and_forces(double &u, Cartesian *f)
{
  if (!any_site_has_q && !any_site_has_lj)
    return;
  TIMESTART("ElecLJ::add_to_energy_and_forces");
  coordinates_changed();
  if (any_bci && !(dyn && have_solved_for_bci))
    solve_for_bci(Cartesian(0,0,0));
  calc_phi_evec_lj();
  for (int i = 0; i < nsite; i++)
    site[i]->apply_force(msys->atom, f, q[i]*evec[i] + lj_force[i]);
  elec_energy = 0;
  if (any_site_has_q) {
    elec_energy = (q*phi)/2;
    if (any_bci) {
      HTabIterator<BondChargeIncrement,Int2> b;
      for (b.init(bci); b.ok(); b.next())
	elec_energy += b.val().self_energy();
    }
  }
  u += elec_energy + lj_energy;
  TIMESTOP("ElecLJ::add_to_energy_and_forces");
}

void ElecLJ::update_q()
{
  int i;
  for (i = 0; i < nsite; i++)
    q[i] = site[i]->q0;
  if (any_bci) {
    HTabIterator<BondChargeIncrement,Int2> b;
    for (b.init(bci); b.ok(); b.next()) {
      q[b.key().a] += b.val().q0 + b.val().dq;
      q[b.key().b] -= b.val().q0 + b.val().dq;
    }
  }
}

void ElecLJ::solve_for_bci(const Cartesian &external_field)
{
  int it;
  bool converged = false;
  HTabIterator<BondChargeIncrement,Int2> b;
  for (it = 0; !converged && it < max_iterations; it++) {
    converged = true;
    calc_phi();
    if (external_field.magnitude() > 0)
      for (int i = 0; i < nsite; i++)
	phi[i] -= external_field * r[i];
    for (b.init(bci); b.ok(); b.next()) {
      double &bq = b.val().dq;
      if (is_almost_zero(b.val().k)) {
	bq = 0;
	continue;
      }
      double last = bq;
      bq = (phi[b.key().b] - phi[b.key().a])/b.val().k;
      if (is_not_a_number(bq))
	die("ElecLJ::solve_for_bci: NaN detected");
      if (fabs(bq) > tolerance) {
        if (fabs(last-bq)/fabs(bq) > tolerance)
          converged = false;
      } else {
        if (fabs(last-bq) > tolerance)
          converged = false;
      }
    }
    update_q();
  }
  if (!converged)
    die("ElecLJ::solve_for_bci: did not converge to tolerance of %f within %d iterations", 
	tolerance, max_iterations); 
  if (verbose >= 4)
    Out() << "Converged to tolerance of " << tolerance << " within "
	  << it << " iterations.\n" << flush;
  for (b.init(bci); b.ok(); b.next())
    b.val().velocity = 0;
  have_solved_for_bci = true;
}

void ElecLJ::integrate_positions(double dt)
{
  HTabIterator<BondChargeIncrement,Int2> b;
  for (b.init(bci); b.ok(); b.next())
    b.val().dq += b.val().velocity * dt;
}

void ElecLJ::integrate_velocities(double dt)
{
  HTabIterator<BondChargeIncrement,Int2> b;
  for (b.init(bci); b.ok(); b.next())
    if (is_almost_zero(b.val().k))
      b.val().velocity = 0;
    else
      b.val().velocity += (phi[b.key().b] - phi[b.key().a] - b.val().k * b.val().dq)/b.val().mass() * dt;
}

void ElecLJ::zero_bci()
{
  HTabIterator<BondChargeIncrement,Int2> b;
  for (b.init(bci); b.ok(); b.next())
    b.val().dq = 0;
}

Tensor ElecLJ::polarizability() const
{
  if (!any_bci)
    return Tensor(0);
  ElecLJ &that = *(ElecLJ *) this;
  that.zero_bci();
  that.solve_for_bci(Cartesian(0,0,0));
  const Cartesian mu0 = dipole_moment();
  Tensor a;
  that.zero_bci();
  that.solve_for_bci(Cartesian(1,0,0));
  a.col(0) = dipole_moment() - mu0;
  that.zero_bci();
  that.solve_for_bci(Cartesian(0,1,0));
  a.col(1) = dipole_moment() - mu0;
  that.zero_bci();
  that.solve_for_bci(Cartesian(0,0,1));
  a.col(2) = dipole_moment() - mu0;
  that.zero_bci();
  that.solve_for_bci(Cartesian(0,0,0));
  return a;
}

void ElecLJ::calc_phi_evec_lj()
{
  TIMESTART("ElecLJ::calc_phi_evec_lj");
  const BoundaryConditions &bc = *msys->boundary_conditions;
  evec.zero();
  phi.zero();
  lj_force.zero();
  lj_energy = 0;
  const Tensor a = bc.lattice_vectors();
  if (include_rspace) {
    const Int2 *g = group;
    TIMESTART("ElecLJ::RSpace");
    if (bc.type == non_periodic)
      RSpaceNonPeriodic(nsite, ngroup, (const int (*)[2]) g, r,
			group_excluded, excluded,
			site_mask, q, phi, evec, rscreen, 
			sig, eps, geometric_combining_for_sigma,
			&lj_energy, lj_force,
			one_four, elec_1_4_scale_factor, lj_1_4_scale_factor);
    else
      RSpaceEwald(nsite, ngroup, (const int (*)[2]) g, r, bc.type, (const tensor_t *) &a,
		  nlist->neighbors(), nlist->number_of_neighbors(),
		  group_excluded, excluded,
		  sq(rc-smoothing_width), sq(rc), 
		  site_mask, q, phi, evec, eta, rspace_chebyshev, rscreen,
		  sig, eps, geometric_combining_for_sigma,
		  &lj_energy, lj_force,
		  one_four, elec_1_4_scale_factor, lj_1_4_scale_factor);
    TIMESTOP("ElecLJ::RSpace");
  }
  if (bc.type != non_periodic && any_site_has_q) {
    if (include_kspace) {
      if (use_p3m) {
	TIMESTART("p3m_field");
	p3m_field(p3m, r, q, phi, evec);
	TIMESTOP("p3m_field");
      } else {
	const Tensor latv = bc.reciprocal_lattice_vectors();
	kspace_field(kspace, r, q, (const tensor_t *) &latv, bc.volume(), eta, phi, evec);
      }
    }
    if (include_self)
      calc_self_phi();
  }
  calc_phi_one_three();
  recalculate_lj_correction();
  TIMESTOP("ElecLJ::calc_phi_evec_lj");
}

void ElecLJ::calc_self_phi()
{
  const double a = eta * M_2_SQRTPI;
  int i;
  double qsum = 0;
  for (i = 0; i < nsite; i++) {
    qsum += q[i];
    phi[i] -= a * q[i];
  }
  if (!is_almost_zero(qsum)) {
    /* Boguzs, Cheatham, Brooks JCP 108, 7070 (1998) eq 3 */
    const double b = -M_PI/(sq(eta)*msys->boundary_conditions->volume()) * qsum;
    for (i = 0; i < nsite; i++)
      phi[i] += b;
  }
}

void ElecLJ::calc_phi()
{
  const BoundaryConditions &bc = *msys->boundary_conditions;
  const Tensor a = bc.lattice_vectors();
  phi.zero();
  if (include_rspace) {
    const Int2 *g = group;
    TIMESTART("ElecLJ::RSpacePhi");
    if (bc.type == non_periodic)
      RSpaceNonPeriodicPhi(nsite, ngroup, (const int (*)[2]) g, r,
			   group_excluded, excluded, site_mask, 
			   q, phi, rscreen, one_four, elec_1_4_scale_factor);
    else {
      RSpaceEwaldPhi(nsite, ngroup, (const int (*)[2]) g, r,
		     bc.type, (const tensor_t *) &a,
		     nlist->neighbors(), nlist->number_of_neighbors(),
		     group_excluded, excluded,
		     sq(rc-smoothing_width), sq(rc), site_mask, q, phi,
		     eta, rspace_chebyshev, rscreen, one_four, elec_1_4_scale_factor);
    }
    TIMESTOP("ElecLJ::RSpacePhi");
  }
  if (bc.type != non_periodic && any_site_has_q) {
    if (include_kspace) {
      if (use_p3m) {
	p3m_field(p3m, r, q, phi, 0);
      } else {
	const Tensor latv = bc.reciprocal_lattice_vectors();
	kspace_phi(kspace, r, q, (const tensor_t *) &latv, bc.volume(), eta, phi);
      }
    }
    if (include_self)
      calc_self_phi();
  }
  calc_phi_one_three();
}

void ElecLJ::calc_phi_one_three()
{
  HTabIterator<OneThree,Int2> k;
  for (k.init(one_three); k.ok(); k.next()) {
    const int i = k.key().a;
    const int j = k.key().b;
    phi[i] += q[j] * k.val().param;
    phi[j] += q[i] * k.val().param;
  }
}

void ElecLJ::write() const
{
  const BoundaryConditions &bc = *msys->boundary_conditions;
  if (verbose) {
    if (any_site_has_q)
      Out() << "Electrostatic energy: " << elec_energy << " kcal/mol\n";
    if (any_site_has_lj) {
      Out() << "LJ energy: " << lj_energy << " kcal/mol\n";
      if (include_lj_correction)
	Out() << "LJ correction: " << lj_correction << " kcal/mol\n";
    }
    if (any_bci) {
      double ke = 0;
      ConstHTabIterator<BondChargeIncrement,Int2> b;
      int nbci = 0;
      for (b.init(bci); b.ok(); b.next()) {
	ke += b.val().kinetic_energy();
	nbci++;
      }
      Out() << "Extended Lagrangian kinetic energy: " << ke << " kcal/mol\n"
	    << "Extended Lagrangian temperature: " 
	    << (nbci > 0 ? kcal_mol_to_K(2*ke/nbci) : 0.0) << " K\n";
    }
    if (any_site_has_q && bc.type != non_periodic)
      Out() << "Static dielectric constant at " << kcal_mol_to_K(kT) << " K: " 
	    << static_dielectric_constant() << "\n";
  }
  if (verbose >= 2 && bc.type != non_periodic)
    Out() << "Optical dielectric constant: " << optical_dielectric_constant() << "\n";
  if (verbose >= 3) {
    const Tensor a = polarizability();
    Tensor evec;
    Cartesian eval;
    Diagonalize(a,eval,evec);
    Out() << "Polarizability (A^3): " << a << eval << "\n";
  }
  if (verbose >= 4) {
    OutPrintf("%-5s %-10s %12s %12s %12s %12s %8s %8s\n","sym","type","x","y","z","charge (e)","sig","eps");
    int i;
    for (i = 0; i < nsite; i++) {
      const char *sym = " ";
      if (dynamic_cast<AtomSite *>(site[i]))
	sym = msys->atom[(dynamic_cast<AtomSite *>(site[i]))->iatom].symbol;
      else if (dynamic_cast<MSite *>(site[i]))
	sym = "M";
      else if (dynamic_cast<Sp3Site *>(site[i]))
	sym = "sp3";
      else if (dynamic_cast<Sp2Site *>(site[i]))
	sym = "sp2";
      OutPrintf("%-5s %-10s ", sym, (const char *) site[i]->type1);
      OutPrintf("%12.8f %12.8f %12.8f ", r[i].x, r[i].y, r[i].z);
      if (site_mask[i] & SITE_HAS_Q)
	OutPrintf("%12.8f ", charge_unit_to_e(q[i]));
      else
	OutPrintf("%12s ", " ");
      if (site_mask[i] & SITE_HAS_LJ)
	OutPrintf("%8.4f %8.4f ", sig[i], eps[i]);
      else
	OutPrintf("%8s %8s ", " ", " ");
      Out() << "\n";
    }
    Out() << "\n";
    if (any_bci_params) {
      Out() << "Bond-charge increments:\n";
      OutPrintf("%5s %-8s %-4s %5s %-8s %12s %12s %12s\n", 
                "i", "type","bond", "j", "type","q0 (e)","k","dq (e)");
      for (i = 0; i < nsite; i++) {
	LstItr<int> j;
	for (j.init(neighbors[i]); j.ok(); j.next())
	  if (i < j()) {
	    const char *p = bond_type_for_sites(msys->atom,i,j());
	    OutPrintf("%5d %-8s %-4s %5d %-8s", 
		      i, (const char *) site[i]->type1, p ? p : "",
		      j(), (const char *) site[j()]->type1);
            const BondChargeIncrement *b = bci.get(Int2(i,j()));
            if (b)
              OutPrintf("%12.8f %12.8f %12.8f", 
                        charge_unit_to_e(b->q0), 
                        b->k, charge_unit_to_e(b->dq));
            Out() << "\n";
	  }
      }
      Out() << "\n";
    }
    const Cartesian mu = dipole_moment().map(dipole_unit_to_debye);
    Out() << "Net charge: " << charge_unit_to_e(net_charge()) << " e\n"
	  << "Dipole moment (D): " << mu << "\n"
	  << "Dipole moment magnitude: " << mu.magnitude() << " D\n";
    Quadrupole Q = quadrupole_moment().map(dipole_unit_to_debye);
    Out() << "Quadrupole moment (D A):\n" << Q << "\n";
    Q.make_traceless();
    Out()  << "Traceless quadrupole moment (D A): " << Q << "\n";
  }
  if (verbose >= 5) {
    Octopole O = octopole_moment().map(dipole_unit_to_debye);
    Out() << "Octopole moment (D A):\n" << O << "\n";
    O.make_traceless();
    Out() << "Traceless octopole moment (D A):\n" << O << "\n";
    Hexadecapole H = hexadecapole_moment().map(dipole_unit_to_debye);
    Out() << "Hexadecapole moment (D A):\n" << H << "\n";
    H.make_traceless();
    Out() << "Traceless hexadecapole moment (D A):\n" << H << "\n";
    H.make_traceless();
  }
  if (verbose >= 6) {
    Out() << *this;
    Out() << "sites:\n";
    for (int i = 0; i < nsite; i++)
      site[i]->write(Out());
    Out() << "excluded:\n" << flush;
    PairBinary_show(excluded);
    Out() << "group_excluded:\n" << flush;
    PairBinary_show(group_excluded);
    if (one_four) {
      Out() << "one_four:\n" << flush;
      PairBinary_show(one_four);
    }
    if (nlist) {
      Out() << "Neighbor list:\n" << flush;
      nlist->show_self();
    }
    if (kspace)
      kspace_show(kspace);
  }
}

double ElecLJ::static_dielectric_constant() const
{
  return 1 + 4*M_PI*dipole_moment().sq() / 
    (3*msys->boundary_conditions->volume()*kT);
}

double ElecLJ::optical_dielectric_constant() const
{
  return 1 + 4*M_PI*polarizability().trace() / 
    (3*msys->boundary_conditions->volume());
}

Cartesian ElecLJ::dipole_moment() const
{
  Cartesian mu(0,0,0);
  if (any_site_has_q)
    for (int i = 0; i < nsite; i++)
      mu += q[i]*r[i];
  return mu;
}

/* Same convention as Gaussian */
Quadrupole ElecLJ::quadrupole_moment() const
{
  Quadrupole Q;
  Q.zero();
  if (any_site_has_q)
    for (int i = 0; i < nsite; i++) {
      const Cartesian s = r[i];
      Q.XX += q[i] * s.x * s.x;
      Q.YY += q[i] * s.y * s.y;
      Q.ZZ += q[i] * s.z * s.z;
      Q.XY += q[i] * s.x * s.y;
      Q.XZ += q[i] * s.x * s.z;
      Q.YZ += q[i] * s.y * s.z;
    }
  return Q;
}

Octopole ElecLJ::octopole_moment() const
{
  Octopole o;
  o.zero();
  if (any_site_has_q)
    for (int i = 0; i < nsite; i++) {
      const Cartesian s = r[i];
      o.XXX += q[i] * s.x * s.x * s.x;
      o.YYY += q[i] * s.y * s.y * s.y;
      o.ZZZ += q[i] * s.z * s.z * s.z;
      o.XYY += q[i] * s.x * s.y * s.y;
      o.XXY += q[i] * s.x * s.x * s.y;
      o.XXZ += q[i] * s.x * s.x * s.z;
      o.XZZ += q[i] * s.x * s.z * s.z;
      o.YZZ += q[i] * s.y * s.z * s.z;
      o.YYZ += q[i] * s.y * s.y * s.z;
      o.XYZ += q[i] * s.x * s.y * s.z;
    }
  return o;
}

Hexadecapole ElecLJ::hexadecapole_moment() const
{
  Hexadecapole h;
  h.zero();
  if (any_site_has_q)
    for (int i = 0; i < nsite; i++) {
      const Cartesian s = r[i];
      h.XXXX += q[i] * s.x * s.x * s.x * s.x;
      h.YYYY += q[i] * s.y * s.y * s.y * s.y;
      h.ZZZZ += q[i] * s.z * s.z * s.z * s.z;
      h.XXXY += q[i] * s.x * s.x * s.x * s.y;
      h.XXXZ += q[i] * s.x * s.x * s.x * s.z;
      h.YYYX += q[i] * s.y * s.y * s.y * s.x;
      h.YYYZ += q[i] * s.y * s.y * s.y * s.z;
      h.ZZZX += q[i] * s.z * s.z * s.z * s.x;
      h.ZZZY += q[i] * s.z * s.z * s.z * s.y;
      h.XXYY += q[i] * s.x * s.x * s.y * s.y;
      h.XXZZ += q[i] * s.x * s.x * s.z * s.z;
      h.YYZZ += q[i] * s.y * s.y * s.z * s.z;
      h.XXYZ += q[i] * s.x * s.x * s.y * s.z;
      h.YYXZ += q[i] * s.y * s.y * s.x * s.z;
      h.ZZXY += q[i] * s.z * s.z * s.x * s.y;
    }
  return h;
}

double ElecLJ::net_charge() const
{
  double totq = 0;
  if (any_site_has_q)
    for (int i = 0; i < nsite; i++)
    totq += q[i];
  return totq;
}

void ElecLJ::recalculate_lj_correction()
{
  if (!include_lj_correction || msys->boundary_conditions->type == non_periodic) {
    lj_correction = 0;
    return;
  }
  const double rhi = rc;
  const double rlo = rhi - smoothing_width;
  const double rlo2 = rlo*rlo, rlo3 = rlo2*rlo, rlo4 = rlo3*rlo;
  const double rlo5 = rlo4*rlo, rlo6 = rlo5*rlo, rlo7 = rlo6*rlo, rlo9 = rlo7*rlo2;
  const double rhi2 = rhi*rhi, rhi3 = rhi2*rhi, rhi4 = rhi3*rhi, rhi5 = rhi4*rhi;
  const double rhi6 = rhi5*rhi, rhi7 = rhi6*rhi;
  lj_correction = 2*M_PI/msys->boundary_conditions->volume() *
    (totA/(9*rlo9) + totB/(3*rlo3) -
     (rhi - rlo)/(63*rhi3*rlo9*pow(rhi + rlo,5)) *
     (3*totB*rhi3*rlo6*(7*rhi4 + 42*rhi3*rlo + 112*rhi2*rlo2 + 150*rhi*rlo3 + 25*rlo4) + 
      totA*(7*rhi7 + 42*rhi6*rlo + 112*rhi5*rlo2 + 
	    182*rhi4*rlo3 + 217*rhi3*rlo4 + 224*rhi2*rlo5 + 192*rhi*rlo6 + 32*rlo7)));
  lj_energy += lj_correction;
}

void ElecLJ::scale(double s)
{
  if (p3m)
    p3m_scale(p3m,s);
  if (include_rspace && any_site_has_q) {
    eta *= 1/s;
    if (rspace_chebyshev)
      free(rspace_chebyshev);
    rspace_chebyshev = RSpaceChebyshev_new(rc, eta, 0, 0);
  }
  if (nlist)
    nlist->set_boundary_conditions(*msys->boundary_conditions);
}

void ElecLJ::write_electrostatic_potential_at_gridpoints(const CVec gp)
{
  DVec p;
  electrostatic_potential_at_gridpoints(gp,p);
  OutPrintf("Electrostatic potential:\n");
  OutPrintf("%15s %15s %15s %15s\n", "x (A)" , "y (A)", "z (A)", "ESP (a.u.)");
  OutPrintf("---------------------------------------------------------------\n");
  for (int i = 0; i < gp.size(); i++)
    OutPrintf("%15.8f %15.8f %15.8f %15.8f\n", gp[i].x, gp[i].y, gp[i].z, 
	      potential_unit_to_hartree_e(p[i]));
  OutPrintf("\n\n");
}

void ElecLJ::electrostatic_potential_at_gridpoints(const CVec gp, DVec &p)
{
  int i;
  const BoundaryConditions &bc = *msys->boundary_conditions;
  coordinates_changed();
  if (any_bci)
    solve_for_bci(Cartesian(0,0,0));
  const int ngrid = gp.size();
  const int ntot = ngrid+nsite;
  p.resize(ntot);
  p.zero();
  if (bc.type == non_periodic) {
    RSpaceNonPeriodicPhiAtGridpoints(nsite, ngrid, r, gp, site_mask, q, p, rscreen);
  } else {
    DVec qtot(ngrid+nsite, 0.0);
    CVec rtot(ngrid+nsite);
    for (i = 0; i < ngrid; i++)
      rtot[i] = gp[i];
    for (i = 0; i < nsite; i++) {
      qtot[i+ngrid] = q[i];
      rtot[i+ngrid] = r[i];
    }
    if (include_rspace) {
      const Tensor a = bc.lattice_vectors();
      RSpaceEwaldPhiAtGridpoints(nsite, ngrid, r, gp, bc.type, (const tensor_t *) &a,
				 sq(rc-smoothing_width), sq(rc), site_mask,
				 q, p, eta, rspace_chebyshev, rscreen);
    }
    if (include_kspace) {
      struct kspace *kspace_tmp = kspace_new(kspace_cutoff, ntot);
      const Tensor a = bc.reciprocal_lattice_vectors();
      kspace_phi(kspace_tmp, rtot, qtot, (const tensor_t *) &a, bc.volume(), eta, p);
      kspace_delete(kspace_tmp);
    }
    if (include_self) {
      double qsum = 0;
      for (i = 0; i < nsite; i++)
	qsum += q[i];
      if (!is_almost_zero(qsum)) {
	/* Boguzs, Cheatham, Brooks JCP 108, 7070 (1998) eq 3 */
	const double b = -M_PI/(sq(eta)*msys->boundary_conditions->volume()) * qsum;
	for (i = 0; i < ngrid; i++)
	  p[i] += b;
      }
    }
  }
}

double ElecLJ::charge_on_atom(int i) const
{
  return q[site_on_atom[i].atom->index];
}

/* Sk = sum q exp(-ikr) */
Complex ElecLJ::structure_factor(const Cartesian &k) const
{
  Complex s(0,0);
  for (int i = 0; i < nsite; i++)
    s += q[i]*expi(-k*r[i]);
  return s;
}
