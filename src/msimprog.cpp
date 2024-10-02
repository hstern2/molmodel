#include <cstdlib>
#include <sstream>
#include <cstring>
#include "msimprog.hpp"
#include "mdyn.hpp"
#include "checkgr.hpp"
#include "random.h"
#include "mmin.hpp"
#include "prepro.hpp"
#include "units.h"
#include "geograd.hpp"
#include "linalg.hpp"
#include "waterintra.hpp"
#include "covalent.hpp"
#include "contact.hpp"
#include "eleclj.hpp"
#include "intra.hpp"
#include "harmonic.hpp"
#include "espfit.hpp"
#include "bcifit.hpp"
#include "ppot.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
#endif
  OutInit();
  const char *dir = getenv("MSIM_DIR");
  if (!dir || strlen(dir) == 0)
    die("msim: must set environment variable MSIM_DIR");
  if (strstr(dir,"/param"))
    AddToSearchPath(getenv("MSIM_DIR"));
  else
    AddToSearchPath(Str(dir)+"/param");
  PreprocessedStream s;
  try {
    s.init(argc,argv);
    MSimProgram msim;
    s >> msim;
  } catch (FatalException) {
    s.write_current_file_and_line(Out());
  }
  Atom::finalize();
  OutFinalize();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}

FEP::FEP() : 
  start_lambda(0), 
  end_lambda(0), 
  perturbed_msys(0), 
  perturbed_potential(0)
{ }

FEP::~FEP()
{
  if (perturbed_msys)
    delete perturbed_msys;
  if (perturbed_potential)
    delete perturbed_potential;
  perturbed_msys = 0;
  perturbed_potential = 0;
}

void FEP::init_perturbed(MSys *msys, MSimPotential *p, Constraint *c)
{
  Out() << "Performing free energy perturbation, start lambda = " 
	<< start_lambda << ", "
	<< "end lambda = " << end_lambda
	<< "\n";
  if (strlen(modify) == 0)
    die("MSimProgram::free_energy_perturbation: expecting 'modify'");
  perturbed_msys = new MSys(*msys);
  istringstream ss((const char *) modify);
  ss >> *perturbed_msys;
  insist(msys->has_same_topology_as(*perturbed_msys));
  p->init(msys);
  p->types_changed(*perturbed_msys);
  p->set_lambda(start_lambda);
  perturbed_potential = new MSimPotential(*p);
  perturbed_potential->init(msys);
  perturbed_potential->types_changed(*perturbed_msys);
  if (c) {
    c->init(msys);
    c->add_constraints(*perturbed_msys);
  }
  f.resize(msys->atom.size());
}

void FEP::update(const MSys &)
{
  ustart = 0;
  perturbed_potential->set_lambda(start_lambda);
  perturbed_potential->add_to_energy_and_forces(ustart,f);
  uend = 0;
  perturbed_potential->set_lambda(end_lambda);
  perturbed_potential->add_to_energy_and_forces(uend,f);
  dU = uend - ustart;
}

void FEP::write() const
{
  Out() << "FEP potential energy for lambda=" << start_lambda << ": " << ustart << " kcal/mol\n"
	<< "FEP potential energy for lambda=" << end_lambda << ": " << uend << " kcal/mol\n"
	<< "FEP energy difference: " << dU << " kcal/mol\n";
}

void MSimProgram::dynamics(istream &s)
{
  insist(!(free_energy_perturbation && weak_couple));
  if (free_energy_perturbation) {
    free_energy_perturbation->init_perturbed(&molsys,potential,constraint);
    callback.push(CallbackPtr(free_energy_perturbation));
  } else {
    potential->init(&molsys);
    if (constraint)
      constraint->init(&molsys);
  }
  MDyn dyn(*potential,constraint,callback);
  if (weak_couple) {
    weak_couple->wc_init(potential->eleclj_, &dyn);
    callback.push(CallbackPtr(weak_couple));
  }
  s >> dyn;
  if (free_energy_perturbation) {
    callback.pop();
    delete free_energy_perturbation;
    free_energy_perturbation = 0;
  }
  if (weak_couple) {
    callback.pop();
    delete weak_couple;
    weak_couple = 0;
  }
}

void MSimProgram::minimize(istream &s)
{
  potential->init(&molsys);
  if (constraint)
    constraint->init(&molsys);
  MMin min(molsys,*potential,constraint);
  s >> min;
}

void MSimProgram::check_gradients(istream &s)
{
  potential->init(&molsys);
  if (constraint)
    constraint->init(&molsys);
  CheckGradients c(molsys,*potential,constraint);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> c;
  c.run();
}

void MSimProgram::energy()
{
  potential->init(&molsys);
  double u = 0;
  CVec f(molsys.atom.size());
  f.zero();
  potential->add_to_energy_and_forces(u,f);
  Out() << "Boundary conditions: " << molsys.boundary_conditions << "\n"
	<< "Energy: " << u << " kcal/mol\n";
  const int nmol = molsys.molecule.size();
  if (nmol > 0)
    Out() << "Energy per molecule: " << u/nmol << " kcal/mol\n";
  potential->write();
  Out() << "\n" << flush;
}

void MSimProgram::energies_from_xyz_trajectory(Str fname)
{
  potential->init(&molsys);
  CVec f(molsys.atom.size());
  FILE *fxyz = FileSearch(fname);
  Out() << "Energies from xyz trajectory file '" << fname << "':\n";
  while (!end_of_file(fxyz)) {
    molsys.read_as_xyz(fxyz);
    double u = 0;
    potential->add_to_energy_and_forces(u,f);
    Out() << "Boundary conditions: " << molsys.boundary_conditions << "\n"
	  << "Energy: " << u << " kcal/mol\n";
    potential->write();
    Out() << "\n";
  }
  Out() << "\n" << flush;
}

void MSimProgram::zmatrices_from_xyz_trajectory(ZMatrix zmat, Str fname)
{
  CVec f(molsys.atom.size());
  FILE *fxyz = FileSearch(fname);
  Out() << "Z-matrices from xyz trajectory file '" << fname << "':\n";
  while (!end_of_file(fxyz)) {
    molsys.read_as_xyz(fxyz);
    molsys.show_as_zmatrix(zmat);
  }
  Out() << "\n" << flush;
}

void MSimProgram::forces()
{
  potential->init(&molsys);
  double u = 0;
  CVec f(molsys.atom.size());
  f.zero();
  potential->add_to_energy_and_forces(u,f);
  potential->write();
  Out() << "Total energy: " << u << " kcal/mol\n" 
	<< "Total energy per molecule: " << u/molsys.molecule.size() 
	<< " kcal/mol\n" 
	<< "Forces: " << f
	<< "Sum of forces: " << f.sum() << "\n\n"
	<< flush;
}

void MSimProgram::binding_energy()
{
  potential->init(&molsys);
  double u = 0;
  CVec f(molsys.atom.size());
  potential->add_to_energy_and_forces(u,f);
  potential->write();
  Out() << "Energy of oligomer: " << u << " kcal/mol\n";
  for (int i = 0; i < molsys.molecule.size(); i++) {
    MSys m;
    m.create_from_molecule(molsys,i);
    potential->init(&m);
    double u0 = 0;
    potential->add_to_energy_and_forces(u0,f);
    potential->write();
    Out() << "Energy of monomer " << i << ": " << u0 << " kcal/mol\n";
    u -= u0;
  }
  for (int i = 0; i < molsys.molecule.size(); i++)
    for (int j = i+1; j < molsys.molecule.size(); j++) {
      MSys m;
      m.create_from_molecule(molsys,i);
      MSys mtmp;
      mtmp.create_from_molecule(molsys,j);
      m.append(mtmp);
      potential->init(&m);
      double utmp = 0;
      potential->add_to_energy_and_forces(utmp,f);
      potential->write();
      Out() << "Energy of dimer " << i << "," << j << ": " << utmp << " kcal/mol\n\n";
    }
  Out() << "Binding energy: " << u << " kcal/mol\n\n";
}

void MSimProgram::finite_difference_forces()
{
  const double eps = 1e-5;
  potential->init(&molsys);
  double u0 = 0;
  CVec f(molsys.atom.size());
  f.zero();
  potential->add_to_energy_and_forces(u0,f);
  CVec fsave = f.copy();
  for (int i = 0; i < molsys.atom.size(); i++) {
    double u;
    Cartesian ftmp;
    molsys.atom[i].position.x += eps;
    u = 0;
    potential->add_to_energy_and_forces(u,f);
    ftmp.x = (u0-u)/eps;
    molsys.atom[i].position.x -= eps;
    molsys.atom[i].position.y += eps;
    u = 0;
    potential->add_to_energy_and_forces(u,f);
    ftmp.y = (u0-u)/eps;
    molsys.atom[i].position.y -= eps;
    molsys.atom[i].position.z += eps;
    u = 0;
    potential->add_to_energy_and_forces(u,f);
    ftmp.z = (u0-u)/eps;
    molsys.atom[i].position.z -= eps;
    OutPrintf("%15.8f %15.8f %15.8f\n", fsave[i].x, fsave[i].y, fsave[i].z);
    OutPrintf("%15.8f %15.8f %15.8f\n\n", ftmp.x, ftmp.y, ftmp.z);
    Out() << flush;
  }
}

void MSimProgram::random_seed(int n)
{
  Out() << "Setting random seed to " << n << "\n\n";
  RandomSeed(n);
}

void MSimProgram::rms_force_error(istream &s)
{
  MSimPotential *p1 = 0, *p2 = 0;
  int ntrial;
  s >> ntrial;
  if (!s)
    die("MSimProgram::rms_force_error: expecting number of trials, then two potentials");
  if (ntrial <= 0)
    die("MSimProgram::rms_force_error: number of trials must be > 0");
  s >> p1 >> p2;
  if (!s)
    die("MSimProgram::rms_force_error: expecting number of trials, then two potentials");
  p1->init(&molsys);
  p2->init(&molsys);
  const int n = molsys.atom.size();
  if (n <= 0)
    die("MSimProgram::rms_force_error: number of atoms must be > 0");
  CVec f1(n), f2(n);
  double frms = 0;
  for (int i = 0; i < ntrial; i++) {
    if (ntrial > 1)
      molsys.randomize_positions();
    double u1 = 0;
    f1.zero();
    p1->add_to_energy_and_forces(u1,f1);
    if (ntrial == 1)
      Out() << "Potential energy 1: " << u1 << " kcal/mol\n";
    double u2 = 0;
    f2.zero();
    p2->add_to_energy_and_forces(u2,f2);
    if (ntrial == 1)
      Out() << "Potential energy 2: " << u2 << " kcal/mol\n";
    frms += (f1 -= f2).sq()/n;
  }
  frms = sqrt(frms/ntrial);
  Out() << "RMS error in forces: " << frms << " kcal/(mol A)\n";
}

void MSimProgram::rms_electric_field_error(istream &s)
{
  int n;
  s >> n;
  if (!s || n < 2)
    die("MSimProgram::rms_electric_field: expecting n >= 2");
  ElecLJ e1, e2;
  s >> e1 >> e2;
  if (!s)
    die("MSimProgram::rms_electric_field: expecting two sets of ElecLJ params");
  MSys m;
  LstStr ls;
  ls.add("H");
  int i;
  for (i = 0; i < n-1; i++)
    ls.add(";Lj");
  m.create_from_topology(ls);
  delete m.boundary_conditions;
  m.boundary_conditions = BoundaryConditions::new_cubic(1);
  e1.fixed_charge["H"] = charge_unit_to_e(1.0);
  e2.fixed_charge["H"] = charge_unit_to_e(1.0);
  e1.init(&m);
  e2.init(&m);
  double rms = 0;
  CVec f(n);
  double u;
  for (i = 0; i < n; i++) {
    m.randomize_positions();
    e1.add_to_energy_and_forces(u,f);
    e2.add_to_energy_and_forces(u,f);
    for (int j = 1; j < n; j++)
      rms += (e1.electric_field()[j] - e2.electric_field()[j]).sq();
  }
  Out() << "RMS error in electric field per unit charge: " 
	<< sqrt(rms/(n*(n-1))) << " kcal/(mol A)"
	<< "\n";
}

void MSimProgram::test_scale(double a)
{
  if (a <= 0)
    die("MSimProgram::test_scale: expecting scale factor > 0");
  potential->init(&molsys);
  CVec f(molsys.atom.size());
  double u = 0;
  f.zero();
  potential->add_to_energy_and_forces(u,f);
  Out() << "Boundary conditions: " << molsys.boundary_conditions << "\n"
	<< "Energy: " << u << " kcal/mol\n";
  potential->write();
  molsys.scale(a);
  potential->scale(a);
  u = 0;
  f.zero();
  potential->add_to_energy_and_forces(u,f);
  Out() << "Boundary conditions: " << molsys.boundary_conditions << "\n"
	<< "Energy: " << u << " kcal/mol\n";
  potential->write();
  Out() << "\n" << flush;
}

void MSimProgram::write_torsion_profile(Str fname, int i, int j, int k, int l)
{
  Atom *at = molsys.atom;
  const double kT = K_to_kcal_mol(298.15);
  const double r = at[i].position.distance(at[j].position);
  const double t = Angle(at[i].position,at[j].position,at[k].position);
  const Cartesian save = at[i].position;
  potential->init(&molsys);
  double usum = 0;
  CVec f(molsys.atom.size());
  const int n = 100;
  double *uarr = new double[n];
  int ia;
  for (ia = 0; ia < n; ia++) {
    double phi = -M_PI + ia*(2*M_PI)/n;
    at[i].position = ZLocation(at[j].position,at[k].position,at[l].position,r,t,phi);
    double u = 0;
    potential->add_to_energy_and_forces(u,f);
    uarr[ia] = exp(-u/kT);
    usum += uarr[ia];
  }
  ostream *ss = FileStream(fname);
  for (ia = 0; ia < n; ia++)
    *ss << -M_PI + ia*(2*M_PI)/n << " " << n*uarr[ia]/(usum*2*M_PI) << "\n";
  delete ss;
  delete[] uarr;
  at[i].position = save;
}

void MSimProgram::normal_modes(double displacement)
{
  Out() << "Calculating normal modes...\n" << flush;
  potential->init(&molsys);
  const int n = molsys.atom.size();
  Atom *atom = molsys.atom;
  /*
   * Calculate mass-weighted Hessian matrix by 
   * finite difference of gradient vector
   */
  Out() << "Computing mass-weighted Hessian by finite difference...\n" << flush;
  DMat m(3*n,3*n);
  CVec f(n), f0(n);
  double u = 0;
  int i;
  for (i = 0; i < n; i++) {
    insist(!atom[i].is_dummy());
    int j;
    /* x */
    atom[i].position.x -= displacement;
    f0.zero();
    potential->add_to_energy_and_forces(u,f0);
    atom[i].position.x += 2*displacement;
    f.zero();
    potential->add_to_energy_and_forces(u,f);
    for (j = 0; j < n; j++) {
      Cartesian tmp = (f0[j] - f[j])/(2*displacement*sqrt(atom[i].mass()*atom[j].mass()));
      m(3*i,3*j) = tmp.x;
      m(3*i,3*j+1) = tmp.y;
      m(3*i,3*j+2) = tmp.z;
    }
    atom[i].position.x -= displacement;
    /* y */
    atom[i].position.y -= displacement;
    f0.zero();
    potential->add_to_energy_and_forces(u,f0);
    atom[i].position.y += 2*displacement;
    f.zero();
    potential->add_to_energy_and_forces(u,f);
    for (j = 0; j < n; j++) {
      Cartesian tmp = (f0[j] - f[j])/(2*displacement*sqrt(atom[i].mass()*atom[j].mass()));
      m(3*i+1,3*j) = tmp.x;
      m(3*i+1,3*j+1) = tmp.y;
      m(3*i+1,3*j+2) = tmp.z;
    }
    atom[i].position.y -= displacement;
    /* z */
    atom[i].position.z -= displacement;
    f0.zero();
    potential->add_to_energy_and_forces(u,f0);
    atom[i].position.z += 2*displacement;
    f.zero();
    potential->add_to_energy_and_forces(u,f);
    for (j = 0; j < n; j++) {
      Cartesian tmp = (f0[j] - f[j])/(2*displacement*sqrt(atom[i].mass()*atom[j].mass()));
      m(3*i+2,3*j) = tmp.x;
      m(3*i+2,3*j+1) = tmp.y;
      m(3*i+2,3*j+2) = tmp.z;
    }
    atom[i].position.z -= displacement;
  }
  Out() << "Mass-weighted Hessian matrix:\n";
  for (int i = 0; i < 3*n; i++) {
    for (int j = 0; j < 3*n; j++)
      OutPrintf("%12.8f ", m(i,j));
    Out() << "\n";
  }
  Out() << "Diagonalizing...\n" << flush;
  DVec eval(3*n);
  DMat evec(3*n,3*n);
  Diagonalize(m, eval, evec);
  Out()  << "Eigenvectors: " << evec
	 << "Frequencies:\n";
  const double hbar = hbar_to_action_unit(1.0);
  double zero_point_energy = 0;
  for (i = 0; i < 3*n; i++) {
    double freq;
    char im[5];
    if (eval[i] >= 0) {
      freq = sqrt(eval[i]);
      snprintf(im,5," ");
      zero_point_energy += 0.5*hbar*freq;
    } else {
      freq = sqrt(-eval[i]);
      snprintf(im,5," i ");
    }
    Out() << freq << (const char *) im << "psec-1 ("
	  << angular_velocity_unit_to_cm1(freq) << " cm-1)\n"
	  << flush;
  }
  Out() << "Zero-point energy: " << zero_point_energy << " kcal/mol\n\n" << flush;
}

MSimPotential::MSimPotential() : 
  eleclj_(0), intra_(0), contact_(0),
  covalent_(0), water_intra_(0), harmonic_(0), pair_potential_(0)
{ }

MSimPotential::MSimPotential(const MSimPotential &m) :
  eleclj_(0), intra_(0), contact_(0),
  covalent_(0), water_intra_(0), harmonic_(0), pair_potential_(0)
{
  if (m.eleclj_)
    pp.add(eleclj_ = new ElecLJ(*m.eleclj_));
  if (m.intra_)
    pp.add(intra_ = new IntramolecularPotential(*m.intra_));
  if (m.contact_)
    pp.add(contact_ = new Contact(*m.contact_));
  if (m.covalent_)
    pp.add(covalent_ = new Covalent(*m.covalent_));
  if (m.water_intra_)
    pp.add(water_intra_ = new WaterIntra(*m.water_intra_));
  if (m.harmonic_)
    pp.add(harmonic_ = new Harmonic(*m.harmonic_));
  if (m.pair_potential_)
    pp.add(pair_potential_ = new PairPotential(*m.pair_potential_));
}

MSimPotential::~MSimPotential()
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    delete p();
}

void MSimPotential::init(MSys *m)
{
  Potential::init(m);
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->init(m);
}

void MSimPotential::add_to_energy_and_forces(double &u, Cartesian *f)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->add_to_energy_and_forces(u,f);
}

void MSimPotential::add_to_hessian(double **h)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->add_to_hessian(h);
}

void MSimPotential::write() const
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->write();
}

void MSimPotential::is_running_dynamics(bool b)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->is_running_dynamics(b);
}

void MSimPotential::temperature(double t)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->temperature(t);
}

bool MSimPotential::is_translationally_invariant() const
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    if (!p()->is_translationally_invariant())
      return false;
  return true;
}

bool MSimPotential::is_rotationally_invariant() const
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    if (!p()->is_rotationally_invariant())
      return false;
  return true;
}

void MSimPotential::integrate_positions(double dt)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->integrate_positions(dt);
}

void MSimPotential::integrate_velocities(double dt)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->integrate_velocities(dt);
}

void MSimPotential::scale(double s)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->scale(s);
}

void MSimPotential::types_changed()
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->types_changed();
}

void MSimPotential::types_changed(const MSys &m)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->types_changed(m);
}

void MSimPotential::set_lambda(double l)
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    p()->set_lambda(l);
}

bool MSimPotential::can_calculate_hessian() const
{
  LstItr<Potential *> p;
  for (p.init(pp); p.ok(); p.next())
    if (!p()->can_calculate_hessian())
      return false;
  return true;
}

void MSimPotential::eleclj(istream &s)
{
  if (!eleclj_)
    pp.add(eleclj_ = new ElecLJ);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *eleclj_;
}

void MSimPotential::intra(istream &s)
{
  if (!intra_)
    pp.add(intra_ = new IntramolecularPotential);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *intra_;
}

void MSimPotential::contact(istream &s)
{
  if (!contact_)
    pp.add(contact_ = new Contact);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *contact_;
}

void MSimPotential::covalent(istream &s)
{
  if (!covalent_)
    pp.add(covalent_ = new Covalent);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *covalent_;
}

void MSimPotential::harmonic(istream &s)
{
  if (!harmonic_)
    pp.add(harmonic_ = new Harmonic);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *harmonic_;
}

void MSimPotential::pair_potential(istream &s)
{
  if (!pair_potential_)
    pp.add(pair_potential_ = new PairPotential);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *pair_potential_;
}

void MSimPotential::water_intra(istream &s)
{
  if (!water_intra_)
    pp.add(water_intra_ = new WaterIntra);
  char bracket = 0;
  s >> bracket;
  s.putback(bracket);
  if (bracket == '{')
    s >> *water_intra_;
}

void MSimProgram::write_electrostatic_potential_at_gridpoints(CVec gp)
{
  if (!potential->eleclj_)
    die("MSimPotential::write_electrostatic_potential_at_gridpoints: no eleclj in potential");    
  potential->init(&molsys);
  double u = 0;
  CVec f(molsys.atom.size());
  f.zero();
  potential->add_to_energy_and_forces(u,f);
  potential->eleclj_->write_electrostatic_potential_at_gridpoints(gp);
}

void MSimProgram::fit_charges_to_electrostatic_potential(istream &s)
{
  ESPFit esp(molsys);
  s >> esp;
  esp.fit();
}

void MSimProgram::fit_bond_charge_increments(istream &s)
{
  BCIFit b;
  s >> b;
}

void MSimProgram::extract_dimers(double r, Str dir)
{
  const BoundaryConditions &bc = *molsys.boundary_conditions;
  char cmd[1024];
  snprintf(cmd, 1024, "rm -rf %s ; mkdir %s", 
	  (const char *) dir,
	  (const char *) dir);
  ::system(cmd);
  int ndimer = 0;
  for (int i = 0; i < molsys.molecule.size(); i++) {
    const Cartesian cmi = molsys.molecule_center_of_mass(i);
    for (int j = i+1; j < molsys.molecule.size(); j++) {
      const Cartesian cmj = molsys.molecule_center_of_mass(j);
      const double rij = bc.minimum_image_distance(cmi,cmj);
      if (rij > r)
	continue;
      Lst<int> ij;
      ij.add(i);
      ij.add(j);
      MSys dimer;
      dimer.create_from_molecules(molsys,ij);
      dimer.translate(-cmi);
      dimer.map_to_central_box();
      delete dimer.boundary_conditions;
      dimer.boundary_conditions = BoundaryConditions::new_non_periodic();
      potential->init(&dimer);
      double u = 0;
      CVec f(dimer.atom.size());
      potential->add_to_energy_and_forces(u,f);
      char fname[1024];
      snprintf(fname, 1024, "%s/%s%04d.xyz", 
	      (const char *) dir, 
	      (const char *) dir,
	      ndimer++);
      char cmt[1024];
      snprintf(cmt, 1024, "Distance: %f A  Energy: %f kcal/mol", rij, u);
      dimer.write_as_xyz(Str(fname), Str(cmt));
    }
  }
}

WeakCouple::WeakCouple() 
  : p(0), d(0),
    target_density(UNDEF_VAL), 
    target_potential_energy_per_molecule(UNDEF_VAL),
    target_optical_dielectric_constant(UNDEF_VAL),
    sigma_adjustment_factor(1e-2),
    epsilon_adjustment_factor(1e-3),
    bci_adjustment_factor(1e-3)
{ }

WeakCouple::~WeakCouple()
{
  p = 0;
  d = 0;
}

void WeakCouple::wc_init(ElecLJ *p_, MDyn *d_)
{
  if (are_approximately_equal(target_density,UNDEF_VAL))
    die("WeakCouple::init: need to specify target_density");
  if (are_approximately_equal(target_potential_energy_per_molecule,UNDEF_VAL))
    die("WeakCouple::init: need to specify target_potential_energy_per_molecule");
  insist(target_density > 0);
  insist(sigma_adjustment_factor >= 0);
  insist(epsilon_adjustment_factor >= 0);
  p = p_;
  d = d_;
  if (!p)
    die("WeakCouple::init: no ElecLJ potential specified");
  insist(d);
  if (d->msys.boundary_conditions->type == non_periodic)
    die("WeakCouple::init: must have periodic boundary conditions");
  LstItr<Str> sig, eps, bci;
  for (sig.init(sigma); sig.ok(); sig.next())
    if (!p->sigma.exists(sig()))
      die("WeakCouple::init_potential: no sigma for '%s'", (const char *) sig());
  for (eps.init(epsilon); eps.ok(); eps.next())
    if (!p->epsilon.exists(eps()))
      die("WeakCouple::init_potential: no epsilon for '%s'", (const char *) eps());
  for (bci.init(bond_charge_increment); bci.ok(); bci.next())
    if (!p->bond_charge_increment.exists(bci()))
      die("WeakCouple::init_potential: no bond charge increment for '%s'", (const char *) bci());
}

void WeakCouple::update(const MSys &msys)
{
  LstItr<Str> sig, eps, bci;
  if (sigma.size() > 0) {
    const double density = density_unit_to_g_cm3(msys.density());
    for (sig.init(sigma); sig.ok(); sig.next())
      p->sigma[sig()] *= 1 + sigma_adjustment_factor * (density - target_density)/target_density;
  }
  if (epsilon.size() > 0) {
    const double u = d->u / msys.molecule.size();
    for (eps.init(epsilon); eps.ok(); eps.next())
      p->epsilon[eps()] *= 1 - epsilon_adjustment_factor * (u - target_potential_energy_per_molecule)/
	target_potential_energy_per_molecule;
  }
  if (bond_charge_increment.size() > 0) {
    const double e = p->optical_dielectric_constant();
    for (bci.init(bond_charge_increment); bci.ok(); bci.next())
      p->bond_charge_increment[bci()].k *= 1 + bci_adjustment_factor * 
	(e - target_optical_dielectric_constant)/target_optical_dielectric_constant;
  }
  p->types_changed();
}
  
void WeakCouple::write() const
{
  LstItr<Str> sig, eps, bci;
  for (sig.init(sigma); sig.ok(); sig.next())
    Out() << "Sigma for " << sig() << ": " << p->sigma[sig()] << "\n";
  for (eps.init(epsilon); eps.ok(); eps.next())
    Out() << "Epsilon for " << eps() << ": " << p->epsilon[eps()] << "\n";
  for (bci.init(bond_charge_increment); bci.ok(); bci.next())
    Out() << "Bond charge increment for " << bci() << ": " << p->bond_charge_increment[bci()].k << "\n";
}
