#ifndef MSIMPROG_H
#define MSIMPROG_H

#include "msys.hpp"
#include "prog.hpp"
#include "pot.hpp"
#include "constraint.hpp"
#include "callback.hpp"
#include "blockstr.hpp"

class ElecLJ;
class IntramolecularPotential;
class Contact;
class Covalent;
class WaterIntra;
class Harmonic;
class PairPotential;
class MDyn;


typedef Lst<Str> LstStr;

class MSimPotential : public Potential
{

public:

  ElecLJ *eleclj_;
  IntramolecularPotential *intra_;
  Contact *contact_;
  Covalent *covalent_;
  WaterIntra *water_intra_;
  Harmonic *harmonic_;
  PairPotential *pair_potential_;
  
  Lst<Potential *> pp;

  MSimPotential();
  MSimPotential(const MSimPotential &);
  MSimPotential & operator=(const MSimPotential &);
  ~MSimPotential();

  void eleclj(istream &); // in
  void intra(istream &); // in
  void covalent(istream &); // in
  void harmonic(istream &); // in
  void contact(istream &); // in
  void water_intra(istream &); // in
  void pair_potential(istream &); // in
  
  void init(MSys *);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  void is_running_dynamics(bool);
  void temperature(double);
  bool is_translationally_invariant() const;
  bool is_rotationally_invariant() const;
  void integrate_positions(double dt);
  void integrate_velocities(double dt);
  void scale(double);
  void types_changed();
  void types_changed(const MSys &);
  void set_lambda(double);
  bool can_calculate_hessian() const;
  void add_to_hessian(double **);
  
  friend istream & operator>>(istream &s, MSimPotential *&p)
  {
    if (p)
      delete p;
    return s >> *(p = new MSimPotential);
  }

  classIO(MSimPotential);
};

class FEP : public Callback // in
{
public:
  BlockStr modify; // in
  double start_lambda, end_lambda; // in
  double ustart, uend, dU;
  CVec f;
  MSys *perturbed_msys;
  MSimPotential *perturbed_potential;
  FEP();
  ~FEP();
  void init_perturbed(MSys *, MSimPotential *, Constraint *);
  friend istream & operator>>(istream &s, FEP *&p)
  {
    if (p)
      delete p;
    return s >> *(p = new FEP);
  }
  classIO(FEP);
private:
  void write() const;
  void update(const MSys &);
};

class WeakCouple : public Callback // in
{
public:
  ElecLJ *p;
  MDyn *d;
  double target_density, target_potential_energy_per_molecule, target_optical_dielectric_constant; // io
  double sigma_adjustment_factor, epsilon_adjustment_factor, bci_adjustment_factor; // io
  LstStr sigma, epsilon, bond_charge_increment; // io
  void wc_init(ElecLJ *, MDyn *);
  WeakCouple();
  ~WeakCouple();
  friend istream & operator>>(istream &s, WeakCouple *&w)
  {
    if (w)
      delete w;
    return s >> *(w = new WeakCouple);
  }
  classIO(WeakCouple);
private:
  void write() const;
  void update(const MSys &);
};

class MSimProgram : public Program // in
{

public:

  MSys molsys; // in
  MSimPotential* potential; // in
  CallbackList callback; // in
  Constraint* constraint; // in
  FEP* free_energy_perturbation; // in
  WeakCouple* weak_couple; // in
  void random_seed(int n); // in
  void energy(); // in
  void forces(); // in
  void finite_difference_forces(); // in
  void test_scale(double a); // in
  void dynamics(istream &); // in
  void check_gradients(istream &); // in
  void minimize(istream &); // in
  void rms_force_error(istream &); // in
  void rms_electric_field_error(istream &); // in
  void write_torsion_profile(Str fname, int i, int j, int k, int l); // in
  void normal_modes(double displacement); // in
  void binding_energy(); // in
  void energies_from_xyz_trajectory(Str fname); // in
  void zmatrices_from_xyz_trajectory(ZMatrix zmat, Str fname); // in
  void write_electrostatic_potential_at_gridpoints(CVec gp); // in
  void fit_charges_to_electrostatic_potential(istream &); // in
  void fit_bond_charge_increments(istream &); // in
  void extract_dimers(double r, Str dir); // in
  
  MSimProgram() : 
    potential(new MSimPotential), 
    constraint(0),
    free_energy_perturbation(0),
    weak_couple(0)
  { }
  ~MSimProgram()
  { 
    if (potential)
      delete potential; 
    if (constraint)
      delete constraint;
    if (free_energy_perturbation)
      delete free_energy_perturbation;
    if (weak_couple)
      delete weak_couple;
    potential = 0;
    constraint = 0;
    free_energy_perturbation = 0;
    weak_couple = 0;
    LstItr<CallbackPtr> p;
    for (p.init(callback); p.ok(); p.next())
      delete p().p;
  }
  classIO(MSimProgram);
};

#endif /* MSIMPROG_H */
