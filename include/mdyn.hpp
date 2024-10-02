#ifndef MDYN_H
#define MDYN_H

#include "cvec.hpp"
#include "str.hpp"
#include "lst.hpp"
#include "boo.hpp"
#include "msys.hpp"
#include "range.hpp"
#include "callback.hpp"
#include "zmatrix.hpp"

class Potential;
class Constraint;

class MDyn
{

  friend class WeakCouple;
  
public:
  
  MDyn(Potential &, Constraint *, CallbackList &);

  int steps, write_interval, verbose; // in
  double timestep, kT, pressure, pH; // in
  Boolean constant_temperature, constant_pressure, constant_pH; // in
  Boolean replica_exchange, path_integrals; // in
  double hbar; // in
  
  int reinit_velocities_interval; // in
  int volume_move_interval; // in
  double volume_move_factor; // in
  int protonation_state_move_interval; // in
  int protonation_state_move_steps; // in
  int protonation_state_move_interpolation_order; // in
  int replica_exchange_interval; // in
  
  void temperature(double T); // in
  Str restart_file, xyz_trajectory_file; // in
  Boolean show_as_xyz, show_geometry, show_forces; // in
  ZMatrix show_as_zmatrix; // in 
 
  void init_velocities(); // in
  void init_velocities(const Vec<int>);
  void run(); // in
  void read_frame(FILE *f, double &time_elapsed);
  void read_xyz_trajectory(Str traj); // in

  void freeze_atoms(LstRange r); // in
  void freeze_everything_except_for_types(LstStr l); // in

  double time_elapsed() const;
  int degrees_of_freedom() const;

private:

  MSys &msys;
  Potential &potential;
  Lst<CallbackPtr> &callback;
  Constraint *constraint;
  int istep, mpi_size, mpi_rank;
  double u, upath;
  CVec f, rtmp, vtmp, rnext, rprev;
  Vec<bool> is_frozen;
  double energy_offset;
  int nreptry, nrepaccept, npmovetry, npmoveaccept;
  MSys tmpmsys;
  
  classIO(MDyn);
  void write();
  void callback_init(), callback_update(), callback_write();
  void make_volume_move(), make_protonation_state_move();
  int how_many_frozen() const;
  void zero_velocities_of_frozen_atoms();
  double probability_of_randomly_placing_atom(int i);
  void randomly_place_atom(int i);
  void try_replica_exchange();
  void path_integral_energy_and_forces(Cartesian *);
  double average_bead_energy();
  
};

#endif /* MDYN_H */
