#include "units.h"
#include "mdyn.hpp"
#include "pot.hpp"
#include "constraint.hpp"
#include "random.h"
#include "callback.hpp"
#include "timing.h"
#include "cgmin.h"
#include "geograd.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

MDyn::MDyn(Potential &p, Constraint *c, CallbackList &a) :
  steps(10), write_interval(1), verbose(0),
  timestep(0.001), 
  kT(K_to_kcal_mol(298.15)),
  pressure(1.01325),
  pH(7.0),
  constant_temperature(false),
  constant_pressure(false),
  constant_pH(false),
  replica_exchange(false),
  path_integrals(false),
  hbar(1.0),
  reinit_velocities_interval(500),
  volume_move_interval(50),
  volume_move_factor(1e-3),
  protonation_state_move_interval(5000),
  protonation_state_move_steps(5000),
  protonation_state_move_interpolation_order(1),
  replica_exchange_interval(50),
  msys(*p.msys), potential(p), 
  callback(a), constraint(c), istep(0), u(0), f(msys.atom.size()),
  is_frozen(msys.atom.size(), false), energy_offset(0), 
  nreptry(0), nrepaccept(0), npmovetry(0), npmoveaccept(0)
{ }

int MDyn::how_many_frozen() const
{
  int n = 0;
  for (int i = 0; i < is_frozen.size(); i++)
    if (is_frozen[i])
      n++;
  return n;
}

static bool do_every(int interval)
{
  return interval > 0 && interval * UniformRandom() < 1;
}

void MDyn::run()
{
  Out() << "Starting dynamics...\n"
	<< "Steps: " << steps << "\n"
	<< "Timestep: " << timestep << " psec\n"
	<< "Constraints: " << (constraint ? constraint->how_many() : 0) << "\n"
	<< "Frozen atoms: " << how_many_frozen() << "\n"
	<< "Degrees of freedom: " << degrees_of_freedom() << "\n"
	<< "Equilibration temperature: " << kcal_mol_to_K(kT) << " K\n"
	<< "Equilibration pressure: " << pressure << " bar\n"
	<< "Equilibration kT: " << kT << " kcal/mol\n";
  if (constant_temperature && reinit_velocities_interval > 0)
    Out() << "Constant-temperature simulation:\n"
	  << "Reinitalizing velocites from Boltzmann distribution every "
	  << reinit_velocities_interval << " steps on average.\n"; 
  if (constant_pressure && volume_move_interval > 0) {
    if (msys.boundary_conditions->type == non_periodic)
      die("MDyn::run: boundary conditions must be periodic "
	  "for stochastic volume moves\n");
    Out() << "Constant-pressure simulation:\n"
	  << "Making volume moves every " 
	  << volume_move_interval << " steps on average.\n"
	  << "Volume move scale factor: " 
	  << volume_move_factor << "\n";
  }
  if (constant_pH && protonation_state_move_interval > 0 && 
      msys.acid.size() > 0)
    Out() << "Constant-pH simulation at pH " << pH << "\n"
	  << "Making protonation state moves every "
	  << protonation_state_move_interval << " steps on average.\n"
	  << "Each protonation state move is "
	  << protonation_state_move_steps << " steps.\n"
	  << "Order of polynomial for interpolation: "
	  << protonation_state_move_interpolation_order << "\n";
  if (replica_exchange && constant_pH)
    die("MDyn::run: cannot run replica exchange and constant pH");
  if (path_integrals && constant_pH)
    die("MDyn::run: cannot run path-integral MD and constant pH");
  if (replica_exchange && path_integrals)
    die("MDyn::run: cannot run replica exchange and path integral MD");
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
  if (replica_exchange) {
#ifdef USE_MPI
    Out() << "Running replica exchange with " << mpi_size << " replicas.\n"
	  << "This is replica " << mpi_rank << ".\n";
#else
    die("MDyn::run: cannot do replica exchange "
	"because this executable was not compiled with MPI");
#endif 
  }
  if (path_integrals) {
#ifdef USE_MPI
    Out() << "Running path-integral MD with " << mpi_size << " beads.\n"
	  << "This is bead " << mpi_rank << ".\n";
#else
    die("MDyn::run: cannot run path-integral MD "
	"because this executable was not compiled with MPI");
#endif 
  }
  if (show_as_zmatrix.size() > 0)
    show_as_zmatrix.check(msys.atom.size());
  potential.is_running_dynamics(true);
  Out() << "\n" << flush;
  if (constraint)
    constraint->constrain_positions(0);
  u = 0;
  f.zero();
  potential.add_to_energy_and_forces(u,f);
  if (path_integrals)
    path_integral_energy_and_forces(f);
  istep = 0;
  energy_offset = 0;
  nreptry = nrepaccept = npmovetry = npmoveaccept = 0;
  LstItrMod<CallbackPtr> c;
  for (c.init(callback); c.ok(); c.next()) {
    c()->f = f;
    c()->init(msys,write_interval);
    c()->maybe_update(0, 0, msys);
  }
  write();
  ostream *xyz_traj = 0;
  if (strlen(xyz_trajectory_file) > 0) {
    xyz_traj = FileStream(xyz_trajectory_file);
    msys.write_as_xyz(*xyz_traj,0.0);
  }
  while (istep < steps) {
    TIMESTART("Timestep");
    msys.integrate_velocities(f, timestep/2);
    zero_velocities_of_frozen_atoms();
    potential.integrate_velocities(timestep/2);
    msys.integrate_positions(timestep);
    potential.integrate_positions(timestep);
    if (constraint)
      constraint->constrain_positions(timestep);
    u = 0;
    f.zero();
    potential.add_to_energy_and_forces(u,f);
    if (path_integrals)
      path_integral_energy_and_forces(f);
    msys.integrate_velocities(f, timestep/2);
    zero_velocities_of_frozen_atoms();
    potential.integrate_velocities(timestep/2);
    if (constraint)
      constraint->constrain_velocities();
    if (constant_temperature && do_every(reinit_velocities_interval)) {
      energy_offset += msys.kinetic_energy();
      init_velocities();
      energy_offset -= msys.kinetic_energy();
    }
    if (constant_pressure) {
      bool should = do_every(volume_move_interval);
#ifdef USE_MPI
      if (path_integrals)
	MPI_Bcast(&should, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if (should) {
	energy_offset += u;
	make_volume_move();
	energy_offset -= u;
      }
    }
    if (constant_pH && do_every(protonation_state_move_interval)) {
      energy_offset += u + msys.kinetic_energy();
      make_protonation_state_move();
      energy_offset -= u + msys.kinetic_energy();
    }
    if (replica_exchange && do_every(replica_exchange_interval))
      try_replica_exchange();
    istep++;
    LstItrMod<CallbackPtr> c;
    for (c.init(callback); c.ok(); c.next())
      c()->maybe_update(time_elapsed(), istep, msys);
    if (is_not_a_number(u)) {
      write();
      die("Dynamics: NaN detected.");
    }
    if (fopen("stopdyn","r")) {
      write();
      die("Dynamics: 'stopdyn' detected.");
    }
    if (write_interval > 0 && istep % write_interval == 0) {
      write();
      if (xyz_traj)
	msys.write_as_xyz(*xyz_traj, time_elapsed());
    }
    TIMESTOP("Timestep");
  }
  if (xyz_traj)
    delete xyz_traj;
  potential.is_running_dynamics(false);
}

void MDyn::read_frame(FILE *f, double &time_elapsed)
{
  time_elapsed = 0;
  fpos_t pos;
  fgetpos(f, &pos);
  char buf[1024];
  int i;
  for (i = 0; i < 2; i++)
    if (!fgets(buf, 1024, f))
      break;
  if (i == 2) {
    char *tmp = strstr(buf, "Time elapsed:");
    if (tmp)
      sscanf(tmp, "Time elapsed: %lf", &time_elapsed);
  }
  fsetpos(f, &pos);
  msys.read_as_xyz(f);
}

void MDyn::read_xyz_trajectory(Str traj)
{
  FILE *ftraj = FileSearch(traj);
  Out() << "Reading trajectory " << traj << "\n"
	<< "Starting dynamics...\n"
	<< "Constraints: " << (constraint ? constraint->how_many() : 0) << "\n"
	<< "Frozen atoms: " << how_many_frozen() << "\n"
	<< "Degrees of freedom: " << degrees_of_freedom() << "\n"
	<< "\n" << flush;
  const double timestep_save = timestep;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
  if (path_integrals) {
#ifdef USE_MPI
    Out() << "Running path-integral MD with " << mpi_size << " beads.\n"
	  << "This is bead " << mpi_rank << ".\n";
#else
    die("MDyn::run: cannot run path-integral MD "
	"because this executable was not compiled with MPI");
#endif 
  }
  for (istep = 0; !end_of_file(ftraj); istep++) {
    LstItrMod<CallbackPtr> c;
    double t;
    read_frame(ftraj, t);
    if (istep == 0) {
      potential.init(&msys);
      if (constraint)
	constraint->init(&msys);
      for (c.init(callback); c.ok(); c.next()) {
	c()->update_interval = 1;
	c()->init(msys,write_interval);
      }
      timestep = 0;
    } else {
      timestep = t/istep;
    }
    u = 0;
    f.zero();
    potential.add_to_energy_and_forces(u,f);
    if (path_integrals)
      path_integral_energy_and_forces(f);
    for (c.init(callback); c.ok(); c.next())
      c()->maybe_update(t, istep, msys);
    write();
  }
  timestep = timestep_save;
  fclose(ftraj);
}

#ifdef USE_MPI

static void send_int(int proc, int buf)
{
  MPI_Send(&buf, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
}

static int receive_int(int proc)
{
  int buf;
  MPI_Status s;
  MPI_Recv(&buf, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &s);
  return buf;
}

static void send_double(int proc, double buf)
{
  MPI_Send(&buf, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
}

static double receive_double(int proc)
{
  double buf;
  MPI_Status s;
  MPI_Recv(&buf, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &s);
  return buf;
}

static void exchange(int proc, CVec r)
{
  MPI_Status s;
  MPI_Sendrecv_replace((Cartesian *) r, 3*r.size(), MPI_DOUBLE, proc, 0, proc, 0, 
		       MPI_COMM_WORLD, &s);
}

#endif

void MDyn::try_replica_exchange()
{
#ifdef USE_MPI
  if (mpi_rank < mpi_size-1) {
    const double kT0 = kT;
    const double u0 = u;
    const double kT1 = receive_double(mpi_rank+1);
    const double u1 = receive_double(mpi_rank+1);
    const double du = -u1/kT0 - u0/kT1 + u0/kT0 + u1/kT1;
    if (du < 0 || UniformRandom() < exp(-du)) {
      send_int(mpi_rank+1, 1);
      rtmp.resize(msys.atom.size());
      msys.copy_positions_to(rtmp);
      exchange(mpi_rank+1, rtmp);
      msys.copy_positions_from(rtmp);
      nrepaccept++;
    } else {
      send_int(mpi_rank+1, 0);
    }
  }
  if (mpi_rank > 0) {
    send_double(mpi_rank-1, kT);
    send_double(mpi_rank-1, u);
    if (receive_int(mpi_rank-1)) {
      rtmp.resize(msys.atom.size());
      msys.copy_positions_to(rtmp);
      exchange(mpi_rank-1, rtmp);
      msys.copy_positions_from(rtmp);
    }
  }
  nreptry++;
#else
  die("MDyn::try_replica_exchange: this executable was not compiled with MPI");
#endif
}
	     
void MDyn::write()
{
  const double kin = msys.kinetic_energy();
  const int ndof = degrees_of_freedom();
  const double tmpkT = ndof > 0 ? 2.0*kin/ndof : 0;
  double conserved = kin + u + energy_offset;
#ifdef USE_MPI
  if (path_integrals) {
    conserved += upath;
    double etot;
    MPI_Reduce(&conserved, &etot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&etot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    conserved = etot;
  }
#endif
  Out() << "Timestep: " << istep << "\n"
	<< "Time elapsed: " << time_elapsed() << " psec\n"
	<< "Total energy: " << kin + u << " kcal/mol\n"
	<< "Potential energy: " << u << " kcal/mol\n"
	<< "Potential energy per molecule: " << u/msys.molecule.size() << " kcal/mol\n"
	<< "Kinetic energy: " << kin << " kcal/mol\n"
	<< "Conserved energy: " << conserved << " kcal/mol\n"
	<< "kT: " << tmpkT << " kcal/mol\n"
	<< "Temperature: " << kcal_mol_to_K(tmpkT) << " K\n"
	<< "Boundary conditions: " << msys.boundary_conditions << "\n";
  if (msys.boundary_conditions->type != non_periodic) {
    const double v = msys.boundary_conditions->volume();
    Out() << "Volume: " << v << " A^3\n"
	  << "Molar volume: " << A3_to_L_mol(v/msys.molecule.size()) << " L/mol\n"
	  << "Density: " << density_unit_to_g_cm3(msys.density()) << " g/cm^3\n"
	  << "Number density: " << msys.number_density() << "/A^3\n";
  }
  Out() << "Center of mass (A): " << msys.center_of_mass() << "\n"
	<< "Linear momentum (kcal/mol psec/A): " << msys.linear_momentum() << "\n"
	<< "Angular momentum (kcal/mol psec): " << msys.angular_momentum() << "\n";
  if (path_integrals)
    Out() << "Path-integral energy: " << upath << " kcal/mol\n";
  if (replica_exchange)
    Out() << "Replica exchanges attempted between this process and next higher: "
	  << nreptry << "\n"
	  << "Fraction of replica exchanges accepted between this process and next higher: "
	  << (nreptry > 0 ? double(nrepaccept)/double(nreptry) : 0.0)
	  << "\n";
  if (constant_pH && npmovetry > 0)
    Out() << "Number of protonation state moves attempted: "
	  << npmovetry << "\n"
	  << "Number of protonation state moves accepted: "
	  << npmoveaccept << "\n"
	  << "Fraction of protonation state moves accepted: "
	  << double(npmoveaccept)/double(npmovetry) << "\n"
	  << "Fraction of protonation state moves accepted per step: "
	  << double(npmoveaccept)/double(npmovetry*protonation_state_move_steps) << "\n";
  if (show_geometry)
    msys.show_geometry();
  if (show_as_xyz)
    msys.show_as_xyz();
  if (show_as_zmatrix.size() > 0) {
    Out() << "ZMatrix: ";
    msys.show_as_zmatrix(show_as_zmatrix);
  }
  if (show_forces)
    Out() << "Forces: " << f << "\n";
  potential.write();
  LstItr<CallbackPtr> c;
  for (c.init(callback); c.ok(); c.next())
    c()->maybe_write(time_elapsed());
  Out() << "\n\n" << flush;
  if (strlen(restart_file) > 0)
    msys.write_as_xyz(restart_file);
}

double MDyn::time_elapsed() const { return istep*timestep; }

int MDyn::degrees_of_freedom() const
{
  const int nat = msys.atom.size();
  const int nfrozen = how_many_frozen();
  int n = 3*nat;
  if (nfrozen == 0 && potential.is_translationally_invariant())
    n -= 3;
  if (nfrozen == 0 && 
      msys.boundary_conditions->type == non_periodic &&
      potential.is_rotationally_invariant()) {
    if (nat == 2)
      n -= 2;
    else
      n -= 3;
  }
  if (constraint)
    n -= constraint->how_many();
  n -= 3*nfrozen;
  return n > 0 ? n : 0;
}

void MDyn::freeze_atoms(LstRange r)
{
  LstItr<Range> i;
  for (i.init(r); i.ok(); i.next()) {
    if (i().lo < 0 || i().hi >= msys.atom.size())
      die("MDyn::freeze_atoms: indices %d..%d are out of range",
	  i().lo, i().hi);
    for (int j = i().lo; j <= i().hi; j++)
      is_frozen[j] = true;
  }
}

void MDyn::freeze_everything_except_for_types(LstStr l)
{
  LstItr<Str> j;
  HTab<bool> in;
  for (j.init(l); j.ok(); j.next())
    in[j()] = true;
  for (int i = 0; i < msys.atom.size(); i++)
    is_frozen[i] = !in.exists(msys.atom[i].type);
}

void MDyn::zero_velocities_of_frozen_atoms()
{
  for (int i = 0; i < msys.atom.size(); i++)
    if (is_frozen[i])
      msys.atom[i].velocity.zero();
}

void MDyn::init_velocities()
{
  bool any_frozen = false;
  for (int i = 0; i < msys.atom.size(); i++) {
    Atom &at = msys.atom[i];
    if (is_frozen[i]) {
      any_frozen = true;
      at.velocity.zero();
      continue;
    }
    const double a = sqrt(kT/at.mass());
    at.velocity.x = a*GaussianRandom();
    at.velocity.y = a*GaussianRandom();
    at.velocity.z = a*GaussianRandom();
  }
  if (constraint)
    constraint->constrain_velocities();
  if (!any_frozen && potential.is_translationally_invariant())
    msys.rezero_linear_momentum();
  if (!any_frozen && potential.is_rotationally_invariant())
    msys.rezero_angular_momentum();
}

void MDyn::temperature(double T)
{
  kT = K_to_kcal_mol(T);
}

double MDyn::average_bead_energy()
{
#ifdef USE_MPI
  double utot = 0;
  MPI_Reduce(&u, &utot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return utot/mpi_size;
#else
  return 0;
#endif
}

void MDyn::make_volume_move()
{
  if (verbose)
    Out() << "Attempting volume move...\n";
  const double p = bar_to_pressure_unit(pressure);
  const double v0 = msys.boundary_conditions->volume();
  double h0 = 0;
  double a = exp(volume_move_factor * GaussianRandom());
  if (path_integrals) {
#ifdef USE_MPI
    h0 = average_bead_energy() + p*v0;
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    msys.scale_centroids(a);
#else
    die("MDyn::make_volume_move: executable not compiled with MPI");
#endif
  } else {
    h0 = u + p*v0;
    msys.scale(a);
  }
  potential.scale(a);
  const double v = msys.boundary_conditions->volume();
  u = 0;
  f.zero();
  potential.add_to_energy_and_forces(u,f);
  bool should_reject = false;
  if (path_integrals) {
#ifdef USE_MPI
    path_integral_energy_and_forces(f);
    const double dh = (average_bead_energy() + p*v - h0)/kT 
      - msys.molecule.size()*log(v/v0);
    should_reject = dh > 0 && exp(-dh) < UniformRandom();
    MPI_Bcast(&should_reject, 1, MPI_INT, 0, MPI_COMM_WORLD);
#else
    die("MDyn::make_volume_move: executable not compiled with MPI");
#endif
  } else {
    const double dh = (u + p*v - h0)/kT 
      - msys.molecule.size()*log(v/v0);
    should_reject = dh > 0 && exp(-dh) < UniformRandom();
  }
  if (should_reject) {
    msys.scale(1/a);
    potential.scale(1/a);
    u = 0;
    f.zero();
    potential.add_to_energy_and_forces(u,f);
    if (path_integrals)
      path_integral_energy_and_forces(f);
    if (verbose)
      Out() << "Volume move rejected.\n" << flush;
  } else {
    if (verbose)
      Out() << "Volume move accepted.\n" << flush;
  }
}

static double interpolate(double x, int n)
{
  double s = 0;
  for (int k = 0; k <= n-1; k++)
    s += combinations(n-1+k,k)*pow(1-x,k);
  s *= pow(x,n);
  return s;
}

void MDyn::make_protonation_state_move()
{
  TIMESTART("MDyn::protonation_state_move");
  const double log10pH = log(10.0) * pH;
  rtmp.resize(msys.atom.size());
  vtmp.resize(msys.atom.size());
  for (int i = 0; i < msys.acid.size(); i++) {
    tmpmsys = msys;
    AcidicGroup &a = msys.acid[i];
    const int nstate = a.desc->number_of_states();
    if (nstate == 1)
      continue;
    const int sold = a.state;
    const ProtonationState &oldstate = a.desc->state[sold];
    const double e0 = u + msys.kinetic_energy();
    msys.copy_positions_to(rtmp);
    msys.copy_velocities_to(vtmp);
    /* Pick a new protonation state at random
       (distinct from the old, and differing 
       by at most one proton */
    int snew, dN;
    do {
      snew = (int) floor(nstate * UniformRandom());
      dN = a.desc->state[snew].number_of_protons() -
	a.desc->state[sold].number_of_protons();
    } while (snew == sold || dN < -1 || dN > 1);
    if (verbose) {
      Out() << "Attempting protonation state move for group "
	    << i << " from state " << sold << " to state " << snew << ".\n"
	    << "Before attempted move:\n"
	    << "Total energy: " << e0 << " kcal/mol\n";
      potential.write();
      Out() << "\n" << flush;
    }
    tmpmsys.change_protonation_state(i,snew);
    tmpmsys.assign_properties_and_types();
    potential.types_changed(tmpmsys);
    potential.set_lambda(0);
    const bool reverse_velocities = UniformRandom() < 0.5;
    if (reverse_velocities)
      msys.reverse_velocities();
    /* Do MD steps */
    for (int k = 1; k <= protonation_state_move_steps; k++) {
      TIMESTART("Timestep");
      msys.integrate_velocities(f, timestep/2);
      zero_velocities_of_frozen_atoms();
      potential.integrate_velocities(timestep/2);
      msys.integrate_positions(timestep);
      potential.integrate_positions(timestep);
      if (constraint)
	constraint->constrain_positions(timestep);
      potential.set_lambda(interpolate((double) k/(double) protonation_state_move_steps,
				       protonation_state_move_interpolation_order));
      u = 0;
      f.zero();
      potential.add_to_energy_and_forces(u,f);
      msys.integrate_velocities(f, timestep/2);
      zero_velocities_of_frozen_atoms();
      potential.integrate_velocities(timestep/2);
      if (constraint)
	constraint->constrain_velocities();
      TIMESTOP("Timestep");
    }
    if (reverse_velocities)
      msys.reverse_velocities();
    const double e = u + msys.kinetic_energy();
    if (verbose) {
      Out() << "After attempted move:\n"
	    << "Total energy: " << e << " kcal/mol\n";
      potential.write();
      Out() << "\n" << flush;
    }
    npmovetry++;
    /* Determine if move is to be accepted */
    const ProtonationState &newstate = a.desc->state[snew];
    const double x = 
      (e - e0)/kT + log10pH*dN + log((double) newstate.degeneracy) - 
      log((double) oldstate.degeneracy);
    if (x <= 0 || exp(-x) > UniformRandom()) {
      if (verbose)
	Out() << "Protonation state move accepted.\n\n";
      msys.change_protonation_state(i,snew);
      potential.types_changed();
      npmoveaccept++;
    } else {
      if (verbose)
	Out() << "Protonation state move rejected.\n\n";
      tmpmsys.change_protonation_state(i,sold);
      potential.set_lambda(0);
      msys.copy_positions_from(rtmp);
      msys.copy_velocities_from(vtmp);
      if (constraint)
	constraint->reset();
      u = 0;
      f.zero();
      potential.add_to_energy_and_forces(u,f);
    }
  }
  TIMESTOP("MDyn::protonation_state_move");
}

void MDyn::path_integral_energy_and_forces(Cartesian *f)
{
#ifdef USE_MPI
  const double omega2 = mpi_size * sq(kT/hbar_to_action_unit(hbar));
  const int natom = msys.atom.size();
  const int prev = mymod(mpi_rank-1,mpi_size);
  const int next = mymod(mpi_rank+1,mpi_size);
  rnext.resize(natom);
  rprev.resize(natom);
  rtmp.resize(natom);
  msys.copy_positions_to(rtmp);
  MPI_Status status;
  /* Send to previous, and receive from next */
  MPI_Sendrecv(rtmp, 3*natom, MPI_DOUBLE, prev, 0,
	       rnext, 3*natom, MPI_DOUBLE, next, 0,
	       MPI_COMM_WORLD, &status);
  /* Send to next, and receive from previous */
  MPI_Sendrecv(rtmp, 3*natom, MPI_DOUBLE, next, 0,
	       rprev, 3*natom, MPI_DOUBLE, prev, 0,
	       MPI_COMM_WORLD, &status);
  upath = 0;
  for (int i = 0; i < natom; i++) {
    const double k = msys.atom[i].mass()*omega2;
    upath += 0.5*k*(rtmp[i] - rprev[i]).sq();
    f[i] -= k*(2*rtmp[i] - rprev[i] - rnext[i]);
  }
#endif
}
