#ifndef POT_H
#define POT_H

#include "msys.hpp"

class Potential
{
public:
  Potential() : msys(0) { }
  virtual ~Potential() { }
  virtual void init(MSys *m) { msys = m; }
  virtual void add_to_energy_and_forces(double &, Cartesian *) { }
  virtual void write() const { }
  virtual bool is_translationally_invariant() const { return true; }
  virtual bool is_rotationally_invariant() const
  { return msys->boundary_conditions->type == non_periodic; }
  virtual void integrate_positions(double) { }
  virtual void integrate_velocities(double) { }
  virtual void scale(double) { }
  virtual void types_changed() { }
  virtual void types_changed(const MSys &perturbed) { }
  virtual void set_lambda(double lambda) { }
  virtual void is_running_dynamics(bool) { }
  virtual void temperature(double) { }
  virtual bool can_calculate_hessian() const { return false; }
  virtual void add_to_hessian(double **) { }
  MSys *msys;
};

#endif /* POT_H */
