#include <cstring>
#include "mmin.hpp"
#include "msys.hpp"
#include "pot.hpp"
#include "cgmin.h"
#include "constraint.hpp"

MMin::MMin(MSys &m, Potential &p_, Constraint *c_) : 
  maximum_evaluations(100000), 
  verbose(1), 
  tolerance(1e-3),
  linesearch_tolerance(0.5),
  msys(m), 
  p(p_),
  c(c_)
{ }

static double calc_fr(int n, const double *x, double *r, void *user)
{
  MMin *t = (MMin *) user;
  t->msys.copy_positions_from((const Cartesian *) x);
  double u = 0;
  memset(r, 0, n*sizeof(double));
  t->p.add_to_energy_and_forces(u, (Cartesian *) r);
  if (t->c)
    t->c->add_to_energy_and_forces(u, (Cartesian *) r);
  return u;
}

void MMin::run()
{
  const int n = 3*msys.atom.size();
  double *x = new double[n];
  double *r = new double[n];
  double *work = new double[4*n];
  msys.copy_positions_to((Cartesian *) x);
  conjugate_gradient_minimize(n, x, r, 0, tolerance, linesearch_tolerance,
			      maximum_evaluations, verbose, calc_fr, 0, this, work);
  p.write();
  delete[] x;
  delete[] r;
  delete[] work;
}
