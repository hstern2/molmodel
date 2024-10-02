#include "checkgr.hpp"
#include "random.h"
#include "cvec.hpp"
#include "out.hpp"
#include "pot.hpp"
#include "msys.hpp"
#include "fns.h"
#include "constraint.hpp"
#include "array.h"

static void randomize(Cartesian &c, double displacement)
{
  c.x = displacement*GaussianRandom();
  c.y = displacement*GaussianRandom();
  c.z = displacement*GaussianRandom();
}

static void randomize(CVec dx, double displacement)
{
  for (int i = 0; i < dx.size(); i++)
    randomize(dx[i], displacement);
}

static void multiply(int n, double *a, const double **m, const double *b)
{
  for (int i = 0; i < n; i++) {
    a[i] = 0;
    for (int j = 0; j < n; j++)
      a[i] += m[i][j] * b[j];
  }
}

void CheckGradients::run()
{
  const int n = msys.atom.size();
  Out() << "Checking gradients with finite difference...\n"
	<< "Displacement: " << displacement << "\n";
  CVec x(n), x0(n), dx(n), f(n);
  msys.copy_positions_to(x0);
  double u0 = 0;
  f.zero();
  p.add_to_energy_and_forces(u0,f);
  if (c)
    c->add_to_energy_and_forces(u0,f);
  OutPrintf("\n%8s%17s%17s%17s\n", "Trial", "U(x+dx) - U(x)", "F dx", "% error");
  Out() << "-----------------------------------------------------------\n";
  for (int itrial = 1; itrial <= trials; itrial++) {
    randomize(dx, displacement);
    x.copy(x0);
    x += dx;
    msys.copy_positions_from(x);
    double u = 0;
    f.zero();
    p.add_to_energy_and_forces(u,f);
    if (c)
      c->add_to_energy_and_forces(u,f);
    const double du = u - u0;
    const double fdx = f*dx;
    OutPrintf("%8d%17.8g%17.8g%17.8f\n", itrial, du, -fdx, 100.0*fabs((du+fdx)/(fdx+machine_epsilon())));
  }
  if (p.can_calculate_hessian()) {
    msys.copy_positions_from(x0);
    double **h = double_array2D_new(3*n,3*n);
    memset(h[0], 0, 9*n*n*sizeof(double));
    p.add_to_hessian(h);
    if (verbose) {
      Out() << "Hessian:\n";
      for (int i = 0; i < 3*n; i++) {
	for (int j = 0; j < 3*n; j++)
	  OutPrintf("%12.8f ", h[i][j]);
	Out() << "\n";
      }
    }
    Out() << "\nChecking Hessian:\n";
    OutPrintf("\n%8s%17s%17s%17s\n", "Trial", "dU", "dx1 H dx2", "% error");
    Out() << "-----------------------------------------------------------\n";
    CVec dy(n);
    CVec Hdx(n);
    for (int itrial = 1; itrial <= trials; itrial++) {
      randomize(dx, 100*displacement);
      randomize(dy, 100*displacement);
      int i;
      for (i = 0; i < n; i++)
	x[i] = x0[i] + dx[i];
      msys.copy_positions_from(x);
      double udx = 0;
      p.add_to_energy_and_forces(udx,f);
      x.copy(x0);
      for (i = 0; i < n; i++)
	x[i] = x0[i] + dy[i];
      msys.copy_positions_from(x);
      double udy = 0;
      p.add_to_energy_and_forces(udy,f);
      for (i = 0; i < n; i++)
	x[i] = x0[i] + dx[i] + dy[i];
      msys.copy_positions_from(x);
      double udxdy = 0;
      p.add_to_energy_and_forces(udxdy,f);
      multiply(3*n, Hdx, (const double **) h, dx);
      const double dxHdy = Hdx*dy;
      const double dU = udxdy + u0 - udx - udy;
      OutPrintf("%8d%17.8g%17.8g%17.8f\n", itrial, dU, dxHdy, 100.0*fabs((dU-dxHdy)/(dxHdy+small_val())));
    }
    double_array2D_delete(h);
  }
  msys.copy_positions_from(x0);
  Out() << "\n" << flush;
}
