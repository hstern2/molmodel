#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "cheby.h"

static void chebft(double a, double b, double coeff[], double (*f)(double))
{
  int k,j;
  double g[50], c[4];
  const double fac = 2.0/50;
  const double bma=0.5*(b-a);
  const double bpa=0.5*(b+a);

  for (k = 0; k < 50; k++)
    g[k]=f(cos(M_PI*(k+0.5)/50) * bma + bpa);

  for (j = 0; j < 4; j++) {
    double sum = 0;
    for (k = 0; k < 50; k++)
      sum += g[k]*cos(M_PI*j*(k+0.5)/50);
    c[j] = fac*sum;
  }
  /* coefficients for function */
  coeff[0] = c[0]/2 - c[1] + c[2] - c[3];
  coeff[1] = 2*c[1] - 8*c[2] + 18*c[3];
  coeff[2] = 8*c[2] - 48*c[3];
  coeff[3] = 32*c[3];
  /* coefficients for -2 * derivative */
  coeff[4] = -2*coeff[1]/(b-a);
  coeff[5] = -2*2*coeff[2]/(b-a);
  coeff[6] = -2*3*coeff[3]/(b-a);
}

void Chebyshev_init(struct Chebyshev *t, double min, double max, double (*f)(double))
{
  int i;
  const double dx = (max-min)/NCHEBY;
  assert(max > min);
  t->min = min;
  t->max = max;
  t->inv_dx = 1.0/dx;
  for (i = 0; i < NCHEBY; i++)
    chebft(i*dx + min, (i+1)*dx + min, t->coeff[i], f);
}

void Chebyshev_test(const struct Chebyshev *t, double (*f)(double),
		    double *maxerr, double *xerr)
{
  int i;
  *maxerr = *xerr = 0;
  for (i = 0; i < NCHEBY-1; i++) {
    const double x = (i+drand48())/t->inv_dx + t->min;
    double a, da;
    Chebyshev_evaluate(t, x, &a, &da);
    a = fabs(a - f(x));
    if (a > *maxerr) {
      *xerr = x;
      *maxerr = a;
    }
  }
}

#ifdef SIMPLE_EXAMPLE

int main()
{
  double x, min = 0, max = 10, mrel = 0, mabs = 0, maxerr, xerr;
  struct Chebyshev c;
  Chebyshev_init(&c, min, max, sin);
  Chebyshev_test(&c, sin, &maxerr, &xerr);
  printf("Max error: %g at x=%g\n", maxerr, xerr);
  for (x = min; x < max; x += (max-min)/103.1415) {
    double f, df, fexact, neg2dfexact, rel, abs;
    fexact = sin(x);
    neg2dfexact = -2*cos(x);
    Chebyshev_evaluate(&c, x, &f, &df);
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f\n", x, fexact, f, neg2dfexact, df);
    abs = fabs(f - fexact);
    rel = f > 0 ? abs/f : abs;
    if (abs > mabs)
      mabs = abs;
    if (rel > mrel)
      mrel = rel;
  }
  printf("# Maximum absolute error: %12.8g\n", mabs);
  printf("# Maximum relative error: %12.8g\n", mrel);
  return 1;
}

#endif
