/* $Id: cheby.h,v 1.4 2005/08/25 23:34:14 hstern Exp $ */

#ifndef CHEBY_H
#define CHEBY_H

/***
 * Approximate a function by a lookup table
 * of third-order Chebyshev polynomials
 * (see Numerical Recipes section 5.8)
 ***/

#define NCHEBY 512 /* size of lookup table */

#ifdef __cplusplus
extern "C" {
#endif

struct Chebyshev { double min, max, inv_dx, coeff[NCHEBY][7]; };

/* Approximate a function from min to max */
void Chebyshev_init(struct Chebyshev *t, double min, double max, double (*f)(double));

void Chebyshev_test(const struct Chebyshev *t, double (*f)(double), 
		      double *maxerr, double *xerr);

/* x must be in range [min,max] */
static void Chebyshev_evaluate(const struct Chebyshev *t, double x, 
			       double *f, double *n2df)
{
  unsigned i;
  const double *c;
  x -= t->min;
  x *= t->inv_dx;
  i = (unsigned) x;
  c = t->coeff[i];
  x -= i;
  *f = ((c[3]*x + c[2])*x + c[1])*x + c[0];
  *n2df = (c[6]*x + c[5])*x + c[4];
}

#ifdef __cplusplus
}
#endif

#endif /* CHEBY_H */
