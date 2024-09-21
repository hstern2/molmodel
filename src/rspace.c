#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "rspace.h"
#include "pbc.h"
#include "pbin.h"

static double sq(double x) { return x*x; }



#define EWCHEB /* Chebyshev approximation for erfc */
#ifdef EWCHEB

#include "cheby.h"

static double eta_;
static double pair_function_for_chebyshev(double r2)
{ 
  const double r = sqrt(r2);
  return erfc(eta_*r)/r;
}

struct Chebyshev *RSpaceChebyshev_new(double rhi, double eta, double *maxerr, double *xerr)
{
  struct Chebyshev *t = malloc(sizeof(struct Chebyshev));
  assert(rhi > 0);
  assert(eta > 0);
  const double r2max = 1.2*sq(rhi);
  const double r2min = r2max > 4 ? 4 : 0.9*r2max;
  eta_ = eta;
  Chebyshev_init(t, r2min, r2max, pair_function_for_chebyshev);
  if (maxerr && xerr)
    Chebyshev_test(t, pair_function_for_chebyshev, maxerr, xerr);
  return t;
}

#else /* EWCHEB */

struct Chebyshev *RSpaceChebyshev_new(double rhi, double eta, double *maxerr, double *xerr)
{
  *maxerr = *xerr = 0;
  assert(rhi > 0);
  assert(eta > 0);
  return 0;
}

#endif /* EWCHEB */

/* Long-range switching function */
static void switching_fn(double r2, double rlo2, double rhi2, double *f, double *_2df)
{
  const double a = 1.0/(rhi2 - rlo2);
  const double z = (r2 - rlo2)*a;
  *f = ((-6*z + 15)*z - 10)*z*z*z + 1;
  *_2df = ((60*z - 120)*z + 60)*z*z*a;
}

/* Short-range Coulombic screening function */
static inline void screening_function(double r2, double rscreen, double *b0, double *b1)
{
  *b0 = -0.5*r2*(*b1 = 1.0/(rscreen*rscreen*rscreen)) + 1.5/rscreen;
}

static void coulomb_non_periodic(int i, int j,
				 double x, double y, double z,
				 double r2, double qi, double qj,
				 double *phi, cart_t *evec,
				 double rscreen, int are_smoothing,
				 double smooth, double *tmpsmooth,
				 int are_scaling, double scale_factor)
{
  double b0, b1, qib1, qjb1;
  if (r2 < rscreen*rscreen) {
    assert(!are_smoothing);
    screening_function(r2,rscreen,&b0,&b1);
  } else {
    const double invr = 1.0/sqrt(r2);
    b0 = invr;
    b1 = invr*invr*invr;
  }
  if (are_smoothing) {
    *tmpsmooth += qi*qj*b0;
    b0 *= smooth;
    b1 *= smooth;
  }
  if (are_scaling) {
    b0 *= scale_factor;
    b1 *= scale_factor;
  }
  qib1 = qi*b1;
  qjb1 = qj*b1;
  phi[i] += qj*b0;
  evec[i].x += qjb1*x;
  evec[i].y += qjb1*y;
  evec[i].z += qjb1*z;
  phi[j] += qi*b0;
  evec[j].x -= qib1*x;
  evec[j].y -= qib1*y;
  evec[j].z -= qib1*z;
}

static void ewald(int i, int j,
		  double x, double y, double z,
		  double r2, double qi, double qj,
		  double *phi, cart_t *evec,
		  double eta,
		  const struct Chebyshev *rspace_chebyshev,
		  double rscreen, int are_smoothing,
		  double smooth, double *tmpsmooth,
		  int are_scaling, double scale_factor)
{
  double b0, b1, qib1, qjb1;
#ifdef EWCHEB
  if (rspace_chebyshev->min < r2 && r2 < rspace_chebyshev->max) {
    Chebyshev_evaluate(rspace_chebyshev, r2, &b0, &b1);
  } else {
#endif
    const double rm = sqrt(r2);
    const double etar = eta*rm;
    const double eta2r2 = sq(etar);
    const double expeta2r2 = exp(-eta2r2);
    b0 = erfc(etar)/rm;
    b1 = (b0 + eta*M_2_SQRTPI*expeta2r2)/r2;
#ifdef EWCHEB
  }
#endif
  if (r2 < rscreen*rscreen) {
    assert(!are_smoothing);
    const double invr = 1.0/sqrt(r2);
    double f, df;
    screening_function(r2,rscreen,&f,&df);
    b0 += f - invr;
    b1 += df - invr*invr*invr;
  }
  if (are_smoothing) {
    *tmpsmooth += qi*qj*b0;
    b0 *= smooth;
    b1 *= smooth;
  }
  if (are_scaling) {
    b0 *= scale_factor;
    b1 *= scale_factor;
  }
  qib1 = qi*b1;
  qjb1 = qj*b1;
  phi[i] += qj*b0;
  evec[i].x += qjb1*x;
  evec[i].y += qjb1*y;
  evec[i].z += qjb1*z;
  phi[j] += qi*b0;
  evec[j].x -= qib1*x;
  evec[j].y -= qib1*y;
  evec[j].z -= qib1*z;
}

static void coulomb_non_periodic_phi(int i, int j,
				     double x, double y, double z,
				     double r2, double qi, double qj,
				     double *phi,
				     double rscreen, int are_smoothing,
				     double smooth, int are_scaling,
				     double scale_factor)
{
  double b0, b1;
  if (r2 < rscreen*rscreen) {
    assert(!are_smoothing);
    screening_function(r2,rscreen,&b0,&b1);
  } else {
    const double invr = 1.0/sqrt(r2);
    b0 = invr;
  }
  if (are_smoothing)
    b0 *= smooth;
  if (are_scaling)
    b0 *= scale_factor;
  phi[i] += qj*b0;
  phi[j] += qi*b0;
}

static void ewald_phi(int i, int j,
		      double x, double y, double z,
		      double r2, double qi, double qj,
		      double *phi,
		      double eta,
		      const struct Chebyshev *rspace_chebyshev,
		      double rscreen, int are_smoothing,
		      double smooth, int are_scaling,
		      double scale_factor)
{
  double b0, b1;
#ifdef EWCHEB
  if (rspace_chebyshev->min < r2 && r2 < rspace_chebyshev->max) {
    Chebyshev_evaluate(rspace_chebyshev, r2, &b0, &b1);
  } else {
#endif
    const double rm = sqrt(r2);
    const double etar = eta*rm;
    b0 = erfc(etar)/rm;
#ifdef EWCHEB
  }
#endif
  if (r2 < rscreen*rscreen) {
    const double invr = 1.0/sqrt(r2);
    double f, df;
    assert(!are_smoothing);
    screening_function(r2,rscreen,&f,&df);
    b0 += f - invr;
  }
  if (are_smoothing)
    b0 *= smooth;
  if (are_scaling)
    b0 *= scale_factor;
  phi[i] += qj*b0;
  phi[j] += qi*b0;
}

static void ewald_excluded(int i, int j,
			     double x, double y, double z,
			     double r2,
			     double qi, double qj,
			     double *phi, cart_t *evec,
			     double eta)
{
  const double rm = sqrt(r2);
  const double etar = eta*rm;
  const double eta2r2 = sq(etar);
  const double expeta2r2 = exp(-eta2r2);
  const double b0 = -erf(etar)/rm;
  const double b1 = (b0 + eta*M_2_SQRTPI*expeta2r2)/r2;
  const double qib1 = qi*b1;
  const double qjb1 = qj*b1;
  phi[i] += qj*b0;
  evec[i].x += qjb1*x;
  evec[i].y += qjb1*y;
  evec[i].z += qjb1*z;
  phi[j] += qi*b0;
  evec[j].x -= qib1*x;
  evec[j].y -= qib1*y;
  evec[j].z -= qib1*z;
}

static void ewald_excluded_phi(int i, int j,
				 double x, double y, double z,
				 double r2,
				 double qi, double qj,
				 double *phi,
				 double eta)
{
  const double rm = sqrt(r2);
  const double etar = eta*rm;
  const double b0 = -erf(etar)/rm;
  phi[i] += qj*b0;
  phi[j] += qi*b0;
}

static void lj(int i, int j, 
	       double x, double y, double z, double r2,
	       const double *sig, const double *eps,
	       int geometric_combining_for_sigma,
	       double *u, cart_t *f,
	       int are_smoothing, double smooth, 
	       double *tmpsmooth, int are_scaling,
	       double scale_factor)
{
  const double s2 = geometric_combining_for_sigma ? sig[i]*sig[j] : sq(0.5*(sig[i]+sig[j]));
  const double s6 = s2*s2*s2;
  const double b = 4*sqrt(eps[i]*eps[j])*s6;
  const double a = b*s6;
  const double invr2 = 1.0/r2;
  const double invr6 = invr2*invr2*invr2;
  const double binvr6 = b*invr6;
  const double ainvr12 = a*invr6*invr6;
  double du = 12*ainvr12*invr2 - 6*binvr6*invr2;
  double utmp = ainvr12 - binvr6;
  double cx, cy, cz;
  if (are_smoothing) {
    *tmpsmooth += utmp;
    utmp *= smooth;
    du *= smooth;
  }
  if (are_scaling) {
    utmp *= scale_factor;
    du *= scale_factor;
  }
  *u += utmp;
  cx = du*x;
  f[i].x += cx;
  f[j].x -= cx;
  cy = du*y;
  f[i].y += cy;
  f[j].y -= cy;
  cz = du*z;
  f[i].z += cz;
  f[j].z -= cz;
}

static void interact(int i, int j,
		     double x, double y, double z, 
		     double r2, double qi, double qj,
		     double *phi, cart_t *evec,
		     int mask, const double *rscreen,
		     const double *sig, const double *eps,
		     int geometric_combining_for_sigma,
		     double *lj_ener, cart_t *lj_force,
		     int are_smoothing, double smooth, 
		     double *tmpsmooth,
		     int are_scaling_elec,
		     double elec_scale_factor,
		     int are_scaling_lj,
		     double lj_scale_factor)
{
  if (mask & SITE_HAS_Q)
    coulomb_non_periodic(i, j, x, y, z, r2, qi, qj, phi, evec, 
			 rscreen[i]+rscreen[j], are_smoothing, 
			 smooth, tmpsmooth, are_scaling_elec, elec_scale_factor);
  if (mask & SITE_HAS_LJ)
    lj(i, j, x, y, z, r2, sig, eps, geometric_combining_for_sigma, lj_ener, lj_force,
       are_smoothing, smooth, tmpsmooth, are_scaling_lj, lj_scale_factor);
}

static void interact_ewald(int i, int j,
			   double x, double y, double z, 
			   double r2, double qi, double qj,
			   double *phi, cart_t *evec,
			   double eta, int mask, 
			   const struct Chebyshev *rspace_chebyshev,
			   const double *rscreen,
			   const double *sig, const double *eps,
			   int geometric_combining_for_sigma,
			   double *lj_ener, cart_t *lj_force,
			   int are_smoothing, double smooth, 
			   double *tmpsmooth,
			   int are_scaling_elec,
			   double elec_scale_factor,
			   int are_scaling_lj,
			   double lj_scale_factor)
{
  if (mask & SITE_HAS_Q)
    ewald(i, j, x, y, z, r2, qi, qj, phi, evec, 
	  eta, rspace_chebyshev, rscreen[i]+rscreen[j], are_smoothing, 
	  smooth, tmpsmooth, are_scaling_elec, elec_scale_factor);
  if (mask & SITE_HAS_LJ)
    lj(i, j, x, y, z, r2, sig, eps, geometric_combining_for_sigma, lj_ener, lj_force,
       are_smoothing, smooth, tmpsmooth, are_scaling_lj, lj_scale_factor);
}

static void displacement_to_central_box(double x, double y, double z,
					double *dx, double *dy, double *dz,
					enum bcondtype_t bctype,
					double negbx, double invbx,
					double negby, double invby,
					double negbz, double invbz,
					double negtwobx, double invtwobx,
					const tensor_t *latv, 
					const tensor_t *invlatv)
{
  switch (bctype) {
  case cubic:
    cubic_displacement_to_central_box(x,y,z,dx,dy,dz,negbx,invbx);
    return;
  case orthorhombic:
    orthorhombic_displacement_to_central_box(x,y,z,dx,dy,dz,
					     negbx,negby,negbz,
					     invbx,invby,invbz);
    return;
  case bcc:
    bcc_displacement_to_central_box(x,y,z,dx,dy,dz,negtwobx,invtwobx);
    return;
  case fcc:
    fcc_displacement_to_central_box(x,y,z,dx,dy,dz,negbx,invbx);
    return;
  case triclinic:
    triclinic_displacement_to_central_box(x,y,z,dx,dy,dz,latv,invlatv);
    return;
  default:
    assert(0);
  }
}

void RSpaceEwald(int nsite, int ngroup, const int (*group)[2], const cart_t *r,
		 enum bcondtype_t bctype, const tensor_t *latv,
		 int *const *neighbors, const int *number_of_neighbors,
		 const struct PairBinary *group_excluded,
		 const struct PairBinary *excluded,
		 double rlo2, double rhi2,
		 const int *site_mask, const double *q, 
		 double *phi, cart_t *evec,
		 double eta, const struct Chebyshev *rspace_chebyshev,
		 const double *rscreen,
		 const double *sig, const double *eps,
		 int geometric_combining_for_sigma,
		 double *lj_ener, cart_t *lj_force,
		 const struct PairBinary *one_four,
		 double elec_1_4_scale_factor, 
		 double lj_1_4_scale_factor)
{
  const double negbx = -latv->xx;
  const double invbx = 1/latv->xx;
  const double negby = -latv->yy;
  const double invby = 1/latv->yy;
  const double negbz = -latv->zz;
  const double invbz = 1/latv->zz;
  const double negtwobx = -2*latv->xx;
  const double invtwobx = 1/(2*latv->xx);
  int i, ig;
  tensor_t invlatv;
  assert(bctype != non_periodic);
  assert(rhi2 >= rlo2);
  assert(rlo2 > 0);
  inverse(&invlatv, latv);
  for (ig = 0; ig < ngroup; ig++) {
    const int istart = group[ig][0];
    const int iend = group[ig][1];
    const int *jgp;
    const int *jgend = &neighbors[ig][number_of_neighbors[ig]];
    for (jgp = neighbors[ig]; jgp < jgend; jgp++) {
      const int jg = *jgp;
      const int jstart = group[jg][0];
      const int jend = group[jg][1];
      const int group_excl = PairBinary_exists(group_excluded,ig,jg);
      double x = r[istart].x - r[jstart].x;
      double y = r[istart].y - r[jstart].y;
      double z = r[istart].z - r[jstart].z;
      double dx, dy, dz, r2;
      double smooth=1, dsmooth=0, tmpsmooth=0, dsmx=0, dsmy=0, dsmz=0;
      int are_smoothing = 0;
      displacement_to_central_box(x,y,z,&dx,&dy,&dz,bctype,
				  negbx,invbx,negby,invby,negbz,invbz,
				  negtwobx,invtwobx,
				  latv,&invlatv);
      x += dx;
      y += dy;
      z += dz;
      r2 = x*x + y*y + z*z;
      if (r2 > rhi2)
	continue;
      if (r2 > rlo2) {
	are_smoothing = 1;
	switching_fn(r2, rlo2, rhi2, &smooth, &dsmooth);
	dsmx = dsmooth*x;
	dsmy = dsmooth*y;
	dsmz = dsmooth*z;
      }
      if (!group_excl || !PairBinary_exists(excluded,istart,jstart))
	interact_ewald(istart, jstart, x, y, z, r2, q[istart], q[jstart],
		       phi, evec, eta, site_mask[istart] & site_mask[jstart],
		       rspace_chebyshev, rscreen,
		       sig, eps, geometric_combining_for_sigma, 
		       lj_ener, lj_force,
		       are_smoothing, smooth, &tmpsmooth,
		       0, 1, 0, 1);
      for (i = istart; i < iend; i++) {
	const double xi = r[i].x + dx;
	const double yi = r[i].y + dy;
	const double zi = r[i].z + dz;
	const int imask = site_mask[i];
	const double qi = q[i];
	int j;
	for (j = (i == istart ? jstart+1 : jstart); j < jend; j++) {
	  const int mask = imask & site_mask[j];
	  if (!mask)
	    continue;
	  if (group_excl && PairBinary_exists(excluded,i,j))
	    continue;
	  x = xi - r[j].x;
	  y = yi - r[j].y;
	  z = zi - r[j].z;
	  r2 = x*x + y*y + z*z;
	  interact_ewald(i, j, x, y, z, r2, qi, q[j],
			 phi, evec, eta, mask, rspace_chebyshev, rscreen, 
			 sig, eps, geometric_combining_for_sigma,
			 lj_ener, lj_force,
			 are_smoothing, smooth, &tmpsmooth, 
			 0, 1, 0, 1);
	}
      }
      if (are_smoothing) {
	dsmx *= tmpsmooth;
	dsmy *= tmpsmooth;
	dsmz *= tmpsmooth;
	lj_force[istart].x += dsmx;
	lj_force[istart].y += dsmy;
	lj_force[istart].z += dsmz;
	lj_force[jstart].x -= dsmx;
	lj_force[jstart].y -= dsmy;
	lj_force[jstart].z -= dsmz;
      }
    }
  }
  /* Excluded */
  for (i = 0; i < nsite; i++) {
    int *const *excl = PairBinary_elements(excluded);
    const int *nexcl = PairBinary_number_of_elements(excluded);
    if (site_mask[i] & SITE_HAS_Q) {
      const double xi = r[i].x;
      const double yi = r[i].y;
      const double zi = r[i].z;
      const int *jp, *jend = &excl[i][nexcl[i]];
      const double qi = q[i];
      for (jp = excl[i]; jp < jend; jp++) {
	const int j = *jp;
	if (site_mask[j] & SITE_HAS_Q) {
	  double x = xi - r[j].x;
	  double y = yi - r[j].y;
	  double z = zi - r[j].z;
	  double dx, dy, dz;
	  displacement_to_central_box(x,y,z,&dx,&dy,&dz,bctype,
				      negbx,invbx,negby,invby,negbz,invbz,
				      negtwobx,invtwobx,
				      latv,&invlatv);
	  x += dx;
	  y += dy;
	  z += dz;
	  const double r2 = x*x + y*y + z*z;
	  ewald_excluded(i, j, x, y, z, r2, qi, q[j], phi, evec, eta);
	}
      }
    }
  }
  /* One-four */
  if (one_four) {
    const int are_scaling_elec = fabs(elec_1_4_scale_factor-1) > 1e-10;
    const int are_scaling_lj = fabs(lj_1_4_scale_factor-1) > 1e-10;
    for (i = 0; i < nsite; i++) {
      int *const *of = PairBinary_elements(one_four);
      const int *nof = PairBinary_number_of_elements(one_four);
      const double xi = r[i].x;
      const double yi = r[i].y;
      const double zi = r[i].z;
      const double qi = q[i];
      const int *jp, *jend = &of[i][nof[i]];
      const int imask = site_mask[i];
      for (jp = of[i]; jp < jend; jp++) {
	const int j = *jp;
	double x = xi - r[j].x;
	double y = yi - r[j].y;
	double z = zi - r[j].z;
	double dx, dy, dz;
	displacement_to_central_box(x,y,z,&dx,&dy,&dz,bctype,
				    negbx,invbx,negby,invby,negbz,invbz,
				    negtwobx,invtwobx,
				    latv,&invlatv);
	x += dx;
	y += dy;
	z += dz;
	const double r2 = x*x + y*y + z*z;
	interact(i, j, x, y, z, r2, qi, q[j],
		 phi, evec, imask & site_mask[j], 
		 rscreen, sig, eps, geometric_combining_for_sigma,
		 lj_ener, lj_force, 0, 0, 0, 
		 are_scaling_elec, elec_1_4_scale_factor - 1, 
		 are_scaling_lj, lj_1_4_scale_factor - 1);
      }
    }
  }
}

void RSpaceNonPeriodic(int nsite, int ngroup, const int (*group)[2], const cart_t *r,
		       const struct PairBinary *group_excluded,
		       const struct PairBinary *excluded,
		       const int *site_mask, const double *q, 
		       double *phi, cart_t *evec,
		       const double *rscreen,
		       const double *sig, const double *eps,
		       int geometric_combining_for_sigma,
		       double *lj_ener, cart_t *lj_force,
		       const struct PairBinary *one_four,
		       double elec_1_4_scale_factor,
		       double lj_1_4_scale_factor)
{
  int i, ig;
  for (ig = 0; ig < ngroup; ig++) {
    const int istart = group[ig][0];
    const int iend = group[ig][1];
    int jg;
    for (jg = ig+1; jg < ngroup; jg++) {
      const int jstart = group[jg][0];
      const int jend = group[jg][1];
      const int group_excl = PairBinary_exists(group_excluded,ig,jg);
      double x = r[istart].x - r[jstart].x;
      double y = r[istart].y - r[jstart].y;
      double z = r[istart].z - r[jstart].z;
      double r2 = x*x + y*y + z*z;
      if (!group_excl || !PairBinary_exists(excluded,istart,jstart))
	interact(istart, jstart, x, y, z, r2, q[istart], q[jstart],
		 phi, evec, site_mask[istart] & site_mask[jstart],
		 rscreen, sig, eps, geometric_combining_for_sigma,
		 lj_ener, lj_force, 0, 0, 0,
		 0, 1, 0, 1);
      for (i = istart; i < iend; i++) {
	const double xi = r[i].x;
	const double yi = r[i].y;
	const double zi = r[i].z;
	const int imask = site_mask[i];
	const double qi = q[i];
	int j;
	for (j = (i == istart ? jstart+1 : jstart); j < jend; j++) {
	  const int mask = imask & site_mask[j];
	  if (!mask)
	    continue;
	  if (group_excl && PairBinary_exists(excluded,i,j))
	    continue;
	  x = xi - r[j].x;
	  y = yi - r[j].y;
	  z = zi - r[j].z;
	  r2 = x*x + y*y + z*z;
	  interact(i, j, x, y, z, r2, qi, q[j],
		   phi, evec, mask, rscreen, 
		   sig, eps, geometric_combining_for_sigma,
		   lj_ener, lj_force, 0, 0, 0,
		   0, 1, 0, 1);
	}
      }
    }
  }
  /* One-four */
  if (one_four) {
    const int are_scaling_elec = fabs(elec_1_4_scale_factor-1) > 1e-10;
    const int are_scaling_lj = fabs(lj_1_4_scale_factor-1) > 1e-10;
    for (i = 0; i < nsite; i++) {
      int *const *of = PairBinary_elements(one_four);
      const int *nof = PairBinary_number_of_elements(one_four);
      const double xi = r[i].x;
      const double yi = r[i].y;
      const double zi = r[i].z;
      const double qi = q[i];
      const int *jp, *jend = &of[i][nof[i]];
      const int imask = site_mask[i];
      for (jp = of[i]; jp < jend; jp++) {
	const int j = *jp;
	const double x = xi - r[j].x;
	const double y = yi - r[j].y;
	const double z = zi - r[j].z;
	const double r2 = x*x + y*y + z*z;
	interact(i, j, x, y, z, r2, qi, q[j],
		 phi, evec, imask & site_mask[j], 
		 rscreen, sig, eps, geometric_combining_for_sigma,
		 lj_ener, lj_force,
		 0, 0, 0, 
		 are_scaling_elec, elec_1_4_scale_factor - 1, 
		 are_scaling_lj, lj_1_4_scale_factor - 1);
      }
    }
  }
}

void RSpaceEwaldPhi(int nsite, int ngroup, const int (*group)[2], const cart_t *r,
		    enum bcondtype_t bctype, const tensor_t *latv,
		    int *const *neighbors, const int *number_of_neighbors,
		    const struct PairBinary *group_excluded,
		    const struct PairBinary *excluded,
		    double rlo2, double rhi2,
		    const int *site_mask, const double *q, 
		    double *phi,
		    double eta, const struct Chebyshev *rspace_chebyshev,
		    const double *rscreen,
		    const struct PairBinary *one_four,
		    double elec_1_4_scale_factor)
{
  const double negbx = -latv->xx;
  const double invbx = 1/latv->xx;
  const double negby = -latv->yy;
  const double invby = 1/latv->yy;
  const double negbz = -latv->zz;
  const double invbz = 1/latv->zz;
  const double negtwobx = -2*latv->xx;
  const double invtwobx = 1/(2*latv->xx);
  int i, ig;
  tensor_t invlatv;
  assert(bctype != non_periodic);
  assert(rhi2 >= rlo2);
  assert(rlo2 > 0);
  inverse(&invlatv, latv);
  for (ig = 0; ig < ngroup; ig++) {
    const int istart = group[ig][0];
    const int iend = group[ig][1];
    const int *jgp;
    const int *jgend = &neighbors[ig][number_of_neighbors[ig]];
    for (jgp = neighbors[ig]; jgp < jgend; jgp++) {
      const int jg = *jgp;
      const int jstart = group[jg][0];
      const int jend = group[jg][1];
      const int group_excl = PairBinary_exists(group_excluded,ig,jg);
      double x = r[istart].x - r[jstart].x;
      double y = r[istart].y - r[jstart].y;
      double z = r[istart].z - r[jstart].z;
      double dx, dy, dz, r2;
      double smooth=1, dsmooth=0;
      int are_smoothing = 0;
      displacement_to_central_box(x,y,z,&dx,&dy,&dz,bctype,
				  negbx,invbx,negby,invby,negbz,invbz,
				  negtwobx,invtwobx,latv,&invlatv);
      x += dx;
      y += dy;
      z += dz;
      r2 = x*x + y*y + z*z;
      if (r2 > rhi2)
	continue;
      if (r2 > rlo2) {
	are_smoothing = 1;
	switching_fn(r2, rlo2, rhi2, &smooth, &dsmooth);
      }
      if ((!group_excl || !PairBinary_exists(excluded,istart,jstart))
	  && (site_mask[istart] & SITE_HAS_Q)
	  && (site_mask[jstart] & SITE_HAS_Q))
	ewald_phi(istart, jstart, x, y, z, r2, q[istart], q[jstart],
		  phi, eta, rspace_chebyshev, 
		  rscreen[istart]+rscreen[jstart], are_smoothing, smooth, 0, 1);
      for (i = istart; i < iend; i++) {
	const double xi = r[i].x + dx;
	const double yi = r[i].y + dy;
	const double zi = r[i].z + dz;
	const double qi = q[i];
	int j;
	if (!(site_mask[i] & SITE_HAS_Q))
	  continue;
	for (j = (i == istart ? jstart+1 : jstart); j < jend; j++) {
	  if (group_excl && PairBinary_exists(excluded,i,j))
	    continue;
	  if (!(site_mask[j] & SITE_HAS_Q))
	    continue;
	  x = xi - r[j].x;
	  y = yi - r[j].y;
	  z = zi - r[j].z;
	  r2 = x*x + y*y + z*z;
	  ewald_phi(i, j, x, y, z, r2, qi, q[j],
		    phi, eta, rspace_chebyshev,
		    rscreen[i]+rscreen[j], are_smoothing, smooth, 0, 1);
	}
      }
    }
  }
  for (i = 0; i < nsite; i++) {
    int *const *excl = PairBinary_elements(excluded);
    const int *nexcl = PairBinary_number_of_elements(excluded);
    if (site_mask[i] & SITE_HAS_Q) {
      const double xi = r[i].x;
      const double yi = r[i].y;
      const double zi = r[i].z;
      const int *jp, *jend = &excl[i][nexcl[i]];
      const double qi = q[i];
      for (jp = excl[i]; jp < jend; jp++) {
	const int j = *jp;
	if (site_mask[j] & SITE_HAS_Q) {
	  const double x = xi - r[j].x;
	  const double y = yi - r[j].y;
	  const double z = zi - r[j].z;
	  const double r2 = x*x + y*y + z*z;
	  ewald_excluded_phi(i, j, x, y, z, r2, qi, q[j], phi, eta);
	}
      }
    }
  }
  /* One-four */
  if (one_four) {
    const int are_scaling_elec = fabs(elec_1_4_scale_factor-1) > 1e-10;
    for (i = 0; i < nsite; i++)
      if (site_mask[i] & SITE_HAS_Q) {
	int *const *of = PairBinary_elements(one_four);
	const int *nof = PairBinary_number_of_elements(one_four);
	const double xi = r[i].x;
	const double yi = r[i].y;
	const double zi = r[i].z;
	const double qi = q[i];
	const int *jp, *jend = &of[i][nof[i]];
	for (jp = of[i]; jp < jend; jp++) {
	  const int j = *jp;
	  if (site_mask[j] & SITE_HAS_Q) {
	    const double x = xi - r[j].x;
	    const double y = yi - r[j].y;
	    const double z = zi - r[j].z;
	    const double r2 = x*x + y*y + z*z;
	    coulomb_non_periodic_phi(i, j, x, y, z, r2, qi, q[j],
				     phi, rscreen[i]+rscreen[j], 0, 0,
				     are_scaling_elec, elec_1_4_scale_factor - 1);
	  }
	}
      }
  }
}

void RSpaceNonPeriodicPhi(int nsite, int ngroup, const int (*group)[2], const cart_t *r,
			  const struct PairBinary *group_excluded,
			  const struct PairBinary *excluded,
			  const int *site_mask, const double *q, 
			  double *phi, const double *rscreen,
			  const struct PairBinary *one_four,
			  double elec_1_4_scale_factor)
{
  int i, ig;
  for (ig = 0; ig < ngroup; ig++) {
    const int istart = group[ig][0];
    const int iend = group[ig][1];
    int jg;
    for (jg = ig+1; jg < ngroup; jg++) {
      const int jstart = group[jg][0];
      const int jend = group[jg][1];
      const int group_excl = PairBinary_exists(group_excluded,ig,jg);
      double x = r[istart].x - r[jstart].x;
      double y = r[istart].y - r[jstart].y;
      double z = r[istart].z - r[jstart].z;
      double r2 = x*x + y*y + z*z;
      if (!group_excl || !PairBinary_exists(excluded,istart,jstart))
	coulomb_non_periodic_phi(istart, jstart, x, y, z, r2, q[istart], q[jstart],
				 phi, rscreen[istart]+rscreen[jstart], 0, 0, 0, 1);
      for (i = istart; i < iend; i++) {
	const double xi = r[i].x;
	const double yi = r[i].y;
	const double zi = r[i].z;
	const int imask = site_mask[i];
	const double qi = q[i];
	int j;
	for (j = (i == istart ? jstart+1 : jstart); j < jend; j++) {
	  const int mask = imask & site_mask[j];
	  if (!mask)
	    continue;
	  if (group_excl && PairBinary_exists(excluded,i,j))
	    continue;
	  x = xi - r[j].x;
	  y = yi - r[j].y;
	  z = zi - r[j].z;
	  r2 = x*x + y*y + z*z;
	  coulomb_non_periodic_phi(i, j, x, y, z, r2, qi, q[j],
				   phi, rscreen[i]+rscreen[j], 0, 0, 0, 1);
	}
      }
    }
  }
  if (one_four) {
    const int are_scaling_elec = fabs(elec_1_4_scale_factor-1) > 1;
    for (i = 0; i < nsite; i++)
      if (site_mask[i] & SITE_HAS_Q) {
	int *const *of = PairBinary_elements(one_four);
	const int *nof = PairBinary_number_of_elements(one_four);
	const double xi = r[i].x;
	const double yi = r[i].y;
	const double zi = r[i].z;
	const double qi = q[i];
	const int *jp, *jend = &of[i][nof[i]];
	for (jp = of[i]; jp < jend; jp++) {
	  const int j = *jp;	  
	  if (site_mask[j] & SITE_HAS_Q) {
	    const double x = xi - r[j].x;
	    const double y = yi - r[j].y;
	    const double z = zi - r[j].z;
	    const double r2 = x*x + y*y + z*z;
	    coulomb_non_periodic_phi(i, j, x, y, z, r2, qi, q[j],
				     phi, rscreen[i]+rscreen[j], 0, 0,
				     are_scaling_elec, elec_1_4_scale_factor - 1);
	  }
	}
      }
  }
}

void RSpaceNonPeriodicPhiAtGridpoints(int nsite, int ngrid, const cart_t *rsite,
				      const cart_t *rgrid,
				      const int *site_mask,
				      const double *q,
				      double *phi, const double *rscreen)
{
  int i;
  for (i = 0; i < nsite; i++) {
    if (!site_mask[i])
      continue;
    const double xi = rsite[i].x;
    const double yi = rsite[i].y;
    const double zi = rsite[i].z;
    const double qi = q[i];
    const double rscr = rscreen[i];
    const double rscr2 = rscr*rscr;
    int j;
    for (j = 0; j < ngrid; j++) {
      double b0, b1;
      const double x = xi - rgrid[j].x;
      const double y = yi - rgrid[j].y;
      const double z = zi - rgrid[j].z;
      const double r2 = x*x + y*y + z*z;
      if (r2 < rscr2) {
	screening_function(r2,rscr,&b0,&b1);
      } else {
	const double invr = 1.0/sqrt(r2);
	b0 = invr;
      }
      phi[j] += qi*b0;
    }
  }
}

void RSpaceEwaldPhiAtGridpoints(int nsite, int ngrid, const cart_t *rsite,
				const cart_t *rgrid, enum bcondtype_t bctype,
				const tensor_t *latv,
				double rlo2, double rhi2,
				const int *site_mask,
				const double *q,
				double *phi,
				double eta, const struct Chebyshev *rspace_chebyshev,
				const double *rscreen)
{
  const double negbx = -latv->xx;
  const double invbx = 1/latv->xx;
  const double negby = -latv->yy;
  const double invby = 1/latv->yy;
  const double negbz = -latv->zz;
  const double invbz = 1/latv->zz;
  const double negtwobx = -2*latv->xx;
  const double invtwobx = 1/(2*latv->xx);
  int i, j;
  tensor_t invlatv;
  assert(bctype != non_periodic);
  assert(rhi2 >= rlo2);
  assert(rlo2 > 0);
  inverse(&invlatv, latv);
  for (i = 0; i < nsite; i++) {
    if (!site_mask[i])
      continue;
    const double xi = rsite[i].x;
    const double yi = rsite[i].y;
    const double zi = rsite[i].z;
    const double qi = q[i];
    const double rscr = rscreen[i];
    const double rscr2 = rscr*rscr;
    for (j = 0; j < ngrid; j++) {
      double x = xi - rgrid[j].x;
      double y = yi - rgrid[j].y;
      double z = zi - rgrid[j].z;
      double dx, dy, dz, smooth=1, dsmooth=0;
      int are_smoothing = 0;
      displacement_to_central_box(x,y,z,&dx,&dy,&dz,bctype,
				  negbx,invbx,negby,invby,negbz,invbz,
				  negtwobx,invtwobx,latv,&invlatv);
      x += dx;
      y += dy;
      z += dz;
      const double r2 = x*x + y*y + z*z;
      if (r2 > rhi2)
	continue;
      if (r2 > rlo2) {
	are_smoothing = 1;
	switching_fn(r2, rlo2, rhi2, &smooth, &dsmooth);
      }
      double b0, b1;
#ifdef EWCHEB
      if (rspace_chebyshev->min < r2 && r2 < rspace_chebyshev->max) {
	Chebyshev_evaluate(rspace_chebyshev, r2, &b0, &b1);
      } else {
#endif
	const double rm = sqrt(r2);
	const double etar = eta*rm;
	b0 = erfc(etar)/rm;
#ifdef EWCHEB
      }
#endif
      if (r2 < rscr2) {
	const double invr = 1.0/sqrt(r2);
	double f, df;
	screening_function(r2,rscr,&f,&df);
	b0 += f - invr;
      }
      if (are_smoothing)
	b0 *= smooth;
      phi[j] += qi*b0;
    }
  }
}
