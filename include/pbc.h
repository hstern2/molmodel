#ifndef PBC_H
#define PBC_H

#include "cart.h"

static inline double round_to_nearest_int(double x) { return floor(x+0.5); }

static inline void cubic_map_to_central_box(double *x, double *y, double *z, 
					    double boxl, double invboxl)
{ 
  *x -= boxl * round_to_nearest_int((*x)*invboxl);
  *y -= boxl * round_to_nearest_int((*y)*invboxl);
  *z -= boxl * round_to_nearest_int((*z)*invboxl);
}

static inline void cubic_displacement_to_central_box(double x, double y, double z, 
						     double *dx, double *dy, double *dz,
						     double negboxl, double invboxl)
{ 
  *dx = negboxl * round_to_nearest_int(x*invboxl);
  *dy = negboxl * round_to_nearest_int(y*invboxl);
  *dz = negboxl * round_to_nearest_int(z*invboxl);
}

static inline void cubic_nearest_lattice_point(double x, double y, double z,
					       double boxl, int *i, int *j, int *k)
{
  *i = (int) round_to_nearest_int(x/boxl);
  *j = (int) round_to_nearest_int(y/boxl);
  *k = (int) round_to_nearest_int(z/boxl);
}

static inline void orthorhombic_map_to_central_box(double *x, double *y, double *z,
						   double bx, double by, double bz,
						   double invbx, double invby, double invbz)
{
  *x -= bx * round_to_nearest_int((*x)*invbx);
  *y -= by * round_to_nearest_int((*y)*invby);
  *z -= bz * round_to_nearest_int((*z)*invbz);
}

static inline void orthorhombic_displacement_to_central_box(double x, double y, double z, 
							    double *dx, double *dy, double *dz,
							    double negbx, double negby, double negbz,
							    double invbx, double invby, double invbz)
{ 
  *dx = negbx * round_to_nearest_int(x*invbx);
  *dy = negby * round_to_nearest_int(y*invby);
  *dz = negbz * round_to_nearest_int(z*invbz);
}

static inline void orthorhombic_nearest_lattice_point(double x, double y, double z,
						      double bx, double by, double bz,
						      int *i, int *j, int *k)
{
  *i = (int) round_to_nearest_int(x/bx);
  *j = (int) round_to_nearest_int(y/by);
  *k = (int) round_to_nearest_int(z/bz);
}

static inline void triclinic_map_to_central_box(double *x, double *y, double *z,
						const tensor_t *latv,
						const tensor_t *invlatv)
{
  cart_t a;
  a.x = *x;
  a.y = *y;
  a.z = *z;
  right_operate(invlatv,&a);
  a.x -= round_to_nearest_int(a.x);
  a.y -= round_to_nearest_int(a.y);
  a.z -= round_to_nearest_int(a.z);
  right_operate(latv,&a);
  *x = a.x;
  *y = a.y;
  *z = a.z;
}

static inline void triclinic_displacement_to_central_box(double x, double y, double z, 
							 double *dx, double *dy, double *dz,
							 const tensor_t *latv,
							 const tensor_t *invlatv)
{ 
  const double x0 = x, y0 = y, z0 = z;
  triclinic_map_to_central_box(&x,&y,&z,latv,invlatv);
  *dx = x-x0;
  *dy = y-y0;
  *dz = z-z0;
}

static inline void triclinic_nearest_lattice_point(double x, double y, double z,
						   const tensor_t *invlatv,
						   int *i, int *j, int *k)
{
  cart_t a;
  a.x = x;
  a.y = y;
  a.z = z;
  right_operate(invlatv,&a);
  *i = (int) round_to_nearest_int(a.x);
  *j = (int) round_to_nearest_int(a.y);
  *k = (int) round_to_nearest_int(a.z);
}

static inline void bcc_map_to_central_box(double *x, double *y, double *z,
					  double side, double invside)
{
  double tx, ty, tz;
  *x *= invside;
  *y *= invside;
  *z *= invside;
  tx = *x; 
  ty = *y; 
  tz = *z;
  *x -= round_to_nearest_int(*x);
  *y -= round_to_nearest_int(*y);
  *z -= round_to_nearest_int(*z);
  if (fabs(*x)+fabs(*y)+fabs(*z) > 0.75) {
    *x = tx - floor(tx) - 0.5;
    *y = ty - floor(ty) - 0.5;
    *z = tz - floor(tz) - 0.5;
  }
  *x *= side;
  *y *= side;
  *z *= side;
}

static inline void bcc_displacement_to_central_box(double x, double y, double z,
						   double *dx, double *dy, double *dz,
						   double negside, double invside)
{
  x *= invside;
  y *= invside;
  z *= invside;
  *dx = round_to_nearest_int(x);
  *dy = round_to_nearest_int(y);
  *dz = round_to_nearest_int(z);
  if (fabs(x-(*dx))+fabs(y-(*dy))+fabs(z-(*dz)) > 0.75) {
    *dx = floor(x) + 0.5;
    *dy = floor(y) + 0.5;
    *dz = floor(z) + 0.5;
  }
  *dx *= negside;
  *dy *= negside;
  *dz *= negside;
}

static inline void bcc_nearest_lattice_point(double x, double y, double z, double side,
					     int *i, int *j, int *k)
{
  const double sx = x/side, sy = y/side, sz = z/side;
  const double dx = round_to_nearest_int(sx);
  const double dy = round_to_nearest_int(sy);
  const double dz = round_to_nearest_int(sz);
  if (fabs(sx-dx)+fabs(sy-dy)+fabs(sz-dz) < 0.75) {
    const int ix = (int) dx;
    const int iy = (int) dy;
    const int iz = (int) dz;
    *i = ix + iy;
    *j = iz - ix;
    *k = iz - iy;
  } else {
    const int ix = (int) floor(sx);
    const int iy = (int) floor(sy);
    const int iz = (int) floor(sz);
    *i = ix + iy + 1;
    *j = iz - ix;
    *k = iz - iy;
  }
}

static inline void fcc_map_to_central_box(double *x, double *y, double *z,
					  double halfside, double invhalfside)
{
  double qx = round_to_nearest_int(*x *= invhalfside);
  double qy = round_to_nearest_int(*y *= invhalfside);
  double qz = round_to_nearest_int(*z *= invhalfside);
  int px = (int) qx;
  int py = (int) qy;
  int pz = (int) qz;
  if ((px + py + pz) % 2) {
    const double x_px = fabs((*x)-qx);
    const double y_py = fabs((*y)-qy);
    const double z_pz = fabs((*z)-qz);
    if (x_px > y_py && x_px > z_pz)
      qx = 2*((int) floor(*x)) + 1 - px;
    else if (y_py > z_pz)
      qy = 2*((int) floor(*y)) + 1 - py;
    else
      qz = 2*((int) floor(*z)) + 1 - pz;
  }
  *x = halfside*((*x) - qx);
  *y = halfside*((*y) - qy);
  *z = halfside*((*z) - qz);
}

static inline void fcc_displacement_to_central_box(double x, double y, double z,
						   double *dx, double *dy, double *dz,
						   double neghalfside, double invhalfside)
{
  double qx = round_to_nearest_int(x *= invhalfside);
  double qy = round_to_nearest_int(y *= invhalfside);
  double qz = round_to_nearest_int(z *= invhalfside);
  int px = (int) qx;
  int py = (int) qy;
  int pz = (int) qz;
  if ((px + py + pz) % 2) {
    const double x_px = fabs(x-qx);
    const double y_py = fabs(y-qy);
    const double z_pz = fabs(z-qz);
    if (x_px > y_py && x_px > z_pz)
      qx = 2*((int) floor(x)) + 1 - px;
    else if (y_py > z_pz)
      qy = 2*((int) floor(y)) + 1 - py;
    else
      qz = 2*((int) floor(z)) + 1 - pz;
  }
  *dx = neghalfside*qx;
  *dy = neghalfside*qy;
  *dz = neghalfside*qz;
}

static inline void fcc_nearest_lattice_point(double x, double y, double z, double halfside,
					     int *i, int *j, int *k)
{
  double qx = round_to_nearest_int(x /= halfside);
  double qy = round_to_nearest_int(y /= halfside);
  double qz = round_to_nearest_int(z /= halfside);
  int px = (int) qx;
  int py = (int) qy;
  int pz = (int) qz;
  if ((px + py + pz) % 2) {
    const double x_px = fabs(x-qx);
    const double y_py = fabs(y-qy);
    const double z_pz = fabs(z-qz);
    if (x_px > y_py && x_px > z_pz)
      px = 2*((int) floor(x)) + 1 - px;
    else if (y_py > z_pz)
      py = 2*((int) floor(y)) + 1 - py;
    else
      pz = 2*((int) floor(z)) + 1 - pz;
  }
  *i = (px + py - pz)/2;
  *j = (py - px + pz)/2;
  *k = (px - py + pz)/2;
}

#endif /* PBC_H */
