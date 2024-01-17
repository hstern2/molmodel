/* $Id: pppm.c,v 1.15 2008/09/24 14:07:38 hstern Exp $ */

/*
 * Copyright (c) 2008 Harry A. Stern
 * Copyright (c) 2008 University of Rochester
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */


#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_multimin.h>

#include "pppm.h"
#include "array.h"

struct p3m
{
  int nsite; /* number of sites */
  int nalias; /* limit for alias sum.  suggested value: 2 */
  int ik_differentiate; /* differentiate on k-space mesh (requires 3 extra FFT's) */
  int assignment_order; /* must be an integer > 0 */
  double grid_spacing; /* suggested grid spacing */
  tensor_t latv, invlatv; /* lattice vectors */
  icart_t ngrid; /* number of grid points in each dimension */
  cart_t ***kvec; /* grid-based k-space wave vectors */
  double ***green_hat; /* grid-based reciprocal-space Green function */
  double ***qphigrid; /* grid-based real-space charge, potential */
  cart_t **weight; /* weights for each particle for assigned grid points */
  cart_t **ndw; /* negative derivatives of weights */
  icart_t **index; /* for each particle, indices of nearest grid points */
  complex_t ***qphigrid_hat; /* grid-based reciprocal-space charge, potential */
  double ***evecgrid_x, ***evecgrid_y, ***evecgrid_z;
  complex_t ***evecgrid_hat_x, ***evecgrid_hat_y, ***evecgrid_hat_z;
  fftw_plan fft_q, fft_phi; /* forwards and backwards FFT for potential */
  fftw_plan fft_evec_x, fft_evec_y, fft_evec_z; /* FFTs for field */
};

/* Return n mod m; works for negative n */
static int mymod(int n, int m)
{
  int k;
  if (n >= 0)
    return n % m;
  k = -n % m;
  if (k > 0)
    return m - k;
  else
    return 0;
}

static void die_at(const char *s, const char *file, int line)
{
  printf("%s at %s: %d\n", s, file, line);
  abort();
}

#define insist(x) if (!(x)) die_at("insist failed: " #x,__FILE__,__LINE__)

static double sq(double x) { return x*x; }

static double sincpi(double x) { 
  const double y  = M_PI*x;
  return fabs(y) < 1e-16 ? 1 - sq(y)/6 : sin(y)/y;
}

#define MAX_ASSIGNMENT_ORDER 7
#define POLY2(x,a,b) (a+b*x)
#define POLY3(x,a,b,c) (a+POLY2(x,b,c)*x)
#define POLY4(x,a,b,c,d) (a+POLY3(x,b,c,d)*x)
#define POLY5(x,a,b,c,d,e) (a+POLY4(x,b,c,d,e)*x)
#define POLY6(x,a,b,c,d,e,f) (a+POLY5(x,b,c,d,e,f)*x)
#define POLY7(x,a,b,c,d,e,f,g) (a+POLY6(x,b,c,d,e,f,g)*x)

static void weight(const cart_t *r, cart_t *w, int order)
{
  switch (order) {
  case 1:
    w[0].x = 1;
    w[0].y = 1;
    w[0].z = 1;
    return;
  case 2:
    w[0].x = POLY2(r->x,1,-2)/2;
    w[0].y = POLY2(r->y,1,-2)/2;
    w[0].z = POLY2(r->z,1,-2)/2;
    w[1].x = POLY2(r->x,1,2)/2;
    w[1].y = POLY2(r->y,1,2)/2;
    w[1].z = POLY2(r->z,1,2)/2;
    return;
  case 3:
    w[0].x = POLY3(r->x,1,-4,4)/8;
    w[0].y = POLY3(r->y,1,-4,4)/8;
    w[0].z = POLY3(r->z,1,-4,4)/8;
    w[1].x = POLY3(r->x,3,0,-4)/4;
    w[1].y = POLY3(r->y,3,0,-4)/4;
    w[1].z = POLY3(r->z,3,0,-4)/4;
    w[2].x = POLY3(r->x,1,4,4)/8;
    w[2].y = POLY3(r->y,1,4,4)/8;
    w[2].z = POLY3(r->z,1,4,4)/8;
    return;
  case 4:
    w[0].x = POLY4(r->x,1,-6,12,-8)/48;
    w[0].y = POLY4(r->y,1,-6,12,-8)/48;
    w[0].z = POLY4(r->z,1,-6,12,-8)/48;
    w[1].x = POLY4(r->x,23,-30,-12,24)/48;
    w[1].y = POLY4(r->y,23,-30,-12,24)/48;
    w[1].z = POLY4(r->z,23,-30,-12,24)/48;
    w[2].x = POLY4(r->x,23,30,-12,-24)/48;
    w[2].y = POLY4(r->y,23,30,-12,-24)/48;
    w[2].z = POLY4(r->z,23,30,-12,-24)/48;
    w[3].x = POLY4(r->x,1,6,12,8)/48;
    w[3].y = POLY4(r->y,1,6,12,8)/48;
    w[3].z = POLY4(r->z,1,6,12,8)/48;
    return;
  case 5:
    w[0].x = POLY5(r->x,1,-8,24,-32,16)/384;
    w[0].y = POLY5(r->y,1,-8,24,-32,16)/384;
    w[0].z = POLY5(r->z,1,-8,24,-32,16)/384;
    w[1].x = POLY5(r->x,19,-44,24,16,-16)/96;
    w[1].y = POLY5(r->y,19,-44,24,16,-16)/96;
    w[1].z = POLY5(r->z,19,-44,24,16,-16)/96;
    w[2].x = POLY5(r->x,115,0,-120,0,48)/192;
    w[2].y = POLY5(r->y,115,0,-120,0,48)/192;
    w[2].z = POLY5(r->z,115,0,-120,0,48)/192;
    w[3].x = POLY5(r->x,19,44,24,-16,-16)/96;
    w[3].y = POLY5(r->y,19,44,24,-16,-16)/96;
    w[3].z = POLY5(r->z,19,44,24,-16,-16)/96;
    w[4].x = POLY5(r->x,1,8,24,32,16)/384;
    w[4].y = POLY5(r->y,1,8,24,32,16)/384;
    w[4].z = POLY5(r->z,1,8,24,32,16)/384;
    return;
  case 6:
    w[0].x = POLY6(r->x,1,-10,40,-80,+80,-32)/3840;
    w[0].y = POLY6(r->y,1,-10,40,-80,+80,-32)/3840;
    w[0].z = POLY6(r->z,1,-10,40,-80,+80,-32)/3840;
    w[1].x = POLY6(r->x,237,-750,840,-240,-240,160)/3840;
    w[1].y = POLY6(r->y,237,-750,840,-240,-240,160)/3840;
    w[1].z = POLY6(r->z,237,-750,840,-240,-240,160)/3840;
    w[2].x = POLY6(r->x,841,-770,-440,560,80,-160)/1920;
    w[2].y = POLY6(r->y,841,-770,-440,560,80,-160)/1920;
    w[2].z = POLY6(r->z,841,-770,-440,560,80,-160)/1920;
    w[3].x = POLY6(r->x,841,770,-440,-560,80,160)/1920;
    w[3].y = POLY6(r->y,841,770,-440,-560,80,160)/1920;
    w[3].z = POLY6(r->z,841,770,-440,-560,80,160)/1920;
    w[4].x = POLY6(r->x,237,750,840,240,-240,-160)/3840;
    w[4].y = POLY6(r->y,237,750,840,240,-240,-160)/3840;
    w[4].z = POLY6(r->z,237,750,840,240,-240,-160)/3840;
    w[5].x = POLY6(r->x,1,10,40,80,80,32)/3840;
    w[5].y = POLY6(r->y,1,10,40,80,80,32)/3840;
    w[5].z = POLY6(r->z,1,10,40,80,80,32)/3840;
    return;
  case 7:
    w[0].x = POLY7(r->x,1,-12,60,-160,240,-192,64)/46080;
    w[0].y = POLY7(r->y,1,-12,60,-160,240,-192,64)/46080;
    w[0].z = POLY7(r->z,1,-12,60,-160,240,-192,64)/46080;
    w[1].x = POLY7(r->x,361,-1416,2220,-1600,240,384,-192)/23040;
    w[1].y = POLY7(r->y,361,-1416,2220,-1600,240,384,-192)/23040;
    w[1].z = POLY7(r->z,361,-1416,2220,-1600,240,384,-192)/23040;
    w[2].x = POLY7(r->x,10543,-17340,4740,6880,-4080,-960,960)/46080;
    w[2].y = POLY7(r->y,10543,-17340,4740,6880,-4080,-960,960)/46080;
    w[2].z = POLY7(r->z,10543,-17340,4740,6880,-4080,-960,960)/46080;
    w[3].x = POLY7(r->x,5887,0,-4620,0,1680,0,-320)/11520;
    w[3].y = POLY7(r->y,5887,0,-4620,0,1680,0,-320)/11520;
    w[3].z = POLY7(r->z,5887,0,-4620,0,1680,0,-320)/11520;
    w[4].x = POLY7(r->x,10543,17340,4740,-6880,-4080,960,960)/46080;    
    w[4].y = POLY7(r->y,10543,17340,4740,-6880,-4080,960,960)/46080;    
    w[4].z = POLY7(r->z,10543,17340,4740,-6880,-4080,960,960)/46080;    
    w[5].x = POLY7(r->x,361,1416,2220,1600,240,-384,-192)/23040;
    w[5].y = POLY7(r->y,361,1416,2220,1600,240,-384,-192)/23040;
    w[5].z = POLY7(r->z,361,1416,2220,1600,240,-384,-192)/23040;
    w[6].x = POLY7(r->x,1,12,60,160,240,192,64)/46080;
    w[6].y = POLY7(r->y,1,12,60,160,240,192,64)/46080;
    w[6].z = POLY7(r->z,1,12,60,160,240,192,64)/46080;
  default:
    assert(0);
  }
}

static void neg_weight_deriv(const cart_t *r, cart_t *w, int order)
{
  switch (order) {
  case 1:
    insist(0);
    return;
  case 2:
    w[0].x = 1;
    w[0].y = 1;
    w[0].z = 1;
    w[1].x = -1;
    w[1].y = -1;
    w[1].z = -1;
    return;
  case 3:
    w[0].x = -POLY2(r->x,-4,8)/8;
    w[0].y = -POLY2(r->y,-4,8)/8;
    w[0].z = -POLY2(r->z,-4,8)/8;
    w[1].x = -POLY2(r->x,0,-8)/4;
    w[1].y = -POLY2(r->y,0,-8)/4;
    w[1].z = -POLY2(r->z,0,-8)/4;
    w[2].x = -POLY2(r->x,4,8)/8;
    w[2].y = -POLY2(r->y,4,8)/8;
    w[2].z = -POLY2(r->z,4,8)/8;
    return;
  case 4:
    w[0].x = -POLY3(r->x,-6,24,-24)/48;
    w[0].y = -POLY3(r->y,-6,24,-24)/48;
    w[0].z = -POLY3(r->z,-6,24,-24)/48;
    w[1].x = -POLY3(r->x,-30,-24,72)/48;
    w[1].y = -POLY3(r->y,-30,-24,72)/48;
    w[1].z = -POLY3(r->z,-30,-24,72)/48;
    w[2].x = -POLY3(r->x,30,-24,-72)/48;
    w[2].y = -POLY3(r->y,30,-24,-72)/48;
    w[2].z = -POLY3(r->z,30,-24,-72)/48;
    w[3].x = -POLY3(r->x,6,24,24)/48;
    w[3].y = -POLY3(r->y,6,24,24)/48;
    w[3].z = -POLY3(r->z,6,24,24)/48;
    return;
  case 5:
    w[0].x = -POLY4(r->x,-8,48,-96,64)/384;
    w[0].y = -POLY4(r->y,-8,48,-96,64)/384;
    w[0].z = -POLY4(r->z,-8,48,-96,64)/384;
    w[1].x = -POLY4(r->x,-44,48,48,-64)/96;
    w[1].y = -POLY4(r->y,-44,48,48,-64)/96;
    w[1].z = -POLY4(r->z,-44,48,48,-64)/96;
    w[2].x = -POLY4(r->x,0,-240,0,192)/192;
    w[2].y = -POLY4(r->y,0,-240,0,192)/192;
    w[2].z = -POLY4(r->z,0,-240,0,192)/192;
    w[3].x = -POLY4(r->x,44,48,-48,-64)/96;
    w[3].y = -POLY4(r->y,44,48,-48,-64)/96;
    w[3].z = -POLY4(r->z,44,48,-48,-64)/96;
    w[4].x = -POLY4(r->x,8,48,96,64)/384;
    w[4].y = -POLY4(r->y,8,48,96,64)/384;
    w[4].z = -POLY4(r->z,8,48,96,64)/384;
    return;
  case 6:
    w[0].x = -POLY5(r->x,-10,80,-240,320,-160)/3840;
    w[0].y = -POLY5(r->y,-10,80,-240,320,-160)/3840;
    w[0].z = -POLY5(r->z,-10,80,-240,320,-160)/3840;
    w[1].x = -POLY5(r->x,-750,1680,-720,-960,800)/3840;
    w[1].y = -POLY5(r->y,-750,1680,-720,-960,800)/3840;
    w[1].z = -POLY5(r->z,-750,1680,-720,-960,800)/3840;
    w[2].x = -POLY5(r->x,-770,-880,1680,320,-800)/1920;
    w[2].y = -POLY5(r->y,-770,-880,1680,320,-800)/1920;
    w[2].z = -POLY5(r->z,-770,-880,1680,320,-800)/1920;
    w[3].x = -POLY5(r->x,770,-880,-1680,320,800)/1920;
    w[3].y = -POLY5(r->y,770,-880,-1680,320,800)/1920;
    w[3].z = -POLY5(r->z,770,-880,-1680,320,800)/1920;
    w[4].x = -POLY5(r->x,750,1680,720,-960,-800)/3840;
    w[4].y = -POLY5(r->y,750,1680,720,-960,-800)/3840;
    w[4].z = -POLY5(r->z,750,1680,720,-960,-800)/3840;
    w[5].x = -POLY5(r->x,10,80,240,320,160)/3840;
    w[5].y = -POLY5(r->y,10,80,240,320,160)/3840;
    w[5].z = -POLY5(r->z,10,80,240,320,160)/3840;
    return;
  case 7:
    w[0].x = -POLY6(r->x,-12,120,-480,960,-960,384)/46080;
    w[0].y = -POLY6(r->y,-12,120,-480,960,-960,384)/46080;
    w[0].z = -POLY6(r->z,-12,120,-480,960,-960,384)/46080;
    w[1].x = -POLY6(r->x,-1416,4440,-4800,960,1920,-1152)/23040;
    w[1].y = -POLY6(r->y,-1416,4440,-4800,960,1920,-1152)/23040;
    w[1].z = -POLY6(r->z,-1416,4440,-4800,960,1920,-1152)/23040;
    w[2].x = -POLY6(r->x,-17340,9480,20640,-16320,-4800,5760)/46080;
    w[2].y = -POLY6(r->y,-17340,9480,20640,-16320,-4800,5760)/46080;
    w[2].z = -POLY6(r->z,-17340,9480,20640,-16320,-4800,5760)/46080;
    w[3].x = -POLY6(r->x,0,-9240,0,6720,0,-1920)/11520;
    w[3].y = -POLY6(r->y,0,-9240,0,6720,0,-1920)/11520;
    w[3].z = -POLY6(r->z,0,-9240,0,6720,0,-1920)/11520;
    w[4].x = -POLY6(r->x,17340,9480,-20640,-16320,4800,5760)/46080;
    w[4].y = -POLY6(r->y,17340,9480,-20640,-16320,4800,5760)/46080;
    w[4].z = -POLY6(r->z,17340,9480,-20640,-16320,4800,5760)/46080;
    w[5].x = -POLY6(r->x,1416,4440,4800,960,-1920,-1152)/23040;
    w[5].y = -POLY6(r->y,1416,4440,4800,960,-1920,-1152)/23040;
    w[5].z = -POLY6(r->z,1416,4440,4800,960,-1920,-1152)/23040;
    w[6].x = -POLY6(r->x,12,120,480,960,960,384)/46080;
    w[6].y = -POLY6(r->y,12,120,480,960,960,384)/46080;
    w[6].z = -POLY6(r->z,12,120,480,960,960,384)/46080;
  default:
    assert(0);
  }
}

static int multiple_of_2357_near(double x)
{
  int p2, p23, p235, p2357, best = 1;
  double bestdiff = fabs(x-1);
  for (p2 = 1; p2 <= x; p2 *= 2)
    for (p23 = p2; p23 <= x; p23 *= 3)
      for (p235 = p23; p235 <= x; p235 *= 5)
	for (p2357 = p235; p2357 <= x; p2357 *= 7) {
	  const double diff = fabs(x - p2357);
	  if (fabs(diff) < bestdiff) {
	    bestdiff = diff;
	    best = p2357;
	  }
	}
  return best;
}

p3m_t p3m_copy(const p3m_t that, int nsite)
{
  int ngrid[3];
  p3m_t p = p3m_new(nsite, 
		    (double *) &that->latv,
		    that->grid_spacing,
		    that->nalias,
		    that->ik_differentiate,
		    that->assignment_order,
		    ngrid);
  assert(ngrid[0] == that->ngrid.x);
  assert(ngrid[1] == that->ngrid.y);
  assert(ngrid[2] == that->ngrid.z);
  const int n = ngrid[0]*ngrid[1]*(ngrid[2]/2+1);
  memcpy(p->kvec, that->kvec, n*sizeof(cart_t));
  memcpy(p->green_hat, that->green_hat, n*sizeof(double));
  return p;
} 

void p3m_init(p3m_t that, double eta, double *chi)
{
  const int nx = that->ngrid.x;
  const int ny = that->ngrid.y;
  const int nz = that->ngrid.z;
  const double fourpi = 4*M_PI;
  const double inv_4eta2 = 1/sq(2*eta);
  const double inv_v = 1.0/fabs(det(&that->latv));
  double ***green_hat = that->green_hat;
  const int nalias = that->nalias;
  const int order = that->assignment_order;
  int ix, iy, iz;
  *chi = 0;
  tensor_t reciprocal_latv;
  cart_t *uv2 = cart_t_array_new(2*nalias+1) + nalias;
  transpose(&reciprocal_latv, &that->invlatv);
  scalar_multiply(&reciprocal_latv, 2*M_PI);
  /* Compute grid-based reciprocal-space Green function */
  for (ix = 0; ix < nx; ix++) {
    cart_t j;
    j.x = ix > nx/2 ? ix - nx : ix;
    for (iy = 0; iy < ny; iy++) {
      j.y = iy > ny/2 ? iy - ny : iy;
      for (iz = 0; iz <= nz/2; iz++) {
	cart_t kvec, p;
	double kvec2, numer = 0, denom1 = 0, denom2 = 0, denom = 0;
	int q, px, py, pz;
	j.z = iz;
	multiply(&kvec, &reciprocal_latv, &j);
	that->kvec[ix][iy][iz] = kvec;
	if (ix == 0 && iy == 0 && iz == 0) {
	  green_hat[0][0][0] = 0;
	  continue;
	}
	kvec2 = cart_sq(&kvec);
	for (q = -nalias; q <= nalias; q++) {
	  uv2[q].x = pow(sincpi(j.x/nx+q), 2*order);
	  uv2[q].y = pow(sincpi(j.y/ny+q), 2*order);
	  uv2[q].z = pow(sincpi(j.z/nz+q), 2*order);
	}
	for (px = -nalias; px <= nalias; px++) {
	  p.x = j.x + px*nx;
	  for (py = -nalias; py <= nalias; py++) {
	    p.y = j.y + py*ny;
	    for (pz = -nalias; pz <= nalias; pz++) {
	      double kp2, u2, g, g2;
	      cart_t kp;	      
	      p.z = j.z + pz*nz;
	      if (fabs(p.x) < 1e-16 && 
		  fabs(p.y) < 1e-16 &&
		  fabs(p.z) < 1e-16)
		continue;
	      multiply(&kp, &reciprocal_latv, &p);
	      kp2 = cart_sq(&kp);
	      u2 = uv2[px].x * uv2[py].y * uv2[pz].z;
	      g = fourpi/kp2 * exp(-kp2*inv_4eta2);
	      g2 = g*g;
	      if (that->ik_differentiate) {
		numer += g*u2 * dot(&kvec,&kp);
		denom1 += u2;
	      } else {
		numer += g*kp2*u2;
		denom1 += u2;
		denom2 += u2*kp2;
	      }
	      *chi += g2*kp2;
	    }
	  }
	}
	if (that->ik_differentiate)
	  denom = sq(denom1) * kvec2;
	else
	  denom = denom1 * denom2;
	green_hat[ix][iy][iz] = inv_v * numer/denom;
	*chi -= sq(numer)/denom;
      }
    }
  }
  cart_t_array_delete(uv2 - nalias);
  *chi = pow(inv_v,1.0/3.0) * sqrt(2*fabs(*chi)); /* factor of 2 from +/- k */
}

p3m_t p3m_new(int nsite, 
	      const double primitive_translation_vectors[9],
	      double grid_spacing, 
	      int nalias, 
	      int ik_differentiate, 
	      int assignment_order, 
	      int ngrid[3])
{
  int nx, ny, nz;
  double a1, a2, a3;
  p3m_t that = malloc(sizeof(struct p3m));
  memset(that, 0, sizeof(struct p3m));
  insist(nsite >= 0);
  insist(assignment_order > 1);
  insist(assignment_order <= MAX_ASSIGNMENT_ORDER);
  insist(assignment_order >= 2 || ik_differentiate);
  insist(nalias > 0);
  insist(grid_spacing > 0);
  that->nsite = nsite;
  that->nalias = nalias;
  that->grid_spacing = grid_spacing;
  that->assignment_order = assignment_order;
  that->ik_differentiate = ik_differentiate;
  that->weight = cart_t_array2D_new(nsite, assignment_order);
  that->index = icart_t_array2D_new(nsite, assignment_order);
  that->kvec = 0;
  const tensor_t *lattice_vectors = (const tensor_t *) primitive_translation_vectors;
  if (!ik_differentiate)
    that->ndw = cart_t_array2D_new(nsite, assignment_order);
  a1 = sqrt(sq(lattice_vectors->xx) + sq(lattice_vectors->yx) + sq(lattice_vectors->zx));
  a2 = sqrt(sq(lattice_vectors->xy) + sq(lattice_vectors->yy) + sq(lattice_vectors->zy));
  a3 = sqrt(sq(lattice_vectors->xz) + sq(lattice_vectors->yz) + sq(lattice_vectors->zz));
  that->latv = *lattice_vectors;
  inverse(&that->invlatv, lattice_vectors);
  nx = that->ngrid.x = multiple_of_2357_near(floor(a1/grid_spacing + 0.5));
  ny = that->ngrid.y = multiple_of_2357_near(floor(a2/grid_spacing + 0.5));
  nz = that->ngrid.z = multiple_of_2357_near(floor(a3/grid_spacing + 0.5));
  that->qphigrid = double_array3D_new(nx, ny, nz);
  that->qphigrid_hat = complex_t_array3D_new(nx, ny, nz/2 + 1);
  that->fft_q = fftw_plan_dft_r2c_3d(nx, ny, nz, 
				     that->qphigrid[0][0],
				     (fftw_complex *) that->qphigrid_hat[0][0], 
				     FFTW_MEASURE|FFTW_DESTROY_INPUT);
  that->fft_phi = fftw_plan_dft_c2r_3d(nx, ny, nz,
				       (fftw_complex *) that->qphigrid_hat[0][0], 
				       that->qphigrid[0][0], 
				       FFTW_MEASURE|FFTW_DESTROY_INPUT);
  that->kvec = cart_t_array3D_new(nx, ny, nz/2 + 1);
  that->green_hat = double_array3D_new(nx, ny, nz/2 + 1);
  if (that->ik_differentiate) {
    that->evecgrid_x = double_array3D_new(nx, ny, nz);
    that->evecgrid_y = double_array3D_new(nx, ny, nz);
    that->evecgrid_z = double_array3D_new(nx, ny, nz);
    that->evecgrid_hat_x = complex_t_array3D_new(nx, ny, nz/2 + 1);
    that->evecgrid_hat_y = complex_t_array3D_new(nx, ny, nz/2 + 1);
    that->evecgrid_hat_z = complex_t_array3D_new(nx, ny, nz/2 + 1);
    that->fft_evec_x = fftw_plan_dft_c2r_3d(nx, ny, nz,
					    (fftw_complex *) that->evecgrid_hat_x[0][0], 
					    that->evecgrid_x[0][0], 
					    FFTW_MEASURE|FFTW_DESTROY_INPUT);
    that->fft_evec_y = fftw_plan_dft_c2r_3d(nx, ny, nz,
					    (fftw_complex *) that->evecgrid_hat_y[0][0], 
					    that->evecgrid_y[0][0], 
					    FFTW_MEASURE|FFTW_DESTROY_INPUT);
    that->fft_evec_z = fftw_plan_dft_c2r_3d(nx, ny, nz,
					    (fftw_complex *) that->evecgrid_hat_z[0][0], 
					    that->evecgrid_z[0][0], 
					    FFTW_MEASURE|FFTW_DESTROY_INPUT);
  }
  ngrid[0] = that->ngrid.x;
  ngrid[1] = that->ngrid.y;
  ngrid[2] = that->ngrid.z;
  return that;
}

void p3m_delete(p3m_t that)
{
  cart_t_array2D_delete(that->weight);
  icart_t_array2D_delete(that->index);
  if (!that->ik_differentiate)
    cart_t_array2D_delete(that->ndw);
  if (that->kvec) {
    double_array3D_delete(that->qphigrid);
    complex_t_array3D_delete(that->qphigrid_hat);
    fftw_destroy_plan(that->fft_q);
    fftw_destroy_plan(that->fft_phi);
    cart_t_array3D_delete(that->kvec);
    double_array3D_delete(that->green_hat);
    if (that->ik_differentiate) {
      double_array3D_delete(that->evecgrid_x);
      double_array3D_delete(that->evecgrid_y);
      double_array3D_delete(that->evecgrid_z);
      complex_t_array3D_delete(that->evecgrid_hat_x);
      complex_t_array3D_delete(that->evecgrid_hat_y);
      complex_t_array3D_delete(that->evecgrid_hat_z);
      fftw_destroy_plan(that->fft_evec_x);
      fftw_destroy_plan(that->fft_evec_y);
      fftw_destroy_plan(that->fft_evec_z);
    }
  }
  free(that);
}

void p3m_field(p3m_t that, 
	       const double *rp, 
	       const double *q, 
	       double *phi, 
	       double *evecp)
{
  const cart_t *r = (const cart_t *) rp;
  cart_t *evec = (cart_t *) evecp;
  const int nx = that->ngrid.x;
  const int ny = that->ngrid.y;
  const int nz = that->ngrid.z;
  const int nsite = that->nsite;
  const int ik_differentiate = that->ik_differentiate;
  const int order = that->assignment_order;
  double ***qphigrid = that->qphigrid;
  complex_t ***qphigrid_hat = that->qphigrid_hat;
  double ***green_hat = 0;
  double ***ex = that->evecgrid_x;
  double ***ey = that->evecgrid_y;
  double ***ez = that->evecgrid_z;
  tensor_t invlatv = that->invlatv;
  int i, j, k;
  const int odd = order % 2;
  const int istart = -((order+1)/2 - 1);
  green_hat = that->green_hat;
  invlatv.xx *= nx;  invlatv.xy *= nx;  invlatv.xz *= nx;
  invlatv.yx *= ny;  invlatv.yy *= ny;  invlatv.yz *= ny;
  invlatv.zx *= nz;  invlatv.zy *= nz;  invlatv.zz *= nz;
  /* Compute weights and assign charges to grid */
  memset(qphigrid[0][0], 0, nx*ny*nz*sizeof(double));
  for (i = 0; i < nsite; i++) {
    icart_t *index = that->index[i];
    cart_t s;
    int p, jx, jy, jz;
    cart_t *w = that->weight[i];
    multiply(&s, &invlatv, &r[i]);
    if (odd) {
      jx = (int) floor(s.x + 0.5);
      jy = (int) floor(s.y + 0.5);
      jz = (int) floor(s.z + 0.5);
      s.x -= jx;
      s.y -= jy;
      s.z -= jz;
    } else {
      jx = (int) floor(s.x);
      jy = (int) floor(s.y);
      jz = (int) floor(s.z);
      s.x -= jx + 0.5;
      s.y -= jy + 0.5;
      s.z -= jz + 0.5;
    }
    weight(&s, w, order);
    if (evec && !ik_differentiate)
      neg_weight_deriv(&s, that->ndw[i], order);
    for (p = 0; p < order; p++) {
      const int m = istart+p;
      index[p].x = mymod(jx+m, nx);
      index[p].y = mymod(jy+m, ny);
      index[p].z = mymod(jz+m, nz);
    }
    for (jx = 0; jx < order; jx++) {
      for (jy = 0; jy < order; jy++) {
	const double qwxwy = w[jx].x*w[jy].y*q[i];
	for (jz = 0; jz < order; jz++)
	  qphigrid[index[jx].x][index[jy].y][index[jz].z] += qwxwy * w[jz].z;
      }
    }
  }
  /* Transform grid charges to reciprocal space */
  fftw_execute(that->fft_q);
  /* Multiply by Green function to get grid potential in reciprocal space */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k <= nz/2; k++) {
	const double a = green_hat[i][j][k];
	complex_t *b = &qphigrid_hat[i][j][k];
	b->re *= a;
	b->im *= a;
      }
  if (evec && ik_differentiate) {
    /* Multiply by -ik to get grid-based field */
    cart_t ***kvec = that->kvec;
    complex_t ***ehx = that->evecgrid_hat_x;
    complex_t ***ehy = that->evecgrid_hat_y;
    complex_t ***ehz = that->evecgrid_hat_z;
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
	for (k = 0; k <= nz/2; k++) {
	  const cart_t kv = kvec[i][j][k];
	  const complex_t a = qphigrid_hat[i][j][k];
	  ehx[i][j][k].re = kv.x * a.im;
	  ehx[i][j][k].im = -kv.x * a.re;
	  ehy[i][j][k].re = kv.y * a.im;
	  ehy[i][j][k].im = -kv.y * a.re;
	  ehz[i][j][k].re = kv.z * a.im;
	  ehz[i][j][k].im = -kv.z * a.re;
	}
    /* Transform grid-based field back to real space */
    fftw_execute(that->fft_evec_x);
    fftw_execute(that->fft_evec_y);
    fftw_execute(that->fft_evec_z);
  }
  /* Transform grid potential back to real space */
  fftw_execute(that->fft_phi);
  /* Interpolate grid potential (and field) back to particles */
  for (i = 0; i < nsite; i++) {
    const icart_t *index = that->index[i];
    const cart_t *wei = that->weight[i];
    int jx, jy, jz;
    const cart_t *ndw = (evec && !ik_differentiate ? that->ndw[i] : 0);
    for (jx = 0; jx < order; jx++) {
      const int px = index[jx].x;
      const double Wx = wei[jx].x;
      for (jy = 0; jy < order; jy++) {
	const int py = index[jy].y;
	const double Wy = wei[jy].y;
	const double WxWy = Wx*Wy;
	for (jz = 0; jz < order; jz++) {
	  const int pz = index[jz].z;
	  const double Wz = wei[jz].z;
	  const double w = WxWy*Wz;
	  const double p = qphigrid[px][py][pz];
	  phi[i] += w*p;
	  if (evec) {
	    if (ik_differentiate) {
	      evec[i].x += w * ex[px][py][pz];
	      evec[i].y += w * ey[px][py][pz];
	      evec[i].z += w * ez[px][py][pz];
	    } else {
	      cart_t e;
	      e.x = ndw[jx].x * Wy * Wz * p;
	      e.y = Wx * ndw[jy].y * Wz * p;
	      e.z = WxWy * ndw[jz].z * p;
	      left_operate(&e, &invlatv);
	      evec[i].x += e.x;
	      evec[i].y += e.y;
	      evec[i].z += e.z;
	    }
	  }
	}
      }
    }
  }
}

void p3m_scale(p3m_t t, double s)
{
  const int nx = t->ngrid.x, ny = t->ngrid.y, nz = t->ngrid.z;
  int ix, iy, iz;
  t->latv.xx *= s; t->latv.xy *= s; t->latv.xz *= s;
  t->latv.yx *= s; t->latv.yy *= s; t->latv.yz *= s;
  t->latv.zx *= s; t->latv.zy *= s; t->latv.zz *= s;
  t->invlatv.xx /= s; t->invlatv.xy /= s; t->invlatv.xz /= s;
  t->invlatv.yx /= s; t->invlatv.yy /= s; t->invlatv.yz /= s;
  t->invlatv.zx /= s; t->invlatv.zy /= s; t->invlatv.zz /= s;
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
      cart_t *k = t->kvec[ix][iy];
      double *g = t->green_hat[ix][iy];
      for (iz = 0; iz <= nz/2; iz++) {
	k[iz].x /= s;
	k[iz].y /= s;
	k[iz].z /= s;
	g[iz] /= s;
      }
    }
}


static double tmp_v = 0;
static double tmp_rc = 0;
static double *tmp_chi_rspace = 0;
static double *tmp_chi_p3m = 0;

static double p3m_optimize_help(const gsl_vector *x, void *that)
{
  double eta = gsl_vector_get(x, 0);
  if (eta < 0)
    return 1e8;
  *tmp_chi_rspace = 2/sqrt(tmp_rc/pow(tmp_v,1.0/3.0)) * exp(-sq(tmp_rc*eta));
  p3m_init((p3m_t) that, eta, tmp_chi_p3m);
  return sq(*tmp_chi_rspace) + sq(*tmp_chi_p3m);
}

void p3m_optimize(p3m_t that, 
		  double rc,
		  double *eta, 
		  double *chi_rspace, 
		  double *chi_p3m)
{
  const double tol = 1e-4;
  const int maxiter = 100;
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *fmin = gsl_multimin_fminimizer_alloc(T, 1);
  gsl_multimin_function f;
  gsl_vector *x = gsl_vector_alloc(1);
  gsl_vector *step = gsl_vector_alloc(1);
  int iter = 0;
  f.n = 1;
  f.f = &p3m_optimize_help;
  f.params = (void *) that;
  tmp_v = fabs(det(&that->latv));
  tmp_rc = rc;
  tmp_chi_rspace = chi_rspace;
  tmp_chi_p3m = chi_p3m;
  gsl_vector_set(x, 0, 0.1/that->grid_spacing);
  gsl_vector_set(step, 0, 0.1);
  gsl_multimin_fminimizer_set(fmin, &f, x, step);
  for (iter = 1; iter <= maxiter; iter++)
    if (gsl_multimin_fminimizer_iterate(fmin) ||
	gsl_multimin_test_size(gsl_multimin_fminimizer_size(fmin), tol) ==
	GSL_SUCCESS)
      break;
  *eta = gsl_vector_get(fmin->x,0);
  gsl_multimin_fminimizer_free(fmin);
  gsl_vector_free(x);
  gsl_vector_free(step);
}
