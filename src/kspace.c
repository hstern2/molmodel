#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "kspace.h"
#include "array.h"

struct kspace
{
  int nsite, nk, nmax; /* number of sites, wavevectors, max index for integral reciprocal lattice */
  icart_t *ikvec; /* integral reciprocal lattice */
  cart_t *kvec, *ctmp, *stmp; /* wavevectors, workspace */
  double *ak, **skr, **ckr, *sQ, *cQ; /* constants, sin(k*r), cos(k*r), q sin(k*r), q cos(k*r) */
};

int kspace_number_of_wavevectors(const struct kspace *that)
{
  return that->nk;
}

static double sq(const cart_t *r) { return r->x * r->x + r->y * r->y + r->z * r->z; }
static int isq(const icart_t *i) { return i->x * i->x + i->y * i->y + i->z * i->z; }

static int icart_t_cmp(const void *i, const void *j)
{
  const int i2 = isq((const icart_t *) i);
  const int j2 = isq((const icart_t *) j);
  if (i2 < j2)
    return 1;
  else if (j2 < i2)
    return -1;
  else
    return 0;
}

struct kspace *kspace_new(double cutoff, int nsite)
{
  int nk = 0, nx, ny, nz;
  struct kspace *that = (struct kspace *) malloc(sizeof(struct kspace));
  const double nsqmax = cutoff*cutoff + 1e-8;
  const int nmax = (int) floor(cutoff);
  for (nx = 0; nx <= nmax; nx++)
    for (ny = (nx == 0) ? 0 : -nmax; ny <= nmax; ny++)
      for (nz = (nx == 0 && ny == 0) ? 1 : -nmax; nz <= nmax; nz++)
	if (nx*nx + ny*ny + nz*nz <= nsqmax)
	  nk++;
  that->ikvec = (icart_t *) malloc(nk * sizeof(icart_t));
  that->nk = nk;
  nk = 0;
  for (nx = 0; nx <= nmax; nx++)
    for (ny = (nx == 0) ? 0 : -nmax; ny <= nmax; ny++)
      for (nz = (nx == 0 && ny == 0) ? 1 : -nmax; nz <= nmax; nz++)
	if (nx*nx + ny*ny + nz*nz <= nsqmax) {
	  that->ikvec[nk].x = nx;
	  that->ikvec[nk].y = ny;
	  that->ikvec[nk].z = nz;
	  nk++;
	}
  qsort(that->ikvec, nk, sizeof(icart_t), icart_t_cmp);
  that->kvec = (cart_t *) malloc(nk * sizeof(cart_t));
  that->ak = (double *) malloc(nk * sizeof(double));
  that->skr = double_array2D_new(nsite, nk);
  that->ckr = double_array2D_new(nsite, nk);
  that->ctmp = ((cart_t *) malloc((2*nmax + 1) * sizeof(cart_t))) + nmax;
  that->stmp = ((cart_t *) malloc((2*nmax + 1) * sizeof(cart_t))) + nmax;
  that->sQ = (double *) malloc(nk * sizeof(double));
  that->cQ = (double *) malloc(nk * sizeof(double));
  that->nmax = nmax;
  that->nk = nk;
  that->nsite = nsite;
  return that;
}

void kspace_delete(struct kspace *that)
{
  free(that->ikvec);
  free(that->kvec);
  free(that->ak);
  double_array2D_delete(that->skr);
  double_array2D_delete(that->ckr);
  free(that->ctmp - that->nmax);
  free(that->stmp - that->nmax);
  free(that->sQ);
  free(that->cQ);
  free(that);
}

void kspace_field(struct kspace *that, const cart_t *r, const double *q, 
		  const tensor_t *latv, double volume, double eta,
		  double *phi, cart_t *evec)
{
  const icart_t *ikvec = that->ikvec;
  cart_t *kvec = that->kvec, *ctmp = that->ctmp, *stmp = that->stmp;
  double *ak = that->ak, **ckr = that->ckr, **skr = that->skr, *sQ = that->sQ, *cQ = that->cQ;
  int i, k, nk = that->nk, nmax = that->nmax, nsite = that->nsite;
  const double c = 0.25/(eta*eta), fourpi_V = 4*M_PI/volume;
  for (k = 0; k < nk; k++) {
    double tmp;
    kvec[k].x = latv->xx*ikvec[k].x + latv->xy*ikvec[k].y + latv->xz*ikvec[k].z;
    kvec[k].y = latv->yx*ikvec[k].x + latv->yy*ikvec[k].y + latv->yz*ikvec[k].z;
    kvec[k].z = latv->zx*ikvec[k].x + latv->zy*ikvec[k].y + latv->zz*ikvec[k].z;
    tmp = sq(&kvec[k]);
    ak[k] = fourpi_V*exp(-tmp*c)/tmp;
  }
  for (i = 0; i < nsite; i++) {
    cart_t tmp, ck, sk;
    tmp.x = latv->xx*r[i].x + latv->yx*r[i].y + latv->zx*r[i].z;
    tmp.y = latv->xy*r[i].x + latv->yy*r[i].y + latv->zy*r[i].z;
    tmp.z = latv->xz*r[i].x + latv->yz*r[i].y + latv->zz*r[i].z;
    ck.x = cos(tmp.x);
    ck.y = cos(tmp.y);
    ck.z = cos(tmp.z);
    sk.x = sin(tmp.x);
    sk.y = sin(tmp.y);
    sk.z = sin(tmp.z);
    ctmp[0].x = ctmp[0].y = ctmp[0].z = 1;
    stmp[0].x = stmp[0].y = stmp[0].z = 0;
    for (k = 1; k <= nmax; k++) {
      ctmp[k].x = ck.x*ctmp[k-1].x - sk.x*stmp[k-1].x;
      stmp[k].x = sk.x*ctmp[k-1].x + ck.x*stmp[k-1].x;
      ctmp[-k].x = ctmp[k].x;
      stmp[-k].x = -stmp[k].x;
      ctmp[k].y = ck.y*ctmp[k-1].y - sk.y*stmp[k-1].y;
      stmp[k].y = sk.y*ctmp[k-1].y + ck.y*stmp[k-1].y;
      ctmp[-k].y = ctmp[k].y;
      stmp[-k].y = -stmp[k].y;
      ctmp[k].z = ck.z*ctmp[k-1].z - sk.z*stmp[k-1].z;
      stmp[k].z = sk.z*ctmp[k-1].z + ck.z*stmp[k-1].z;
      ctmp[-k].z = ctmp[k].z;
      stmp[-k].z = -stmp[k].z;
    }
    for (k = 0; k < nk; k++) {
      const int ix = ikvec[k].x;
      const int iy = ikvec[k].y;
      const int iz = ikvec[k].z;
      const double cxi = ctmp[ix].x;
      const double sxi = stmp[ix].x;
      const double cyi = ctmp[iy].y;
      const double syi = stmp[iy].y;
      const double czi = ctmp[iz].z;
      const double szi = stmp[iz].z;
      const double ccss = cxi*cyi - sxi*syi;
      const double sccs = sxi*cyi + cxi*syi;
      ckr[i][k] = ccss*czi - sccs*szi;
      skr[i][k] = sccs*czi + ccss*szi;
    }
  }
  for (k = 0; k < nk; k++)
    sQ[k] = cQ[k] = 0;
  for (i = 0; i < nsite; i++) {
    const double qi = q[i], *skri = skr[i], *ckri = ckr[i];
    for (k = 0; k < nk; k++) {
      sQ[k] += qi*skri[k];
      cQ[k] += qi*ckri[k];
    }
  }
  for (i = 0; i < nsite; i++) {
    double gx = 0, gy = 0, gz = 0, p = 0;
    const double *skri = skr[i], *ckri = ckr[i];
    for (k = 0; k < nk; k++) {
      const double two_a = 2*ak[k];
      const double s = skri[k];
      const double c = ckri[k];
      const double tmp = two_a*(cQ[k]*s - sQ[k]*c);
      gx += tmp*kvec[k].x;
      gy += tmp*kvec[k].y;
      gz += tmp*kvec[k].z;
      p += two_a*(sQ[k]*s + cQ[k]*c);
    }
    phi[i] += p;
    evec[i].x += gx;
    evec[i].y += gy;
    evec[i].z += gz;
  }
}

void kspace_phi(struct kspace *that, const cart_t *r, const double *q, 
		const tensor_t *latv, double volume, double eta, double *phi)
{
  const icart_t *ikvec = that->ikvec;
  cart_t *kvec = that->kvec, *ctmp = that->ctmp, *stmp = that->stmp;
  double *ak = that->ak, **ckr = that->ckr, **skr = that->skr, *sQ = that->sQ, *cQ = that->cQ;
  int i, k, nk = that->nk, nmax = that->nmax, nsite = that->nsite;
  const double c = 0.25/(eta*eta), fourpi_V = 4*M_PI/volume;
  for (k = 0; k < nk; k++) {
    double tmp;
    kvec[k].x = latv->xx*ikvec[k].x + latv->xy*ikvec[k].y + latv->xz*ikvec[k].z;
    kvec[k].y = latv->yx*ikvec[k].x + latv->yy*ikvec[k].y + latv->yz*ikvec[k].z;
    kvec[k].z = latv->zx*ikvec[k].x + latv->zy*ikvec[k].y + latv->zz*ikvec[k].z;
    tmp = sq(&kvec[k]);
    ak[k] = fourpi_V*exp(-tmp*c)/tmp;
  }
  for (i = 0; i < nsite; i++) {
    cart_t tmp, ck, sk;
    tmp.x = latv->xx*r[i].x + latv->yx*r[i].y + latv->zx*r[i].z;
    tmp.y = latv->xy*r[i].x + latv->yy*r[i].y + latv->zy*r[i].z;
    tmp.z = latv->xz*r[i].x + latv->yz*r[i].y + latv->zz*r[i].z;
    ck.x = cos(tmp.x);
    ck.y = cos(tmp.y);
    ck.z = cos(tmp.z);
    sk.x = sin(tmp.x);
    sk.y = sin(tmp.y);
    sk.z = sin(tmp.z);
    ctmp[0].x = ctmp[0].y = ctmp[0].z = 1;
    stmp[0].x = stmp[0].y = stmp[0].z = 0;
    for (k = 1; k <= nmax; k++) {
      ctmp[k].x = ck.x*ctmp[k-1].x - sk.x*stmp[k-1].x;
      stmp[k].x = sk.x*ctmp[k-1].x + ck.x*stmp[k-1].x;
      ctmp[-k].x = ctmp[k].x;
      stmp[-k].x = -stmp[k].x;
      ctmp[k].y = ck.y*ctmp[k-1].y - sk.y*stmp[k-1].y;
      stmp[k].y = sk.y*ctmp[k-1].y + ck.y*stmp[k-1].y;
      ctmp[-k].y = ctmp[k].y;
      stmp[-k].y = -stmp[k].y;
      ctmp[k].z = ck.z*ctmp[k-1].z - sk.z*stmp[k-1].z;
      stmp[k].z = sk.z*ctmp[k-1].z + ck.z*stmp[k-1].z;
      ctmp[-k].z = ctmp[k].z;
      stmp[-k].z = -stmp[k].z;
    }
    for (k = 0; k < nk; k++) {
      const int ix = ikvec[k].x;
      const int iy = ikvec[k].y;
      const int iz = ikvec[k].z;
      const double cxi = ctmp[ix].x;
      const double sxi = stmp[ix].x;
      const double cyi = ctmp[iy].y;
      const double syi = stmp[iy].y;
      const double czi = ctmp[iz].z;
      const double szi = stmp[iz].z;
      const double ccss = cxi*cyi - sxi*syi;
      const double sccs = sxi*cyi + cxi*syi;
      ckr[i][k] = ccss*czi - sccs*szi;
      skr[i][k] = sccs*czi + ccss*szi;
    }
  }
  for (k = 0; k < nk; k++)
    sQ[k] = cQ[k] = 0;
  for (i = 0; i < nsite; i++) {
    const double qi = q[i], *skri = skr[i], *ckri = ckr[i];
    for (k = 0; k < nk; k++) {
      sQ[k] += qi*skri[k];
      cQ[k] += qi*ckri[k];
    }
  }
  for (i = 0; i < nsite; i++) {
    double p = 0;
    const double *skri = skr[i], *ckri = ckr[i];
    for (k = 0; k < nk; k++)
      p += 2*ak[k]*(sQ[k]*skri[k] + cQ[k]*ckri[k]);
    phi[i] += p;
  }
}

static void icart_show(const char *name, int n, const icart_t *i)
{
  int k;
  printf("%s:\n", name);
  for (k = 0; k < n; k++)
    printf("%d %d %d\n", i[k].x, i[k].y, i[k].z);
  printf("\n");
}

static void cart_show(const char *name, int n, const cart_t *r)
{
  int k;
  printf("%s:\n", name);
  for (k = 0; k < n; k++)
    printf("%18.12g %18.12g %18.12g\n", r[k].x, r[k].y, r[k].z);
  printf("\n");
}

static void double_show(const char *name, int n, const double *d)
{
  int k;
  printf("%s:\n", name);
  for (k = 0; k < n; k++)
    printf("%18.12g\n", d[k]);
  printf("\n");
}

void kspace_show(const struct kspace *t)
{
  printf("nsite: %d  nk: %d  nmax: %d\n", t->nsite, t->nk, t->nmax);
  icart_show("ikvec", t->nk, t->ikvec);
  cart_show("kvec", t->nk, t->kvec);
  cart_show("ctmp", 2*t->nmax+1, t->ctmp - t->nmax);
  cart_show("stmp", 2*t->nmax+1, t->stmp - t->nmax);
  double_show("ak", t->nk, t->ak);
  double_show("sQ", t->nk, t->sQ);
  double_show("cQ", t->nk, t->cQ);

}
