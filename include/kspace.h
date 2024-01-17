#ifndef KSPACE_H
#define KSPACE_H

#include "cart.h"

#ifdef __cplusplus
extern "C" {
#endif

  struct kspace;
  struct kspace *kspace_new(double cutoff, int nsite);
  void kspace_delete(struct kspace *);
  void kspace_show(const struct kspace *);
  void kspace_field(struct kspace *, const cart_t *r, const double *q, 
		    const tensor_t *latv, double volume, double eta, double *phi, cart_t *evec);
  void kspace_phi(struct kspace *, const cart_t *r, const double *q, 
		  const tensor_t *latv, double volume, double eta, double *phi);
  int kspace_number_of_wavevectors(const struct kspace *);

#ifdef __cplusplus
}
#endif

#endif /* KSPACE_H */
