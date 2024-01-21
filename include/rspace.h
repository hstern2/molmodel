/* $Id: rspace.h,v 1.17 2010/09/13 18:13:44 hstern Exp $ */

#ifndef RSPACE_H
#define RSPACE_H

#include "cart.h"
#include "bcondtype.h"

#define SITE_HAS_Q 1
#define SITE_HAS_LJ 2

#ifdef __cplusplus
extern "C" {
#endif

  struct Chebyshev;
  struct PairBinary;
  
  struct Chebyshev *RSpaceChebyshev_new(double rhi, double eta, double *maxerr, double *xerr);

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
		   double lj_1_4_scale_factor);
  
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
			 double lj_1_4_scale_factor);
  
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
		      double elec_1_4_scale_factor);

  void RSpaceEwaldPhiAtGridpoints(int nsite, int ngrid, const cart_t *rsite,
				  const cart_t *rgrid, enum bcondtype_t bctype,
				  const tensor_t *latv,
				  double rlo2, double rhi2,
				  const int *site_mask,
				  const double *q,
				  double *phi,
				  double eta, const struct Chebyshev *rspace_chebyshev,
				  const double *rscreen);
  
  void RSpaceNonPeriodicPhi(int nsite, int ngroup, const int (*group)[2], const cart_t *r,
			    const struct PairBinary *group_excluded,
			    const struct PairBinary *excluded,
			    const int *site_mask, const double *q, 
			    double *phi, const double *rscreen,
			    const struct PairBinary *one_four,
			    double elec_1_4_scale_factor);
  
  void RSpaceNonPeriodicPhiAtGridpoints(int nsite, int ngrid, const cart_t *rsite,
					const cart_t *rgrid,
					const int *site_mask,
					const double *q,
					double *phi, const double *rscreen);
  
#ifdef __cplusplus
}
#endif

#endif /* RSPACE_H */
