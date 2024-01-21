#ifndef PPPM_H
#define PPPM_H

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct p3m *p3m_t;
  
  /* p3m_new:
     create new P3M object */
  p3m_t p3m_new(int nsite, /* number of sites */
		const double primitive_translation_vectors[9],
		double grid_spacing, /* suggested grid spacing */
		int nalias, /* limit of alias sum.  suggested value: 2 */
		int ik_differentiate, /* if nonzero, differentiate by multiplying by ik_n  
					 This requires three extra FFTs.
					 Otherwise, differentiate in real space,
					 using gradient of assignment function */
		int assignment_order, /* must be an integer >= 1 (and >= 2
				         if real-space differentiation is used) */
		int ngrid[3]);         /* output -- actual number of grid points 
					  in each dimension */
  
  /* p3m_copy:
     Return new P3M object with same parameters, 
     but different number of sites */
  p3m_t p3m_copy(const p3m_t, int nsite);
  
  /* p3m_init:
     Initialize (calculate Green function coefficients)
     eta is the Ewald screening parameter 
     output: chi, the dimensionless RMS error in the forces,
     chi^2 = V^{-2/3} \int_V\int_V |E_P3M(r1,r2) - E_Ewald(r1,r2)| ^2 dr1 dr2 */
  void p3m_init(p3m_t, double eta, double *chi);
  
  /* p3m_optimize:
     Given a real-space cutoff (rc),
     find the optimal value for eta 
     Output: eta and dimensionless RMS errors 
     for real-space and P3M */
  void p3m_optimize(p3m_t, 
		    double rc, 
		    double *eta, 
		    double *chi_rspace, 
		    double *chi_p3m);
  
  /* p3m_field:
     Compute potential and electric field and add them to phi and evec.
     Note, phi and evec are added to, not replaced. 
     If evec is null, only the potential is computed. */
  void p3m_field(p3m_t,
		 const double *r, /* site coordinates (3 nsite) */
		 const double *q, /* site charges (nsite) */
		 double *phi, /* potential (nsite) */
		 double *evec); /* electric field (3 nsite) */
  
  /* p3m_scale:
     Scale boundary conditions by a factor of s */
  void p3m_scale(p3m_t, double s);
  
  /* p3m_delete:
     Free memory */
  void p3m_delete(p3m_t);
  
#ifdef __cplusplus
}
#endif

#endif /* PPPM_H */
