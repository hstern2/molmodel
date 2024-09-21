#include "pppm.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NSITE 5

static double my_random()
{
  static unsigned my_random_seed = 314159265;
  my_random_seed = (1103515245*my_random_seed + 12345);
  return ((my_random_seed/65536) % 32768)/32769.0;
}

int main()
{
  double ptv[9];              /* primitive translation vectors */
  double q[NSITE];            /* Charges for sites */
  double r[3*NSITE];          /* x,y,z coordinates for sites */
  int ngrid[3];               /* actual number of grid points that will be used */
  double rc;                  /* real-space cutoff for Ewald sum */
  double eta;                 /* Ewald screening parameter */
  double chi_rspace;          /* real-space error */
  double chi_p3m;             /* reciprocal-space error */
  double grid_spacing = 0.1;  /* suggested grid spacing */
  const int nalias = 2;       /* limit for alias sum */
  int ik_differentiate;       /* whether or not to do differentiation in k-space,
				 or real space */
  int assignment_order;       /* P3M assignment order */
  double phi[NSITE];          /* electrostatic potential at each site */
  double evec[3*NSITE];       /* electric field at each site */

  int i;

  /* Set primitive translation vectors to a simple cubic box of side 1.0 */
  ptv[0] = 1.0; ptv[1] = 0.0; ptv[2] = 0.0;
  ptv[3] = 0.0; ptv[4] = 1.0; ptv[5] = 0.0;
  ptv[6] = 0.0; ptv[7] = 0.0; ptv[8] = 1.0;

  rc = 0.4; /* set real-space cutoff */
  ik_differentiate = 1;  /* do differentiation in k space */
  assignment_order = 3;

  /* Create the object */
  p3m_t p3m = p3m_new(NSITE, ptv, grid_spacing, nalias,
		      ik_differentiate, assignment_order, ngrid);
  
  /* Find optimal value for eta, given current parameters */
  eta = 5.6; /* initial guess */
  p3m_optimize(p3m, rc, &eta, &chi_rspace, &chi_p3m);

  printf("Using %s-space differentiation\n", ik_differentiate ? "reciprocal" : "real");
  printf("Number of grid points: %d %d %d\n", ngrid[0], ngrid[1], ngrid[2]);
  printf("Assignment order: %d\n", assignment_order);
  printf("Screening parameter: %f\n", eta);
  printf("Real-space cutoff: %f\n", rc);
  printf("Estimated dimensionless real-space force error: %f\n", chi_rspace);
  printf("Estimated dimensionless P3M RMS force error: %f\n", chi_p3m);
  
  /* Assign random charges from -1 to 1 */
  printf("\n");
  for (i = 0; i < NSITE; i++) {
    q[i] = 2*my_random() - 1;
    printf("Charge of site %d: %f\n", i, q[i]);
  }
  printf("\n");
  
  /* Assign random coordinates from -1 to 1 */
  for (i = 0; i < NSITE; i++) {
    r[3*i] = 2*my_random() - 1;
    r[3*i+1] = 2*my_random() - 1;
    r[3*i+2] = 2*my_random() - 1;
    printf("Coordinates of site %d: %f %f %f\n", i, r[3*i], r[3*i+1], r[3*i+2]);
  }
  printf("\n");
  
  memset(phi, 0, NSITE*sizeof(double));
  memset(evec, 0, 3*NSITE*sizeof(double));

  /* Calculate k-space contribution to electrostatic potential and field */
  p3m_field(p3m, r, q, phi, evec);
  /* Write potential and field */
  for (i = 0; i < NSITE; i++)
    printf("Electrostatic potential at site %d: %f\n", i, phi[i]);
  printf("\n");
  for (i = 0; i < NSITE; i++)
    printf("Electric field at site %d: %f %f %f\n", i, 
            evec[3*i], evec[3*i+1], evec[3*i+2]);
  printf("\n");

  p3m_delete(p3m);

}
