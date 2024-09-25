#include "rmsd.hpp"
#include "cvec.hpp"
#include "lapack.hpp"
#include "fns.h"

/***
 * Calculates the RMSD between two arrays of Cartesians.
 * Rotates and translates the second array to obtain
 * maximal coincidence with the first.
 ***/
double RootMeanSquareDistance(const CVec &c1, CVec &c2)
{
  const int n = c1.size();
  insist(c2.size() == n);
  if (n <= 1)
    return 0;
  /* Get centroids */
  const Cartesian centroid1 = c1.average();
  const Cartesian centroid2 = c2.average();
  const Tensor rotation_matrix = RotationMatrixToAlign(c1,centroid1,c2,centroid2);
  /* Rotate and translate */
  int i;
  for (i = 0; i < n; i++)
    c2[i] = rotation_matrix*(c2[i]-centroid2) + centroid1;
  /* Calculate and return RMSD */
  double rmsd = 0;
  for (i = 0; i < n; i++)
    rmsd += (c1[i]-c2[i]).sq();
  return sqrt(rmsd/n);
}

Tensor RotationMatrixToAlign(const CVec &c1, const Cartesian &centroid1,
			     const CVec &c2, const Cartesian &centroid2)
{
  /* Calculate rotation matrix */
  const int n = c1.size();
  insist(c2.size() == n);
  if (n <= 1)
    return Tensor(1);
  Tensor u(0), vt;
  Cartesian s;
  for (int i = 0; i < n; i++)
    u += Tensor(c1[i] - centroid1, c2[i] - centroid2); /* outer product */
  SingularValueDecomposition(u, s, vt);
  Tensor rotation_matrix = u*vt;
  if (rotation_matrix.determinant() < 0) {
    if (s.x < s.y && s.x < s.z)
      vt.set_row1(-vt.row(0));
    else if (s.y < s.x && s.y < s.z)
      vt.set_row2(-vt.row(1));
    else
      vt.set_row3(-vt.row(2));
    rotation_matrix = u*vt;
  }
  if (fabs(rotation_matrix.determinant() - 1.0) > 1e-8)
    die("RotationMatrixToAlign: rotation matrix determinant is %.15f",
	rotation_matrix.determinant());
  return rotation_matrix;
}
