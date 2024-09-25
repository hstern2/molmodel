#ifndef RMSD_H
#define RMSD_H

#include "cvec.hpp"

/***
 * Calculates the RMSD between two arrays of Cartesians.
 * Rotates and translates the second array to obtain
 * maximal coincidence with the first.
 ***/
double RootMeanSquareDistance(const CVec &c1, CVec &c2);

/***
 * Return the rotation matrix to align the second array with the first
 ***/
Tensor RotationMatrixToAlign(const CVec &c1, const Cartesian &centroid1,
			     const CVec &c2, const Cartesian &centroid2);
  
#endif /* RMSD_H */
