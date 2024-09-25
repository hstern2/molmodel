#ifndef LAPACK_H
#define LAPACK_H

#include "rmat.hpp"
#include "coord.hpp"

/***
 * Solve A x = b
 */
void LinearSolve(const RMat &a, RVec &x, const RVec &b);

/***
 * Solve A x = b where A is symmetric
 ***/
void SymmetricLinearSolve(const RMat &a, RVec &x, const RVec &b);
void SymmetricLinearSolve(const Tensor &a, Cartesian &x, const Cartesian &b);

/***
 * Multiply a times b 
 ***/
RMat operator*(const RMat &a, const RMat &b);

/*** 
 * Minimize |A x - b|^2  
 * Destroys A
 ***/
void LeastSquaresFit(RMat &a, RVec &x, RVec &b);

/***
 * Minimize |A x - c|^2 subject to the constraint B x = d 
 * Destroys A, B, c, d.
 ***/
void ConstrainedLeastSquaresFit(RMat &a, RVec &x, RVec &c,
				RMat &b, RVec &d, double *rms = 0);

/***
 * Computes singular value decomposition of A, and stores 
 * the factorization in A ( = u), s, and vt.  Destroys A!!
 * A(m,n) with m >= n;  s[n];  vt[n]
 ***/
void SingularValueDecomposition(RMat &a, RVec &s, RMat &vt);
void SingularValueDecomposition(Tensor &a, Cartesian &s, Tensor &vt);
void InvertSVD(RMat &a);
void InvertSymmetricTensor(Tensor &a);

/***
 * Computes eigenvalues and eigenvectors of a symmetric matrix A
 ***/
void Diagonalize(const RMat &a, RVec &eval, RMat &evec);
void Diagonalize(const Tensor &a, Cartesian &eval, Tensor &evec);

/***
 * Orthonormalize columns via singular value decomposition
 ***/
void Orthonormalize(RMat &);

#endif /* LAPACK_H */
