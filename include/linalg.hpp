#ifndef LINALG_H
#define LINALG_H

#include "dmat.hpp"
#include "dvec.hpp"
#include "coord.hpp"

/**
 * solve Ax = b
 * A is matrix on input, factorized matrix on output
 * x is right hand side (b) on input, solution on output
 **/
void LinearSolve(DMat A, DVec x);

/**
 * find x to minimize |Ax-b|^2 for an overdetermined system
 * A is an m by n matrix on input, factorized matrix on output
 * x is a vector of length m on input containing right hand side,
 * on output, the first n entries of x are the solution.
 **/
void LeastSquares(DMat A, DVec x);


/***
 * find x to minimize |Ax-c|^2 subject to constraint Bx = d
 * A is an m by n matrix
 * B is a p by n  matrix
 * c is a vector of length m
 * d is a vector of length p
 * x is a vector of length n
 ***/
void ConstrainedLeastSquares(DMat A, DMat B, DVec c, DVec d, DVec x);

/***
 * find x to minimize |Ax-c|^2
 * A is an m by n matrix
 * c is a vector of length m
 * x is a vector of length n
 * zero out singular values if they are smaller than tol
 ***/
void SVDLeastSquares(DMat A, DVec x, DVec c, double tol, bool verbose = false);

/***
 * Computes singular value decomposition of A, and stores 
 * the factorization in A ( = u), s, and vt.  Destroys A!!
 * A(m,n) with m >= n;  s[n];  vt[n]
 ***/
void SingularValueDecomposition(DMat a, DVec s, DMat vt);
void SingularValueDecomposition(Tensor &a, Cartesian &s, Tensor &vt);
void InvertSVD(DMat a);
void InvertSymmetricTensor(Tensor &a);

/***
 * Computes eigenvalues and eigenvectors of a symmetric matrix A
 * A = evec . eval . evec^t
 * Eigenvectors are the columns of evec
 ***/
void Diagonalize(const DMat a, DVec eval, DMat evec);
void Diagonalize(const Tensor &a, Cartesian &eval, Tensor &evec);

/*** Computes eigenvalues and eigenvectors of a general matrix A
 * A = evec . eval . evec^-1
 ***/
void Diagonalize(DMat a, DVec eval_re, DVec eval_im, DMat evec);
  
/*** 
 * Computes determinant of a symmetric matrix A
 ***/
double Determinant(const DMat a);

/***
 * Computes log of the determinant of a positive-definite matrix A
 ***/
double LogDeterminant(const DMat a);

/***
 * Orthonormalize columns via singular value decomposition
 ***/
void Orthonormalize(DMat);

#endif /* LINALG_H */
