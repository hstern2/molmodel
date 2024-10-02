#include "linalg.hpp"
#include "fmacro.h"
#include "fns.h"

extern "C" 
{
  void FORT(dgesv)(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  void FORT(dsysv)(const char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, 
	    double *b, int *ldb, double *work, int *lwork, int *info);
  void FORT(dgels)(const char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b,
	    int *ldb, double *work, int *lwork, int *info);
  void FORT(dgglse)(int *m, int *n, int *p, double *a, int *lda, double *b, int *ldb,
	     double *c, double *d, double *x, double *work, int *lwork, int *info);
  void FORT(dgesvd)(const char *jobu, const char *jobvt, int *m, int *n, double *a, int *lda,
	     double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
  void FORT(dsyev)(const char *jobz, const char *uplo, const int *n, double *a, const int *lda, 
		   double *w, double *work, const int *lwork, int *info);
  void FORT(dgeev)(const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda,
		   double *wr, double *wi, double *vl, const int *ldvl, double *vr, const int *ldvr,
		   double *work, const int *lwork, int *info);
}

static const char jobV[128] = "V";
static const char jobU[128] = "U";
static const char jobN[128] = "N";
static const char jobO[128] = "O";
static const char jobS[128] = "S";

void LinearSolve(DMat A, DVec x)
{
  int n = A.rows();
  int nrhs = 1;
  Vec<int> ipiv(n);
  int info;
  assert(A.columns() == n && x.size() == n);
  FORT(dgesv)(&n, &nrhs, A, &n, ipiv, x, &n, &info);
}

void LeastSquares(DMat A, DVec x)
{
  int m = A.rows();
  int n = A.columns();
  int nrhs = 1;
  assert(m == x.size());
  assert(m >= n);
  int info = 0;
  int lwork = 4*m*m;
  DVec work(lwork);
  FORT(dgels)(jobN, &m, &n, &nrhs, A, &m, x, &m, work, &lwork, &info);
}

void SingularValueDecomposition(DMat a, DVec s, DMat vt)
{
  int m = a.rows();
  int n = a.columns();
  if (m < n)
    die("lapack.C: SVD is underdetermined");
  assert(n > 0);
  int lda = (m > 1 ? m : 1);
  assert(s.size() == n);
  s.zero();
  int ldu = 1;
  int ldvt = n;
  assert(vt.rows() == n);
  assert(vt.columns() == n);
  vt.zero();
  int lwork = 6*m;
  DVec work(lwork, 0.0);
  int info = 0;
  FORT(dgesvd)(jobO, jobS, &m, &n, a, &lda, s, 0, &ldu, vt, &ldvt, work, &lwork, &info);
}

void SingularValueDecomposition(Tensor &a, Cartesian &s, Tensor &vt)
{
  int m = 3;
  int n = 3;
  int lda = 3;
  int ldu = 1;
  int ldvt = 3;
  int lwork = 256;
  DVec work(lwork, 0.0);
  int info = 0;
  s.zero();
  vt.zero();
  FORT(dgesvd)(jobO, jobS, &m, &n, a, &lda, s, 0, &ldu, vt, &ldvt, work, &lwork, &info);
}

void InvertSVD(DMat a)
{
  assert(a.is_square());
  int n = a.rows();
  DMat acopy(a.copy());
  DVec s(n);
  DMat vt(n,n);
  SingularValueDecomposition(acopy, s, vt);
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++) {
      a(i,k) = 0;
      for (int j = 0; j < n; j++)
	if (s[j] > 1e-9)
	  a(i,k) += vt(j,i) * acopy(k,j) / s[j];
    }
}

void InvertSymmetricTensor(Tensor &a)
{
  Cartesian eval;
  Tensor evec;
  Diagonalize(a,eval,evec);
  Tensor s(0);
  s.xx = fabs(eval.x) > 1e-9 ? 1/eval.x : 0;
  s.yy = fabs(eval.y) > 1e-9 ? 1/eval.y : 0;
  s.zz = fabs(eval.z) > 1e-9 ? 1/eval.z : 0;
  a = evec * s * evec.transpose();
}

void Diagonalize(const DMat a, DVec eval, DMat evec)
{
  insist(a.is_square());
  insist(evec.is_conformable_with(a));
  const int n = a.rows();
  insist(eval.size() == n);
  const int lda = (n > 1 ? n : 1);
  evec.copy(a);
  const int lwork = 3*lda*lda;
  DVec work(lwork);
  int info;
  FORT(dsyev)(jobV, jobU, &n, evec, &lda, eval, work, &lwork, &info);
}

void Diagonalize(DMat a, DVec eval_re, DVec eval_im, DMat evec)
{
  insist(a.is_square());
  insist(evec.is_conformable_with(a));
  const int n = a.rows();
  insist(eval_re.size() == n);
  insist(eval_im.size() == n);
  const int lwork = 10*n;
  DVec work(lwork);
  int info;
  FORT(dgeev)(jobN, jobV, &n, a, &n, eval_re, eval_im, 0, &n, evec, &n, work, &lwork, &info);
}

double LogDeterminant(const DMat a)
{
  const int n = a.rows();
  insist(a.columns() == n);
  DMat evec(n,n);
  DVec eval(n);
  Diagonalize(a,eval,evec);
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += log(eval[i]);
  return sum;
}


double Determinant(const DMat a)
{
  insist(a.is_square());
  const int n = a.rows();
  DVec eval(n);
  DMat evec(n,n);
  Diagonalize(a, eval, evec);
  double x = 1;
  for (int i = 0; i < n; i++)
    x *= eval[i];
  return x;
}

void Diagonalize(const Tensor &a, Cartesian &eval, Tensor &evec)
{
  int three = 3, lwork = 27;
  evec = a;
  DVec work(lwork);
  int info;
  FORT(dsyev)(jobV, jobU, &three, evec, &three, eval, work, &lwork, &info);
}

void Orthonormalize(DMat a)
{
  assert(a.is_square());
  const int n = a.rows();
  DMat vt(n,n);
  DVec s(n);
  SingularValueDecomposition(a, s, vt);
  a.apply(zero_if_almost_zero);
}

void ConstrainedLeastSquares(DMat A, DMat B, DVec c, DVec d, DVec x)
{
  int m = A.rows();
  int n = A.columns();
  int p = B.rows();
  if (m == 0 || n == 0 || p == 0)
    return;
  int lwork = m+n+p;
  DVec work(lwork);
  insist(B.columns() == n);
  insist(c.size() == m);
  insist(d.size() == p);
  insist(x.size() == n);
  insist(0 <= p && p <= n && n <= m+p);
  int info;
  FORT(dgglse)(&m, &n, &p, A, &m, B, &p, c, d, x, work, &lwork, &info);
}

void SVDLeastSquares(DMat A, DVec x, DVec c, double tol, bool verbose)
{
  int m = A.rows();
  int n = A.columns();
  insist(m >= n);
  insist(x.size() == n);
  insist(c.size() == m);
  DVec s(n);
  DMat vt(n,n);
  SingularValueDecomposition(A,s,vt);
  if (verbose)
    Out() << "SVDLeastSquares...\n"
	  << "Tolerance: " << tol << "\n";
  DVec utc = A.transpose() * c;
  for (int i = 0; i < n; i++)
    if (s[i] < tol) {
      if (verbose)
	Out() << "Zeroing singular value " << s[i] << "\n";
      utc[i] = 0;
    } else {
      if (verbose)
	Out() << "Including singular value " << s[i] << "\n";
      utc[i] /= s[i];
    }
  x.copy(vt.transpose()*utc);
  if (verbose)
    Out() << "Done.\n\n" << flush;
}
