#include "lapack.hpp"
#include "fns.h"
#include "fmacro.h"
#include "out.hpp"

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
  void FORT(dsyev)(const char *jobz, const char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
}

static const char jobV[128] = "V";
static const char jobU[128] = "U";
static const char jobN[128] = "N";
static const char jobO[128] = "O";
static const char jobS[128] = "S";

void LinearSolve(const RMat &a, RVec &x, const RVec &b)
{
  RMat aa(a.copy());
  int n = aa.rows();
  int nrhs = 1;
  Vec<int> ipiv(n);
  int info;
  x = b.copy();
  assert(aa.columns() == n && x.size() == n);
  FORT(dgesv)(&n, &nrhs, aa, &n, ipiv, x, &n, &info);
  assert(info == 0);
}

void SymmetricLinearSolve(const RMat &a, RVec &x, const RVec &b)
{
  RMat aa(a.copy());
  int n = aa.rows();
  int nrhs = 1;
  Vec<int> ipiv(n);
  int lwork = n*n;
  RVec work(lwork);
  int info;
  x = b.copy();
  assert(aa.columns() == n && x.size() == n);
  FORT(dsysv)(jobU, &n, &nrhs, aa, &n, ipiv, x, &n, work, &lwork, &info);
  assert(info == 0);
}

void SymmetricLinearSolve(const Tensor &a, Cartesian &x, const Cartesian &b)
{
  Tensor aa(a);
  int n = 3;
  int nrhs = 1;
  int ipiv[3];
  int lwork = 9;
  double work[9];
  int info;
  x = b;
  FORT(dsysv)(jobU, &n, &nrhs, aa, &n, ipiv, x, &n, work, &lwork, &info);
  assert(info == 0);
}

void LeastSquaresFit(RMat &a, RVec &x, RVec &b)
{
  int m = a.rows();
  int n = a.columns();
  int nrhs = 1;
  int lda = (m > 1 ? m : 1);
  int ldb = (lda > n ? lda : n);
  int lwork = 2*ldb*ldb;
  RVec work(lwork);
  int info;
  assert(m >= 0);
  assert(n >= 0);
  assert(b.size() == m);
  x.resize(ldb);
  for (int i = 0; i < m; i++)
    x[i] = b[i];
  FORT(dgels)(jobN, &m, &n, &nrhs, a, &lda, x, &ldb, work, &lwork, &info);
  assert(info == 0);
  x.resize(n);
}

void ConstrainedLeastSquaresFit(RMat &a, RVec &x, RVec &c,
				RMat &b, RVec &d, double *rms)
{
  int m = a.rows();
  int n = a.columns();
  int p = b.rows();
  int lda = (m > 1 ? m : 1);
  int ldb = (p > 1 ? p : 1);
  int lwork = 10*(m+n+p);
  RVec work(lwork);
  int info;
  x.resize(n);
  assert(m >= 0);
  assert(n >= 0);
  assert(b.columns() == n);
  assert(c.size() == m);
  assert(d.size() == p);
  assert(0 <= p && p <= n && n <= m+p);
  FORT(dgglse)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
  assert(info == 0);
  if (rms) {
    *rms = 0;
    for (int i = n - p; i < m; i++)
      *rms += sq(c[i]);
    *rms = sqrt(*rms/m);
  }
}

void SingularValueDecomposition(RMat &a, RVec &s, RMat &vt)
{
  int m = a.rows();
  int n = a.columns();
  if (m < n)
    die("lapack.C: SVD is underdetermined");
  assert(n > 0);
  int lda = (m > 1 ? m : 1);
  s.resize(n);
  s.zero();
  int ldu = 1;
  int ldvt = n;
  vt.resize(n,n);
  vt.zero();
  int lwork = 6*m;
  RVec work(lwork, 0.0);
  int info = 0;
  FORT(dgesvd)(jobO, jobS, &m, &n, a, &lda, s, 0, &ldu, vt, &ldvt, work, &lwork, &info);
  assert(info == 0);
}

void SingularValueDecomposition(Tensor &a, Cartesian &s, Tensor &vt)
{
  int m = 3;
  int n = 3;
  int lda = 3;
  int ldu = 1;
  int ldvt = 3;
  int lwork = 256;
  RVec work(lwork, 0.0);
  int info = 0;
  s.zero();
  vt.zero();
  FORT(dgesvd)(jobO, jobS, &m, &n, a, &lda, s, 0, &ldu, vt, &ldvt, work, &lwork, &info);
  assert(info == 0);
}

void InvertSVD(RMat &a)
{
  assert(a.is_square());
  int n = a.rows();
  RMat acopy(a.copy());
  RVec s;
  RMat vt;
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

void Diagonalize(const RMat &a, RVec &eval, RMat &evec)
{
  assert(a.is_square());
  evec = a.copy();
  int n = a.rows();
  int lda = (n > 1 ? n : 1);
  eval.resize(n);
  int lwork(3*lda*lda);
  RVec work(lwork);
  int info;
  FORT(dsyev)(jobV, jobU, &n, evec, &lda, eval, work, &lwork, &info);
  assert(info == 0);
}

void Diagonalize(const Tensor &a, Cartesian &eval, Tensor &evec)
{
  int three = 3, lwork = 27;
  evec = a;
  RVec work(lwork);
  int info;
  FORT(dsyev)(jobV, jobU, &three, evec, &three, eval, work, &lwork, &info);
  assert(info == 0);
}

void Orthonormalize(RMat &a)
{
  RMat vt;
  RVec s;
  SingularValueDecomposition(a, s, vt);
  a.apply(zero_if_almost_zero);
}
