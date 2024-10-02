#ifndef DMAT_H
#define DMAT_H

#include <cmath>
#include "dvec.hpp"
#include "array.h"

class DMat
{
  struct MatRep
  {
    int ref, n, m;
    double **p;
    MatRep(int nn, int mm) : ref(1), n(nn), m(mm)
    { p = double_array2D_new(m,n); }
    ~MatRep() 
    { double_array2D_delete(p); }
    double & operator()(int i, int j)
    {
      assert(0 <= i && i < n && 0 <= j && j < m);
      return p[j][i];
    }
  };

  void destroy()
  {
    if (--rep->ref == 0) {
      delete rep;
      rep = 0;
    }
  }

public:

  MatRep *rep;

  int rows() const
  { return rep->n; }

  int columns() const
  { return rep->m; }

  double & operator()(int i, int j)
  { return (*rep)(i,j); }

  const double & operator()(int i, int j) const
  { return (*rep)(i,j); }

  operator double*() { return rep->p ? *rep->p : 0; }

  operator const double*() const { return rep->p ? *rep->p : 0; }

  DMat(int n = 0, int m = 0) : rep(new MatRep(n,m)) { }

  DMat(const DMat &mat) : rep(mat.rep) { rep->ref++; }

  DMat(int n, int m, const double &t) : rep(new MatRep(n,m))
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) = t;
  }

  ~DMat()
  { destroy(); }

  DMat & operator=(const DMat &mat)
  { 
    mat.rep->ref++;
    destroy();
    rep = mat.rep;
    return *this;
  }

  void resize(int n, int m)
  {
    if (n == rows() && m == columns())
      return;
    DMat mat(n, m);
    for (int i = 0; i < rows() && i < mat.rows(); i++)
      for (int j = 0; j < columns() && j < mat.columns(); j++)
	mat(i,j) = (*this)(i,j);
    *this = mat;
  }

  DMat transpose() const
  {
    DMat mat(columns(), rows());
    for (int i = 0; i < rows(); i++) {
      for (int j = 0; j < columns(); j++)
	mat(j,i) = (*this)(i,j);
    }
    return mat;
  }

  friend DMat HorizontalPaste(const DMat m1, const DMat m2)
  {
    assert(m1.rows() == m2.rows());
    DMat mat(m1.rows(), m1.columns() + m2.columns());
    for (int i = 0; i < m1.rows(); i++) {
      int j;
      for (j = 0; j < m1.columns(); j++)
	mat(i,j) = m1(i,j);
      for (j = 0; j < m2.columns(); j++)
	mat(i,j+m1.columns()) = m2(i,j);
    }
    return mat;
  }

  friend DMat VerticalPaste(const DMat m1, const DMat m2)
  {
    assert(m1.columns() == m2.columns());
    DMat mat(m1.rows() + m2.rows(), m1.columns());
    int i;
    for (i = 0; i < m1.rows(); i++) {
      for (int j = 0; j < m1.columns(); j++)
	mat(i,j) = m1(i,j);
    }
    for (i = 0; i < m2.rows(); i++) {
      for (int j = 0; j < m1.columns(); j++)
	mat(i+m1.rows(), j) = m2(i,j);
    }
    return mat;
  }

  DMat & copy(const DMat &mat)
  {
    assert(is_conformable_with(mat));
    memcpy((double *) (*this), (const double *) mat, rows()*columns()*sizeof(double));
    return *this;
  }

  DMat copy() const
  {
    DMat r(rows(), columns());
    r.copy(*this);
    return r;
  }

  int is_square() const { return rows() == columns(); }

  int is_conformable_with(const DMat &mat) const
  { return rows() == mat.rows() && columns() == mat.columns(); }

  int is_diagonally_dominant() const
  {
    for (int i = 0; i < rows(); i++) {
      double sum = 0;
      for (int j = 0; j < columns(); j++)
	if (j != i)
	  sum += fabs((*this)(i,j));
      if (sum >= (*this)(i,i))
	return 0;
    }
    return 1;
  }

  double minimum() const
  {
    if (rows() == 0 || columns() == 0)
      return 0;
    double min = (*this)(0,0);
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	if ((*this)(i,j) < min)
	  min = (*this)(i,j);
    return min;
  }

  double maximum() const
  {
    if (rows() == 0 || columns() == 0)
      return 0;
    double max = (*this)(0,0);
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	if ((*this)(i,j) > max)
	  max = (*this)(i,j);
    return max;
  }

  double minimum_in_row(int i) const
  {
    assert(0 <= i && i < rows());
    if (columns() == 0)
      return 0;
    double min = (*this)(i,0);
    for (int j = 1; j < columns(); j++)
      if ((*this)(i,j) < min)
	min = (*this)(i,j);
    return min;
  }

  double maximum_in_row(int i) const
  {
    assert(0 <= i && i < rows());
    if (columns() == 0)
      return 0;
    double max = (*this)(i,0);
    for (int j = 1; j < columns(); j++)
      if ((*this)(i,j) > max)
	max = (*this)(i,j);
    return max;
  }

  double minimum_in_column(int j) const
  {
    assert(0 <= j && j < columns());
    if (rows() == 0)
      return 0;
    double min = (*this)(0,j);
    for (int i = 1; i < rows(); i++)
      if ((*this)(i,j) < min)
	min = (*this)(i,j);
    return min;
  }

  double maximum_in_column(int j) const
  {
    assert(0 <= j && j < columns());
    if (rows() == 0)
      return 0;
    double max = (*this)(0,j);
    for (int i = 1; i < rows(); i++)
      if ((*this)(i,j) > max)
	max = (*this)(i,j);
    return max;
  }

  DMat & apply(double (*f)(double))
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) = f((*this)(i,j));
    return *this;
  }

  DMat & operator+=(const DMat &t)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) += t(i,j);
    return *this;
  }

  DMat & operator-=(const DMat &t)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) -= t(i,j);
    return *this;
  }

  DMat & operator*=(double r)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) *= r;
    return *this;
  }

  DMat operator-() const
  { 
    DMat m(rows(), columns());
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	m(i,j) = -(*this)(i,j);
    return m;
  }

  DVec row(int i)
  {
    DVec r(columns());
    for (int j = 0; j < columns(); j++)
      r[j] = (*this)(i,j);
    return r;
  }

  DVec column(int j)
  {
    DVec c(rows());
    for (int i = 0; i < rows(); i++)
      c[i] = (*this)(i,j);
    return c;
  }

  void zero()
  { memset((double *) (*this), 0, rows()*columns()*sizeof(double)); }

  static DMat diagonal(const DVec v)
  {
    const int n = v.size();
    DMat m(n,n,0);
    for (int i = 0; i < n; i++)
      m(i,i) = v[i];
    return m;
  }
  
  friend ostream & operator<<(ostream &s, const DMat &m)
  {
    s << m.rows() << " " << m.columns() << "\n";
    for (int i = 0; i < m.rows(); i++) {
      for (int j = 0; j < m.columns(); j++)
	s << m(i,j) << " ";
      s << "\n";
    }
    s << "\n";
    return s;
  }

  friend istream & operator>>(istream &s, DMat &mat)
  {
    int n, m;
    s >> n >> m;
    if (n != mat.rows() || m != mat.columns())
      mat = DMat(n,m);
    for (int i = 0; i < mat.rows(); i++)
      for (int j = 0; j < mat.columns(); j++)
	s >> mat(i,j);
    return s;
  }
  
  void show_self() { Out() << *this; }
  
  friend DMat operator+(const DMat &a, const DMat &b)
  {
    DMat c = a.copy();
    return c += b;
  }
  
  friend DMat operator-(const DMat &a, const DMat &b)
  {
    DMat c = a.copy();
    return c -= b;
  }
  
  friend DMat operator*(double a, const DMat &b)
  {
    DMat c = b.copy();
    return c *= a;
  }
  
  friend DMat operator*(const DMat &a, double b)
  {
    DMat c = a.copy();
    return c *= b;
  }
  
  friend DMat operator*(const DMat &a, const DMat &b)
  {
    assert(a.columns() == b.rows());
    DMat c(a.rows(), b.columns());
    for (int i = 0; i < a.rows(); i++)
      for (int j = 0; j < b.columns(); j++) {
	double x = 0;
	for (int k = 0; k < a.columns(); k++)
	  x += a(i,k)*b(k,j);
	c(i,j) = x;
      }
    return c;
  }
  
  friend DVec operator*(const DVec &a, const DMat &b)
  {
    assert(a.size() == b.rows());
    DVec c(b.columns());
    for (int j = 0; j < b.columns(); j++) {
      double x = 0;
      for (int k = 0; k < a.size(); k++)
	x += a[k]*b(k,j);
      c[j] = x;
    }
    return c;
  }
  
  friend DVec operator*(const DMat &a, const DVec &b)
  {
    assert(a.columns() == b.size());
    DVec c(a.rows());
    for (int i = 0; i < a.rows(); i++) {
      double x = 0;
      for (int k = 0; k < b.size(); k++)
	x += a(i,k)*b[k];
      c[i] = x;
    }
    return c;
  }
};

#endif /* DMAT_H */
