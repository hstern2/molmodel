#ifndef RMAT_H
#define RMAT_H

#include <cstring>
#include <cmath>
#include "mat.hpp"
#include "vec.hpp"

class RMat;

class RVec : public Vec<double>
{
public:
  RVec(int sz = 0) : Vec<double>(sz) { } 
  RVec(const RVec &v) : Vec<double>(v) { }
  RVec(const Vec<double> &v) : Vec<double>(v) { }
  RVec(int sz, const double &t) : Vec<double>(sz,t) { }
  RVec & operator=(const RVec &v)
  { 
    Vec<double>::operator=(v);
    return *this;
  }
  RVec & apply(double (*f)(double))
  { 
    for (int i = 0; i < size(); i++)
      (*this)[i] = f((*this)[i]);
    return *this;
  }
  RVec map(double (*f)(double)) const
  {
    RVec vec(size());
    for (int i = 0; i < size(); i++)
      vec[i] = f((*this)[i]);
    return vec;
  }
  RVec & operator+=(const RVec &t)
  { 
    assert(is_conformable_with(t));
    for (int i = 0; i < size(); i++)
      (*this)[i] += t[i];
    return *this;
  }
  RVec & operator-=(const RVec &t)
  { 
    assert(is_conformable_with(t));
    for (int i = 0; i < size(); i++)
      (*this)[i] -= t[i];
    return *this;
  }
  RVec & operator*=(double r)
  { 
    for (int i = 0; i < size(); i++)
      (*this)[i] *= r;
    return *this;
  }
  RVec & operator/=(double r)
  { 
    for (int i = 0; i < size(); i++)
      (*this)[i] /= r;
    return *this;
  }
  RVec operator-() const
  { 
    RVec v(size());
    for (int i = 0; i < size(); i++)
      v[i] = -(*this)[i];
    return v;
  }
  RVec & copy(const RVec &vec)
  {
    assert(is_conformable_with(vec));
    memcpy((double *) *this, (const double *) vec, size()*sizeof(double));
    return *this;
  }
  RVec & copy(const double *vec)
  { 
    memcpy((double *) *this, vec, size()*sizeof(double));
    return *this;
  }
  const RVec & copy_to(double *vec) const
  { 
    memcpy(vec, (const double *) *this, size()*sizeof(double));
    return *this;
  }
  RVec & copy(const RVec &vec, int start, int end)
  {
    assert(size() >= end - start + 1);
    assert(0 <= start && start < vec.size());
    assert(0 <= end && end < vec.size());
    memcpy(*this, &vec[start], (end - start + 1) * sizeof(double));
    return *this;
  }
  RVec copy() const
  {
    RVec r(size());
    r.copy(*this);
    return r;
  }
  void zero()
  { memset(*this, 0, size()*sizeof(double)); }
  void garbage()
  { memset(*this, 0xFF, size()*sizeof(double)); }
  double sum()
  {
    double s = 0.0;
    for (int i = 0; i < size(); i++)
      s += (*this)[i];
    return s;
  }
  double sq() const
  {
    double c = 0;
    for (int i = 0; i < size(); i++)
      c += (*this)[i] * (*this)[i];
    return c;
  }
  double minimum() const
  {
    if (size() == 0)
      return 0;
    double m = (*this)[0];
    for (int i = 1; i < size(); i++)
      if ((*this)[i] < m)
	m = (*this)[i];
    return m;
  }
  double maximum() const
  {
    if (size() == 0)
      return 0;
    double m = (*this)[0];
    for (int i = 1; i < size(); i++)
      if ((*this)[i] > m)
	m = (*this)[i];
    return m;
  }
};

inline RVec operator+(const RVec &a, const RVec &b)
{
  RVec c(a.copy());
  return c += b;
}

inline RVec operator-(const RVec &a, const RVec &b)
{
  RVec c(a.copy());
  return c -= b;
}

inline RVec operator*(double a, const RVec &b)
{
  RVec c(b.copy());
  return c *= a;
}

inline RVec operator*(const RVec &a, double b)
{
  RVec c(a.copy());
  return c *= b;
}

inline double operator*(const RVec &a, const RVec &b)
{
  assert(a.is_conformable_with(b));
  double c = 0;
  for (int i = 0; i < a.size(); i++)
    c += a[i]*b[i];
  return c;
}

class RMat : public Mat<double>
{
public:
  RMat(int n = 0, int m = 0) : Mat<double>(n, m) { }
  RMat(const RMat &mat) : Mat<double>(mat) { }
  RMat(const Mat<double> &mat) : Mat<double>(mat) { }
  RMat(int n, int m, const double &t) : Mat<double>(n, m, t) { }
  RMat & operator=(const RMat &mat)
  { 
    Mat<double>::operator=(mat); 
    return *this;
  }
  RMat transpose() const 
  { return RMat(Mat<double>::transpose()); }
  RMat & copy(const RMat &mat)
  {
    assert(is_conformable_with(mat));
    memcpy(rep->p, mat.rep->p, rows()*columns()*sizeof(double));
    return *this;
  }
  RMat copy() const
  {
    RMat r(rows(), columns());
    r.copy(*this);
    return r;
  }
  int is_symmetric() const
  {
    if (!is_square())
      return 0;
    for (int i = 0; i < rows(); i++)
      for (int j = i+1; j < columns(); j++)
	if ((*this)(i,j) != (*this)(j,i))
	  return 0;
    return 1;
  }
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
  RMat & apply(double (*f)(double))
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) = f((*this)(i,j));
    return *this;
  }
  RMat & operator+=(const RMat &t)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) += t(i,j);
    return *this;
  }
  RMat & operator-=(const RMat &t)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) -= t(i,j);
    return *this;
  }
  RMat & operator*=(double r)
  {
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	(*this)(i,j) *= r;
    return *this;
  }
  RMat operator-() const
  { 
    RMat m(rows(), columns());
    for (int i = 0; i < rows(); i++)
      for (int j = 0; j < columns(); j++)
	m(i,j) = -(*this)(i,j);
    return m;
  }
  void zero()
  { memset(rep->p, 0, rows()*columns()*sizeof(double)); }
  void garbage()
  { memset(rep->p, 0xFF, rows()*columns()*sizeof(double)); }
};

inline RMat operator+(const RMat &a, const RMat &b)
{
  RMat c = a.copy();
  return c += b;
}

inline RMat operator-(const RMat &a, const RMat &b)
{
  RMat c = a.copy();
  return c -= b;
}

inline RMat operator*(double a, const RMat &b)
{
  RMat c = b.copy();
  return c *= a;
}

inline RMat operator*(const RMat &a, double b)
{
  RMat c = a.copy();
  return c *= b;
}

inline RMat operator*(const RMat &a, const RMat &b)
{
  assert(a.columns() == b.rows());
  RMat c(a.rows(), b.columns());
  for (int i = 0; i < a.rows(); i++)
    for (int j = 0; j < b.columns(); j++) {
      double x = 0;
      for (int k = 0; k < a.columns(); k++)
	x += a(i,k)*b(k,j);
      c(i,j) = x;
    }
  return c;
}

inline RVec operator*(const RVec &a, const RMat &b)
{
  assert(a.size() == b.rows());
  RVec c(b.columns());
  for (int j = 0; j < b.columns(); j++) {
    double x = 0;
    for (int k = 0; k < a.size(); k++)
      x += a[k]*b(k,j);
    c[j] = x;
  }
  return c;
}

inline RVec operator*(const RMat &a, const RVec &b)
{
  assert(a.columns() == b.size());
  RVec c(a.rows());
  for (int i = 0; i < a.rows(); i++) {
    double x = 0;
    for (int k = 0; k < b.size(); k++)
      x += a(i,k)*b[k];
    c[i] = x;
  }
  return c;
}

inline ostream & operator<<(ostream &s, const RVec &v)
{
  s << v.size() << "\n";
  for (int i = 0; i < v.size(); i++)
    s << v[i] << " ";
  s << "\n";
  return s;
}

inline ostream & operator<<(ostream &s, const RMat &m)
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

#endif /* RMAT_H */
