#ifndef DVEC_H
#define DVEC_H

#include <cstring>
#include "vec.hpp"

class DVec : public Vec<double>
{
public:
  DVec(int sz = 0) : Vec<double>(sz) { } 
  DVec(const DVec &v) : Vec<double>(v) { }
  DVec(int sz, const double &t) : Vec<double>(sz,t) { }
  DVec & operator=(const DVec &v) 
  { return (DVec &) Vec<double>::operator=(v); }
  DVec & apply(double (*f)(double))
  { 
    for (int i = 0; i < size(); i++)
      (*this)[i] = f((*this)[i]);
    return *this;
  }
  DVec map(double (*f)(double)) const
  {
    DVec v(size());
    for (int i = 0; i < size(); i++)
      v[i] = f((*this)[i]);
    return v;
  }
  operator double *() { return (double *) p; }
  operator const double *() const { return (const double *) p; }
  operator const double *() { return (const double *) p; }
  DVec & copy(const double *v)
  {
    memcpy((double *) *this, v, size()*sizeof(double));
    return *this;
  }
  DVec copy() const
  {
    DVec v(size());
    v.copy(*this);
    return v;
  }
  void zero()
  { memset((double *) *this, 0, size()*sizeof(double)); }
  friend double operator*(const DVec a, const DVec b)
  {
    double c = 0.0;
    assert(a.is_conformable_with(b));
    for (int i = 0; i < a.size(); i++)
      c += a[i]*b[i];
    return c;
  }
  double sq() const { return (*this)*(*this); }
  DVec & operator+=(const DVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] += a[i];
    return *this;
  }
  DVec & operator-=(const DVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] -= a[i];
    return *this;
  }
  DVec & operator*=(double x)
  {
    for (int i = 0; i < size(); i++)
      (*this)[i] *= x;
    return *this;
  }
  double sum() const
  {
    double c = 0;
    for (int i = 0; i < size(); i++)
      c += (*this)[i];
    return c;
  }
  double average() const
  { return sum()/size(); }
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

#endif /* DVEC_H */
