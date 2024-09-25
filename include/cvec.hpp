#ifndef CVEC_H
#define CVEC_H

#include <cstring>
#include "vec.hpp"
#include "coord.hpp"

class CVec : public Vec<Cartesian>
{
public:
  CVec(int sz = 0) : Vec<Cartesian>(sz) { } 
  CVec(const CVec &v) : Vec<Cartesian>(v) { }
  CVec(int sz, const Cartesian &t) : Vec<Cartesian>(sz,t) { }
  CVec(int sz, Cartesian *q) : Vec<Cartesian>(sz,q) { }
  CVec & operator=(const CVec &v) 
  { return (CVec &) Vec<Cartesian>::operator=(v); }
  CVec & apply(double (*f)(double))
  { 
    for (int i = 0; i < size(); i++)
	(*this)[i].apply(f);
    return *this;
  }
  CVec map(Cartesian (Cartesian::*f)() const) const
  {
    CVec v(size());
    for (int i = 0; i < size(); i++)
      v[i] = ((*this)[i].*f)();
    return v;
  }

  operator cart_t *() { return (cart_t *) p; }
  operator const cart_t *() const { return (const cart_t *) p; }
  operator const cart_t *() { return (const cart_t *) p; }

  operator double *() { return (double *) p; }
  operator const double *() const { return (const double *) p; }
  operator const double *() { return (const double *) p; }

  operator void *() { return (void *) p; }
  operator const void *() const { return (const void *) p; }
  operator const void *() { return (const void *) p; }

  CVec as_unit_vector() const
  { return map(&Cartesian::as_unit_vector); }
  CVec & copy(const Cartesian *v)
  {
    memcpy((Cartesian *) *this, v, size()*sizeof(Cartesian));
    return *this;
  }
  CVec copy() const
  {
    CVec v(size());
    v.copy(*this);
    return v;
  }
  void zero()
  { memset((Cartesian *) *this, 0, size()*sizeof(Cartesian)); }
  friend double operator*(const CVec a, const CVec b)
  {
    double c = 0.0;
    assert(a.is_conformable_with(b));
    for (int i = 0; i < a.size(); i++)
      c += a[i]*b[i];
    return c;
  }
  double sq() const { return (*this)*(*this); }
  CVec & operator+=(const CVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] += a[i];
    return *this;
  }
  CVec & operator-=(const CVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] -= a[i];
    return *this;
  }
  CVec & operator*=(double x)
  {
    for (int i = 0; i < size(); i++)
      (*this)[i] *= x;
    return *this;
  }
  Cartesian sum() const
  {
    Cartesian c(0,0,0);
    for (int i = 0; i < size(); i++)
      c += (*this)[i];
    return c;
  }
  Cartesian average() const
  { return sum()/size(); }
  Tensor outer_product(const CVec bv) const
  {
    assert(is_conformable_with(bv));
    Tensor t(0);
    const int n = size();
    const Cartesian *a = *this;
    const Cartesian *b = bv;
    for (int i = 0; i < n; i++) {
      t.xx += a[i].x * b[i].x;
      t.xy += a[i].x * b[i].y;
      t.xz += a[i].x * b[i].z;
      t.yx += a[i].y * b[i].x;
      t.yy += a[i].y * b[i].y;
      t.yz += a[i].y * b[i].z;
      t.zx += a[i].z * b[i].x;
      t.zy += a[i].z * b[i].y;
      t.zz += a[i].z * b[i].z;
    }
    return t;
  }
};
  
#endif /* CVEC_H */
