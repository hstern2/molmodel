#ifndef CMPLX_H
#define CMPLX_H

#include <cmath>
#include "cart.h"
#include "out.hpp"

class Complex : public complex_t
{
public:
  Complex() { }
  Complex(double r, double i) { re = r; im = i; }
  Complex & operator+=(const Complex &c)
  { re += c.re; im += c.im; return *this; }
  Complex & operator-=(const Complex &c)
  { re -= c.re; im -= c.im; return *this; }
  Complex & operator*=(double r)
  { re *= r; im *= r; return *this; }
  Complex & operator/=(double r)
  { re /= r; im /= r; return *this; }
  Complex operator-() const
  { return Complex(-re,-im); }
  Complex & apply(double (*f)(double))
  { re = f(re); im = f(im); return *this; }
  Complex map(double (*f)(double)) const
  { return Complex(f(re),f(im)); }
  friend Complex operator+(const Complex &c1, const Complex &c2)
  { return Complex(c1.re + c2.re, c1.im + c2.im); }
  friend Complex operator-(const Complex &c1, const Complex &c2)
  { return Complex(c1.re - c2.re, c1.im - c2.im); }
  friend Complex operator+(const Complex &c1, double c2)
  { return Complex(c1.re + c2, c1.im); }
  friend Complex operator-(const Complex &c1, double c2)
  { return Complex(c1.re - c2, c1.im); }
  friend Complex operator+(double c2, const Complex &c1)
  { return Complex(c2 + c1.re, c1.im); }
  friend Complex operator-(double c2, const Complex &c1)
  { return Complex(c2 - c1.re, c1.im); }
  friend Complex operator*(const Complex &a, const Complex &b)
  { return Complex(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re); }
  friend Complex operator*(double r, const Complex &c)
  { return Complex(r*c.re, r*c.im); }
  friend Complex operator*(const Complex &c, double r)
  { return r*c; }
  friend Complex operator/(const Complex &c, double r)
  { return Complex(c.re/r, c.im/r); }
};

inline static Complex expi(const double &t)
{
  return Complex(cos(t),sin(t));
}

inline static ostream & operator<<(ostream &s, const Complex &c) 
{ 
  return s << c.re << " " << c.im;
}

inline static istream & operator>>(istream &s, Complex &c)
{ 
  return s >> c.re >> c.im;
}


#endif /* CMPLX_H */
