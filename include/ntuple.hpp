#ifndef NTUPLE_H
#define NTUPLE_H

#include "str.hpp"
#include "hash.h"

template<class T> struct Pair
{
  Pair() { }
  Pair(T a_, T b_) : a(a_), b(b_) { }
  T a, b;
};

template<class T> inline bool operator==(const Pair<T> &i, const Pair<T> &j)
{ return i.a == j.a && i.b == j.b; }

template<class T> inline bool operator<(const Pair<T> &i, const Pair<T> &j)
{ 
  if (i.a < j.a)
    return true;
  if (i.a > j.a)
    return false;
  if (i.b < j.b)
    return true;
  return false;
}

template<class T> inline bool operator>(const Pair<T> &i, const Pair<T> &j)
{ 
  if (i.a > j.a)
    return true;
  if (i.a < j.a)
    return false;
  if (i.b > j.b)
    return true;
  return false;
}

template<class T> inline istream & operator>>(istream &s, Pair<T> &p)
{ return s >> p.a >> p.b; }

template<class T> inline ostream & operator<<(ostream &s, const Pair<T> &p)
{ return s << p.a << " " << p.b << " "; }

template<class T> inline unsigned hash_val(const Pair<T> &key)
{ return (unsigned) (1223*hash_val(key.a) + hash_val(key.b)); }

template<class T> struct Triplet
{
  Triplet() { }
  Triplet(T a_, T b_, T c_) : a(a_), b(b_), c(c_) { }
  T a, b, c;
};

template<class T> inline bool operator==(const Triplet<T> &i, const Triplet<T> &j)
{ return i.a == j.a && i.b == j.b && i.c == j.c; }

template<class T> inline bool operator<(const Triplet<T> &i, const Triplet<T> &j)
{ 
  if (i.a < j.a)
    return true;
  if (i.a > j.a)
    return false;
  if (i.b < j.b)
    return true;
  if (i.b > j.b)
    return false;
  if (i.c < j.c)
    return true;
  return false;
}

template<class T> inline bool operator>(const Triplet<T> &i, const Triplet<T> &j)
{ 
  if (i.a > j.a)
    return true;
  if (i.a < j.a)
    return false;
  if (i.b > j.b)
    return true;
  if (i.b < j.b)
    return false;
  if (i.c > j.c)
    return true;
  return false;
}

template<class T> inline istream & operator>>(istream &s, Triplet<T> &p)
{ return s >> p.a >> p.b >> p.c; }

template<class T> inline ostream & operator<<(ostream &s, const Triplet<T> &p)
{ return s << p.a << " " << p.b << " " << p.c << " "; }

template<class T> inline unsigned hash_val(const Triplet<T> &key)
{ return (unsigned) (1223*(1223*hash_val(key.a) + hash_val(key.b) + hash_val(key.c))); }

template<class T> struct Quadruplet
{
  Quadruplet() { }
  Quadruplet(T a_, T b_, T c_, T d_) : a(a_), b(b_), c(c_), d(d_) { }
  T a, b, c, d;
};

template<class T> inline istream & operator>>(istream &s, Quadruplet<T> &c) 
{ return s >> c.a >> c.b >> c.c >> c.d; }

template<class T> inline ostream & operator<<(ostream &s, const Quadruplet<T> &c)
{ return s << c.a << " " << c.b << " " << c.c << " " << c.d << " "; }

template<class T> inline bool operator==(const Quadruplet<T> &i, const Quadruplet<T> &j)
{ return i.a == j.a && i.b == j.b && i.c == j.c && i.d == j.d; }

template<class T> inline bool operator<(const Quadruplet<T> &i, const Quadruplet<T> &j)
{ 
  if (i.a < j.a)
    return true;
  if (i.a > j.a)
    return false;
  if (i.b < j.b)
    return true;
  if (i.b > j.b)
    return false;
  if (i.c < j.c)
    return true;
  if (i.c > j.c)
    return false;
  if (i.d < j.d)
    return true;
  return false;
}

template<class T> inline bool operator>(const Quadruplet<T> &i, const Quadruplet<T> &j)
{ 
  if (i.a > j.a)
    return true;
  if (i.a < j.a)
    return false;
  if (i.b > j.b)
    return true;
  if (i.b < j.b)
    return false;
  if (i.c > j.c)
    return true;
  if (i.c < j.c)
    return false;
  if (i.d > j.d)
    return true;
  return false;
}

template<class T> inline unsigned hash_val(const Quadruplet<T> &key)
{ return (unsigned) (1223*(1223*(1223*hash_val(key.a) + hash_val(key.b)) + hash_val(key.c)) + hash_val(key.d)); }

typedef Pair<int> Int2;
typedef Triplet<int> Int3;
typedef Quadruplet<int> Int4;

typedef Pair<double> Double2;
typedef Triplet<double> Double3;
typedef Quadruplet<double> Double4;

typedef Pair<Str> Str2;
typedef Triplet<Str> Str3;
typedef Quadruplet<Str> Str4;

#endif /* NTUPLE_H */
