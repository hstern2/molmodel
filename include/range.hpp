#ifndef RANGE_H
#define RANGE_H

#include "out.hpp"

struct Range
{
  int lo, hi;
  int size() const { return hi - lo + 1; }
  friend istream & operator>>(istream &s, Range &c)
  {
    char a, b;
    s >> c.lo >> a >> b >> c.hi;
    if (!(s && a == '.' && b == '.'))
      die("Range: expecting n1..n2");
    if (c.lo > c.hi) {
      int tmp = c.lo;
      c.lo = c.hi;
      c.hi = tmp;
    }
    return s;
  }
  friend ostream & operator<<(ostream &s, Range &c)
  {
    return s << c.lo << ".." << c.hi;
  }
};

typedef Lst<Range> LstRange;

#endif /* RANGE_H */
