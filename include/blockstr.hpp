#ifndef BLOCKSTR_H
#define BLOCKSTR_H

#include "str.hpp"

class BlockStr : public Str
{
public:
  BlockStr() { }
  BlockStr(const Str &s) : Str(s) { }
  BlockStr(const char *s) : Str(s) { }
  BlockStr & operator=(const Str &s) { 
    ((Str &) (*this)) = s;
    return *this;
  }
  BlockStr & operator=(const char *v) { 
    ((Str &) (*this)) = v;
    return *this;
  }
  friend istream & operator>>(istream &c, BlockStr &s)
  {
    char buf[64];
    char *b = buf;
    c >> *b;
    if (*b != '{')
      die("Str::read_block: expecting '{'");
    b++;
    int n = 1;
    s = "";
    while (c.get(*b) && n > 0) {
      if (*b == '{')
	n++;
      else if (*b == '}')
	n--;
      b++;
      if (b - buf >= 64) {
	*b = '\0';
	s += buf;
	b = buf;
      }
    }
    *b = '\0';
    s += buf;
    return c;
  }
};

#endif /* BLOCKSTR_H */
