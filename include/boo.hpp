#ifndef BOO_H
#define BOO_H

#include <cstdlib>
#include "out.hpp"
#include "str.hpp"

class Boolean
{
public:
  Boolean(bool v = false) : val(v) { }
  Boolean(const Boolean &b) : val(b.val) { }
  Boolean & operator=(const Boolean &b)
  {
    val = b.val;
    return *this;
  }
  Boolean & operator=(bool i)
  {
    val = i;
    return *this;
  }
  operator bool() const { return val; }
private:
  bool val;
  friend istream & operator>>(istream &s, Boolean &b)
  {
    Str buf;
    s >> buf;
    const char *c = buf;
    if (!strcasecmp(c,"y") || !strcasecmp(c,"yes") ||
	!strcasecmp(c,"t") || !strcasecmp(c,"true") ||
	!strcmp(c,"1"))
      b = 1;
    else if (!strcasecmp(c,"n") || !strcasecmp(c,"no") ||
	     !strcasecmp(c,"f") || !strcasecmp(c,"false") ||
	     !strcmp(c,"0"))
      b = 0;
    else
      die("Boolean: expecting 'yes'/'true' or 'no'/'false': read %s", c);
    return s;
  }
  friend ostream & operator<<(ostream &s, const Boolean &b)
  { return s << (b ? "yes" : "no"); }
};

#endif /* BOO_H */
