#ifndef STR_H
#define STR_H

#include <cstring>
#include <cstdlib>
#include <cctype>
#include "out.hpp"

class Str
{
  struct StrRep
  {
    int ref;
    char *p;
    StrRep(const char *s) : ref(1), p(new char[strlen(s)+1]) { strcpy(p,s); }
    ~StrRep() { delete[] p; p = 0; }
  };
  StrRep *rep;
  void destroy() 
  { 
    if (--rep->ref == 0) { 
      delete rep; 
      rep = 0; 
    } 
  }
public:
  operator char*() { return rep->p; }
  operator const char*() const { return rep->p; }
  operator const char*() { return rep->p; }
  ~Str() { destroy(); }
  Str() : rep(new StrRep("\0")) { }
  Str(const Str &s) : rep(s.rep) { rep->ref++; }
  Str(const char *s) : rep(new StrRep(s)) { }
  Str & operator=(const Str &s)
  {
    s.rep->ref++;
    destroy();
    rep = s.rep;
    return *this;
  }
  
  Str & operator=(const char *v) { return *this = Str(v); }
  Str copy() const { return Str((const char *) (*this)); }

  friend Str operator+(const Str s1, const char *s2)
  {
    char *tmp = new char[strlen(s1)+strlen(s2)+1];
    strcpy(tmp,s1);
    strcat(tmp,s2);
    Str s(tmp);
    delete[] tmp;
    return s;
  }
  
  friend Str operator+(const char *s1, const Str s2)
  {
    char *tmp = new char[strlen(s1)+strlen(s2)+1];
    strcpy(tmp,s1);
    strcat(tmp,s2);
    Str s(tmp);
    delete[] tmp;
    return s;
  }

  friend Str operator+(const Str s1, const Str s2)
  {
    return s1 + (const char *) s2;
  }
 
#define STRBUFSIZE 64 
  friend istream & operator>>(istream &c, Str &s)
  {
    char q, buf[STRBUFSIZE];
    char *b = buf;
    c >> q;
    s = "";
    if (q == '"') { // quoted string
      while(c.get(*b) && *b != '"') {
	b++;
	if (b - buf + 1 >= STRBUFSIZE) {
	  *b = '\0';
	  s += buf;
	  b = buf;
	}
      }
    } else {        // unquoted string
      c.putback(q);
      while(c.get(*b) && !isspace(*b) && *b != '{' && *b != '}') {
	b++;
	if (b - buf + 1 >= STRBUFSIZE) {
	  *b = '\0';
	  s += buf;
	  b = buf;
	}
      }
      if (c)
	c.putback(*b);
    }
    *b = '\0';
    s += buf;
    return c;
  }
#undef STRBUFSIZE  
  friend ostream & operator<<(ostream &os, const Str &s)
  { 
    const char *c = s;
    while (*c && !isspace(*c))
      c++;
    if (isspace(*c))
      return os << '"' << (const char *) s << '"';
    else
      return os << (const char *) s;
  }

  Str & operator+=(const char *s)
  {
    return *this = *this + s;
  }
  
  
  friend bool operator==(const Str &a, const Str &b) { return !strcmp(a,b); }
  friend bool operator!=(const Str &a, const Str &b) { return strcmp(a,b); }

  friend bool operator<(const Str &s1, const Str &s2)
  { return strcmp(s1, s2) < 0; }
  
  friend bool operator>(const Str &s1, const Str &s2)
  { return strcmp(s1, s2) > 0; }

  friend const Str &max(const Str &s1, const Str &s2)
  { return s1 > s2 ? s1 : s2; }

  friend const Str &min(const Str &s1, const Str &s2)
  { return s1 < s2 ? s1 : s2; }

};

#endif /* STR_H */
