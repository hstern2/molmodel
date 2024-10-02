#ifndef HIST_H
#define HIST_H

#include "out.hpp"
#include "boo.hpp"
#include "vec.hpp"

class Histogram
{
public:
  int nbin; // in
  double lo; // in
  double hi; // in
  double mean() const;
  double variance() const;
  Boolean is_dynamic; // in
  Histogram();
  Histogram(const Histogram &);
  void init();
  void update(double);
  void write(const char *, const char * = 0) const;
  void write(ostream &) const;
  /* For input/output */
  friend istream & operator>>(istream &, Histogram &);
  int read_field(istream &, const char *);
  void write_field_names() const;
private:
  int ncount;
  Vec<int> v;
  void update_v(double);
};

#endif /* HIST_H */
