#ifndef APATTERN_H
#define APATTERN_H

#include "vec.hpp"
#include "hset.hpp"
#include "htab.hpp"

class AtomTest;
class Atom;

class AtomPattern
{
public:
  int ntest;
  AtomPattern(const char *pattern);
  AtomPattern(const AtomPattern &); /* dummy routine */
  AtomPattern &operator=(const AtomPattern &); /* dummy routine */
  virtual ~AtomPattern();
  void run(Vec<Atom> &, const HTab<HSet<int> > &atoms_with_property);
private:
  virtual void succeed(const int *path, Atom *atom);
  void search(int a, int i, HSet<int> &seen, int *path, Atom *atom);
  Vec<AtomTest *> test;
  Vec<int> next;
};

#endif /* APATTERN_H */
