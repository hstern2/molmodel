#ifndef MMIN_H
#define MMIN_H

#include "out.hpp"

struct MSys;
class Potential;
class Constraint;

class MMin
{

public:

  MMin(MSys &, Potential &, Constraint *);
  int maximum_evaluations, verbose; // in
  double tolerance, linesearch_tolerance; // in
  void run(); // in

  MSys &msys;
  Potential &p;
  Constraint *c;

  classIO(MMin);
  
};

#endif /* MMIN_H */
