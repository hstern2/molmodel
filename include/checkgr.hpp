#ifndef CHECKGR_H
#define CHECKGR_H

#include "pot.hpp"

class Constraint;

class CheckGradients
{
public:
  CheckGradients(MSys &m, Potential &p_, Constraint *c_) : 
    verbose(0), trials(10), displacement(1e-8), 
    msys(m), p(p_), c(c_) { }
  int verbose, trials; // in
  double displacement; // in
  void run();
private:
  MSys &msys;
  Potential &p;
  Constraint *c;
  classIO(CheckGradients);
};

#endif /* CHECKGR_H */
