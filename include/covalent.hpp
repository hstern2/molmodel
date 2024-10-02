#ifndef COVALENT_H
#define COVALENT_H

#include "pot.hpp"

class Covalent : public Potential
{
public:
  HTab<double,Str2> param; // io
  int verbose; // io
  void init(MSys *);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  void types_changed();
  void types_changed(const MSys &perturbed);
  void set_lambda(double lambda);
  Covalent() : verbose(1), u1(0), u2(0), u(0) { }
  classIO(Covalent);
private:
  double u1, u2, u;
};

#endif /* COVALENT_H */
