#ifndef PPOT_H
#define PPOT_H

#include "pot.hpp"
#include "htab.hpp"
#include "array.h"
#include "cvec.hpp"
  
struct PairPotentialParam
{
  double c2, c4, c6, c8;  // io
  bool matches(const PairPotentialParam &) const;
  classIO(PairPotentialParam);
};

extern "C" { 
  struct PairBinary;
}
class NeighborList;

class PairPotential : public Potential
{
public:
  PairPotential();
  ~PairPotential();
  HTab<PairPotentialParam,Str2> covalent_param, non_covalent_param; // io
  void init(MSys *m);
  void add_to_energy_and_forces(double &, Cartesian *);
  void add_to_hessian(double **);
  bool can_calculate_hessian() const;
  void write() const;
private:
  CVec r;
  double rc, ener;
  Vec<int> type;
  const void ***covp, ***noncovp;
  NeighborList *nlist;
  struct PairBinary *covalent;
  void add_to_energy_and_forces(double &, Cartesian *, int, int,
				const BoundaryConditions &, double);
  void add_to_hessian(double **h, int, int, const BoundaryConditions &, double);
  void pair_function(int, int, double, double &, double &);
  void pair_function(int, int, double, double &, double &, double &);
  classIO(PairPotential);
};

#endif /* PPOT_H */
