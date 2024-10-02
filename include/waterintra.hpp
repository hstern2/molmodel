#ifndef WATERINTRA_H
#define WATERINTRA_H

#include "pot.hpp"

/***
 * Intramolecular potential energy surface for water
 * from Partridge and Schwenke, JCP 106, 4618 (1996)
 ***/
class WaterIntra : public Potential
{
public:
  WaterIntra();
  void init(MSys *m);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  int verbose; // io
private:
  double ener;
  Vec<Int3> waters; // out
  classIO(WaterIntra);
};

#endif /* WATERINTRA_H */
