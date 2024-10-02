#ifndef HARMONIC_H
#define HARMONIC_H

#include "pot.hpp"
#include "str.hpp"
#include "dmat.hpp"
#include "cvec.hpp"
#include "zmatrix.hpp"
#include "htab.hpp"
#include "ntuple.hpp"

typedef HSet<int> HSetInt;

class Harmonic : public Potential
{
  Lst<Potential *> pp;  
public:
  Harmonic();
  ~Harmonic();

  Boolean ignore_dihedrals_with_variance_greater_than_pi2_over_4; // io
  ZMatrix zmatrix; // io
  HSetInt ignore_dihedral; // io 
  double temperature; // in
  int verbose;  // io
  double E0; // io
  void init(MSys *);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  classIO(Harmonic); 
  DMat covariance_matrix, harmonic_matrix; // in
  DVec average_position; // in
private:
  void detect_dihedrals_to_ignore();
  DMat jacobian_matrix;
  DVec delta_position;
  int nfreedom; 
  double u; // just the energy for Harmonic potential
  void cleanup();
};
#endif /* HARMONIC_H */

