#ifndef ZMATRIX_H
#define ZMATRIX_H

#include "vec.hpp"

struct ZMatrixEntry
{
  int i, a, b, c; // in
  double r, theta, phi; // in
  ZMatrixEntry() : i(-1), a(-1), b(-1), c(-1), r(0), theta(0), phi(0) { }
  classIO(ZMatrixEntry);
};

class Atom;
class DMat;
class DVec;

struct ZMatrix : public Vec<ZMatrixEntry>
{
  void check(int n) const;
  void to_cartesian(Atom *) const;
  void from_cartesian(const Atom *);
  int number_of_bonds() const;
  int number_of_angles() const;
  int number_of_dihedrals() const;
  int degrees_of_freedom() const;
  double jacobian() const; /* based on current values of r, theta, phi */
  void check_ignore_dihedral(const HSet<int> &ignore_dihedral) const;
  void to_internal(DVec &, const HSet<int> &ignore_dihedral) const;
  void from_internal(const DVec &, const HSet<int> &ignore_dihedral);
  void harmonic_partition_function(const DVec &avgpos,
				   const DMat &cov,
				   const HSet<int> &ignore_dihedral,
				   double &logz0,
				   double &logz2,
				   double &logz4);
  
  friend ostream & operator<<(ostream &, const ZMatrix &);
};

#endif /* ZMATRIX_H */
