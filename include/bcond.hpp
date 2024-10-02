#ifndef BCOND_H
#define BCOND_H

#include "coord.hpp"
#include "bcondtype.h"
#include "fns.h"

enumIO(bcondtype_t);

class BoundaryConditions
{
  
public:
  
  static BoundaryConditions *new_non_periodic();
  static BoundaryConditions *new_cubic(double);
  static BoundaryConditions *new_orthorhombic(Cartesian);
  static BoundaryConditions *new_bcc(double);
  static BoundaryConditions *new_fcc(double);
  static BoundaryConditions *new_triclinic(const Tensor &);
  static BoundaryConditions *new_triclinic(double a, double b, double c,
					   double alpha, double beta, double gamma);
  static BoundaryConditions *new_from_cell(double a, double b, double c,
					   double alpha, double beta, double gamma);
  
  const bcondtype_t type;

  BoundaryConditions(bcondtype_t t) : type(t) { }
  virtual ~BoundaryConditions();
  virtual BoundaryConditions *copy() const = 0;
  virtual void map_to_central_box(Cartesian &) const = 0;
  virtual void nearest_lattice_point(const Cartesian &, int &i, int &j, int &k) const = 0;
  virtual double volume() const = 0;
  virtual double min_diameter() const = 0; /* diameter of inscribed sphere */
  virtual double max_diameter() const = 0; /* diameter of circumscribed sphere */
  virtual void scale(double) = 0;
  virtual Tensor lattice_vectors() const = 0;
  Tensor reciprocal_lattice_vectors() const;
  virtual void write(ostream &) const = 0;

  Cartesian & minimum_image_displacement(Cartesian &ab, const Cartesian &a, const Cartesian &b) const
  {
    ab.x = a.x - b.x;
    ab.y = a.y - b.y;
    ab.z = a.z - b.z;
    map_to_central_box(ab);
    return ab;
  }
  
  double minimum_image_distance(const Cartesian &a, const Cartesian &b) const
  {
    Cartesian r;
    return minimum_image_displacement(r,a,b).magnitude();
  }
  
  double square_minimum_image_distance(const Cartesian &a, const Cartesian &b) const
  {
    Cartesian r;
    return minimum_image_displacement(r,a,b).sq();
  }
  
  int is_in_central_box(const Cartesian &a) const
  {
    Cartesian r = a;
    map_to_central_box(r);
    r -= a;
    return is_almost_zero(r.sq());
  }
  
};

istream & operator>>(istream &, BoundaryConditions *&);
ostream & operator<<(ostream &, const BoundaryConditions *);

#endif /* BCOND_H */
