#include "bcond.hpp"
#include "pbc.h"
#include "str.hpp"
#include "units.h" 
#include "geograd.hpp"

class NonPeriodic : public BoundaryConditions
{
public:
  NonPeriodic() : BoundaryConditions(non_periodic) { }
  BoundaryConditions *copy() const { return new NonPeriodic; }
  void map_to_central_box(Cartesian &) const { }
  void nearest_lattice_point(const Cartesian &, int &, int &, int &) const { }
  double volume() const { return 0; }
  double min_diameter() const { return 0; }
  double max_diameter() const { return 0; }
  void scale(double) { }
  Tensor lattice_vectors() const { return Tensor(0); }
  void write(ostream &s) const { s << "non_periodic"; }
};

class Cubic : public BoundaryConditions
{
public:
  Cubic(double b) : BoundaryConditions(cubic), boxl(b) { }
  BoundaryConditions *copy() const { return new Cubic(boxl); }
  void map_to_central_box(Cartesian &c) const 
  { 
    cubic_map_to_central_box(&c.x, &c.y, &c.z, boxl, 1/boxl);
  }
  void nearest_lattice_point(const Cartesian &c, int &i, int &j, int &k) const
  {
    cubic_nearest_lattice_point(c.x, c.y, c.z, boxl, &i, &j, &k);
  }
  double volume() const { return cube(boxl); }
  double min_diameter() const { return boxl; }
  double max_diameter() const { return sqrt(3.0) * boxl; }
  void scale(double s) { boxl *= s; }
  Tensor lattice_vectors() const { return Tensor(boxl); }
  void write(ostream &s) const { s << "cubic " << min_diameter(); }
private:
  double boxl;
};

class Orthorhombic : public BoundaryConditions
{
public:
  Orthorhombic(Cartesian b) : BoundaryConditions(orthorhombic), boxl(b) { }
  BoundaryConditions *copy() const { return new Orthorhombic(boxl); }
  void map_to_central_box(Cartesian &c) const 
  { 
    orthorhombic_map_to_central_box(&c.x, &c.y, &c.z, boxl.x, boxl.y, boxl.z,
				    1/boxl.x, 1/boxl.y, 1/boxl.z);
  }
  void nearest_lattice_point(const Cartesian &c, int &i, int &j, int &k) const
  {
    orthorhombic_nearest_lattice_point(c.x, c.y, c.z, boxl.x, boxl.y, boxl.z, &i, &j, &k);
  }
  double volume() const { return boxl.x*boxl.y*boxl.z; }
  double min_diameter() const { return min(min(boxl.x,boxl.y),boxl.z); }
  double max_diameter() const { return boxl.magnitude(); }
  void scale(double s) { boxl *= s; }
  Tensor lattice_vectors() const { return Tensor(boxl); }
  void write(ostream &s) const { s << "orthorhombic " << boxl; }
private:
  Cartesian boxl;
};

class BCC : public BoundaryConditions
{
public:
  BCC(double boxl) : BoundaryConditions(bcc), side(boxl/(0.5*sqrt(3.0))) { }
  BoundaryConditions *copy() const { return new BCC(min_diameter()); }
  void map_to_central_box(Cartesian &c) const 
  { 
    bcc_map_to_central_box(&c.x, &c.y, &c.z, side, 1/side);
  }
  void nearest_lattice_point(const Cartesian &c, int &i, int &j, int &k) const
  {
    bcc_nearest_lattice_point(c.x, c.y, c.z, side, &i, &j, &k);
  }
  double volume() const { return 0.5*cube(side); }
  double min_diameter() const { return 0.5*sqrt(3.0)*side; }
  double max_diameter() const { return 0.5*sqrt(5.0)*side; }
  void scale(double s) { side *= s; }
  Tensor lattice_vectors() const { return (0.5*side)*Tensor(1,-1, 1, 1, 1,-1, 1, 1, 1); }
  void write(ostream &s) const { s << "bcc " << min_diameter(); }
private:
  double side; // length of side of containing cube
};

class FCC : public BoundaryConditions
{
public:
  FCC(double b) : BoundaryConditions(fcc), halfside(b/sqrt(2.0)) { }
  BoundaryConditions *copy() const { return new FCC(min_diameter()); }
  void map_to_central_box(Cartesian &c) const
  {
    fcc_map_to_central_box(&c.x, &c.y, &c.z, halfside, 1/halfside);
  }
  void nearest_lattice_point(const Cartesian &c, int &i, int &j, int &k) const
  {
    fcc_nearest_lattice_point(c.x, c.y, c.z, halfside, &i, &j, &k);
  }
  double volume() const { return 0.25*cube(2*halfside); }
  double min_diameter() const { return sqrt(2.0)*halfside; }
  double max_diameter() const { return 2*halfside; }
  void scale(double s) { halfside *= s; }
  Tensor lattice_vectors() const { return halfside*Tensor(1,0,1,1,1,0,0,1,1); }
  void write(ostream &s) const { s << "fcc " << min_diameter(); }
private:
  double halfside; // half of length of containing cube
};

class Triclinic : public BoundaryConditions
{
public:
  Triclinic(double a, double b, double c, double alpha, double beta, double gamma) : BoundaryConditions(triclinic)
  {
    const double cosa = cos(alpha);
    const double cosb = cos(beta);
    const double cosg = cos(gamma);
    const double sing = sin(gamma);
    const double tmp = 1 - sq(cosa) - sq(cosb) - sq(cosg) + 2*cosa*cosb*cosg;
    insist(tmp > 0);
    box_vectors = Tensor(a, b*cosg, c*cosb,
			 0, b*sing, c*(cosa-cosb*cosg)/sing,
			 0,      0, c*sqrt(tmp)/sing);
    inverse_box_vectors = box_vectors.inverse();
  }
  Triclinic(const Tensor &h) : BoundaryConditions(triclinic), box_vectors(h), inverse_box_vectors(h.inverse()) { }
  BoundaryConditions *copy() const { return new Triclinic(a(),b(),c(),alpha(),beta(),gamma()); }
  void map_to_central_box(Cartesian &c) const
  { 
    triclinic_map_to_central_box(&c.x,&c.y,&c.z,
				 (const tensor_t *) &box_vectors,
				 (const tensor_t *) &inverse_box_vectors);
  }
  void nearest_lattice_point(const Cartesian &c, int &i, int &j, int &k) const
  {
    triclinic_nearest_lattice_point(c.x,c.y,c.z,
				    (const tensor_t *) &inverse_box_vectors,
				    &i,&j,&k);
  }
  double volume() const { return fabs(box_vectors.determinant()); }
  double min_diameter() const 
  { 
    const double a01 = box_vectors.col(0).cross(box_vectors.col(1)).magnitude();
    const double a02 = box_vectors.col(0).cross(box_vectors.col(2)).magnitude();
    const double a12 = box_vectors.col(1).cross(box_vectors.col(2)).magnitude();
    return volume()/max(max(a01,a02),a12);
  }
  double max_diameter() const 
  {
    const double d1 = (box_vectors.col(0)+box_vectors.col(1)+box_vectors.col(2)).magnitude();
    const double d2 = (-box_vectors.col(0)+box_vectors.col(1)+box_vectors.col(2)).magnitude();
    const double d3 = (box_vectors.col(0)-box_vectors.col(1)+box_vectors.col(2)).magnitude();
    const double d4 = (box_vectors.col(0)+box_vectors.col(1)-box_vectors.col(2)).magnitude();
    return max(max(d1,d2),max(d3,d4));
  }
  void scale(double s) 
  { 
    box_vectors *= s; 
    inverse_box_vectors *= 1/s;
  }
  Tensor lattice_vectors() const { return box_vectors; }
  double a() const { return box_vectors.col(0).magnitude(); }
  double b() const { return box_vectors.col(1).magnitude(); }
  double c() const { return box_vectors.col(2).magnitude(); }
  double alpha() const 
  { return Angle(box_vectors.col(1), Cartesian(0,0,0), box_vectors.col(2)); }
  double beta() const
  { return Angle(box_vectors.col(0), Cartesian(0,0,0), box_vectors.col(2)); }
  double gamma() const
  { return Angle(box_vectors.col(0), Cartesian(0,0,0), box_vectors.col(1)); }
  void write(ostream &s) const 
  { 
    s << "triclinic " << a() << " " << b() << " " << c() << "  "
      << radians_to_degrees(alpha()) << " "
      << radians_to_degrees(beta()) << " "
      << radians_to_degrees(gamma());
  }
private:
  Tensor box_vectors, inverse_box_vectors;
};

BoundaryConditions::~BoundaryConditions() { }

Tensor BoundaryConditions::reciprocal_lattice_vectors() const
{ return (2*M_PI)*lattice_vectors().inverse().transpose(); }

BoundaryConditions *BoundaryConditions::new_non_periodic() 
{ return new NonPeriodic; }

BoundaryConditions *BoundaryConditions::new_cubic(double boxl) 
{ return new Cubic(boxl); }

BoundaryConditions *BoundaryConditions::new_orthorhombic(Cartesian boxl) 
{ return new Orthorhombic(boxl); }
  
BoundaryConditions *BoundaryConditions::new_bcc(double boxl) 
{ return new BCC(boxl); }

BoundaryConditions *BoundaryConditions::new_fcc(double boxl) 
{ return new FCC(boxl); }

BoundaryConditions *BoundaryConditions::new_triclinic(const Tensor &h) 
{ return new Triclinic(h); }

BoundaryConditions *BoundaryConditions::new_triclinic(double a, double b, double c,
						      double alpha, double beta, double gamma) 
{ 
  return new Triclinic(a,b,c,
		       degrees_to_radians(alpha),
		       degrees_to_radians(beta),
		       degrees_to_radians(gamma)); 
}

BoundaryConditions *BoundaryConditions::new_from_cell(double a, double b, double c,
						      double alpha, double beta, double gamma) 
{ 
  if (are_approximately_equal(alpha,90) && 
      are_approximately_equal(beta,90) && 
      are_approximately_equal(gamma,90))
    if (are_approximately_equal(a,b) && are_approximately_equal(b,c))
      return BoundaryConditions::new_cubic(a);
    else
      return BoundaryConditions::new_orthorhombic(Cartesian(a,b,c));
  else
    return BoundaryConditions::new_triclinic(a,b,c,alpha,beta,gamma);
}

istream & operator>>(istream &s, BoundaryConditions *&bc)
{
  if (bc)
    delete bc;
  bcondtype_t t;
  s >> t;
  double boxl, a, b, c, alpha, beta, gamma;
  Cartesian boxd;
  switch (t) {
  case non_periodic:
    bc = BoundaryConditions::new_non_periodic();
    break;
  case cubic:
    s >> boxl;
    if (!s)
      die("BoundaryConditions: error reading Cubic: expecting box length");
    bc = BoundaryConditions::new_cubic(boxl);
    break;
  case orthorhombic:    
    s >> boxd;
    if (!s)
      die("BoundaryConditions: error reading Orthorhombic: expecting box dimensions");
    bc = BoundaryConditions::new_orthorhombic(boxd);
    break;
  case bcc:
    s >> boxl;
    if (!s)
      die("BoundaryConditions: error reading BCC: expecting box length");
    bc = BoundaryConditions::new_bcc(boxl);
    break;
  case fcc:
    s >> boxl;
    if (!s)
      die("BoundaryConditions: error reading FCC: expecting box length");
    bc = BoundaryConditions::new_fcc(boxl);
    break;
  case triclinic:
    s >> a >> b >> c >> alpha >> beta >> gamma;
    if (!s)
      die("BoundaryConditions: error reading Triclinic: expecting a,b,c,alpha,beta,gamma");
    bc = BoundaryConditions::new_from_cell(a,b,c,alpha,beta,gamma);
    break;
  }
  return s;
}

ostream & operator<<(ostream &s, const BoundaryConditions *bc)
{
  bc->write(s);
  return s;
}
