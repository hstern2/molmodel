#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "msys.hpp"
#include "htab.hpp"
#include "cvec.hpp"

class Cnstr;
class ConstraintGroup;
class Potential;

struct CnstrPtr
{
  Cnstr * const p;
  CnstrPtr(Cnstr *p_) : p(p_) { }
  friend unsigned hash_val(const CnstrPtr &key);
  friend bool operator==(const CnstrPtr &a, const CnstrPtr &b);
};

class Constraint
{

public:

  HTab<double,Str2> bond_length; // in
  HTab<double,Str3> angle; // in
  HTab<double,Int4> dihedral; // in
  HTab<double,Str4> sp2_proton; // in
  HTab<double,Int2> distance; // in
  double tolerance; // in
  int max_iterations, verbose; // in

  Constraint();
  ~Constraint();
  void init(MSys *);
  void reset();
  void add_constraints(const MSys &);
  int how_many() const;
  int how_many_on_atom(int) const;
  void constrain_positions(double dt);
  void constrain_velocities();
  void write() const;
  void add_to_energy_and_forces(double &, Cartesian *);
  friend istream & operator>>(istream &s, Constraint *&p)
  {
    if (p)
      delete p;
    return s >> *(p = new Constraint);
  }
  classIO(Constraint);
  
private:

  void cleanup(), add_constraint(Cnstr *);
  MSys *msys;
  HSet<CnstrPtr> cnstr;
  Lst<ConstraintGroup *> group;
  Vec<Lst<Cnstr *> > cnstr_on_atom;
  Vec<ConstraintGroup *> group_on_atom;
  CVec r0;
  int ncnstr;
    
};

#endif /* CONSTRAINT_H */

