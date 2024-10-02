#ifndef INTRA_H
#define INTRA_H

#include "out.hpp"
#include "ntuple.hpp"
#include "htab.hpp"
#include "pot.hpp"
#include "str.hpp"
#include "hset.hpp"
#include "cvec.hpp"
#include "fns.h"

/* u = k(r-r0)^2 */
struct StretchParam
{ 
  double k, r0; // io
  StretchParam() : k(UNDEF_VAL), r0(UNDEF_VAL) { }
  classIO(StretchParam);
};

/* u = k(theta-theta0)^2 */
struct BendParam
{ 
  double k, theta0; // io
  BendParam() : k(UNDEF_VAL), theta0(UNDEF_VAL) { }
  classIO(BendParam);
};

/* u = sum_n v_n[1 +/- cos(n phi)] */
struct TorsionParam
{
  double v, sign;
  int n;
  TorsionParam() : v(UNDEF_VAL), sign(UNDEF_VAL) { }
  friend istream & operator>>(istream &s, TorsionParam &p)
  { 
    s >> p.v >> p.sign >> p.n; 
    assert(p.n > 0);
    return s;
  }
  friend ostream & operator<<(ostream &s, const TorsionParam &p)
  { return s << p.v << " " << p.sign << " " << p.n << " "; }
};

/***
 * Improper torsions are given by the same functional form
 * as ordinary torsions.  Following the AMBER convention, 
 * the third atom of four is the trigonal planar atom 
 * on which the improper torsion is centered.
 ***/

typedef Lst<TorsionParam> LstTorsionParam;

/***
 * Dihedral restraint: u = k/2 (phi-phi0)^2
 * phi0 in degrees, k in kcal/(mol degrees^2)
 ***/
struct DihedralRestraint
{
  Int4 atoms; // io
  double phi0, k; // io
  Boolean phi0_from_initial_geometry; // io
  void add_to_energy_and_forces(const Atom *, double &u, Cartesian *f) const;
  void init(int n, const Atom *);
  DihedralRestraint() : phi0(0), k(-1), phi0_from_initial_geometry(false) { }
  classIO(DihedralRestraint);
};

struct Stretch
{
  const int i, j;
  StretchParam p, p1, p2;
  double energy(const Atom *) const;
  void add_to_energy_and_forces(const Atom *, double &u, Cartesian *f) const;
  void set_lambda(double lambda);
  Stretch(int i_, int j_) : i(i_), j(j_) { }
  friend ostream & operator<<(ostream &s, const Stretch &c)
  { 
    s << c.i << " " << c.j << "  " << c.p.k << " " << c.p.r0;
    return s;
  }
};

struct Bend
{
  const int i, j, k;
  BendParam p, p1, p2;
  double energy(const Atom *) const;
  void add_to_energy_and_forces(const Atom *, double &u, Cartesian *f) const;
  void set_lambda(double lambda);
  Bend(int i_, int j_, int k_) : i(i_), j(j_), k(k_) { }
  friend ostream & operator<<(ostream &s, const Bend &c)
  { 
    s << c.i << " " << c.j << " " << c.k << "  " << c.p.k << " " << c.p.theta0;
    return s;
  }
};

struct Torsion
{
  const int i, j, k, l;
  LstTorsionParam p, p1, p2;
  double energy(const Atom *) const;
  void add_to_energy_and_forces(const Atom *, double &u, Cartesian *f) const;
  void set_lambda(double lambda);
  Torsion(int i_, int j_, int k_, int l_)
    : i(i_), j(j_), k(k_), l(l_) { }
  friend ostream & operator<<(ostream &s, const Torsion &c)
  { 
    s << c.i << " " << c.j << " " << c.k << " " << c.l << "   " << c.p;
    return s;
  }
};

class Atom;

class IntramolecularPotential : public Potential
{

public:

  HTab<StretchParam,Str2> stretch; // io
  HTab<BendParam,Str3> bend; // io
  HTab<LstTorsionParam,Str4> torsion, improper; // io
  HTab<Str> general_type; // io
  Lst<DihedralRestraint> dihedral_restraint; // io
  int verbose; // io
  Boolean include_stretch, include_bend, include_torsion, include_improper; // io
  
  IntramolecularPotential() 
    : verbose(1), include_stretch(true), include_bend(true),
      include_torsion(true), include_improper(true)
  { }

  void init(MSys *);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  void write_details() const;

  void types_changed();
  void types_changed(const MSys &perturbed);
  void set_lambda(double lambda);

  const StretchParam *stretch_for(const Str, const Str) const;
  const BendParam *bend_for(const Str, const Str, const Str) const;
  const LstTorsionParam *torsion_for(const Str, const Str, const Str, const Str) const;
  const LstTorsionParam *improper_for(const Str, const Str, const Str, const Str) const;

private:

  double ustretch, ubend, utorsion, uimproper, udihedral_restraint;
  bool any_stretch, any_bend, any_torsion, any_improper, any_dihedral_restraint;
  CVec ftmp;
  Lst<Stretch> stretch_list; // out
  Lst<Bend> bend_list; // out
  Lst<Torsion> torsion_list, improper_list; // out

  Str code_for_type(const Str) const;

  void update_stretch(Stretch &, HTab<Int2,Str2> &undef) const;
  void update_bend(Bend &, HTab<Int3,Str3> &undef) const;
  void update_torsion(Torsion &, HTab<Int4,Str4> &undef) const;
  void update_improper(Torsion &, HTab<Int4,Str4> &undef) const;

  void update_stretch(Stretch &, HTab<Int2,Str2> &undef, const MSys &perturbed) const;
  void update_bend(Bend &, HTab<Int3,Str3> &undef, const MSys &perturbed) const;
  void update_torsion(Torsion &, HTab<Int4,Str4> &undef, const MSys &perturbed) const;
  void update_improper(Torsion &, HTab<Int4,Str4> &undef, const MSys &perturbed) const;

  classIO(IntramolecularPotential);

};

#endif /* INTRA_H */
