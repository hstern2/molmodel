#ifndef ATOM_H
#define ATOM_H

#include "str.hpp"
#include "hset.hpp"
#include "lst.hpp"
#include "pdb.hpp"
#include "coord.hpp"
#include "fns.h"
#include "random.h"
#include "cmplx.hpp"

typedef HSet<Str> HSetStr;

struct Element
{
  int atomic_number; // io
  double mass, covalent_radius, vdw_radius; // io

  classIO(Element);
};

ostream & operator<<(ostream &, const Element *);

class Atom
{

public:

  Atom();
  Atom(Str sym);
  Atom(Str sym, Str t);
  Atom(const PDBAtomRecord &, const Cartesian &);

  Cartesian position, velocity; // out
  Str symbol, type; // out
  LstInt neighbors; // out
  LstStr bond_types; // out
  PDBAtomRecord pdb; // out
  HSetStr property; // out
  const Element* element; // out

  bool is_hydrogen() const;
  bool is_dummy() const;
  const double &mass() const { return element->mass; }
  const double &covalent_radius() const { return element->covalent_radius; }
  const double &vdw_radius() const { return element->vdw_radius; }
  const int &atomic_number() const { return element->atomic_number; }
  Complex nuclear_structure_factor(const Cartesian &k) const;

  void copy_from(const Atom &);
  void write_as_xyz(ostream &s) const;
  void write_as_mol2(ostream &s) const;
  void create_from_tinker(istream &s);
  Cartesian linear_momentum() const { return mass()*velocity; }
  double kinetic_energy() const { return 0.5*mass()*velocity.sq(); }
  const char *bond_type_for_neighbor(int j) const;
  static void finalize();
  static void set_mass(Str sym, double mass);
  classIO(Atom);
};

#endif /* ATOM_H */
