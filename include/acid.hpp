#ifndef ACID_H
#define ACID_H

#include "atom.hpp"
#include "vec.hpp"
#include "htab.hpp"

struct ProtonationState
{
  Vec<Boolean> is_protonated; // io
  int degeneracy; // io
  int number_of_protons() const; /* number of protons for this particular state */
  ProtonationState() : degeneracy(1) { }
  classIO(ProtonationState);
};

class AcidicGroupDesc;
struct MSys;
class Atom;
class Constraint;

struct AcidicGroup
{
  const AcidicGroupDesc *desc;
  Vec<int> index; // out
  Vec<Vec<Atom> > atom; // out
  int state; // out
  int number_of_protons() const;
  AcidicGroup();
  AcidicGroup(const AcidicGroupDesc *, int state, const Vec<int> index, const Atom *);
  AcidicGroup(const AcidicGroup &);
  AcidicGroup &operator=(const AcidicGroup &);
  void reindex(const int *iold, const int *inew);
  classIO(AcidicGroup);
};

class AcidicGroupDesc
{
public:
  Str pattern; // io
  Vec<int> protonation_site; // io    indices of atom to which protons are attached
  Vec<ProtonationState> state; // io
  classIO(AcidicGroupDesc);
  void check() const; // in
  void find_groups(Lst<AcidicGroup> &group, Vec<Atom> atom, Lst<Atom> &newatom, 
		   Lst<int> &bonded_to, const HTab<HSet<int> > &atoms_with_property);
  int number_of_states() const { return state.size(); }
};

#endif /* ACID_H */
