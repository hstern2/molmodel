#include "acid.hpp"
#include "msys.hpp"
#include "apattern.hpp"
#include "geograd.hpp"
#include "random.h"
#include "constraint.hpp"

static int bonded_proton_not_in_path(const Atom *at, int i, const HSet<int> &p)
{
  LstItr<int> n;
  for (n.init(at[i].neighbors); n.ok(); n.next())
    if ((at[n()].is_dummy() || at[n()].is_hydrogen()) && !p.exists(n()))
      return n();
  return -1;
}

int ProtonationState::number_of_protons() const
{
  int n = 0;
  for (int i = 0; i < is_protonated.size(); i++)
    if (is_protonated[i])
      n++;
  return n;
}

void AcidicGroupDesc::check() const
{
  int i;
  if (state.size() == 0)
    die("AcidicGroupDesc::check: no protonation states for acid group %s", (const char *) pattern);
  int min = state[0].number_of_protons();
  int max = state[0].number_of_protons();
  for (i = 1; i < state.size(); i++) {
    if (state[i].number_of_protons() < min)
      min = state[i].number_of_protons();
    if (state[i].number_of_protons() > max)
      max = state[i].number_of_protons();
  }
  for (i = min; i <= max; i++) {
    int j;
    for (j = 0; j < state.size(); j++)
      if (i == state[j].number_of_protons())
	break;
    if (j == state.size())
      die("AcidicGroupDesc::check: states with %d protons and %d protons, "
	  "but not %d protons, for acidic group %s", min, max, i, (const char *) pattern);
  }
}

void AcidicGroupDesc::find_groups(Lst<AcidicGroup> &group, Vec<Atom> atom, 
				  Lst<Atom> &newatom, Lst<int> &bonded_to,
				  const HTab<HSet<int> > &atoms_with_property)
{
  struct AcidPattern : public AtomPattern
  {
    AcidicGroupDesc &a;
    Lst<AcidicGroup> &group;
    Lst<Atom> &newatom;
    Lst<int> &bonded_to;
    HSet<HSet<int> > paths;
    AcidPattern(AcidicGroupDesc &a_, Lst<AcidicGroup> &g_, Lst<Atom> &newatom_, Lst<int> &bonded_to_) : 
      AtomPattern(a_.pattern), a(a_), group(g_), newatom(newatom_), bonded_to(bonded_to_)
    { 
      const int nprotons = a.protonation_site.size();
      for (int i = 0; i < a.state.size(); i++) {
	const ProtonationState &p = a.state[i];
	if (p.is_protonated.size() != nprotons)
	  die("AcidPattern: size of list is_protonated (%d) "
	      "must be equal to number of protonated atoms (%d)", 
	      p.is_protonated.size(), nprotons);
      }
    }
    void succeed(const int *path, Atom *at)
    { 
      HSet<int> p(ntest);
      p.add(ntest,path);
      if (paths.exists(p))
	return; /* seen a permutation of this path before */
      const int np = a.protonation_site.size(); /* number of labile protons */
      Vec<int> in(ntest+np);
      for (int is = 0; is < a.state.size(); is++) {
	Lst<Atom> newatom_tmp;
	Lst<int> bonded_to_tmp;
	const ProtonationState &ps = a.state[is];
	int i;
	for (i = 0; i < np; i++) {
	  const int ia = path[a.protonation_site[i]];
	  const int k = bonded_proton_not_in_path(at,ia,p);
	  if (k >= 0) {
	    if (ps.is_protonated[i] == at[k].is_dummy())
	      break; /* protonation state doesn't match */
	    in[ntest+i] = k; /* so far, so good */
	  } else {
	    if (ps.is_protonated[i])
	      break; /* protonation state doesn't match */
	    /* need to add a new atom */
	    Atom atmp;
	    atmp.position = at[ia].position;
	    atmp.position.x += GaussianRandom();
	    atmp.position.y += GaussianRandom();
	    atmp.position.z += GaussianRandom();
	    newatom_tmp.add(atmp);
	    bonded_to_tmp.add(ia);
	    in[ntest+i] = -(newatom.size() + newatom_tmp.size());
	  }
	}
	if (i < np)
	  continue; /* protonation state doesn't match */
	/* This protonation state matches */
	for (i = 0; i < ntest; i++)
	  in[i] = path[i];
	group.add(AcidicGroup(&a,is,in,at));
	paths.add(p);
	newatom.add(newatom_tmp);
	bonded_to.add(bonded_to_tmp);
	return;
      }
      die("AcidPattern: cannot find matching protonation state for pattern %s", (const char *) a.pattern);
    }
  };
  check();
  AcidPattern a(*this,group,newatom,bonded_to);
  a.run(atom,atoms_with_property);
}

AcidicGroup::AcidicGroup() : desc(0), state(-1) { }

AcidicGroup::AcidicGroup(const AcidicGroupDesc *d, int s, const Vec<int> i, const Atom *a)
  : desc(d), index(i), state(s)
{ 
  const int nat = index.size();
  const int nstate = desc->number_of_states();
  const int np = desc->protonation_site.size();
  atom.resize(nstate);
  for (s = 0; s < nstate; s++) {
    atom[s].resize(nat);
    for (int j = 0; j < nat; j++) {
      const int ip = j - nat + np;
      if (ip >= 0) {
	insist(ip < np);
	atom[s][j] = desc->state[s].is_protonated[ip] ?
	  Atom("H") : Atom();
      } else {
	atom[s][j].copy_from(a[index[j]]);
      }
    }
  }
}

AcidicGroup::AcidicGroup(const AcidicGroup &a) 
  : desc(a.desc), index(a.index.copy()), 
    atom(a.atom.copy()), state(a.state)
{ }

AcidicGroup & AcidicGroup::operator=(const AcidicGroup &a)
{ 
  desc = a.desc;
  index = a.index.copy();
  atom = a.atom.copy();
  state = a.state;
  return *this;
}

void AcidicGroup::reindex(const int *iold, const int *inew)
{
  for (int i = 0; i < index.size(); i++) {
    int &p = index[i];
    if (p >= 0)
      p = iold[p];
    else
      p = inew[-(p+1)];
  }
}

int AcidicGroup::number_of_protons() const
{
  return desc->state[state].number_of_protons();
}

