#include "intra.hpp"
#include "str.hpp"
#include "msys.hpp"
#include "geograd.hpp"
#include "timing.h"
#include "units.h"

static bool matches(const LstTorsionParam &a,
		   const LstTorsionParam &b)
{
  if (a.size() != b.size())
    return 0;
  LstItr<TorsionParam> pi, pj;
  for (pi.init(a), pj.init(b); pi.ok(); pi.next(), pj.next())
    if (!are_approximately_equal(pi().v,pj().v) || 
	!are_approximately_equal(pi().sign,pj().sign) || 
	pi().n != pj().n)
      return false;
  return true;
}

const StretchParam *IntramolecularPotential::stretch_for(const Str si, const Str sj) const
{
  const Str2 pstr(si,sj);
  const Str2 qstr(sj,si);
  const StretchParam *p = stretch.get(pstr);
  const StretchParam *q = stretch.get(qstr);
  if (p) {
    if (q && !(are_approximately_equal(p->k,q->k) && are_approximately_equal(p->r0,q->r0)))
      Out() << "*** conflicting stretches:\n"
	    << "*** " << pstr << " " << *p << "\n"
	    << "*** " << qstr << " " << *q << "\n"
	    << flush;
    return p;
  }
  return q;
}

const BendParam *IntramolecularPotential::bend_for(const Str si, const Str sj, const Str sk) const
{
  const Str3 pstr(si,sj,sk);
  const Str3 qstr(sk,sj,si);
  const BendParam *p = bend.get(pstr);
  const BendParam *q = bend.get(qstr);
  if (p) {
    if (q && !(are_approximately_equal(p->k,q->k) && are_approximately_equal(p->theta0,q->theta0)))
      Out() << "*** conflicting bends:\n"
	    << "*** " << pstr << " " << *p << "\n"
	    << "*** " << qstr << " " << *q << "\n"
	    << flush;
    return p;
  }
  return q;
}

const LstTorsionParam *IntramolecularPotential::torsion_for(const Str si, const Str sj, const Str sk, const Str sl) const
{
  const Str4 pstr(si,sj,sk,sl);
  const Str4 qstr(sl,sk,sj,si);
  const LstTorsionParam *p = torsion.get(pstr);
  const LstTorsionParam *q = torsion.get(qstr);
  if (p) {
    if (q && !matches(*p,*q))
      Out() << "*** conflicting torsions:\n"
	    << "*** " << pstr << " " << *p << "\n"
	    << "*** " << qstr << " " << *q << "\n"
	    << flush;
    return p;
  }
  return q;
}

const LstTorsionParam *IntramolecularPotential::improper_for(const Str si, const Str sj, const Str sk, const Str sl) const
{
  return improper.get(Str4(si,sj,sk,sl));
}

Str IntramolecularPotential::code_for_type(const Str t) const
{
  Str *k = general_type.get(t);
  return k ? *k : t;
}

void IntramolecularPotential::update_stretch(Stretch &s, HTab<Int2,Str2> &undef) const
{
  const Atom &ai = msys->atom[s.i];
  const Atom &aj = msys->atom[s.j];
  const Str si = code_for_type(ai.type);
  const Str sj = code_for_type(aj.type);
  const StretchParam *p = stretch_for(si,sj);
  if (p) {
    s.p = *p;
    return;
  }
  s.p.k = s.p.r0 = 0;
  if (!ai.is_dummy() && !aj.is_dummy())
    undef[Str2(si,sj)] = Int2(s.i,s.j);
}

void IntramolecularPotential::update_stretch(Stretch &s, HTab<Int2,Str2> &undef,
					     const MSys &perturbed) const
{
  {
    const Atom &ai = msys->atom[s.i];
    const Atom &aj = msys->atom[s.j];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const StretchParam *p = stretch_for(si,sj);
    if (p) {
      s.p1 = *p;
    } else {
      s.p1.k = s.p1.r0 = 0;
      if (!ai.is_dummy() && !aj.is_dummy())
	undef[Str2(si,sj)] = Int2(s.i,s.j);
    }
  }
  {
    const Atom &ai = perturbed.atom[s.i];
    const Atom &aj = perturbed.atom[s.j];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const StretchParam *p = stretch_for(si,sj);
    if (p) {
      s.p2 = *p;
    } else {
      s.p2.k = s.p2.r0 = 0;
      if (!ai.is_dummy() && !aj.is_dummy())
	undef[Str2(si,sj)] = Int2(s.i,s.j);
    }
  }
}

void IntramolecularPotential::update_bend(Bend &b, HTab<Int3,Str3> &undef) const
{
  const Atom &ai = msys->atom[b.i];
  const Atom &aj = msys->atom[b.j];
  const Atom &ak = msys->atom[b.k];
  const Str si = code_for_type(ai.type);
  const Str sj = code_for_type(aj.type);
  const Str sk = code_for_type(ak.type);  
  const BendParam *p = bend_for(si,sj,sk);
  if (p) {
    b.p = *p;
  } else {
    b.p.k = b.p.theta0 = 0;
    if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy())
      undef[Str3(si,sj,sk)] = Int3(b.i,b.j,b.k);
  }
}

void IntramolecularPotential::update_bend(Bend &b, HTab<Int3,Str3> &undef,
					  const MSys &perturbed) const
{
  {
    const Atom &ai = msys->atom[b.i];
    const Atom &aj = msys->atom[b.j];
    const Atom &ak = msys->atom[b.k];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);  
    const BendParam *p = bend_for(si,sj,sk);
    if (p) {
      b.p1 = *p;
    } else {
      b.p1.k = b.p1.theta0 = 0;
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy())
	undef[Str3(si,sj,sk)] = Int3(b.i,b.j,b.k);
    }
  }
  {
    const Atom &ai = perturbed.atom[b.i];
    const Atom &aj = perturbed.atom[b.j];
    const Atom &ak = perturbed.atom[b.k];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);  
    const BendParam *p = bend_for(si,sj,sk);
    if (p) {
      b.p2 = *p;
    } else {
      b.p2.k = b.p2.theta0 = 0;
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy())
	undef[Str3(si,sj,sk)] = Int3(b.i,b.j,b.k);
    }
  }
}

void IntramolecularPotential::update_torsion(Torsion &t, HTab<Int4,Str4> &undef) const
{
  const Atom &ai = msys->atom[t.i];
  const Atom &aj = msys->atom[t.j];
  const Atom &ak = msys->atom[t.k];
  const Atom &al = msys->atom[t.l];
  const Str si = code_for_type(ai.type);
  const Str sj = code_for_type(aj.type);
  const Str sk = code_for_type(ak.type);
  const Str sl = code_for_type(al.type);
  const LstTorsionParam *p = torsion_for(si,sj,sk,sl);
  if (p) {
    t.p = *p;
  } else {
    t.p.remove_all();
    if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
      undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
  }
}

void IntramolecularPotential::update_torsion(Torsion &t, HTab<Int4,Str4> &undef,
					     const MSys &perturbed) const
{
  t.p.remove_all();
  {
    const Atom &ai = msys->atom[t.i];
    const Atom &aj = msys->atom[t.j];
    const Atom &ak = msys->atom[t.k];
    const Atom &al = msys->atom[t.l];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);
    const Str sl = code_for_type(al.type);
    const LstTorsionParam *p = torsion_for(si,sj,sk,sl);
    if (p) {
      t.p.add(t.p1 = *p);
    } else {
      t.p1.remove_all();
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
	undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
    }
  }
  {
    const Atom &ai = perturbed.atom[t.i];
    const Atom &aj = perturbed.atom[t.j];
    const Atom &ak = perturbed.atom[t.k];
    const Atom &al = perturbed.atom[t.l];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);
    const Str sl = code_for_type(al.type);
    const LstTorsionParam *p = torsion_for(si,sj,sk,sl);
    if (p) {
      t.p.add(t.p2 = *p);
    } else {
      t.p2.remove_all();
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
	undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
    }
  }
}

void IntramolecularPotential::update_improper(Torsion &t, HTab<Int4,Str4> &undef) const
{
  const Atom *a = msys->atom;
  const Atom &ai = a[t.i];
  const Atom &aj = a[t.j];
  const Atom &ak = a[t.k];
  const Atom &al = a[t.l];
  const Str si = code_for_type(ai.type);
  const Str sj = code_for_type(aj.type);
  const Str sk = code_for_type(ak.type);
  const Str sl = code_for_type(al.type);
  const int aeqb = !strcmp(si,sj);
  const int beqc = !strcmp(sj,sl);
  const int aeqc = !strcmp(si,sl);
  const double scale = aeqb && beqc ? 1.0/6.0 : (aeqb || beqc || aeqc ? 0.5 : 1);
  const LstTorsionParam *p = improper_for(si,sj,sk,sl);
  if (p) {
    t.p = *p;
    LstItrMod<TorsionParam> tp;
    for (tp.init(t.p); tp.ok(); tp.next())
      tp().v *= scale;
  } else {
    t.p.remove_all();
    if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
      undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
  }
}

void IntramolecularPotential::update_improper(Torsion &t, HTab<Int4,Str4> &undef,
					      const MSys &perturbed) const
{
  t.p.remove_all();
  t.p1.remove_all();
  t.p2.remove_all();
  {
    const Atom *a = msys->atom;
    const Atom &ai = a[t.i];
    const Atom &aj = a[t.j];
    const Atom &ak = a[t.k];
    const Atom &al = a[t.l];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);
    const Str sl = code_for_type(al.type);
    const int aeqb = !strcmp(si,sj);
    const int beqc = !strcmp(sj,sl);
    const int aeqc = !strcmp(si,sl);
    const double scale = aeqb && beqc ? 1.0/6.0 : (aeqb || beqc || aeqc ? 0.5 : 1);
    const LstTorsionParam *p = improper_for(si,sj,sk,sl);
    if (p) {
      t.p1 = *p;
      LstItrMod<TorsionParam> tp;
      for (tp.init(t.p1); tp.ok(); tp.next())
	tp().v *= scale;
      t.p.add(t.p1);
    } else {
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
	undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
    }
  }
  {
    const Atom *a = perturbed.atom;
    const Atom &ai = a[t.i];
    const Atom &aj = a[t.j];
    const Atom &ak = a[t.k];
    const Atom &al = a[t.l];
    const Str si = code_for_type(ai.type);
    const Str sj = code_for_type(aj.type);
    const Str sk = code_for_type(ak.type);
    const Str sl = code_for_type(al.type);
    const int aeqb = !strcmp(si,sj);
    const int beqc = !strcmp(sj,sl);
    const int aeqc = !strcmp(si,sl);
    const double scale = aeqb && beqc ? 1.0/6.0 : (aeqb || beqc || aeqc ? 0.5 : 1);
    const LstTorsionParam *p = improper_for(si,sj,sk,sl);
    if (p) {
      t.p2 = *p;
      LstItrMod<TorsionParam> tp;
      for (tp.init(t.p2); tp.ok(); tp.next())
	tp().v *= scale;
      t.p.add(t.p2);
    } else {
      if (!ai.is_dummy() && !aj.is_dummy() && !ak.is_dummy() && !al.is_dummy())
	undef[Str4(si,sj,sk,sl)] = Int4(t.i,t.j,t.k,t.l);
    }
  }
}

void IntramolecularPotential::types_changed()
{
  TIMESTART("Intra::types_changed");
  HTab<Int2,Str2> undef_stretch;
  HTab<Int3,Str3> undef_bend;
  HTab<Int4,Str4> undef_torsion, undef_improper;
  LstItrMod<Stretch> s;
  for (s.init(stretch_list); s.ok(); s.next())
    update_stretch(s(),undef_stretch);
  LstItrMod<Bend> b;
  for (b.init(bend_list); b.ok(); b.next())
    update_bend(b(),undef_bend);
  LstItrMod<Torsion> t;
  for (t.init(torsion_list); t.ok(); t.next())
    update_torsion(t(),undef_torsion);
  for (t.init(improper_list); t.ok(); t.next())
    update_improper(t(),undef_improper);
  const Atom *at = msys->atom;
  if (undef_stretch.size() > 0) {
    Out() << "*** undefined stretch parameters:\n";
    OutPrintf("%5s %3s %8s %8s    %5s %3s %8s %8s\n", 
	      "i", "sym", "type", "general", 
	      "j", "sym", "type", "general");
    HTabIterator<Int2,Str2> u;
    for (u.init(undef_stretch); u.ok(); u.next()) {
      const int i = u.val().a;
      const int j = u.val().b;
      OutPrintf("%5d %3s %8s %8s    %5d %3s %8s %8s\n", 
		i, (const char *) at[i].symbol, (const char *) at[i].type, (const char *) u.key().a,
		j, (const char *) at[j].symbol, (const char *) at[j].type, (const char *) u.key().b);
    }
    Out() << "\n";
  }
  if (undef_bend.size() > 0) {
    Out() << "*** undefined bend parameters:\n";
    OutPrintf("%5s %3s %8s %8s    %5s %3s %8s %8s    %5s %3s %8s %8s\n", 
	      "i", "sym", "type", "general", 
	      "j", "sym", "type", "general", 
	      "k", "sym", "type", "general");
    HTabIterator<Int3,Str3> u;
    for (u.init(undef_bend); u.ok(); u.next()) {
      const int i = u.val().a;
      const int j = u.val().b;
      const int k = u.val().c;
      OutPrintf("%5d %3s %8s %8s    %5d %3s %8s %8s    %5d %3s %8s %8s\n", 
		i, (const char *) at[i].symbol, (const char *) at[i].type, (const char *) u.key().a,
		j, (const char *) at[j].symbol, (const char *) at[j].type, (const char *) u.key().b,
		k, (const char *) at[k].symbol, (const char *) at[k].type, (const char *) u.key().c);
    }
    Out() << "\n";
  }
  if (undef_torsion.size() > 0)
    Out() << "*** no torsion for " << undef_torsion << "\n";
  if (verbose >= 2 && undef_improper.size() > 0)
    Out() << "*** no improper torsion for " << undef_improper << "\n";
  if (verbose >= 2)
    write_details();
  if (verbose >= 3)
    Out() << *this << flush;
  Out() << flush;
  TIMESTOP("Intra::types_changed");
}

void IntramolecularPotential::types_changed(const MSys &perturbed)
{
  TIMESTART("Intra::types_changed_perturbed");
  insist(msys->has_same_topology_as(perturbed));
  HTab<Int2,Str2> undef_stretch;
  HTab<Int3,Str3> undef_bend;
  HTab<Int4,Str4> undef_torsion, undef_improper;
  LstItrMod<Stretch> s;
  for (s.init(stretch_list); s.ok(); s.next())
    update_stretch(s(),undef_stretch,perturbed);
  LstItrMod<Bend> b;
  for (b.init(bend_list); b.ok(); b.next())
    update_bend(b(),undef_bend,perturbed);
  LstItrMod<Torsion> t;
  for (t.init(torsion_list); t.ok(); t.next())
    update_torsion(t(),undef_torsion,perturbed);
  for (t.init(improper_list); t.ok(); t.next())
    update_improper(t(),undef_improper,perturbed);
  if (undef_stretch.size() > 0)
    Out() << "*** no stretch for " << undef_stretch << "\n";
  if (undef_bend.size() > 0)
    Out() << "*** no bend for " << undef_bend << "\n";
  if (undef_torsion.size() > 0)
    Out() << "*** no torsion for " << undef_torsion << "\n";
  if (verbose >= 2 && undef_improper.size() > 0)
    Out() << "*** no improper torsion for " << undef_improper << "\n";
  if (verbose >= 2)
    write_details();
  if (verbose >= 3)
    Out() << *this << flush;
  Out() << flush;
  TIMESTOP("Intra::types_changed_perturbed");
}

void IntramolecularPotential::set_lambda(double lambda)
{
  TIMESTART("Intra::set_lambda");
  LstItrMod<Stretch> s;
  if (include_stretch)
    for (s.init(stretch_list); s.ok(); s.next())
      s().set_lambda(lambda);
  LstItrMod<Bend> b;
  if (include_bend)
    for (b.init(bend_list); b.ok(); b.next())
      b().set_lambda(lambda);
  LstItrMod<Torsion> t;
  if (include_torsion)
    for (t.init(torsion_list); t.ok(); t.next())
      t().set_lambda(lambda);
  if (include_improper)
    for (t.init(improper_list); t.ok(); t.next())
      t().set_lambda(lambda);
  TIMESTOP("Intra::set_lambda");
}

void IntramolecularPotential::init(MSys *m)
{
  Potential::init(m);
  ustretch = ubend = utorsion = uimproper = udihedral_restraint = 0;
  stretch_list.remove_all();
  bend_list.remove_all();
  torsion_list.remove_all();
  improper_list.remove_all();
  any_stretch = stretch.size() > 0;
  any_bend = bend.size() > 0;
  any_torsion = torsion.size() > 0;
  any_improper = improper.size() > 0;
  any_dihedral_restraint = dihedral_restraint.size() > 0;
  if (!any_stretch && !any_bend && !any_torsion && !any_improper && !any_dihedral_restraint)
    return;
  const Atom *atom = msys->atom;
  const int natom = msys->atom.size();
  ftmp.resize(natom);
  int i;
  for (i = 0; i < natom; i++) {
    const Atom &ai = atom[i];
    if (any_stretch) {
      LstItr<int> n;
      for (n.init(ai.neighbors); n.ok(); n.next())
	if (n() > i)
	  stretch_list.add(Stretch(i,n()));
    }
    if (any_bend) {
      LstPairItr<int> n;
      for (n.init(ai.neighbors); n.ok(); n.next()) {
	bend_list.add(Bend(n.i(),i,n.j()));
      }
    }
    if (any_improper && ai.neighbors.size() == 3) {
      LstItr<int> n;
      n.init(ai.neighbors);
      const int n1 = n();
      n.next();
      const int n2 = n();
      n.next();
      const int n3 = n();
      improper_list.add(Torsion(n1, n2, i, n3));
      improper_list.add(Torsion(n1, n3, i, n2));
      improper_list.add(Torsion(n2, n1, i, n3));
      improper_list.add(Torsion(n2, n3, i, n1));
      improper_list.add(Torsion(n3, n1, i, n2));
      improper_list.add(Torsion(n3, n2, i, n1));
    }
  }
  if (any_torsion) {
    Lst<Int2> bonds;
    LstItr<int> ni, nl;
    for (i = 0; i < natom; i++)
      for (ni.init(atom[i].neighbors); ni.ok(); ni.next())
	if (i < ni())
	  bonds.add(Int2(i,ni()));
    LstItr<Int2> b;
    for (b.init(bonds); b.ok(); b.next()) {
      const int j = b().a, k = b().b;
      for (ni.init(atom[j].neighbors); ni.ok(); ni.next())
	if (ni() != k)
	  for (nl.init(atom[k].neighbors); nl.ok(); nl.next())
	    if (nl() != j && nl() != ni())
	      torsion_list.add(Torsion(ni(),j,k,nl()));
    }
  }
  LstItrMod<DihedralRestraint> dh;
  for (dh.init(dihedral_restraint); dh.ok(); dh.next())
    dh().init(msys->atom.size(), msys->atom);
  types_changed();
}

void Stretch::add_to_energy_and_forces(const Atom *a, double &u, Cartesian *f) const
{
  const Cartesian rij = a[i].position - a[j].position;
  const double rmag = rij.magnitude();
  const double rr0 = rmag - p.r0;
  u += p.k*sq(rr0);
  const double du = 2*p.k*rr0;
  const double du_r = du/rmag;
  Cartesian fij(rij);
  fij *= -du_r;
  f[i] += fij;
  f[j] -= fij;
}

double Stretch::energy(const Atom *a) const
{
  return p.k*sq((a[i].position - a[j].position).magnitude() - p.r0);
}

void Stretch::set_lambda(double lambda)
{
  p.k = (1-lambda)*p1.k + lambda*p2.k;
  p.r0 = (1-lambda)*p1.r0 + lambda*p2.r0;
}

void Bend::add_to_energy_and_forces(const Atom *a, double &u, Cartesian *f) const
{
  Cartesian g1, g2, g3;
  const double tt0 = AngleGradient(a[i].position,a[j].position,a[k].position,g1,g2,g3)
    - degrees_to_radians(p.theta0);
  u += p.k*sq(tt0);
  const double du = 2*p.k*tt0;
  f[i] -= du*g1;
  f[j] -= du*g2;
  f[k] -= du*g3;
}

double Bend::energy(const Atom *a) const
{
  return p.k*sq(Angle(a[i].position,a[j].position,a[k].position) - degrees_to_radians(p.theta0));
}

void Bend::set_lambda(double lambda)
{
  p.k = (1-lambda)*p1.k + lambda*p2.k;
  p.theta0 = (1-lambda)*p1.theta0 + lambda*p2.theta0;
}

void Torsion::add_to_energy_and_forces(const Atom *a, double &u, Cartesian *f) const
{
  if (p.size() == 0)
    return;
  Cartesian g1, g2, g3, g4;
  const double phi = DihedralGradient(a[i].position,a[j].position,a[k].position,a[l].position,
				      g1, g2, g3, g4);
  double du = 0;
  LstItr<TorsionParam> ip;
  for (ip.init(p); ip.ok(); ip.next()) {
    const double a = ip().n*phi;
    const double v = ip().v;
    u += v*(1+ip().sign*cos(a));
    du += v*(-ip().sign*ip().n*sin(a));
  }
  f[i] -= du*g1;
  f[j] -= du*g2;
  f[k] -= du*g3;
  f[l] -= du*g4;
}

double Torsion::energy(const Atom *a) const
{
  if (p.size() == 0)
    return 0;
  const double phi = Dihedral(a[i].position,a[j].position,a[k].position,a[l].position);
  double u = 0;
  LstItr<TorsionParam> ip;
  for (ip.init(p); ip.ok(); ip.next())
    u += ip().v*(1+ip().sign*cos(ip().n*phi));
  return u;
}

void Torsion::set_lambda(double lambda)
{
  LstItrMod<TorsionParam> ip;
  LstItr<TorsionParam> iq;
  for (ip.init(p), iq.init(p1); iq.ok(); ip.next(), iq.next())
    ip().v = (1-lambda)*iq().v;
  for (iq.init(p2); iq.ok(); ip.next(), iq.next())
    ip().v = lambda*iq().v;
  insist(!ip.ok());
}

void DihedralRestraint::add_to_energy_and_forces(const Atom *a, double &u, Cartesian *f) const
{
  Cartesian g1, g2, g3, g4;
  const double phi = DihedralGradient(a[atoms.a].position,
                                      a[atoms.b].position,
                                      a[atoms.c].position,
                                      a[atoms.d].position,
                                      g1, g2, g3, g4);
  const double dphi = periodic(phi - degrees_to_radians(phi0), 2*M_PI);
  const double tmpk = k*sq(180/M_PI);
  u += 0.5*tmpk*sq(dphi);
  const double du = tmpk*dphi;
  f[atoms.a] -= du*g1;
  f[atoms.b] -= du*g2;
  f[atoms.c] -= du*g3;
  f[atoms.d] -= du*g4;
}

void DihedralRestraint::init(int n, const Atom *a)
{
  if (atoms.a < 0 || atoms.a >= n)
    die("DihedralRestraint::check: index %d is out of range", atoms.a);
  if (atoms.b < 0 || atoms.b >= n)
    die("DihedralRestraint::check: index %d is out of range", atoms.b);
  if (atoms.c < 0 || atoms.c >= n)
    die("DihedralRestraint::check: index %d is out of range", atoms.c);
  if (atoms.d < 0 || atoms.d >= n)
    die("DihedralRestraint::check: index %d is out of range", atoms.d);
  if (k < 0)
    die("DihedralRestraint::check: must set force constant 'k' to a number > 0");
  if (phi0_from_initial_geometry)
    phi0 = radians_to_degrees(Dihedral(a[atoms.a].position, a[atoms.b].position, a[atoms.c].position, a[atoms.d].position));
}

void IntramolecularPotential::add_to_energy_and_forces(double &u, Cartesian *f)
{
  if (!any_stretch && !any_bend && !any_torsion && !any_improper && !any_dihedral_restraint)
    return;
  TIMESTART("Intra::add_to_energy_and_forces");
  ustretch = 0;
  LstItr<Stretch> s;
  if (include_stretch)
    for (s.init(stretch_list); s.ok(); s.next())
      s().add_to_energy_and_forces(msys->atom, ustretch, f);
  ubend = 0;
  LstItr<Bend> b;
  if (include_bend)
    for (b.init(bend_list); b.ok(); b.next())
      b().add_to_energy_and_forces(msys->atom, ubend, f);
  utorsion = 0;
  LstItr<Torsion> t;
  if (include_torsion)
    for (t.init(torsion_list); t.ok(); t.next())
      t().add_to_energy_and_forces(msys->atom, utorsion, f);
  uimproper = 0;
  if (include_improper)
    for (t.init(improper_list); t.ok(); t.next())
      t().add_to_energy_and_forces(msys->atom, uimproper, f);
  udihedral_restraint = 0;
  LstItr<DihedralRestraint> dh;
  for (dh.init(dihedral_restraint); dh.ok(); dh.next())
    dh().add_to_energy_and_forces(msys->atom, udihedral_restraint, f);
  u += ustretch + ubend + utorsion + uimproper + udihedral_restraint;
  TIMESTOP("Intra::add_to_energy_and_forces");
}

void IntramolecularPotential::write() const
{
  if (verbose) {
    if (any_stretch)
      Out() << "Stretch energy: " << ustretch << " kcal/mol\n";
    if (any_bend)
      Out() << "Bend energy: " << ubend << " kcal/mol\n";
    if (any_torsion)
      Out() << "Torsion energy: " << utorsion << " kcal/mol\n";
    if (any_improper)
      Out() << "Improper torsion energy: " << uimproper << " kcal/mol\n";
    if (any_dihedral_restraint)
      Out() << "Dihedral restraint energy: " << udihedral_restraint << " kcal/mol\n";
  }
  if (verbose >= 2)
    write_details();
  if (verbose >= 3)
    Out() << *this << flush;
}

void IntramolecularPotential::write_details() const
{
  const Atom *a = msys->atom;
  Out() << "Atoms:\n";
  for (int i = 0; i < msys->atom.size(); i++) {
    const Atom &ai = a[i];
    Out() << ai.symbol << ":" << ai.type << ":" << code_for_type(ai.type) << "\n";
  }
  Out() << "Stretches (type1,type2,i,j,forceconstant,r0,r,energy)\n";
  LstItr<Stretch> s;
  for (s.init(stretch_list); s.ok(); s.next()) {
    const Atom &ai = a[s().i], &aj = a[s().j];
    Out() << ai.symbol << ":" << ai.type << ":" << code_for_type(ai.type) << "  "
	  << aj.symbol << ":" << aj.type << ":" << code_for_type(aj.type) << "  "
	  << s() << " " << ai.position.distance(aj.position) << " " << s().energy(a) << "\n";
  }
  Out() << "Bends (type1,type2,type3,i,j,k,forceconstant,theta0,theta,energy)\n";
  LstItr<Bend> b;
  for (b.init(bend_list); b.ok(); b.next()) {
    const Atom &ai = a[b().i], &aj = a[b().j], &ak = a[b().k];
    Out() << ai.symbol << ":" << ai.type << ":" << code_for_type(ai.type) << "  "
	  << aj.symbol << ":" << aj.type << ":" << code_for_type(aj.type) << "  "
	  << ak.symbol << ":" << ak.type << ":" << code_for_type(ak.type) << "  "
	  << b() << " " 
	  << radians_to_degrees(Angle(ai.position,aj.position,ak.position))
	  << " "
	  << b().energy(a) << "\n";
  }
  Out() << "Torsions:\n";
  LstItr<Torsion> t;
  for (t.init(torsion_list); t.ok(); t.next()) {
    const Atom &ai = a[t().i], &aj = a[t().j], &ak = a[t().k], &al = a[t().l];
    Out() << ai.symbol << ":" << ai.type << ":" << code_for_type(ai.type) << "  "
	  << aj.symbol << ":" << aj.type << ":" << code_for_type(aj.type) << "  "
	  << ak.symbol << ":" << ak.type << ":" << code_for_type(ak.type) << "  "
	    << al.symbol << ":" << al.type << ":" << code_for_type(al.type) << "  "
	  << t() << " " << t().energy(a) << "\n";
  }
  Out() << "Impropers:\n";
  for (t.init(improper_list); t.ok(); t.next()) {
    const Atom &ai = a[t().i], &aj = a[t().j], &ak = a[t().k], &al = a[t().l];
    Out() << ai.symbol << ":" << ai.type << ":" << code_for_type(ai.type) << "  "
	  << aj.symbol << ":" << aj.type << ":" << code_for_type(aj.type) << "  "
	  << ak.symbol << ":" << ak.type << ":" << code_for_type(ak.type) << "  "
	  << al.symbol << ":" << al.type << ":" << code_for_type(al.type) << "  "
	  << t() << " " << t().energy(a) << "\n";
  }
  if (dihedral_restraint.size() > 0) {
    Out() << "Dihedral restraints (i,j,k,l,phi,phi0)\n";
    LstItr<DihedralRestraint> dh;
    for (dh.init(dihedral_restraint); dh.ok(); dh.next()) {
      const int i = dh().atoms.a;
      const int j = dh().atoms.b;
      const int k = dh().atoms.c;
      const int l = dh().atoms.d;
      OutPrintf("%4d %4d %4d %4d  %12.8f  %12.8f\n", i, j, k, l,
		radians_to_degrees(Dihedral(a[i].position, a[j].position, a[k].position, a[l].position)),
		dh().phi0);
    }
    Out() << "\n" << flush;
  }
}
