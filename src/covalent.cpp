#include "covalent.hpp"

void Covalent::init(MSys *msys)
{
  Potential::init(msys);
  u1 = u2 = u = 0;
  types_changed();
}

void Covalent::add_to_energy_and_forces(double &utmp, Cartesian *)
{
  utmp += u;
}

void Covalent::write() const
{
  if (verbose)
    Out() << "Covalent energy: " << u << " kcal/mol\n";
}

void Covalent::types_changed()
{
  u1 = u2 = u = 0;
  const Atom *at = msys->atom;
  const int nat = msys->atom.size();
  for (int i = 0; i < nat; i++) {
    LstItr<int> j;
    for (j.init(at[i].neighbors); j.ok(); j.next()) {
      if (i >= j())
	continue;
      double *p = param.get(Str2(at[i].type,at[j()].type));
      if (!p)
	p = param.get(Str2(at[j()].type,at[i].type));
      if (p)
	u1 += *p;
    }
  }
  u2 = u = u1;
}

void Covalent::types_changed(const MSys &pert)
{
  u1 = u2 = u = 0;
  const Atom *at = msys->atom;
  const Atom *pat = pert.atom;
  const int nat = msys->atom.size();
  insist(pert.atom.size() == nat);
  for (int i = 0; i < nat; i++) {
    LstItr<int> j;
    for (j.init(at[i].neighbors); j.ok(); j.next()) {
      if (i >= j())
	continue;
      double *p = 0;
      p = param.get(Str2(at[i].type,at[j()].type));
      if (!p)
	p = param.get(Str2(at[j()].type,at[i].type));
      if (p)
	u1 += *p;
    }
    for (j.init(pat[i].neighbors); j.ok(); j.next()) {
      if (i >= j())
	continue;
      double *p = 0;
      p = param.get(Str2(pat[i].type,pat[j()].type));
      if (!p)
	p = param.get(Str2(pat[j()].type,pat[i].type));
      if (p)
	u2 += *p;
    }
  }
  set_lambda(0);
}

void Covalent::set_lambda(double lambda)
{
  u = (1-lambda)*u1 + lambda*u2;
}
