#include "geohist.hpp"
#include "out.hpp"
#include "msys.hpp"
#include "units.h"
#include "geograd.hpp"

void BondStatistics::init(const MSys &msys)
{
  if (strlen(types.a) == 0)
    die("BondStatistics::init: need to specify types.a");
  if (strlen(types.b) == 0)
    die("BondStatistics::init: need to specify types.b");
  lst.remove_all();
  histogram.init();
  LstItr<int> n;
  for (int i = 0; i < msys.atom.size(); i++) {
    const Str si = msys.atom[i].type;
    for (n.init(msys.atom[i].neighbors); n.ok(); n.next()) {
      const Str sj = msys.atom[n()].type;
      if ((!strcmp(si,types.a) && !strcmp(sj,types.b)) ||
	  (!strcmp(si,types.b) && !strcmp(sj,types.a)))
	lst.add(Int2(i,n()));
    }
  }
  if (lst.size() == 0)
    Out() << "*** BondStatistics::init: there are no "
	  << types.a << "-" << types.b << " bonds in this molecular system.\n"
	  << flush;
}

void BondStatistics::update(const MSys &msys)
{
  LstItr<Int2> w;
  for (w.init(lst); w.ok(); w.next()) {
    const Atom &a = msys.atom[w().a];
    const Atom &b = msys.atom[w().b];
    if (!a.is_dummy() && !b.is_dummy())
      histogram.update(a.position.distance(b.position));
  }
}

void BondStatistics::write() const
{
  Out() << "Average " << types.a << "-" << types.b << " bond length: "
	<< histogram.mean() << " A\n";
  histogram.write(file);
}

void AngleStatistics::init(const MSys &msys)
{
  if (strlen(types.a) == 0)
    die("AngleStatistics::init: need to specify types.a");
  if (strlen(types.b) == 0)
    die("AngleStatistics::init: need to specify types.b");
  if (strlen(types.c) == 0)
    die("AngleStatistics::init: need to specify types.c");
  lst.remove_all();
  histogram.init();
  for (int i = 0; i < msys.atom.size(); i++) {
    const Str sj = msys.atom[i].type;
    if (strcmp(sj,types.b))
      continue;
    LstPairItr<int> n;
    for (n.init(msys.atom[i].neighbors); n.ok(); n.next()) {
      const Str si = msys.atom[n.i()].type;
      const Str sk = msys.atom[n.j()].type;
      if ((!strcmp(si,types.a) && !strcmp(sk,types.c)) ||
	  (!strcmp(si,types.c) && !strcmp(sk,types.a)))
	lst.add(Int3(n.i(),i,n.j()));
    }
  }
  if (lst.size() == 0)
    Out() << "*** AngleStatistics::init: there are no "
	  << types.a << "-" << types.b << "-" << types.c << " angles in this molecular system.\n"
	  << flush;
}

void AngleStatistics::update(const MSys &msys)
{
  LstItr<Int3> w;
  for (w.init(lst); w.ok(); w.next()) {
    const Atom &a = msys.atom[w().a];
    const Atom &b = msys.atom[w().b];
    const Atom &c = msys.atom[w().c];
    if (!a.is_dummy() && !b.is_dummy() && !c.is_dummy())
      histogram.update(radians_to_degrees(Angle(a.position,b.position,c.position)));
  }
}

void AngleStatistics::write() const
{
  Out() << "Average " << types.a << "-" << types.b << "-" << types.c << " bond angle: "
	<< histogram.mean() << " degrees\n";
  histogram.write(file);
}

void DihedralStatistics::init(const MSys &msys)
{
  if (strlen(types.a) == 0)
    die("DihedralStatistics::init: need to specify types.a");
  if (strlen(types.b) == 0)
    die("DihedralStatistics::init: need to specify types.b");
  if (strlen(types.c) == 0)
    die("DihedralStatistics::init: need to specify types.c");
  if (strlen(types.d) == 0)
    die("DihedralStatistics::init: need to specify types.d");
  lst.remove_all();
  histogram.init();
  Lst<Int2> bonds;
  LstItr<int> ni, nl;
  const int natom = msys.atom.size();
  const Atom *atom = msys.atom;
  for (int i = 0; i < natom; i++) {
    for (ni.init(atom[i].neighbors); ni.ok(); ni.next())
      if (i < ni() && ((!strcmp(atom[i].type, types.b) && !strcmp(atom[ni()].type, types.c)) ||
		       (!strcmp(atom[i].type, types.c) && !strcmp(atom[ni()].type, types.b))))
	bonds.add(Int2(i,ni()));
  }
  LstItr<Int2> b;
  for (b.init(bonds); b.ok(); b.next()) {
    const int j = b().a, k = b().b;
    for (ni.init(atom[j].neighbors); ni.ok(); ni.next())
      if (ni() != k)
	for (nl.init(atom[b().b].neighbors); nl.ok(); nl.next())
	  if (nl() != j && nl() != ni())
	    if ((!strcmp(atom[ni()].type, types.a) &&
		 !strcmp(atom[j].type, types.b) &&
		 !strcmp(atom[k].type, types.c) &&
		 !strcmp(atom[nl()].type, types.d)) ||
		(!strcmp(atom[ni()].type, types.d) &&
		 !strcmp(atom[j].type, types.c) &&
		 !strcmp(atom[k].type, types.b) &&
		 !strcmp(atom[nl()].type, types.a)))
	      lst.add(Int4(nl(),k,j,ni()));
  }
  if (lst.size() == 0)
    Out() << "*** DihedralStatistics::init: there are no "
	  << types.a << "-" << types.b << "-" << types.c << "-" << types.d << " dihedrals in this molecular system.\n"
	  << flush;
}

void DihedralStatistics::update(const MSys &msys)
{
  LstItr<Int4> w;
  for (w.init(lst); w.ok(); w.next()) {
    const Atom &a = msys.atom[w().a];
    const Atom &b = msys.atom[w().b];
    const Atom &c = msys.atom[w().c];
    const Atom &d = msys.atom[w().d];
    if (!a.is_dummy() && !b.is_dummy() && !c.is_dummy() && !d.is_dummy())
      histogram.update(radians_to_degrees(Dihedral(a.position,b.position,c.position,d.position)));
  }
}

void DihedralStatistics::write() const
{
  Out() << "Average " << types.a << "-" << types.b << "-" << types.c << "-" << types.d << " dihedral angle: "
	<< histogram.mean() << " degrees\n";
  histogram.write(file);
}

void GeometryStatistics::init(const MSys &msys, int write_interval)
{
  Callback::init(msys,write_interval);
  LstItrMod<BondStatistics> b;
  for (b.init(bond); b.ok(); b.next())
    b().init(msys);
  LstItrMod<AngleStatistics> a;
  for (a.init(angle); a.ok(); a.next())
    a().init(msys);
  LstItrMod<DihedralStatistics> c;
  for (c.init(dihedral); c.ok(); c.next())
    c().init(msys);
}

void GeometryStatistics::update(const MSys &msys) 
{
  LstItrMod<BondStatistics> b;
  for (b.init(bond); b.ok(); b.next())
    b().update(msys);
  LstItrMod<AngleStatistics> a;
  for (a.init(angle); a.ok(); a.next())
    a().update(msys);
  LstItrMod<DihedralStatistics> c;
  for (c.init(dihedral); c.ok(); c.next())
    c().update(msys);
}

void GeometryStatistics::write() const
{
  LstItr<BondStatistics> b;
  for (b.init(bond); b.ok(); b.next())
    b().write();
  LstItr<AngleStatistics> a;
  for (a.init(angle); a.ok(); a.next())
    a().write();
  LstItr<DihedralStatistics> c;
  for (c.init(dihedral); c.ok(); c.next())
    c().write();
}
