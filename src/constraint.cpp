#include "constraint.hpp"
#include "cvec.hpp"
#include "linalg.hpp"
#include "units.h"
#include "timing.h"
#include "geograd.hpp"
#include "pot.hpp"

class Cnstr
{
public:
  const int natom;
  int i[4];
  Cnstr(int n) : natom(n) { }
  virtual void write() const = 0;
  virtual void update(const Atom *) = 0;
  virtual void update(const Cartesian *) = 0;
  virtual double value() const = 0;
  virtual Cartesian deriv(int k) const = 0;
  virtual unsigned hash_val() const = 0;
  virtual bool matches(const Cnstr *) const = 0;
  virtual bool params_match(const Cnstr *) const = 0;
  virtual ~Cnstr() { }
};

unsigned hash_val(const CnstrPtr &key) { return key.p->hash_val(); }
bool operator==(const CnstrPtr &a, const CnstrPtr &b)
{ return a.p->matches(b.p); }

class LengthCnstr : public Cnstr
{
  friend class Constraint;
private:
  const double len;
  Cartesian r;
public:
  LengthCnstr(int i0, int i1, double l) : Cnstr(2), len(l)
  { 
    i[0] = i0;
    i[1] = i1;
  }
  unsigned hash_val() const { return 31*i[0] + i[1]; }
  bool matches(const Cnstr *c) const 
  {
    const LengthCnstr *cc = dynamic_cast<const LengthCnstr *>(c);
    return cc && i[0] == cc->i[0] && cc->i[1];
  }
  bool params_match(const Cnstr *c) const 
  {
    const LengthCnstr *cc = dynamic_cast<const LengthCnstr *>(c);
    return cc && are_approximately_equal(len,cc->len);
  }
  void write() const { OutPrintf("length constraint %d-%d: %f\n", i[0],i[1],len); }
  void update(const Atom *a)
  { r = a[i[0]].position - a[i[1]].position; }
  void update(const Cartesian *a)
  { r = a[i[0]] - a[i[1]]; }
  double value() const
  { return 0.5*(r.sq() - len*len); }
  Cartesian deriv(int k) const
  {
    if (k == i[0])
      return r;
    if (k == i[1])
      return -r;
    return Cartesian(0,0,0);
  }
};

class AngleCnstr : public Cnstr
{
  friend class Constraint;
private:
  const double t0;
  double t;
  Cartesian g[3];
public:
  AngleCnstr(int i0, int i1, int i2, const double t0) 
    : Cnstr(3), t0(degrees_to_radians(t0))
  { 
    i[0] = i0;
    i[1] = i1;
    i[2] = i2;
  }
  unsigned hash_val() const { return 31*(31*i[0] + i[1])+i[2]; }
  bool matches(const Cnstr *c) const 
  {
    const AngleCnstr *cc = dynamic_cast<const AngleCnstr *>(c);
    return cc && i[0] == cc->i[0] && i[1] == cc->i[1] && i[2] == cc->i[2];
  }
  bool params_match(const Cnstr *c) const 
  {
    const AngleCnstr *cc = dynamic_cast<const AngleCnstr *>(c);
    return cc && are_approximately_equal(t0,cc->t0);
  }
  void write() const { OutPrintf("angle constraint %d-%d-%d: %f\n", i[0],i[1],i[2],
				 radians_to_degrees(t0)); }
  void update(const Atom *a)
  { t = AngleGradient(a[i[0]].position,a[i[1]].position,a[i[2]].position,g[0],g[1],g[2]); }
  void update(const Cartesian *a)
  { t = AngleGradient(a[i[0]],a[i[1]],a[i[2]],g[0],g[1],g[2]); }
  double value() const
  { return t - t0; }
  Cartesian deriv(int k) const
  {
    if (k == i[0])
      return g[0];
    if (k == i[1])
      return g[1];
    if (k == i[2])
      return g[2];
    return Cartesian(0,0,0);
  }
};

class DihedralCnstr : public Cnstr
{
  friend class Constraint;
private:
  const double phi0;
  double phi;
  Cartesian g[4];
public:
  DihedralCnstr(int i0, int i1, int i2, int i3, const double phi0_) 
    : Cnstr(4), phi0(degrees_to_radians(phi0_))
  { 
    i[0] = i0;
    i[1] = i1;
    i[2] = i2;
    i[3] = i3;
  }
  unsigned hash_val() const { return 31*(31*(31*i[0] + i[1])+i[2])+i[3]; }
  bool matches(const Cnstr *c) const 
  {
    const DihedralCnstr *cc = dynamic_cast<const DihedralCnstr *>(c);
    return cc && i[0] == cc->i[0] && i[1] == cc->i[1] && i[2] == cc->i[2] && i[3] == cc->i[3];
  }
  bool params_match(const Cnstr *c) const 
  {
    const DihedralCnstr *cc = dynamic_cast<const DihedralCnstr *>(c);
    return cc && are_approximately_equal(phi0,cc->phi0);
  }
  void write() const { OutPrintf("dihedral constraint %d-%d-%d-%d: %f\n", i[0],i[1],i[2],i[3],
				 radians_to_degrees(phi0)); }
  void update(const Atom *a)
  { 
    phi = DihedralGradient(a[i[0]].position,a[i[1]].position,a[i[2]].position,a[i[3]].position,
			   g[0],g[1],g[2],g[3]); 
  }
  void update(const Cartesian *a)
  { phi = DihedralGradient(a[i[0]],a[i[1]],a[i[2]],a[i[3]],g[0],g[1],g[2],g[3]); }
  double value() const
  { return phi - phi0; }
  Cartesian deriv(int k) const
  {
    if (k == i[0])
      return g[0];
    if (k == i[1])
      return g[1];
    if (k == i[2])
      return g[2];
    if (k == i[3])
      return g[3];
    return Cartesian(0,0,0);
  }
};

struct Sp2ProtonCnstrHelp
{
  Cartesian r12_last, r32_last;
  Cartesian uhat;
  Tensor g[4];
  Sp2ProtonCnstrHelp()
  {
    g[0] = Tensor(1);
  }
  void update(const Cartesian r12, const Cartesian r32) 
  {
    if (r12 == r12_last && r32 == r32_last)
      return;
    r12_last = r12;
    r32_last = r32;
    const Cartesian r12hat = r12.as_unit_vector();
    const Cartesian r32hat = r32.as_unit_vector();
    const Cartesian u = r12hat + r32hat;
    uhat = u.as_unit_vector();
    const Tensor duhat_du = (1.0/u.magnitude()) * (g[0] - Tensor(uhat,uhat));
    const Tensor dr12hat_dr12 = (1.0/r12.magnitude()) * (g[0] - Tensor(r12hat,r12hat));
    const Tensor dr32hat_dr32 = (1.0/r32.magnitude()) * (g[0] - Tensor(r32hat,r32hat));
    g[1] = duhat_du*dr12hat_dr12;
    g[3] = duhat_du*dr32hat_dr32;
    g[2] = -g[1]-g[3];
  }
};

class Sp2ProtonCnstr : public Cnstr
{
private:
  double len;
  const int icoord;
  Sp2ProtonCnstrHelp *help;
  double r0, r2;
public:
  Sp2ProtonCnstr(int i0, int i1, int i2, int i3, int icoord_, Sp2ProtonCnstrHelp *h, double l) 
    : Cnstr(4), len(l), icoord(icoord_), help(h)
  { 
    i[0] = i0;
    i[1] = i1;
    i[2] = i2;
    i[3] = i3;
  }
  ~Sp2ProtonCnstr() { if (icoord == 0) delete help; }
  unsigned hash_val() const { return 31*(31*(31*i[0] + i[1])+i[2])+i[3]; }
  bool matches(const Cnstr *c) const 
  {
    const Sp2ProtonCnstr *cc = dynamic_cast<const Sp2ProtonCnstr *>(c);
    return cc && i[0] == cc->i[0] && i[1] == cc->i[1] && i[2] == cc->i[2] && 
      i[3] == cc->i[3] && icoord == cc->icoord;
  }
  bool params_match(const Cnstr *c) const 
  {
    const Sp2ProtonCnstr *cc = dynamic_cast<const Sp2ProtonCnstr *>(c);
    return cc && are_approximately_equal(len,cc->len);
  }
  void write() const { 
    const char xyz[4] = "xyz";
    OutPrintf("sp2 constraint %d-%d-%d-%d (%c): %f (current value: %f)\n", i[0],i[1],i[2],i[3],
	      xyz[icoord],len,value());
  }
  void update(const Atom *a) { 
    help->update(a[i[1]].position - a[i[2]].position,
		 a[i[3]].position - a[i[2]].position);
    r0 = a[i[0]].position.coord(icoord);
    r2 = a[i[2]].position.coord(icoord);
  }
  void update(const Cartesian *a) { 
    help->update(a[i[1]]-a[i[2]],a[i[3]]-a[i[2]]); 
    r0 = a[i[0]].coord(icoord);
    r2 = a[i[2]].coord(icoord);
  }
  void place(Atom *a) const {
    help->update(a[i[1]].position - a[i[2]].position,
		 a[i[3]].position - a[i[2]].position);
    a[i[0]].position.coord(icoord) = a[i[2]].position.coord(icoord) - len*help->uhat.coord(icoord);
  }
  double value() const { return r0 + len*help->uhat.coord(icoord) - r2; }
  Cartesian deriv(int k) const
  {
    if (k == i[0])
      return help->g[0].row(icoord);
    if (k == i[1])
      return len*help->g[1].row(icoord);
    if (k == i[2])
      return len*help->g[2].row(icoord) - help->g[0].row(icoord);
    if (k == i[3])
      return len*help->g[3].row(icoord);
    return Cartesian(0,0,0);
  }
};

class ConstraintGroup
{
  friend class Constraint;
  int nc, ni;
  Vec<Cnstr *> cnstr;
  Vec<int> index;
  Vec<CVec> sci0, sci;
  DVec x;
  DMat s;
public:
  ConstraintGroup(const Lst<Cnstr *> &clst) : 
    nc(clst.size()), sci0(nc), sci(nc), x(nc), s(nc,nc)
  {
    LstItr<Cnstr *> c;
    HSet<int> iat;
    int k;
    for (c.init(clst); c.ok(); c.next())
      for (k = 0; k < c()->natom; k++)
	iat.add(c()->i[k]);
    cnstr.copy_from_list(clst);
    index.copy_from_set(iat);
    ni = index.size();
    for (k = 0; k < nc; k++) {
      sci0[k] = CVec(ni);
      sci[k] = CVec(ni);
    }
  }
  void write() const
  {
    insist(nc == cnstr.size());
    for (int i = 0; i < nc; i++)
      cnstr[i]->write();
    Out() << "\n";
  }
  void constrain_positions(const CVec r0, Atom *a,
			   double dt, double tol,
			   int maxit, int verbose)
  {
    int i, c;
    for (c = 0; c < nc; c++) {
      cnstr[c]->update(r0);
      for (i = 0; i < ni; i++)
	sci0[c][i] = cnstr[c]->deriv(index[i])/a[index[i]].mass();
    }
    for (int iter = 0; iter < maxit; iter++) {
      bool is_satisfied = true;
      for (c = 0; c < nc; c++) {
	cnstr[c]->update(a);
	x[c] = -cnstr[c]->value();
	if (fabs(x[c]) > tol)
	  is_satisfied = false;
      }
      if (is_satisfied) {
	if (verbose)
	  Out() << "ConstraintGroup::constrain_positions: satisfied tolerance of "
		<< tol << " A " << " in " << iter+1 << " iterations\n" << flush;
	return;
      }
      for (c = 0; c < nc; c++)
	for (i = 0; i < ni; i++)
	  sci[c][i] = cnstr[c]->deriv(index[i]);
      for (c = 0; c < nc; c++)
	for (int cb = 0; cb < nc; cb++)
	  s(c,cb) = sci[c]*sci0[cb];
      LinearSolve(s,x);
      const double invdt = dt > 0 ? 1/dt : 0;
      for (c = 0; c < nc; c++)
	for (i = 0; i < ni; i++) {
	  Atom &ai = a[index[i]];
	  const Cartesian dx = x[c] * sci0[c][i];
	  ai.position += dx;
	  ai.velocity += dx*invdt;
	}
    }
  }
  void constrain_velocities(Vec<Atom> a)
  {
    int i, c;
    for (c = 0; c < nc; c++) {
      x[c] = 0;
      cnstr[c]->update(a);
      for (i = 0; i < ni; i++) {
	const Atom &ai = a[index[i]];
	sci[c][i] = cnstr[c]->deriv(index[i]);
	sci0[c][i] = sci[c][i]/ai.mass();
	x[c] -= sci[c][i]*ai.velocity;
      }
    }
    for (c = 0; c < nc; c++)
      for (int cb = c; cb < nc; cb++)
	s(c,cb) = s(cb,c) = sci[c]*sci0[cb];
    LinearSolve(s, x);
    for (c = 0; c < nc; c++)
      for (i = 0; i < ni; i++) {
	Atom &ai = a[index[i]];
	ai.velocity += x[c] * sci0[c][i];
      }
  }
  void add_to_energy_and_forces(double &u, Cartesian *f, const Atom *a)
  {
    for (int c = 0; c < nc; c++) {
      cnstr[c]->update(a);
      const double val = 50*cnstr[c]->value();
      u += sq(val);
      for (int i = 0; i < ni; i++)
	f[index[i]] -= 100*val*cnstr[c]->deriv(index[i]);
    }
  }
};

Constraint::Constraint() : tolerance(1e-8), max_iterations(10000), verbose(0), msys(0), ncnstr(0) { }
Constraint::~Constraint() { cleanup(); }

void Constraint::add_constraint(Cnstr *c)
{
  CnstrPtr p(c);
  CnstrPtr *pold = cnstr.get(p);
  if (!pold) {
    cnstr.add(p);
    return;
  }
  Cnstr *cold = pold->p;
  if (!c->params_match(cold)) {
    Out() << "Constraint::add_constraint: incompatible constraints\n";
    c->write();
    cold->write();
    die("");
  }
  delete c;
} 

void Constraint::add_constraints(const MSys &msys_)
{
  const int natom = msys_.atom.size();
  insist(msys->atom.size() == natom);
  const Atom *a = msys_.atom;
  for (int i = 0; i < natom; i++) {
    LstItr<int> j;
    for (j.init(a[i].neighbors); j.ok(); j.next()) {
      const double *len = 0;
      if (i < j() &&
	  ((len = bond_length.get(Str2(a[i].type,a[j()].type))) != 0 ||
	   (len = bond_length.get(Str2(a[j()].type,a[i].type))) != 0))
	add_constraint(new LengthCnstr(i,j(),*len));
    }
    if (!strcmp(a[i].symbol,"H") && 
	a[i].neighbors.size() == 1 &&
	a[a[i].neighbors.first()].neighbors.size() == 3) {
      const int n2 = a[i].neighbors.first();
      int n1 = -1, n3 = -1;
      LstPairItr<int> j;
      for (j.init(a[n2].neighbors); j.ok(); j.next())
	if (j.i() != i && j.j() != i) {
	  n1 = j.i();
	  n3 = j.j();
	  break;
	}
      insist(n1 >= 0 && n3 >= 0);
      const double *len = 0;
      if (((len = sp2_proton.get(Str4(a[i].type,a[n1].type,a[n2].type,a[n3].type))) != 0) ||
	  ((len = sp2_proton.get(Str4(a[i].type,a[n3].type,a[n2].type,a[n1].type))) != 0)) {
	Sp2ProtonCnstrHelp *h = new Sp2ProtonCnstrHelp;
	for (int icoord = 0; icoord < 3; icoord++)
	  add_constraint(new Sp2ProtonCnstr(i,n1,n2,n3,icoord,h,*len));
      }
    }
    LstPairItr<int> n;
    double *t0 = 0;
    int nc = 0;
    const int nmax = 2*a[i].neighbors.size() - 3;
    for (n.init(a[i].neighbors); n.ok(); n.next()) {
      Str3 sij(a[n.i()].type,a[i].type,a[n.j()].type);
      Str3 sji(a[n.j()].type,a[i].type,a[n.i()].type);
      if ((t0 = angle.get(sij)) != 0 || (t0 = angle.get(sji)) != 0) {
	add_constraint(new AngleCnstr(n.i(),i,n.j(),*t0));
	if (++nc == nmax)
	  break;
      }
    }
  }
  HTabIterator<double,Int4> ti;
  for (ti.init(dihedral); ti.ok(); ti.next()) {
    if (ti.key().a < 0 || ti.key().a >= natom)
      die("Constraint::add_constraint: index %d is out of range", ti.key().a);
    if (ti.key().b < 0 || ti.key().b >= natom)
      die("Constraint::add_constraint: index %d is out of range", ti.key().a);
    if (ti.key().c < 0 || ti.key().c >= natom)
      die("Constraint::add_constraint: index %d is out of range", ti.key().a);
    if (ti.key().d < 0 || ti.key().d >= natom)
      die("Constraint::add_constraint: index %d is out of range", ti.key().a);
    add_constraint(new DihedralCnstr(ti.key().a,ti.key().b,ti.key().c,ti.key().d,ti.val()));
  }
  HTabIterator<double,Int2> di;
  for (di.init(distance); di.ok(); di.next()) {
    if (di.key().a < 0 || di.key().a >= natom) 
      die("Constraint::add_constraint: index %d is out of range for distance constraint", ti.key().a);
    if (di.key().b < 0 || di.key().b >= natom) 
      die("Constraint::add_constraint: index %d is out of range for distance constraint", ti.key().a);
    add_constraint(new LengthCnstr(di.key().a, di.key().b, di.val()));
  }
  cnstr_on_atom = Vec<Lst<Cnstr *> >(natom);
  HSetItr<CnstrPtr> c;
  for (c.init(cnstr); c.ok(); c.next()) {
    for (int i = 0; i < c().p->natom; i++)
      cnstr_on_atom[c().p->i[i]].add(c().p);
  }
  LstItrMod<ConstraintGroup *> g;
  for (g.init(group); g.ok(); g.next()) {
    if (g())
      delete g();
    g() = 0;
  }
  group.remove_all();
  HSet<void *> is_in_group;
  ncnstr = 0;
  for (c.init(cnstr); c.ok(); c.next()) {
    struct tmp { 
      void f(Cnstr *cn, HSet<void *> &in_group, Lst<Cnstr *> &g, Lst<Cnstr *> *on_atom)
      {
	if (!in_group.exists(cn)) {
	  in_group.add(cn);
	  g.add(cn);
	  for (int i = 0; i < cn->natom; i++) {
	    LstItr<Cnstr *> d;
	    for (d.init(on_atom[cn->i[i]]); d.ok(); d.next())
	      f(d(),in_group,g,on_atom);
	  }
	}
      }
    };
    Lst<Cnstr *> grp;
    tmp t;
    t.f(c().p,is_in_group,grp,cnstr_on_atom);
    if (grp.size() > 0)
      group.add(new ConstraintGroup(grp));
    ncnstr++;
  }
  group_on_atom = Vec<ConstraintGroup *>(natom);
  group_on_atom.set(0);
  for (g.init(group); g.ok(); g.next())
    for (int i = 0; i < g()->index.size(); i++) {
      const int k = g()->index[i];
      insist(group_on_atom[k] == 0 || group_on_atom[k] == g());
      group_on_atom[k] = g();
    }
}

void Constraint::init(MSys *msys_)
{
  cleanup();
  msys = msys_;
  add_constraints(*msys);
  /* Add constraints for each protonation state */
  for (int i = 0; i < msys->acid.size(); i++) {
    AcidicGroup &a = msys->acid[i];
    const int stmp = a.state;
    for (int s = 0; s < a.desc->number_of_states(); s++)
      if (s != stmp) {
	msys->change_protonation_state(i,s);
	add_constraints(*msys);
      }
    msys->change_protonation_state(i,stmp);
  }
  r0 = CVec(msys->atom.size());
  reset();
  if (verbose)
    write();
}

void Constraint::reset()
{
  msys->copy_positions_to(r0);
}

void Constraint::write() const
{
  LstItr<ConstraintGroup *> g;
  for (g.init(group); g.ok(); g.next())
    g()->write();
}

void Constraint::constrain_positions(double dt)
{
  TIMESTART("constrain_positions");
  LstItrMod<ConstraintGroup *> g;
  if (verbose > 1) {
    Out() << "Before constraining positions:\n";
    msys->show_as_xyz();
  }
  for (g.init(group); g.ok(); g.next())
    g()->constrain_positions(r0,msys->atom,dt,tolerance,max_iterations,verbose);
  msys->copy_positions_to(r0);
  if (verbose)
    write();
  if (verbose > 1) {
    Out() << "After constraining positions:\n";
    msys->show_as_xyz();
  }
  TIMESTOP("constrain_positions");
}

void Constraint::constrain_velocities()
{
  TIMESTART("constrain_velocities");
  LstItrMod<ConstraintGroup *> g;
  for (g.init(group); g.ok(); g.next())
    g()->constrain_velocities(msys->atom);
  TIMESTOP("constrain_velocities");
}

void Constraint::add_to_energy_and_forces(double &u, Cartesian *f)
{
  LstItrMod<ConstraintGroup *> g;
  for (g.init(group); g.ok(); g.next())
    g()->add_to_energy_and_forces(u,f,msys->atom);
}

void Constraint::cleanup()
{
  HSetItr<CnstrPtr> c;
  for (c.init(cnstr); c.ok(); c.next())
    delete c().p;
  cnstr.remove_all();
  LstItrMod<ConstraintGroup *> g;
  for (g.init(group); g.ok(); g.next())
    delete g();
  group.remove_all();
  msys = 0;
  ncnstr = 0;
}

int Constraint::how_many() const
{
  return ncnstr;
}

int Constraint::how_many_on_atom(int i) const
{
  LstItr<Cnstr *> c;
  int n = 0;
  for (c.init(cnstr_on_atom[i]); c.ok(); c.next())
    n++;
  return n;
}
