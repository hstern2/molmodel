#ifndef ELECLJ_H
#define ELECLJ_H

#include "pot.hpp"
#include "str.hpp"
#include "htab.hpp"
#include "cvec.hpp"
#include "dvec.hpp"
#include "boo.hpp"
#include "pppm.h"
#include "ntuple.hpp"
#include "units.h"
#include "cmpole.hpp"

extern "C" {
  struct Chebyshev;
  struct PairBinary;
  struct kspace;
  typedef struct p3m *p3m_t;
}

struct ElecLJConvergenceException { };

struct MSiteSpec
{ 
  Str type; // io
  double r; // io
  MSiteSpec() : type("X"), r(0.128012065) /* default: TIP4P */ { }
  classIO(MSiteSpec);
};

struct Sp3SiteSpec
{
  Str type; // io
  double r, theta; // io
  Sp3SiteSpec() : type("X"), r(0.7), theta(125.26439) { }
  classIO(Sp3SiteSpec);
};

struct Sp2SiteSpec
{
  Str type; // io
  double r, theta; // io
  Sp2SiteSpec() : type("X"), r(0.7), theta(120) { }
  classIO(Sp2SiteSpec);
};

struct Site
{
  int index; // out
  Str type1, type2; // out
  double q0, q01, q02, rscreen1, rscreen2, sig1, sig2, eps1, eps2; // out
  Site() : q0(0), q01(0), q02(0), rscreen1(0), rscreen2(0), sig1(0), sig2(0), eps1(0), eps2(0) { }
  void rezero() { q0 = q01 = q02 = rscreen1 = rscreen2 = sig1 = sig2 = eps1 = eps2 = 0; }
  virtual ~Site() { }
  virtual Cartesian position(const Atom *) const = 0;
  virtual void apply_force(const Atom *, Cartesian *, const Cartesian &) const = 0;
  virtual void set_lambda(double lambda)
  { q0 = (1-lambda)*q01 + lambda*q02; }
  classIO(Site);
  virtual void write(ostream &) const = 0;
};

struct AtomSite : public Site // out
{
  const int iatom; // out
  AtomSite(int iat) : iatom(iat) { }
  Cartesian position(const Atom *a) const { return a[iatom].position; }
  void apply_force(const Atom *, Cartesian *f, const Cartesian &sf) const 
  { f[iatom] += sf; }
  classIO(AtomSite);
  void write(ostream &s) const { s << *this; }
};
  
struct MSite : public Site // out
{
  const int a, b, c; // out
  double r1, r2, r; // out
  MSite(int a_, int b_, int c_) : a(a_), b(b_), c(c_) { }
  void init(const MSiteSpec *m1, const MSiteSpec *m2)
  {
    r1 = r2 = 0;
    if (m1) {
      type1 = m1->type;
      r1 = m1->r;
    }
    if (m2) {
      r2 = m2->r;
      type2 = m2->type;
    }
    if (!m1)
      r1 = r2;
    if (!m2)
      r2 = r1;
  }
  Cartesian position(const Atom *at) const 
  {
    return at[b].position + 
      r * ((at[a].position - at[b].position)+(at[c].position-at[b].position));
  }
  void apply_force(const Atom *, Cartesian *f, const Cartesian &sf) const
  {
    f[a] += r*sf;
    f[b] += (1 - 2*r) * sf;
    f[c] += r*sf;
  }
  void set_lambda(double lambda)
  { 
    Site::set_lambda(lambda);
    r = (1-lambda)*r1 + lambda*r2;
  }
  classIO(MSite);
  void write(ostream &s) const { s << *this; }
};

struct Sp3Site : public Site
{
  const int a, b, c; // out
  const double sign; // out
  double rct1, rst1, rct2, rst2, rct, rst; // out
  Sp3Site(int a_, int b_, int c_, double si) 
    : a(a_), b(b_), c(c_), sign(si) { }
  void init(const Sp3SiteSpec *s1, const Sp3SiteSpec *s2)
  {
    rct1 = rst1 = rct2 = rst2 = 0;
    if (s1) {
      type1 = s1->type;
      rct1 = s1->r*cos(degrees_to_radians(sign*s1->theta));
      rst1 = s1->r*sin(degrees_to_radians(sign*s1->theta));
    }
    if (s2) {
      type2 = s2->type;
      rct2 = s2->r*cos(degrees_to_radians(sign*s2->theta));
      rst2 = s2->r*sin(degrees_to_radians(sign*s2->theta));
    }
    if (!s1) {
      rct1 = rct2;
      rst1 = rst2;
    }
    if (!s2) {
      rct2 = rct1;
      rst2 = rst1;
    }
  }
  Cartesian position(const Atom *at) const
  {
    const Cartesian ab = at[a].position - at[b].position;
    const Cartesian abhat = ab.as_unit_vector();
    const Cartesian cb = at[c].position - at[b].position;
    const Cartesian cbhat = cb.as_unit_vector();
    const Cartesian uhat = (abhat+cbhat).as_unit_vector();
    const Cartesian vhat = uhat.cross(cbhat).as_unit_vector();
    return at[b].position + rct*uhat + rst*vhat;
  }
  void apply_force(const Atom *at, Cartesian *f, const Cartesian &sf) const
  {
    const Cartesian ab = at[a].position - at[b].position;
    const Cartesian abhat = ab.as_unit_vector();
    const Cartesian cb = at[c].position - at[b].position;
    const Cartesian cbhat = cb.as_unit_vector();
    const Cartesian u = abhat+cbhat;
    const Cartesian uhat = u.as_unit_vector();
    const Cartesian v = uhat.cross(cbhat);
    const Cartesian vhat = v.as_unit_vector();
    const Tensor unit(1);
    const Tensor duhat_du = (unit - Tensor(uhat,uhat))/u.magnitude();
    const Tensor dvhat_dv = (unit - Tensor(vhat,vhat))/v.magnitude();
    const Tensor dabhat_dab = (unit - Tensor(abhat,abhat))/ab.magnitude();
    const Tensor dcbhat_dcb = (unit - Tensor(cbhat,cbhat))/cb.magnitude();
    const Tensor duab = duhat_du * dabhat_dab;
    const Tensor ducb = duhat_du * dcbhat_dcb;
    Tensor tmp;
    tmp.col(0) = duab.row(0).cross(cbhat);
    tmp.col(1) = duab.row(1).cross(cbhat);
    tmp.col(2) = duab.row(2).cross(cbhat);
    const Cartesian fa = sf*(rct*duab + rst*dvhat_dv*tmp);
    tmp.col(0) = ducb.row(0).cross(cbhat) + uhat.cross(dcbhat_dcb.row(0));
    tmp.col(1) = ducb.row(1).cross(cbhat) + uhat.cross(dcbhat_dcb.row(1));
    tmp.col(2) = ducb.row(2).cross(cbhat) + uhat.cross(dcbhat_dcb.row(2));
    const Cartesian fc = sf*(rct*ducb + rst*dvhat_dv*tmp);
    f[a] += fa;
    f[b] += sf - fa - fc;
    f[c] += fc;
  }
  void set_lambda(double lambda)
  { 
    Site::set_lambda(lambda);
    rct = (1-lambda)*rct1 + lambda*rct2;
    rst = (1-lambda)*rst1 + lambda*rst2;
  }
  classIO(Sp3Site);
  void write(ostream &s) const { s << *this; }
};

struct Sp2Site : public Site
{
  const int a, b, c, i; // out
  const double sign; // out
  double rct1, rst1, rct2, rst2, rct, rst; // out
  Sp2Site(int a_, int b_, int c_, int i_, double si) : 
    a(a_), b(b_), c(c_), i(i_), sign(si)
  { }
  void init(const Sp2SiteSpec *s1, const Sp2SiteSpec *s2)
  {
    rct1 = rst1 = rct2 = rst2 = 0;
    if (s1) {
      type1 = s1->type;
      rct1 = s1->r*cos(degrees_to_radians(sign*s1->theta));
      rst1 = s1->r*sin(degrees_to_radians(sign*s1->theta));
    }
    if (s2) {
      type2 = s2->type;
      rct2 = s2->r*cos(degrees_to_radians(sign*s2->theta));
      rst2 = s2->r*sin(degrees_to_radians(sign*s2->theta));
    }
    if (!s1) {
      rct1 = rct2;
      rst1 = rst2;
    }
    if (!s2) {
      rct2 = rct1;
      rst2 = rst1;
    }
  }
  Cartesian position(const Atom *at) const
  {
    const Cartesian u = at[i].position - at[a].position;
    const Cartesian uhat = u.as_unit_vector();
    const Cartesian cb = at[c].position - at[b].position;
    const Cartesian v = cb - (cb*uhat)*uhat;
    const Cartesian vhat = v.as_unit_vector();
    return at[a].position + rct*uhat + rst*vhat;
  }
  void apply_force(const Atom *at, Cartesian *f, const Cartesian &sf) const
  {
    const Cartesian u = at[i].position - at[a].position;
    const Cartesian uhat = u.as_unit_vector();
    const Cartesian cb = at[c].position - at[b].position;
    const Cartesian v = cb - (cb*uhat)*uhat;
    const Cartesian vhat = v.as_unit_vector();
    const Tensor unit(1);
    const Tensor uhatuhat = unit - Tensor(uhat,uhat);
    const Tensor dvhat_dv = (unit - Tensor(vhat,vhat))/v.magnitude();
    const Cartesian fi = sf*(Tensor(rct) - rst*dvhat_dv*(Tensor(uhat,cb)+Tensor(uhat*cb))) * 
      uhatuhat/u.magnitude();
    const Cartesian fc = sf*rst*dvhat_dv*uhatuhat;
    f[i] += fi;
    f[a] += sf - fi;
    f[b] -= fc;
    f[c] += fc;
  }
  void set_lambda(double lambda)
  { 
    Site::set_lambda(lambda);
    rct = (1-lambda)*rct1 + lambda*rct2;
    rst = (1-lambda)*rst1 + lambda*rst2;
  }
  classIO(Sp2Site);
  void write(ostream &s) const { s << *this; }
};

struct Sites
{
  AtomSite *atom;
  MSite *m;
  Sp3Site *sp3a, *sp3b;
  Sp2Site *sp2a, *sp2b;
  Sites() : atom(0), m(0), sp3a(0), sp3b(0), sp2a(0), sp2b(0) { }
  Sites(const Sites &);
  int how_many() const;
  void add_to(Site **, int &) const;
  void add_virtual_neighbors(const Sites *, Lst<int> *) const;
};

struct BondChargeIncrement
{
  double k1, k2, k; // out
  double q01, q02, q0; // out
  double omega2, velocity; // out
  double dq; // out
  double mass() const { return k/omega2; }
  double self_energy() const { return 0.5*k*sq(dq); }
  double kinetic_energy() const { return 0.5*mass()*sq(velocity); }
  classIO(BondChargeIncrement);
};

struct OneThree
{
  double p1, p2, param; // out
  classIO(OneThree);
};

class Atom;
class NeighborList;

struct BondChargeIncrementParam
{
  double k, q0; // io
  BondChargeIncrementParam() : k(0), q0(0) { }
  classIO(BondChargeIncrementParam);
};

class ElecLJ : public Potential
{

  friend struct BCIFit;

public:

  HTab<MSiteSpec,Str3> msite; // io
  HTab<Sp3SiteSpec,Str3> sp3; // io
  HTab<Sp2SiteSpec,Str2> sp2; // io
  HTab<double> fixed_charge, screening_radius; // io
  HTab<BondChargeIncrementParam> bond_charge_increment; // io  
  HTab<double,Str3> one_three_interaction; // io
  HTab<double> sigma, epsilon; // io
  HTab<Str> general_type; // io
  double elec_1_4_scale_factor, lj_1_4_scale_factor; // io
  double rspace_cutoff, smoothing_width, neighbor_list_width; // io
  double kspace_cutoff, ewald_screening; // io
  double tolerance, bond_charge_increment_frequency; // io
  double p3m_grid_spacing, p3m_grid_density; // io    if set, density overrides spacing
  int verbose, max_iterations, p3m_alias_sum_limit, p3m_assignment_order; // io
  Boolean reorder_sites; // io
  Boolean geometric_combining_for_sigma, include_lj_correction; // io
  Boolean include_rspace, include_kspace, include_self; // io
  Boolean use_p3m, p3m_ik_differentiate; // io
  double kT; // out
  
  ElecLJ();
  ElecLJ(const ElecLJ &);
  ~ElecLJ();
  
  void init(MSys *m);
  void types_changed();
  void types_changed(const MSys &perturbed);
  void set_lambda(double lambda);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  void integrate_positions(double dt);
  void integrate_velocities(double dt);

  Cartesian dipole_moment() const;
  Quadrupole quadrupole_moment() const;
  Octopole octopole_moment() const;
  Hexadecapole hexadecapole_moment() const;

  double net_charge() const;
  Tensor polarizability() const;
  double static_dielectric_constant() const;
  double optical_dielectric_constant() const;
  const Cartesian *electric_field() const { return evec; }
  void scale(double s);
  void is_running_dynamics(bool b)  { dyn = b; }
  void temperature(double t)  { kT = K_to_kcal_mol(t); }
  void write_electrostatic_potential_at_gridpoints(const CVec gp);
  void electrostatic_potential_at_gridpoints(const CVec gp, DVec &p);
  Complex structure_factor(const Cartesian &k) const;
  double charge_on_atom(int i) const;
  
private:
  
  int nsite, nbci, ngroup; // out
  Boolean any_fixed_charge_params, any_bci_params, any_lj_params, any_virtual_sites; // out
  Boolean any_site_has_q, any_site_has_lj, any_bci, have_solved_for_bci; // out
  Boolean dyn; // out
  double totA, totB, lj_correction, elec_energy, lj_energy, rc, nlw, eta; // out
  DVec q, sig, eps, phi, rscreen; // out
  CVec evec, r, gr, lj_force; // out
  Vec<Site *> site;
  Vec<Sites> site_on_atom;
  HTab<BondChargeIncrement,Int2> bci; // out
  Vec<int> site_mask, lj_type; // out
  Vec<Lst<int> > neighbors; // out
  NeighborList *nlist; // non-bonded group neighbors
  HTab<OneThree,Int2> one_three; // out
  Vec<Int2> group; // out
  Vec<int> group_for_site; // out
  struct PairBinary *excluded, *group_excluded, *one_four;
  struct Chebyshev *rspace_chebyshev;
  struct kspace *kspace;
  p3m_t p3m;

  void init_sites();
  void types_changed(bool all, const Vec<int> index);
  void coordinates_changed();
  void calc_phi_one_three(), calc_phi(), calc_phi_nonperiodic(), calc_phi_cubic();
  void calc_phi_evec_lj(), calc_self_phi();
  void solve_for_bci(const Cartesian &external_field);
  void zero_bci();
  void update_q(), recalculate_lj_correction(), cleanup();

  const char *bond_type_for_sites(const Atom *, int, int) const;
  
  classIO(ElecLJ);

};

#endif /* ELECLJ_H */
