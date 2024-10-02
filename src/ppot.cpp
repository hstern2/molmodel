#include "ppot.hpp"
#include "fns.h"
#include "nbrlist.hpp"
#include "pbin.h"

bool PairPotentialParam::matches(const PairPotentialParam &p) const
{
  return are_approximately_equal(c2, p.c2) &&
    are_approximately_equal(c4, p.c4) &&
    are_approximately_equal(c6, p.c6) &&
    are_approximately_equal(c8, p.c8);
}

PairPotential::PairPotential() : ener(0), covp(0), noncovp(0), nlist(0), covalent(0) { }

PairPotential::~PairPotential()
{
  if (covp)
    const_ptr_t_array2D_delete(covp);
  if (noncovp)
    const_ptr_t_array2D_delete(noncovp);
  if (nlist)
    delete nlist;
  covp = noncovp = 0;
  nlist = 0;
}

void PairPotential::init(MSys *m)
{
  Potential::init(m);
  const BoundaryConditions &bc = *m->boundary_conditions;
  const int nat = m->atom.size();
  int ntype = 0;
  HTab<int> ti;
  type.resize(nat);
  if (covalent)
    PairBinary_delete(covalent);
  covalent = PairBinary_new(nat);
  for (int i = 0; i < nat; i++) {
    const char *t = m->atom[i].type;
    if (!ti.exists(t))
      ti[t] = ntype++;
    type[i] = ti[t];
    LstItr<int> j;
    for (j.init(m->atom[i].neighbors); j.ok(); j.next())
      PairBinary_add(covalent,i,j());
    LstPairItr<int> k;
    for (k.init(m->atom[i].neighbors); k.ok(); k.next())
      PairBinary_add(covalent,k.i(),k.j());
  }
  PairBinary_init(covalent);
  r.resize(nat);
  m->copy_positions_to(r);
  if (bc.type != non_periodic) {
    if (nlist)
      delete nlist;
    rc = 0.45 * bc.min_diameter();
    const double nlw = 0.499 * bc.min_diameter() - rc;
    nlist = new NeighborList(nat, rc, nlw);
    nlist->set_boundary_conditions(bc);
  } else {
    rc = -1;
  }
  if (noncovp)
    const_ptr_t_array2D_delete(noncovp);
  noncovp = const_ptr_t_array2D_new(ntype,ntype);
  if (covp)
    const_ptr_t_array2D_delete(covp);
  covp = const_ptr_t_array2D_new(ntype,ntype);
  ConstHTabIterator<int> j, k;
  for (j.init(ti); j.ok(); j.next())
    for (k.init(ti); k.ok(); k.next()) {
      const Str2 k1(j.key(), k.key());
      const Str2 k2(k.key(), j.key());
      /* Non-covalent parameters */
      const PairPotentialParam *p1 = non_covalent_param.get(k1);
      const PairPotentialParam *p2 = non_covalent_param.get(k2);
      if (!(p1 || p2))
	Out() << "*** PairPotential::init: no non-covalent parameters for types " 
	      << j.key() << " " << k.key() << "\n";
      if (p1 && p2 && !p1->matches(*p2))
	Out() << "*** PairPotential::init: conflicting non-covalent parameters for types " 
	      << j.key() << " " << k.key() << "\n";
      const int v1 = j.val();
      const int v2 = k.val();
      noncovp[v1][v2] = noncovp[v2][v1] = p1 ? p1 : p2;
      /* Covalent parameteres */
      p1 = covalent_param.get(k1);
      p2 = covalent_param.get(k2);
      if (!(p1 || p2))
	Out() << "*** PairPotential::init: no covalent parameters for types " 
	      << j.key() << " " << k.key() << "\n";
      if (p1 && p2 && !p1->matches(*p2))
	Out() << "*** PairPotential::init: conflicting covalent parameters for types " 
	      << j.key() << " " << k.key() << "\n";
      covp[v1][v2] = covp[v2][v1] = p1 ? p1 : p2;
    }
}

void PairPotential::pair_function(int i, int j, double r2, double &u, double &du)
{
  const int ti = type[i];
  const int tj = type[j];
  const double r4 = r2*r2;
  const double r6 = r4*r2;
  const double r8 = r6*r2;
  const double r10 = r8*r2;
  const PairPotentialParam *p = (const PairPotentialParam *)
    (PairBinary_exists(covalent,i,j) ? covp : noncovp)[ti][tj];
  u = p->c2/r2 + p->c4/r4 + p->c6/r6 + p->c8/r8;
  du = -p->c2/r4 - 2*p->c4/r6 - 3*p->c6/r8 - 4*p->c8/r10;
}

void PairPotential::pair_function(int i, int j, double r2, double &u, double &du, double &d2u)
{
  const int ti = type[i];
  const int tj = type[j];
  const double r4 = r2*r2;
  const double r6 = r4*r2;
  const double r8 = r6*r2;
  const double r10 = r8*r2;
  const double r12 = r10*r2;
  const PairPotentialParam *p = (const PairPotentialParam *)
    (PairBinary_exists(covalent,i,j) ? covp : noncovp)[ti][tj];
  u = p->c2/r2 + p->c4/r4 + p->c6/r6 + p->c8/r8;
  du = -p->c2/r4 - 2*p->c4/r6 - 3*p->c6/r8 - 4*p->c8/r10;
  d2u = 2*p->c2/r6 + 6*p->c4/r8 + 12*p->c6/r10 + 20*p->c8/r12;
}

void PairPotential::add_to_energy_and_forces(double &ener, Cartesian *f, int i, int j,
					     const BoundaryConditions &bc, double rc2)
{
  Cartesian rij;
  bc.minimum_image_displacement(rij,r[i],r[j]);
  const double r2 = rij.sq();
  if (rc2 > 0 && r2 > rc2)
    return;
  double u, du;
  pair_function(i,j,r2,u,du);
  ener += u;
  rij *= 2*du;
  f[i] -= rij;
  f[j] += rij;
}

static void add_to_hessian(double **h, int i, int j, const Tensor &t)
{
  i *= 3;
  j *= 3;
  h[i][j]     += t.xx;
  h[i][j+1]   += t.xy;
  h[i][j+2]   += t.xz;
  h[i+1][j]   += t.yx;
  h[i+1][j+1] += t.yy;
  h[i+1][j+2] += t.yz;
  h[i+2][j]   += t.zx;
  h[i+2][j+1] += t.zy;
  h[i+2][j+2] += t.zz;
}

void PairPotential::add_to_hessian(double **h, int i, int j,
				   const BoundaryConditions &bc, double rc2)
{
  Cartesian rij;
  bc.minimum_image_displacement(rij,r[i],r[j]);
  const double r2 = rij.sq();
  if (rc2 > 0 && r2 > rc2)
    return;
  double u, du, d2u;
  pair_function(i,j,r2,u,du,d2u);
  Tensor t = 4*d2u*Tensor(rij,rij) + 2*Tensor(du);
  ::add_to_hessian(h,i,i,t);
  ::add_to_hessian(h,j,j,t);
  t *= -1.0;
  ::add_to_hessian(h,i,j,t);
  ::add_to_hessian(h,j,i,t);
}

void PairPotential::add_to_energy_and_forces(double &u, Cartesian *f)
{
  const BoundaryConditions &bc = *msys->boundary_conditions;
  const int nat = msys->atom.size();
  msys->copy_positions_to(r);
  ener = 0;
  if (nlist) {
    nlist->coordinates_changed(r);
    const double rc2 = sq(rc);
    for (int i = 0; i < nat; i++) {
      const int nj = nlist->number_of_neighbors()[i];
      const int *j = nlist->neighbors()[i];
      for (int k = 0; k < nj; k++)
	add_to_energy_and_forces(ener,f,i,j[k],bc,rc2);
    }
  } else {
    for (int i = 0; i < nat; i++)
      for (int j = i+1; j < nat; j++)
	add_to_energy_and_forces(ener,f,i,j,bc,-1);
  }
  u += ener;
}

void PairPotential::write() const
{
  Out() << "Pair potential energy: " << ener << " kcal/mol\n";
}

bool PairPotential::can_calculate_hessian() const
{
  return true;
}

void PairPotential::add_to_hessian(double **h)
{
  const BoundaryConditions &bc = *msys->boundary_conditions;
  const int nat = msys->atom.size();
  msys->copy_positions_to(r);
  if (nlist) {
    nlist->coordinates_changed(r);
    const double rc2 = sq(rc);
    for (int i = 0; i < nat; i++) {
      const int nj = nlist->number_of_neighbors()[i];
      const int *j = nlist->neighbors()[i];
      for (int k = 0; k < nj; k++)
	add_to_hessian(h,i,j[k],bc,rc2);
    }
  } else {
    for (int i = 0; i < nat; i++)
      for (int j = i+1; j < nat; j++)
	add_to_hessian(h,i,j,bc,-1);
  }
}
