#include "nbrlist.hpp"
#include "bcond.hpp"
#include "pbin.h"
#include "timing.h"
#include "fns.h"

#define PARTICLES_PER_CELL 2

NeighborList::NeighborList(int n_, double cutoff_, double width_)
  : first_time(1),
    n(n_), ncell((int) floor(cbrt(n/PARTICLES_PER_CELL+3))),
    ncell3(ncell*ncell*ncell),
    cutoff(cutoff_), width(width_), 
    last(n), bc(0), bccell(0),
    p(PairBinary_new(n)), cell(PairBinary_new(ncell3)),
    cell_head(new int[ncell3]), cell_list(new int[n])
{ }

NeighborList::~NeighborList()
{
  if (bc)
    delete bc;
  if (bccell)
    delete bccell;
  PairBinary_delete(p);
  PairBinary_delete(cell);
  delete[] cell_head;
  delete[] cell_list;
}

int * const * NeighborList::neighbors() const
{
  return PairBinary_elements(p);
}

const int *NeighborList::number_of_neighbors() const
{
  return PairBinary_number_of_elements(p);
}

const int *NeighborList::start(int i) const
{
  return PairBinary_start(p,i);
}

const int *NeighborList::end(int i) const
{
  return PairBinary_end(p,i);
}

void NeighborList::coordinates_changed(const CVec r)
{
  TIMESTART("NeighborList::coordinates_changed");
  insist(bc);
  insist(r.size() == n);
  insist(last.size() == n);
  insist(cutoff > 0);
  if (first_time) {
    first_time = 0;
    TIMESTOP("NeighborList::coordinates_changed");
    rebuild(r);
    return;
  }
  double largest = 0, second_largest = 0;
  for (int i = 0; i < n; i++) {
    const double x = r[i].distance(last[i]);
    if (x > largest) {
      second_largest = largest;
      largest = x;
    } else if (x > second_largest) {
      second_largest = x;
    }
  }
  TIMESTOP("NeighborList::coordinates_changed");
  if (second_largest + largest > width)
    rebuild(r);
}

static inline int cell_index(int ix, int iy, int iz, int n)
{
  return (mymod(ix,n)*n + mymod(iy,n))*n + mymod(iz,n);
}

void NeighborList::set_boundary_conditions(const BoundaryConditions &bc_)
{
  TIMESTART("NeighborList::set_boundary_conditions");
  /* Create a neighbor list for the cells themselves */
  if (bc)
    delete bc;
  if (bccell)
    delete bccell;
  bc = bc_.copy();
  if (bc->type == non_periodic)
    return;
  insist(cutoff > 0);
  bccell = bc_.copy();
  bccell->scale(1.0/ncell);
  PairBinary_clear(cell);
  const Tensor bclat = bccell->lattice_vectors();
  const double c2 = sq(cutoff + width + bccell->max_diameter());
  BoundaryConditions *bctri = BoundaryConditions::new_triclinic(bclat);
  const int nmax = min((int) ceil((cutoff + width + bccell->max_diameter())/bctri->min_diameter()), ncell/2+1);
  delete bctri;
  bctri = 0;
  for (int jx = 0; jx <= nmax; jx++)
    for (int jy = (jx == 0 ? 0 : -nmax); jy <= nmax; jy++)
      for (int jz = (jx == 0 && jy == 0 ? 1 : -nmax); jz <= nmax; jz++) {
	Cartesian jr = bclat*Cartesian(jx,jy,jz);
	bccell->map_to_central_box(jr);
	if (jr.sq() < c2)
	  for (int ix = 0; ix < ncell; ix++)
	    for (int iy = 0; iy < ncell; iy++)
	      for (int iz = 0; iz < ncell; iz++) {
		const int ia = cell_index(ix,iy,iz,ncell);
		const int ib = cell_index(ix+jx,iy+jy,iz+jz,ncell);
		if (ia != ib)
		  PairBinary_add(cell, ia, ib);
	      }
      }
  PairBinary_init(cell);
  first_time = 1;
  TIMESTOP("NeighborList::set_boundary_conditions");
}

inline void NeighborList::add_neighbor(int i, int j)
{
  PairBinary_add(p,i,j);
}

void NeighborList::rebuild(const Cartesian *r)
{
  TIMESTART("NeighborList::rebuild");
  PairBinary_clear(p);
  for (int i = 0; i < ncell3; i++)
    cell_head[i] = -1;
  if (bc->type == non_periodic)
    rebuild_nonperiodic(r);
  else
    rebuild_periodic(r);
  last.copy(r);
  TIMESTOP("NeighborList::rebuild");
}

void NeighborList::rebuild_nonperiodic(const Cartesian *r)
{
  if (n == 0)
    return;
  int i;
  insist(cutoff > 0);
  const double c2 = sq(cutoff+width);
  Cartesian rmin(r[0]), rmax(r[0]);
  for (i = 0; i < n; i++) {
    const Cartesian &a = r[i];
    if (a.x < rmin.x)
      rmin.x = a.x;
    else if (a.x > rmax.x)
      rmax.x = a.x;
    if (a.y < rmin.y)
      rmin.y = a.y;
    else if (a.y > rmax.y)
      rmax.y = a.y;
    if (a.z < rmin.z)
      rmin.z = a.z;
    else if (a.z > rmax.z)
      rmax.z = a.z;
  }
  rmin -= Cartesian(1e-2,1e-2,1e-2);
  rmax += Cartesian(1e-2,1e-2,1e-2);
  Cartesian d = rmax - rmin;
  const Cartesian lcell(d.x/ncell, d.y/ncell, d.z/ncell);
  const int ncell2 = ncell*ncell;
  /* Assign particles to cells */
  for (i = 0; i < n; i++) {
    const Cartesian &a = r[i];
    const int icell = (int)
      (floor((a.x - rmin.x)/lcell.x) + 
       floor((a.y - rmin.y)/lcell.y) * ncell +
       floor((a.z - rmin.z)/lcell.z) * ncell2);
    cell_list[i] = cell_head[icell];
    cell_head[icell] = i;
  }
  /* Loop over all pairs of cells */
  const int wcellx = min((int) ceil((cutoff+width)/lcell.x), ncell/2+1);
  const int wcelly = min((int) ceil((cutoff+width)/lcell.y), ncell/2+1);
  const int wcellz = min((int) ceil((cutoff+width)/lcell.z), ncell/2+1);
  const double cell_cutoff2 = sq(cutoff + width + lcell.magnitude());
  for (int ijx = 0; ijx <= wcellx; ijx++) {
    const double dx2 = sq(ijx*lcell.x);
    for (int ijy = (ijx == 0) ? 0 : -wcelly; ijy <= wcelly; ijy++) {
      const double dx2dy2 = dx2 + sq(ijy*lcell.y);
      for (int ijz = (ijx == 0 && ijy == 0) ? 0 : -wcellz; ijz <= wcellz; ijz++) {
	if (dx2dy2 + sq(ijz*lcell.z) > cell_cutoff2)
	  continue;
	const int same_cell = ijx == 0 && ijy == 0 && ijz == 0;
	for (int ix = 0; ix < ncell; ix++) {
	  const int jx = ix + ijx;
	  if (0 <= jx && jx < ncell)
	    for (int iy = 0; iy < ncell; iy++) {
	      const int jy = iy + ijy;
	      if (0 <= jy && jy < ncell)
		for (int iz = 0; iz < ncell; iz++) {
		  const int jz = iz + ijz;
		  if (0 <= jz && jz < ncell) {
		    const int icell = ix + iy*ncell + iz*ncell2;
		    const int jcell = jx + jy*ncell + jz*ncell2;
		    /* Loop over all pairs of particles within this pair of cells */
		    for (int i = cell_head[icell]; i != -1; i = cell_list[i])
		      for (int j = same_cell ? cell_list[i] : cell_head[jcell]; j != -1; j = cell_list[j]) {
			const double x = r[i].x - r[j].x, y = r[i].y - r[j].y, z = r[i].z - r[j].z;
			if (x*x + y*y + z*z < c2)
			  add_neighbor(i,j);
		      }
		  }
		}
	    }
	}
      }
    }
  }
}

void NeighborList::rebuild_periodic(const Cartesian *r)
{
  insist(cutoff > 0);
  const double c2 = sq(cutoff + width);
  /* Assign particles to cells */
  int i;
  for (i = 0; i < n; i++) {
    int ix, iy, iz;
    bccell->nearest_lattice_point(r[i],ix,iy,iz);
    const int icell = cell_index(ix,iy,iz,ncell);
    cell_list[i] = cell_head[icell];
    cell_head[icell] = i;
  }
  /* Loop over pairs of particles in neighboring cells and assign to neighbor list */
  const int *nel = PairBinary_number_of_elements(cell);
  int * const *el = PairBinary_elements(cell);
  for (i = 0; i < ncell3; i++) {
    const int nj = nel[i];
    const int *j = el[i];
    for (int k = cell_head[i]; k != -1; k = cell_list[k]) {
      int l;
      for (l = cell_list[k]; l != -1; l = cell_list[l])
	add_neighbor(k,l); /* particles in the same cell */
      for (int m = 0; m < nj; m++)
	for (l = cell_head[j[m]]; l != -1; l = cell_list[l])
	  if (bc->square_minimum_image_distance(r[k],r[l]) < c2)
	    add_neighbor(k,l); /* particles in neighboring cells */
    }
  }
}

void NeighborList::show_self() const
{
  Out() << "Cell neighbor list:\n" << flush;
  PairBinary_show(cell);
  Out() << "Neighbor list:\n" << flush;
  PairBinary_show(p);
}
