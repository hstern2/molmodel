#ifndef NBRLIST_H
#define NBRLIST_H

#include "cvec.hpp"

class BoundaryConditions;

class NeighborList
{

public:
  
  NeighborList(int n, double cutoff, double width);
  ~NeighborList();
  
  void set_boundary_conditions(const BoundaryConditions &);
  void coordinates_changed(const CVec r);
  void show_self() const;
  
  int * const * neighbors() const; /* return a pointer to list of neighbors for atom i */
  const int *number_of_neighbors() const;
  const int *start(int i) const;
  const int *end(int i) const;
  
private:
  
  int first_time;
  const int n, ncell, ncell3;
  const double cutoff, width;
  CVec last;
  BoundaryConditions *bc, *bccell;
  struct PairBinary *p, *cell;
  int *cell_head, *cell_list;
  
  inline void add_neighbor(int i, int j);
  void rebuild(const Cartesian *r);
  void rebuild_periodic(const Cartesian *r);
  void rebuild_nonperiodic(const Cartesian *r);
  
};

#endif /* NBRLIST_H */
