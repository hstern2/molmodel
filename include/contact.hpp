#ifndef CONTACT_H
#define CONTACT_H

#include "pot.hpp"
#include "cvec.hpp"
#include "ptab.hpp"
#include "dvec.hpp"
#include "callback.hpp"

class NeighborList;
struct PairBinary;

class Contact : public Potential
{
public:
  HTab<double> number_of_contacts; // io
  Str average_ncontact_file; // io
  MSys native_structure; // io
  int verbose, wait; // io
  double contact_radius, contact_width, stability, ncontact_weight; // io
  void init(MSys *);
  void add_to_energy_and_forces(double &, Cartesian *);
  void write() const;
  Contact();
  ~Contact();
  classIO(Contact);
private:
  CVec r;
  NeighborList *nlist;
  struct PairBinary *partial_contact;
  PTab<double> native;
  DVec ncontact, ncontact0, dncontact, tmpncontact, sumcontact, sumcontact2;
  int nnative, nstat, nstep;
  double ener_native, ener_ncontact;
  void cleanup(), update_avg();
};

#endif /* CONTACT_H */
