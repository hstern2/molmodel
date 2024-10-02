#include "contact.hpp"
#include "nbrlist.hpp"
#include "timing.h"
#include "fns.h"
#include "pbin.h"
#include "random.h"

Contact::Contact() : 
  verbose(1),
  wait(0),
  contact_radius(8),
  contact_width(1),
  stability(1000),
  ncontact_weight(10),
  nlist(0),
  partial_contact(0),
  nnative(0),
  nstat(0),
  nstep(0),
  ener_native(0),
  ener_ncontact(0)
{ }

Contact::~Contact() { cleanup(); }

void Contact::cleanup()
{
  if (nlist)
    delete nlist;
  if (partial_contact)
    PairBinary_delete(partial_contact);
  native.remove_all();
  nlist = 0;
  partial_contact = 0;
}

void Contact::init(MSys *m)
{
  cleanup();
  Potential::init(m);
  const int n = native_structure.atom.size();
  if (n == 0)
    return;
  insist(msys->atom.size() == n);
  r.resize(n);
  ncontact.resize(n);
  ncontact0.resize(n);
  tmpncontact.resize(n);
  dncontact.resize(n);
  sumcontact.resize(n);
  sumcontact.zero();
  sumcontact2.resize(n);
  sumcontact2.zero();
  native_structure.copy_positions_to(r);
  const BoundaryConditions &bc = *native_structure.boundary_conditions;
  nlist = new NeighborList(n, contact_radius + contact_width, 1.0);
  nlist->set_boundary_conditions(bc);
  nlist->coordinates_changed(r);
  partial_contact = PairBinary_new(n);
  nnative = 0;
  nstat = 0;
  nstep = 0;
  for (int i = 0; i < n; i++) {
    double *p = 0;
    if ((p = number_of_contacts.get(msys->atom[i].type)) != 0) {
      ncontact0[i] = *p;
    } else {
      ncontact0[i] = -1;
    }
    tmpncontact[i] = ncontact0[i];
    const int nj = nlist->number_of_neighbors()[i];
    const int *j = nlist->neighbors()[i];
    for (int k = 0; k < nj; k++) {
      const double rij = bc.minimum_image_distance(r[i], r[j[k]]);
      if (rij < contact_radius) {
	native(i,j[k]) = rij;
	nnative++;
      }
    }
  }
  Out() << "Initializing contact potential...\n"
	<< "Number of native contacts: " << nnative << "\n"
	<< "Stability: " << stability << " kcal/mol\n"
	<< "Weight factor for number of contacts: " << ncontact_weight << " kcal/mol\n"
	<< "Contact radius: " << contact_radius << " A\n"
	<< "Contact width: " << contact_width << " A\n"
	<< "\n" << flush;
}

/* Long-range switching function */
static void switching_fn(double r2, double rlo2, double rhi2, double &f, double &_2df)
{
  const double a = 1.0/(rhi2 - rlo2);
  const double z = (r2 - rlo2)*a;
  f = ((-6*z + 15)*z - 10)*z*z*z + 1;
  _2df = ((60*z - 120)*z + 60)*z*z*a;
}

void Contact::add_to_energy_and_forces(double &u, Cartesian *f)
{
  TIMESTART("Contact::energy_and_forces");
  msys->copy_positions_to(r);
  const BoundaryConditions &bc = *msys->boundary_conditions;
  PTabItr<double> p;
  ener_native = 0;
  ener_ncontact = 0;
  const double a = stability/nnative;
  nstep++;
  /* Native contacts */
  for (p.init(native); p.ok(); p.next()) {
    const int i = p.i();
    const int j = p.j();
    Cartesian rij;
    bc.minimum_image_displacement(rij,r[i],r[j]);
    const double r2 = rij.sq();
    const double rm = sqrt(r2);
    const double r0 = p.val();
    const double rm_r0 = rm - r0;
    if (fabs(rm_r0) < contact_width) {
      const double x = rm_r0/contact_width;
      const double x2 = x*x;
      ener_native += 2*x2 - x2*x2 - 1;
      const Cartesian fij = (-4*a*x*(1-x2)/(contact_width*rm)) * rij;
      f[i] += fij;
      f[j] -= fij;
    }
  }
  ener_native *= a;
  const int n = r.size();
  int i;
  /* Number of contacts */
  ncontact.zero();
  nlist->coordinates_changed(r);
  const double cr2 = sq(contact_radius);
  const double crw2 = sq(contact_radius+contact_width);
  PairBinary_clear(partial_contact);
  for (i = 0; i < n; i++) {
    const int nj = nlist->number_of_neighbors()[i];
    const int *j = nlist->neighbors()[i];
    for (int k = 0; k < nj; k++) {
      const int jk = j[k];
      const double r2 = bc.square_minimum_image_distance(r[i], r[jk]);
      if (r2 < crw2) {
	if (r2 < cr2) {
	  ncontact[i] += 1.0;
	  ncontact[jk] += 1.0;
	} else {
	  double fn, _2dfn;
	  switching_fn(r2, cr2, crw2, fn, _2dfn);
	  ncontact[i] += fn;
	  ncontact[jk] += fn;
	  PairBinary_quick_add(partial_contact, i, jk);
	}
      }
    }
  }
  for (i = 0; i < n; i++) {
    if (ncontact0[i] >= 0) {
      const double tmp = tmpncontact[i] + 1e-4*GaussianRandom();
      if (tmp >= 0 && UniformRandom() < exp((tmpncontact[i] - tmp)/ncontact0[i]))
	tmpncontact[i] = tmp;
      dncontact[i] = ncontact[i] - tmpncontact[i];
    } else {
      dncontact[i] = 0;
    }
  }
  for (i = 0; i < n; i++) {
    ener_ncontact += ncontact_weight*sq(dncontact[i]);
    const int nj = PairBinary_number_of_elements(partial_contact)[i];
    const int *j = PairBinary_elements(partial_contact)[i];
    for (int k = 0; k < nj; k++) {
      const int jk = j[k];
      Cartesian rij;
      bc.minimum_image_displacement(rij, r[i], r[jk]);
      const double r2 = rij.sq();
      insist(cr2 <= r2 && r2 <= crw2);
      double fn, _2dfn;
      switching_fn(r2, cr2, crw2, fn, _2dfn);
      const Cartesian fij = (2*ncontact_weight*(dncontact[i]+dncontact[jk])*_2dfn)*rij;
      f[i] += fij;
      f[jk] -= fij;
    }
  }
  u += ener_native + ener_ncontact;
  TIMESTOP("Contact::energy_and_forces");
}

void Contact::update_avg()
{
  for (int i = 0; i < ncontact.size(); i++) {
    sumcontact[i] += ncontact[i];
    sumcontact2[i] += sq(ncontact[i]);
  }
  nstat++;
}

void Contact::write() const
{
  Out() << "Native contact energy: " << ener_native << " kcal/mol\n"
	<< "Restraint to number of contacts: " << ener_ncontact << " kcal/mol\n";
  if (verbose > 1)
    Out() << "Number of contacts:\n" << ncontact
	  << "Deviation from constrained number of contacts:\n" << dncontact;
  Out() << flush;
  if (strlen(average_ncontact_file) > 0 && nstep >= wait) {
    ((Contact *) this)->update_avg();
    ostream *s = FileStream(average_ncontact_file);
    *s << "# avg stddev restrained\n";
    for (int i = 0; i < sumcontact.size(); i++) 
      if (ncontact0[i] > 0) {
	const double avg = sumcontact[i]/nstat;
	*s << i << " " << avg << " " 
	   << sqrt(sumcontact2[i]/nstat - sq(avg))
	   << " " << ncontact0[i]
	   << "  # " << msys->atom[i].type
	   << "\n";
      }
    delete s;
  }
}
