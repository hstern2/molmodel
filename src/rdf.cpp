#include "rdf.hpp"
#include "msys.hpp"
#include "fns.h"
#include "timing.h"

void RDF::init(const MSys &msys, int write_interval)
{
  Callback::init(msys,write_interval);
  if (msys.boundary_conditions->type == non_periodic)
    die("RDF::init: boundary conditions must be periodic");
  LstItrMod<RDFPair> p;
  for (p.init(pair); p.ok(); p.next())
    p().init(msys);
}

void RDF::update(const MSys &msys)
{
  TIMESTART("RDF::update");
  LstItrMod<RDFPair> p;
  for (p.init(pair); p.ok(); p.next())
    p().update(msys);
  TIMESTOP("RDF::update");
}

void RDF::write() const
{
  LstItr<RDFPair> p;
  for (p.init(pair); p.ok(); p.next())
    p().write();
}

RDFPair::RDFPair() : 
  number_of_bins(200), peak_location(0), peak_height(0),
  include_intramolecular_contribution(0)
{ }

void RDFPair::init(const MSys &msys)
{
  b1 = Vec<Boolean>(msys.atom.size(), Boolean(0));
  b2 = Vec<Boolean>(msys.atom.size(), Boolean(0));
  n1 = 0;
  n2 = 0;
  Boolean at_least_one1(0);
  Boolean at_least_one2(0);
  int any1 = (strlen(type1) == 0);
  int any2 = (strlen(type2) == 0);
  for (int i = 0; i < msys.atom.size(); i++) {
    if (any1 || !strcmp(msys.atom[i].type, type1)) {
      n1++;
      b1[i] = at_least_one1 = Boolean(1);
    }
    if (any2 || !strcmp(msys.atom[i].type, type2)) {
      n2++;
      b2[i] = at_least_one2 = Boolean(1);
    }
  }
  if (!at_least_one1)
    Out() << "*** RDFPair::init: no atoms of type " << type1 << "\n";
  if (!at_least_one2)
    Out() << "*** RDFPair::init: no atoms of type " << type2 << "\n";
  bin = Vec<int>(number_of_bins, 0);
  ntot = 0;
  dr = 0.5*msys.boundary_conditions->min_diameter()/number_of_bins;
  if (is_almost_zero(dr))
    die("RDFPair::init(): maximum is zero");
  // vol is the average reciprocal volume
  volsum = vol = 0; 
  nvol = 0;
}

void RDFPair::count(int iat, int jat, Vec<Atom> a, const BoundaryConditions &bc)
{
  if (b1[jat] || b2[jat]) {
    const int index = int(floor(bc.minimum_image_distance(a[iat].position, a[jat].position)/dr));
    if (index < number_of_bins) {
      if (b1[iat] && b2[jat])
	bin[index]++;
      if (b2[iat] && b1[jat])
	bin[index]++;
    }
  }
}

void RDFPair::update(const MSys &msys)
{
  const BoundaryConditions &bc = *msys.boundary_conditions;
  ntot++;
  for (int imol = 0; imol < msys.molecule.size(); imol++)
    for (int iat = msys.molecule[imol].a; iat < msys.molecule[imol].b; iat++)
      if (b1[iat] || b2[iat]) {
	for (int jmol = imol+1; jmol < msys.molecule.size(); jmol++)
          for (int jat = msys.molecule[jmol].a; jat < msys.molecule[jmol].b; jat++)
	    count(iat, jat, msys.atom, bc);
	if (include_intramolecular_contribution)
	  for (int jat = iat+1; jat < msys.molecule[imol].b; jat++)
	    count(iat, jat, msys.atom, bc);
      }
  volsum += 1.0/bc.volume();
  nvol++;
  vol = volsum/nvol; // average reciprocal volume
  /* Find the highest peak */
  peak_location = 0;
  peak_height = 0;
  double norm = 4.0 * M_PI * n2 * n1 * ntot * vol / 3.0;
  for (int i = 0; i < number_of_bins; i++) {
    double r1 = i*dr;
    double r2 = r1 + dr;    
    double nbin = bin[i];
    double z = nbin/(norm*(cube(r2)-cube(r1)));
    if (z > peak_height) {
      peak_height = z;
      peak_location = (r1+r2)/2.0;
    }
  }
}

void RDFPair::write() const
{
  if (strlen(file) > 0)
    write_as_gnuplot();
  if (peak_location > 0)
    Out() << "RDF peak location for " << type1 << ", " << type2 
	  << ": " << peak_location << " angstrom\n" << flush
	  << "RDF peak height for " << type1 << ", " << type2 
	  << ": " << peak_height << "\n" << flush;
}

void RDFPair::write_as_gnuplot() const
{
  FILE *f;
  f = fopen(file, "w");
  if (!f)
    die("RDFPair::write_as_gnuplot(): cannot create %s", (const char *) file);
  fprintf(f, "# RDF for types %s and %s\n", (const char *) type1, (const char *) type2);
  fprintf(f, "# Peak location: %16.8f angstrom\n", peak_location);
  fprintf(f, "# Peak height: %16.8f angstrom\n", peak_height);
  fprintf(f, "# %16s%16s\n", "r (angstrom)", "g(r)");
  fprintf(f, "  %16.8f%16.8f\n", 0.0, 0.0);
  // vol is the average reciprocal volume
  double c = 4.0 * M_PI * n2 * n1 * ntot * vol / 3.0;
  for (int i = 0; i < number_of_bins; i++) {
    double r1 = i*dr;
    double r2 = r1 + dr;
    double nbin = bin[i];
    fprintf(f, "  %16.8f%16.8f\n", (r1+r2)/2.0, nbin/(c*(cube(r2)-cube(r1))));
  }
  fclose(f);
}
