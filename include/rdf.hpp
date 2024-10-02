#ifndef RDF_H
#define RDF_H

#include "callback.hpp"
#include "msys.hpp"
#include "boo.hpp"

class MDyn;

class RDFPair 
{
public:
  RDFPair();
  void init(const MSys &msys);
  void count(int iat, int jat, Vec<Atom> a, const BoundaryConditions &bc);
  void update(const MSys &msys);
  void write() const;
  void write_as_gnuplot() const;
private:
  int number_of_bins; // io
  double dr, volsum, vol, peak_location, peak_height;
  Boolean include_intramolecular_contribution; // io
  Vec<Boolean> b1, b2;
  Vec<int> bin;
  int n1, n2, ntot, nvol;
  Str type1, type2, file; // io
  classIO(RDFPair);
};

class RDF : public Callback // io
{
public:
  void init(const MSys &msys, int write_interval);
  void update(const MSys &msys);
  void write() const;
private:
  Lst<RDFPair> pair; // io
  classIO(RDF);
};

#endif /* RDF_H */
