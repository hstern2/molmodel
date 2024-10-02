#ifndef AVGMF_H
#define AVGMF_H

#include "msys.hpp"
#include "callback.hpp"
#include "cvec.hpp"

class AvgMolecularForce : public Callback // io
{
public:
  void init(const MSys &msys, int write_interval);
  void update(const MSys &msys);
  void write() const;
  MSys molsys; // io
private:
  CVec fsum;
  int nmol, n;
  classIO(AvgMolecularForce);
};

#endif /* AVGMF_H */
