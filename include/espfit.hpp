#ifndef ESPFIT_H
#define ESPFIT_H

#include "msys.hpp"
#include "cvec.hpp"
#include "dvec.hpp"

struct ESPFit
{
  ESPFit(MSys &msys);
  MSys &msys; // io
  double net_charge; // io
  CVec gridpoints; // io
  DVec potential; // io
  Str write_charges_to_file; // io
  void fit();
  
  classIO(ESPFit);
};

#endif /* ESPFIT_H */
