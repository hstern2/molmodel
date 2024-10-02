#ifndef BCIFIT_H
#define BCIFIT_H

#include "msys.hpp"
#include "eleclj.hpp"
#include <gsl/gsl_multimin.h>

/***
 * S(k) = int rho(r) exp(-ikr) d^3r
 * should be given in atomic units
 ***/
struct StructureFactor
{
  Cartesian k; // io   /* A.U. */
  Complex Sk; // io  /* A.U. */
  classIO(StructureFactor);
};

typedef Lst<StructureFactor> LstStructureFactor;

struct MolSpec
{
  MolSpec();

  int number_of_symmetry_copies; // io
  MSys msys; // io
  DVec charges; // io
  Cartesian dipole; // io        /* Debye */
  Quadrupole quadrupole; // io   /* Debye Angstrom */
  Octopole octopole; // io   /* Debye Angstrom^2 */
  Hexadecapole hexadecapole; // io   /* Debye Angstrom^3 */
  LstStructureFactor electronic_structure_factors; // io   /* A.U. */
  Tensor polarizability; // io    /* Angstrom^3 */
  ElecLJ elec; // out

  classIO(MolSpec);
};

typedef Lst<MolSpec> LstMolSpec;

struct BCIFit
{
  BCIFit();
  LstMolSpec molecule_specification; // io
  
  int verbose; // io
  double kspace_cutoff, elec_1_4_scale_factor; // io
  HTab<double,Str3> one_three_interaction; // io
  HTab<double> fixed_charge, screening_radius; // io
  HTab<BondChargeIncrementParam> bci_to_fit, bci_to_fix; // io
  HTab<MSiteSpec,Str3> msite; // io
  HTab<Sp3SiteSpec,Str3> sp3; // io
  HTab<Sp2SiteSpec,Str2> sp2; // io
  double atomic_charge_weight; // io
  double dipole_moment_weight; // io
  double quadrupole_moment_weight; // io
  double octopole_moment_weight; // io
  double hexadecapole_moment_weight; // io
  double structure_factor_weight; // io
  double svd_tolerance; // io
  Boolean use_traceless_multipoles; // io
  Str write_parameters_to; // io

  HTab<BondChargeIncrementParam> bcitab;

  double tolerance, bci_minimum; // io
  int max_iterations; // io

  void fit_bci(); // in
  void fit(); // in
  double help(const gsl_vector *);
  classIO(BCIFit);
};

#endif /* BCIFIT_H */
