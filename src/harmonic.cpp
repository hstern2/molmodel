#include "harmonic.hpp"
#include "out.hpp"
#include "timing.h"
#include "linalg.hpp"
#include "units.h"
#include "geograd.hpp"

Harmonic::Harmonic():
  ignore_dihedrals_with_variance_greater_than_pi2_over_4(true),
  temperature(298.15),  
  verbose(1),
  E0(0),
  nfreedom(0),
  u(0)
{ }

Harmonic::~Harmonic() 
{
  cleanup();
}

void Harmonic::cleanup()
{
  nfreedom = 0;
  u = 0;
}

void Harmonic::add_to_energy_and_forces(double &utmp, Cartesian *f )
{
  if (covariance_matrix.rows() == 0)
    return;
  TIMESTART("Harmonic::energy_and_forces");
  const Atom *a = msys->atom;
  /* Compute delta/jacobian factor for each internal coordinate */
  double phi;
  int n = 0, m = 0;
  Cartesian a0,a1,a2,a3;
  
  /* Loop over all zmatrix entries,
     keep track of which one (with 'm')
  */
  for (m = 0; m < zmatrix.size(); m++) {
    if (m>0){
      delta_position[n] = a[zmatrix[m].i].position.distance(a[zmatrix[m].a].position) - average_position[n]; 
      a0 = a[zmatrix[m].i].position - a[zmatrix[m].a].position;
      a0.scale_to_unit_magnitude();
      jacobian_matrix (n, 3*zmatrix[m].i) = a0.x;
      jacobian_matrix (n, 3*zmatrix[m].i+1) = a0.y;
      jacobian_matrix (n, 3*zmatrix[m].i+2) = a0.z;
      jacobian_matrix (n, 3*zmatrix[m].a) = -a0.x;
      jacobian_matrix (n, 3*zmatrix[m].a+1) = -a0.y;
      jacobian_matrix (n, 3*zmatrix[m].a+2) = -a0.z;
      n++;
     }
    if (m>1){
      delta_position[n] = AngleGradient(a[zmatrix[m].i].position,
					a[zmatrix[m].a].position,
					a[zmatrix[m].b].position,a0,a1,a2) - average_position[n];
      jacobian_matrix (n, 3*zmatrix[m].i) = a0.x;
      jacobian_matrix (n, 3*zmatrix[m].i+1) = a0.y;
      jacobian_matrix (n, 3*zmatrix[m].i+2) = a0.z;
      jacobian_matrix (n, 3*zmatrix[m].a) = a1.x;
      jacobian_matrix (n, 3*zmatrix[m].a+1) = a1.y;
      jacobian_matrix (n, 3*zmatrix[m].a+2) = a1.z;
      jacobian_matrix (n, 3*zmatrix[m].b) = a2.x;
      jacobian_matrix (n, 3*zmatrix[m].b+1) = a2.y;
      jacobian_matrix (n, 3*zmatrix[m].b+2) = a2.z;
      n++;
     }
  
    if (m>2 && !ignore_dihedral.exists(m)) {
      phi =  DihedralGradient(a[zmatrix[m].i].position,
			      a[zmatrix[m].a].position,
			      a[zmatrix[m].b].position,
			      a[zmatrix[m].c].position,a0,a1,a2,a3);
      jacobian_matrix (n, 3*zmatrix[m].i) = a0.x;
      jacobian_matrix (n, 3*zmatrix[m].i+1) = a0.y;
      jacobian_matrix (n, 3*zmatrix[m].i+2) = a0.z;
      jacobian_matrix (n, 3*zmatrix[m].a) = a1.x;
      jacobian_matrix (n, 3*zmatrix[m].a+1) = a1.y;
      jacobian_matrix (n, 3*zmatrix[m].a+2) = a1.z;
      jacobian_matrix (n, 3*zmatrix[m].b) = a2.x;
      jacobian_matrix (n, 3*zmatrix[m].b+1) = a2.y;
      jacobian_matrix (n, 3*zmatrix[m].b+2) = a2.z;
      jacobian_matrix (n, 3*zmatrix[m].c) = a3.x;
      jacobian_matrix (n, 3*zmatrix[m].c+1) = a3.y;
      jacobian_matrix (n, 3*zmatrix[m].c+2) = a3.z;
      delta_position[n] = periodic(phi - average_position[n], 2*M_PI);
      n++;
    }
  }
  
  utmp += (u = E0 + 0.5 * (delta_position * harmonic_matrix * delta_position));
  
  DVec tmp = delta_position * harmonic_matrix * jacobian_matrix;
  const double *tmp2 = tmp;
  for (int i = 0; i < msys->atom.size(); i++)
    f[i] -= ((const Cartesian *) tmp2)[i];
  
  TIMESTOP("Harmonic::energy_and_forces");
}

void Harmonic::write() const
{
  if (covariance_matrix.rows() == 0)
    return; // user didn't specify anything
  if (verbose > 2)
    Out() << "Zmatrix coordinates:\n" << zmatrix <<"\n"
	  << "Delta position: " << delta_position << "\n"
	  << "Covariance_matrix: " << covariance_matrix <<"\n";
  
  Out() << "Harmonic energy: " << u << " kcal/mol\n";
}

void Harmonic::detect_dihedrals_to_ignore()
{
  int n = 0, m;
  const int nold = covariance_matrix.rows();
  Vec<bool> keep(nold);
  for (m = 0; m < zmatrix.size(); m++) {
    if (m > 0) {
      keep[n] = true;
      n++; // bond
    }
    if (m > 1) {
      keep[n] = true;
      n++; // angle
    }
    if (m > 2 && !ignore_dihedral.exists(m)) {
      const double sig2 = covariance_matrix(n,n);
      if (sig2 > M_PI*M_PI/4.0) {
	ignore_dihedral.add(m);
	keep[n] = false;
      } else {
	keep[n] = true;
      }
      n++;
    }
  }
  insist(n == nold);
  int nnew = 0;
  for (m = 0; m < nold; m++)
    if (keep[m])
      nnew++;

  DVec tmp(nnew);
  int j = 0;
  for (int i = 0; i < nold; i++)
    if (keep[i]){
      tmp[j] = average_position[i];
      j++;
    }
  insist (j == nnew);
  average_position = tmp;
  DMat tmp2(nnew,nnew);
  int k = 0;
  for (int i = 0; i < nold; i++)
    if (keep[i]) {
      int l = 0;
      for (int j = 0; j < nold; j++) 
	if (keep[j]) {
	  tmp2(k,l) = covariance_matrix(i,j);
	  l++;
	}
      insist(l == nnew);
      k++;
    }
  insist(k == nnew);
  covariance_matrix = tmp2;
}

void Harmonic::init(MSys *m)
{ 
  cleanup();
  Potential::init(m);  
  if (covariance_matrix.rows() == 0)
    return; // user didn't specify anything
  if (ignore_dihedrals_with_variance_greater_than_pi2_over_4)
    detect_dihedrals_to_ignore();
  harmonic_matrix = covariance_matrix.copy();
  InvertSVD(harmonic_matrix);
  harmonic_matrix *= K_to_kcal_mol(temperature);
  nfreedom = covariance_matrix.rows();
  if ( average_position.size() != nfreedom )
    die("Harmonic::init: nfreedom = %d, but average_position contains %d coordinates", 
	nfreedom, average_position.size());
  // here we assume the degrees of freedom we excluded are all dihedrals
  if (ignore_dihedral.size() != 3 * m->atom.size() - 6 - nfreedom)
    die ("Harmonic::init: internal coordinate size exceeds total degree of freedom. covariance matrix size: %d,ignored dihedral is %d,but total freedom is %d \n", nfreedom, ignore_dihedral.size(),3 * m->atom.size() - 6);
  double logz0, logz2, logz4;
  zmatrix.from_internal(average_position,ignore_dihedral);
  zmatrix.harmonic_partition_function(average_position,covariance_matrix,ignore_dihedral,logz0,logz2,logz4);
  const double kT = K_to_kcal_mol(temperature);
  Out() << "Harmonic approximation to free energy...\n"
        << "Temperature: " << temperature << " K\n"
	<< "Energy at average position: " << E0 << " kcal/mol\n"
	<< "With zero-order correction: " << E0 - kT*logz0 << " kcal/mol\n"    
	<< "With second-order correction: " << E0 - kT*(logz0+logz2) << " kcal/mol\n"
	<< "With fourth-order correction: " << E0 - kT*(logz0+logz4) << " kcal/mol\n"
	<< "\n" << flush;
  delta_position.resize(nfreedom);  
  jacobian_matrix.resize(nfreedom,3*msys->atom.size());
  jacobian_matrix.zero();
}
