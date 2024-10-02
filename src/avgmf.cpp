#include "avgmf.hpp"

void AvgMolecularForce::init(const MSys &msys, int write_interval)
{
  Callback::init(msys, write_interval);
  insist(f);
  n = 0;
  nmol = msys.molecule.size();
  fsum = CVec(nmol, Cartesian(0,0,0));
}

void AvgMolecularForce::update(const MSys &msys)
{
  for (int im = 0; im < nmol; im++) {
    Cartesian fs(0,0,0);
    for (int i = msys.molecule[im].a; i < msys.molecule[im].b; i++, n++)
      fs += f[i];
    fsum[im] += fs;
  }
  n++;   
}

void AvgMolecularForce::write() const
{
  if (n == 0)
    return;
  for (int im = 0; im < nmol; im++)
    Out() << "Average force on molecule " << im << ": " << fsum[im]/n << " kcal/(mol A)\n";
}
