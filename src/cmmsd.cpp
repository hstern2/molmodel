#include "cmmsd.hpp"

void CMMeanSquareDistance::init(const MSys &msys, int write_interval)
{
  Callback::init(msys,write_interval);
  const int n = msys.molecule.size();
  c0.resize(n);
  for (int i = 0; i < n; i++)
    c0[i] = msys.molecule_center_of_mass(i);
  msd = 0;
}

void CMMeanSquareDistance::update(const MSys &msys)
{
  const int n = msys.molecule.size();
  msd = 0;
  for (int i = 0; i < n; i++)
    msd += (c0[i] - msys.molecule_center_of_mass(i)).sq();
  if (n > 0)
    msd /= n;
}

void CMMeanSquareDistance::write() const
{
  Out() << "Center of mass mean square distance: " << msd << " A^2\n";
}
