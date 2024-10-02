#include "rmsdcb.hpp"
#include "rmsd.hpp"

void RMSDCallBack::update(const MSys &msys)
{
  insist(msys.atom.size() == molsys.atom.size());
  c1.resize(msys.atom.size());
  c2.resize(msys.atom.size());
  msys.copy_positions_to(c1);
  molsys.copy_positions_to(c2);
  rmsd = RootMeanSquareDistance(c1,c2);
}

void RMSDCallBack::write() const
{
  Out() << "Root mean square distance: " << rmsd << " A\n";
}
