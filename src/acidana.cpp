#include "acidana.hpp"
#include "msys.hpp"
#include "acid.hpp"

void AcidAnalysis::init(const MSys &msys, int write_interval)
{
  Callback::init(msys,write_interval);
  msysp = &msys;
  const int nacid = msys.acid.size();
  istate = Vec<Vec<int> >(nacid);
  for (int i = 0; i < nacid; i++)
    istate[i] = Vec<int>(msys.acid[i].desc->number_of_states(), 0);
  np.resize(0);
}

void AcidAnalysis::update(const MSys &msys)
{
  assert(msysp = &msys);
  int n = 0;
  for (int i = 0; i < msys.acid.size(); i++) {
    istate[i][msys.acid[i].state]++;
    n += msys.acid[i].number_of_protons();
  }
  if (np.size() < n+1)
    np.resize(n+1, 0);
  np[n]++;
}

void AcidAnalysis::write() const
{
  int i;
  const MSys &msys = *msysp;
  for (i = 0; i < msys.acid.size(); i++) {
    Out() << "Protonation state of group " << i << ": " << msys.acid[i].state << "\n";
    const int ns = istate[i].size();
    int k, n=0;
    for (k = 0; k < ns; k++)
      n += istate[i][k];
    if (n > 0)
      for (k = 0; k < ns; k++)
	Out() << "Fraction of time group " << i << " spent in state " << k << ": "
	      << (double) istate[i][k] / (double) n << "\n";
  }
  int ntot = 0;
  for (i = 0; i < np.size(); i++)
    ntot += np[i];
  if (ntot > 0)
    for (i = 0; i < np.size(); i++)
      Out() << "Fraction of time with " << i << " labile protons: " 
	    << (double) np[i] / (double) ntot << "\n";
}
