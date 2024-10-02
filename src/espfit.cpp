#include "eleclj.hpp"
#include "espfit.hpp"
#include "dmat.hpp"
#include "linalg.hpp"

ESPFit::ESPFit(MSys &m) : 
  msys(m), 
  net_charge(0.0)
{ }

void ESPFit::fit()
{
  Out() << "Starting ESPFit::fit...\n"
	<< "Net charge: " << net_charge << " e\n";
  int i;
  const int nat = msys.atom.size();
  Atom *at = msys.atom;
  const int ngrid = gridpoints.size();
  insist(gridpoints.size() == potential.size());
  char buf[256];
  VecStr saved_types(nat);
  for (i = 0; i < nat; i++) {
    snprintf(buf, 256, "%d", i);
    saved_types[i] = at[i].type;
    at[i].type = buf;
  }
  DMat A(ngrid,nat);
  DVec p;
  ElecLJ elec;
  for (i = 0; i < nat; i++) {
    snprintf(buf, 256, "%d", i);
    elec.fixed_charge[buf] = charge_unit_to_e(1.0);
    if (i == 0)
      elec.init(&msys);
    elec.types_changed();
    elec.electrostatic_potential_at_gridpoints(gridpoints,p);
    for (int j = 0; j < ngrid; j++)
      A(j,i) = p[j];
    elec.fixed_charge.remove(buf);
  }
  DVec c(ngrid);
  for (i = 0; i < ngrid; i++)
    c[i] = hartree_e_to_potential_unit(potential[i]);
  DMat B(1,nat,1.0);
  DVec d(1, e_to_charge_unit(net_charge));
  DVec x(nat);
  ConstrainedLeastSquares(A, B, c, d, x);
  for (i = 0; i < nat; i++)
    x[i] = charge_unit_to_e(x[i]);
  if (strlen(write_charges_to_file) > 0) {
    ostream *s = FileStream(write_charges_to_file);
    *s << x;
    delete s;
  }
  OutPrintf("ESP charges on atoms:\n%-5s %12s %12s %12s %12s\n","sym","x","y","z","charge (e)");
  for (i = 0; i < nat; i++) {
    OutPrintf("%-5s %12.8f %12.8f %12.8f %12.8f\n", 
	      (const char *) at[i].symbol,
	      at[i].position.x,
	      at[i].position.y,
	      at[i].position.z,
	      x[i]);
    at[i].type = saved_types[i];
  }
}
