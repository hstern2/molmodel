#include "callback.hpp"
#include "geohist.hpp"
#include "rdf.hpp"
#include "acidana.hpp"
#include "rmsdcb.hpp"
#include "cmmsd.hpp"
#include "avgmf.hpp"

Callback::Callback() : wait(0), update_interval(-1), f(0) { }
Callback::~Callback() { }

void Callback::init(const MSys &m, int write_interval) 
{ 
  if (update_interval <= 0) {
    if (write_interval < 50)
      update_interval = write_interval;
    else
      update_interval = 50;
  }
}

void Callback::maybe_update(double time_elapsed, int istep, const MSys &m)
{
  if (update_interval > 0 && istep % update_interval == 0 && time_elapsed >= wait)
    update(m);
}

void Callback::maybe_write(double time_elapsed) const
{
  if (time_elapsed >= wait)
    write();
}

istream & operator>>(istream &s, CallbackPtr &ptr)
{
  Str t;
  char bracket;
  s >> t >> bracket;
  s.putback(bracket);
  if (ptr.p)
    delete ptr.p;
  if (!strcmp(t, "radial_distribution_function")) {
    RDF *p = new RDF;
    if (bracket == '{')
      s >> *p;
    ptr = p;
  } else if (!strcmp(t, "geometry_statistics")) {
    GeometryStatistics *p = new GeometryStatistics;
    if (bracket == '{')
      s >> *p;
    ptr = p;
  } else if (!strcmp(t, "acid_statistics")) {
    AcidAnalysis *p = new AcidAnalysis;
    if (bracket == '{')
      s >> *p;
    ptr = p;
  } else if (!strcmp(t, "center_of_mass_mean_square_distance")) {
    CMMeanSquareDistance *p = new CMMeanSquareDistance;
    if (bracket == '{')
      s >> *p;
    ptr = p;
  } else if (!strcmp(t, "root_mean_square_distance")) {
    RMSDCallBack *p = new RMSDCallBack;
    s >> *p;
    ptr = p;
  } else if (!strcmp(t, "average_molecular_forces")) {
    AvgMolecularForce *p = new AvgMolecularForce;
    if (bracket == '{')
      s >> *p;
    ptr = p;
  } else {
    die("Callback: error reading callback name: %s\n"
	"Must be one of 'radial_distribution_function',"
	"'geometry_statistics','acid_statistics',"
	"'root_mean_square_distance','average_molecular_forces'",
	(const char *) t);
  }
  return s;
}
