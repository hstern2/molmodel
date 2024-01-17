/* $Id: geoprog.C,v 1.4 2006/02/14 18:09:35 hstern Exp hstern2 $ */

#include "coord.h"
#include "units.h"
#include "geograd.h"
#include "out.h"
#include <cstdlib>
#include <cstring>

const char Usage[] = 
"usage: geo [distance <r1> <r2>]\n"
"           [angle <r1> <r2> <r3>]\n"
"           [dihedral <r1> <r2> <r3> <r4>]\n"
"           [zlocation <a> <b> <c> <r> <theta> <phi>]\n";

void die()
{
  printf(Usage);
  exit(1);
}

int main(int argc, char *argv[])
{
  argc--;
  Cartesian a, b, c, d;
  if (argc < 1)
    die();
  const char *which = *++argv;
  if (argc < 7)
    die();
  a.x = atof(*++argv);
  a.y = atof(*++argv);
  a.z = atof(*++argv);
  b.x = atof(*++argv);
  b.y = atof(*++argv);
  b.z = atof(*++argv);
  if (!strcmp(which,"distance")) {
    printf("%12.8f\n", a.distance(b));
    exit(0);
  }
  if (argc < 10)
    die();
  c.x = atof(*++argv);
  c.y = atof(*++argv);
  c.z = atof(*++argv);
  if (!strcmp(which,"angle")) {
    printf("%12.8f\n", radians_to_degrees(Angle(a,b,c)));
    exit(0);
  }
  if (argc < 13)
    die();
  d.x = atof(*++argv);
  d.y = atof(*++argv);
  d.z = atof(*++argv);
  if (!strcmp(which,"dihedral")) {
    printf("%12.8f\n", radians_to_degrees(Dihedral(a,b,c,d)));
    exit(0);
  }
  if (!strcmp(which,"zlocation") || !strcmp(which,"zlocation")) {
    const Cartesian loc = ZLocation(a,b,c,d.x,
				    degrees_to_radians(d.y),
				    degrees_to_radians(d.z));
    printf("%12.8f %12.8f %12.8f\n", loc.x, loc.y, loc.z);
    exit(0);
  }
  die();
}
