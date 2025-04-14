/* $Id: histoprog.C,v 1.7 2006/08/31 23:22:20 hstern Exp $ */

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "hist.hpp"
#include "lst.hpp"
#include "dmat.hpp"
#include "out.hpp"

int main(int argc, char **argv)
{
  OutInit();
  const char *usage = 
    "usage: hist [-n <nbin>][-l <lo>][-h <hi>]\n\n"
    "Takes a series of numbers on stdin and outputs\n"
    "a histogram on stdout\n\n"
    "Options:\n"
    "-n <nbin>:  number of bins (default 50)\n"
    "-l <lo>: lower boundary of histogram (default lowest value of data)\n"
    "-h <hi>: higher boundary of histogram (default highest value of data)\n";
  int c;
  extern char *optarg;
  int nbin = 50;
  double lo = 0, hi = 0;
  bool lo_undef = true, hi_undef = true;
  while ((c = getopt(argc, argv, "n:l:h:")) != EOF)
    switch (c) {
    case 'n':
      nbin = atoi(optarg);
      break;
    case 'l':
      lo = atof(optarg);
      lo_undef = false;
      break;
    case 'h':
      hi = atof(optarg);
      hi_undef = false;
      break;
    default:
      die(usage);
    }
  Lst<double> da;
  double x;
  while(cin >> x) da.add(x);
  DVec data;
  data.copy_from_list(da);
  if (lo_undef)
    lo = data.minimum() - 1e-8;
  if (hi_undef)
    hi = data.maximum() + 1e-8;
  if (hi < lo)
    die("hist: <hi> is less than <lo>");
  Histogram h;
  h.is_dynamic = false;
  h.nbin = nbin;
  h.lo = lo;
  h.hi = hi;
  h.init();
  for (int i = 0; i < data.size(); i++)
    h.update(data[i]);
  h.write(cout);
}
