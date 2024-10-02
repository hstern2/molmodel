#include <cmath>
#include <cstring>
#include "out.hpp"
#include "hist.hpp"
#include "fns.h"
#include "random.h"

Histogram::Histogram() 
  : nbin(50), lo(0), hi(-1), is_dynamic(true), ncount(0), v(0)
{ }

Histogram::Histogram(const Histogram &h) 
  : nbin(h.nbin), lo(h.lo), hi(h.hi), is_dynamic(h.is_dynamic), ncount(0), v(0)
{ }

void Histogram::init()
{
  ncount = 0;
  v = Vec<int>(nbin, 0);
}

void Histogram::update(double x)
{
  if (hi < lo) {
    if (!is_dynamic)
      die("Histogram: hi (%f) < lo (%f)", hi, lo);
    const double eps = 1e-8*(fabs(x) + 1.0);
    lo = x - eps;
    hi = x + eps;
  }
  if (x < lo || x >= hi) {
    if (!is_dynamic)
      return;
    /* Resize histogram */
    const double eps = 1e-8*(fabs(x) + 1.0);
    double lotmp = lo, hitmp = hi;
    if (x < lo)
      lo = x - eps;
    else
      hi = x + eps;
    Vec<int> vtmp = v;
    v = Vec<int>(nbin, 0);
    ncount = 0;
    for (int i = 0; i < nbin; i++)
      for (int j = 0; j < vtmp[i]; j++)
	update_v((i+UniformRandom())*(hitmp-lotmp)/double(nbin) + lotmp);
  }
  update_v(x);
}

void Histogram::update_v(double x)
{
  v[int(floor(nbin*(x-lo)/(hi-lo)))]++;
  ncount++;
}

void Histogram::write(const char *f, const char *cmt) const
{
  if (strlen(f) > 0) {
    ostream *s = FileStream(f);
    if (cmt)
      *s << cmt;
    write(*s);
    delete s;
  }
}

void Histogram::write(ostream &s) const
{
  if (ncount == 0)
    return;
  char buf[256];
  snprintf(buf, 256, "# gauss(x,%f,%f)\n", mean(), sqrt(variance()));
  s << buf;
  s << "# " << ncount << "\n";
  Vec<double> vc(v.size());
  int i;
  for (i = 0; i < nbin; i++)
    vc[i] = v[i];
  /* Renormalize */
  int j, k;
  for (i = 0; i < nbin && v[i] == 0; i++)
    ;
  for (j = i+1; j < nbin && v[j] == 0; j++)
    ;
  for (k = j+1; k < nbin && v[k] == 0; k++)
    ;
  while (k < nbin && vc[k] > 0) {
    vc[j] /= 0.5*(k-i);
    i = j;
    j = k;
    for (k++; k < nbin && v[k] == 0; k++)
      ;
  }
  vc[j] /= j-i;
  for (int i = 0; i < nbin; i++)
    s << (i+0.5)*(hi-lo)/double(nbin) + lo << " " 
      << vc[i]*nbin/(ncount*(hi-lo)) << "\n";
}

double Histogram::mean() const
{
  if (ncount == 0)
    return 0;
  double sum = 0;
  for (int i = 0; i < nbin; i++)
    sum += v[i] * ((i+0.5)*(hi-lo)/double(nbin) + lo);
  return sum/ncount;
}

double Histogram::variance() const
{
  if (ncount == 0)
    return 0;
  double sum = 0, m = mean();
  for (int i = 0; i < nbin; i++)
    sum += v[i] * sq((i+0.5)*(hi-lo)/double(nbin) + lo - m);
  return sum/ncount;
}
