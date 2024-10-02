#include <stdlib.h>
#include <math.h>
#include "random.h"

void RandomSeed(int s) 
{ 
  srandom(s);
}

double UniformRandom() 
{ 
  return (double) random() / 2147483648.0;
}

static int have_one_saved = 0;
static double saved;

double GaussianRandom()
{
  if (have_one_saved) {
    have_one_saved = 0;
    return saved;
  } else {
    double a = sqrt(-2*log(UniformRandom()));
    double b = 2*M_PI*UniformRandom();
    have_one_saved = 1;
    saved = a*cos(b);
    return a*sin(b);
  }
}

double ExponentialRandom()
{
  return -log(UniformRandom());
}
