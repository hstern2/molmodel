#include "pbin.h"
#include <stdlib.h>
#include <stdio.h>

#define SLACK 64
#define MAKE_CANONICAL(i,j) if (((i+j)%2) ^ (i<j)) { const int tmp = i; i = j; j = tmp; }

struct PairBinary { int n, *nn, *space, **t; };

/***
 * t[i] is array of neighbors.  
 * nn[i] is actual number of neighbors;
 * space[i] is space available.
 ***/
struct PairBinary *PairBinary_new(int n)
{
  int i;
  struct PairBinary *that = (struct PairBinary *) malloc(sizeof(struct PairBinary));
  that->n = n;
  that->nn = (int *) malloc(n * sizeof(int));
  that->space = (int *) malloc(n * sizeof(int));
  that->t = (int **) malloc(n * sizeof(int *));
  for (i = 0; i < n; i++) {
    that->nn[i] = 0;
    that->t[i] = (int *) malloc((that->space[i] = SLACK) * sizeof(int));
  }
  return that;
}

void PairBinary_delete(struct PairBinary *that)
{
  int i;
  for (i = 0; i < that->n; i++)
    free(that->t[i]);
  free(that->t);
  free(that->nn);
  free(that->space);
  free(that);
}

void PairBinary_add(struct PairBinary *that, int i, int j)
{
  int nmax;
  MAKE_CANONICAL(i,j);
  nmax = that->nn[i];
  if (that->space[i] == nmax)
    that->t[i] = (int *) realloc(that->t[i], (that->space[i] = nmax+SLACK) * sizeof(int));
  that->t[i][nmax] = j;
  that->nn[i] = nmax+1;
}

/* Add without checking that i,j order is canonical */
void PairBinary_quick_add(struct PairBinary *that, int i, int j)
{
  int nmax;
  nmax = that->nn[i];
  if (that->space[i] == nmax)
    that->t[i] = (int *) realloc(that->t[i], (that->space[i] = nmax+SLACK) * sizeof(int));
  that->t[i][nmax] = j;
  that->nn[i] = nmax+1;
}

static int compar(const void *a, const void *b)
{
  if (*((int *) a) > *((int *) b))
    return 1;
  else if (*((int *) a) < *((int *) b))
    return -1;
  else
    return 0;
}

void PairBinary_init(struct PairBinary *that)
{
  int i;
  for (i = 0; i < that->n; i++)
    qsort(that->t[i], that->nn[i], sizeof(int), compar);
  /* Remove duplicates */
  for (i = 0; i < that->n; i++) {
    const int n = that->nn[i];
    if (n > 1) {
      int *p = that->t[i];
      int j, k = 0;
      for (j = 1; j < n; j++)
	if (p[j] != p[k])
	  p[++k] = p[j];
      that->nn[i] = k+1;
    }
  }
}

void PairBinary_clear(struct PairBinary *that)
{
  int i;
  for (i = 0; i < that->n; i++)
    that->nn[i] = 0;
}

int * const *PairBinary_elements(const struct PairBinary *that)
{
  return that->t;
}

const int *PairBinary_number_of_elements(const struct PairBinary *that)
{
  return that->nn;
}

const int *PairBinary_start(const struct PairBinary *that, int i)
{
  return that->t[i];
}

const int *PairBinary_end(const struct PairBinary *that, int i)
{
  return &that->t[i][that->nn[i]];
}

int PairBinary_exists(const struct PairBinary *that, int i, int j)
{
  const int *t;
  int low, high;
  MAKE_CANONICAL(i,j);
  t = that->t[i];
  low = 0;
  high = that->nn[i]-1;
  while (low <= high) {
    const int mid = (low+high)/2;
    if (j < t[mid])
      high = mid - 1;
    else if (j > t[mid])
      low = mid + 1;
    else
      return 1;
  }
  return 0;
}

void PairBinary_show(const struct PairBinary *that)
{
  int i;
  for (i = 0; i < that->n; i++) {
    int j;
    printf("%5d: ", i);
    for (j = 0; j < that->nn[i]; j++)
      printf("%d ", that->t[i][j]);
    printf("\n");
  }
  fflush(stdout);
}

#ifdef SIMPLE_EXAMPLE
/*
  should print out
   0: 0 1 3 
   1: 1 2 4 8 
   2: 0 
   3: 1 8 
   4: 
   5: 1 
   6: 
   7: 1 8 
   8: 0 2 6 9 
   9: 
1 1 0 1 0 0 0 0 0 0 
0 1 1 0 1 0 0 0 1 0 
1 0 0 0 0 0 0 0 0 0 
0 1 0 0 0 0 0 0 1 0 
0 0 0 0 0 0 0 0 0 0 
0 1 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 
0 1 0 0 0 0 0 0 1 0 
1 0 1 0 0 0 1 0 0 1 
0 0 0 0 0 0 0 0 0 0 
*/
int main()
{
  int i, j;
  struct PairBinary *p = PairBinary_new(10);
  PairBinary_add(p,0,1);
  PairBinary_add(p,0,3);
  PairBinary_add(p,0,2);
  PairBinary_add(p,0,2);
  PairBinary_add(p,0,2);
  PairBinary_add(p,0,8);
  PairBinary_add(p,0,0);
  PairBinary_add(p,1,3);
  PairBinary_add(p,1,1);
  PairBinary_add(p,1,3);
  PairBinary_add(p,1,2);
  PairBinary_add(p,1,4);  
  PairBinary_add(p,1,7);  
  PairBinary_add(p,1,5);
  PairBinary_add(p,3,8);
  PairBinary_add(p,3,8);
  PairBinary_add(p,3,8);
  PairBinary_add(p,8,2);
  PairBinary_add(p,8,6);
  PairBinary_add(p,8,6);
  PairBinary_add(p,8,0);
  PairBinary_add(p,8,1);
  PairBinary_add(p,8,7);
  PairBinary_add(p,8,9);
  PairBinary_add(p,8,6);
  PairBinary_init(p);
  PairBinary_show(p);
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++)
      printf("%d ", PairBinary_exists(p,i,j));
    printf("\n");
  }
  PairBinary_delete(p);
}
#endif
