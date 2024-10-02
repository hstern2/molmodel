#include <stdlib.h>
#include "array.h"
#include "util.h"

ARRAY_DEFINITION(int)
ARRAY_DEFINITION(double)
ARRAY_DEFINITION(cart_t)
ARRAY_DEFINITION(icart_t)
ARRAY_DEFINITION(complex_t)
ARRAY_DEFINITION(const_ptr_t)

#ifdef SIMPLE_EXAMPLE

int main()
{
  int i, j, k, n, **a, ***b;
  double2_t ***c;
  /* 2D array */
  a = int_array2D_new(4, 2);
  n = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++)
      a[i][j] = n++;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 2; j++)
      printf("%4d", a[i][j]);
    printf("\n");
  }
  printf("\n\n");
  int_array2D_delete(a);
  /* 3D array of int */
  b = int_array3D_new(4, 2, 8);
  n = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < 8; k++)
	b[i][j][k] = n++;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < 8; k++)
	printf("%4d", b[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }
  int_array3D_delete(b);
  /* 3D array of double2_t */
  c = double2_t_array3D_new(4, 2, 8);
  n = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < 8; k++) {
	c[i][j][k][0] = n * 0.1;
	c[i][j][k][1] = n * 0.2;
	n++;
      }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < 8; k++)
	printf("%.1f %.1f  ", c[i][j][k][0], c[i][j][k][1]);
      printf("\n");
    }
    printf("\n");
  }
  double2_t_array3D_delete(c);
}

#endif
