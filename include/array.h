#ifndef ARRAY_H
#define ARRAY_H

#include "util.h"
#include "cart.h"

typedef const void *const_ptr_t;
typedef void *ptr_t;

/* Poor man's template :) */
#define ARRAY_DECLARATION(T) \
  T *T##_array_new(size_t n);\
  void T##_array_delete(T *);\
  T **T##_array2D_new(size_t n1, size_t n2);\
  void T##_array2D_delete(T **);\
  T ***T##_array3D_new(size_t n1, size_t n2, size_t n3);\
  void T##_array3D_delete(T ***)

#define ARRAY_DEFINITION(T) \
\
T * T##_array_new(size_t n) \
{ \
  return (T *) safe_malloc(n*sizeof(T));\
} \
 \
void T##_array_delete(T *p) \
{ \
  if (p) { \
    free(p); \
  } \
} \
 \
T ** T##_array2D_new(size_t n1, size_t n2) \
{ \
  size_t i; \
  T **p; \
  if (n1 <= 0 || n2 <= 0) \
    return 0; \
  p = (T **) safe_malloc(n1*sizeof(T *)); \
  p[0] = (T *) safe_malloc(n1*n2*sizeof(T)); \
  for (i = 1; i < n1; i++) \
    p[i] = p[0] + i*n2; \
  return p; \
} \
 \
void T##_array2D_delete(T **p) \
{ \
  if (p) { \
    free(p[0]); \
    free(p); \
  } \
} \
 \
T ***T##_array3D_new(size_t n1, size_t n2, size_t n3) \
{ \
  size_t i; \
  T ***p = 0; \
  if (n1 <= 0 || n2 <= 0 || n3 <= 0) \
    return 0; \
  p = (T ***) safe_malloc(n1*sizeof(T **)); \
  p[0] = (T **) safe_malloc(n1*n2*sizeof(T *)); \
  p[0][0] = (T *) safe_malloc(n1*n2*n3*sizeof(T)); \
  for (i = 0; i < n1; i++) { \
    size_t j; \
    p[i] = p[0] + i*n2; \
    for (j = 0; j < n2; j++) \
      p[i][j] = p[0][0] + (i*n2+j)*n3; \
  } \
  return p; \
} \
 \
void T##_array3D_delete(T ***p) \
{ \
  if (p) { \
    free(p[0][0]); \
    free(p[0]); \
    free(p); \
  } \
}

#ifdef __cplusplus
extern "C" {
#endif

  ARRAY_DECLARATION(int);
  ARRAY_DECLARATION(double);
  ARRAY_DECLARATION(cart_t);
  ARRAY_DECLARATION(icart_t);
  ARRAY_DECLARATION(complex_t);
  ARRAY_DECLARATION(const_ptr_t);

#ifdef __cplusplus
}
#endif

#endif /* ARRAY_H */
