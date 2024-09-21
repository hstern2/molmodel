#ifndef PBIN_H
#define PBIN_H

#ifdef __cplusplus
extern "C" {
#endif

  struct PairBinary;
  
  struct PairBinary *PairBinary_new(int n);
  void PairBinary_delete(struct PairBinary *);
  void PairBinary_add(struct PairBinary *, int i, int j);
  void PairBinary_quick_add(struct PairBinary *, int i, int j);
  void PairBinary_init(struct PairBinary *);
  void PairBinary_clear(struct PairBinary *);
  const int *PairBinary_number_of_elements(const struct PairBinary *);
  int * const *PairBinary_elements(const struct PairBinary *);
  const int *PairBinary_start(const struct PairBinary *, int i);
  const int *PairBinary_end(const struct PairBinary *, int i);
  int PairBinary_exists(const struct PairBinary *, int i, int j);
  void PairBinary_show(const struct PairBinary *);
  
#ifdef __cplusplus
}
#endif

#endif /* PBIN_H */
