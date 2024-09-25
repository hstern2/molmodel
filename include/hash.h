#ifndef HASH_H
#define HASH_H

static inline unsigned hash_val(int n) { return (unsigned) n; }

static inline unsigned hash_val(void *n) 
{ 
  const char *s = (const char *) &n;
  unsigned i, k = 0;
  for (i = 0; i < sizeof(void *)/sizeof(char); i++)
    k = 31*k + s[i];
  return k;
}

static inline unsigned hash_val(const char *s)
{
  unsigned k;
  for (k = 0; *s; s++)
    k = 31*k + (*s);
  return k;
}

#endif /* HASH_H */
