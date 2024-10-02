#ifndef PREPRO_H
#define PREPRO_H

#include "str.hpp"
#include "htab.hpp"

class PreprocessedStream
{
public:
  PreprocessedStream();
  ~PreprocessedStream();
  operator istream &() { return *s; }
  void init(int argc, char **argv);
  void show_self();
  void cleanup();
  void write_current_file_and_line(ostream &) const;
private:
  void read_fname(const char *fname);
  void read_file(FILE *, const char *fname, bool = false);
  void read_line(char *ln, const char *fname, int nline);
  Str tmpfname, index_fname;
  ostream *ftmp, *findex;
  Str last_fname;
  int last_nline;
  istream *s;
  HTab<Str> var;
#ifdef USE_MPI
  int size, rank; // for MPI
#endif
};

#endif /* PREPRO_H */
