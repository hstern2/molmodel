#ifndef PDB_H
#define PDB_H

#include "out.hpp"
#include "str.hpp"

class Cartesian;

struct PDBAtomRecord
{
  PDBAtomRecord();
  void read(const char *buf, Cartesian &r);
  void write(ostream &, const char *sym, int i, const Cartesian &r) const;
  char record[7], name[5], altLoc[2], resName[4], chainID[2], atomID[32], iCode[2];
  int resSeq;
  int is_heterogen() const;
  int is_alternate_location() const;
  int is_hydrogen() const;
  int is_empty() const;
  Str symbol() const;
  
  friend istream & operator>>(istream &, PDBAtomRecord &);
  friend ostream & operator<<(ostream &, const PDBAtomRecord &);
};

#endif /* PDB_H */
