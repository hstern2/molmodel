#ifndef PROG_H
#define PROG_H

#include "out.hpp"
#include "str.hpp"

class Program
{
public:

  virtual ~Program() { }
  void system(istream &); // in
  void exit(); // in
  virtual void write_header(); // in
  void timing(); // in
  void output_precision(int n); // in
  void output_file(Str f); // in
  void add_to_search_path(Str f); // in

  classIO(Program);

};

#endif /* PROG_H */
