/* $Id: out.h,v 1.15 2009/12/06 18:48:27 hstern Exp $ */

/*
 * Copyright (c) 2008 Harry A. Stern
 * Copyright (c) 2008 University of Rochester
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#ifndef OUT_H
#define OUT_H

#include <cstdio>
#include <iostream>
#include <cassert>
using namespace std;

#define classIO(T) \
  friend istream & operator>>(istream &, T &); \
  int read_field(istream &, const char *); \
  void write_field_names() const; \
  friend ostream & operator<<(ostream &, const T &); \
  void write_fields(ostream &) const; \
  void show_self();

#define enumIO(T) \
  istream & operator>>(istream &, T &); \
  ostream & operator<<(ostream &, const T &)

#define SHOW(x) Out() << #x << ": " << x << "\n"

void OutInit();
void OutFinalize();
void OutSetToFile(const char *);
void OutPrintf(const char *, ...);
void SetPrecision(int);
void AddToSearchPath(const char *);

FILE *FileSearch(const char *);
istream *StreamSearch(const char *);

ostream &Out();
ostream *FileStream(const char *);
ostream *AppendFileStream(const char *);

bool end_of_file(FILE *);
class Str;
void TempName(Str &);
void IndentPush();
void IndentPop();
const char *Indent();

class FatalException { };

#define insist(x) assert(x)
void die(const char *s, ...);

#endif /* OUT_H */
