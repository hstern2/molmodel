#ifndef GEOHIST_H
#define GEOHIST_H

#include "callback.hpp"
#include "ntuple.hpp"
#include "hist.hpp"
#include "str.hpp"
#include "lst.hpp"

class BondStatistics
{
public:
  Str2 types; // in
  Str file; // in
  Histogram histogram; // in
  void init(const MSys &);
  void update(const MSys &);
  void write() const;
  classIO(BondStatistics);
private:
  Lst<Int2> lst;
};

class AngleStatistics
{
public:
  Str3 types; // in
  Str file; // in
  Histogram histogram; // in
  void init(const MSys &);
  void update(const MSys &);
  void write() const;
  classIO(AngleStatistics);
private:
  Lst<Int3> lst;
};

class DihedralStatistics
{
public:
  Str4 types; // in
  Str file; // in
  Histogram histogram; // in
  DihedralStatistics() 
  {
    histogram.lo = -180;
    histogram.hi = 180;
    histogram.is_dynamic = false;
  }
  void init(const MSys &);
  void update(const MSys &);
  void write() const;
  classIO(DihedralStatistics);
private:
  Lst<Int4> lst;
};

class GeometryStatistics : public Callback // io
{
public:
  Lst<BondStatistics> bond; // in
  Lst<AngleStatistics> angle; // in
  Lst<DihedralStatistics> dihedral; // in
  void init(const MSys &, int write_interval);
  void update(const MSys &);
  void write() const;
  classIO(GeometryStatistics);
};

#endif /* GEOHIST_H */
