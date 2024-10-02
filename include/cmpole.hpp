#ifndef CMPOLE_H
#define CMPOLE_H

#include "out.hpp"

struct Quadrupole
{
  double XX, YY, ZZ, XY, XZ, YZ; // io
  void zero() { XX = YY = ZZ = XY = XZ = YZ = 0; }
  void apply(double (*f)(double)) 
  {
    XX = f(XX);
    YY = f(YY);
    ZZ = f(ZZ);
    XY = f(XY);
    XZ = f(XZ);
    YZ = f(YZ);
  }      
  Quadrupole map(double (*f)(double)) const
  {
    Quadrupole that(*this);
    that.apply(f);
    return that;
  }
  void make_traceless()
  {
    const double t = (XX+YY+ZZ)/3;
    XX -= t;
    YY -= t;
    ZZ -= t;
  }
  classIO(Quadrupole);
};

struct Octopole
{
  double XXX, YYY, ZZZ, XYY; // io
  double XXY, XXZ, XZZ, YZZ; // io
  double YYZ, XYZ; // io
  void zero() { XXX = YYY = ZZZ = XYY = XXY = XXZ = XZZ = YZZ = YYZ = XYZ = 0; }
  void apply(double (*f)(double)) 
  {
    XXX = f(XXX);
    YYY = f(YYY);
    ZZZ = f(ZZZ);
    XYY = f(XYY);
    XXY = f(XXY);
    XXZ = f(XXZ);
    XZZ = f(XZZ);
    YZZ = f(YZZ);
    YYZ = f(YYZ);
    XYZ = f(XYZ);
  }
  void make_traceless()
  {
    double t;
    t = (XXX+XYY+XZZ)/3;
    XXX -= t;
    XYY -= t;
    XZZ -= t;
    t = (XXY+YYY+YZZ)/3;
    XXY -= t;
    YYY -= t;
    YZZ -= t;
    t = (XXZ+YYZ+ZZZ)/3;
    XXZ -= t;
    YYZ -= t;
    ZZZ -= t;
  }
  Octopole map(double (*f)(double)) const
  {
    Octopole that(*this);
    that.apply(f);
    return that;
  }
  classIO(Octopole);
};

struct Hexadecapole
{
  double XXXX, YYYY, ZZZZ, XXXY; // io
  double XXXZ, YYYX, YYYZ, ZZZX; // io
  double ZZZY, XXYY, XXZZ, YYZZ; // io
  double XXYZ, YYXZ, ZZXY; // io
  void zero() 
  {
    XXXX = YYYY = ZZZZ = XXXY = 0;
    XXXZ = YYYX = YYYZ = ZZZX = 0;
    ZZZY = XXYY = XXZZ = YYZZ = 0;
    XXYZ = YYXZ = ZZXY = 0;
  }
  void apply(double (*f)(double)) 
  {
    XXXX = f(XXXX);
    YYYY = f(YYYY);
    ZZZZ = f(ZZZZ);
    XXXY = f(XXXY);
    XXXZ = f(XXXZ);
    YYYX = f(YYYX);
    YYYZ = f(YYYZ);
    ZZZX = f(ZZZX);
    ZZZY = f(ZZZY);
    XXYY = f(XXYY);
    XXZZ = f(XXZZ);
    YYZZ = f(YYZZ);
    XXYZ = f(XXYZ);
    YYXZ = f(YYXZ);
    ZZXY = f(ZZXY);
  }
  void make_traceless()
  {
    double t;
    t = (XXXX+XXYY+XXZZ)/3;
    XXXX -= t;
    XXYY -= t;
    XXZZ -= t;
    t = (XXXY+YYYX+ZZXY)/3;
    XXXY -= t;
    YYYX -= t;
    ZZXY -= t;
    t = (XXXZ+YYXZ+ZZZX)/3;
    XXXZ -= t;
    YYXZ -= t;
    ZZZX -= t;
    t = (XXYY+YYYY+YYZZ)/3;
    XXYY -= t;
    YYYY -= t;
    YYZZ -= t;
    t = (XXYZ+YYYZ+ZZZY)/3;
    XXYZ -= t;
    YYYZ -= t;
    ZZZY -= t;
    t = (XXZZ+YYZZ+ZZZZ)/3;
    XXZZ -= t;
    YYZZ -= t;
    ZZZZ -= t;
  }
  Hexadecapole map(double (*f)(double)) const
  {
    Hexadecapole that(*this);
    that.apply(f);
    return that;
  }
  classIO(Hexadecapole);
};

#endif /* CMPOLE_H */
