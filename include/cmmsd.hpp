#ifndef CMMSD_H
#define CMMSD_H

#include "msys.hpp"
#include "callback.hpp"
#include "cvec.hpp"

class CMMeanSquareDistance : public Callback // io
{
public:
  void init(const MSys &, int write_interval);
  void update(const MSys &);
  void write() const;
  MSys molsys; // io
private:
  CVec c0;
  double msd;
  classIO(CMMeanSquareDistance);
};


#endif /* CMMSD_H */
