#ifndef RMSDCB_H
#define RMSDCB_H

#include "msys.hpp"
#include "callback.hpp"
#include "cvec.hpp"

class RMSDCallBack : public Callback // io
{
public:
  void update(const MSys &msys);
  void write() const;
  MSys molsys; // io
private:
  CVec c1, c2;
  double rmsd;
  classIO(RMSDCallBack);
};

#endif /* RMSDCB_H */
