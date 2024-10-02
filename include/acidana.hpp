#ifndef ACIDANA_H
#define ACIDANA_H

#include "callback.hpp"
#include "vec.hpp"

class AcidAnalysis  : public Callback
{
public:
  void init(const MSys &, int write_interval);
  void write() const;
protected:
  void update(const MSys &);
private:
  const MSys *msysp;
  Vec<int> np;
  Vec<Vec<int> > istate;
};

#endif /* ACIDANA_H */
