#ifndef CALLBACK_H
#define CALLBACK_H

#include "out.hpp"
#include "lst.hpp"

struct MSys;
class Cartesian;

class Callback
{
  friend class MDyn;
public:
  double wait; // io
  int update_interval; // io
  Callback();
  virtual ~Callback();
  virtual void init(const MSys &, int write_interval);
  void maybe_update(double time_elapsed, int istep, const MSys &);
  void maybe_write(double time_elapsed) const;
  classIO(Callback);
protected:
  const Cartesian *f;
  virtual void write() const = 0;
  virtual void update(const MSys &) = 0;
};

struct CallbackPtr
{
  Callback *p;
  operator Callback *() { return p; }
  operator const Callback *() const { return p; }
  Callback *operator->() { return p; }
  const Callback *operator->() const { return p; }
  CallbackPtr(Callback *p_ = 0) : p(p_) { }
};

typedef Lst<CallbackPtr> CallbackList;

istream & operator>>(istream &, CallbackPtr &);

#endif /* CALLBACK_H */
