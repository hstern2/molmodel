#ifndef LST_H
#define LST_H

#include "out.hpp"
#include "str.hpp"
#include "boo.hpp"

template<class T> class LstItr;
template<class T> class LstItrMod;
template<class T> class LstPairItr;

template<class T> struct LstRep
{
  T t;
  LstRep *next;
  LstRep(const T &s, LstRep *n = 0) : t(s), next(n) { }
};

template<class T> class Lst
{
  friend class LstItr<T>;
  friend class LstItrMod<T>;
  friend class LstPairItr<T>;
protected:
  LstRep<T> *head, *tail;
  int n;
public:
  Lst() : head(0), tail(0), n(0) { }
  ~Lst() { remove_all(); }
  Lst<T> & operator=(const Lst<T> &lst)
  {
    remove_all();
    add(lst);
    return *this;
  }
  void remove_all()
  { 
    while (head) {
      tail = head->next;
      delete head;
      head = tail;
    }
    n = 0;
  }
  void remove_last()
  {
    for (LstRep<T> *p = head; p && p->next; p = p->next)
      if (p->next == tail) {
	delete p->next;
	p->next = 0;
	tail = p;
	n--;
	return;
      }
  }
  void add(const T &t)
  { 
    n++;
    if (!head)
      tail = head = new LstRep<T>(t);
    else
      tail = tail->next = new LstRep<T>(t);
  }
  void add(int n, const T *t)
  {
    for (int i = 0; i < n; i++)
      add(t[i]);
  }
  void push(const T &t)
  {
    head = new LstRep<T>(t, head);
    n++;
  }
  T pop()
  {
    T t(head->t);
    LstRep<T> *nn = head->next;
    delete head;
    head = nn;
    n--;
    return t;
  }
  T & first() { 
    return head->t; 
  }
  const T & first() const { 
    return head->t; 
  }
  T & last() {
    return tail->t;
  }
  const T & last() const {
    return tail->t;
  }
  int size() const { return n; }
  inline void add(const Lst<T> &l);
  inline Lst(const Lst<T> &);
};

template<class T> class LstItr
{
private:
  LstRep<T> *i;
public:
  LstItr() : i(0) { }
  void init(const Lst<T> &l) 
    { i = l.head; }
  int ok() const { return i != 0; }
  const T & operator()() { return i->t; }
  void next() { i = i->next; }
};

template<class T> class LstItrMod
{
private:
  LstRep<T> *i;
public:
  LstItrMod() : i(0) { }
  void init(Lst<T> &l) 
  { i = l.head; }
  int ok() const { return i != 0; }
  T & operator()() { return i->t; }
  void next() { i = i->next; }
};

template<class T> class LstPairItr
{
private:
  LstRep<T> *p, *q;
public:
  LstPairItr() : p(0), q(0) { }
  void init(const Lst<T> &l)
    { if ((p = l.head) != 0) q = p->next; }
  int ok() const { return p != 0 && q != 0; }
  T & i() { return p->t; }
  T & j() { return q->t; }
  void next() 
    { 
      if ((q = q->next) == 0 && (p = p->next) != 0) 
	q = p->next; 
    }
};

template<class T> inline void Lst<T>::add(const Lst<T> &l)
{
  LstItr<T> i;
  for (i.init(l); i.ok(); i.next())
    add(i());
}

template<class T> inline Lst<T>::Lst(const Lst<T> &l) : head(0), tail(0), n(0) 
{
  add(l);
}

template<class T> inline istream & operator>>(istream &s, Lst<T> &c)
{
  c.remove_all();
  char bracket;
  s >> bracket;
  if (bracket != '{') { // Assume the list contains only one element
    s.putback(bracket);
    T t;
    s >> t;
    c.add(t);
    return s;
  }
  while(s) {
    s >> bracket;
    if (bracket == '}')
      return s;
    s.putback(bracket);
    T t;
    s >> t;
    c.add(t);
  }
  return s;
}

template<class T> inline ostream & operator<<(ostream &s, const Lst<T> &c)
{
  LstItr<T> i;
  s << "{\n";
  IndentPush();
  for (i.init(c); i.ok(); i.next())
    s << Indent() << i() << "\n";
  IndentPop();
  return s << Indent() << "}\n";
}

typedef Lst<int> LstInt;

inline ostream & operator<<(ostream &s, const LstInt &c)
{
  LstItr<int> i;
  s << "{ ";
  for (i.init(c); i.ok(); i.next())
    s << i() << " ";
  return s << "}";
}

typedef Lst<Str> LstStr;

inline ostream & operator<<(ostream &s, const LstStr &c)
{
  LstItr<Str> i;
  s << "{ ";
  for (i.init(c); i.ok(); i.next())
    s << i() << " ";
  return s << "}";
}

typedef Lst<Boolean> LstBoolean;

inline ostream & operator<<(ostream &s, const LstBoolean &c)
{
  LstItr<Boolean> i;
  s << "{ ";
  for (i.init(c); i.ok(); i.next())
    s << i() << " ";
  return s << "}";
}

#endif /* LST_H */
