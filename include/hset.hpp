#ifndef HSET_H
#define HSET_H

#include <cassert>
#include <cstring>
#include "out.hpp"
#include "hash.h"

template<class T> struct HSetItem
{
  T t;
  HSetItem *next;
  HSetItem(const T &t_, HSetItem *next_) : t(t_), next(next_) { }
  ~HSetItem() { if (next) delete next; }
};

template<class T> class HSetItr;

template<class T> class HSet {
  friend class HSet<HSet<T> >;
  friend class HSetItr<T>;
  HSetItem<T> **items;
  int n;
public:
  friend unsigned hash_val(const HSet<T> &key)
  {
    unsigned k = 0;
    for (int i = 0; i < key.n; i++)
      for (HSetItem<T> *j = key.items[i]; j; j = j->next)
	k += hash_val(j->t);
    return k;
  }
  friend bool operator==(const HSet<T> &a, const HSet<T> &b)
  {
    if (a.size() != b.size())
      return false;
    for (int i = 0; i < a.n; i++)
      for (HSetItem<T> *j = a.items[i]; j; j = j->next)
	if (!b.exists(j->t))
	  return false;
    return true;
  }
  friend bool operator!=(const HSet<T> &a, const HSet<T> &b)
  {
    return !(a == b);
  }
  HSet(int sz = 101) : n(sz) 
  {
    items = new HSetItem<T> *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  HSet(const HSet<T> &s) : n(s.n)
  {
    items = new HSetItem<T> *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
    add(s);
  }
  HSet & operator=(const HSet<T> &s)
  {
    remove_all();
    if (n != s.n) {
      n = s.n;
      delete[] items;
      items = new HSetItem<T> *[n];
    }
    add(s);
    return *this;
  }
  void add(const T &t) 
  { 
    const int k = hash_val(t) % n;
    const HSetItem<T> *i;
    for (i = items[k]; i; i = i->next)
      if (t == i->t)
	return;
    items[k] = new HSetItem<T>(t, items[k]);
  }
  void add(const HSet<T> &s)
  {
    for (int i = 0; i < s.n; i++)
      for (const HSetItem<T> *j = s.items[i]; j; j = j->next)
	add(j->t);
  }
  void add(int nt, const T *t)
  {
    for (int i = 0; i < nt; i++)
      add(t[i]);
  }
  void remove(const T &t)
  {
    const int k = hash_val(t) % n;
    HSetItem<T> *i = items[k];
    if (!i)
      return;
    if (t == i->t) {
      items[k] = i->next;
      i->next = 0;
      delete i;
      return;
    }
    while (i->next) {
      if (t == i->next->t) {
	HSetItem<T> *j = i->next;
	i->next = j->next;
	j->next = 0;
	delete j;
	return;
      }
      i = i->next;
    }
  }
  bool exists(const T &t) const
  {
    for (const HSetItem<T> *i = items[hash_val(t) % n]; i; i = i->next)
      if (t == i->t)
	return true;
    return false;
  }
  T *get(const T &t) const
  {
    for (HSetItem<T> *i = items[hash_val(t) % n]; i; i = i->next)
      if (t == i->t)
	return &i->t;
    return 0;
  }
  void remove_all()
  {
    for (int i = 0; i < n; i++)
      if (items[i]) {
	delete items[i];
	items[i] = 0;
      }
  }
  ~HSet()
  {
    remove_all();
    delete[] items;
  }
  int size() const 
  {
    int e = 0;
    for (int i = 0; i < n; i++)
      for (const HSetItem<T> *j = items[i]; j; j = j->next)
	e++;
    return e;
  }
};

template<class T> class HSetItr
{
  const HSet<T> *t;
  int k;
  const HSetItem<T> *i;
public:
  HSetItr() : t(0) { }
  void init(const HSet<T> &tab)
  {
    t = &tab;
    k = 0;
    i = t->items[0];
    if (!i)
      next();
  }
  int ok() const { return k < t->n; }
  const T & operator()() { return i->t; }
  void next()
  { 
    if (i == 0 || i->next == 0)
      for (k++ ; k < t->n && (i = t->items[k]) == 0; k++)
	;
    else
      i = i->next;
  }
};

template<class T> inline istream & operator>>(istream &s, HSet<T> &c)
{
  c.remove_all();
  char bracket;
  s >> bracket;
  if (bracket != '{') { // Assume the set contains only one element
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

template<class T> inline ostream & operator<<(ostream &s, const HSet<T> &c)
{
  HSetItr<T> i;
  s << "{ ";
  for (i.init(c); i.ok(); i.next())
    s << i() << " ";
  s << "}";
  return s;
}

typedef HSet<int> HSetInt;

#endif /* HSET_H */
