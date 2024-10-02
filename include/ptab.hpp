#ifndef PTAB_H
#define PTAB_H

#include <cassert>
#include "out.hpp"

/***
 *  Table of type T indexed by pairs of ints (i,j)
 *  (order of i and j doesn't matter)
 ***/

template<class T> class PTabItrMod;
template<class T> class PTabItr;

template<class T> struct PTabItem 
{
  const int i, j;
  T val;
  PTabItem *next;
  PTabItem(int k1, int k2, const T &ob, PTabItem *n) : 
    i(k1), j(k2), val(ob), next(n) { }
  ~PTabItem() { if (next) delete next; next = 0; }
};

template<class T> class PTab {
  friend class PTabItrMod<T>;
  friend class PTabItr<T>;
  PTabItem<T> **items;
  int n;
  unsigned hash(int k1, int k2) const { return (unsigned) (k1*104729 + k2); }
public:
  T default_value;
  inline const T& operator()(int, int) const;
  inline T& operator()(int, int);
  PTab(int sz = 101) : n(sz) 
  {
    items = new PTabItem<T> *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  PTab(int sz, const T &def) : n(sz), default_value(def)
  {
    items = new PTabItem<T> *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  PTab(const PTab<T> &t) : n(t.n), default_value(t.default_value)
  {
    int i;
    items = new PTabItem<T> *[n];
    for (i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (PTabItem<T> *j = t.items[i]; j; j = j->next)
	(*this)(j->i,j->j) = j->val;
  }
  PTab<T> & operator=(const PTab<T> &t)
  {
    int i;
    remove_all();
    delete[] items;
    n = t.n;
    items = new PTabItem<T> *[n];
    for (i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (PTabItem<T> *j = t.items[i]; j; j = j->next)
	(*this)(j->i,j->j) = j->val;
    return *this;
  }
  ~PTab();
  int size() const;
  bool exists(int, int) const;
  T *get(int, int) const;
  void remove_all();
  void remove(int k1, int k2)
  {
    if (k2 < k1) { const int tmp = k1; k1 = k2; k2 = tmp; }
    const int k = hash(k1,k2) % n;
    PTabItem<T> *i = items[k];
    if (!i)
      return;
    if (k1 == i->i && k2 == i->j) {
      items[k] = i->next;
      i->next = 0;
      delete i;
      return;
    }
    while (i->next) {
      if (k1 == i->next->i && k2 == i->next->j) {
	PTabItem<T> *j = i->next;
	i->next = j->next;
	j->next = 0;
	delete j;
	return;
      }
      i = i->next;
    }
  }
};

template<class T> inline PTab<T>::~PTab()
{
  remove_all();
  delete[] items;
  items = 0;
}

template<class T> inline T& PTab<T>::operator()(int k1, int k2)
{
  if (k2 < k1) { const int tmp = k1; k1 = k2; k2 = tmp; }
  const int k = hash(k1,k2) % n;
  PTabItem<T> *i;
  for (i = items[k]; i && !(k1 == i->i && k2 == i->j); i = i->next)
    ;
  if (i == 0) { /* not found; create a new item */
    i = new PTabItem<T>(k1, k2, default_value, items[k]);
    items[k] = i;
  }
  return i->val;
}

template<class T> inline const T & PTab<T>::operator()(int k1, int k2) const
{
  if (k2 < k1) { const int tmp = k1; k1 = k2; k2 = tmp; }
  const int k = hash(k1,k2) % n;
  for (PTabItem<T> *i = items[k]; i; i = i->next)
    if (k1 == i->i && k2 == i->j)
      return i->val;
  return default_value;
}

template<class T> inline bool PTab<T>::exists(int k1, int k2) const
{
  if (k2 < k1) { const int tmp = k1; k1 = k2; k2 = tmp; }
  const int k = hash(k1,k2) % n;
  for (PTabItem<T> *i = items[k]; i; i = i->next)
    if (k1 == i->i && k2 == i->j)
      return true;
  return false;
}

template<class T> inline T *PTab<T>::get(int k1, int k2) const
{
  if (k2 < k1) { const int tmp = k1; k1 = k2; k2 = tmp; }
  const int k = hash(k1,k2) % n;
  for (PTabItem<T> *i = items[k]; i; i = i->next)
    if (k1 == i->i && k2 == i->j)
      return &i->val;
  return 0;
}

template<class T> inline int PTab<T>::size() const
{
  int e = 0;
  for (int i = 0; i < n; i++)
    for (PTabItem<T> *j = items[i]; j; j = j->next)
      e++;
  return e;
}

template<class T> inline void PTab<T>::remove_all()
{ 
  for (int i = 0; i < n; i++) {
    delete items[i];
    items[i] = 0;
  }
}

template<class T> class PTabItrMod
{
  PTab<T> *t;
  int k;
  PTabItem<T> *it;
public:
  PTabItrMod() : t(0) { }
  void init(PTab<T> &tab)
  {
    t = &tab;
    k = 0;
    it = t->items[0];
    if (!it)
      next();
  }
  bool ok() const
  { return k < t->n; }
  int i() const
  { return it->i; }
  int j() const
  { return it->j; }
  T & val()
  { return it->val; }
  void next()
  { 
    if (it == 0 || it->next == 0)
      for (k++ ; k < t->n && (it = t->items[k]) == 0; k++)
	;
    else
      it = it->next;
  }
};

template<class T> class PTabItr
{
  const PTab<T> *t;
  int k;
  const PTabItem<T> *it;
public:
  PTabItr() : t(0) { }
  void init(const PTab<T> &tab)
  {
    t = &tab;
    k = 0;
    it = t->items[0];
    if (!it)
      next();
  }
  bool ok() const
  { return k < t->n; }
  int i() const
  { return it->i; }
  int j() const
  { return it->j; }
  const T & val()
  { return it->val; }
  void next()
  { 
    if (it == 0 || it->next == 0)
      for (k++ ; k < t->n && (it = t->items[k]) == 0; k++)
	;
    else
      it = it->next;
  }
};

template<class T> inline istream & operator>>(istream &s, PTab<T> &t)
{
  char bracket;
  int k1, k2;
  s >> bracket;
  if (bracket != '{')
    die("PTab: operator>>: expecting '{'");
  while (s) {
    s >> bracket;
    if (bracket == '}')
      return s;
    else
      s.putback(bracket);
    s >> k1 >> k2 >> t(k1,k2);
  }
  die("PTab::operator>>: EOF reached before reading '}'");
  return s;
}

template<class T> inline ostream & operator<<(ostream &s, const PTab<T> &t)
{
  PTabItr<T> ti;
  s << Indent() << "{\n";
  IndentPush();
  for (ti.init(t); ti.ok(); ti.next())
    s << Indent() << ti.i() << " " << ti.j() << " " << ti.val() << "\n";
  IndentPop();
  return s << Indent() << "}\n";
}

#endif /* PTAB_H */
