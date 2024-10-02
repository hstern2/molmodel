#ifndef HTAB_H
#define HTAB_H

#include <cassert>
#include <cstring>
#include "str.hpp"
#include "hash.h"

/***
 *
 *  Hash table of type T indexed by type K
 *
 ***/

template<class T, class K = Str> class HTabIterator;
template<class T, class K = Str> class ConstHTabIterator;

template<class T, class K = Str> struct HTabItem 
{
  K key;
  T val;
  HTabItem *next;
  HTabItem(const K &k, const T &ob, HTabItem *n) : key(k), val(ob), next(n) { }
  ~HTabItem() { if (next) delete next; next = 0; }
};

template<class T, class K = Str> class HTab {
  friend class HTabIterator<T,K>;
  friend class ConstHTabIterator<T,K>;
  HTabItem<T,K> **items;
  int n;
public:
  inline T& operator[](const K &key)
  {
    int k = hash_val(key) % n;
    HTabItem<T,K> *i;
    for (i = items[k]; i && !(key == i->key); i = i->next)
      ;
    if (i == 0) { /* not found; create a new item */
      i = new HTabItem<T,K>(key, T(), items[k]);
      items[k] = i;
    }
    return i->val;
  }
  HTab(int sz = 101) : n(sz) 
  {
    items = new HTabItem<T,K> *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  HTab(const HTab<T,K> &t) : n(t.n)
  {
    int i;
    items = new HTabItem<T,K> *[n];
    for (i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (HTabItem<T,K> *j = t.items[i]; j; j = j->next)
	(*this)[j->key] = j->val;
  }
  HTab<T,K> & operator=(const HTab<T,K> &t)
  {
    int i;
    remove_all();
    delete[] items;
    n = t.n;
    items = new HTabItem<T,K> *[n];
    for (i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (HTabItem<T,K> *j = t.items[i]; j; j = j->next)
	(*this)[j->key] = j->val;
    return *this;
  }
  void remove_all()
  {
    for (int i = 0; i < n; i++) {
      delete items[i];
      items[i] = 0;
    }
  }
  ~HTab()
  {
    remove_all();
    delete[] items;
  }
  int size() const
  {
    int e = 0;
    for (int i = 0; i < n; i++)
      for (HTabItem<T,K> *j = items[i]; j; j = j->next)
	e++;
    return e;
  }
  bool exists(const K &key) const {
    const int k = hash_val(key) % n;
    for (HTabItem<T,K> *i = items[k]; i; i = i->next)
      if (key == i->key)
	return true;
    return false;
  }
  T *get(const K &key) const 
  {
    const int k = hash_val(key) % n;
    for (HTabItem<T,K> *i = items[k]; i; i = i->next)
      if (key == i->key)
	return &i->val;
    return 0;
  }
  void remove(const K &key)
  {
    const int k = hash_val(key) % n;
    HTabItem<T,K> *i = items[k];
    if (!i)
      return;
    if (key == i->key) {
      items[k] = i->next;
      i->next = 0;
      delete i;
      return;
    }
    while (i->next) {
      if (key == i->next->key) {
	HTabItem<T,K> *j = i->next;
	i->next = j->next;
	j->next = 0;
	delete j;
	return;
      }
      i = i->next;
    }
  }
};

template<class T, class K> class HTabIterator
{
  HTab<T,K> *t;
  int k;
  HTabItem<T,K> *i;
public:
  HTabIterator() : t(0) { }
  void init(HTab<T,K> &tab)
  {
    t = &tab;
    k = 0;
    i = t->items[0];
    if (!i)
      next();
  }
  int ok() const
  { return k < t->n; }
  const K &key()
  { return i->key; }
  T & val()
  { return i->val; }
  void next()
  { 
    if (i == 0 || i->next == 0)
      for (k++ ; k < t->n && (i = t->items[k]) == 0; k++)
	;
    else
      i = i->next;
  }
};

template<class T, class K> class ConstHTabIterator
{
  const HTab<T,K> *t;
  int k;
  const HTabItem<T,K> *i;
public:
  ConstHTabIterator() : t(0) { }
  void init(const HTab<T,K> &tab)
  {
    t = &tab;
    k = 0;
    i = t->items[0];
    if (!i)
      next();
  }
  int ok() const { return k < t->n; }
  const K &key() { return i->key; }
  const T &val() { return i->val; }
  void next()
  { 
    if (i == 0 || i->next == 0)
      for (k++ ; k < t->n && (i = t->items[k]) == 0; k++)
	;
    else
      i = i->next;
  }
};

template<class T, class K> inline ostream & operator<<(ostream &s, const HTab<T,K> &t)
{
  ConstHTabIterator<T,K> ti;
  int i;
  int n = t.size();
  const K **keys = new const K*[n];
  for (i = 0, ti.init(t); ti.ok(); i++, ti.next())
    keys[i] = &ti.key();
  /* Bubble sort */
  bool any = true;
  while(any)
    for (any = false, i = 1; i < n; i++)
      if (*keys[i-1] > *keys[i]) {
	const K *tmp = keys[i];
	keys[i] = keys[i-1];
	keys[i-1] = tmp;
	any = true;
      }
  s << "{\n";
  IndentPush();
  for (i = 0; i < n; i++)
    s << Indent() << *keys[i] << " " << *t.get(*keys[i]) << "\n";
  IndentPop();
  delete[] keys;
  keys = 0;
  return s << Indent() << "}\n";
}

template<class T, class K> inline istream & operator>>(istream &s, HTab<T,K> &t)
{
  char bracket;
  s >> bracket;
  if (bracket != '{') {
    s.putback(bracket);
    Str fname;
    s >> fname;
    istream *ftmp = StreamSearch(fname);
    *ftmp >> t;
    delete ftmp;
    return s;
  }
  while (s) {
    s >> bracket;
    if (bracket == '}')
      return s;
    else
      s.putback(bracket);
    K key;
    s >> key;
    s >> t[key];
  }
  die("HTab::operator>>: EOF reached before reading '}'");
  return s;
}

#endif /* HTAB_H */
