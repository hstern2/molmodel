#include <sstream>
#include "apattern.hpp"
#include "atom.hpp"

class AtomTest
{
public:
  virtual ~AtomTest() { }
  virtual int matches(int a, const HSet<int> &seen, const int *path, const Atom *) const = 0;
  virtual const Str *necessary_property() { return 0; }
};

class PropertyTest : public AtomTest
{
public:
  PropertyTest(const char *p_) : p(p_) { }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  { return !seen.exists(a) && atom[a].property.exists(p); }
  const Str *necessary_property() { return &p; }
private:
  const Str p;
};

class NeighborTest : public AtomTest
{
public:
  NeighborTest(int n_) : n(n_) { }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  { 
    if (seen.exists(a))
      return false;
    LstItr<int> j;
    int nn = 0;
    for (j.init(atom[a].neighbors); j.ok(); j.next())
      if (!atom[j()].is_dummy())
	nn++;
    return nn == n;
  }
private:
  const int n;
};

class PropertyNegation : public AtomTest
{
public:
  PropertyNegation(const AtomTest *t_) : t(t_) { }
  ~PropertyNegation() { delete t; }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  { return !seen.exists(a) && !t->matches(a, seen, path, atom); }
private:
  const AtomTest *t;
};

class PropertyDisjunction : public AtomTest
{
public:
  PropertyDisjunction(Lst<AtomTest *> &t) { terms.copy_from_list(t); }
  ~PropertyDisjunction()
  {
    for (int i = 0; i < terms.size(); i++)
      delete terms[i];
  }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  {
    if (seen.exists(a))
      return 0;
    for (int i = 0; i < terms.size(); i++)
      if (terms[i]->matches(a, seen, path, atom))
	return 1;
    return 0;
  }
private:
  Vec<AtomTest *> terms;
};

class PropertyConjunction : public AtomTest
{
public:
  PropertyConjunction(Lst<AtomTest *> &t) { terms.copy_from_list(t); }
  ~PropertyConjunction()
  {
    for (int i = 0; i < terms.size(); i++)
      delete terms[i];
  }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  {
    if (seen.exists(a))
      return 0;
    for (int i = 0; i < terms.size(); i++)
      if (!terms[i]->matches(a, seen, path, atom))
	return 0;
    return 1;
  }
  const Str *necessary_property() 
  {
    for (int i = 0; i < terms.size(); i++) {
      const Str *p = terms[i]->necessary_property();
      if (p)
	return p;
    }
    return 0;
  }
private:
  Vec<AtomTest *> terms;
};

class RingClosure : public AtomTest
{
public:
  RingClosure(int n_) : n(n_) { }
  int matches(int a, const HSet<int> &seen, const int *path, const Atom *atom) const
  { return path[n] == a; }
private:
  const int n;
};

AtomPattern::~AtomPattern() 
{ 
  for (int i = 0; i < ntest; i++)
    delete test[i];
}

static void split(char *s, char c, Lst<char *> &lst)
{
  lst.add(s);
  while(*s++)
    if (*s == c) {
      *s++ = '\0';
      lst.add(s);
    }
}

static AtomTest *new_property_negation(const char *s)
{
  const int ns = strlen(s);
  if (ns == 0)
    die("new_property_negation: expecting a property name");
  if (ns > 1 && s[0] == '!')
    return new PropertyNegation(new_property_negation(&s[1]));
  if (ns == 2 && s[0] == 'x' && isdigit(s[1]))
    return new NeighborTest(atoi(&s[1]));
  return new PropertyTest(s);
}

static AtomTest *new_property_disjunction(char *tmp) 
{
  Lst<char *> lst;
  split(tmp, ',', lst);
  if (lst.size() == 1)
    return new_property_negation(lst.first());
  LstItr<char *> i;
  Lst<AtomTest *> t;
  for (i.init(lst); i.ok(); i.next())
    t.add(new_property_negation(i()));
  return new PropertyDisjunction(t);
}

static AtomTest *new_property_conjunction(char *tmp) 
{
  Lst<char *> lst;
  split(tmp, '&', lst);
  if (lst.size() == 1)
    return new_property_disjunction(lst.first());
  LstItr<char *> i;
  Lst<AtomTest *> t;
  for (i.init(lst); i.ok(); i.next())
    t.add(new_property_disjunction(i()));
  return new PropertyConjunction(t);
}

static Str space_out(const char *tmp, const char *c)
{
  const int n = strlen(tmp);
  int i,nc = 0;
  for (i = 0; i < n; i++)
    if (strchr(c,tmp[i]))
      nc++;
  if (nc == 0)
    return Str(tmp);
  char *tmp2 = new char[n+2*nc+1];
  int j = 0;
  for (i = 0; i < n; i++)
    if (strchr(c,tmp[i])) {
      tmp2[j++] = ' ';
      tmp2[j++] = tmp[i];
      tmp2[j++] = ' ';
    } else {
      tmp2[j++] = tmp[i];
    }
  assert(j == n+2*nc);
  tmp2[j] = 0;
  Str s(tmp2);
  delete[] tmp2;
  return s;
}

AtomPattern::AtomPattern(const char *pat)
{
  /* Find parentheses in pattern */
  Str pattern = space_out(pat,"()@");
  Lst<int> save;
  istringstream s((const char *) pattern);
  ntest = 0;
  int can_deal_with_left_paren = 0;
  Lst<AtomTest *> tmptest;
  Lst<int> tmpnext;
  char c;
  Str tmp;
  int nr;
  int nring[11];
  char cn[2];
  for (int i = 0; i < 11; i++)
    nring[i] = -1;
  while (s >> c)
    switch (c) {
    case '(':
      if (!can_deal_with_left_paren)
	die("AtomPattern::init: too many '(' in pattern '%s'", pat);
      save.add(tmpnext.last());
      can_deal_with_left_paren = 0;
      break;
    case ')':
      if (save.size() == 0)
	die("AtomPattern::init: too many ')' in pattern '%s'", pat);
      tmpnext.last() = save.last();
      save.remove_last();
      can_deal_with_left_paren = 1;
      break;
    case '@':
      s >> cn[0];
      cn[1] = '\0';
      if (isdigit(cn[0])) {
	nr = atoi(cn);
      } else {
	nr = 10;
	s.putback(cn[0]);
      }
      if (ntest == 0)
	die("AtomPattern::init: ring closure must not appear first in pattern '%s'", pat);
      if (nring[nr] == -1)
	nring[nr] = ntest-1;
      else {
	tmptest.add(new RingClosure(nring[nr]));
	tmpnext.add(ntest-1);
	ntest++;
	nring[nr] = -1;
      }
      break;
    default:
      s.putback(c);
      Str tmp;
      s >> tmp;
      tmptest.add(new_property_conjunction(tmp));
      tmpnext.add(ntest++);
      can_deal_with_left_paren = 1;
    }
  test.copy_from_list(tmptest);
  next.copy_from_list(tmpnext);
}

void AtomPattern::search(int a, int i, HSet<int> &seen, int *path, Atom *atom)
{
  const AtomTest &t = *test[i];
  if (t.matches(a, seen, path, atom)) {
    path[i] = a;
    if (i < ntest-1) {
      seen.add(a);
      LstItr<int> n;
      for (n.init(atom[path[next[i]]].neighbors); n.ok(); n.next())
	search(n(), i+1, seen, path, atom);
      seen.remove(a);
    } else {
      succeed(path, atom);
    }
  }
}

void AtomPattern::run(Vec<Atom> &atom, const HTab<HSet<int> > &atoms_with_property)
{
  HSet<int> seen;
  Vec<int> path(ntest);
  if (test.size() == 0)
    return;
  const Str *s = test[0]->necessary_property();
  if (s) {
    const HSet<int> *awp = atoms_with_property.get(*s);
    if (awp) {
      HSetItr<int> i;      
      for (i.init(*awp); i.ok(); i.next())
	search(i(), 0, seen, path, atom);
    }
  } else {
    for (int i = 0; i < atom.size(); i++)
      search(i, 0, seen, path, atom);
  }
}

void AtomPattern::succeed(const int *path, Atom *atom)
{
  Out() << "Path matched: ";
  for (int i = 0; i < ntest; i++)
    Out() << path[i] << " ";
  Out() << "\n";
}
