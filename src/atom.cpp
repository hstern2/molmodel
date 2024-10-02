#include "atom.hpp"
#include "htab.hpp"
#include "units.h"
#include "cmplx.hpp"

static HTab<Element> *element_table = 0;
static HTab<Str> *element_for_tinker_symbol_table = 0;

ostream & operator<<(ostream &s, const Element *e)
{
  if (e)
    return s << *e;
  else
    return s << "{ }\n";
}

void Atom::finalize()
{
  if (element_table)
    delete element_table;
  if (element_for_tinker_symbol_table)
    delete element_for_tinker_symbol_table;
  element_table = 0;
  element_for_tinker_symbol_table = 0;
}

static void init_element_for_tinker_symbol_table()
{
  if (element_for_tinker_symbol_table)
    delete element_for_tinker_symbol_table;
  element_for_tinker_symbol_table = new HTab<Str>;
  istream *s = StreamSearch("_msim_element_for_tinker_symbol");
  *s >> *element_for_tinker_symbol_table;
}

static void init_element_table()
{
  if (element_table)
    delete element_table;
  element_table = new HTab<Element>;
  istream *s = StreamSearch("_msim_elements");
  *s >> *element_table;
  HTabIterator<Element> i;
  for (i.init(*element_table); i.ok(); i.next())
    i.val().mass = g_mol_to_mass_unit(i.val().mass);
  delete s;
}

static const Element *element_for(const char *sym)
{
  if (!element_table)
    init_element_table();
  const Element *e = element_table->get(sym);
  if (!e)
    die("Atom::element_for: unknown element: %s", sym);
  return e;
}

void Atom::set_mass(Str sym, double mass)
{
  if (!element_table)
    init_element_table();
  Element *e = element_table->get(sym);
  if (!e)
    die("Atom::set_mass: unknown element: %s", (const char *) sym);
  e->mass = g_mol_to_mass_unit(mass);
}

static const char *element_for_tinker_symbol(const char *s)
{
  if (!element_for_tinker_symbol_table)
    init_element_for_tinker_symbol_table();
  if (!element_for_tinker_symbol_table->exists(s))
    die("MSys::element_for_tinker_symbol_table: no element for tinker symbol %s", s);
  return (*element_for_tinker_symbol_table)[s];
}

Atom::Atom() : 
  position(0,0,0), 
  velocity(0,0,0), 
  symbol("dummy"),
  type(symbol),
  element(element_for(symbol))
{ }

Atom::Atom(Str sym) : 
  position(0,0,0), 
  velocity(0,0,0),
  symbol(sym),
  type(sym),
  element(element_for(sym))
{ 
  property.add(symbol);
}

Atom::Atom(Str sym, Str t) :
  position(0,0,0),
  velocity(0,0,0),
  symbol(sym),
  type(t),
  element(element_for(sym))
{ 
  property.add(symbol);
}

Atom::Atom(const PDBAtomRecord &p, const Cartesian &r) : 
  position(r), 
  velocity(0,0,0), 
  symbol(p.symbol()),
  type(symbol),
  pdb(p), 
  element(element_for(symbol))
{
  property.add(symbol);
}

void Atom::copy_from(const Atom &a)
{
  symbol = a.symbol;
  type = a.type;
  property = a.property;
  pdb = a.pdb;
  element = a.element;
}

bool Atom::is_hydrogen() const { return !strcmp(symbol,"H"); }
bool Atom::is_dummy() const { return !strcmp(symbol,"dummy"); }

void Atom::write_as_xyz(ostream &s) const
{
  s << symbol << "  " << position << " "
    << (strlen(type) > 0 ? (const char *) type : "*")
    << "  " << velocity << "\n";
}

void Atom::write_as_mol2(ostream &s) const
{
  s << symbol << "  " << position << " "
    << (strlen(type) > 0 ? (const char *) type : "*")
    << "\n";
}

void Atom::create_from_tinker(istream &s)
{
  int j;
  char tsym[64];
  s >> j >> tsym >> position >> type;
  symbol = element_for_tinker_symbol(tsym);
  element = element_for(symbol);
  property.remove_all();
  property.add(symbol);
  neighbors.remove_all();
  bond_types.remove_all();
  while (s >> j)
    neighbors.add(j-1);
}

const char *Atom::bond_type_for_neighbor(int j) const
{
  LstItr<int> n;
  LstItr<Str> t;
  for (n.init(neighbors), t.init(bond_types); n.ok(); n.next(), t.next()) {
    if (!t.ok())
      break;
    if (n() == j)
      return t();
  }
  return 0;
}

Complex Atom::nuclear_structure_factor(const Cartesian &k) const
{
  return e_to_charge_unit(atomic_number()) * expi(-k*position);
}
