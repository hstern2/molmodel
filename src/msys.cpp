#include <sstream>
#include <cctype>
#include "msys.hpp"
#include "htab.hpp"
#include "units.h"
#include "linalg.hpp"
#include "random.h"
#include "cvec.hpp"
#include "nbrlist.hpp"
#include "pbin.h"
#include "timing.h"
#include "euler.hpp"
#include "apattern.hpp"
#include "geograd.hpp"
#include "rmsd.hpp"
#include "fns.h"
#include "dvec.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif

MSys::MSys() : boundary_conditions(BoundaryConditions::new_non_periodic()) { }
MSys::MSys(const BoundaryConditions &bc) : boundary_conditions(bc.copy()) { }

MSys::MSys(const MSys &msys) : 
  boundary_conditions(msys.boundary_conditions->copy()),
  atom(msys.atom.copy()), 
  molecule(msys.molecule.copy()), 
  acidic_group(msys.acidic_group.copy()),
  properties(msys.properties), 
  types(msys.types),
  acid(msys.acid.copy())
{ }

MSys & MSys::operator=(const MSys &msys)
{
  if (boundary_conditions != msys.boundary_conditions) {
    delete boundary_conditions;
    boundary_conditions = msys.boundary_conditions->copy();
  }
  atom.copy(msys.atom);
  molecule.copy(msys.molecule);
  acidic_group.copy(msys.acidic_group);
  properties.copy(msys.properties);
  types.copy(msys.types);
  acid.copy(msys.acid);
  return *this;
}

MSys::~MSys() 
{
  delete boundary_conditions;
}

void MSys::append(const MSys &m)
{
  const int nat = atom.size(), nmol = molecule.size();
  int i;
  atom.resize(nat + m.atom.size());
  LstItrMod<int> j;
  for (i = 0; i < m.atom.size(); i++) {
    atom[i+nat] = m.atom[i];
    for (j.init(atom[i+nat].neighbors); j.ok(); j.next())
      j() += nat;
  }
  molecule.resize(nmol + m.molecule.size());
  for (i = 0; i < m.molecule.size(); i++)
    molecule[i+nmol] = Int2(m.molecule[i].a + nat, m.molecule[i].b + nat);
  init_acid();
}

void MSys::create_from_topology(istream &s)
{
  LstStr l;
  s >> l;
  LstStr tmp;
  LstItr<Str> t;
  for (t.init(l); t.ok(); t.next()) {
    istringstream ss((const char *) t());
    Str buf;
    while (ss) {
      ss >> buf;
      tmp.add(buf);
    }
  }
  create_from_topology(tmp);
}

void MSys::create_protein_from_sequence(Str seq)
{
  LstStr top;
  istream *f = StreamSearch("_msim_amino_acid_templates");
  HTab<LstStr> aa;
  LstStr start, end;
  *f >> aa >> start >> end;
  insist(start.size() > 0);
  insist(end.size() > 0);
  HSetStr starts, ends;
  LstItr<Str> stmp;
  for (stmp.init(start); stmp.ok(); stmp.next())
    starts.add(stmp());
  for (stmp.init(end); stmp.ok(); stmp.next())
    ends.add(stmp());
  bool default_start = true;
  bool default_end = true;
  char key[256];
  for (const char *c = seq; *c; c++) {
    if (!default_end)
      die("MSys::create_protein_from_sequence: residue '%s' "
	  "must only appear at end of sequence", key);
    if (isspace(*c))
      continue;
    if (isupper(*c)) {
      *key = *c;
      key[1] = '\0';
    } else {
      int i;
      for (i = 0; i < 255; i++, c++) {
	if (*c == '\0' || isspace(*c) || isupper(*c))
	  break;
	key[i] = *c;
      }
      if (isupper(*c))
	c--;
      insist(i < 256);
      key[i] = '\0';
    }
    if (starts.exists(key)) {
      if (top.size() != 0)
	die("MSys::residue '%s' must only appear "
	    "at start of sequence", key);
      default_start = false;
    }
    if (ends.exists(key)) {
      default_end = false;
    }
    LstStr *t = aa.get(key);
    if (!t)
      die("MSys::create_protein_from_sequence: unknown residue: '%s'", key);
    top.add(aa[key]);
  }
  if (default_start && start.size() > 0) {
    LstStr *t = aa.get(start.first());
    if (!t)
      die("MSys::create_protein_from_sequence: unknown residue: '%s'", key);
    LstStr tmp(*t);
    tmp.add(top);
    top = tmp;
  }
  if (default_end && end.size() > 0) {
    LstStr *t = aa.get(end.first());
    if (!t)
      die("MSys::create_protein_from_sequence: unknown residue: '%s'", key);
    top.add(*t);
  }
  create_from_topology(top);
}

void MSys::create_from_topology(const LstStr &ls)
{
  Lst<Atom> a;
  Lst<Int2> tmpn;
  Lst<Int2> g;
  LstItr<Str> s;
  LstInt istack;
  int i = -1, ilast = -1, istartmol = 0, j;
  HTab<int> ring;
  int first = -1;
  for (s.init(ls); s.ok(); s.next()) {
    Str tmp(s().copy());
    char *c = tmp;
    while (*c) {
      char tok = 0;
      char *ctmp = strpbrk(c,"()@;");
      if (ctmp) {
	tok = *ctmp;
	*ctmp = 0;
      }
      if (*c) {
	a.add(Atom(c));
	i++;
	tmpn.add(Int2(ilast,i));
	ilast = i;
	if (first == -1 && istack.size() == 0)
	  first = i;
      }
      if (!tok)
	break;
      c = ctmp+1;
      char jstr[64];
      switch (tok) {
      case '(': // start branch
	istack.push(ilast);
	continue;
      case ')': // end branch
	if (istack.size() == 0)
	  die("MSys::create_from_topology: unmatched ')'");
	ilast = istack.pop();
	continue;
      case '@': // ring closure
	if (isdigit(*c)) {
	  j = atoi(c);
	  insist(j >= 0);
	  while (isdigit(*c))
	    c++;
	} else {
	  j = -1;
	}
	snprintf(jstr,64,"%d",j);
	if (ring.exists(jstr)) {
	  tmpn.add(Int2(ring[jstr],ilast));
	  ring.remove(jstr);
	} else {
	  ring[jstr] = ilast;
	}
	continue;
      case ';': // molecule separator
	if (i+1 > istartmol) {
	  g.add(Int2(istartmol,i+1));
	  istartmol = i+1;
	}
	LstItrMod<Int2> ntmp;
	for (ntmp.init(tmpn); ntmp.ok(); ntmp.next()) {
	  if (ntmp().a == -1)
	    ntmp().a = first;
	  if (ntmp().b == -1)
	    ntmp().b = first;
	}
	first = -1;
	ilast = -1;
	istack.remove_all();
	continue;
      }
    }
  }
  if (ring.size() > 0)
    die("MSys::create_from_topology: all '@' must be paired");
  if (i+1 > istartmol)
    g.add(Int2(istartmol,i+1));
  atom.copy_from_list(a);
  LstItrMod<Int2> n;
  HSet<Int2> seen;
  for (n.init(tmpn); n.ok(); n.next()) {
    if (n().a == -1)
      n().a = first;
    if (n().b == -1)
      n().b = first;
    if (n().a == -1 || n().b == -1 || n().a == n().b || seen.exists(n()))
      continue;
    atom[n().a].neighbors.add(n().b);
    atom[n().b].neighbors.add(n().a);
    seen.add(n());
    seen.add(Int2(n().b,n().a));
  }
  molecule.copy_from_list(g);
  random_atomic_displacements(1);
  assign_properties_and_types();
  init_acid();
}

void MSys::replicate(int n)
{
  MSys m(*boundary_conditions);
  m.atom.resize(n * atom.size());
  m.molecule.resize(n * molecule.size());
  int ja = 0, jg = 0;
  for (int i = 0; i < n; i++) {
    for (int ig = 0; ig < molecule.size(); ig++)
      m.molecule[jg++] = Int2(molecule[ig].a + ja, molecule[ig].b + ja);
    const int offset = ja;
    for (int ia = 0; ia < atom.size(); ia++) {
      m.atom[ja] = atom[ia];
      LstItrMod<int> n;
      for (n.init(m.atom[ja].neighbors); n.ok(); n.next())
	n() += offset;
      ja++;
    }
  }
  *this = m;
}

void MSys::create_from_molecules(const MSys &m, const Lst<int> &imol)
{
  MSys mtmp, mtot(*m.boundary_conditions);
  LstItr<int> i;
  for (i.init(imol); i.ok(); i.next()) {
    mtmp.create_from_molecule(m,i());
    mtot.append(mtmp);
  }
  *this = mtot;
}

void MSys::create_from_molecule(const MSys &m, int imol)
{
  insist(0 <= imol && imol < m.molecule.size());
  delete boundary_conditions;
  boundary_conditions = m.boundary_conditions->copy();
  const int a = m.molecule[imol].a;
  const int b = m.molecule[imol].b;
  atom.resize(b-a);
  molecule.resize(1);
  molecule[0].a = 0;
  molecule[0].b = b-a;
  properties = m.properties;
  types = m.types;
  acid = m.acid.copy();
  int i, j;
  for (i = a, j = 0; i < b; i++, j++) {
    atom[j] = m.atom[i];
    LstItrMod<int> n;
    for (n.init(atom[j].neighbors); n.ok(); n.next())
      n() -= a;
  }
}

void MSys::set_density(double d)
{
  d = g_cm3_to_density_unit(d);
  if (boundary_conditions->type == non_periodic)
    die("MSys::set_density: boundary conditions must be periodic");
  scale(cbrt(density()/d));
}

void MSys::set_number_density(double d)
{
  if (boundary_conditions->type == non_periodic)
    die("MSys::set_number_density: boundary conditions must be periodic");
  scale(cbrt(number_density()/d));
}

void MSys::set_volume(double v)
{
  if (boundary_conditions->type == non_periodic)
    die("MSys::set_volume: boundary conditions must be periodic");
  scale(cbrt(v/boundary_conditions->volume()));
}

void MSys::set_molar_volume(double v)
{
  if (boundary_conditions->type == non_periodic)
    die("MSys::set_molar_volume: boundary conditions must be periodic");
  const int nmol = molecule.size();
  if (nmol == 0)
    die("MSys::set_molar_volume: no molecules");
  set_volume(L_mol_to_A3(v) * molecule.size());
}

double MSys::density() const
{
  if (boundary_conditions->type == non_periodic)
    return 0;
  double m = 0;
  for (int i = 0; i < atom.size(); i++)
    m += atom[i].mass();
  return m/boundary_conditions->volume();
}

double MSys::number_density() const
{
  return molecule.size()/boundary_conditions->volume();
}

Cartesian MSys::center_of_mass() const
{
  Cartesian cm(0,0,0);
  double m = 0;
  for (int i = 0; i < atom.size(); i++) {
    m += atom[i].mass();
    cm += atom[i].mass() * atom[i].position;
  }
  if (m > 0)
    cm /= m;
  return cm;
}

double MSys::mass() const
{
  double m = 0;
  for (int i = 0; i < atom.size(); i++)
    m += atom[i].mass();
  return m;
}

void MSys::translate_atoms_by_random_lattice_vectors()
{
  if (boundary_conditions->type == non_periodic)
    return;
  Tensor a = boundary_conditions->lattice_vectors();
  for (int i = 0; i < atom.size(); i++) {
    Cartesian j;
    j.x = round_to_nearest_integer(2*GaussianRandom());
    j.y = round_to_nearest_integer(2*GaussianRandom());
    j.z = round_to_nearest_integer(2*GaussianRandom());
    atom[i].position += a*j;
  }
}

void MSys::translate(const Cartesian &r)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].position += r;
}

void MSys::rezero_center_of_mass()
{
  translate(-center_of_mass());
}

void MSys::rezero_linear_momentum()
{
  const Cartesian vcm = linear_momentum()/mass();
  for (int i = 0; i < atom.size(); i++)
    atom[i].velocity -= vcm;
}

void MSys::scale(double x)
{
  for (int i = 0; i < molecule.size(); i++)
    translate_molecule(i, (x-1)*molecule_center_of_mass(i));
  boundary_conditions->scale(x);
}

void MSys::scale_centroids(double x)
{
#ifdef USE_MPI
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (int i = 0; i < molecule.size(); i++) {
    Cartesian rcm = molecule_center_of_mass(i);
    Cartesian rcm_centroid;
    MPI_Reduce(&rcm, &rcm_centroid, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    rcm_centroid /= size;
    MPI_Bcast(&rcm_centroid, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    translate_molecule(i, (x-1)*rcm_centroid);
  }
  boundary_conditions->scale(x);
#else
  die("MSys::scale_centroids: this executable was not compiled with MPI");
#endif
}

void MSys::translate_molecule(int j, const Cartesian &r)
{
  for (int i = molecule[j].a; i < molecule[j].b; i++)
    atom[i].position += r;
}

void MSys::translate_atom(int j, const Cartesian &r)
{
  atom[j].position += r;
}

Cartesian MSys::molecule_center_of_mass(int j) const
{
  Cartesian cm(0,0,0);
  double m = 0;
  for (int i = molecule[j].a; i < molecule[j].b; i++) {
    m += atom[i].mass();
    cm += atom[i].mass() * atom[i].position;
  }
  if (m > 0)
    cm /= m;
  return cm;
}

double MSys::molecule_mass(int j) const
{
  double m = 0;
  for (int i = molecule[j].a; i < molecule[j].b; i++)
    m += atom[i].mass();
  return m;
}

void MSys::tile(int n)
{
  if (n <= 0)
    die("MSys::tile: n must be > 0");
  if (boundary_conditions->type == non_periodic)
    die("MSys::tile: boundary conditions must be periodic");
  const Tensor a = boundary_conditions->lattice_vectors();
  const int ng = molecule.size(), n2 = n*n, n3 = n*n2;
  replicate(n3);
  for (int i0 = 0; i0 < n; i0++)
    for (int i1 = 0; i1 < n; i1++)
      for (int i2 = 0; i2 < n; i2++)
	for (int i = 0; i < ng; i++)
	  translate_molecule(i0*n2*ng + i1*n*ng + i2*ng + i, a*Cartesian(i0,i1,i2));
  boundary_conditions->scale(n);
}

void MSys::cubic_lattice(int n, double d)
{
  delete boundary_conditions;
  boundary_conditions = BoundaryConditions::new_cubic(1);
  if (is_perfect_cube(n)) {
    /* Simple cubic */
    tile((int) round_to_nearest_integer(cbrt(n)));
  } else if (is_perfect_cube(0.25*n)) {
    /* Face-centered cubic */
    Cartesian a[4] = { Cartesian(0,0,0), 
		       Cartesian(0,0.5,0.5), 
		       Cartesian(0.5,0,0.5),
		       Cartesian(0.5,0.5,0) };
    replicate(4);
    for (int i = 0; i < molecule.size(); i++)
      translate_molecule(i, a[i%4]);
    tile((int) round_to_nearest_integer(cbrt(0.25*n)));
  } else if (is_perfect_cube(0.5*n)) {
    /* Body-centered cubic */
    Cartesian a[2] = { Cartesian(0,0,0), 
		       Cartesian(0.5,0.5,0.5) };
    replicate(2);
    for (int i = 0; i < molecule.size(); i++)
      translate_molecule(i, a[i%2]);
    tile((int) round_to_nearest_integer(cbrt(0.5*n)));
  } else {
    die("MSys::cubic_lattice: n must be equal to\n"
	"N (for simple cubic),\n"
	"2N (for body-centered cubic), or\n"
	"4N (for face-centered cubic),\n"
	"where N is a perfect cube.");
  }
  set_density(d);
}

void MSys::create_from_pdb(CreateFromPDB p)
{
  Lst<Atom> alst;
  FILE *f = FileSearch(p.file);
  char buf[1024];
  while (fgets(buf, 1024, f)) {
    if (strlen(buf) >= 1024)
      die("Msys::create_from_pdb: error reading %s", (const char *) p.file);
    if (!strncmp(buf, "ATOM", 4) || !strncmp(buf,"HETATM",6)) {
      PDBAtomRecord pdb;
      Cartesian r;
      pdb.read(buf,r);
      if (p.ignore_alternate_locations && pdb.is_alternate_location())
	continue;
      if (p.ignore_hydrogens && pdb.is_hydrogen())
	continue;
      if (pdb.is_heterogen() && !p.include_heterogen.exists(pdb.resName))
	continue;
      alst.add(Atom(pdb,r));
    }
  }
  atom.copy_from_list(alst);
  find_neighbors();
  find_molecules();
  assign_properties_and_types();
  init_acid();
}

void MSys::write_as_pdb(Str fname) const
{
  ostream *s = FileStream(fname);
  write_as_pdb(*s);
  delete s;
}

void MSys::write_as_pdb(ostream &s) const
{
  int k = 1;
  for (int i = 0; i < atom.size(); i++)
    if (!atom[i].is_dummy() && !atom[i].pdb.is_heterogen())
      atom[i].pdb.write(s, atom[i].symbol, k++, atom[i].position);
  for (int i = 0; i < atom.size(); i++)
    if (!atom[i].is_dummy() && atom[i].pdb.is_heterogen())
      atom[i].pdb.write(s, atom[i].symbol, k++, atom[i].position);
  s << "END\n" << flush;
}

void MSys::create_from_xyz(Str s)
{
  FILE *f = FileSearch(s);
  char buf[256];
  int n;
  if (!fgets(buf, 256, f) || buf[strlen(buf)-1] != '\n')
    die("MSys::create_from_xyz: error reading xyz file");
  if (sscanf(buf,"%d",&n) != 1)
    die("MSys::create_from_xyz: could not read number of atoms");
  fclose(f);
  atom = Vec<Atom>(n);
  read_as_xyz(s);
  find_neighbors();
  find_molecules();
  assign_properties_and_types();
  init_acid();
}

static void look_for(const char *what, char *buf, int nbuf, FILE *f, const char *fn)
{
  while (fgets(buf,nbuf,f))
    if (!strncmp(buf, what, strlen(what)))
      return;
  die("%s: cannot find '%s'", fn, what);
}

static bool is_all_white_space(const char *buf)
{
  for ( ; *buf; buf++)
    if (!isspace(*buf))
      return false;
  return true;
}

struct SymTerm
{
  double a;
  int i;
  friend ostream & operator<<(ostream &s, const SymTerm &c) { return s << c.a << " " << c.i; }
};

struct SymOp
{
  Vec<Lst<SymTerm> > terms;
  SymOp(char *s) : terms(3)
  {
    int i = 0;
    for (char *c = strtok(s,","); i < 3; i++, c = strtok(0,",")) {
      while (c) {
	SymTerm t;
	t.a = 0;
	t.i = 0;
	int n, d;
	char sn[4], xyz[4];
	if (sscanf(c, "%d/%d", &n, &d) == 2) {
	  t.a = (double) n / (double) d;
	  t.i = -1;
	} else if ((n = sscanf(c, "%[+-]%[xyz]", sn, xyz)) == 2 ||
		   (n = sscanf(c, "%[xyz]", xyz)) == 1) {
	  t.a = n == 2 && *sn == '-' ? -1.0 : 1.0;
	  if (*xyz == 'x')
	    t.i = 0;
	  else if (*xyz == 'y')
	    t.i = 1;
	  else if (*xyz == 'z')
	    t.i = 2;
	  else
	    die("SymOp: cannot parse %s", c);
	} else {
	  die("SymOpt: cannot parse %s", c);
	}
	terms[i].add(t);
	c++;
	c = strpbrk(c, "+-");
      }
    };
  }
  Cartesian new_position(const Cartesian &c0) const
  {
    const double *r0 = (const double *) &c0;
    Cartesian c(0,0,0);
    double *r = (double *) &c;
    for (int i = 0; i < 3; i++) {
      LstItr<SymTerm> t;
      for (t.init(terms[i]); t.ok(); t.next())
	r[i] += t().a * (t().i == -1 ? 1.0 : r0[t().i]);
    }
    return c;
  }
  friend ostream & operator<<(ostream &s, const SymOp &c) { return s << c.terms; }  
};

void MSys::create_from_cif(Str fname)
{
  const char *fn = fname;
  FILE *f = FileSearch(fn);
  char buf[1024];
  look_for("_symmetry_equiv_pos_as_xyz", buf, 1024, f, fn);
  Lst<SymOp> symm;
  while (fgets(buf, 1024, f)) {
    int d;
    if (!strncmp(buf, "_cell_length_a", 14))
      goto found_cell_length_a;
    char sbuf[1024];
    if (sscanf(buf,"%d %s",&d,sbuf) != 2)
      die("MSys::create_from_cif: %s: error reading symmetry lines", fn);
    symm.add(SymOp(sbuf));
  }
  die("MSys::create_from_cif: %s: cannot find '_cell_length_a'", fn);
 found_cell_length_a:
  double a=0, b=0, c=0, alpha=0, beta=0, gamma=0;
  if (sscanf(buf,"_cell_length_a %lf", &a) != 1)
    die("MSys::create_from_cif: %s: error reading '_cell_length_a'", fn);
  if (!(fgets(buf, 1024, f) && sscanf(buf,"_cell_length_b %lf", &b) == 1))
    die("MSys::create_from_cif: %s: cannot find '_cell_length_b'", fn);
  if (!(fgets(buf, 1024, f) && sscanf(buf,"_cell_length_c %lf", &c) == 1))
    die("MSys::create_from_cif: %s: cannot find '_cell_length_c'", fn);
  if (!(fgets(buf, 1024, f) && sscanf(buf,"_cell_angle_alpha %lf", &alpha) == 1))
    die("MSys::create_from_cif: %s: cannot find '_cell_angle_alpha'", fn);
  if (!(fgets(buf, 1024, f) && sscanf(buf,"_cell_angle_beta %lf", &beta) == 1))
    die("MSys::create_from_cif: %s: cannot find '_cell_angle_beta'", fn);
  if (!(fgets(buf, 1024, f) && sscanf(buf,"_cell_angle_gamma %lf", &gamma) == 1))
    die("MSys::create_from_cif: %s: cannot find '_cell_angle_gamma'", fn);
  delete boundary_conditions;
  boundary_conditions = BoundaryConditions::new_from_cell(a,b,c,alpha,beta,gamma);
  look_for("_cell_formula_units_Z", buf, 1024, f, fn);
  int formula_units;
  if (sscanf(buf, "_cell_formula_units_Z %d", &formula_units) != 1)
    die("MSys::create_from_cif: %s: error reading '_cell_formula_units_Z'", fn);    
  look_for("_atom_site_fract_z", buf, 1024, f, fn);
  Lst<Atom> alst;
  while (fgets(buf, 1024, f) && !is_all_white_space(buf)) {
    char sym[256];
    double x, y, z;
    if (sscanf(buf, "%*s %s %lf %lf %lf", sym, &x, &y, &z) != 4)
      die("MSys::create_from_cif: %s: error reading coordinate line %s", fn, buf);
    Atom a(sym);
    a.position.x = x;
    a.position.y = y;
    a.position.z = z;
    alst.add(a);
  }
  atom = Vec<Atom>(formula_units * alst.size());
  LstItr<SymOp> o;
  o.init(symm);
  int n = 0;
  const Tensor h = boundary_conditions->lattice_vectors();
  for (int i = 0; i < formula_units; i++, o.next()) {
    insist(o.ok());
    LstItr<Atom> a;
    for (a.init(alst); a.ok(); a.next(), n++) {
      atom[n] = a();
      atom[n].position = h * o().new_position(a().position);
    }
  }
  find_neighbors();
  find_molecules();
  assign_properties_and_types();
  init_acid();
}

void MSys::create_from_csd(Str s)
{
  FILE *f = FileSearch(s);
  char buf[256];
  /* First line is a comment - skip it */
  if (!fgets(buf, 256, f) || buf[strlen(buf)-1] != '\n')
    die("MSys::create_from_xyz: error reading csd file '%s'", (const char *) s);
  /* Read CELL description */
  while (fgets(buf, 256, f) && strncmp(buf, "CELL", 4))
    ;
  if (strncmp(buf, "CELL", 4))
    die("MSys::create_from_csd: cannot find 'CELL' in file %s", (const char *) s);
  double a, b, c, alpha, beta, gamma;
  if (sscanf(buf,"CELL %lf%lf%lf%lf%lf%lf",&a,&b,&c,&alpha,&beta,&gamma) != 6)
    die("MSys::create_from_csd: expecting a,b,c,alpha,beta,gamma "
	"in 'CELL' record in file %s", (const char *) s);
  delete boundary_conditions;
  const Tensor bc = boundary_conditions->lattice_vectors();
  Lst<Tensor> rot;
  Lst<Cartesian> trans, xyz;
  Lst<Str> sym;
  /* Read SYMM or atom records */
  while (fgets(buf, 256, f)) 
    if (buf[strlen(buf)-1] != '\n')
      die("MSys::create_from_xyz: error reading csd file '%s'", (const char *) s);
    else if (!strncmp(buf, "SYMM", 4)) {
      Tensor r;
      Cartesian t;
      if (sscanf(buf, "SYMM %lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	 	 &r.xx,&r.xy,&r.xz,&t.x,
		 &r.yx,&r.yy,&r.yz,&t.y,
		 &r.zx,&r.zy,&r.zz,&t.z) != 12)
	die("MSys::create_from_csd: expecting 12 numbers "
	    "in 'SYMM' record in file %s", (const char *) s);
      rot.add(r);
      trans.add(t);
    } else {
      char sy[256];
      double x,y,z;
      if (sscanf(buf, "%s%lf%lf%lf",sy,&x,&y,&z) != 4)
	die("MSys::create_from_csd: expecting a string and 3 numbers "
	    "in atom record in file %s", (const char *) s);
      for (char *ctmp = sy + strlen(sy) - 1; ctmp >= sy; ctmp--)
	if (isdigit(*ctmp))
	  *ctmp = '\0';
      sym.add(sy);
      xyz.add(Cartesian(x,y,z));
    }
  fclose(f);
  atom = Vec<Atom>(rot.size() * xyz.size());
  LstItr<Tensor> r;
  LstItr<Cartesian> t, x;
  LstItr<Str> sy;
  int n = 0;
  for (r.init(rot), t.init(trans); r.ok(); r.next(), t.next())
    for (x.init(xyz), sy.init(sym); x.ok(); x.next(), sy.next(), n++) {
      atom[n] = Atom(sy());
      atom[n].position = bc*(r()*x() + t());
    }
  insist(n == atom.size());
  find_neighbors();
  find_molecules();
  assign_properties_and_types();
  init_acid();
}

void MSys::create_from_tinker(Str s)
{
  FILE *f = FileSearch(s);
  char buf[512];
  if (!fgets(buf,512,f) || buf[strlen(buf)-1] != '\n')
    die("MSys::create_from_tinker: error reading Tinker file");
  int n;
  if (sscanf(buf,"%d",&n) != 1)
    die("MSys::create_from_tinker: could not read number of atoms");
  if (n <= 0)
    return;
  atom = Vec<Atom>(n);
  for (int i = 0; i < n; i++) {
    if (!fgets(buf,512,f) || buf[strlen(buf)-1] != '\n')
      die("MSys::create_from_tinker: error reading Tinker file");
    istringstream ss(buf);
    atom[i].create_from_tinker(ss);
  }
  check_neighbors();
  find_molecules();
  init_acid();
}

void MSys::read_as_xyz(Str fname)
{
  FILE *f = FileSearch(fname);
  read_as_xyz(f);
  fclose(f);
}

void MSys::read_as_xyz(FILE *f)
{
  char buf[256];
  int n;
  const bool alldum = contains_only_dummy_atoms();
  if (!fgets(buf, 256, f) || buf[strlen(buf)-1] != '\n')
    die("MSys::read_as_xyz: error reading line 1 of xyz file");
  if (sscanf(buf,"%d",&n) != 1)
    die("MSys::read_as_xyz: could not read number of atoms");
  if (n != atom.size())
    die("MSys::read_as_xyz: number of atoms in molsys (%d) "
	"does not match number of atoms in xyz file (%d)",
	atom.size(), n);
  if (!fgets(buf, 256, f) || buf[strlen(buf)-1] != '\n')
    die("MSys::read_as_xyz: error reading line 2 of xyz file");
  char *tmp = strstr(buf,"boundary_conditions");
  if (tmp) {
    istringstream s(tmp + strlen("boundary_conditions"));
    s >> boundary_conditions;
  }
  for (int i = 0; i < n; i++) {
    if (!fgets(buf, 256, f) || buf[strlen(buf)-1] != '\n')
      die("MSys::read_as_xyz: error reading line %d of xyz file", i+2);
    Cartesian r, v;
    char sym[256], type[256];
    const int m = sscanf(buf,"%s%lf%lf%lf%s%lf%lf%lf",
			 sym,&r.x,&r.y,&r.z,type,&v.x,&v.y,&v.z);
    if (!(m == 4 || m == 5 || m == 8))
      die("MSys::read_as_xyz: error reading line %d of xyz file", i+2);
    if (!alldum) {
      if (strcmp(sym, atom[i].symbol))
	die("MSys::read_as_xyz: symbols (%s,%s) do not match for atom %d",
	    sym, (const char *) atom[i].symbol, i+1);
      if (m >= 5 && strcmp(type,"*") && strcmp(type, atom[i].type))
	die("MSys::read_as_xyz: types (%s,%s) do not match for atom %d",
	    type, (const char *) atom[i].type, i+1);
    } else {
      atom[i] = Atom(sym);
    }
    atom[i].position = r;
    if (m == 8)
      atom[i].velocity = v;
  }
}

void MSys::write_as_mol2(Str f) const
{
  ostream *stmp = FileStream(f);
  write_as_mol2(*stmp);
  delete stmp;
}

void MSys::write_as_mol2(ostream &s) const
{
  s << "@<TRIPOS>ATOM\n";
  int i;
  for (i = 0; i < atom.size(); i++) {
    s << i+1 << " ";
    atom[i].write_as_mol2(s);
  }
  s << "@<TRIPOS>BOND\n";
  int ib = 1; 
  for (i = 0; i < atom.size(); i++) {
    LstItr<int> n;
    LstItr<Str> bt;
    for (bt.init(atom[i].bond_types),n.init(atom[i].neighbors); n.ok(); n.next()) { 
      if (n() > i)
        s << ib++ << " " << i+1 << " " << n()+1 << " " << (bt.ok() ? (const char *) bt() : "1") << "\n";
      if (bt.ok())
        bt.next();
    } 
  }
  if (boundary_conditions->type != non_periodic) {
    s << "\n@<TRIPOS>CRYSIN\n";
    const Tensor h = boundary_conditions->lattice_vectors();
    const Cartesian zero(0,0,0);
    const double a = h.col(0).magnitude();
    const double b = h.col(1).magnitude();
    const double c = h.col(2).magnitude();
    const double alpha = Angle(h.col(1), zero, h.col(2));
    const double beta = Angle(h.col(0), zero, h.col(2));
    const double gamma = Angle(h.col(0), zero, h.col(1)); 
    s << a << " " << b << " " << c << "  ";
    s << radians_to_degrees(alpha) << " "
      << radians_to_degrees(beta) << " "
      << radians_to_degrees(gamma) << "\n";
  }
}

void MSys::write_as_xyz(ostream &s) const
{
  s << atom.size() << "\nboundary_conditions " << boundary_conditions << "\n";
  for (int i = 0; i < atom.size(); i++)
    atom[i].write_as_xyz(s);
  s << flush;
}

void MSys::write_as_xyz(ostream &s, double t) const
{
  char cmt[1024];
  snprintf(cmt, 1024, "Time elapsed: %f psec", t);
  write_as_xyz(s,cmt);
}

void MSys::write_as_xyz(ostream &s, const char *cmt) const
{
  s << atom.size() << "\nboundary_conditions " << boundary_conditions
    << "  " << cmt << "\n";
  for (int i = 0; i < atom.size(); i++)
    atom[i].write_as_xyz(s);
  s << flush;
}

void MSys::show_center_of_mass() const
{
  Out() << "Center of mass (A): " << center_of_mass() << "\n" << flush;
}

void MSys::show_volume() const
{
  Out() << "Volume: " << boundary_conditions->volume() << " A^3\n" << flush;
}

void MSys::show_lattice_vectors() const
{
  Out() << "Lattice vectors: " << boundary_conditions->lattice_vectors() << "\n" << flush;
}

void MSys::show_reciprocal_lattice_vectors() const
{
  Out() << "Reciprocal lattice vectors: " << boundary_conditions->reciprocal_lattice_vectors() << "\n" << flush;
}

void MSys::show_properties() const
{
  for (int i = 0; i < atom.size(); i++) {
    OutPrintf("%-5s ", (const char *) atom[i].symbol);
    HSetItr<Str> p;
    for (p.init(atom[i].property); p.ok(); p.next())
      if (strcmp(p(),atom[i].symbol))
	Out() << p() << " ";
    Out() << "\n";
  }
}

void MSys::show_acid() const
{
  Out() << acid.size() << " acidic groups...\n" << acid << "\n" << flush;
}

void MSys::show_mass() const
{
  Out() << "Masses:\n";
  for (int i = 0; i < atom.size(); i++)
    Out() << atom[i].symbol << " " << mass_unit_to_g_mol(atom[i].mass()) << "\n";
  Out() << "\n" << flush;
}

void MSys::show_geometry() const
{
  Out() << "Bond lengths (A):\n";
  int i;
  Lst<Int2> bonds;
  for (i = 0; i < atom.size(); i++) {
    LstItr<int> n;
    for (n.init(atom[i].neighbors); n.ok(); n.next())
      if (i < n()) {
	Out() << atom[i].symbol<< i << "-" 
	      << atom[n()].symbol<< n() << " "
	      << atom[i].position.distance(atom[n()].position)
	      << "\n";
	bonds.add(Int2(i,n()));
      }
  }
  Out() << "Angles (degrees):\n";
  for (i = 0; i < atom.size(); i++) {
    LstPairItr<int> n;
    for (n.init(atom[i].neighbors); n.ok(); n.next())
      Out() << atom[n.i()].symbol<< n.i() << "-"
	    << atom[i].symbol<< i << "-"
	    << atom[n.j()].symbol << n.j() << " "
	    << radians_to_degrees(Angle(atom[n.i()].position,atom[i].position,atom[n.j()].position))
	    << "\n";
  }
  Out() << "Dihedrals (degrees):\n";
  LstItr<int> ni, nl;
  LstItr<Int2> b;
  for (b.init(bonds); b.ok(); b.next()) {
    const int j = b().a, k = b().b;
    for (ni.init(atom[j].neighbors); ni.ok(); ni.next())
      if (ni() != k)
	for (nl.init(atom[k].neighbors); nl.ok(); nl.next())
	  if (nl() != j && nl() != ni())
	    Out() << atom[ni()].symbol<< ni() << "-"
		  << atom[j].symbol<< j << "-"
		  << atom[k].symbol<< k << "-"
		  << atom[nl()].symbol << nl() << " "
		  << radians_to_degrees(Dihedral(atom[ni()].position,
						 atom[j].position,
						 atom[k].position,
						 atom[nl()].position))
		  << "\n";
  }
}

void MSys::geometry_rmsd(Str xyzfile)
{
  const int n = atom.size();
  CVec rsave(n);
  copy_positions_to(rsave);
  read_as_xyz(xyzfile);
  double bsum = 0, asum = 0, dsum = 0;
  int bn = 0, an = 0, dn = 0;
  int i;
  Lst<Int2> bonds;
  for (i = 0; i < atom.size(); i++) {
    LstItr<int> n;
    for (n.init(atom[i].neighbors); n.ok(); n.next())
      if (i < n()) {
	bsum += sq(atom[i].position.distance(atom[n()].position) -
		   rsave[i].distance(rsave[n()]));
	bn++;
	bonds.add(Int2(i,n()));
      }
  }
  Out() << "Bond length RMSD: " << (bn > 0 ? sqrt(bsum/bn) : 0.0) << " A\n";
  for (i = 0; i < atom.size(); i++) {
    LstPairItr<int> n;
    for (n.init(atom[i].neighbors); n.ok(); n.next()) {
      asum += sq(radians_to_degrees(Angle(atom[n.i()].position,
					  atom[i].position,
					  atom[n.j()].position)) -
		 radians_to_degrees(Angle(rsave[n.i()],
					  rsave[i],
					  rsave[n.j()])));
      an++;
    }
  }
  Out() << "Angle RMSD: " << (an > 0 ? sqrt(asum/an) : 0.0) << " degrees\n";  
  LstItr<int> ni, nl;
  LstItr<Int2> b;
  for (b.init(bonds); b.ok(); b.next()) {
    const int j = b().a, k = b().b;
    for (ni.init(atom[j].neighbors); ni.ok(); ni.next())
      if (ni() != k)
	for (nl.init(atom[k].neighbors); nl.ok(); nl.next())
	  if (nl() != j && nl() != ni()) {
	    dsum += sq(radians_to_degrees(Dihedral(atom[ni()].position,atom[j].position,
						   atom[k].position,atom[nl()].position)) -
		       radians_to_degrees(Dihedral(rsave[ni()],rsave[j],
						   rsave[k],rsave[nl()])));
	    dn++;
	  }
  }
  Out() << "Dihedral RMSD: " << (dn > 0 ? sqrt(dsum/dn) : 0.0) << " degrees\n";    
  copy_positions_from(rsave);
}

void MSys::write_as_xyz(Str f) const
{
  ostream *s = FileStream(f);
  write_as_xyz(*s);
  delete s;
}

void MSys::write_as_xyz(Str f, Str cmt) const
{
  ostream *s = FileStream(f);
  write_as_xyz(*s,cmt);
  delete s;
}

void MSys::append_as_xyz(Str f) const
{
  ostream *s = AppendFileStream(f);
  write_as_xyz(*s);
  delete s;
}

#define LPAREN -1
#define RPAREN -2

static void write_atom_as_topology(int i, int ilast, const Atom *atom, 
				   Lst<int> &tok, int *seen, Vec<HSet<int> > &extrabond)
{
  if (seen[i]) {
    extrabond[i].add(ilast);
    extrabond[ilast].add(i);
    return;
  }
  tok.add(i);
  seen[i] = 1;
  LstItr<int> n;
  int nn = 0;
  for (n.init(atom[i].neighbors); n.ok(); n.next())
    if (n() != ilast)
      nn++;
  for (n.init(atom[i].neighbors); n.ok(); n.next())
    if (n() != ilast) {
      if (nn > 1)
	tok.add(LPAREN);
      write_atom_as_topology(n(),i,atom,tok,seen,extrabond);
      if (nn > 1)
	tok.add(RPAREN);
      nn--;
    }
}

void MSys::write_as_topology(ostream &s) const
{
  const int nmol = molecule.size();
  const int nat = atom.size();
  Vec<int> seen(nat, 0);
  Vec<HSet<int> > extrabond(nat);
  Vec<HSet<int> > ringmark(nat);
  HSet<int> used_ringmark;
  s << "{" << (nmol > 1 ? "\n" : " ");
  for (int i = 0; i < nmol; i++) {
    Lst<int> tok;
    for (int j = molecule[i].a; j < molecule[i].b; j++) {
      write_atom_as_topology(j,-1,atom,tok,seen,extrabond);
      break;
    }
    if (tok.size() == 0)
      continue;
    LstItr<int> t;
    int sp = 0, nparen = 0;
    for (t.init(tok); t.ok(); t.next()) {
      if (t() == LPAREN) {
	nparen++;
      } else if (t() == RPAREN) {
	if (nparen > 0)
	  nparen--;
	else {
	  s << ")";
	  sp = 0;
	}
      } else {
	while (nparen > 0) {
	  s << "(";
	  nparen--;
	  sp = 0;
	}
	if (sp)
	  s << " ";
	s << atom[t()].symbol;
	sp = 1;
	HSetItr<int> r;
	for (r.init(extrabond[t()]); r.ok(); r.next()) {
	  int k = 0;
	  while (used_ringmark.exists(k))
	    k++;
	  s << "@";
	  if (k > 0)
	    s << k;
	  used_ringmark.add(k);
	  ringmark[r()].add(k);
	  extrabond[r()].remove(t());
	}
	for (r.init(ringmark[t()]); r.ok(); r.next()) {
	  s << "@";
	  if (r() > 0)
	    s << r();
	  used_ringmark.remove(r());
	}
      }
    }
    if (i < nmol-1)
      s << ";\n";
  }
  s << (nmol > 1 ? "\n" : " ") << "}\n";
}

void MSys::write_as_topology(Str f) const
{
  ostream *s = FileStream(f);
  write_as_topology(*s);
  delete s;
}

void MSys::show_as_topology() const
{
  write_as_topology(Out());
}

void MSys::show_neighbors() const
{
  for (int i = 0; i < atom.size(); i++) {
    Out() << i << " " << atom[i].symbol << ":" << atom[i].type << " ";
    LstItr<int> j;
    for (j.init(atom[i].neighbors); j.ok(); j.next())
      Out() << j() << " ";
    Out() << "\n";
  }
}

void MSys::show_as_xyz() const
{
  write_as_xyz(Out());
}

void MSys::show_as_pdb() const
{
  write_as_pdb(Out());
}

void MSys::copy_positions_from(const Cartesian *r)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = r[i];
}

void MSys::copy_positions_to(Cartesian *r) const
{
  for (int i = 0; i < atom.size(); i++)
    r[i] = atom[i].position;
}

void MSys::copy_positions_from(const MSys &m)
{
  insist(atom.size() == m.atom.size());
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = m.atom[i].position;
}

void MSys::copy_velocities_from(const Cartesian *r)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].velocity = r[i];
}

void MSys::copy_velocities_to(Cartesian *r) const
{
  for (int i = 0; i < atom.size(); i++)
    r[i] = atom[i].velocity;
}

double MSys::kinetic_energy() const
{
  double k = 0;
  for (int i = 0; i < atom.size(); i++)
    k += atom[i].mass() * atom[i].velocity.sq();
  return k/2;
}

void MSys::integrate_velocities(const Cartesian *f, double dt)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].velocity += f[i] * (dt/atom[i].mass());
}

void MSys::integrate_positions(double dt)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].position += atom[i].velocity * dt;
}

Cartesian MSys::linear_momentum() const
{
  Cartesian P(0,0,0);
  for (int i = 0; i < atom.size(); i++)
    P += atom[i].linear_momentum();
  return P;
}

Cartesian MSys::angular_momentum() const
{
  Cartesian L(0,0,0);
  const Cartesian cm = center_of_mass();
  for (int i = 0; i < atom.size(); i++)
    L += atom[i].mass()*(atom[i].position-cm).cross(atom[i].velocity);
  return L;
}

void MSys::align_principal_axes_with_xyz()
{
  Tensor t;
  Cartesian c;
  Diagonalize(inertia_tensor(), c, t);
  t = t.transpose();
  if (t.determinant() < 0)
    t *= -1.0;
  const Cartesian cm = center_of_mass();
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = t*(atom[i].position - cm) + cm;
}

Tensor MSys::inertia_tensor() const
{
  Tensor d(0);
  const Cartesian cm = center_of_mass();
  for (int i = 0; i < atom.size(); i++) {
    const Cartesian s = atom[i].position - cm;
    d += atom[i].mass() * (Tensor(s.sq()) - Tensor(s,s));
  }
  return d;
}

void MSys::rezero_angular_momentum()
{
  Tensor I = inertia_tensor();
  InvertSymmetricTensor(I);
  const Cartesian omega = angular_momentum() * I; // angular velocity
  const Cartesian cm = center_of_mass();
  for (int i = 0; i < atom.size(); i++)
    atom[i].velocity -= omega.cross(atom[i].position-cm);
}

void MSys::map_to_central_box()
{
  for (int i = 0; i < molecule.size(); i++) {
    const Cartesian cm = molecule_center_of_mass(i);
    Cartesian r = cm;
    boundary_conditions->map_to_central_box(r);
    translate_molecule(i, r-cm);
  }
}

void MSys::randomize_orientations()
{
  double t;
  for (int i = 0; i < molecule.size(); i++) {
    t = UniformRandom() * 2*M_PI;
    Tensor u1(cos(t), sin(t), 0,
	      -sin(t), cos(t), 0,
	      0,       0,    1);
    t = UniformRandom() * 2*M_PI;
    Tensor u2(cos(t), 0, sin(t),
	      0,      1, 0,
	      -sin(t),0, cos(t));
    t = UniformRandom() * 2*M_PI;
    Tensor u3(1,      0,  0,
	      0,    cos(t), sin(t),
	      0,   -sin(t), cos(t));
    rotate_molecule(i, u1*u2*u3);
  }
}

void MSys::rotate_molecule(int i, double phi, double theta, double psi)
{
  Tensor r;
  EulerAnglesToRotationMatrix(r, phi, theta, psi);
  rotate_molecule(i,r);
}

void MSys::rotate_molecule(int i, const Tensor &u)
{
  const Cartesian cm = molecule_center_of_mass(i);
  translate_molecule(i,-cm);
  for (int j = molecule[i].a; j < molecule[i].b; j++)
    atom[j].position = u * atom[j].position;
  translate_molecule(i,cm);
}

void MSys::rotate(double phi, double theta, double psi)
{
  Tensor r;
  EulerAnglesToRotationMatrix(r, phi, theta, psi);
  const Cartesian cm = center_of_mass();
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = r*(atom[i].position - cm) + cm;
}

void MSys::apply_rotation_matrix(const Tensor r)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = r*atom[i].position;
}

void MSys::apply_symmetry_operations(LstRotTrans ops)
{
  const int nat = atom.size();
  replicate(ops.size());
  LstItr<RotTrans> o;
  int i = 0;
  for (o.init(ops); o.ok(); o.next()) {
    const Tensor &r = o().rotation_matrix;
    const Cartesian &t = o().translation;
    for (int j = 0; j < nat; i++, j++)
      atom[i].position = r*atom[i].position + t;
  }
}

void MSys::randomize_positions(unsigned seed)
{
  for (int i = 0; i < atom.size(); i++) {
    double c[3];
    for (int j = 0; j < 3; j++) {
      seed = (1103515245*seed + 12345);
      c[j] = ((seed/65536) % 32768) * boundary_conditions->min_diameter()/32769;
    }
    atom[i].position.x = c[0];
    atom[i].position.y = c[1];
    atom[i].position.z = c[2];
  }
}

void MSys::randomize_positions()
{
  const double d = boundary_conditions->min_diameter();
  for (int i = 0; i < atom.size(); i++)
    atom[i].position = d*Cartesian(UniformRandom() - 0.5,
				   UniformRandom() - 0.5,
				   UniformRandom() - 0.5);
}

void MSys::random_displacements(double r)
{
  for (int i = 0; i < molecule.size(); i++)
    translate_molecule(i, r*Cartesian(GaussianRandom(),GaussianRandom(),GaussianRandom()));
}

void MSys::random_atomic_displacements(double r)
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].position += r*Cartesian(GaussianRandom(),GaussianRandom(),GaussianRandom());
}

void MSys::check_types()
{
  for (int i = 0; i < atom.size(); i++)
    if (!strcmp(atom[i].symbol, atom[i].type))
      Out() << "*** atom " << i << " does not have a type assigned:\n" << atom[i] << "\n";
}

void MSys::check_distances(double r)
{
  const double r2 = r;
  for (int i = 0; i < atom.size(); i++)
    for (int j = i+1; j < atom.size(); j++) {
      const Cartesian a = atom[i].position;
      const Cartesian b = atom[j].position;
      const double c2 = boundary_conditions->square_minimum_image_distance(a,b);
      if (c2 < r2)
	Out() << i << "," << j << ": " << sqrt(c2) << "\n" << flush;
    }
}

void MSys::arrange()
{
  const int nmol = molecule.size();
  const int nl = (int) (is_perfect_cube(nmol) ?
			floor(cbrt(nmol) + 0.5) : 
			ceil(cbrt(nmol)));
  const int nl3 = nl*nl*nl;
  const int nskip = nl3 - nmol;
  Vec<int> skip(nl3, 0);
  int i, j, k;
  for (i = 0; i < nskip; i++) {
    int iskip = (int) (UniformRandom() * nl3);
    while (skip[iskip%nl3])
      iskip++;
    skip[iskip%nl3] = 1;
  }
  int imol = 0;
  const Tensor a = boundary_conditions->lattice_vectors();
  for (i = 0; i < nl; i++)
    for (j = 0; j < nl; j++)
      for (k = 0; k < nl; k++) {
	if (skip[(i*nl+j)*nl+k])
	  continue;
	translate_molecule(imol, a*Cartesian((double) i / (double) nl,
					     (double) j / (double) nl,
					     (double) k / (double) nl) - 
			   molecule_center_of_mass(imol));
	imol++;
      }
  insist(imol == nmol);
  rezero_center_of_mass();
}

void MSys::show_distances_within(double w)
{
  const int n = atom.size();
  const BoundaryConditions &bc = *boundary_conditions;
  Out() << "Close contacts:\n";
  CVec r(n);
  DVec c(n);
  int i;
  for (i = 0; i < n; i++) {
    r[i] = atom[i].position;
    c[i] = atom[i].vdw_radius();
  }
  NeighborList nlist(n, w, 0);
  nlist.set_boundary_conditions(bc);
  nlist.coordinates_changed(r);
  for (i = 0; i < n; i++) {
    const int nj = nlist.number_of_neighbors()[i];
    const int *j = nlist.neighbors()[i];
    for (int k = 0; k < nj; k++) {
      const double r2 = bc.square_minimum_image_distance(r[i], r[j[k]]);
      if (r2 < w*w) {
	LstItr<int> q;
	for (q.init(atom[i].neighbors); q.ok(); q.next())
	  if (q() == j[k])
	    break;
	if (!q.ok())
	  Out() << sqrt(r2) << " A:\n" 
		<< i << " " << atom[i]
		<< j[k] << " " << atom[j[k]] << "\n";
      }
    }
  }
  Out() << flush;
}

static HSet<HSet<int> > clusters(const struct PairBinary *p, int n, int nc)
{
  HSet<HSet<int> > c;
  if (nc > 2) {
    HSet<HSet<int> > cm = clusters(p,n,nc-1);
    Vec<Lst<const HSet<int> *> > containing(n);
    HSetItr<HSet<int> > cmi;
    for (cmi.init(cm); cmi.ok(); cmi.next()) {
      HSetItr<int> j;
      for (j.init(cmi()); j.ok(); j.next())
	containing[j()].add(&cmi());
    }
    for (int i = 0; i < n; i++) {
      const int nj = PairBinary_number_of_elements(p)[i];
      const int *j = PairBinary_elements(p)[i];
      for (int k = 0; k < nj; k++) {
	LstItr<const HSet<int> *> cmi;
	for (cmi.init(containing[j[k]]); cmi.ok(); cmi.next()) {
	  if (cmi()->exists(i))
	    continue;
	  HSetItr<int> q;
	  HSet<int> cn;
	  for (q.init(*cmi()); q.ok(); q.next())
	    if (!PairBinary_exists(p,i,q()))
	      goto next_cluster;
	  cn.add(*cmi());
	  cn.add(i);
	  c.add(cn);
	next_cluster:
	  continue;
	}
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      const int nj = PairBinary_number_of_elements(p)[i];
      const int *j = PairBinary_elements(p)[i];
      for (int k = 0; k < nj; k++) {
	HSet<int> c2;
	c2.add(i);
	c2.add(j[k]);
	c.add(c2);
      }
    }
  }
  return c;
}

void MSys::show_clusters_as_xyz(int n, double w)
{
  const int nmol = molecule.size();
  const BoundaryConditions &bc = *boundary_conditions;
  CVec r(nmol);
  int i;
  for (i = 0; i < nmol; i++)
    r[i] = molecule_center_of_mass(i);
  NeighborList nlist(nmol,w,0);
  nlist.set_boundary_conditions(bc);
  nlist.coordinates_changed(r);
  struct PairBinary *nbr = PairBinary_new(nmol);
  const double w2 = w*w;
  for (i = 0; i < nmol; i++) {
    const int nj = nlist.number_of_neighbors()[i];
    const int *j = nlist.neighbors()[i];
    for (int k = 0; k < nj; k++)
      if (bc.square_minimum_image_distance(r[i],r[j[k]]) < w2)
	PairBinary_add(nbr, i, j[k]);
  }
  PairBinary_init(nbr);
  const HSet<HSet<int> > c = clusters(nbr,nmol,n);
  HSetItr<HSet<int> > ci;
  for (ci.init(c); ci.ok(); ci.next()) {
    Lst<int> lci;
    HSetItr<int> j;
    for (j.init(ci()); j.ok(); j.next())
      lci.add(j());
    MSys m;
    m.create_from_molecules(*this,lci);
    m.rezero_center_of_mass();
    m.map_to_central_box();
    m.rezero_center_of_mass();
    m.show_as_xyz();
  }
  Out() << flush;
}

void MSys::show_intermolecular_distances_within(double w)
{
  const int n = atom.size();
  const BoundaryConditions &bc = *boundary_conditions;
  Out() << "Close contacts:\n";
  CVec r(n);
  DVec c(n);
  int i;
  for (i = 0; i < n; i++) {
    r[i] = atom[i].position;
    c[i] = atom[i].vdw_radius();
  }
  NeighborList nlist(n, w, 0);
  nlist.set_boundary_conditions(bc);
  nlist.coordinates_changed(r);
  for (i = 0; i < n; i++) {
    const int nj = nlist.number_of_neighbors()[i];
    const int *j = nlist.neighbors()[i];
    for (int k = 0; k < nj; k++) {
      const double r2 = bc.square_minimum_image_distance(r[i], r[j[k]]);
      if (r2 < w*w) {
	LstItr<int> q;
	for (q.init(atom[i].neighbors); q.ok(); q.next())
	  if (q() == j[k])
	    break;
	if (are_atoms_in_same_molecule(i,j[k]))
	  continue;
	if (!q.ok())
	  Out() << sqrt(r2) << " A:\n" 
		<< i << " " << atom[i]
		<< j[k] << " " << atom[j[k]] << "\n";
      }
    }
  }
  Out() << flush;
}

void MSys::assign_bonds(LstInt2 bonds)
{
  for (int i = 0; i < atom.size(); i++) {
    atom[i].neighbors.remove_all();
    atom[i].bond_types.remove_all();
  }
  connect(bonds);
}

void MSys::connect(LstInt2 bonds)
{
  LstItrMod<Int2> b;
  HSet<Int2> seen;
  LstItr<int> n;
  for (int i = 0; i < atom.size(); i++)
    for (n.init(atom[i].neighbors); n.ok(); n.next())
      if (i < n())
	seen.add(Int2(i,n()));
  for (b.init(bonds); b.ok(); b.next()) {
    if (b().a > b().b) {
      const int tmp = b().b;
      b().b = b().a;
      b().a = tmp;
    }
    insist(0 <= b().a && b().a < atom.size());
    insist(0 <= b().b && b().b < atom.size());
    insist(b().a < b().b);
    if (seen.exists(b()))
      continue;
    seen.add(b());
    atom[b().a].neighbors.add(b().b);
    atom[b().b].neighbors.add(b().a);
  }
  find_molecules();
}

void MSys::find_neighbors()
{
  TIMESTART("MSys::find_neighbors");
  int i;
  const BoundaryConditions &bc = *boundary_conditions;
  const int n = atom.size();
  CVec r(n);
  DVec c(n);
  for (i = 0; i < n; i++) {
    r[i] = atom[i].position;
    c[i] = atom[i].covalent_radius();
    atom[i].neighbors.remove_all();
    atom[i].bond_types.remove_all();
  }
  NeighborList nlist(n, 6, 0);
  nlist.set_boundary_conditions(bc);
  nlist.coordinates_changed(r);
  for (i = 0; i < n; i++) {
    const int nj = nlist.number_of_neighbors()[i];
    const int *j = nlist.neighbors()[i];
    for (int k = 0; k < nj; k++)
      if (bc.square_minimum_image_distance(r[i], r[j[k]]) < sq(1.2*(c[i]+c[j[k]]))) {
	atom[i].neighbors.add(j[k]);
	atom[j[k]].neighbors.add(i);
      }
  }
  TIMESTOP("MSys::find_neighbors");
}

static void build_property_index(HTabHSetInt &atoms_with_property, const Vec<Atom> atom)
{
  HSetItr<Str> p;
  for (int i = 0; i < atom.size(); i++)
    for (p.init(atom[i].property); p.ok(); p.next())
      atoms_with_property[p()].add(i);
}

void MSys::match(Str pat)
{
  HTabHSetInt atoms_with_property;
  build_property_index(atoms_with_property, atom);
  AtomPattern(pat).run(atom,atoms_with_property);
}

void MSys::assign_properties_and_types()
{
  if (properties.size() == 0 && types.size() == 0)
    return;
  struct AssignProp : public AtomPattern
  {
    const VecStr &prop;
    HTabHSetInt &atoms_with_property;
    AssignProp(const char *p, const VecStr &prop_,
	       HTab<HSet<int> > &awp) : 
      AtomPattern(p), prop(prop_), atoms_with_property(awp)
    {
      if (ntest != prop.size())
	die("MSys::assign_property: number of properties (%d)\n"
	    "must match number of atom tests (%d)\nfor pattern '%s'", prop.size(), ntest, p);
    }
    void succeed(const int *path, Atom *at)
    {
      for (int i = 0; i < ntest; i++)
	if (strcmp(prop[i],"*")) {
	  at[path[i]].property.add(prop[i]);
	  atoms_with_property[prop[i]].add(path[i]);
	}
    }
  };
  struct AssignType : public AtomPattern
  {
    const VecStr &type;
    AssignType(const char *p, const VecStr &type_) : AtomPattern(p), type(type_) 
    {
      if (ntest != type.size())
	die("MSys::assign_type: number of types must match number of atom tests");
    }
    void succeed(const int *path, Atom *at)
    {
      for (int i = 0; i < ntest; i++)
	if (strcmp(type[i],"*"))
	  at[path[i]].type = type[i];
    }
  };
  HTabHSetInt atoms_with_property;
  int i;
  for (i = 0; i < atom.size(); i++) {
    atom[i].property.remove_all();
    atom[i].property.add(atom[i].symbol);
    atom[i].type = atom[i].symbol;
  }
  HSet<HSet<int> > rings;
  find_rings(rings, 10);
  HSetItr<HSet<int> > r;
  Vec<int> nring(atom.size(), 0);
  char tmp[256];
  for (r.init(rings); r.ok(); r.next()) {
    snprintf(tmp, 256, "ring%d", r().size());
    HSetItr<int> a;
    for (a.init(r()); a.ok(); a.next()) {
      atom[a()].property.add("ring");
      atom[a()].property.add(tmp);
      nring[a()]++;
    }
  }
  for (i = 0; i < atom.size(); i++)
    if (nring[i] == 1)
      atom[i].property.add("contained_in_exactly_1_ring");
    else if (nring[i] > 1)
      atom[i].property.add("contained_in_more_than_1_ring");
  build_property_index(atoms_with_property, atom);
  for (i = 0; i < properties.size(); i++) {
    AssignProp p(properties[i].str, properties[i].vec_str, atoms_with_property);
    p.run(atom,atoms_with_property);
  }
  for (i = 0; i < types.size(); i++) {
    AssignType t(types[i].str, types[i].vec_str);
    t.run(atom,atoms_with_property);
  }
}

void MSys::solvate(istream &stmp)
{
  TIMESTART("MSys::solvate");
  Solvate c;
  char bracket;
  stmp >> bracket;
  stmp.putback(bracket);
  if (bracket == '{')
    stmp >> c;
  map_to_central_box();
  rezero_center_of_mass();
  double diameter = 0;
  int i, imol;
  for (i = 0; i < atom.size(); i++) {
    const double r2 = atom[i].position.sq();
    if (r2 > diameter)
      diameter = r2;
  }
  diameter = 2 * c.box_factor * sqrt(diameter);
  Vec<Str> solvent_boxes;
  if (strlen(c.solvent_box) > 0) {
    solvent_boxes.resize(1);
    solvent_boxes[0] = c.solvent_box;
  } else {
    istream *s = StreamSearch(c.default_solvent_boxes);
    *s >> solvent_boxes;
    delete s;
  }
  const int nbox = solvent_boxes.size();
  if (nbox == 0)
    die("Msys::solvate: no solvent boxes found: specify solvent_box or create default_solvent_boxes file");
  int ibest = -1;
  double leftover = 1e10;
  for (i = 0; i < nbox; i++) {
    FILE *f = FileSearch(solvent_boxes[i]);
    char buf[256];
    BoundaryConditions *bc = 0;
    while (fgets(buf,256,f)) {
      if (buf[strlen(buf)-1] != '\n')
	die("MSys::solvate: error reading solvent box file");
      char *tmp = strstr(buf,"boundary_conditions");
      if (tmp) {
	istringstream s(tmp + strlen("boundary_conditions"));
	s >> bc;
	if (bc->type == non_periodic)
	  die("MSys::solvate: non-periodic solvent box: %s", (const char *) solvent_boxes[i]);
	break;
      }
    }
    fclose(f);
    if (!bc)
      die("MSys::solvate: could not find boundary_conditions in %s", (const char *) solvent_boxes[i]);
    if (c.boundary_conditions_type == non_periodic ||
	c.boundary_conditions_type == bc->type) {
      const double left = bc->min_diameter()*ceil(diameter/bc->min_diameter()) - diameter;
      if (left < leftover) {
	leftover = left;
	ibest = i;
      }
    }
    delete bc;
  }
  if (ibest == -1) {
    Out() << "Msys::solvate: no solvent box found with "
	  << c.boundary_conditions_type << " boundary conditions\n"
	  << flush;
    die("");
  }
  MSys solvent_box;
  solvent_box.create_from_xyz(solvent_boxes[ibest]);
  insist(solvent_box.boundary_conditions->type != non_periodic);
  solvent_box.tile((int) ceil(diameter/solvent_box.boundary_conditions->min_diameter()));
  insist(solvent_box.boundary_conditions->min_diameter() >= diameter);
  solvent_box.map_to_central_box();
  solvent_box.rezero_center_of_mass();
  const int nsolvent = solvent_box.molecule.size();
  const int nsolute = atom.size();
  const int n = nsolvent+nsolute;
  CVec r(n);
  DVec vdw(n, 0.0);
  for (imol = 0; imol < nsolvent; imol++) {
    r[imol] = solvent_box.molecule_center_of_mass(imol);
    for (i = solvent_box.molecule[imol].a; i < solvent_box.molecule[imol].b; i++) {
      const double v = solvent_box.atom[i].vdw_radius();
      if (v > vdw[imol])
	vdw[imol] = v;
    }
  }
  for (i = 0; i < nsolute; i++) {
    const Atom &a = atom[i];
    r[i+nsolvent] = a.position;
    vdw[i+nsolvent] = a.vdw_radius();
  }
  delete boundary_conditions;
  boundary_conditions = solvent_box.boundary_conditions->copy();
  NeighborList nlist(n, 5, 0);
  nlist.set_boundary_conditions(*boundary_conditions);
  nlist.coordinates_changed(r);
  Vec<int> included(nsolvent, 1);
  for (i = 0; i < n; i++) {
    const int nj = nlist.number_of_neighbors()[i];
    const int *ji = nlist.neighbors()[i];
    for (int k = 0; k < nj; k++) {
      const int j = ji[k];
      if (i < nsolvent && j < nsolvent)
	      continue;
      if (i >= nsolvent && j >= nsolvent)
	      continue;
      if (boundary_conditions->square_minimum_image_distance(r[i], r[j]) < sq(0.8*(vdw[i]+vdw[j]))) {
      	if (i < nsolvent)
      	  included[i] = 0;
      	else
      	  included[j] = 0;
      }
    }
  }
  int nsolvatom = 0, nsolvmol = 0;
  for (imol = 0; imol < nsolvent; imol++)
    if (included[imol]) {
      nsolvmol++;
      nsolvatom += solvent_box.molecule[imol].b - solvent_box.molecule[imol].a;
    }
  MSys mtmp;
  mtmp.properties = properties;
  mtmp.types = types;
  mtmp.atom.resize(nsolvatom);
  mtmp.molecule.resize(nsolvmol);
  int iat = 0, im = 0;
  for (imol = 0; imol < nsolvent; imol++)
    if (included[imol]) {
      const int iat0 = iat;
      const int i0 = solvent_box.molecule[imol].a;
      for (i = i0; i < solvent_box.molecule[imol].b; i++) {
	mtmp.atom[iat] = solvent_box.atom[i];
	LstItrMod<int> nei;
	for (nei.init(mtmp.atom[iat].neighbors); nei.ok(); nei.next())
	  nei() += iat0 - i0;
	iat++;
      }
      mtmp.molecule[im++] = Int2(iat0,iat);
    }
  insist(im == nsolvmol);
  insist(iat == nsolvatom);
  mtmp.assign_properties_and_types();
  append(mtmp);
  Out() << "Solvated with " << nsolvmol << " molecules ("
	<< nsolvatom << " atoms) from " 
	<< solvent_boxes[ibest] << "\n" << flush;
  TIMESTOP("MSys::solvate");
}

struct MSys_insert_t 
{ 
  const Atom *a; 
  int b; 
};

static int MSys_insert_cmp(const void *a, const void *b)
{
  const MSys_insert_t *atmp = (const MSys_insert_t *) a;
  const MSys_insert_t *btmp = (const MSys_insert_t *) b;
  if (atmp->b > btmp->b)
    return 1;
  else if (atmp->b < btmp->b)
    return -1;
  else
    return 0;
}

void MSys::insert(const Lst<Atom> &newatomlst, const Lst<int> &bonded_to, int *iold, int *inew)
{
  const int natom = atom.size();
  const int nnew = newatomlst.size();
  insist(bonded_to.size() == nnew);
  int i;
  if (nnew == 0) {
    for (i = 0; i < natom; i++)
      iold[i] = i;
    return;
  }
  MSys_insert_t *tmp = new MSys_insert_t[nnew];
  LstItr<Atom> a;
  LstItr<int> ib;
  int k = 0;
  for (a.init(newatomlst), ib.init(bonded_to); a.ok(); a.next(), ib.next(), k++) {
    tmp[k].a = &a();
    tmp[k].b = ib();
  }
  qsort((void *) tmp, nnew, sizeof(MSys_insert_t), MSys_insert_cmp);
  Vec<Atom> newatom(natom+nnew);
  for (i = k = 0; i < natom; i++) {
    newatom[i+k] = atom[i];
    iold[i] = i+k;
    while (k < nnew && i == tmp[k].b) {
      newatom[i+k+1] = *tmp[k].a;
      inew[k] = i+k+1;
      k++;
    }
  }
  insist(k == nnew);
  for (i = 0; i < natom+nnew; i++) {
    LstItrMod<int> ni;
    for (ni.init(newatom[i].neighbors); ni.ok(); ni.next())
      ni() = iold[ni()];
  }
  for (i = 0; i < nnew; i++) {
    const int p = inew[i];
    const int q = iold[tmp[i].b];
    newatom[p].neighbors.add(q);
    newatom[q].neighbors.add(p);
  }
  delete[] tmp;
  atom = newatom;
  find_molecules();
}

bool MSys::is_homogeneous() const
{
  if (molecule.size() < 2)
    return 1;
  const int a = molecule[0].a;
  const int b = molecule[0].b;
  insist(a == 0);
  for (int imol = 1; imol < molecule.size(); imol++) {
    if (molecule[imol].b - molecule[imol].a != b)
      return 0;
    int i, j;
    for (i = 0, j = molecule[imol].a; i < b; i++,j++) {
      if (strcmp(atom[i].symbol, atom[j].symbol))
	return 0;
      if (strcmp(atom[i].type, atom[j].type))
	return 0;
    }
  }
  return 1;
}

void MSys::write_solvation_shell_as_xyz(Str f) const
{
  const int n = atom.size();
  const int nmol = molecule.size();
  CVec r(n);
  copy_positions_to(r);
  NeighborList nlist(n,3.0,0); /* solvation shell: waters within 3.0 A */
  nlist.set_boundary_conditions(*boundary_conditions);
  nlist.coordinates_changed(r);
  Lst<const Atom *> a;
  Vec<int> imol(n);
  int i;
  for (i = 0; i < nmol; i++)
    for (int j = molecule[i].a; j < molecule[i].b; j++)
      imol[j] = i;
  for (i = 0; i < nmol; i++) {
    if (is_molecule_water(i)) {
      int j;
      for (j = molecule[i].a; j < molecule[i].b; j++) {
	const int nj = nlist.number_of_neighbors()[j];
	const int *jp = nlist.neighbors()[j];
	int k;
	for (k = 0; k < nj; k++)
	  if (!is_molecule_water(imol[jp[k]]))
	    break;
	if (k < nj)
	  break;
      }
      if (j == molecule[i].b)
	continue;
    }
    for (int j = molecule[i].a; j < molecule[i].b; j++)
      if (!atom[j].is_dummy())
	a.add(&atom[j]);
  }
  ostream *ss = FileStream(f);
  *ss << a.size() << "\nboundary_conditions " << boundary_conditions << "\n";
  LstItr<const Atom *> ai;
  for (ai.init(a); ai.ok(); ai.next())
    ai()->write_as_xyz(*ss);
  delete ss;
}


bool MSys::is_molecule_water(int imol) const
{
  const int a = molecule[imol].a;
  return molecule[imol].b - a == 3 &&
    ((!strcmp(atom[a].symbol,"O") && !strcmp(atom[a+1].symbol,"H") && !strcmp(atom[a+2].symbol,"H")) ||
     (!strcmp(atom[a].symbol,"H") && !strcmp(atom[a+1].symbol,"O") && !strcmp(atom[a+2].symbol,"H")) ||
     (!strcmp(atom[a].symbol,"H") && !strcmp(atom[a+1].symbol,"H") && !strcmp(atom[a+2].symbol,"O")));
}

bool MSys::has_same_topology_as(const MSys &m) const
{
  if (atom.size() != m.atom.size())
    return false;
  int i;
  LstItr<int> j, k;
  for (i = 0; i < atom.size(); i++) {
    if (atom[i].neighbors.size() != m.atom[i].neighbors.size())
      return false;
    for (j.init(atom[i].neighbors), k.init(m.atom[i].neighbors); j.ok(); j.next(), k.next())
      if (j() != k())
	return false;
  }
  return true;
}

static void color_atoms(const Atom *a, int i, int *imol, int im)
{
  if (imol[i] != -1)
    return;
  imol[i] = im;
  LstItr<int> n;
  for (n.init(a[i].neighbors); n.ok(); n.next())
    color_atoms(a,n(),imol,im);
}

void MSys::find_molecules()
{
  const int natom = atom.size();
  Vec<int> imol(natom, -1);
  int i, im = 0;
  for (i = 0; i < natom; i++)
    if (imol[i] == -1)
      color_atoms(atom, i, imol, im++);
  Vec<Lst<int> > a(im);
  for (i = 0; i < natom; i++)
    a[imol[i]].add(i);
  Lst<int> itot;
  for (i = 0; i < im; i++)
    itot.add(a[i]);
  insist(itot.size() == natom);
  Vec<int> index(natom);
  LstItr<int> j;
  int k;
  for (k = 0, j.init(itot); j.ok(); k++, j.next())
    index[j()] = k;
  Vec<Atom> atmp(natom);
  for (k = 0, j.init(itot); j.ok(); k++, j.next()) {
    atmp[k] = atom[j()];
    atmp[k].neighbors.remove_all();
    LstItr<int> n;
    for (n.init(atom[j()].neighbors); n.ok(); n.next())
      atmp[k].neighbors.add(index[n()]);
  }
  atom = atmp;
  molecule.resize(im);
  k = 0;
  for (i = 0; i < im; i++) {
    molecule[i].a = k;
    k += a[i].size();
    molecule[i].b = k;
  }
  insist(k == natom);
}

int MSys::number_of_water_molecules() const
{
  int n = 0;
  for (int i = 0; i < molecule.size(); i++)
    if (is_molecule_water(i))
      n++;
  return n;
}

int MSys::number_of_hydrogens() const
{
  int n = 0;
  for (int i = 0; i < atom.size(); i++)
    if (atom[i].is_hydrogen())
      n++;
  return n;
}

void MSys::init_acid()
{
  Lst<AcidicGroup> g;
  Lst<Atom> newatom;
  Lst<int> bonded_to;
  HTabHSetInt atoms_with_property;
  build_property_index(atoms_with_property, atom);
  int i;
  for (i = 0; i < acidic_group.size(); i++)
    acidic_group[i].find_groups(g,atom,newatom,bonded_to,atoms_with_property);
  Vec<int> iold(atom.size()), inew(newatom.size());
  insert(newatom, bonded_to, iold, inew);
  acid.copy_from_list(g);
  for (i = 0; i < acid.size(); i++)
    acid[i].reindex(iold,inew);
}

void MSys::change_protonation_state(int iacid, int state)
{
  if (iacid < 0 || iacid >= acid.size())
    die("MSys::change_protonation_state: acid index %d is out of range\n"
	"(only %d acidic groups)", iacid, acid.size());
  AcidicGroup &a = acid[iacid];
  a.state = state;
  for (int k = 0; k < a.index.size(); k++)
    atom[a.index[k]].copy_from(a.atom[state][k]);
  assign_properties_and_types();
} 

void MSys::set_velocities(LstCartesian c)
{
  const int n = atom.size();
  if (c.size() != n)
    die("MSys::set_velocities: %d velocities, %d atoms", 
	c.size(), n);
  LstItr<Cartesian> r;
  r.init(c);
  for (int i = 0; i < atom.size(); i++) {
    insist(r.ok());
    atom[i].velocity = r();
    r.next();
  }
  insist(!r.ok());
}

void MSys::set_positions(LstCartesian c)
{
  const int n = atom.size();
  if (c.size() != n)
    die("MSys::set_positions: %d positions, %d atoms", 
	c.size(), n);
  LstItr<Cartesian> r;
  r.init(c);
  for (int i = 0; i < atom.size(); i++) {
    insist(r.ok());
    atom[i].position = r();
    r.next();
  }
  insist(!r.ok());
}

void MSys::set_types(LstStr c)
{
  const int n = atom.size();
  if (c.size() != n)
    die("MSys::set_types: %d types, %d atoms", 
	c.size(), n);
  LstItr<Str> r;
  r.init(c);
  for (int i = 0; i < atom.size(); i++) {
    insist(r.ok());
    atom[i].type = r();
    r.next();
  }
  insist(!r.ok());
}

void MSys::assign_individual_types()
{
  for (int i = 0; i < atom.size(); i++) {
    char buf[256];
    snprintf(buf, 256, "%s%d",(const char *) atom[i].symbol, i);
    atom[i].type = buf;
  }
}

void MSys::align_with(MSys m)
{
  insist(m.atom.size() == atom.size());
  CVec c1(atom.size()), c2(atom.size());
  m.copy_positions_to(c1);
  copy_positions_to(c2);
  Out() << "RMSD: " << RootMeanSquareDistance(c1,c2) << " A\n";
  copy_positions_from(c2);
}

void MSys::set_mass(Str sym, double m)
{
  Atom::set_mass(sym,m);
}

void MSys::show_dihedral(Int4 i) const
{
  if (i.a < 0 || i.a >= atom.size())
    die("MSys::show_dihedral: first index out of range");
  if (i.b < 0 || i.b >= atom.size())
    die("MSys::show_dihedral: first index out of range");
  if (i.c < 0 || i.c >= atom.size())
    die("MSys::show_dihedral: first index out of range");
  if (i.d < 0 || i.d >= atom.size())
    die("MSys::show_dihedral: first index out of range");
  OutPrintf("Dihedral %d-%d-%d-%d (%s-%s-%s-%s): %f degrees\n",
	    i.a,i.b,i.c,i.d,
	    (const char *) atom[i.a].type,
	    (const char *) atom[i.b].type,
	    (const char *) atom[i.c].type,
	    (const char *) atom[i.d].type,
	    radians_to_degrees(Dihedral(atom[i.a].position,atom[i.b].position,
					atom[i.c].position,atom[i.d].position)));
}

void MSys::subset(const LstInt &i)
{
  MSys m;
  delete m.boundary_conditions;
  m.boundary_conditions = boundary_conditions->copy();
  m.atom.resize(i.size());
  Vec<int> newindex(atom.size(), -1);
  LstItr<int> j;
  int k;
  for (k = 0, j.init(i); j.ok(); k++, j.next()) {
    if (j() < 0 || j() >= atom.size())
      die("MSys::subset: index %d is out of range", j());
    newindex[j()] = k;
  }
  insist(k == i.size());
  for (k = 0, j.init(i); j.ok(); k++, j.next()) {
    if (j() < 0 || j() >= atom.size())
      die("MSys::subset: index %d is out of range", j());
    m.atom[k] = atom[j()];
    m.atom[k].neighbors.remove_all();
    m.atom[k].bond_types.remove_all();
    LstItr<int> n;
    LstItr<Str> bt;
    for (n.init(atom[j()].neighbors),bt.init(atom[j()].bond_types); n.ok(); n.next()) {
      const int nn = newindex[n()];
      insist(nn >= 0 && nn < i.size());
      m.atom[k].neighbors.add(nn);
      if (bt.ok()) {
	m.atom[k].bond_types.add(bt());
	bt.next();
      }
    }
  }
  *this = m;
  find_molecules();
  init_acid();
}

void MSys::remove_water()
{
  Lst<int> j;
  for (int i = 0; i < molecule.size(); i++)
    if (!is_molecule_water(i))
      for (int k = molecule[i].a; k < molecule[i].b; k++)
	j.add(k);
  subset(j);
}

void MSys::remove_water_and_ions()
{
  Lst<int> j;
  for (int i = 0; i < molecule.size(); i++)
    if (!is_molecule_water(i) && molecule[i].b - molecule[i].a > 1)
      for (int k = molecule[i].a; k < molecule[i].b; k++)
	j.add(k);
  subset(j);
}

void MSys::show_molalities() const
{
  Out() << "Molalities of ions:\n";
  const int nwater = number_of_water_molecules();
  if (nwater == 0)
    return;
  HTab<int> n;
  double water_mass = -1;
  for (int i = 0; i < molecule.size(); i++)
    if (molecule[i].b - molecule[i].a == 1)
      n[atom[molecule[i].a].symbol]++;
    else if (water_mass < 0 && is_molecule_water(i))
      water_mass = molecule_mass(i);
  insist(water_mass > 0);
  water_mass = mass_unit_to_g_mol(water_mass) * 1e-3;
  HTabIterator<int> ni;
  for (ni.init(n); ni.ok(); ni.next())
    Out() << ni.key() << ": " << (ni.val()/(nwater*water_mass)) << " molal\n";
  Out() << "\n\n" << flush;
}

void MSys::replace_waters_with_atoms(int n, Str sym, Str type)
{
  Lst<int> tmp;
  int i;
  for (i = 0; i < molecule.size(); i++)
    if (is_molecule_water(i))
      tmp.add(i);
  if (n > tmp.size())
    die("MSys::replace_waters_with_ions: not enough waters to change");
  Vec<int> waters;
  waters.copy_from_list(tmp);
  HSet<int> replace;
  while (replace.size() < n)
    replace.add(waters[(int) floor(UniformRandom() * waters.size())]);
  insist(replace.size() == n);
  tmp.remove_all();
  HSetItr<int> r;
  Lst<Cartesian> water_cm;
  for (r.init(replace); r.ok(); r.next())
    water_cm.add(molecule_center_of_mass(r()));
  for (i = 0; i < molecule.size(); i++)
    if (!replace.exists(i))
      for (int k = molecule[i].a; k < molecule[i].b; k++)
	tmp.add(k);
    else
      insist(is_molecule_water(i));
  subset(tmp);
  MSys m;
  m.atom.resize(water_cm.size());
  LstItr<Cartesian> wcm;
  for (i = 0, wcm.init(water_cm); wcm.ok(); i++, wcm.next()) {
    if (strlen(type) > 0 && strcmp(type,"*"))
      m.atom[i] = Atom(sym,type);
    else
      m.atom[i] = Atom(sym);
    m.atom[i].position = wcm();
  }
  m.find_molecules();
  append(m);
}

void MSys::more_types(istream &st)
{
  Lst<Str_VecStr> t;
  st >> t;
  types.add_from_list(t);
}

void MSys::more_properties(istream &st)
{
  Lst<Str_VecStr> p;
  st >> p;
  properties.add_from_list(p);
}

void MSys::reverse_velocities()
{
  for (int i = 0; i < atom.size(); i++)
    atom[i].velocity *= -1;
}

void MSys::set_coordinates_from_zmatrix(const ZMatrix z)
{
  z.check(atom.size());
  z.to_cartesian(atom);
}

void MSys::show_as_zmatrix(ZMatrix z) const
{
  z.check(atom.size());
  z.from_cartesian(atom);
  Out() << z << flush;
}

bool MSys::contains_only_dummy_atoms() const
{
  for (int i = 0; i < atom.size(); i++)
    if (!atom[i].is_dummy())
      return false;
  return true;
}

void MSys::create_from_mol2(Str s)
{
  const char *fname = s;
  FILE *f = FileSearch(fname);
  char buf[256];
  int nline = 0;
  while (fgets(buf, 256, f)) {
    nline++;
    if (!strncmp(buf, "@<TRIPOS>ATOM", 13))
      goto found_atom;
  }
  die("Mol::create_from_mol2: no @<TRIPOS>ATOM in '%s',", fname);
 found_atom:
  Lst<Atom> alst;
  while (fgets(buf, 256, f)) {
    nline++;
    if (!strncmp(buf, "@<TRIPOS>BOND", 13))
      goto found_bond;
    int n;
    char sym[256], t[256];
    double x, y, z;
    if (sscanf(buf, "%d %s %lf %lf %lf %s", &n, sym, &x, &y, &z, t) != 6)
      die("Mol::create_from_mol2: error reading coordinate line %s:%d\n%s", 
	  fname, nline, buf);
    if (islower(sym[1]))
      sym[2] = '\0';
    else
      sym[1] = '\0';
    Atom a(sym,t);
    a.position.x = x;
    a.position.y = y;
    a.position.z = z;
    alst.add(a);
  }
  die("Mol::create_from_mol2: no @<TRIPOS>BOND in '%s',", fname);
 found_bond:
  atom.copy_from_list(alst);
  while (fgets(buf, 256, f)) {
    nline++;
    if (is_all_white_space(buf) || 
	!strncmp(buf, "@<TRIPOS>SUBSTRUCTURE", 21))
      break;
    int n, i, j;
    char t[256];
    if (sscanf(buf, "%d %d %d %s", &n, &i, &j, t) != 4)
      die("Mol::create_from_mol2: error reading bond line %s:%d\n%s", 
	  fname, nline, buf);
    i--;
    j--;
    atom[i].neighbors.add(j);
    atom[j].neighbors.add(i);
    atom[i].bond_types.add(t);
    atom[j].bond_types.add(t);
  }
  while (fgets(buf, 256, f)) {
    nline++;
    if (!strncmp(buf, "@<TRIPOS>CRYSIN", 15)) {
      const bool ok = fgets(buf, 256, f);
      nline++;
      if (!ok)
	die("Mol::create_from_mol2: error reading CRYSIN line %s:%d\n%s",
	    fname, nline, buf);
      double a,b,c,alpha,beta,gamma;
      if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", 
		 &a, &b, &c, &alpha, &beta, &gamma) != 6)
	die("Mol::create_from_mol2: error reading CRYSIN line %s:%d\n%s", 
	    fname, nline, buf);
      delete boundary_conditions;
      boundary_conditions = BoundaryConditions::new_from_cell(a,b,c,alpha,beta,gamma);
      break;
    }
  }
  find_molecules();
  assign_properties_and_types();
  init_acid();
}

void MSys::add_kinetic_energy_to_atoms(double ke, Range r)
{
  if (r.size() < 2)
    die("MSys::add_kinetic_energy_to_atoms: range must contain at least 2 atoms");
  insist(0 <= r.lo && r.hi < atom.size());
  double m = 0;
  Cartesian p(0,0,0);
  int i;
  for (i = r.lo; i <= r.hi; i++) {
    m += atom[i].mass();
    p += atom[i].linear_momentum();
  }
  const Cartesian vcm = p/m;
  double old_ke = 0;
  for (i = r.lo; i <= r.hi; i++) {
    atom[i].velocity -= vcm;
    old_ke += atom[i].kinetic_energy();
  }
  if (is_almost_zero(old_ke))
    die("MSys::add_kinetic_energy_to_atoms: atoms must have some kinetic energy to start with");
  const double z = sqrt((ke+old_ke)/old_ke); 
  for (int i = r.lo; i <= r.hi; i++) {
    atom[i].velocity *= z;
    atom[i].velocity += vcm;
  }
}

void MSys::show_kinetic_energy() const
{
  Out() << "Kinetic energy: " << kinetic_energy() << " kcal/mol\n" << flush;
}

Complex MSys::nuclear_structure_factor(const Cartesian &k) const
{
  Complex s(0,0);
  for (int i = 0; i < atom.size(); i++)
    s += atom[i].nuclear_structure_factor(k);
  return s;
}

void MSys::show_nuclear_structure_factor(const Cartesian &k) const
{
  Out() << "Nuclear structure factor for " << k << ": " << nuclear_structure_factor(k) << "\n" << flush;
}

void MSys::check_neighbors() const
{
  bool err = false;
  for (int i = 0; i < atom.size(); i++) {
    LstItr<int> j;
    for (j.init(atom[i].neighbors); j.ok(); j.next()) {
      LstItr<int> k;
      for (k.init(atom[j()].neighbors); k.ok(); k.next())
	if (k() == i)
	  goto found;
      err = true;
      OutPrintf("*** MSys::check_neighbors: atom %d is a neighbor of atom %d, but not the reverse\n", i, j());
    found:
      continue;
    }
  }
  if (err)
    die("");
}

int MSys::molecule_containing_atom(int i) const
{
  for (int j = 0; j < molecule.size(); j++) 
    if (molecule[j].a <= i && i < molecule[j].b)
      return j;
  die("MSys::molecule_containing_atom: %d out of range", i);
  return -1; // should never get here
}

bool MSys::are_atoms_in_same_molecule(int i, int j) const
{
  return molecule_containing_atom(i) == molecule_containing_atom(j);
}

static double round_to_nearest_twelveth(double x)
{ 
  return round_to_nearest_integer(12.0*x)/12.0; 
}

void MSys::correct_symmetry(int nsymm)
{
  const int nat = atom.size();
  if (nat == 0)
    return;
  if (nsymm <= 0)
    die("MSys::correct_symmetry: nsymm = %d must be greater than 0", nsymm);
  if (nat % nsymm != 0)
    die("MSys::correct_symmetry: number of atoms must be a multiple of nsymm = %d", nsymm);
  if (boundary_conditions->type == non_periodic)
    die("MSys::correct_symmetry: boundary conditions must be periodic");
  const int n1 = nat/nsymm;
  int i;
  for (i = 0; i < nat; i++)
    if (strcmp(atom[i].symbol,atom[i%n1].symbol))
      die("MSys::correct_symmetry: atom symbols do not match for atoms %d (%s), %d (%s)",
	  i, (const char *) atom[i].symbol,
	  i%n1, (const char*) atom[i%n1].symbol);
    else if (strcmp(atom[i].type,atom[i%n1].type))
      die("MSys::correct_symmetry: atom types do not match for atoms %d (%s), %d (%s)",
	  i, (const char *) atom[i].type,
	  i%n1, (const char*) atom[i%n1].type);
  const Tensor a = boundary_conditions->lattice_vectors();
  const Tensor ainv = a.inverse();
  CVec r(nat); // crystallographic coordinates
  for (i = 0; i < nat; i++)
    r[i] = ainv * atom[i].position;
  Cartesian sum1(0,0,0);
  for (i = 0; i < n1; i++)
    sum1 += r[i];
  for (i = 1; i < nsymm; i++) {
    int j;
    Cartesian sum2(0,0,0);
    for (j = 0; j < n1; j++)
      sum2 += r[i*n1+j];
    Tensor u(0);
    for (j = 0; j < n1; j++)
      u += Tensor(r[j], r[i*n1+j]);
    u -= 1.0/double(n1) * Tensor(sum1,sum2);
    Tensor vt;
    Cartesian s;
    SingularValueDecomposition(u, s, vt);
    Tensor rotation = u*vt;
    Cartesian translation(0,0,0);
    for (j = 0; j < n1; j++)
      translation += r[i*n1+j] - rotation*r[j];
    translation /= n1;
    translation.apply(round_to_nearest_twelveth);
    rotation.apply(round_to_nearest_twelveth);
    for (j = 0; j < n1; j++)
      atom[i*n1+j].position = a*(rotation*r[j] + translation);
  }
}

static void ring_search(const int a, 
			const int i, 
			const int maxsize,
			int *path, 
			HSet<int> &seen, 
			Atom *atom, 
			HSet<HSet<int> > &rings)
{
  if (i > maxsize)
    return;
  if (i > 2 && a == path[0]) {
    HSet<int> r;
    for (int j = 0; j < i; j++)
      r.add(path[j]);
    rings.add(r);
    return;
  }
  if (seen.exists(a))
    return;
  path[i] = a;
  seen.add(a);
  LstItr<int> n;
  for (n.init(atom[a].neighbors); n.ok(); n.next())
    ring_search(n(), i+1, maxsize, path, seen, atom, rings);
  seen.remove(a);
}

void MSys::find_rings(HSet<HSet<int> > &rings, const int maxsize)
{
  rings.remove_all();
  HSet<int> seen;
  Vec<int> path(maxsize+1);
  for (int a = 0; a < atom.size(); a++)
    ring_search(a, 0, maxsize, path, seen, atom, rings);
}

void MSys::show_rings(int maxsize)
{
  HSet<HSet<int> > rings;
  find_rings(rings, maxsize);
  Out() << rings << "\n" << flush;
}

void MSys::show_inertia_tensor() const
{
  Out() << "Inertia tensor: " << inertia_tensor() << "\n";
}

void MSys::replace_with_fragments(istream &s)
{
  Vec<Fragment> f;
  s >> f;
  replace_with_fragments(f.size(), f);
}
  
void MSys::replace_with_fragments(int n, Fragment *f)
{
  struct FindFragment : public AtomPattern
  {
    Lst<Vec<int> > &p;
    FindFragment(const char *pattern, Lst<Vec<int> > &q) : AtomPattern(pattern), p(q) { }
    HSet<HSet<int> > seen;
    void succeed(const int *path, Atom *) 
    {
      Vec<int> pi(ntest);
      HSet<int> ps;
      for (int i = 0; i < ntest; i++) {
	pi[i] = path[i];
	ps.add(path[i]);
      }
      if (seen.exists(ps))
	return;
      p.add(pi);
      seen.add(ps);
    }
  };
  MSys m(*boundary_conditions);
  HTabHSetInt atoms_with_property;
  build_property_index(atoms_with_property, atom);
  for (int i = 0; i < n; i++) {
    Lst<Vec<int> > imsys, ifrag;
    FindFragment(f[i].pattern,imsys).run(atom, atoms_with_property);
    HTabHSetInt fragment_atoms_with_property;
    build_property_index(fragment_atoms_with_property, f[i].molsys.atom);
    FindFragment(f[i].pattern,ifrag).run(f[i].molsys.atom, fragment_atoms_with_property);
    if (ifrag.size() < 1)
      die("MSys::replace_with_fragments: pattern '%s' should match fragment at least once",
	  (const char *) f[i].pattern);
    const Vec<int> &iv = ifrag.first();
    CVec c2(iv.size());
    int k;
    for (k = 0; k < iv.size(); k++)
      c2[k] = f[i].molsys.atom[iv[k]].position;
    const Cartesian centroid2 = c2.average();
    LstItr<Vec<int> > im;
    for (im.init(imsys); im.ok(); im.next()) {
      CVec c1(im().size());
      for (k = 0; k < im().size(); k++)
	c1[k] = atom[im()[k]].position;
      const Cartesian centroid1 = c1.average();
      const Tensor rotation_matrix = RotationMatrixToAlign(c1,centroid1,c2,centroid2);
      MSys mtmp(f[i].molsys);
      for (k = 0; k < mtmp.atom.size(); k++) {
	Cartesian &r = mtmp.atom[k].position;
	r = rotation_matrix*(r - centroid2) + centroid1;
      }
      m.append(mtmp);
    }
  }
  *this = m;
}

void MSys::cap_methyl_groups()
{
  const double len = 1.1; // C-H bond length
  const double ang = 2*atan(sqrt(2)); // tetrahedral angle
  int a;
  Lst<Atom> h;
  Lst<int> bonded_to;
  for (a = 0; a < atom.size(); a++) {
    if (strcmp(atom[a].symbol, "C"))
      continue;
    const int nn = atom[a].neighbors.size();
    if (nn == 1) {
      const int b = atom[a].neighbors.first();
      if (atom[b].neighbors.size() < 2)
	continue;
      int c = -1;
      LstItr<int> ctmp;
      for (ctmp.init(atom[b].neighbors); ctmp.ok(); ctmp.next())
	if (!strcmp(atom[ctmp()].symbol, "H")) {
	  c = ctmp();
	  break;
	}
      if (c == -1)
	for (ctmp.init(atom[b].neighbors); ctmp.ok(); ctmp.next())
	  if (!strcmp(atom[ctmp()].symbol, "O")) {
	    c = ctmp();
	    break;
	  }
      if (c == -1)
	for (ctmp.init(atom[b].neighbors); ctmp.ok(); ctmp.next())
	  if (ctmp() != a) {
	    c = ctmp();
	    break;
	  }
      insist(c >= 0);
      Atom a1("H"), a2("H"), a3("H");
      if (!atom[a].pdb.is_empty()) {
	a1.pdb = atom[a].pdb;
	a2.pdb = atom[a].pdb;
	a3.pdb = atom[a].pdb;
	snprintf(a1.pdb.name, 3, " H");
	snprintf(a2.pdb.name, 3, " H");
	snprintf(a3.pdb.name, 3, " H");
      }
      a1.position = ZLocation(atom[a].position, atom[b].position, atom[c].position,
			      len, ang, M_PI/3);
      a2.position = ZLocation(atom[a].position, atom[b].position, atom[c].position,
			      len, ang, M_PI);
      a3.position = ZLocation(atom[a].position, atom[b].position, atom[c].position,
			      len, ang, -M_PI/3);
      h.add(a1);
      h.add(a2);
      h.add(a3);
      bonded_to.add(a);
      bonded_to.add(a);
      bonded_to.add(a);
    } else if (nn == 2) {
      Atom a1("H"), a2("H");
      if (!atom[a].pdb.is_empty()) {
	a1.pdb = atom[a].pdb;
	a2.pdb = atom[a].pdb;
	snprintf(a1.pdb.name, 3, " H");
	snprintf(a2.pdb.name, 3, " H");
      }
      const Cartesian b = atom[atom[a].neighbors.first()].position;
      const Cartesian c = atom[atom[a].neighbors.last()].position;
      a1.position = ZLocation(atom[a].position, b, c, len, ang, M_PI/2);
      a2.position = ZLocation(atom[a].position, b, c, len, ang, -M_PI/2);
      h.add(a1);
      h.add(a2);
      bonded_to.add(a);
      bonded_to.add(a);
    }
  }
  if (h.size() > 0) {
    Vec<int> iold(atom.size()), inew(h.size());
    insert(h, bonded_to, iold, inew);
  }
}

void MSys::rotate_methyl_group(int ic, double phi)
{
  phi = degrees_to_radians(phi);
  if (strcmp(atom[ic].symbol,"C"))
    die("MSys::rotate_methyl_group: atom %d is not a methyl carbon", ic);
  if (atom[ic].neighbors.size() != 4)
    die("MSys::rotate_methyl_group: atom %d does not have four neighbors", ic);
  int nh = 0, h[3], a = -1;
  LstItr<int> n;
  for (n.init(atom[ic].neighbors); n.ok(); n.next())
    if (atom[n()].is_hydrogen())
      h[nh++] = n();
    else
      a = n();
  if (nh != 3)
    die("MSys::rotate_methyl_group: atom %d does not have three hydrogen neighbors", ic);
  insist(a >= 0);
  const Cartesian rC = atom[ic].position;
  const Cartesian ra = atom[a].position;
  const Cartesian rb = 2*ra - atom[h[0]].position;
  for (int i = 0; i < nh; i++) {
    const Cartesian rh = atom[h[i]].position;
    atom[h[i]].position = ZLocation(rC,ra,rb,rh.distance(rC),Angle(rh,rC,ra),
				    Dihedral(rh,rC,ra,rb) + phi);
    
  }
}

void MSys::protonate_monovalent_oxygens_making_hydrogen_bonds_to(MSys m)
{
  const BoundaryConditions &bc = *boundary_conditions;
  const double len = 0.96; // OH bond length
  const double hbond_cutoff = 3.0; // cutoff for H bond
  Lst<Atom> h;
  Lst<int> bonded_to;
  for (int i = 0; i < atom.size(); i++) {
    const Atom &a = atom[i];
    if (strcmp(a.symbol,"O") || a.neighbors.size() != 1)
      continue;
    for (int j = 0; j < m.atom.size(); j++) {
      const Atom &ma = m.atom[j];
      if (strcmp(ma.symbol,"O") && strcmp(ma.symbol,"N"))
	continue;
      if (m.is_atom_bonded_to_hydrogen(j))
	continue;
      Cartesian r;
      bc.minimum_image_displacement(r, ma.position, a.position);
      if (r.sq() > sq(hbond_cutoff))
	continue;
      Atom newH("H");
      newH.position = a.position + len * r.as_unit_vector();
      if (!atom[i].pdb.is_empty()) {
	newH.pdb = a.pdb;
	snprintf(newH.pdb.name, 3, " H");
      }
      h.add(newH);
      bonded_to.add(i);
    }
  }
  if (h.size() > 0) {
    Vec<int> iold(atom.size()), inew(h.size());
    insert(h, bonded_to, iold, inew);
  }
}

void MSys::reorient_hydrogen_bonds_to(MSys m)
{
  const double hbond_cutoff = 3.0; // cutoff for H bond
  const BoundaryConditions &bc = *boundary_conditions;
  for (int i = 0; i < atom.size(); i++) {
    const Atom &a = atom[i];
    if (strcmp(a.symbol,"O"))
      continue;
    if (a.neighbors.size() != 2)
      continue;
    Atom *h = &atom[a.neighbors.first()];
    if (!h->is_hydrogen())
      h = &atom[a.neighbors.last()];
    if (!h->is_hydrogen())
      continue;
    for (int j = 0; j < m.atom.size(); j++) {
      const Atom &ma = m.atom[j];
      if (strcmp(ma.symbol,"O") && strcmp(ma.symbol,"N"))
	continue;
      Cartesian r;
      bc.minimum_image_displacement(r, ma.position, a.position);
      if (r.sq() > sq(hbond_cutoff))
	continue;
      LstItr<int> n;
      for (n.init(ma.neighbors); n.ok(); n.next())
	if (m.atom[n()].is_hydrogen() &&
	    bc.square_minimum_image_distance(a.position, 
					     m.atom[n()].position) < sq(2.2))
	  break;
      if (!n.ok())
	h->position = a.position + 
	  bc.minimum_image_distance(a.position,h->position) * r.as_unit_vector();
    }
  }
}

bool MSys::is_atom_bonded_to_hydrogen(int i) const
{
  LstItr<int> k;
  for (k.init(atom[i].neighbors); k.ok(); k.next())
    if (atom[k()].is_hydrogen())
      return true;
  return false;
}

