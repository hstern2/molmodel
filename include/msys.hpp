#ifndef MSYS_H
#define MSYS_H

#include "atom.hpp"
#include "bcond.hpp"
#include "boo.hpp"
#include "hset.hpp"
#include "acid.hpp"
#include "htab.hpp"
#include "zmatrix.hpp"
#include "range.hpp"
#include "ntuple.hpp"

struct StrInt 
{
  Str s;
  int i;
  friend ostream & operator<<(ostream &s, const StrInt &c) { return s << c.s << " " << c.i; }
  friend istream & operator>>(istream &s, StrInt &c) { return s >> c.s >> c.i; }
};

struct RotTrans
{
  Tensor rotation_matrix; // io
  Cartesian translation; // io
  classIO(RotTrans);
};

typedef Lst<StrInt> LstStrInt;
typedef HTab<HSet<int> > HTabHSetInt;
typedef Lst<Cartesian> LstCartesian;
typedef Lst<Int2> LstInt2;
typedef Lst<RotTrans> LstRotTrans;

struct PDBAtomRecord;

struct Solvate
{
  Str solvent_box, default_solvent_boxes; // in
  double box_factor; // in
  bcondtype_t boundary_conditions_type; // in
  
  Solvate() : 
    default_solvent_boxes("_msim_solvent_boxes"), 
    box_factor(1.25), 
    boundary_conditions_type(non_periodic) /* don't care */
  { }
  
  classIO(Solvate);
};

struct CreateFromPDB
{
  Str file; // in
  Boolean ignore_alternate_locations, ignore_hydrogens; // in
  HSet<Str> include_heterogen; // in
  CreateFromPDB() 
    : ignore_alternate_locations(1), 
      ignore_hydrogens(0)
  { }
  classIO(CreateFromPDB);
};

struct Str_VecStr 
{
  Str str;
  VecStr vec_str;
  Str_VecStr &operator=(const Str_VecStr &v)
  {
    str = v.str;
    vec_str = v.vec_str;
    return *this;
  }
  friend istream & operator>>(istream &s, Str_VecStr &c)
  { return s >> c.str >> c.vec_str; }
  friend ostream & operator<<(ostream &s, const Str_VecStr &c)
  { return s << c.str << "  " << c.vec_str; }
};

struct Fragment;

struct MSys
{

public:

  BoundaryConditions* boundary_conditions; // io
  Vec<Atom> atom; // out
  Vec<Int2> molecule; // out
  Vec<AcidicGroupDesc> acidic_group; // io
  Vec<Str_VecStr> properties, types; // io
  Vec<AcidicGroup> acid;

  MSys();
  MSys(const BoundaryConditions &bc);
  MSys(const MSys &);
  ~MSys();
  MSys & operator=(const MSys &);
  
  void show_self() const; // in
  void init_acid();
  void more_types(istream &); // in
  void more_properties(istream &); // in

  void assign_properties_and_types(); // in
  void match(Str pat); // in
  
  void replace_with_fragments(int n, Fragment *f);
  void replace_with_fragments(istream &); // in
  void replace_waters_with_atoms(int n, Str sym, Str type); // in
  void append(const MSys &m); // in
  void create_protein_from_sequence(Str seq); // in
  void create_from_topology(istream &); // in
  void create_from_topology(const LstStr &ls);
  void create_from_tinker(Str f); // in
  void create_from_molecule(const MSys &m, int imol);
  void create_from_molecules(const MSys &m, const Lst<int> &imol);
  void replicate(int n); // in
  void set_mass(Str sym, double m); // in
  void set_density(double d); // in
  void set_number_density(double d); // in
  void set_volume(double v); // in
  void set_molar_volume(double v); // in   L/mol
  void translate_atoms_by_random_lattice_vectors(); // in
  void translate(const Cartesian &r); // in
  void rezero_center_of_mass(); // in
  void rezero_linear_momentum(); // in
  void rezero_angular_momentum(); // in
  void scale(double d); // in
  void scale_centroids(double d); // in
  void translate_molecule(int i, const Cartesian &r); // in
  void translate_atom(int i, const Cartesian &r); // in
  void rotate(double phi, double theta, double psi); // in
  void rotate_molecule(int i, double phi, double theta, double psi); // in
  void apply_rotation_matrix(const Tensor r); // in
  void apply_symmetry_operations(LstRotTrans ops); // in

  void tile(int n); // in
  void cubic_lattice(int n, double d); // in
  void create_from_cif(Str f); // in
  void create_from_csd(Str f); // in
  void create_from_xyz(Str f); // in
  void create_from_mol2(Str f); // in
  void read_as_xyz(FILE *f);
  void read_as_xyz(Str fname); // in
  void write_as_xyz(Str f) const; // in
  void write_as_xyz(Str f, Str cmt) const; // in
  void write_as_xyz(ostream &) const;
  void write_as_xyz(ostream &, const char *cmt) const;
  void write_as_xyz(ostream &, double t) const;
  void write_as_mol2(Str f) const; // in
  void write_as_mol2(ostream &) const;
  void write_solvation_shell_as_xyz(Str f) const; // in
  void append_as_xyz(Str f) const; // in
  void create_from_pdb(CreateFromPDB p); // in
  void write_as_pdb(ostream &) const;
  void write_as_pdb(Str f) const; // in
  void show_center_of_mass() const; // in
  void show_clusters_as_xyz(int n, double w); // in
  void show_as_xyz() const; // in
  void show_geometry() const; // in
  void show_as_pdb() const; // in
  void show_properties() const; // in  
  void show_acid() const; // in
  void show_volume() const; // in
  void show_lattice_vectors() const; // in
  void show_reciprocal_lattice_vectors() const; // in
  void show_dihedral(Int4 i) const; // in
  void write_as_topology(ostream &) const;
  void write_as_topology(Str f) const; // in
  void show_as_topology() const; // in
  void show_neighbors() const; // in
  void map_to_central_box(); // in
  void randomize_positions(unsigned seed); // in
  void randomize_positions();
  void randomize_orientations(); // in
  void random_displacements(double r); // in
  void random_atomic_displacements(double r); // in
  void rotate_molecule(int i, const Tensor &);
  void copy_positions_from(const Cartesian *);
  void copy_positions_to(Cartesian *) const;
  void copy_velocities_from(const Cartesian *);
  void copy_velocities_to(Cartesian *) const;
  void copy_positions_from(const MSys &);
  void geometry_rmsd(Str xyzfile); // in
  double density() const;
  double number_density() const;
  Cartesian center_of_mass() const;
  Cartesian molecule_center_of_mass(int i) const;
  double molecule_mass(int i) const;
  void change_protonation_state(int iacid, int state); // in
  void arrange(); // in
  void check_types(); // in
  void check_distances(double r); // in
  void check_neighbors() const; // in
  void assign_bonds(LstInt2 bonds); // in
  void assign_individual_types(); // in
  void connect(LstInt2 bonds); // in
  void find_neighbors(); // in
  void find_molecules(); // in
  void subset(const LstInt &i);
  void solvate(istream &); // in
  void remove_water(); // in
  void remove_water_and_ions(); // in
  bool is_homogeneous() const;
  bool is_molecule_water(int) const;
  bool has_same_topology_as(const MSys &) const;
  bool contains_only_dummy_atoms() const;
  double kinetic_energy() const;
  void show_kinetic_energy() const; // in
  Cartesian linear_momentum() const;
  Cartesian angular_momentum() const;
  Tensor inertia_tensor() const;
  void show_inertia_tensor() const; // in
  double mass() const;
  void reverse_velocities();
  void integrate_velocities(const Cartesian *f, double dt);
  void integrate_positions(double dt);
  int number_of_hydrogens() const;
  int number_of_water_molecules() const;
  void show_molalities() const; // in
  void cap_methyl_groups(); // in
  void rotate_methyl_group(int ic, double phi); // in

  void protonate_monovalent_oxygens_making_hydrogen_bonds_to(MSys m); // in
  void reorient_hydrogen_bonds_to(MSys m); // in

  void show_mass() const; // in
  void set_types(LstStr c); // in
  void set_positions(LstCartesian c); // in
  void set_velocities(LstCartesian c); // in
  void show_distances_within(double w); // in
  void show_intermolecular_distances_within(double w); // in
  void align_with(MSys m); // in
  void align_principal_axes_with_xyz(); // in
  void set_coordinates_from_zmatrix(ZMatrix z); // in
  void show_as_zmatrix(ZMatrix z) const; // in
  void add_kinetic_energy_to_atoms(double ke, Range r); // in 
  Complex nuclear_structure_factor(const Cartesian &k) const;
  void show_nuclear_structure_factor(const Cartesian &k) const; // in
  int molecule_containing_atom(int i) const;
  bool are_atoms_in_same_molecule(int i, int j) const;
  void correct_symmetry(int nsymm); // in
  void find_rings(HSet<HSet<int> > &, const int maxsize);
  void show_rings(int maxsize); // in
  bool is_atom_bonded_to_hydrogen(int i) const;
  
  classIO(MSys);

private:

  void insert(const Lst<Atom> &newatom, const Lst<int> &bonded_to, int *iold, int *inew);
  
};

struct Fragment
{
  MSys molsys; // in
  Str pattern; // in
  classIO(Fragment);
};

#endif
