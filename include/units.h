/*
 * Copyright (c) 2008 Harry A. Stern
 * Copyright (c) 2008 University of Rochester
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */


/* $Id: units.h,v 1.27 2010/09/13 18:16:56 hstern Exp $ */

#ifndef UNITS_H
#define UNITS_H

#include "fns.h"

/***
 *
 * Base units:
 *
 * energy    kcal/mol
 * length    angstrom
 * time      psec
 *
 * All other internal units are derived from these.
 * For example, the unit of charge is sqrt(A kcal/mol)
 * and the unit of mass is (kcal/mol) (psec^2/A^2)
 *
 * The following constants are taken from NIST:
 * physics.nist.gov/cuu/Constants
 *
 ***/

#define BOHR_RADIUS 0.5291772108   /* 10^-10 m      */
#define AVOGADRO 6.0221415         /* 10^23 mol^-1  */
#define HARTREE 4.35974417         /* 10^-18 J */
#define CALORIE 4.184              /* J */
#define BOLTZMANN 1.3806505        /* 10^-23 J/K */
#define SPEED_OF_LIGHT 2.99792548  /* 10^8 m/s */
#define HBAR 1.05457168            /* 10^-34 J s */

inline double radians_to_degrees(double x)
{
  return x * (180/M_PI);
}

inline double degrees_to_radians(double x)
{
  return x * (M_PI/180);
}

inline double angstrom_to_bohr(double x)
{
  return x / BOHR_RADIUS;
}

inline double bohr_to_angstrom(double x)
{
  return x * BOHR_RADIUS;
}

inline double angstrom2_to_bohr2(double x)
{
  return x / sq(BOHR_RADIUS);
}

inline double bohr2_to_angstrom2(double x)
{
  return x * sq(BOHR_RADIUS);
}

inline double angstrom3_to_bohr3(double x)
{
  return x / cube(BOHR_RADIUS);
}

inline double bohr3_to_angstrom3(double x)
{
  return x * cube(BOHR_RADIUS);
}

inline double hartree_to_kcal_mol(double x)
{
  return x * (HARTREE*AVOGADRO*100/CALORIE);
}

inline double kcal_mol_to_hartree(double x)
{
  return x / (HARTREE*AVOGADRO*100/CALORIE);
}

inline double hartree_bohr_to_kcal_mol_A(double x)
{
  return x * ((HARTREE*AVOGADRO*100)/(BOHR_RADIUS*CALORIE));
}

inline double kcal_mol_A_to_hartree_bohr(double x)
{
  return x / ((HARTREE*AVOGADRO*100)/(BOHR_RADIUS*CALORIE));
}

inline double e_to_charge_unit(double x)
{
  return x * sqrt(BOHR_RADIUS*HARTREE*AVOGADRO*100/CALORIE);
}

inline double charge_unit_to_e(double x)
{
  return x / sqrt(BOHR_RADIUS*HARTREE*AVOGADRO*100/CALORIE);
}

inline double hartree_e_to_potential_unit(double x)
{
  return e_to_charge_unit(angstrom_to_bohr(x));
}

inline double potential_unit_to_hartree_e(double x)
{
  return bohr_to_angstrom(charge_unit_to_e(x));
}

inline double kcal_mol_e_to_potential_unit(double x)
{
  return hartree_e_to_potential_unit(kcal_mol_to_hartree(x));
}

inline double potential_unit_to_kcal_mol_e(double x)
{
  return hartree_to_kcal_mol(potential_unit_to_hartree_e(x));
}

inline double potential_unit_to_e_A(double x)
{
  return charge_unit_to_e(x);
}

inline double e_A_to_potential_unit(double x)
{
  return e_to_charge_unit(x);
}

inline double dipole_unit_to_debye(double x)
{
  return x * sqrt(CALORIE/(AVOGADRO*10));
}

inline double debye_to_dipole_unit(double x)
{
  return x / sqrt(CALORIE/(AVOGADRO*10));
}

inline double quadrupole_unit_to_au(double x)
{
  return charge_unit_to_e(angstrom2_to_bohr2(x));
}

inline double au_to_quadrupole_unit(double x)
{
  return bohr2_to_angstrom2(e_to_charge_unit(x));
}

inline double dipole_unit2_to_debye2(double x)
{
  return dipole_unit_to_debye(dipole_unit_to_debye(x));
}

inline double debye2_to_dipole_unit2(double x)
{
  return debye_to_dipole_unit(debye_to_dipole_unit(x));
}

inline double mass_unit_to_g_mol(double x)
{
  return x * (100*CALORIE);
}

inline double g_mol_to_mass_unit(double x)
{
  return x / (100*CALORIE);
}

inline double angstrom3_to_cm3_mol(double x)
{
  return x * (AVOGADRO*1e-1);
}

inline double cm3_mol_to_angstrom3(double x)
{
  return x / (AVOGADRO*1e-1);
}

inline double density_unit_to_g_cm3(double x)
{
  return mass_unit_to_g_mol(cm3_mol_to_angstrom3(x));
}

inline double g_cm3_to_density_unit(double x)
{
  return angstrom3_to_cm3_mol(g_mol_to_mass_unit(x));
}

inline double pressure_unit_to_bar(double x)
{
  return x * (CALORIE*1e5/AVOGADRO);
}

inline double bar_to_pressure_unit(double x)
{
  return x / (CALORIE*1e5/AVOGADRO);
}

inline double K_to_kcal_mol(double x)
{
  return x * (BOLTZMANN*AVOGADRO*1e-3/CALORIE);
}

inline double kcal_mol_to_K(double x)
{
  return x / (BOLTZMANN*AVOGADRO*1e-3/CALORIE);
}

inline double angstrom2_psec_to_m2_sec(double x)
{
  return x * 1e-8;
}

inline double m2_sec_to_angstrom2_psec(double x)
{
  return x / 1e-8;
}

inline double angular_velocity_unit_to_cm1(double x)
{
  return x / (2*M_PI*SPEED_OF_LIGHT*1e-2);
}

inline double cm1_to_angular_velocity_unit(double x)
{
  return x * (2*M_PI*SPEED_OF_LIGHT*1e-2);
}

inline double concentration_unit_to_molar(double x)
{
  return x * (1e4/AVOGADRO);
}

inline double molar_to_concentration_unit(double x)
{
  return x / (1e4/AVOGADRO);
}

inline double hbar_to_action_unit(double x)
{
  return x * (HBAR*AVOGADRO*1e-2/CALORIE);
}

inline double action_unit_to_hbar(double x)
{
  return x / (HBAR*AVOGADRO*1e-2/CALORIE);
}

inline double c_to_velocity_unit(double x)
{
  return x * (SPEED_OF_LIGHT*1e6);
}

inline double velocity_unit_to_c(double x)
{
  return x / (SPEED_OF_LIGHT*1e6);
}

inline double nm_to_kcal_mol(double x)
{
  return (0.1*2*M_PI*hbar_to_action_unit(1.0)*c_to_velocity_unit(1.0)) / x;
}

inline double kcal_mol_to_nm(double x)
{
  return nm_to_kcal_mol(x);
}

inline double A3_to_L_mol(double x)
{
  return x * AVOGADRO * 1e-4;
}

inline double L_mol_to_A3(double x)
{
  return x / AVOGADRO * 1e-4;
}

#endif /* UNITS_H */
