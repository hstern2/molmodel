/* $Id: geograd.h,v 1.6 2009/04/27 18:39:44 hstern Exp $ */

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

#ifndef GEOGRAD_H
#define GEOGRAD_H

class Cartesian;
class Tensor;

double Angle(const Cartesian &a, const Cartesian &b, const Cartesian &c);

double Dihedral(const Cartesian &a, const Cartesian &b, 
		const Cartesian &c, const Cartesian &d);

/* Distance from r0 to plane given by r1, r2, r3 */
double NormalDistance(const Cartesian &r0, const Cartesian &r1, 
		      const Cartesian &r2, const Cartesian &r3);

/***
 * Return theta and derivatives
 ***/
double AngleGradient(const Cartesian &r1, const Cartesian &r2, const Cartesian &r3,
		     Cartesian &g1, Cartesian &g2, Cartesian &g3);

/***
 * Return phi and derivatives
 ***/
double DihedralGradient(const Cartesian &r1, const Cartesian &r2,
			const Cartesian &r3, const Cartesian &r4,
			Cartesian &g1, Cartesian &g2, 
			Cartesian &g3, Cartesian &g4);

/* Distance from r0 to plane given by r1, r2, r3 and derivatives */
double NormalDistanceGradient(const Cartesian &r0, const Cartesian &r1, 
			      const Cartesian &r2, const Cartesian &r3,
			      Cartesian &g0, Cartesian &g1, Cartesian &g2, Cartesian &g3);

/***
 * Return z, where z is a distance r from a, 
 * the angle z-a-b is given by theta,
 * and the dihedral z-a-b-c is given by phi
 ***/
Cartesian ZLocation(const Cartesian &a, const Cartesian &b, const Cartesian &c,
		    double r, double theta, double phi);

#endif /* GEOGRAD_H */
