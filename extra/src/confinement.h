/* header file for confine.c --
 * confinement interaction to the particles
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#ifndef	_CONFINE_H_
#define	_CONFINE_H_


#include <libstokes-core.h> // struct stokes


struct confinement {
  int type; /* 0 : sphere
	     *       R : radius of cavity centered at (0, 0, 0)
	     * 1 : sphere + hole
	     *       R : radius of cavity centered at (0, 0, 0)
	     *       r : radius of the hole at (0, 0, 1) direction
	     * 2 : cylinder
	     *       r       : radius of the cylinder
	     *       x, y, z : direction vector of the cylinder
	     *       the cylinder center goes through (0, 0, 0) and (x, y, z).
	     * 3 : dumbbell
	     *       R  : left cavity radius centered at (center1, 0, 0)
	     *       R2 : right cavity radius centered at (center2, 0, 0)
	     *       L  : length of the cylinder
	     *       r  : cylinder radius
	     *       the origin is at the center of the cylinder
	     * 4 : hex2d
	     *       R : cavity radius
	     *       r : cylinder radius
	     *       L : lattice spacing
	     * 5 == porous
	     *       R : particle radius
	     *       L : lattice spacing in x (2R for touching case)
	     * 6 == pillars
	     *       R : radius of cylinder (pillar)
	     *       L : lattice spacing
	     *       r : half width of the channel (in z direction)
	     *           negative value means no confinement in z direction
	     */
  // primary parameters
  double R;
  double r;
  double x, y, z;
  double R2;
  double L;

  /* LJ parameters, where the potential is given in DIMENSIONAL FORM as 
   *  U(r) = e ((r0/r)^12 - 2(r0/r)^6)
   * Note that r0 = 2^{1/6} sigma is the distance of potential minimum.
   * The repulsive force is evaluated by taking r = r0 at the contact point.
   * The parameters e and r0 should be dimensionless numbers defined by 
   *         e = hat(e) for flag_LJ == 0, or
   *         e = hat(e) / (peclet * hat(r0)) for flag_LJ == 1,
   *         r0 = hat(r0) = "r0" / length
   * where
   *         hat(e)  = "e" / kT
   *         hat(r0) = "r0" / length
   * and "e" and "r0" are dimensional numbers.
   * The dimensionless force (scalar part) is then given by
   *  hat(F) = 12 e (r0/r)^7 ((r0/r)^6 - 1)
   * where e is in the form of flag_LJ == 1.
   */
  int flag_LJ;
  double e;
  double r0;

  // derived parameters
  double theta; // angle of the cylinder opening from the sphere center
  double theta2; // angle of the cylinder opening on the right cavity

  // for dumbbell
  double center1; // x component of the center of left cavity
  double center2; // x component of the center of right cavity

  // for hex2d
  double Lx; // = 0.5 * L;
  double Ly; // = sqrt (3.0) * Lx;

  // for porous
  //         Lx = L
  //         Ly = L * sqrt (3.0)
  double Lz; // = L * sqrt (6.0) * 2.0 / 3.0;
};


/* initialize struct confinement
 * INPUT
 *  type : confinement type parameter (0 - 4)
 *         0 == sphere
 *               R : radius of cavity centered at (0, 0, 0)
 *         1 == sphere + hole
 *               R : radius of cavity centered at (0, 0, 0)
 *               r : radius of the hole at (0, 0, 1) direction
 *         2 == cylinder
 *               r       : radius of the cylinder
 *               x, y, z : direction vector of the cylinder
 *               the cylinder center goes through (0, 0, 0) and (x, y, z).
 *         3 == dumbbell
 *               R  : left cavity radius centered at (center1, 0, 0)
 *               R2 : right cavity radius centered at (center2, 0, 0)
 *               L  : length of the cylinder
 *               r  : cylinder radius
 *               the origin is at the center of the cylinder
 *         4 == hex2d
 *               R : cavity radius
 *               r : cylinder radius
 *               L : lattice spacing
 *         5 == porous
 *               R : particle radius
 *               L : lattice spacing in x (2R for touching case)
 *  R       : cavity radius
 *  r       : cylinder radius
 *  x, y, z : cylinder direction (for type 2, cylinder)
 *  R2      : right cavity radius (for type 3, dumbbell)
 *  L       : lattice spacing (for type 4, hex2d)
 *  Lennard-Jones parameters
 *  flag_LJ : if 0, e = hat(e) = "e" / kT
 *            if 1, e = hat(e) / (peclet * hat(r0))
 *  e  : the dimensionless number defined by
 *         e = hat(e) for flag_LJ == 0, or
 *         e = hat(e) / (peclet * hat(r0)) for flag_LJ == 1,
 *       where
 *         hat(e)  = "e" / kT
 *         hat(r0) = "r0" / length
 *       note that "e" and "r0" are dimensional numbers.
 *  r0 : the dimensionless number defined by
 *         r0 = hat(r0) = "r0" / length
 *       note that "r0" is a dimensional number.
 */
struct confinement *
CF_init (int type,
	 double R,
	 double r,
	 double x, double y, double z,
	 double R2,
	 double L,
	 int flag_LJ,
	 double e,
	 double r0);

void
CF_free (struct confinement *cf);


/* set LJ parameters for run
 * INPUT
 *  cf     : struct confinement
 *  peclet : peclet number
 * OUTPUT
 *  cf->e  : defined as 
 *             cf->e = hat(e)/(peclet * hat(r0))
 *           where
 *             hat(e)  = "e" / kT
 *             hat(r0) = "r0" / length
 *             "e" and "r0" are dimensional numbers
 *           therefore, cf->flag-LJ is set by 1.
 */
void
CF_set (struct confinement *cf,
	double peclet);

/* wrapper for confinement forces
 * INPUT
 *  c          : struct confinement
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_calc_force (struct confinement *cf,
	       struct stokes *sys,
	       double *f,
	       int flag_add);


#endif /* !_CONFINE_H_ */
