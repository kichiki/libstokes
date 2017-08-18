/* header file for confinement-guile.c --
 * guile interface for struct confinement
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: confinement-guile.h,v 1.1 2008/05/24 05:39:15 kichiki Exp $
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
#ifndef	_CONFINEMENT_GUILE_H_
#define	_CONFINEMENT_GUILE_H_


#include "confinement.h" // struct confinement

/* get confinement from SCM
 * in SCM, confinement is given by one of these:
 * for spherical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *  ))
 * for spherical confinement with a hole,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere+hole"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *    1.0  ;; radius of the hole at (0, 0, 1) direction
 *  ))
 * for cylindrical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "cylinder"    ;; the cylinder center goes through (0, 0, 0) and (x, y, z).
 *    10.0          ;; radius of the cylinder
 *    1.0  0.0  0.0 ;; direction vector (x, y, z) of the cylinder
 *  ))
 * for dumbbell confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "dumbbell" ;; the origin is at the center of the cylinder
 *    10.0       ;; left cavity radius centered at (center1, 0, 0)
 *    10.0       ;; right cavity radius centered at (center2, 0, 0)
 *    2.0        ;; length of the cylinder
 *    1.0        ;; cylinder radius
 *  ))
 * for 2D hexagonal confinement with cylinder pipe,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "hex2d"
 *    10.0    ;; cavity radius
 *    1.0     ;; cylinder radius
 *    12.0    ;; lattice spacing
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "confinement".
 * OUTPUT
 *  returned value : struct confinement
 *                   if NULL is returned, it failed (not defined)
 */
struct confinement *
CF_guile_get (const char *var);


#endif /* !_CONFINEMENT_GUILE_H_ */
