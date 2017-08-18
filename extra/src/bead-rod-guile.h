/* header file for bead-rod-guile.c --
 * guile interface for struct BeadRod
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
#ifndef	_BEAD_ROD_GUILE_H_
#define	_BEAD_ROD_GUILE_H_

#include <libstokes-core.h> // struct stokes


/* get constraints from SCM and set struct BeadRod
 * "constraints" is given in SCM as 
 *  (define constraints '(
 *   ; system parameters
 *   1.0e-6    ; 1) tolerance
 *   "NITSOL"  ; 2) scheme for solving nonlinear equations
 *             ;    "linear" for iterative scheme in linear approximation
 *             ;    "NITSOL" for Newton-GMRES scheme by NITSOL library
 *   ; the following is for each constraint
 *   (         ; 3) constraint type 1
 *    5.0      ; 3-1) distance [nm]
 *    (        ; 3-2) list of particle-pairs
 *     (0 1)
 *     (1 2)
 *     (2 3)
 *   ))
 *   (         ; 4) constraint type 2
 *    10.0     ; 4-1) distance [nm]
 *    (        ; 4-2) list of particle-pairs
 *     (3 4)
 *     (4 5)
 *   ))
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "constraints".
 *  sys : struct stokes
 *  length : unit length in the simulation
 * OUTPUT
 *  returned value : struct BeadRod
 *                   if NULL is returned, it failed (not defined)
 */
struct BeadRod *
BeadRod_guile_get (const char *var,
		   struct stokes *sys,
		   double length);


#endif /* !_BEAD_ROD_GUILE_H_ */
