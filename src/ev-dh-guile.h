/* header file for ev-dh-guile.c --
 * guile interface for struct EV_DH.
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-guile.h,v 1.5 2008/11/01 05:46:40 kichiki Exp $
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
#ifndef	_EV_DH_GUILE_H_
#define	_EV_DH_GUILE_H_


/* get ev_dh from SCM
 * in SCM, "ev-dh" are given by something like
 *  (define ev-dh '(
 *    ; system parameters
 *    1.0e-6   ; 1) epsilon for the cut-off distance of EV_DH interaction
 *    298.0    ; 2) temperature [K]
 *    80.0     ; 3) dielectric constant of the solution
 *    3.07     ; 4) Debye length [nm]
 *    1        ; 5) flag_grid (0 == particle-particle loop, 1 == grid loop)
 *    (        ; 6) list of DH types
 *     (; DH type 1
 *      2.43    ; 1) nu [e/nm]
 *      5.00    ; 2) l0 [nm]
 *      (0 1 2) ; 3) list of particles
 *     )
 *     (; DH type 2
 *      2.00    ; 1) nu [e/nm]
 *      4.00    ; 2) l0 [nm]
 *      (3 4)   ; 3) list of particles
 *     )
 *    )
 *  ))
 * INPUT
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  length : the characteristic length (dimensional, usually in nm)
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_DH
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_DH *
EV_DH_guile_get (const char *var,
		 double length, double peclet, int np);


#endif /* !_EV_DH_GUILE_H_ */
