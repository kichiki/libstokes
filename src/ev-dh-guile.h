/* header file for ev-dh-guile.c --
 * guile interface for struct EV_DH.
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-guile.h,v 1.1 2008/04/26 04:34:07 kichiki Exp $
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
 * in SCM, angles are given by something like
 *  (define ev-dh '(
 *    ; system parameters
 *    4.0      ; 1) max distance for EV_DH interaction [nm]
 *    298.0    ; 2) temperature [K]
 *    80.0     ; 3) dielectric constant of the solution
 *    3.07     ; 4) Debye length [nm]
 *    (        ; 5) list of chain types
 *     (; chain type 1
 *      2.43    ; 1) nu [e/nm]
 *      5.00    ; 2) l0 [nm]
 *      (0 1 2) ; 3) list of particles
 *     )
 *     (; chain type 2
 *      2.00    ; 1) nu [e/nm]
 *      4.00    ; 2) l0 [nm]
 *      (3 4)   ; 3) list of particles
 *     )
 *    )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-dh".
 *  a   : the characteristic length [nm]
 *  pe  : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_DH
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_DH *
EV_DH_guile_get (const char *var,
		 double a, double pe, int np);


#endif /* !_EV_DH_GUILE_H_ */
