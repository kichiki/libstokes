/* header file for ev-LJ-guile.c --
 * guile interface for struct EV_LJ.
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-LJ-guile.h,v 1.2 2008/06/13 03:11:39 kichiki Exp $
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
#ifndef	_EV_LJ_GUILE_H_
#define	_EV_LJ_GUILE_H_


/* get ev-LJ from SCM
 * in SCM, "ev-LJ" are given by something like
 *  (define ev-LJ '(
 *   (; LJ type 1
 *    10.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     0 1 2
 *    )
 *   )
 *   (; LJ type 2
 *    8.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    2.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     3 4
 *    )
 *   )
 *  ))
 * INPUT
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_LJ
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_LJ *
EV_LJ_guile_get (const char *var, int np);


#endif /* !_EV_LJ_GUILE_H_ */
