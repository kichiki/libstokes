/* header file for excluded-volume-guile.c --
 * guile interface for struct EV
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume-guile.h,v 1.4 2008/07/25 22:18:12 kichiki Exp $
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
#ifndef	_EXCLUDED_VOLUME_GUILE_H_
#define	_EXCLUDED_VOLUME_GUILE_H_


#include "excluded-volume.h" // struct EV

/* get ev from SCM and set struct EV
 * in SCM, ev are given by something like
 *  (define ev '(
 *   5.0     ; max distance [nm] (or in the same dimension of "length")
 *   ( ; for EV type 1
 *    0.0012 ; v [nm^3] (or in the same dimension of "length")
 *    0      ; fene flag. if fene == 0, (p1,p2) = (A^{sp},L_{s})
 *    1.0    ; p1 = A^{sp}, scaled spring const
 *    2.1    ; p2 = L_{s} / length, scaled max extension
 *    (0 1 2); list of particles belongs to the EV parameters
 *   )
 *   ( ; for EV type 2
 *    0.002  ; v [nm^3] (or in the same dimension of "length")
 *    1      ; fene flag. if fene == 1, (p1,p2) = (N_{K,s},b_{K})
 *    19.8   ; p1 = N_{K,s}, the Kuhn steps for a spring
 *    106.0  ; p2 = b_{K} [nm], the Kuhn length
 *    (3 4)  ; list of particles belongs to the EV parameters
 *   )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-v".
 *  np     : number of particles (beads)
 *  length : unit length given by "length" in SCM (dimensional value)
 *  peclet : peclet number
 * OUTPUT
 *  returned value : struct EV
 *                   if NULL is returned, it failed (not defined)
 */
struct EV *
EV_guile_get (const char *var,
	      int np,
	      double length, double peclet);


#endif /* !_EXCLUDED_VOLUME_GUILE_H_ */
