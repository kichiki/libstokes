/* header file for excluded-volume-guile.c --
 * guile interface for struct EV
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume-guile.h,v 1.2 2008/05/24 05:42:50 kichiki Exp $
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


#include "bonds.h" // struct bonds
#include "excluded-volume.h" // struct EV

/* get ev-v from SCM and set struct EV
 * in SCM, ev-v is a list of parameter v [nm^3] or [micro m^3]
 * (depending on the dimension of the parameter "length")
 * for each spring:
 *  (define ev-v '(
 *   0.0012 ; for the spring 1
 *   0.002  ; for the spring 2
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-v".
 *  bonds : struct bonds
 *  length : unit length given by "length" in SCM (dimensional value)
 *  peclet : peclet number
 *  ev_r2  : square of max distance for EV interaction
 *  np     : number of particles (beads)
 * OUTPUT
 *  returned value : struct EV
 *                   if NULL is returned, it failed (not defined)
 */
struct EV *
EV_guile_get (const char *var,
	      const struct bonds *bonds,
	      double length, double peclet,
	      double ev_r2, int np);


#endif /* !_EXCLUDED_VOLUME_GUILE_H_ */
