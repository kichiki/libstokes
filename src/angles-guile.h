/* header file for angles-guile.c --
 * guile interface for struct angles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: angles-guile.h,v 1.3 2008/05/24 05:46:24 kichiki Exp $
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
#ifndef	_ANGLES_GUILE_H_
#define	_ANGLES_GUILE_H_


#include "angles.h" // struct angles

/* get angles from SCM
 * in SCM, angles are given by something like
 *  (define angles '(
 *    (; angle type 1
 *     10.0    ; 1) constant (k^{angle})
 *     0.0     ; 2) angle in degree (theta_0)
 *     0       ; 3) scale flag (0 == scaled)
 *             ;    in this case, the above value for k is just used.
 *     ((0 1 2); 4) list of triplets
 *      (1 2 3)
 *      (2 3 4)
 *     )
 *    )
 *    (; angle type 2
 *     20.0    ; 1) constant (k^{angle})
 *     90.0    ; 2) angle in degree (theta_0)
 *     1       ; 3) scale flag (1 == not scaled yet)
 *             ;    in this case, the potential is given by 
 *             ;    (k/2) * kT * (theta - theta_0)^2
 *     ((3 4 5); 4) list of triplets
 *      (4 5 6)
 *     )
 *    )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "angles".
 * OUTPUT
 *  returned value : struct angles
 *                   if NULL is returned, it failed (not defined)
 */
struct angles *
angles_guile_get (const char *var);


#endif /* !_ANGLES_GUILE_H_ */
