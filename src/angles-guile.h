/* header file for angles-guile.c --
 * guile interface for struct angles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: angles-guile.h,v 1.1 2008/04/12 18:18:42 kichiki Exp $
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
 *     180.0   ; 2) angle in degree (theta_0)
 *     ((0 1 2); 3) list of triplets
 *      (1 2 3)
 *      (2 3 4)
 *     )
 *    )
 *    (; angle type 2
 *     20.0    ; 1) constant (k^{angle})
 *     90.0    ; 2) angle in degree (theta_0)
 *     ((3 4 5); 3) list of triplets
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
guile_get_angles (const char *var);


#endif /* !_ANGLES_GUILE_H_ */
