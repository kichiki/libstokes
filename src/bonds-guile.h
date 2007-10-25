/* header file for bonds-guile.c --
 * guile interface for struct bonds
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-guile.h,v 1.2 2007/10/25 05:55:20 kichiki Exp $
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
#ifndef	_BONDS_GUILE_H_
#define	_BONDS_GUILE_H_


/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     0         ; 1) spring type
 *     (         ; 2) spring parameters (list with 3 elements)
 *      0        ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
 *      1.0      ;    p1   = A^{sp}, scaled spring constant  (for fene == 0)
 *      2.1)     ;    p2   = L_{s} / a, scaled max extension (for fene == 0)
 *     ((0 1)    ; 3) list of pairs
 *      (1 2)
 *      (2 3)))
 *    (; bond 2
 *     2         ; 1) spring type
 *     (         ; 2) spring parameters (list with 3 elements)
 *      1        ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
 *      19.8     ;    p1 = N_{K,s}, the Kuhn steps for a spring (for fene = 1)
 *      106.0)   ;    p2 = b_{K} [nm], the Kuhn length          (for fene = 1)
 *     ((4 5)    ; 3) list of pairs
 *      (5 6)
 *      (6 7)))
 *   ))
 * OUTPUT
 *  returned value : struct bonds
 *                   if NULL is returned, it failed (not defined)
 */
struct bonds *
guile_get_bonds (const char * var);


#endif /* !_BONDS_GUILE_H_ */
