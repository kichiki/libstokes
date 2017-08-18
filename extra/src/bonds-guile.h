/* header file for bonds-guile.c --
 * guile interface for struct BONDS
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-guile.h,v 1.11 2008/07/25 22:20:35 kichiki Exp $
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
 *     0       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      0      ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
 *      1.0    ;    p1   = A^{sp}, scaled spring constant
 *      2.1)   ;    p2   = L_{s} / length, scaled max extension
 *     ((0 1)  ; 3) list of pairs
 *      (1 2)
 *      (2 3))
 *      -1)    ; 4) number of exclusion for lubrication
 *             ;    negative means all particles in the chain is excluded.
 *    (; bond 2
 *     2       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      1      ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
 *      19.8   ;    p1 = N_{K,s}, the Kuhn steps for a spring
 *      106.0) ;    p2 = b_{K} [nm], the Kuhn length
 *             ;    note that, for dWLC (type == 6),
 *             ;    (p1, p2) = (k, r0 [nm]), where the potential is
 *             ;    (k/2) * (kT / r0^2) * (r-r0)^2
 *     ((4 5)  ; 3) list of pairs
 *      (5 6)
 *      (6 7))
 *      1)     ; 4) number of exclusion for lubrication
 *    (; bond 3
 *     7       ; 1) spring type (FENE-Fraenkel)
 *     (       ; 2) spring parameters (list with 4 elements)
 *      0      ;    fene = 0 means (p1, p2, p3) = (H, r0 [nm], tol)
 *      1.0e6  ;    p1 = H, the spring constant
 *      0.5    ;    p2 = r0 [nm], the natural length of the spring
 *      0.01)  ;    p3 = tol, the tolerance parameter "s"
 *             ;    note that, for FENE-Fraenkel (type == 7),
 *             ;    the scalar part of the force is
 *             ;    fr = H * (r/hat(r0) - 1.0) / (1 - ((1-r/hat(r0))/tol)^2)
 *             ;    where hat(r0) = r0 / L0 (L0 is given by "length" [nm])
 *     ((8 9)  ; 3) list of pairs
 *      (9 10))
 *      0)     ; 4) number of exclusion for lubrication
 *   ))
 * where spring types are
 *   0 : Hookean spring (Asp * (r - Ls))
 *   1 : wormlike chain (WLC)
 *   2 : inverse Langevin chain (ILC)
 *   3 : Cohen's Pade approximation
 *   4 : Warner spring
 *   5 : Hookean spring (Asp * r / Ls)
 *   6 : Hookean spring for dWLC
 *   7 : FENE-Fraenkel
 * OUTPUT
 *  returned value : struct BONDS
 *                   if NULL is returned, it failed (not defined)
 */
struct BONDS *
BONDS_guile_get (const char * var);


#endif /* !_BONDS_GUILE_H_ */
