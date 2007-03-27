/* header file for bonds-guile.c --
 * guile interface for struct bonds
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-guile.h,v 1.1 2007/03/27 00:54:20 kichiki Exp $
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
 *     1.0 ; spring const
 *     2.1 ; natural distance
 *     ((0 1) ; list of pairs
 *      (1 2)
 *      (2 3)))
 *    (; bond 2
 *     1.0 ; spring const
 *     2.5 ; natural distance
 *     ((4 5) ; list of pairs
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
