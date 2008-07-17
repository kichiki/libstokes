/* header file for check-solve-cubic.c --
 * test code for solve_cubic() with GSL routine poly_solve_cubic()
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-solve-cubic.h,v 1.1 2008/07/17 03:08:23 kichiki Exp $
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
#ifndef	_CHECK_SOLVE_CUBIC_H_
#define	_CHECK_SOLVE_CUBIC_H_


int
check_solve_cubic (int verbose, double tiny);


#endif /* !_CHECK_SOLVE_CUBIC_H_ */
