/* header file foe dnaupd_c.c --
 * C wrappers for ARPACK's dnaupd_()
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dnaupd_c.h,v 1.1 2007/09/29 20:20:47 kichiki Exp $
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
#ifndef	_DNAUPD_C_H_
#define	_DNAUPD_C_H_


/* obtain min and max of real part eigenvalues by dnaupd_()
 * INPUT
 *  n : dimension of the matrix
 *  atimes (n, x, b, user_data) : routine to calc A.x and return b[]
 *  user_data : pointer to be passed to solver and atimes routines
 *  eps : required precision
 * OUTPUT
 *  l[2] : l[0] = min
 *         l[1] = max
 */
void dnaupd_wrap_min_max (int n, double *l,
			  void (*atimes)
			  (int, const double *, double *, void *),
			  void *user_data,
			  double eps);


#endif /* !_DNAUPD_C_H_ */
