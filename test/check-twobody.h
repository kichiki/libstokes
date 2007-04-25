/* header file for check-twobody.c --
 * test code for twobody.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-twobody.h,v 1.1 2007/04/25 05:56:32 kichiki Exp $
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
#ifndef	_CHECK_TWOBODY_H_
#define	_CHECK_TWOBODY_H_


/* check twobody_lub()
 */
int
check_twobody_lub (int version, double r, double a1, double a2, int nmax,
		   int verbose, double tiny);


/* check scalars_minv_f_poly() with scalar_minv_f() for equal sphere
 */
int
check_twobody_scalars_res_with_equal (double r, int nmax,
				      int verbose, double tiny);


#endif /* !_CHECK_TWOBODY_H_ */
