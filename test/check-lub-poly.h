/* header file for check-lub-poly.c --
 * test code for lubrication for polydisperse systems
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lub-poly.h,v 1.1 2007/04/14 00:37:16 kichiki Exp $
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
#ifndef	_CHECK_LUB_POLY_H_
#define	_CHECK_LUB_POLY_H_


/* check calc_lub_fts_2b_poly() with a1=a2=a
 * comparing with calc_lub_fts_2b()
 */
int
check_lub_fts_2b_poly (double r, int verbose, double tiny);

/* check matrix_lub_fts_2b_poly() with a1=a2=a
 * comparing with matrix_lub_fts_2b()
 */
int
check_matrix_lub_fts_2b_poly (double r, int verbose, double tiny);

/* compare calc_lub_fts_2b_poly() and matrix_lub_fts_2b_poly()
 */
int
check_atimes_matrix_lub_fts_2b_poly (double r, double a1, double a2,
				     int verbose, double tiny);


#endif /* !_CHECK_LUB_POLY_H_ */
