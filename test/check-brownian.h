/* header file for check-brownian.c --
 * test code for brownian.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-brownian.h,v 1.2 2007/11/04 00:17:31 kichiki Exp $
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
#ifndef	_CHECK_BROWNIAN_H_
#define	_CHECK_BROWNIAN_H_


int
check_cheb_minv (int n,
		 int verbose, double tiny);

int
check_cheb_lub (int n,
		int verbose, double tiny);


/* compare minv routines in the matrix and atimes versions
 *  BD_matrix_minv_FU() and BT_atimes_mob_FU()
 */
int
check_minv_FU (int verbose, double tiny);

/* compare lub routines in the matrix and atimes versions
 *  BD_matrix_lub_FU() and BT_atimes_lub_FU()
 */
int
check_lub_FU (int verbose, double tiny);

int
benchmark_BD_minv_FU_in_FTS (int np, int verbose, double tiny);

int
check_inv_by_submatrices (int n1, int n2, int verbose, double tiny);


#endif /* !_CHECK_BROWNIAN_H_ */
