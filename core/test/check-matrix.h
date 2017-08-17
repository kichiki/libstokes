/* header file for check-matrix.c --
 * test code for matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-matrix.h,v 1.1 2007/11/03 22:55:16 kichiki Exp $
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
#ifndef	_CHECK_MATRIX_H_
#define	_CHECK_MATRIX_H_


int
benchmark_dgemm (int n, int verbose, double tiny);

int
benchmark_mul_matrices (int n, int verbose, double tiny);


#endif /* !_CHECK_MATRIX_H_ */
