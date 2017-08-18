/* header file for check-solve-gen-linear --
 * test code for solve_gen_linear() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-solve-gen-linear.h,v 1.2 2007/09/30 04:00:50 kichiki Exp $
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
#ifndef	_CHECK_SOLVE_GEN_LINEAR_H_
#define	_CHECK_SOLVE_GEN_LINEAR_H_


/* check for some local routines
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_split_merge (int n1, int n2,
		   int verbose, double tiny);


/* check for local test implementation inverse_by_sub()
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_inverse_by_sub (int n1, int n2,
		      int verbose, double tiny);


/* 
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_solve_gen_linear (int n1, int n2,
			int verbose, double tiny);


#endif /* !_CHECK_SOLVE_GEN_LINEAR_H_ */
