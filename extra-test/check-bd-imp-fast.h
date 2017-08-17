/* header file for check-bd-imp-fast.c --
 * test code for bd-imp-fast.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd-imp-fast.h,v 1.1 2008/07/27 00:59:25 kichiki Exp $
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
#ifndef	_CHECK_BD_IMP_FAST_H_
#define	_CHECK_BD_IMP_FAST_H_


/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_rhs (int type, int np, int flag_noHI,
		  int verbose, double tiny);

/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_solve_cubic (int type, int np, int flag_noHI,
			  int verbose, double tiny);

/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_solve (int type, int np, int flag_noHI, double dt,
		    int verbose, double tiny);


#endif /* !_CHECK_BD_IMP_FAST_H_ */
