/* header file for check-bead-rod.c --
 * test code for bead-rod.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bead-rod.h,v 1.1 2008/06/13 03:12:48 kichiki Exp $
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
#ifndef	_CHECK_BEAD_ROD_H_
#define	_CHECK_BEAD_ROD_H_

#include <bead-rod.h>

struct BeadRod *
BeadRod_init_for_test (int nc,
		       double dt, double zeta,
		       double dr);


int
check_BeadRod_constraint_displacement (int n,
				       int verbose, double tiny);

int
check_BeadRod_solve_iter_gamma (int n, double eps,
				int verbose, double tiny);

int
check_BeadRod_solve_gamma_by_NITSOL (int n, double eps,
				     int verbose, double tiny);

#endif /* !_CHECK_BEAD_ROD_H_ */
