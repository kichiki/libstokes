/* header file for check-grid.c --
 * test for grid.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-grid.h,v 1.2 2008/11/01 05:52:55 kichiki Exp $
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
#ifndef	_CHECK_GRID_H_
#define	_CHECK_GRID_H_


/*
 * INPUT
 */
int
check_GRID_ixyz_to_in_to_ixyz (int verbose, double tiny);


/*
 * INPUT
 *  np : number of particles
 *  flag_conf : 0 == random configuration
 *              1 == regular array configuration
 */
int
check_GRID_init_all_by_cutoff (int np,
			       int flag_conf,
			       int verbose, double tiny);


#endif /* !_CHECK_GRID_H_ */
