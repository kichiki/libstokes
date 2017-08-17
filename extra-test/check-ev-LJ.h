/* header file for check-ev-LJ.c --
 * test code for ev-LJ.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-LJ.h,v 1.1 2008/05/24 06:06:33 kichiki Exp $
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
#ifndef	_CHECK_EV_LJ_H_
#define	_CHECK_EV_LJ_H_


/* check reading SCM script
 */
int
check_EV_LJ_calc_force (double r,
			int verbose, double tiny);


#endif /* !_CHECK_EV_LJ_H_ */
