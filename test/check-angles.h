/* header file for check-angles.c --
 * test code for angles.c and angles-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-angles.h,v 1.1 2008/04/12 19:16:56 kichiki Exp $
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
#ifndef	_CHECK_ANGLES_H_
#define	_CHECK_ANGLES_H_


/* check reading SCM script
 */
int
check_guile_get_angles (int verbose, double tiny);

/* check angle_calc_force()
 * INPUT
 *  k  : 
 *  t0 : (degree)
 */
int
check_angles_calc_force (double k, double t0,
			 int verbose, double tiny);


#endif /* !_CHECK_ANGLES_H_ */
