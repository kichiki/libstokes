/* header file for check.c --
 * utility routines for check
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check.h,v 1.2 2007/12/01 18:23:19 kichiki Exp $
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
#ifndef	_CHECK_H_
#define	_CHECK_H_



/* compare x and y
 */
int
compare (double x, double y, char *label,
	 int verbose, double tiny);

/* compare x and y and keep the max error
 */
int
compare_max (double x, double y, char *label,
	     int verbose, double tiny,
	     double *max);


#endif /* !_CHECK_H_ */
