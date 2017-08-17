/* header file for check-bd.c --
 * test code for Brownian dynamics scheme
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd.h,v 1.1 2007/09/29 20:28:17 kichiki Exp $
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

#ifndef	_CHECK_BD_H_
#define	_CHECK_BD_H_


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_bd (int ncheb, int np, double r,
	  int verbose, double tiny);


#endif /* !_CHECK_BD_H_ */
