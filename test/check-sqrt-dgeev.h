/* header file for check-sqrt-dgeev.c --
 * test code for BD_sqrt_by_dgeev() in brownian.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-sqrt-dgeev.h,v 1.1 2007/11/01 04:57:51 kichiki Exp $
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
#ifndef	_CHECK_SQRT_DGEEV_H_
#define	_CHECK_SQRT_DGEEV_H_


int
check_BD_sqrt_by_dgeev (int n,
			int verbose, double tiny);


#endif /* !_CHECK_SQRT_DGEEV_H_ */
