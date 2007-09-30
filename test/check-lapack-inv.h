/* header file for check-lapack-inv.c --
 * test code for lapack_inv_() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lapack-inv.h,v 1.2 2007/09/30 03:59:15 kichiki Exp $
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
#ifndef	_CHECK_LAPACK_INV_H_
#define	_CHECK_LAPACK_INV_H_


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_lapack_inv_ (int n, int verbose, double tiny);


#endif /* !_CHECK_LAPACK_INV_H_ */
