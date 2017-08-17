/* header file foe check-ewald.c --
 * test code for ewald.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ewald.h,v 1.1 2007/04/20 02:00:00 kichiki Exp $
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
#ifndef	_CHECK_EWALD_H_
#define	_CHECK_EWALD_H_


/** check routines **/

/* compare atimes and matrix processes for ewald_3all
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  phi       : volume fraction, that is, phi = (4/3)pi a^3/l^3
 *  ewald_tr  :
 *  ewald_eps :
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_ewald_3all_atimes_matrix_SC (int version,
				   double phi,
				   double ewald_tr, double ewald_eps,
				   int verbose, double tiny);


#endif /* !_CHECK_EWALD_H_ */
