/* header file for check-ewald-shear.c --
 * test code for ewald.c in the new shear mode
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ewald-shear.h,v 1.1 2007/12/22 18:24:47 kichiki Exp $
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
#ifndef	_CHECK_EWALD_SHEAR_H_
#define	_CHECK_EWALD_SHEAR_H_


/* compare atimes_3all in periodic B.C. with shear-mode 1 or 2
 * based on the fact that for the simple-cubic lattice with 
 * the size (2L, L, L) containing 2 particles, the shift with 
 * (+/-)L in shear B.C. is the same without the shift.
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  shear_mode: 1 for (x = flow dir, y = grad dir)
 *              2 for (x = flow dir, z = grad dir)
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
check_atimes_3all_ewald_shear (int version,
			       int shear_mode,
			       double phi,
			       double ewald_tr, double ewald_eps,
			       int verbose, double tiny);


#endif /* !_CHECK_EWALD_SHEAR_H_ */
