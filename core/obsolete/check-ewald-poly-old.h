/* header file for check-ewald-poly.c --
 * test code for polydisperse code in ewald.c
 * Copyright (C) 2007-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_CHECK_EWALD_POLY_OLD_H_
#define	_CHECK_EWALD_POLY_OLD_H_


/* check make_matrix_mob_ewald_3all() for mono and poly(a=1)
 * with SC config with N=1
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
check_make_matrix_mob_ewald_3all_poly_SC_1
(int version,
 double phi,
 double ewald_tr, double ewald_eps,
 int verbose, double tiny);

/* check make_matrix_mob_ewald_3all() for mono and poly(a=1)
 * with SC config with N=2
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  dir       : direction of the config, 0 (x), 1 (y), 2(z).
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
check_make_matrix_mob_ewald_3all_poly_SC_2
(int version,
 int dir,
 double phi,
 double ewald_tr, double ewald_eps,
 int verbose, double tiny);


#endif /* !_CHECK_EWALD_POLY_OLD_H_ */
