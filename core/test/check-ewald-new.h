/* header file for check-ewald-new.c --
 * test code for polydisperse bugs, regarding to
 * non-ewald-new.c and ewald-new.c
 * Copyright (C) 2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_CHECK_EWALD_NEW_H_
#define	_CHECK_EWALD_NEW_H_


/** check routines **/

// compare with monodisperse systems
int
check_atimes_ewald_3all_new_0
(int version,
 int verbose, double tiny);

// compare with polydisperse systems with a = 1.0 (monodisperse)
int
check_atimes_ewald_3all_new_1
(int version,
 int verbose, double tiny);

// compare with polydisperse systems
// scaled by the factor "scale"
// F, T, S and U, O, E are properly scaled for comparison
// compare with polydisperse systems with a = 1.0 (monodisperse)
int
check_atimes_ewald_3all_new_2
(int version,
 double scale,
 double phi,
 int verbose, double tiny);


#endif /* !_CHECK_EWALD_NEW_H_ */
