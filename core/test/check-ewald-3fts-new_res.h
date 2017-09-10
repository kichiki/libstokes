/* header file for check-ewald-3fts-new_res.c --
 * test code for ewald-3fts-new.c
 * natural resistance problems
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
#ifndef	_CHECK_EWALD_3FTS_NEW_RES_H_
#define	_CHECK_EWALD_3FTS_NEW_RES_H_


/** check routines **/

// compare with new (mono) and new (poly a=1)
int
check_solve_res_3fts_b
(int verbose, double tiny);

// compare with mono and scaled poly (a = "scale")
int
check_solve_res_3fts_b2
(double scale,
 int verbose, double tiny);


// ** peroidic systems **

// compare with new (mono) and new (poly a=1)
int
check_solve_res_3fts_d
(int verbose, double tiny);

// compare with mono and scaled poly (a = "scale")
int
check_solve_res_3fts_d2
(double scale,
 int verbose, double tiny);


#endif /* !_CHECK_EWALD_3FTS_NEW_RES_H_ */
