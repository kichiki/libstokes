/* header file for check-mob-fts.c --
 * test code for mob-fts problem
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-mob-fts.h,v 1.1 2007/12/22 18:24:00 kichiki Exp $
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
#ifndef	_CHECK_MOB_FTS_H_
#define	_CHECK_MOB_FTS_H_


/* compare solve_mob_3fts_matrix() with solve_mob_3fts_matrix_0()
 * and solve_mob_3fts_0().
 */
int
check_mob_fts (int verbose, double tiny);

/* compare solve_mob_lub_3fts_matrix()
 * with solve_mob_lub_3fts_matrix_0()
 * and solve_mob_lub_3fts_0().
 */
int
check_mob_lub_fts (int verbose, double tiny);


#endif /* !_CHECK_MOB_FTS_H_ */
