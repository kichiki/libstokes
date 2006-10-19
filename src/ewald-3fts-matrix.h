/* header file for ewald-3fts-matrix.h --
 * Ewald summation technique with FTS version -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts-matrix.h,v 2.5 2006/10/19 04:15:49 ichiki Exp $
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
#ifndef	_EWALD_3FTS_MATRIX_H_
#define	_EWALD_3FTS_MATRIX_H_


/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_ewald_3fts_matrix (struct stokes * sys,
			     const double *u, const double *o, const double *e,
			     double *f, double *t, double *s);

/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_lub_ewald_3fts_matrix (struct stokes * sys,
				 const double *u, const double *o,
				 const double *e,
				 double *f, double *t, double *s);

/* solve natural mobility problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
solve_mob_ewald_3fts_matrix (struct stokes * sys,
			     const double *f, const double *t, const double *e,
			     double *u, double *o, double *s);

/* solve natural mobility problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
solve_mob_lub_ewald_3fts_matrix (struct stokes * sys,
				 const double *f, const double *t,
				 const double *e,
				 double *u, double *o, double *s);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_ewald_3fts_matrix (struct stokes * sys,
			     const double *f, const double *t,
			     const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_lub_ewald_3fts_matrix (struct stokes * sys,
				 const double *f, const double *t,
				 const double *e,
				 const double *uf, const double *of,
				 const double *ef,
				 double *u, double *o, double *s,
				 double *ff, double *tf, double *sf);

#endif /* !_EWALD_3FTS_MATRIX_H_ */

