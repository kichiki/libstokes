/* header file for ewald-3f-matrix.c --
 * Ewald summation technique with F version -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f-matrix.h,v 2.4 2006/10/05 04:52:26 ichiki Exp $
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
#ifndef	_EWALD_3F_MATRIX_H_
#define	_EWALD_3F_MATRIX_H_


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_ewald_3f_matrix (struct stokes * sys,
			  const double *u,
			  double *f);

/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *u,
			      double *f);

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_ewald_3f_matrix (struct stokes * sys,
			  const double *f,
			  double *u);

/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *f,
			      double *u);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
calc_mob_fix_ewald_3f_matrix (struct stokes * sys,
			      const double *f, const double *uf,
			      double *u, double *ff);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_3f_matrix (struct stokes * sys,
				  const double *f, const double *uf,
				  double *u, double *ff);

#endif /* !_EWALD_3F_MATRIX_H_ */
