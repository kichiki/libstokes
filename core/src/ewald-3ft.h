/* header file for 'ewald-3ft.c' --
 * Solvers for 3 dimensional FT version problems
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft.h,v 4.7 2007/10/27 03:47:42 kichiki Exp $
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
#ifndef	_EWALD_3FT_H_
#define	_EWALD_3FT_H_


/** natural resistance problem **/
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t);

/* solve natural resistance problem in FT version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft_0 (struct stokes * sys,
		 const double *u, const double *o,
		 double *f, double *t);


/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       double *u, double *o);


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf);


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t);

/* solve natural resistance problem with lubrication in FT version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft_0 (struct stokes * sys,
		     const double *u, const double *o,
		     double *f, double *t);


/** mob_lub_3ft **/
/* solve natural mobility problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   double *u, double *o);


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf);


#endif /* !_EWALD_3FT_H_ */
