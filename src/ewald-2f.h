/* header file for ewald-2f.c --
 * Ewald summation technique under 2D
 * this is a wrapper package for ewald-3f.c
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2f.h,v 1.3 2006/10/19 04:16:32 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 2D configuration
 * periodic boundary condition in 3 direction
 * F version
 * non-dimension formulation
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
#ifndef	_EWALD_2F_H_
#define	_EWALD_2F_H_


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_ewald_2f (struct stokes * sys,
		    const double *u,
		    double *f);

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
solve_mob_ewald_2f (struct stokes * sys,
		    const double *f,
		    double *u);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_ewald_2f (struct stokes * sys,
		    const double *f,
		    const double *uf,
		    double *u,
		    double *ff);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_lub_ewald_2f (struct stokes * sys,
			const double *u,
			double *f);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_lub_ewald_2f (struct stokes * sys,
			const double *f,
			const double *uf,
			double *u,
			double *ff);

#endif /* !_EWALD_2F_H_ */
