/* header file for 'ewald-3f.c' --
 * Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f.h,v 4.3 2006/10/05 19:23:59 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 3D configuration
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
#ifndef	_EWALD_3F_H_
#define	_EWALD_3F_H_


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *  f [np * 3] :
 */
void
calc_res_ewald_3f (struct stokes * sys,
		   const double *u,
		   double *f);

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 * OUTPUT
 *  u [np * 3] :
 */
void
calc_mob_ewald_3f (struct stokes * sys,
		   const double *f,
		   double *u);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_fix_ewald_3f (struct stokes * sys,
		       const double *f,
		       const double *uf,
		       double *u,
		       double *ff);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_lub_ewald_3f (struct stokes * sys,
		       const double *u,
		       double *f);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_3f (struct stokes * sys,
			   const double *f,
			   const double *uf,
			   double *u,
			   double *ff);

#endif /* !_EWALD_3F_H_ */
