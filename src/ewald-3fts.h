/* header file for 'ewald-3fts.c' --
 * Solvers for 3 dimensional FTS version problems
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts.h,v 5.10 2007/12/22 17:36:29 kichiki Exp $
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
#ifndef	_EWALD_3FTS_H_
#define	_EWALD_3FTS_H_


/** natural resistance problem **/
/* solve natural resistance problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts_0 (struct stokes * sys,
		  const double *u, const double *o, const double *e,
		  double *f, double *t, double *s);
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts (struct stokes * sys,
		const double *u, const double *o, const double *e,
		double *f, double *t, double *s);


/** natural mobility problem **/
/* solve natural mobility problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys  : system parameters
 *  iter : struct iter (if NULL is given, use sys->it for the solver)
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_3fts_0 (struct stokes *sys, struct iter *iter,
		  const double *f, const double *t, const double *e,
		  double *u, double *o, double *s);
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  s [np * 5] :
 */
void
solve_mob_3fts (struct stokes * sys,
		const double *f, const double *t, const double *e,
		double *u, double *o, double *s);


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  iter : struct iter (if NULL is given, use sys->it for the solver)
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : in the fluid-rest frame
 *  of [nf * 3] : in the fluid-rest frame
 *  ef [nf * 5] : in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts_0 (struct stokes *sys, struct iter *iter,
		  const double *f, const double *t, const double *e,
		  const double *uf, const double *of, const double *ef,
		  double *u, double *o, double *s,
		  double *ff, double *tf, double *sf);
/* solve natural mobility problem with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts (struct stokes * sys,
		const double *f, const double *t, const double *e,
		const double *uf, const double *of, const double *ef,
		double *u, double *o, double *s,
		double *ff, double *tf, double *sf);


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_lub_3fts_0 (struct stokes * sys,
		      const double *u, const double *o, const double *e,
		      double *f, double *t, double *s);
/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] : in the labo frame.
 *   o [np * 3] : in the labo frame.
 *   e [np * 5] : in the labo frame.
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_lub_3fts (struct stokes * sys,
		    const double *u, const double *o, const double *e,
		    double *f, double *t, double *s);


/** mob_lub_3fts **/
/* solve natural mobility problem with lubrication in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_lub_3fts_0 (struct stokes * sys,
		      const double *f, const double *t, const double *e,
		      double *u, double *o, double *s);
/* solve natural mobility problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] : in the labo frame.
 * OUTPUT
 *   u [np * 3] : in the labo frame.
 *   o [np * 3] : in the labo frame.
 *   s [np * 5] :
 */
void
solve_mob_lub_3fts (struct stokes * sys,
		    const double *f, const double *t, const double *e,
		    double *u, double *o, double *s);


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  of [nf * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  ef [nf * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_lub_3fts_0 (struct stokes * sys,
		      const double *f, const double *t, const double *e,
		      const double *uf, const double *of,
		      const double *ef,
		      double *u, double *o, double *s,
		      double *ff, double *tf, double *sf);
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_lub_3fts (struct stokes * sys,
		    const double *f, const double *t, const double *e,
		    const double *uf, const double *of,
		    const double *ef,
		    double *u, double *o, double *s,
		    double *ff, double *tf, double *sf);


#endif /* !_EWALD_3FTS_H_ */
