/* header file for noHI.c --
 * Solvers for no-hydrodynamic interaction problems
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: noHI.h,v 1.1 2008/04/29 03:21:26 kichiki Exp $
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
#ifndef	_NOHI_H_
#define	_NOHI_H_

#include "stokes.h" // struct stokes

/* solve natural mobility problem with fixed particles in F version
 * without HI
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  ff [nf * 3] :
 */
void
solve_mix_3f_noHI (struct stokes * sys,
		   const double *f,
		   const double *uf,
		   double *u,
		   double *ff);

/* solve natural mobility problem with fixed particles in FT version
 * without HI
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  o [nm * 3] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 */
void
solve_mix_3ft_noHI (struct stokes * sys,
		    const double *f, const double *t,
		    const double *uf, const double *of,
		    double *u, double *o,
		    double *ff, double *tf);

/* solve natural mobility problem with fixed particles in FTS version
 * without HI
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
solve_mix_3fts_noHI (struct stokes * sys,
		     const double *f, const double *t, const double *e,
		     const double *uf, const double *of, const double *ef,
		     double *u, double *o, double *s,
		     double *ff, double *tf, double *sf);


#endif /* !_NOHI_H_ */
