/* header file for 'ewald-3ft.c' --
 * Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft.h,v 4.1 2006/09/26 01:05:28 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 3D configuration
 * periodic boundary condition in 3 direction
 * FT version
 * non-dimension formulation
 */
#ifndef	_EWALD_3FT_H_
#define	_EWALD_3FT_H_


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FT version
 * INPUT
 *  n := np * 11
 *  x [n * 6] : FT
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n * 6] : UO
 */
void
atimes_ewald_3ft (int n, double *x, double *y, void * user_data);

/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_ewald_3ft (struct stokes * sys,
		    double *u, double *o,
		    double *f, double *t);

/** natural mobility problem **/
/* solve natural mobility problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
calc_mob_ewald_3ft (struct stokes * sys,
		    double *f, double *t,
		    double *u, double *o);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * under Ewald sum
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
calc_mob_fix_ewald_3ft (struct stokes * sys,
			double *f, double *t,
			double *uf, double *of,
			double *u, double *o,
			double *ff, double *tf);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_lub_ewald_3ft (struct stokes * sys,
			double *u, double *o,
			double *f, double *t);

/** natural mobility problem with lubrication with fixed particles **/
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
calc_mob_lub_fix_ewald_3ft (struct stokes * sys,
			    double *f, double *t,
			    double *uf, double *of,
			    double *u, double *o,
			    double *ff, double *tf);

#endif /* !_EWALD_3FT_H_ */