/* header file for ewald-2ft.c --
 * Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2ft.h,v 1.1 2006/09/26 05:40:01 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 2D configuration
 * periodic boundary condition in 3 direction
 * FT version
 * non-dimension formulation
 */
#ifndef	_EWALD_2FT_H_
#define	_EWALD_2FT_H_


/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
calc_res_ewald_2ft (struct stokes * sys,
		    const double *u, const double *o,
		    double *f, double *t);

/** natural mobility problem **/
/* solve natural mobility problem in FT version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 */
void
calc_mob_ewald_2ft (struct stokes * sys,
		    const double *f, const double *t3,
		    double *u, double *o);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
calc_mob_fix_ewald_2ft (struct stokes * sys,
			const double *f, const double *t3,
			const double *uf, const double *of,
			double *u, double *o,
			double *ff, double *tf);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
calc_res_lub_ewald_2ft (struct stokes * sys,
			const double *u, const double *o,
			double *f, double *t);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_2ft (struct stokes * sys,
			    const double *f, const double *t3,
			    const double *uf, const double *of,
			    double *u, double *o,
			    double *ff, double *tf);

#endif /* !_EWALD_2FT_H_ */
