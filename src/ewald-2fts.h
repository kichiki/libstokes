/* header file for ewald-2fts.c --
 * Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2fts.h,v 1.1 2006/09/26 05:41:25 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 2D configuration
 * periodic boundary condition in 3 direction
 * FTS version
 * non-dimension formulation
 */
#ifndef	_EWALD_2FTS_H_
#define	_EWALD_2FTS_H_


/** natural resistance problem **/
/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_ewald_2fts (struct stokes * sys,
		     const double *u, const double *o, const double *e,
		     double *f, double *t, double *s);

/** natural mobility problem **/
/* solve natural mobility problem in FTS version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_2fts (struct stokes * sys,
		     const double *f, const double *t3, const double *e,
		     double *u, double *o, double *s);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   e [nm * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   ef [nf * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_fix_ewald_2fts (struct stokes * sys,
			 const double *f, const double *t3, const double *e,
			 const double *uf, const double *of, const double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_lub_ewald_2fts (struct stokes * sys,
			 const double *u, const double *o, const double *e,
			 double *f, double *t, double *s);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   e [nm * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   ef [nf * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_lub_fix_ewald_2fts (struct stokes * sys,
			     const double *f, const double *t3,
			     const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf);

#endif /* !_EWALD_2FTS_H_ */