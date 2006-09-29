/* header file for ewald-3fts-matrix.h --
 * Ewald summation technique with FTS version -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts-matrix.h,v 2.2 2006/09/29 03:28:43 ichiki Exp $
 */
#ifndef	_EWALD_3FTS_MATRIX_H_
#define	_EWALD_3FTS_MATRIX_H_


/* make ewald-summed mobility matrix for FTS version
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
void
make_matrix_mob_ewald_3fts (struct stokes * sys, double * mat);

/* make lubrication matrix for FTS version for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
void
make_matrix_lub_ewald_3fts (struct stokes * sys,
			    double * mat);



/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FTS version
 * INPUT
 *  n := np * 11
 *  x [n * 11] : FTS
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3fts_matrix (int n, const double *x, double *y, void * user_data);

/** natural resistance problem **/
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
calc_res_ewald_3fts_matrix (struct stokes * sys,
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
calc_res_lub_ewald_3fts_matrix (struct stokes * sys,
				const double *u, const double *o,
				const double *e,
				double *f, double *t, double *s);

/** natural mobility problem **/
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
calc_mob_ewald_3fts_matrix (struct stokes * sys,
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
calc_mob_lub_ewald_3fts_matrix (struct stokes * sys,
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
calc_mob_fix_ewald_3fts_matrix (struct stokes * sys,
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
calc_mob_lub_fix_ewald_3fts_matrix (struct stokes * sys,
				    const double *f, const double *t,
				    const double *e,
				    const double *uf, const double *of,
				    const double *ef,
				    double *u, double *o, double *s,
				    double *ff, double *tf, double *sf);

#endif /* !_EWALD_3FTS_MATRIX_H_ */

