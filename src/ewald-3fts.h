/* header file for 'ewald-3fts.c' --
 * Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2002 Kengo Ichiki <ichiki@pegasus.me.jhu.edu>
 * $Id: ewald-3fts.h,v 4.1 2002/06/02 19:11:17 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 3D configuration
 * periodic boundary condition in 3 direction,
 * FTS version
 * non-dimension formulation
 */

void
atimes_ewald_3fts (int n, double *x, double *y, void * user_data);

void
calc_res_ewald_3fts (int np,
		     double *u, double *o, double *e,
		     double *f, double *t, double *s);
void
calc_mob_ewald_3fts (int n,
		     double *f, double *t, double *e,
		     double *u, double *o, double *s);
void
calc_mob_fix_ewald_3fts (int np, int nm,
			 double *f, double *t, double *e,
			 double *uf, double *of, double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf);
void
calc_res_lub_ewald_3fts (int np,
			 double *u, double *o, double *e,
			 double *f, double *t, double *s);
void
calc_mob_lub_fix_ewald_3fts (int np, int nm,
			     double *f, double *t, double *e,
			     double *uf, double *of, double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf);
