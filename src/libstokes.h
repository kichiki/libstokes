/* header file for library 'libstokes'
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libstokes.h,v 1.6 2006/10/05 19:25:21 ichiki Exp $
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
#ifndef	_LIBSTOKES_H_
#define	_LIBSTOKES_H_


/***********************************
 ** system parameters             **
 ***********************************/
struct stokes {
  int np; /* number of all particles */
  int nm; /* number of mobile particles */
  double * pos; /* position of particles */

  int version; /* 0 = F, 1 = FT, 2 = FTS */

  /* for ewald codes */
  double rmax;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double zeta,zeta2,zaspi,za2;
  double pivol;

  double lx, ly, lz;
  double llx[27], lly[27], llz[27]; /* for regist and lub */

  // self part
  double self_a;
  double self_c;
  double self_m;

  // table for lattice summation
  int flag_table; // 0 = inactive, 1 = active
  // real space
  int nr; // number of lattice points
  double * rlx;
  double * rly;
  double * rlz;
  // reciprocal space
  int nk; // number of lattice points
  double * ex;
  double * ey;
  double * ez;
  double * k;
  double * k1;
  double * k2;
  double * k3;
  double * ya;
  double * yb;
  double * yc;
  double * yg;
  double * yh;
  double * ym;

  /* for lubrication */
  double lubcut;

  /* for zeta program */
  double cpu1, cpu2, cpu3;

  /* for iterative solvers */
  struct iter * it;
};


/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void);

void
stokes_free (struct stokes * sys);

void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_ll (struct stokes * sys,
	       double lx, double ly, double lz);

void
stokes_set_zeta (struct stokes * sys,
		 double zeta, double cutlim);

double
zeta_by_tratio (struct stokes * sys,
		double tratio);

/** table for lattice summation **/
/* free ewald-table related memory from sys
 */
void
stokes_ewald_table_free (struct stokes * sys);

/* make ewald-summation table
 * INPUT
 *  (struct stokes *) sys :
 *  cutlim   :
 * OUTPUT
 *  (struct stokes *) sys : 
 */
void
stokes_ewald_table_make (struct stokes * sys, double cutlim);


/***********************************
 ** resistance problems           **
 ***********************************/

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

/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_ewald_3f_matrix (struct stokes * sys,
			  const double *u,
			  double *f);
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *u,
			      double *f);

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
		    const double *u, const double *o,
		    double *f, double *t);
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
			const double *u, const double *o,
			double *f, double *t);

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
calc_res_ewald_3ft_matrix (struct stokes * sys,
			   const double *u, const double *o,
			   double *f, double *t);
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
calc_res_lub_ewald_3ft_matrix (struct stokes * sys,
			       const double *u, const double *o,
			       double *f, double *t);

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
calc_res_ewald_3fts (struct stokes * sys,
		     const double *u, const double *o, const double *e,
		     double *f, double *t, double *s);
/* solve natural resistance problem with lubrication
 * in FTS version under Ewald sum
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
calc_res_lub_ewald_3fts (struct stokes * sys,
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

/***********************************
 ** resistance problems           **
 ** monolayer configurations (2D) **
 ***********************************/

/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_ewald_2f (struct stokes * sys,
		   const double *u,
		   double *f);
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_lub_ewald_2f (struct stokes * sys,
		       const double *u,
		       double *f);

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

/***********************************
 ** mobility problems             **
 ***********************************/

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

/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_ewald_3f_matrix (struct stokes * sys,
			  const double *f,
			  double *u);
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
calc_mob_fix_ewald_3f_matrix (struct stokes * sys,
			      const double *f, const double *uf,
			      double *u, double *ff);

/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *f,
			      double *u);
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_3f_matrix (struct stokes * sys,
				  const double *f, const double *uf,
				  double *u, double *ff);


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
		    const double *f, const double *t,
		    double *u, double *o);
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
			const double *f, const double *t,
			const double *uf, const double *of,
			double *u, double *o,
			double *ff, double *tf);
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
			    const double *f, const double *t,
			    const double *uf, const double *of,
			    double *u, double *o,
			    double *ff, double *tf);

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
calc_mob_ewald_3ft_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   double *u, double *o);
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
calc_mob_fix_ewald_3ft_matrix (struct stokes * sys,
			       const double *f, const double *t,
			       const double *uf, const double *of,
			       double *u, double *o,
			       double *ff, double *tf);
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
calc_mob_lub_ewald_3ft_matrix (struct stokes * sys,
			       const double *f, const double *t,
			       double *u, double *o);
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
calc_mob_lub_fix_ewald_3ft_matrix (struct stokes * sys,
				   const double *f, const double *t,
				   const double *uf, const double *of,
				   double *u, double *o,
				   double *ff, double *tf);

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
calc_mob_ewald_3fts (struct stokes * sys,
		     const double *f, const double *t, const double *e,
		     double *u, double *o, double *s);
/* solve natural mobility problem with fixed particles in FTS version
 * under Ewald sum
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
calc_mob_fix_ewald_3fts (struct stokes * sys,
			 const double *f, const double *t, const double *e,
			 const double *uf, const double *of, const double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf);
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
calc_mob_lub_fix_ewald_3fts (struct stokes * sys,
			     const double *f, const double *t, const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf);

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

/***********************************
 ** mobility problems             **
 ** monolayer configurations (2D) **
 ***********************************/

/* solve natural mobility problem in F version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
calc_mob_ewald_2f (struct stokes * sys,
		   const double *f,
		   double *u);
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
calc_mob_fix_ewald_2f (struct stokes * sys,
		       const double *f,
		       const double *uf,
		       double *u,
		       double *ff);
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
calc_mob_lub_fix_ewald_2f (struct stokes * sys,
			   const double *f,
			   const double *uf,
			   double *u,
			   double *ff);

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


/************************************************
 ** atimes forms (unnatural mobility problems) **
 ************************************************/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all (int n, const double *x, double *y, void * user_data);

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix (int n, const double *x,
			  double *y, void * user_data);

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_table (int n, const double *x, double *y, void * user_data);

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_table_matrix (int n, const double *x,
				double *y, void * user_data);

#endif /* !_LIBSTOKES_H_ */
