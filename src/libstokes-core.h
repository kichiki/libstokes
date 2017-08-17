/* header file for library 'libstokes-core'
 * Copyright (C) 1993-2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_LIBSTOKES_CORE_H_
#define	_LIBSTOKES_CORE_H_

#include <stdio.h> // FILE


/***********************************
 ** system parameters             **
 ** from stokes.h                 **
 ***********************************/
struct stokes {
  int version; /* 0 = F, 1 = FT, 2 = FTS  */

  int np;      /* number of all particles  */
  int nm;      /* number of mobile particles  */
  double *pos; /* position of particles  */

  /**
   * parameters for the polydisperse system
   */
  double *a;   /* radius of particles
		* Note : NULL (default) is for monodisperse system  */
  double rmin; /* minimum distance for HI calculation 
		* as in Rotne-Prager for r < 2a.
		* if (r < rmin*(ai+aj)), r = rmin*(ai+aj)
		* (default is 0.0) */
  int twobody_nmax;// max order for the coefficient for twobody_scalars_res()
  int twobody_lub; // 0 (far form) or 1 (lub form) for twobody_scalars_res()
  int *poly_table; // for (i,j), [i*np+j] gives the index of "twobody_f_list"
  struct twobody_f_list *twobody_f_list;

  /**
   * slip parameters
   */
  /* Note: if slip == NULL, the system is treated as the no-slip */
  double *slip;     // slip length
  double *slip_a;   // effective radius for the laplacian terms
  double *slip_G32; // Lambda(3,2) = 1/Lambda(2,3) for a-self.
  double *slip_G30; // Lambda(3,0) = 1/Lambda(0,3) for c-self.
  double *slip_G52; // Lambda(5,2) = 1/Lambda(2,5) for m-self.
  // slip table -- twobody_nmax and twobody_lub below are used, too
  int *slip_table; /* for (i,j), [i*np+j] gives the index
		    * of "twobody_slip_f_list" */
  struct twobody_slip_f_list *twobody_slip_f_list;

  /**
   * imposed flow
   */
  double Ui[3];
  double Oi[3];
  double Ei[5];

  /* auxiliary imposed-flow parameters for simple shear */
  int shear_mode; /* 0 : imposed flow is given by Ui,Oi,Ei.
		   * 1 : x = flow dir, y = grad dir
		   * 2 : x = flow dir, z = grad dir
		   */
  double shear_rate;
  double shear_shift; /* shift of the periodic cell (H.rate.t) at the time t
		       * which is in [0, L)
		       * H : cell size in the grad dir
		       * L : cell size in the flow dir
		       */

  /**
   * periodic parameters
   */
  int periodic; // 0 = non periodic, 1 = periodic

  /* for ewald codes */
  double ewald_eps;
  double rmax2;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double xi, xi2, xiaspi, xia2;
  double pivol;

  double lx, ly, lz;
  int ilx[27], ily[27], ilz[27];
  // note: ll[xyz][i] = l[xyz] * (double)il[xyz][i]

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
  int * rmx;
  int * rmy;
  int * rmz;
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

  /**
   * for lubrication
   */
  double lubmin2;/* square of min distance for lub
		  * (distance is replaced by its root-square)
		  */
  double lubmax; /* max distance for lub;
		  * for the pair beyond this is just ignored.
		  * 0 means no limit for open systems and 
		  * all particles within +/-1 cells in x,y,z for the periodic
		  */
  /* exclusion list for lubrication due to the bonding */
  //struct list_ex *ex_lub;

  /**
   * for zeta program
   */
  double cpu1, cpu2, cpu3;

  /**
   * for iterative solvers
   */
  struct iter * it;
};


/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void);

void
stokes_free (struct stokes * sys);

/* set np and nm and allocate the memory for pos[np*3]
 * also struct list_ex *ex_lub is allocated here
 * (becuase np is necessary for ex_lub).
 */
void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_Ui (struct stokes * sys,
	       double uix, double uiy, double uiz);
void
stokes_set_Oi (struct stokes * sys,
	       double oix, double oiy, double oiz);
void
stokes_set_Ei (struct stokes * sys,
	       double eixx, double eixy, double eixz,
	       double eiyz, double eiyy);

/* set auxiliary imposed flow parameters for simple shear
 * and overwrite the imposed parameters Ui, Oi, Ei.
 * INPUT
 *  shear_mode : 1 (x = flow dir, y = grad dir)
 *               2 (x = flow dir, z = grad dir)
 *  shear_rate :
 * OUTPUT
 *  sys : struct stokes
 */
void
stokes_set_shear (struct stokes *sys,
		  int shear_mode,
		  double shear_rate);
/* get shear_shift at the given time t
 * INPUT
 *  sys : struct stokes
 *  t   : current time 
 *  t0  : time of the reference
 *  s0  : shift at t0
 * OUTPUT
 *  returned value : shift in the range [-lx/2, lx/2)
 */
double
stokes_get_shear_shift (struct stokes *sys,
			double t,
			double t0, double s0);
/* set shear_shift at the given time t
 * INPUT
 *  sys : struct stokes
 *  t   : current time 
 *  t0  : time of the reference
 *  s0  : shift at t0
 * OUTPUT
 *  sys->shear_shift : set and place in the range [-lx/2, lx/2)
 */
void
stokes_set_shear_shift (struct stokes *sys,
			double t,
			double t0, double s0);

void
stokes_set_l (struct stokes * sys,
	      double lx, double ly, double lz);

void
stokes_set_xi (struct stokes * sys,
	       double xi, double ewald_eps);

double
xi_by_tratio (struct stokes * sys,
	      double tratio);


/* set iter param
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *   out   : FILE * to output debug info.
 */
void
stokes_set_iter (struct stokes *sys,
		 const char *solver,
		 int max,
		 int restart,
		 double eps,
		 int debug,
		 FILE *out);

/* set pos for all particles safely by another array
 * INPUT
 *  pos[np*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < np*3)
 */
void
stokes_set_pos (struct stokes *sys,
		const double *pos);

/* set pos for mobile particles safely by another array
 * INPUT
 *  pos[nm*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < nm*3)
 */
void
stokes_set_pos_mobile (struct stokes *sys,
		       const double *pos);

/* set pos for fixed particles safely by another array
 * INPUT
 *  pos[nf*3] : only fixed particles are set, where nf = np - nm
 * OUTPUT
 *  sys->pos[i] for (i = nm*3; i < np*3)
 */
void
stokes_set_pos_fixed (struct stokes *sys,
		      const double *pos);

/* set radius (sys->a[], sys->twobody_f_list, and sys->poly_table).
 * Note that the default setting (sys->a == NULL) is for monodisperse system
 * where a=1 for all particles
 * INPUT
 *  a[np] :
 *  sys->twobody_nmax : define sys->twobody_nmax before calling.
 * OUTPUT
 *  sys->a[np]             :
 *  sys->poly_table[np*np] :
 *  sys->twobody_f_list[]  :
 */
void
stokes_set_radius (struct stokes *sys,
		   const double *a);
/* unset radius (sys->a[], sys->twobody_f_list, and sys->poly_table).
 * that is, the system is treated as a monodisperse system
 * where a=1 for all particles as in the default setting.
 * INPUT
 *  sys                    : struct stokes
 * OUTPUT
 *  sys->a[np]             : freed and set NULL
 *  sys->poly_table[np*np] : freed and set NULL
 *  sys->twobody_f_list[]  : freed and set NULL
 */
void
stokes_unset_radius (struct stokes *sys);


/** slip parameters **/
/* set slip parameters (slip[], slip_a[], slip_G32[], slip_G30[], slip_G52[])
 * Note that the default setting (sys->slip == NULL) is for no-slip system
 * where gamma=0 for all particles
 * INPUT
 *  gamma[np] : slip length
 * OUTPUT
 *  sys->slip[np]     : slip length
 *  sys->slip_a[np]   : effective radius for the laplacian terms
 *  sys->slip_G32[np] : Lambda(3,2) = 1/Lambda(2,3) for a-self.
 *  sys->slip_G30[np] : Lambda(3,0) = 1/Lambda(0,3) for c-self.
 *  sys->slip_G52[np] : Lambda(5,2) = 1/Lambda(2,5) for m-self.
 */
void
stokes_set_slip (struct stokes *sys,
		 const double *gamma);

/* unset slip params (slip[], slip_a[], slip_G32[], slip_G30[], slip_G52[])
 * that is, the system is treated as a no-slip system
 * where gamma=0 for all particles as in the default setting.
 * INPUT
 *  sys                    : struct stokes
 * OUTPUT
 *  sys->slip[np]     : freed and set NULL
 *  sys->slip_a[np]   : freed and set NULL
 *  sys->slip_G32[np] : freed and set NULL
 *  sys->slip_G30[np] : freed and set NULL
 *  sys->slip_G52[np] : freed and set NULL
 */
void
stokes_unset_slip (struct stokes *sys);


/* make a copy of struct stokes s0
 */
struct stokes *
stokes_copy (struct stokes *s0);


/************************************
 ** miscellaneous routines         **
 ************************************/
/* from bench.h */
/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
long
ptime_ms (void);

/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
double
ptime_ms_d (void);

/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
long
ptime_micros (void);


/***********************************
 ** resistance problems           **
 ***********************************/

/** natural resistance problem **/
/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f (struct stokes * sys,
	      const double *u,
	      double *f);

/* solve natural resistance problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_lub_3f (struct stokes * sys,
		  const double *u,
		  double *f);

/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_3f_matrix (struct stokes * sys,
		     const double *u,
		     double *f);

/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_lub_3f_matrix (struct stokes * sys,
			 const double *u,
			 double *f);

/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_3ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t);

/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t);

/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_3ft_matrix (struct stokes * sys,
		      const double *u, const double *o,
		      double *f, double *t);

/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft_matrix (struct stokes * sys,
			  const double *u, const double *o,
			  double *f, double *t);

/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_3fts (struct stokes * sys,
		const double *u, const double *o, const double *e,
		double *f, double *t, double *s);

/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_lub_3fts (struct stokes * sys,
		    const double *u, const double *o, const double *e,
		    double *f, double *t, double *s);

/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_3fts_matrix (struct stokes * sys,
		       const double *u, const double *o, const double *e,
		       double *f, double *t, double *s);

/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_lub_3fts_matrix (struct stokes * sys,
			   const double *u, const double *o,
			   const double *e,
			   double *f, double *t, double *s);

/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_2f (struct stokes * sys,
	      const double *u,
	      double *f);

/* solve natural resistance problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_lub_2f (struct stokes * sys,
		  const double *u,
		  double *f);

/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_2ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t);

/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_lub_2ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t);

/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_2fts (struct stokes * sys,
		const double *u, const double *o, const double *e,
		double *f, double *t, double *s);

/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_lub_2fts (struct stokes * sys,
		    const double *u, const double *o, const double *e,
		    double *f, double *t, double *s);


/***********************************
 ** mobility problems             **
 ***********************************/

/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 * OUTPUT
 *  u [np * 3] :
 */
void
solve_mob_3f (struct stokes * sys,
	      const double *f,
	      double *u);

/* solve natural mobility problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_lub_3f (struct stokes * sys,
		  const double *f,
		  double *u);

/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_3f_matrix (struct stokes * sys,
		     const double *f,
		     double *u);

/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_lub_3f_matrix (struct stokes * sys,
			 const double *f,
			 double *u);

/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       double *u, double *o);

/* solve natural mobility problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   double *u, double *o);

/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft_matrix (struct stokes * sys,
		      const double *f, const double *t,
		      double *u, double *o);

/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft_matrix (struct stokes * sys,
			  const double *f, const double *t,
			  double *u, double *o);

/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_3fts (struct stokes * sys,
		const double *f, const double *t, const double *e,
		double *u, double *o, double *s);

/* solve natural mobility problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_lub_3fts (struct stokes * sys,
		    const double *f, const double *t, const double *e,
		    double *u, double *o, double *s);

/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_3fts_matrix (struct stokes * sys,
		       const double *f, const double *t, const double *e,
		       double *u, double *o, double *s);

/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_lub_3fts_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   const double *e,
			   double *u, double *o, double *s);

/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
solve_mob_2f (struct stokes * sys,
	      const double *f,
	      double *u);

/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 */
void
solve_mob_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       double *u, double *o);

/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
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
solve_mob_2fts (struct stokes * sys,
		const double *f, const double *t3, const double *e,
		double *u, double *o, double *s);


/***********************************
 ** mixed problems                **
 ***********************************/

/* solve natural mobility problem with fixed particles in F version
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_3f (struct stokes * sys,
	      const double *f,
	      const double *uf,
	      double *u,
	      double *ff);

/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_lub_3f (struct stokes * sys,
		  const double *f,
		  const double *uf,
		  double *u,
		  double *ff);

/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3f_matrix (struct stokes * sys,
		     const double *f, const double *uf,
		     double *u, double *ff);

/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_lub_3f_matrix (struct stokes * sys,
			 const double *f, const double *uf,
			 double *u, double *ff);

/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3ft_matrix (struct stokes * sys,
		      const double *f, const double *t,
		      const double *uf, const double *of,
		      double *u, double *o,
		      double *ff, double *tf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3ft_matrix (struct stokes * sys,
			  const double *f, const double *t,
			  const double *uf, const double *of,
			  double *u, double *o,
			  double *ff, double *tf);

/* solve natural mobility problem with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3fts (struct stokes * sys,
		const double *f, const double *t, const double *e,
		const double *uf, const double *of, const double *ef,
		double *u, double *o, double *s,
		double *ff, double *tf, double *sf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3fts (struct stokes * sys,
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
solve_mix_3fts_matrix (struct stokes * sys,
		       const double *f, const double *t,
		       const double *e,
		       const double *uf, const double *of,
		       const double *ef,
		       double *u, double *o, double *s,
		       double *ff, double *tf, double *sf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3fts_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   const double *e,
			   const double *uf, const double *of,
			   const double *ef,
			   double *u, double *o, double *s,
			   double *ff, double *tf, double *sf);

/* solve natural mobility problem with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_2f (struct stokes * sys,
	      const double *f,
	      const double *uf,
	      double *u,
	      double *ff);

/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_lub_2f (struct stokes * sys,
		  const double *f,
		  const double *uf,
		  double *u,
		  double *ff);

/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_2ft (struct stokes * sys,
		   const double *f, const double *t3,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf);

/* solve natural mobility problem with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_2fts (struct stokes * sys,
		const double *f, const double *t3, const double *e,
		const double *uf, const double *of, const double *ef,
		double *u, double *o, double *s,
		double *ff, double *tf, double *sf);

/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_2fts (struct stokes * sys,
		    const double *f, const double *t3,
		    const double *e,
		    const double *uf, const double *of,
		    const double *ef,
		    double *u, double *o, double *s,
		    double *ff, double *tf, double *sf);


/** no-hydrodynamic interaction problems **/
/* from noHI.h */

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
solve_mix_3f_noHI (struct stokes *sys,
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
solve_mix_3ft_noHI (struct stokes *sys,
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
solve_mix_3fts_noHI (struct stokes *sys,
		     const double *f, const double *t, const double *e,
		     const double *uf, const double *of, const double *ef,
		     double *u, double *o, double *s,
		     double *ff, double *tf, double *sf);


#endif /* !_LIBSTOKES_CORE_H_ */
