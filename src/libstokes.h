/* header file for library 'libstokes'
 * Copyright (C) 1993-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libstokes.h,v 1.59 2008/04/29 03:32:39 kichiki Exp $
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
  struct list_ex *ex_lub;

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
 ** NetCDF interface for libstokes **
 ** from stokes-nc.h               **
 ************************************/
struct stokes_nc {
  int id;

  int np; // only MOBILE particles here!!
  int p_dim;
  int p_id;

  int npf; // number of fixed particles (stored separately)
  int pf_dim;
  int pf_id;

  int nvec;
  int vec_dim;
  int vec_id;

  int nstt;
  int stt_dim;
  int stt_id;

  int nquat;
  int quat_dim;
  int quat_id;

  int ntime;
  int time_dim;
  int time_id;

  int l_id;
  int ui0_id;
  int oi0_id;
  int ei0_id;
  int ui_id;
  int oi_id;
  int ei_id;

  int x0_id;
  int u0_id;
  int o0_id;
  int e0_id;
  int f0_id;
  int t0_id;
  int s0_id;

  int xf0_id;
  int uf0_id;
  int of0_id;
  int ef0_id;
  int ff0_id;
  int tf0_id;
  int sf0_id;

  int x_id;
  int u_id;
  int o_id;
  int e_id;
  int f_id;
  int t_id;
  int s_id;
  int q_id;

  int xf_id;
  int uf_id;
  int of_id;
  int ef_id;
  int ff_id;
  int tf_id;
  int sf_id;

  int a_id;
  int af_id;

  /* auxiliary imposed-flow parameters for simple shear */
  int shear_mode_id;
  int shear_rate_id;
  int shear_shift_id;

  /* active/inactive flags : 0 = inactive
   *                         1 = active */
  int flag_ui0;
  int flag_oi0;
  int flag_ei0;
  int flag_ui;
  int flag_oi;
  int flag_ei;

  int flag_x0;
  int flag_u0;
  int flag_o0;
  int flag_e0;
  int flag_f0;
  int flag_t0;
  int flag_s0;

  int flag_xf0;
  int flag_uf0;
  int flag_of0;
  int flag_ef0;
  int flag_ff0;
  int flag_tf0;
  int flag_sf0;

  int flag_x;
  int flag_u;
  int flag_o;
  int flag_e;
  int flag_f;
  int flag_t;
  int flag_s;
  int flag_q;

  int flag_xf;
  int flag_uf;
  int flag_of;
  int flag_ef;
  int flag_ff;
  int flag_tf;
  int flag_sf;

  int flag_a;
  int flag_af;
};


void
stokes_nc_error (int status, const char * message, const char * varname);

void
stokes_nc_print_actives (struct stokes_nc * nc,
			 FILE * out);


/* initialize NetCDF file for libstokes for just x
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  the only activated entry is x.
 */
struct stokes_nc *
stokes_nc_x_init (const char * filename, int np);

/* initialize NetCDF file for libstokes
 * INPUT
 *  np : number of MOBILE particles
 *  nf : number of FIXED particles for "mix" (set 0 for "mob" problem)
 *  version   : 0 = F, 1 = FT, 2 = FTS
 *  flag_poly : 0 = monodisperse
 *              1 = polydisperse (set particle radius)
 *  flag_it   : 0 = constant imposed flow
 *              1 = time-changing imposed flow
 *  shear_mode : 0 == imposed flow is given by Ui,Oi,Ei.
 *               1 == x = flow dir, y = grad dir
 *               2 == x = flow dir, z = grad dir
 *               NOTE: shear_rate and shear_shift[t] are defined only for
 *               shear_mode = 1 or 2.
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are
 *   [F mob]      : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  f0, x, u.
 *   [FT mob]     : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  f0, t0, x, u, o.
 *   [FTS mob]    : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  f0, t0, e0, x, u, o, s.
 *   [F mix]      : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  xf0, f0, uf0, x, u, ff.
 *   [FT mix]     : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  xf0, f0, t0, uf0, of0, x, u, o, ff, tf.
 *   [FTS mix]    : ui0, oi0, ei0 (for it = 0) or ui, oi, ei (for it = 1) and
 *                  xf0, f0, t0, e0, uf0, of0, ef0, x, u, o, s, ff, tf, sf.
 *  for all cases, q is set if flag_Q == 1.
 */
struct stokes_nc *
stokes_nc_init (const char *filename, int np, int nf,
		int version,
		int flag_poly,
		int flag_Q,
		int flag_it,
		int shear_mode);

/* close (and write if necessary) NetCDF file for libstokes
 */
void
stokes_nc_free (struct stokes_nc * nc);


/** set nc data **/
/* set shear_rate (just a scalar)
 */
void
stokes_nc_set_shear_rate (struct stokes_nc *nc,
			  double shear_rate);
/* set l = (lx, ly, lz)
 */
void
stokes_nc_set_l (struct stokes_nc * nc,
		 const double * l);
/* set ui0[vec]
 */
void
stokes_nc_set_ui0 (struct stokes_nc * nc,
		   const double * ui0);
/* set oi[vec]
 */
void
stokes_nc_set_oi0 (struct stokes_nc * nc,
		   const double * oi0);
/* set ei0[stt]
 */
void
stokes_nc_set_ei0 (struct stokes_nc * nc,
		   const double * ei0);
/* set x0[np][vec]
 */
void
stokes_nc_set_x0 (struct stokes_nc * nc,
		  const double * x0);
/* set u0[np][vec]
 */
void
stokes_nc_set_u0 (struct stokes_nc * nc,
		  const double * u0);
/* set o0[np][vec]
 */
void
stokes_nc_set_o0 (struct stokes_nc * nc,
		  const double * o0);
/* set e0[np][stt]
 */
void
stokes_nc_set_e0 (struct stokes_nc * nc,
		  const double * e0);
/* set f0[np][vec]
 */
void
stokes_nc_set_f0 (struct stokes_nc * nc,
		  const double * f0);
/* set t0[np][vec]
 */
void
stokes_nc_set_t0 (struct stokes_nc * nc,
		  const double * t0);
/* set s0[np][stt]
 */
void
stokes_nc_set_s0 (struct stokes_nc * nc,
		  const double * s0);
/* set xf0[npf][vec]
 */
void
stokes_nc_set_xf0 (struct stokes_nc * nc,
		   const double * xf0);
/* set uf0[npf][vec]
 */
void
stokes_nc_set_uf0 (struct stokes_nc * nc,
		   const double * uf0);
/* set of0[npf][vec]
 */
void
stokes_nc_set_of0 (struct stokes_nc * nc,
		   const double * of0);
/* set ef0[npf][stt]
 */
void
stokes_nc_set_ef0 (struct stokes_nc * nc,
		   const double * ef0);
/* set ff0[npf][vec]
 */
void
stokes_nc_set_ff0 (struct stokes_nc * nc,
		   const double * ff0);
/* set tf0[npf][vec]
 */
void
stokes_nc_set_tf0 (struct stokes_nc * nc,
		   const double * tf0);
/* set sf0[npf][stt]
 */
void
stokes_nc_set_sf0 (struct stokes_nc * nc,
		   const double * sf0);
/* set a[np]
 */
void
stokes_nc_set_a (struct stokes_nc * nc,
		 const double * a);
/* set af[npf]
 */
void
stokes_nc_set_af (struct stokes_nc * nc,
		  const double * af);

/* set time (step)
 */
void
stokes_nc_set_time (struct stokes_nc * nc,
		    int step, double time);
/* set shear_shift[time]
 */
void
stokes_nc_set_shear_shift (struct stokes_nc *nc,
			   int step,
			   double shear_shift);
/* set ui[time][vec] at time (step)
 */
void
stokes_nc_set_ui (struct stokes_nc * nc,
		  int step,
		  const double * ui);
/* set oi[time][vec] at time (step)
 */
void
stokes_nc_set_oi (struct stokes_nc * nc,
		  int step,
		  const double * oi);
/* set ei[time][stt] at time (step)
 */
void
stokes_nc_set_ei (struct stokes_nc * nc,
		  int step,
		  const double * ei);
/* set x at time (step)
 */
void
stokes_nc_set_x (struct stokes_nc * nc,
		 int step,
		 const double * x);
/* set q at time (step)
 */
void
stokes_nc_set_q (struct stokes_nc * nc,
		 int step,
		 const double * q);
/* set u at time (step)
 */
void
stokes_nc_set_u (struct stokes_nc * nc,
		 int step,
		 const double * u);
/* set o at time (step)
 */
void
stokes_nc_set_o (struct stokes_nc * nc,
		 int step,
		 const double * o);
/* set e at time (step)
 */
void
stokes_nc_set_e (struct stokes_nc * nc,
		 int step,
		 const double * e);
/* set f at time (step)
 */
void
stokes_nc_set_f (struct stokes_nc * nc,
		 int step,
		 const double * f);
/* set t at time (step)
 */
void
stokes_nc_set_t (struct stokes_nc * nc,
		 int step,
		 const double * t);
/* set s at time (step)
 */
void
stokes_nc_set_s (struct stokes_nc * nc,
		 int step,
		 const double * s);;
/* set xf at time (step)
 */
void
stokes_nc_set_xf (struct stokes_nc * nc,
		  int step,
		  const double * xf);
/* set uf at time (step)
 */
void
stokes_nc_set_uf (struct stokes_nc * nc,
		  int step,
		  const double * uf);
/* set of at time (step)
 */
void
stokes_nc_set_of (struct stokes_nc * nc,
		  int step,
		  const double * of);
/* set ef at time (step)
 */
void
stokes_nc_set_ef (struct stokes_nc * nc,
		  int step,
		  const double * ef);
/* set ff at time (step)
 */
void
stokes_nc_set_ff (struct stokes_nc * nc,
		  int step,
		  const double * ff);
/* set tf at time (step)
 */
void
stokes_nc_set_tf (struct stokes_nc * nc,
		  int step,
		  const double * tf);
/* set sf at time (step)
 */
void
stokes_nc_set_sf (struct stokes_nc * nc,
		  int step,
		  const double * sf);


/** utility routines **/
/* set stokes_nc data by the parameters
 * INPUT
 *  sys        : the following elements are referred.
 *             :   version
 *             :   np, nm
 *             :   a[np] (if NULL, monodisperse mode)
 *             :   periodic
 *  flag_Q     :
 *  Ui, Oi, Ei :
 *  F,  T,  E  :
 *  uf, of, ef : used only for the mix problem
 *  xf         : position of the fixed particles
 *  lat[3]     : used only for the periodic case
 *  shear_mode : 0 == imposed flow is given by Ui,Oi,Ei.
 *               1 == x = flow dir, y = grad dir
 *               2 == x = flow dir, z = grad dir
 *  shear_rate : defined only for shear_mode = 1 or 2.
 */
struct stokes_nc *
stokes_nc_set_by_params (const char *out_file,
			 const struct stokes *sys,
			 int flag_Q,
			 const double *Ui, const double *Oi, const double *Ei,
			 const double *F,  const double *T,  const double *E,
			 const double *uf, const double *of, const double *ef,
			 const double *xf,
			 const double *lat,
			 int shear_mode, double shear_rate);

/* check stokes_nc data with the parameters
 * INPUT
 *  nc         :
 *  sys        : the following elements are referred.
 *             :   version
 *             :   np, nm
 *             :   a[np] (if NULL, monodisperse mode)
 *             :   periodic
 *             :   shear_mode, shear_rate
 *  flag_Q     :
 *  Ui, Oi, Ei :
 *  F,  T,  E  :
 *  uf, of, ef : used only for the mix problem
 *  xf         : position of the fixed particles
 *  lat[3]     : used only for the periodic case
 *  tiny       : small value for the check criteria
 */
int
stokes_nc_check_params (const struct stokes_nc *nc,
			const struct stokes *sys,
			int flag_Q,
			const double *Ui, const double *Oi, const double *Ei,
			const double *F, const double *T, const double *E,
			const double *uf, const double *of, const double *ef,
			const double *xf,
			const double *lat,
			double tiny);


/* from stokes-nc-read.h
 */

/* open stokes_nc file in NC_NOWRITE mode
 * this is for usual analysis
 */
struct stokes_nc *
stokes_nc_open (const char * filename);

/* open stokes_nc file in NC_WRITE mode
 * this is for continuation (appending the results)
 */
struct stokes_nc *
stokes_nc_reopen (const char * filename);

/* read 1d array [vec/stt/np/npf]
 * INPUT
 *  name : either one of them, Ui0, Oi0, Ei0, Ui, Oi, Ei, a, af, l
 * OUTPUT
 *  x[]
 */
void
stokes_nc_get_array1d (const struct stokes_nc * nc,
		       const char * name,
		       double * x);
/* read constant data for particles in 2d array [np/npf][vec/stt]
 */
void
stokes_nc_get_data0 (const struct stokes_nc * nc,
		     const char * name,
		     double * x);
/* read time-dep. particle data at step in 3d array [step][np/npf][vec/stt]
 */
void
stokes_nc_get_data (const struct stokes_nc * nc,
		    const char * name,
		    int step,
		    double * x);

/* read (the whole) time vector
 * INPUT
 *  time[nc->ntime]
 * OUTPUT
 *  time[nc->ntime]
 */
void
stokes_nc_get_time (const struct stokes_nc * nc,
		    double * time);

/* read time at a step
 * INPUT
 *  step
 * OUTPUT
 *  returned value : time[step]
 */
double
stokes_nc_get_time_step (const struct stokes_nc * nc,
			 int step);


/************************************
 ** Guile interface for libstokes  **
 ************************************/
/* from stokes-guile.h */
/* this utility function (original name is does_scm_symbol_exist)
 * written by Michael Gran is taken from
 * http://www.lonelycactus.com/guilebook/c319.html
 */
int
guile_check_symbol (const char *name);

/* check boolean
 * INPUT
 *  var : string of the SCM variable
 * OUTPUT
 *  true  : if var is not nil
 *  false : if var is nil or even not defined
 */
int
guile_get_bool (const char * var);

/*
 * INPUT
 *  var : string of the SCM variable
 *  i0  : default value (for the undefined case)
 * OUTPUT
 */
int
guile_get_int (const char * var, int i0);

/*
 * INPUT
 *  var : string of the SCM variable
 *  d0  : default value (for the undefined case)
 * OUTPUT
 */
double
guile_get_double (const char * var, double d0);

/* get doubles from SCM list or vector with length check
 * OUTPUT
 *  returned value : 0 = failed (not defined)
 *                   1 = success
 */
int
guile_get_doubles (const char * var, int n, double * x);
/* get doubles from SCM list or vector with unknown length
 * OUTPUT
 *  returned value : NULL = failed (not defined)
 */
double *
guile_get_doubles_ (const char * var);
/* get length of SCM list or vector
 * OUTPUT
 *  returned value : length (not defined, 0 is returned)
 */
int
guile_get_length (const char * var);

/*
 */
char *
guile_get_string (const char * var);

/*
 */
FILE *
guile_open_file (const char * var, const char * mode);


/* load scm script file into the libstokes C system
 */
void
guile_load (const char *file);


/************************************
 ** bond-interaction routines      **
 ************************************/
/* from bonds.h */
struct bond_pairs {
  int n;   // number of pairs
  int *ia; // particle a for the pair
  int *ib; // particle b for the pair
};

struct bonds {
  /* table for bond type */
  int n;      // number of bond types
  int *type;  /* type of the spring
	       * 0 : Hookean (p1 = k, spring constant,
	       *              p2 = r0, natural length)
	       * 1 : Wormlike chain (WLC)
	       * 2 : inverse Langevin chain (ILC)
	       * 3 : Cohen's Pade approx for ILC
	       * 4 : Werner spring (approx for ILC)
	       * 5 : another Hookean
	       *     where for these FENE chains,
	       *     p1 = N_{K,s} the Kuhn steps for a spring
	       *     p2 = b_{K}   the Kuhn length [nm]
	       * 6 : discrete Wormlike chain (dWLC), where
	       *     p1 = k, dimensionless spring constant, 
	       *     p2 = r0, the natural length [nm],
	       *     the potential is given by
	       *     U(r) = (k/2) * (kT / r0^2) * (r-r0)^2
	       */
  int *fene;   /* flag for the parameters p1 and p2
		* 0 : (p1, p2) = (A^{sp}, L_{s})
		* 1 : (p1, p2) = (N_{K,s}, b_{K})
		*/
  double *p1;  // the first parameter (k or N_{K,s})
  double *p2;  // the second parameter (r0 or b_{K})

  struct bond_pairs **pairs; // pairs for the bond

  int *nex;    // number of excluded particles in the chain
};

struct list_ex {
  int np;  // total number of particles (must be equal to sys->np)
  int *n;  // n[np] : number of excluded particles for each particles
  int **i; // i[j][k] : k-th particle index to exclude for particle j.
};



void
bond_pairs_free (struct bond_pairs *pairs);

void
bond_pairs_add (struct bond_pairs *pairs,
		int ia, int ib);


/* initialize struct bonds
 * INPUT
 * OUTPUT
 *  returned value : struct bonds
 */
struct bonds *
bonds_init (void);

void
bonds_free (struct bonds *bonds);

/* add a spring into bonds
 * INPUT
 *  bonds  : struct bonds
 *  type   : type of the spring
 *  fene   : 0 == (p1,p2) are (A^{sp}, L_{s})
 *           1 == (p1, p2) = (N_{K,s}, b_{K}) or
 *                (p1, p2) = (k, r0) for dWLC (type == 6).
 *                in the latter case, potential is given by
 *                (k/2) * (kT / r0^2) * (r-r0)^2
 *  p1, p2 : spring parameters
 *  nex    : number of excluded particles in the chain
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2,
		int nex);

/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used
 *          (bonds->fene[i] == 1 is expected), or 
 *          p1 (k) and p2 (r0) are used for dWLC spring (type == 6).
 *            in this case, potential is given by
 *            (k/2) * (kT / r0^2) * (r-r0)^2
 *  a     : length scale in the simulation
 *  pe    : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3a / pe b_{K}
 *  bonds->p2[] := Ls / a = N_{K,s} b_{K} / a 
 *    for dWLC spring, the conversions are given by 
 *  bonds->p1[] := A^{sp} = k / (pe * (r0/a)^2)
 *  bonds->p2[] := Ls / a = r0 / a
 */
void
bonds_set_FENE (struct bonds *bonds,
		double a, double pe);

/*
 * INPUT
 *  bonds      : struct bonds
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
bonds_calc_force (struct bonds *bonds,
		  struct stokes *sys,
		  double *f,
		  int flag_add);

void
bonds_print (FILE *out, struct bonds *bonds);


/**
 * SWIG utility routine
 * For examplean, expected usage in python by SWIG:
 *   n = stokes.bonds_get_pairs_n(bonds, i)
 *   ia = stokes.iarray(n)
 *   ib = stokes.iarray(n)
 *   stokes.bonds_get_pairs(bonds, i, ia, ib)
 * then, you have arrays ia[n] and ib[n].
 */

/* to get the number of pairs for the bond "i"
 */
int
bonds_get_pairs_n (struct bonds *b, int i);
/* 
 * INPUT
 *  b : struct bonds
 *  i : index of the bond
 *  ia[n] : "a" particle index of j-th pair for the bond "i",
 *  ib[n] : "b" particle index of j-th pair for the bond "i",
 *          where j runs from 0 to (n-1) and
 *          n = b->pairs[i]->n is the number of pairs for the bond "i".
 *          before calling, allocate ia and ib with (sizeof(int) * n).
 * OUTPUT
 *  ia[n], ib[n] : 
 */
void
bonds_get_pairs (struct bonds *b, int i,
		 int *ia, int *ib);


/**
 * exclusion list for lubrication due to the bonding
 */
struct list_ex *
list_ex_init (int np);

void
list_ex_add (struct list_ex *ex, int j, int k);

void
list_ex_free (struct list_ex *ex);

struct list_ex *
list_ex_copy (struct list_ex *ex0);

/* construct the excluded list by struct bonds
 */
void
list_ex_set_by_bonds (struct list_ex *ex, const struct bonds *b);

/* check whether j is excluded for i
 * INPUT
 *  ex : struct list_ex
 *  i  : particle now we are considering
 *  j  : particle whether it is in the list or not.
 * OUTPUT
 *  returned value : 0 (false); j is NOT in the excluded list for i.
 *                   1 (true);  j IS in the excluded list for i.
 */
int
list_ex_check (struct list_ex *ex, int i, int j);


/* from bonds-guile.h */
/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     0       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      0      ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
 *      1.0    ;    p1   = A^{sp}, scaled spring constant
 *      2.1)   ;    p2   = L_{s} / a, scaled max extension
 *     ((0 1)  ; 3) list of pairs
 *      (1 2)
 *      (2 3))
 *      -1)    ; 4) number of exclusion for lubrication
 *             ;    negative means all particles in the chain is excluded.
 *    (; bond 2
 *     2       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      1      ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
 *      19.8   ;    p1 = N_{K,s}, the Kuhn steps for a spring
 *      106.0) ;    p2 = b_{K} [nm], the Kuhn length
 *             ;    note that, for dWLC (type == 6),
 *             ;    (p1, p2) = (k, r0 [nm]), where the potential is
 *             ;    (k/2) * (kT / r0^2) * (r-r0)^2
 *     ((4 5)  ; 3) list of pairs
 *      (5 6)
 *      (6 7))
 *       1)    ; 4) number of exclusion for lubrication
 *   ))
 * where spring types are
 *   0 : Hookean spring (Asp * (r - Ls)
 *   1 : wormlike chain (WLC)
 *   2 : inverse Langevin chain (ILC)
 *   3 : Cohen's Pade approximation
 *   4 : Warner spring
 *   5 : Hookean spring (Asp * r / Ls)
 *   6 : Hookean spring for dWLC
 * OUTPUT
 *  returned value : struct bonds
 *                   if NULL is returned, it failed (not defined)
 */
struct bonds *
guile_get_bonds (const char * var);


/* from excluded-volume.h */
struct EV {
  double r2; // square of the max distance for F^{EV}

  /* table for chain type */
  int n;     // number of chain types
  double *l; // characteristic distance = (1/3) N_{K,s} b_{K}^2
  double *A; /* prefactor = (9/2) A^{sp} z', where
	      *   A^{sp} = 3 a / Pe b_{K},
	      *   z' = (N_{K,s}/2 pi)^{3/2} (v/l_s^3)
	      *      = (3 / 2 pi b_{K}^2)^{3/2} v
	      */

  /* table for particles */
  int *ch;   /* chain type for each particle
	      * negative value == no assignement to the chain
	      */
};

/* initialize struct EV
 * INPUT
 *  bonds  : struct bonds (either fene=0 or fene=1 is fine).
 *  a, pe  : parameters for bonds parameters
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring.
 *           the index should correspond to that in bonds.
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] characteristic distance = (1/3) N_{K,s} b_{K}^2,
 *      ev->A[i] prefactor = (9/2) A^{sp} z',
 *    and
 *      A^{sp} = 3 a / Pe b_{K},
 *      z' = (N_{K,s}/2 pi)^{3/2} (v/l_s^3)
 *         = (3 / 2 pi b_{K}^2)^{3/2} v.
 */
struct EV *
EV_init (const struct bonds *bonds, double a, double pe,
	 double r2, const double *v,
	 int np);

void
EV_free (struct EV *ev);


/* retrieve the EV parameters for particles i and j
 * with the combination rules 
 *   l(12) = (l1 + l2) / 2
 *   A(12) = (A1 * A2)^{1/2}
 * INPUT
 *  ev   : struct EV
 *  i, j : particle index (0 ~ np-1)
 * OUTPUT
 *  *l : 
 *  *A : 
 */
void
EV_get_coefficients (struct EV *ev,
		     int i, int j,
		     double *l, double *A);


/*
 * INPUT
 *  ev         : struct EV
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_calc_force (struct EV *ev,
	       struct stokes *sys,
	       double *f,
	       int flag_add);


/************************************
 ** angle-interaction routines     **
 ************************************/
/* from angles.h */
struct angle {
  double k;  // potential factor
  double t0; // natural angle (theta_0) in radian
  int scale; /* flag for the parameters (whether scaled or not)
	      * 0 : (k, t0) is scaled.
	      * 1 : (k, t0) is not scaled yet.
	      */

  // particle indices
  int ia;
  int ib;
  int ic;
};

struct angles {
  /* table for angle type */
  int n;           // number of angles
  struct angle *a; // angles
};


/**
 * struct angles
 */
struct angles *
angles_init (void);

void
angles_free (struct angles *ang);

/*
 * INPUT
 *  ia, ib, ic : particle indices (ib is the center particle)
 *  k  : potential factor
 *  t0 : natural angle (in radian)
 *  scale : flag for scale (0 == k is scaled,
 *                          1 == k is not scaled yet.)
 */
void
angles_add (struct angles *ang,
	    int ia, int ib, int ic,
	    double k, double t0, int scale);

/* scale parameter by the Peclet number
 * INPUT
 *  ang : struct angles
 *        (ang->a [n])->k is scaled as (k / pe)
 *  a     : length scale in the simulation
 *  pe    : peclet number
 * OUTPUT
 *  (ang->a [i])->k : scaled as (k / pe)
 */
void
angles_scale_k (struct angles *ang,
		double a, double pe);

/*
 * INPUT
 *  ang        : struct angles
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
angles_calc_force (struct angles *ang,
		   struct stokes *sys,
		   double *f,
		   int flag_add);


/* from angles-guile.h */
/* get angles from SCM
 * in SCM, angles are given by something like
 *  (define angles '(
 *    (; angle type 1
 *     10.0    ; 1) constant (k^{angle})
 *     0.0     ; 2) angle in degree (theta_0)
 *     0       ; 3) scale flag (0 == scaled)
 *             ;    in this case, the above value for k is just used.
 *     ((0 1 2); 4) list of triplets
 *      (1 2 3)
 *      (2 3 4)
 *     )
 *    )
 *    (; angle type 2
 *     20.0    ; 1) constant (k^{angle})
 *     90.0    ; 2) angle in degree (theta_0)
 *     1       ; 3) scale flag (1 == not scaled yet)
 *             ;    in this case, the potential is given by 
 *             ;    (k/2) * kT * (theta - theta_0)^2
 *     ((3 4 5); 4) list of triplets
 *      (4 5 6)
 *     )
 *    )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "angles".
 * OUTPUT
 *  returned value : struct angles
 *                   if NULL is returned, it failed (not defined)
 */
struct angles *
guile_get_angles (const char *var);


/*************************************
 ** excluded-volume by Debye-Huckel **
 *************************************/
/* from ev-dh.h */
struct EV_DH {
  /* system parameters */
  double r2;    /* square of the max distance for F^{EV}
	         * scaled by the characteristic length */
  double a_sys; /* := (1/pe)(1/kT)(e^2/4pi e0)(1/a) dimensionless number */
  double rd;    /* debye length scaled by the characteristic length */

  /* parameters for each chain */
  /* currently this is implemented particle-wise for simplicity */
  int n;
  double *nu; /* nu := (nu) l0 / e, the dimensionless charge density, where
	       *       (nu) is the line density of charge [e/nm]
	       *       l0   is the bond length [nm]
	       *       e    is the elementary charge [C] */
};


/* initialize struct EV
 * INPUT
 *  a  : characteristic length (in the same dimension for rd below, usually nm)
 *  pe : peclet number
 *  rd : Debye length (in the same dimension for a above, usually nm)
 *  T  : temperature in Kelvin.
 *  e  : dielectric constant of the solution
 *  r2 : square of the max distance for F^{EV_DH}
 *       (in the same dimension squared for a above, usually nm^2)
 *  np : number of particles
 * OUTPUT
 *  returned value : struct EV_DH,
 *      where only the system parameters (a_sys, rd) are defined.
 *      nu[i] is zero cleared.
 */
struct EV_DH *
EV_DH_init (double a, double pe, double rd, double T, double e,
	    double r2,
	    int np);

void
EV_DH_free (struct EV_DH *ev_dh);

/*
 * INPUT
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force (struct EV_DH *ev_dh,
		  struct stokes *sys,
		  double *f,
		  int flag_add);

/* from ev-dh-guile.h */
/* get ev_dh from SCM
 * in SCM, angles are given by something like
 *  (define ev-dh '(
 *    ; system parameters
 *    4.0      ; 1) max distance for EV_DH interaction [nm]
 *    298.0    ; 2) temperature [K]
 *    80.0     ; 3) dielectric constant of the solution
 *    3.07     ; 4) Debye length [nm]
 *    (        ; 5) list of chain types
 *     (; chain type 1
 *      2.43    ; 1) nu [e/nm]
 *      5.00    ; 2) l0 [nm]
 *      (0 1 2) ; 3) list of particles
 *     )
 *     (; chain type 2
 *      2.00    ; 1) nu [e/nm]
 *      4.00    ; 2) l0 [nm]
 *      (3 4)   ; 3) list of particles
 *     )
 *    )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-dh".
 *  a   : the characteristic length [nm]
 *  pe  : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_DH
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_DH *
EV_DH_guile_get (const char *var,
		 double a, double pe, int np);


/************************************
 ** routines for ODE integrator    **
 ************************************/
/* from ode.h */
struct ode_params
{
  struct stokes *sys;
  double *F;
  double *T;
  double *E;
  double *uf;
  double *of;
  double *ef;
  int flag_noHI;
  int flag_mat;
  int flag_lub;
  double st;
  struct bonds *bonds;
  double gamma;

  // auxiliary imposed-flow parameters for simple shear
  double t0; // reference time for s0
  double s0; // cell shift at time t0 (for shear_mode = 1 or 2)
};


/* set the parameters to struct ode_params
 * INPUT
 *  (struct stokes *)sys -- initialize before calling!
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_mat
 *  (int) flag_lub
 *  (double) stokes
 *  (struct bonds *)bonds
 *  (double) gamma (for the bond relaxation scheme)
 * OUTPUT :
 *  (struct ode_params) params
 */
struct ode_params *
ode_params_init (struct stokes *sys,
		 double *F,
		 double *T,
		 double *E,
		 double *uf,
		 double *of,
		 double *ef,
		 int flag_noHI,
		 int flag_lub,
		 int flag_mat,
		 double st,
		 struct bonds *bonds,
		 double gamma);

void 
ode_params_free (struct ode_params *params);

/* set the reference for cell-shift (shear_mode = 1 and 2)
 * INTPUT
 *  t0 : reference time for s0
 *  s0 : cell shift at time t0 (for shear_mode = 1 or 2)
 * OUTPUT
 *  ode->t0, ode->s0 :
 */
void
ode_set_shear_shift_ref (struct ode_params *ode,
			 double t0, double s0);


/* calc dydt for gsl_odeiv ONLY with bond interactions for relaxation
 * this is equivalent to dydt_relax_bond with a constant friction, where
 *  dx/dt = U = f_bond / gamma
 * INPUT
 *  t   : current time
 *  y[] : position of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->bonds     : (struct bonds *)
 *           ode->gamma     : the friction coef
 * OUTPUT
 *  f[] := (d/dt) y (t), that is, the velocity of particles at time t
 */
int
dydt_relax_bond (double t, const double *y, double *f,
		 void *params);

/* calc dydt for gsl_odeiv with hydrodynamics
 * with inertia (scaled Stokes number)
 * INPUT
 *  t   : current time
 *  y[] : position and (real) velocity of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->stokes
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) y(t), where y(t) = (x(t), U(t)).
 *         (d/dt) x(t) = U(t)
 *         (d/dt) U(t) = (1/stokes)(-U(t) + V(t)),
 *         where V(t) := R^-1 . F, the terminal velocity
 */
int
dydt_hydro_st (double t, const double *y, double *dydt,
	       void *params);

/* calc dydt for gsl_odeiv with hydrodynamics
 * without inertia (zero Stokes number)
 * INPUT
 *  t   : current time
 *  y[] : position of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) x(t) = U(t)
 *         where U(t) := R^-1 . F, the terminal velocity
 */
int
dydt_hydro (double t, const double *y, double *dydt,
	    void *params);


/* from ode-quaternion.h */
/* calc dydt for center and angle by gsl_odeiv with hydrodynamics
 * without inertia (zero Stokes number)
 * INPUT
 *  t   : current time
 *  y[nm*3 + nm*4] : center position and quaternion of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) x(t) = U(t)
 *         where U(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro (double t, const double *y, double *dydt,
	      void *params);

/* calc dydt for center and angle by gsl_odeiv with hydrodynamics
 * with inertia (scaled Stokes number)
 * INPUT
 *  t   : current time
 *  y[nm*3 + nm*3 + nm*4] :
 *           position, (real) velocity, and quaternion of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->stokes
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) y(t), where y(t) = (x(t), U(t)).
 *         (d/dt) x(t) = U(t)
 *         (d/dt) U(t) = (1/stokes)(-U(t) + V(t)),
 *         where V(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro_st (double t, const double *y, double *dydt,
		 void *params);


/************************************
 ** routines for Brownian dynamics **
 ************************************/
/* from brownian.h */
struct BD_params
{
  /* note that the following pointers are just pointers, therefore, 
   * you have to take care of them (to free, for example).
   */
  struct stokes *sys;
  struct KIrand *rng;
  double *F;
  double *T;
  double *E;
  double *uf;
  double *of;
  double *ef;

  int flag_noHI;
  int flag_mat;
  int flag_lub;
  int flag_lub_B; // for calc_brownian_force(), lub among mobile particles

  // auxiliary imposed-flow parameters for simple shear
  double t0; // reference time for s0
  double s0; // cell shift at time t0 (for shear_mode = 1 or 2)

  double st; // currently this is just place holders

  struct bonds *bonds;
  double gamma;
  struct EV *ev;
  struct angles *ang;
  struct EV_DH *ev_dh;

  int flag_Q;

  // parameters for Brownian dynamics
  double peclet;
  double eps;

  int n_minv;
  double eig_minv[2];
  double *a_minv;
  int n_lub;
  double eig_lub[2];
  double *a_lub;

  int scheme;  /* 0 : the mid-point algorithm
		* 1 : Banchio-Brady (2003)
		* 2 : Ball-Melrose (1997)
		* 3 : Jendrejack et al (2000)
		* 4 : semi-implicit predictor-corrector
		* note that, for 3 and 4, BD_imp_ode_evolve() should be used.
		*/
  double BB_n; // step parameter for BB03 algorithm

  double rmin;   /* minimum distance for dt-adjustment
		  * if (r < rmin*(ai+aj)), dt is adjusted
		  * rmin == 0 corresponds to "no dt-adjustment"
		  * rmin == 1 corresponds to dt-adjustment for overlap
		  * NOTE: if sys->rmin is defined, dt-adjustment is ignored.
		  */
  double dt_lim; /* lower bound to shrink dt to prevent overlaps
		  * set "dt" if you don't want to adjust dt but just reject
		  */
};


/* set the parameters to struct BD_params
 * INPUT
 *  ** NOTE ** the following pointers are just pointers.
 *             you have to take care of them! (free, for example.)
 *  (struct stokes *)sys -- initialize before calling!
 *  seed : for random number generator
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_noHI
 *  (int) flag_lub
 *  (int) flag_mat
 *        NOTE, flag_lub_B is used for calc_brownian_force() where
 *        lub among ONLY mobile particles are taken.
 *        therefore, check for the existance of mobile pair(s)
 *        not excluded by sys->ex_lub for flag_lub_B.
 *  (double) stokes -- currently this is just a place holder
 *  (struct bonds *)bonds
 *  (double) gamma
 *  (struct EV *)ev
 *  (struct angles *)ang
 *  (struct EV_DH *)ev_dh
 *  (int) flag_Q
 *  (double) peclet
 *  (double) eps
 *  (int) n_minv
 *  (int) n_lub
 *  (int) scheme
 *  (double) BB_n
 *  (double) rmin
 *  (double) dt_lim
 * OUTPUT :
 *  (struct ode_params) params
 */
struct BD_params *
BD_params_init (struct stokes *sys,
		unsigned long seed,
		double *F,
		double *T,
		double *E,
		double *uf,
		double *of,
		double *ef,
		int flag_noHI,
		int flag_lub,
		int flag_mat,
		double st,
		struct bonds *bonds,
		double gamma,
		struct EV *ev,
		struct angles *ang,
		struct EV_DH *ev_dh,
		int flag_Q,
		double peclet,
		double eps,
		int    n_minv,
		int    n_lub,
		int    scheme,
		double BB_n,
		double rmin,
		double dt_lim);

void
BD_params_free (struct BD_params *BD);

/* set the reference for cell-shift (shear_mode = 1 and 2)
 * INTPUT
 *  t0 : reference time for s0
 *  s0 : cell shift at time t0 (for shear_mode = 1 or 2)
 * OUTPUT
 *  BD->t0, BD->s0 :
 */
void
BD_set_shear_shift_ref (struct BD_params *BD,
			double t0, double s0);

/* wrapper for BD_evolve()
 * INPUT
 *  *t    : (input) current time
 *  t_out : output time
 *  *dt   : (input) current inner time-step
 *  y[n]  : (input) current configuration at (*t),
 *          where n = nm*3 for F version
 *                    (y[] in the first (nm*3) elements = velocity),
 *                n = nm*3 + nm*4 with quaternion
 *                    (y[] in the first (nm*3) elements = velocity,
 *                     y[] in the next (nm*4) elements = quaternion).
 * OUTPUT
 *  *t    : (output) output time (= t_out)
 *  *dt   : (output) current (updated, if necessary) inner time-step
 *  y[n]  : (output) updated configuration at t_out
 */
void
BD_ode_evolve (struct BD_params *BD,
	       double *t, double t_out, double *dt,
	       double *y);


/* from bd-imp.h */
#include <gsl/gsl_multiroots.h>
struct BD_imp {
  struct BD_params *BD;
  double dt;
  double *x0;
  double *q0;
  double fact;
  double *z;

  // GSL stuff
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_function *F;
  gsl_multiroot_fsolver  *S;
  gsl_vector *guess;
  int itmax;
  double eps;

  // working area used in BD_imp_JGdP00_func()
  struct FTS *FTS;
  double *pos;
  double *q;

  // for predictor-corrector algorithm
  int flag_PC;
  double *xP;
  double *qP;
  double *uP;
  double *dQdtP;
};


/* initialize struct BD_imp
 * INPUT
 *  BD : struct BD_params
 *       note that BDimp->BD is just a pointer to BD in the argument.
 *  itmax : max of iteration for the root-finding
 *  eps   : tolerance for the root-finding
 */
struct BD_imp *
BD_imp_init (struct BD_params *BD,
	     int itmax, double eps);

void
BD_imp_free (struct BD_imp *BDimp);

/* wrapper for BD_imp_evolve()
 * INPUT
 *  *t    : (input) current time
 *  t_out : output time
 *  *dt   : (input) current inner time-step
 *  y[n]  : (input) current configuration at (*t),
 *          where n = nm*3 for F version
 *                    (y[] in the first (nm*3) elements = velocity),
 *                n = nm*3 + nm*4 with quaternion
 *                    (y[] in the first (nm*3) elements = velocity,
 *                     y[] in the next (nm*4) elements = quaternion).
 * OUTPUT
 *  *t    : (output) output time (= t_out)
 *  *dt   : (output) current (updated, if necessary) inner time-step
 *  y[n]  : (output) updated configuration at t_out
 */
void
BD_imp_ode_evolve (struct BD_imp *BDimp,
		   double *t, double t_out, double *dt,
		   double *y);


/************************************
 ** miscellaneous routines         **
 ************************************/
/* from coll.h */
/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 3] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_particles (struct stokes *sys,
		   const double *x, double *v, double en);

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 3] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 *  x_wall : position of the wall
 *  v_wall : 
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_wall_x (struct stokes *sys,
		const double *x, double *v, double en,
		double x_wall, double v_wall);
void
collide_wall_y (struct stokes *sys,
		const double *x, double *v, double en,
		double y_wall, double v_wall);
void
collide_wall_z (struct stokes *sys,
		const double *x, double *v, double en,
		double z_wall, double v_wall);

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 2] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 2] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 2] : velocity of particles after collisions
 */
void
collide_particles_2d (struct stokes *sys,
		      const double *x, double *v, double en);

/* from periodicity.h */
/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 3] : position to check; the range is [0, l[xyz])
 */
void
check_periodic (struct stokes * sys,
		double *x);

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np * 3] : angle to check; the range is [0, 2\pi)
 */
void
check_angle (struct stokes * sys, double *angle);

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 2] : position to check; the range is [0, l[xy])
 */
void
check_periodic_2d (struct stokes * sys,
		   double *x);
/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np] : angle to check; the range is [0, 2\pi)
 */
void
check_angle_2d (struct stokes * sys, double *angle);

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


#endif /* !_LIBSTOKES_H_ */
