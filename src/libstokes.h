/* header file for library 'libstokes'
 * Copyright (C) 1993-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libstokes.h,v 1.67 2008/07/16 16:43:32 kichiki Exp $
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
 ************************************/
/* from KIrand.h
 */
#define MTRNG_N 624

struct KIrand {
  unsigned long mt[MTRNG_N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */

  /* for KIrand_Gaussian() */
  int Gaussian_has_saved;
  double Gaussian_saved;
};

/* from stokes-nc.h
 */
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

  int nrng;
  int rng_dim;
  int rng_id;
  int flag_rng;

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

  /* states of the random number generator */
  int mt_id;
  int mti_id;
  int mt_Ghs_id;
  int mt_Gs_id;

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
stokes_nc_print_actives (struct stokes_nc *nc,
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
 *  flag_rng : 0 == no states of random number generator
 *             1 == output states of random number generator
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
		int shear_mode,
		int flag_rng);

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
stokes_nc_set_l (struct stokes_nc *nc,
		 const double * l);
/* set ui0[vec]
 */
void
stokes_nc_set_ui0 (struct stokes_nc *nc,
		   const double * ui0);
/* set oi[vec]
 */
void
stokes_nc_set_oi0 (struct stokes_nc *nc,
		   const double * oi0);
/* set ei0[stt]
 */
void
stokes_nc_set_ei0 (struct stokes_nc *nc,
		   const double * ei0);
/* set x0[np][vec]
 */
void
stokes_nc_set_x0 (struct stokes_nc *nc,
		  const double * x0);
/* set u0[np][vec]
 */
void
stokes_nc_set_u0 (struct stokes_nc *nc,
		  const double * u0);
/* set o0[np][vec]
 */
void
stokes_nc_set_o0 (struct stokes_nc *nc,
		  const double * o0);
/* set e0[np][stt]
 */
void
stokes_nc_set_e0 (struct stokes_nc *nc,
		  const double * e0);
/* set f0[np][vec]
 */
void
stokes_nc_set_f0 (struct stokes_nc *nc,
		  const double * f0);
/* set t0[np][vec]
 */
void
stokes_nc_set_t0 (struct stokes_nc *nc,
		  const double * t0);
/* set s0[np][stt]
 */
void
stokes_nc_set_s0 (struct stokes_nc *nc,
		  const double * s0);
/* set xf0[npf][vec]
 */
void
stokes_nc_set_xf0 (struct stokes_nc *nc,
		   const double * xf0);
/* set uf0[npf][vec]
 */
void
stokes_nc_set_uf0 (struct stokes_nc *nc,
		   const double * uf0);
/* set of0[npf][vec]
 */
void
stokes_nc_set_of0 (struct stokes_nc *nc,
		   const double * of0);
/* set ef0[npf][stt]
 */
void
stokes_nc_set_ef0 (struct stokes_nc *nc,
		   const double * ef0);
/* set ff0[npf][vec]
 */
void
stokes_nc_set_ff0 (struct stokes_nc *nc,
		   const double * ff0);
/* set tf0[npf][vec]
 */
void
stokes_nc_set_tf0 (struct stokes_nc *nc,
		   const double * tf0);
/* set sf0[npf][stt]
 */
void
stokes_nc_set_sf0 (struct stokes_nc *nc,
		   const double * sf0);
/* set a[np]
 */
void
stokes_nc_set_a (struct stokes_nc *nc,
		 const double * a);
/* set af[npf]
 */
void
stokes_nc_set_af (struct stokes_nc *nc,
		  const double * af);

/* set time (step)
 */
void
stokes_nc_set_time (struct stokes_nc *nc,
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
stokes_nc_set_ui (struct stokes_nc *nc,
		  int step,
		  const double * ui);
/* set oi[time][vec] at time (step)
 */
void
stokes_nc_set_oi (struct stokes_nc *nc,
		  int step,
		  const double * oi);
/* set ei[time][stt] at time (step)
 */
void
stokes_nc_set_ei (struct stokes_nc *nc,
		  int step,
		  const double * ei);
/* set x at time (step)
 */
void
stokes_nc_set_x (struct stokes_nc *nc,
		 int step,
		 const double * x);
/* set q at time (step)
 */
void
stokes_nc_set_q (struct stokes_nc *nc,
		 int step,
		 const double * q);
/* set u at time (step)
 */
void
stokes_nc_set_u (struct stokes_nc *nc,
		 int step,
		 const double * u);
/* set o at time (step)
 */
void
stokes_nc_set_o (struct stokes_nc *nc,
		 int step,
		 const double * o);
/* set e at time (step)
 */
void
stokes_nc_set_e (struct stokes_nc *nc,
		 int step,
		 const double * e);
/* set f at time (step)
 */
void
stokes_nc_set_f (struct stokes_nc *nc,
		 int step,
		 const double * f);
/* set t at time (step)
 */
void
stokes_nc_set_t (struct stokes_nc *nc,
		 int step,
		 const double * t);
/* set s at time (step)
 */
void
stokes_nc_set_s (struct stokes_nc *nc,
		 int step,
		 const double * s);
/* set xf at time (step)
 */
void
stokes_nc_set_xf (struct stokes_nc *nc,
		  int step,
		  const double * xf);
/* set uf at time (step)
 */
void
stokes_nc_set_uf (struct stokes_nc *nc,
		  int step,
		  const double * uf);
/* set of at time (step)
 */
void
stokes_nc_set_of (struct stokes_nc *nc,
		  int step,
		  const double * of);
/* set ef at time (step)
 */
void
stokes_nc_set_ef (struct stokes_nc *nc,
		  int step,
		  const double * ef);
/* set ff at time (step)
 */
void
stokes_nc_set_ff (struct stokes_nc *nc,
		  int step,
		  const double * ff);
/* set tf at time (step)
 */
void
stokes_nc_set_tf (struct stokes_nc *nc,
		  int step,
		  const double * tf);
/* set sf at time (step)
 */
void
stokes_nc_set_sf (struct stokes_nc *nc,
		  int step,
		  const double * sf);


/* set rng data at time (step)
 */
void
stokes_nc_set_rng (struct stokes_nc *nc,
		   int step,
		   const struct KIrand *rng);


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
 *  flag_rng   :
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
			 int shear_mode, double shear_rate,
			 int flag_rng);

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
stokes_nc_get_array1d (const struct stokes_nc *nc,
		       const char * name,
		       double * x);
/* read constant data for particles in 2d array [np/npf][vec/stt]
 */
void
stokes_nc_get_data0 (const struct stokes_nc *nc,
		     const char * name,
		     double * x);
/* read time-dep. particle data at step in 3d array [step][np/npf][vec/stt]
 */
void
stokes_nc_get_data (const struct stokes_nc *nc,
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
stokes_nc_get_time (const struct stokes_nc *nc,
		    double * time);

/* read time at a step
 * INPUT
 *  step
 * OUTPUT
 *  returned value : time[step]
 */
double
stokes_nc_get_time_step (const struct stokes_nc *nc,
			 int step);

/* read rng data at time (step)
 * INPUT
 *  step
 */
void
stokes_nc_get_rng (struct stokes_nc *nc,
		   int step,
		   struct KIrand *rng);


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
 ** constraints                    **
 ************************************/
/* from bead-rod.h */
struct BeadRod {
  struct stokes *sys;

  int nc; // number of constraints

  // particle indices consisting of the rods
  double *a2;  // a2[nc] : square rod distance
  int *ia; // ia[nc] : particle at one end of the rod
  int *ib; // ib[nc] : particle at the other end of the rod

  int verbose; // 0 == quiet, 1 == verbose

  int scheme; /* 0 == iterative (Liu 1989)
	       * 1 == NITSOL (Newton-GMRES)
	       */
  double eps; // tolerance

  struct NITSOL *nit; // for scheme == 1
  struct iter *it;    // for scheme == 0

  // dynamic properties
  double dt; // time difference between u[] and uu[] below
  double d1; // = (2 dt / zeta)
  double d2; // = (2 dt / zeta)^2
  // note : the linear-term coefficient is (2 * d1)

  double *u;  // u [nc * 3], initial connector at t
  double *uu; // uu[nc * 3], connector without constraint at t+dt
};


/* initialize struct BeadRod
 * INPUT
 *  sys   : struct stokes
 *  nc    : number of constraints
 *  a[nc] : distances for each constraint
 *  ia[nc], ib[nc] : particle indices for each constraint
 */
struct BeadRod *
BeadRod_init (struct stokes *sys,
	      int nc,
	      const double *a,
	      const int *ia,
	      const int *ib);

void
BeadRod_free (struct BeadRod *br);


/* from bead-rod-guile.h */
/* get constraints from SCM and set struct BeadRod
 * "constraints" is given in SCM as 
 *  (define constraints '(
 *   ; system parameters
 *   1.0e-6    ; 1) tolerance
 *   "nitsol"  ; 2) scheme for solving nonlinear equations
 *             ;    "linear" for iterative scheme in linear approximation
 *             ;    "nitsol" for Newton-GMRES scheme by NITSOL library
 *   ; the following is for each constraint
 *   (         ; 3) constraint type 1
 *    5.0      ; 3-1) distance [nm]
 *    (        ; 3-2) list of particle-pairs
 *     (0 1)
 *     (1 2)
 *     (2 3)
 *   ))
 *   (         ; 4) constraint type 2
 *    10.0     ; 4-1) distance [nm]
 *    (        ; 4-2) list of particle-pairs
 *     (3 4)
 *     (4 5)
 *   ))
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "constraints".
 *  sys : struct stokes
 *  length : unit length in the simulation
 * OUTPUT
 *  returned value : struct BeadRod
 *                   if NULL is returned, it failed (not defined)
 */
struct BeadRod *
BeadRod_guile_get (const char *var,
		   struct stokes *sys,
		   double length);


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
  double *p3;  // the third parameter (tol for FENE-Fraenkel)

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
 *  p1, p2, p3 : spring parameters
 *  nex    : number of excluded particles in the chain
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2, double p3,
		int nex);

/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used
 *          (bonds->fene[i] == 1 is expected), or 
 *          p1 (k) and p2 (r0) are used for dWLC spring (type == 6).
 *            in this case, potential is given by
 *            (k/2) * (kT / r0^2) * (r-r0)^2
 *  length : length scale in the simulation
 *  peclet : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3 length / (peclet b_{K})
 *  bonds->p2[] := Ls / length = N_{K,s} b_{K} / length
 *    or, for dWLC spring, the conversions are given by 
 *  bonds->p1[] := A^{sp} = k / (pe * (r0/length)^2)
 *  bonds->p2[] := Ls / length = r0 / length
 */
void
bonds_set_FENE (struct bonds *bonds,
		double length, double peclet);

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
 *      1)     ; 4) number of exclusion for lubrication
 *    (; bond 3
 *     7       ; 1) spring type (FENE-Fraenkel)
 *     (       ; 2) spring parameters (list with 4 elements)
 *      0      ;    fene = 0 means (p1, p2, p3) = (H, r0 [nm], tol)
 *      1.0e6  ;    p1 = H, the spring constant
 *      0.5    ;    p2 = r0 [nm], the natural length of the spring
 *      0.01)  ;    p3 = tol, the tolerance parameter "s"
 *             ;    note that, for FENE-Fraenkel (type == 7),
 *             ;    the scalar part of the force is
 *             ;    fr = H * (r/hat(r0) - 1.0) / (1 - ((1-r/hat(r0))/tol)^2)
 *             ;    where hat(r0) = r0 / L0 (L0 is given by "length" [nm])
 *     ((8 9)  ; 3) list of pairs
 *      (9 10))
 *      1)     ; 4) number of exclusion for lubrication
 *   ))
 * where spring types are
 *   0 : Hookean spring (Asp * (r - Ls))
 *   1 : wormlike chain (WLC)
 *   2 : inverse Langevin chain (ILC)
 *   3 : Cohen's Pade approximation
 *   4 : Warner spring
 *   5 : Hookean spring (Asp * r / Ls)
 *   6 : Hookean spring for dWLC
 *   7 : FENE-Fraenkel
 * OUTPUT
 *  returned value : struct bonds
 *                   if NULL is returned, it failed (not defined)
 */
struct bonds *
bonds_guile_get (const char * var);


/* from excluded-volume.h */
struct EV {
  double r2; // square of the max distance for F^{EV}

  /* table for chain type */
  int n;     // number of chain types
  double *l; /* := sqrt((2/3) hat(ls)^2) */
  double *A; /* := (1/pi)^{3/2} hat(v) N_Ks^2 / (peclet ev->l[i]^5)
	      * where
	      *   ls^2     = N_Ks * b_K^2 / 3 (dimensional, so as b_K)
	      *   hat(ls)  = ls / length      (dimensionless)
	      *   hat(v)   = v / length^3     (dimensionless)
	      * with which the force F^EV is given by 
	      *   F^{EV}_{i} = A * r_{ij} * exp (- r_{ij}^2 / l^2)
	      * (note: r_{ij} is also dimensionless scaled by length)
	      */

  /* table for particles */
  int *ch;   /* chain type for each particle
	      * negative value == no assignement to the chain
	      */
};

/* initialize struct EV
 * INPUT
 *  bonds  : struct bonds (either fene=0 or fene=1 is fine).
 *  length : unit of length given by "length" in SCM (dimensional number)
 *  peclet : peclet number (with respect to "length")
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring.
 *           the index should correspond to that in bonds.
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] = sqrt((2/3) hat(ls)^2)
 *      ev->A[i] = (1/pi)^{3/2} hat(v) N_Ks^2 / (peclet ev->l[i]^5)
 *    where
 *      ls^2     = N_Ks * b_K^2 / 3 (dimensional, so as b_K)
 *      hat(ls)  = ls / length      (dimensionless)
 *      hat(v)   = v / length^3     (dimensionless)
 *    with which the force F^EV is given by 
 *      F^{EV}_{i} = A * r_{ij} * exp (- r_{ij}^2 / l^2)
 *    (note: r_{ij} is also dimensionless scaled by length)
 */
struct EV *
EV_init (const struct bonds *bonds,
	 double length, double peclet,
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


/* from excluded-volume-guile.h */
/* get ev-v from SCM and set struct EV
 * in SCM, ev-v is a list of parameter v [nm^3] or [micro m^3]
 * (depending on the dimension of the parameter "length")
 * for each spring:
 *  (define ev-v '(
 *   0.0012 ; for the spring 1
 *   0.002  ; for the spring 2
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-v".
 *  bonds : struct bonds
 *  length : unit length given by "length" in SCM (dimensional value)
 *  peclet : peclet number
 *  ev_r2  : square of max distance for EV interaction
 *  np     : number of particles (beads)
 * OUTPUT
 *  returned value : struct EV
 *                   if NULL is returned, it failed (not defined)
 */
struct EV *
EV_guile_get (const char *var,
	      const struct bonds *bonds,
	      double length, double peclet,
	      double ev_r2, int np);


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
angles_guile_get (const char *var);


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
 *  length : characteristic length
 *           (in the same dimension for rd below, usually nm)
 *  peclet : peclet number
 *  rd     : Debye length (in the same dimension for a above, usually nm)
 *  T      : temperature in Kelvin.
 *  e      : dielectric constant of the solution (dimensionless number)
 *  eps    : to determine cut-off distance through eps = exp (-r_cutoff / rd) 
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV_DH,
 *      where only the system parameters (a_sys, rd) are defined.
 *      nu[i] is zero cleared.
 */
struct EV_DH *
EV_DH_init (double length, double peclet,
	    double rd, double T, double e,
	    double eps,
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
 * in SCM, "ev-dh" are given by something like
 *  (define ev-dh '(
 *    ; system parameters
 *    1.0e-6   ; 1) epsilon for the cut-off distance of EV_DH interaction
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
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  length : the characteristic length (dimensional, usually in nm)
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_DH
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_DH *
EV_DH_guile_get (const char *var,
		 double length, double peclet, int np);


/**************************************
 ** excluded-volume by Lennard-Jones **
 **************************************/
/* from ev-LJ.h */
struct EV_LJ {
  /* LJ parameters, where the potential is given in DIMENSIONAL FORM as 
   *  U(r) = e ((r0/r)^12 - 2(r0/r)^6)
   * Note that r0 = 2^{1/6} sigma is the distance of potential minimum.
   * The repulsive force is evaluated by taking r = r0 at the contact point.
   * The parameters e and r0 should be dimensionless numbers defined by 
   *         e = hat(e) for flag == 0, or
   *         e = hat(e) / (peclet * hat(r0)) for flag == 1,
   *         r0 = hat(r0) = "r0" / length
   * where
   *         hat(e)  = "e" / kT
   *         hat(r0) = "r0" / length
   * and "e" and "r0" are dimensional numbers.
   * The dimensionless force (scalar part) is then given by
   *  hat(F) = 12 e (r0/r)^7 ((r0/r)^6 - 1)
   * where e is in the form of flag == 1.
   */
  /* parameters for each particle */
  int n;      // number of particles
  int    *flag;
  double *e;  /* always dimensionless
	       * either e/kT              for flag==0
	       * or     e/(kT Pe hat(r0)) for flag==1
	       */
  double *r0; // always dimensionless (scaled by "length")
};


/* initialize struct EV_LJ
 * INPUT
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV_LJ,
 *      where LJ parameters are set by zero.
 */
struct EV_LJ *
EV_LJ_init (int np);

void
EV_LJ_free (struct EV_LJ *ev_LJ);

/* scale LJ parameter e for runs
 * INPUT
 *  ev_LJ  : struct EV_LJ
 *  peclet : peclet number
 */
void
EV_LJ_scale (struct EV_LJ *ev_LJ,
	     double peclet);

/*
 * INPUT
 *  ev_LJ      : struct EV_LJ
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_LJ_calc_force (struct EV_LJ *ev_LJ,
		  struct stokes *sys,
		  double *f,
		  int flag_add);

/* from ev-LJ-guile.h */
/* get ev-LJ from SCM
 * in SCM, "ev-LJ" are given by something like
 *  (define ev-LJ '(
 *   (; LJ type 1
 *    10.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     0 1 2
 *    )
 *   )
 *   (; LJ type 2
 *    8.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    2.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     3 4
 *    )
 *   )
 *  ))
 * INPUT
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_LJ
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_LJ *
EV_LJ_guile_get (const char *var, int np);


/*************************************
 ** confinement forces              **
 *************************************/
/* from confinement.h */
struct confinement {
  int type; /* 0 : sphere
	     *       R : radius of cavity centered at (0, 0, 0)
	     * 1 : sphere + hole
	     *       R : radius of cavity centered at (0, 0, 0)
	     *       r : radius of the hole at (0, 0, 1) direction
	     * 2 : cylinder
	     *       r       : radius of the cylinder
	     *       x, y, z : direction vector of the cylinder
	     *       the cylinder center goes through (0, 0, 0) and (x, y, z).
	     * 3 : dumbbell
	     *       R  : left cavity radius centered at (center1, 0, 0)
	     *       R2 : right cavity radius centered at (center2, 0, 0)
	     *       L  : length of the cylinder
	     *       r  : cylinder radius
	     *       the origin is at the center of the cylinder
	     * 4 : hex2d
	     *       R : cavity radius
	     *       r : cylinder radius
	     *       L : lattice spacing
	     * 5 == porous
	     *       R : particle radius
	     *       L : lattice spacing in x (2R for touching case)
	     */
  // primary parameters
  double R;
  double r;
  double x, y, z;
  double R2;
  double L;

  /* LJ parameters, where the potential is given in DIMENSIONAL FORM as 
   *  U(r) = e ((r0/r)^12 - 2(r0/r)^6)
   * Note that r0 = 2^{1/6} sigma is the distance of potential minimum.
   * The repulsive force is evaluated by taking r = r0 at the contact point.
   * The parameters e and r0 should be dimensionless numbers defined by 
   *         e = hat(e) for flag_LJ == 0, or
   *         e = hat(e) / (peclet * hat(r0)) for flag_LJ == 1,
   *         r0 = hat(r0) = "r0" / length
   * where
   *         hat(e)  = "e" / kT
   *         hat(r0) = "r0" / length
   * and "e" and "r0" are dimensional numbers.
   * The dimensionless force (scalar part) is then given by
   *  hat(F) = 12 e (r0/r)^7 ((r0/r)^6 - 1)
   * where e is in the form of flag_LJ == 1.
   */
  int flag_LJ;
  double e;
  double r0;

  // derived parameters
  double theta; // angle of the cylinder opening from the sphere center
  double theta2; // angle of the cylinder opening on the right cavity

  // for dumbbell
  double center1; // x component of the center of left cavity
  double center2; // x component of the center of right cavity

  // for hex2d
  double Lx; // = 0.5 * L;
  double Ly; // = sqrt (3.0) * Lx;

  // for porous
  //         Lx = L
  //         Ly = L * sqrt (3.0)
  double Lz; // = L * sqrt (6.0) * 2.0 / 3.0;
};


/* initialize struct confinement
 * INPUT
 *  type : confinement type parameter (0 - 4)
 *         0 == sphere
 *               R : radius of cavity centered at (0, 0, 0)
 *         1 == sphere + hole
 *               R : radius of cavity centered at (0, 0, 0)
 *               r : radius of the hole at (0, 0, 1) direction
 *         2 == cylinder
 *               r       : radius of the cylinder
 *               x, y, z : direction vector of the cylinder
 *               the cylinder center goes through (0, 0, 0) and (x, y, z).
 *         3 == dumbbell
 *               R  : left cavity radius centered at (center1, 0, 0)
 *               R2 : right cavity radius centered at (center2, 0, 0)
 *               L  : length of the cylinder
 *               r  : cylinder radius
 *               the origin is at the center of the cylinder
 *         4 == hex2d
 *               R : cavity radius
 *               r : cylinder radius
 *               L : lattice spacing
 *         5 == porous
 *               R : particle radius
 *               L : lattice spacing in x (2R for touching case)
 *  R       : cavity radius
 *  r       : cylinder radius
 *  x, y, z : cylinder direction (for type 2, cylinder)
 *  R2      : right cavity radius (for type 3, dumbbell)
 *  L       : lattice spacing (for type 4, hex2d)
 *  Lennard-Jones parameters
 *  flag_LJ : if 0, e = hat(e) = "e" / kT
 *            if 1, e = hat(e) / (peclet * hat(r0))
 *  e  : the dimensionless number defined by
 *         e = hat(e) for flag_LJ == 0, or
 *         e = hat(e) / (peclet * hat(r0)) for flag_LJ == 1,
 *       where
 *         hat(e)  = "e" / kT
 *         hat(r0) = "r0" / length
 *       note that "e" and "r0" are dimensional numbers.
 *  r0 : the dimensionless number defined by
 *         r0 = hat(r0) = "r0" / length
 *       note that "r0" is a dimensional number.
 */
struct confinement *
CF_init (int type,
	 double R,
	 double r,
	 double x, double y, double z,
	 double R2,
	 double L,
	 int flag_LJ,
	 double e,
	 double r0);

void
CF_free (struct confinement *cf);


/* set LJ parameters for run
 * INPUT
 *  cf     : struct confinement
 *  peclet : peclet number
 * OUTPUT
 *  cf->e  : defined as 
 *             cf->e = hat(e)/(peclet * hat(r0))
 *           where
 *             hat(e)  = "e" / kT
 *             hat(r0) = "r0" / length
 *             "e" and "r0" are dimensional numbers
 *           therefore, cf->flag-LJ is set by 1.
 */
void
CF_set (struct confinement *cf,
	double peclet);

/* wrapper for confinement forces
 * INPUT
 *  c          : struct confinement
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_calc_force (struct confinement *cf,
	       struct stokes *sys,
	       double *f,
	       int flag_add);

/* from confinement-guile.h */
/* get confinement from SCM
 * in SCM, confinement is given by one of these:
 * for spherical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *  ))
 * for spherical confinement with a hole,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere+hole"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *    1.0  ;; radius of the hole at (0, 0, 1) direction
 *  ))
 * for cylindrical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "cylinder"    ;; the cylinder center goes through (0, 0, 0) and (x, y, z).
 *    10.0          ;; radius of the cylinder
 *    1.0  0.0  0.0 ;; direction vector (x, y, z) of the cylinder
 *  ))
 * for dumbbell confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "dumbbell" ;; the origin is at the center of the cylinder
 *    10.0       ;; left cavity radius centered at (center1, 0, 0)
 *    10.0       ;; right cavity radius centered at (center2, 0, 0)
 *    2.0        ;; length of the cylinder
 *    1.0        ;; cylinder radius
 *  ))
 * for 2D hexagonal confinement with cylinder pipe,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "hex2d"
 *    10.0    ;; cavity radius
 *    1.0     ;; cylinder radius
 *    12.0    ;; lattice spacing
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "confinement".
 * OUTPUT
 *  returned value : struct confinement
 *                   if NULL is returned, it failed (not defined)
 */
struct confinement *
CF_guile_get (const char *var);


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
  int flag_lub_B; // for BD_calc_FB(), lub among mobile particles

  // auxiliary imposed-flow parameters for simple shear
  double t0; // reference time for s0
  double s0; // cell shift at time t0 (for shear_mode = 1 or 2)

  double st; // currently this is just place holders

  struct BeadRod *br;
  struct bonds *bonds;
  double gamma;
  struct EV *ev;
  struct angles *ang;
  struct EV_DH *ev_dh;
  struct EV_LJ *ev_LJ;
  struct confinement *cf;

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
 *        NOTE, flag_lub_B is used for BD_calc_FB() where
 *        lub among ONLY mobile particles are taken.
 *        therefore, check for the existance of mobile pair(s)
 *        not excluded by sys->ex_lub for flag_lub_B.
 *  (double) stokes -- currently this is just a place holder
 *  (struct BeadRod *)br
 *  (struct bonds *)bonds
 *  (double) gamma
 *  (struct EV *)ev
 *  (struct angles *)ang
 *  (struct EV_DH *)ev_dh
 *  (struct EV_LJ *)ev_LJ
 *  (struct confinement *)cf
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
		struct BeadRod *br,
		struct bonds *bonds,
		double gamma,
		struct EV *ev,
		struct angles *ang,
		struct EV_DH *ev_dh,
		struct EV_LJ *ev_LJ,
		struct confinement *cf,
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


/* from nitsol_c.h */
struct NITSOL {
  // size of the problem
  int n;

  // tolerance values
  double ftol;
  double stptol;

  // input parameters
  int input[10];

  // work area
  int lrwork;
  double *rwork;

  // parameters to the functions
  double *rpar;
  int *ipar;

  // functions
  void (*f)(int *n, double *xcur, double *fcur,
	    double *rpar, int *ipar, int *itrmf);
  void (*jacv)(int *n, double *xcur, double *fcur,
	       int *ijob, double *v, double *z,
	       double *rpar, int *ipar, int *itrmjv);
  double (*dinpr)(int* N, 
		  double* X, int* incX, 
		  double* Y, int* incY);
  double (*dnorm)(int* N, double* X, int* incX);

  // other parameters
  int iplvl;
  int ipunit;

  double choice1_exp;
  double choice2_exp;
  double choice2_coef;
  double eta_cutoff;
  double etamax;
  double thmin;
  double thmax;
  double etafixed;

  // output parameters
  int iterm;
  int info[6];
};


struct NITSOL *
NITSOL_init (void);

void
NITSOL_free (struct NITSOL *nit);

void
NITSOL_set_n (struct NITSOL *nit,
	      int n);

void
NITSOL_set_GMRES (struct NITSOL *nit,
		  int restart);

void
NITSOL_set_BiCGSTAB (struct NITSOL *nit);

void
NITSOL_set_TFQMR (struct NITSOL *nit);

/*
 * INPUT
 *  p_flag  : 0 == no preconditioning
 *            1 == right preconditioning
 *  j_flag  : 0 == approximated for J*v
 *            1 == analytical J*v
 *  j_order : 1, 2, or 4 for the approximation of J*v
 *            (ignored for j_flag == 1)
 *  jacv  : jacobian and preconditioning function
 *          (ignored for p_flag == 0 and j_flag == 0)
 */
void
NITSOL_set_jacv (struct NITSOL *nit,
		 int p_flag,
		 int j_flag,
		 int j_order,
		 void (*jacv)(int *n, double *xcur, double *fcur,
			      int *ijob, double *v, double *z,
			      double *rpar, int *ipar, int *itrmjv));

void
NITSOL_set_tol (struct NITSOL *nit,
		double ftol, double stptol);

/*
 * i do not know the way of specifying parameters is fine (for choices 1, 2, 3)
 * INPUT
 *  flag : choice 0
 *         choice 1, p1 is for exp (alpha)
 *         choice 2, p1 and p2 are for exp (alpha) and coef (gamma)
 *         choice 3, p1 is for etafixed
 *  p1, p2 : parameters
 */
void
NITSOL_set_forcing (struct NITSOL *nit,
		    int flag,
		    double p1, double p2);

/* set iplvl and ipunit
 * INPUT
 *  iplvl = 0 => no printout
 *        = 1 => iteration numbers and F-norms
 *        = 2 => ... + some stats, step norms, and linear model norms
 *        = 3 => ... + some Krylov solver and backtrack information
 *        = 4 => ... + more Krylov solver and backtrack information
 *  ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
 *           NOTE: If ipunit = 0 on input, then it is set to 6 below.
 */
void
NITSOL_set_iplvl (struct NITSOL *nit,
		  int iplvl, int ipunit);

/* set routines for calculating norm by BLAS routines (ddot_ and dnrm2_)
 */
void
NITSOL_set_norm_by_BLAS (struct NITSOL *nit);


/**
 * parsing functions for the information about NITSOL
 */
void
NITSOL_parse_input (struct NITSOL *nit,
		    FILE *out);

void
NITSOL_parse_info (struct NITSOL *nit,
		   FILE *out);

void
NITSOL_parse_nitinfo (FILE *out);

void
NITSOL_parse_iterm (struct NITSOL *nit,
		    FILE *out);


/**
 * solver routine via struct NITSOL
 */
void
NITSOL_solve (struct NITSOL *nit,
	      double *x);

/* from bd-imp.h */
#include <gsl/gsl_multiroots.h>
struct BD_imp {
  struct BD_params *BD;
  double dt;
  double *x0;
  double *q0;
  double fact;
  double *z;

  int flag_solver; /* 0 == GSL
		    * 1 == NITSOL
		    */

  // GSL stuff
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_function *F;
  gsl_multiroot_fsolver  *S;
  gsl_vector *guess;

  // NITSOL stuff
  struct NITSOL *nit;

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
 *  flag_solver : 0 == GSL solver
 *                1 == NITSOL
 *  itmax : max of iteration for the root-finding
 *  eps   : tolerance for the root-finding
 */
struct BD_imp *
BD_imp_init (struct BD_params *BD,
	     int flag_solver,
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
