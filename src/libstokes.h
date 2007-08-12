/* header file for library 'libstokes'
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libstokes.h,v 1.36 2007/08/12 23:57:53 kichiki Exp $
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
 ***********************************/
struct stokes {
  int np;      /* number of all particles  */
  int nm;      /* number of mobile particles  */
  double *pos; /* position of particles  */
  double *a;   /* radius of particles
		* Note : NULL (default) is for monodisperse system  */

  int version; /* 0 = F, 1 = FT, 2 = FTS  */

  // parameters for the polydisperse system
  int twobody_nmax;// max order for the coefficient for twobody_scalars_res()
  int twobody_lub; // 0 (far form) or 1 (lub form) for twobody_scalars_res()
  int *poly_table; // for (i,j), [i*np+j] gives the index of "twobody_f_list"
  struct twobody_f_list *twobody_f_list;

  /* imposed flow */
  double Ui[3];
  double Oi[3];
  double Ei[5];

  int periodic; // 0 = non periodic, 1 = periodic

  /* for ewald codes */
  double rmax2;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double xi, xi2, xiaspi, xia2;
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
  double lubmin2;/* square of min distance for lub
		  * (distance is replaced by its root-square)
		  */
  double lubmax; /* max distance for lub;
		  * for the pair beyond this is just ignored.
		  * 0 means no limit for open systems and 
		  * all particles within +/-1 cells in x,y,z for the periodic
		  */

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

/* set np and nm and allocate the memory for pos[np*3]
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
stokes_set_iter (struct stokes * sys,
		 const char * solver,
		 int max,
		 int restart,
		 double eps,
		 int debug,
		 FILE * out);

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
		int flag_it);

/* close (and write if necessary) NetCDF file for libstokes
 */
void
stokes_nc_free (struct stokes_nc * nc);


/** set nc data **/
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
 */
struct stokes_nc *
stokes_nc_set_by_params (const char *out_file,
			 const struct stokes *sys,
			 int flag_Q,
			 const double *Ui, const double *Oi, const double *Ei,
			 const double *F,  const double *T,  const double *E,
			 const double *uf, const double *of, const double *ef,
			 const double *xf,
			 const double *lat);

/* check stokes_nc data with the parameters
 * INPUT
 *  nc         :
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
  int ntypes;                // number of bond-types
  double *k;                 // spring constant for the bond-type
  double *r0;                // equilibrium distance for the bond-type
  struct bond_pairs **pairs; // pairs for the bond-type
};


struct bond_pairs *
bond_pairs_init (void);

void
bond_pairs_free (struct bond_pairs *pairs);

void
bond_pairs_add (struct bond_pairs *pairs,
		int ia, int ib);


struct bonds *
bonds_init (void);

void
bonds_free (struct bonds *bonds);

void
bonds_add_type (struct bonds *bonds,
		double k, double r0);

/*
 * INPUT
 *  b          : struct bond
 *  sys        : struct stokes
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
fprint_bonds (FILE *out, struct bonds *bonds);


/* from bonds-guile.h */
/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     1.0 ; spring const
 *     2.1 ; natural distance
 *     ((0 1) ; list of pairs
 *      (1 2)
 *      (2 3)))
 *    (; bond 2
 *     1.0 ; spring const
 *     2.5 ; natural distance
 *     ((4 5) ; list of pairs
 *      (5 6)
 *      (6 7)))
 *   ))
 * OUTPUT
 *  returned value : struct bonds
 *                   if NULL is returned, it failed (not defined)
 */
struct bonds *
guile_get_bonds (const char * var);



/************************************
 ** routines for ODE integrator    **
 ************************************/
/* from ode.h */
struct ode_params
{
  struct stokes *sys;
  double *pos_fixed;
  double *F;
  double *T;
  double *E;
  double *uf;
  double *of;
  double *ef;
  int flag_mat;
  int flag_lub;
  double st;
  struct bonds *bonds;
  double gamma;
};


/* set the parameters to struct ode_params
 * INPUT
 *  (struct stokes *)sys
 *  (double *)pos_fix (the position of fixed particles)
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
		 double *pos_fixed,
		 double *F,
		 double *T,
		 double *E,
		 double *uf,
		 double *of,
		 double *ef,
		 int flag_lub,
		 int flag_mat,
		 double st,
		 struct bonds *bonds,
		 double gamma);


/* calc dydt for gsl_odeiv ONLY with bond interactions for relaxation
 * this is equivalent to dydt_relax_bond with a constant friction, where
 *  dx/dt = U = f_bond / gamma
 * INPUT
 *  t   : current time
 *  y[] : position of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->pos_fixed : 
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
 *           ode->pos_fixed : 
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
 *           ode->pos_fixed : 
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
 *           ode->pos_fixed : 
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
 *           ode->pos_fixed : 
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
