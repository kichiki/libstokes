/* header file for stokes-nc.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.h,v 5.15 2008/10/08 03:25:57 kichiki Exp $
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
#ifndef	_STOKES_NC_H_
#define	_STOKES_NC_H_

#include "stokes.h" // struct stokes

#include "KIrand.h" // struct KIrand

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


/* return stokes version (F/FT/FTS)
 * INPUT
 *  nc : struct stokes_nc *nc
 * OUTPUT
 *  version : returned value
 */
int
stokes_nc_get_version (const struct stokes_nc *nc);


/* get parameters from stokes_nc data
 * for stokes_nc_set_by_params()
 * therefore, struct stokes *sys is incomplete.
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
 *  shear_mode : 0 == imposed flow is given by Ui,Oi,Ei.
 *               1 == x = flow dir, y = grad dir
 *               2 == x = flow dir, z = grad dir
 *  shear_rate : defined only for shear_mode = 1 or 2.
 *  flag_rng   :
 */
void
stokes_nc_get_params (const struct stokes_nc *nc,
		      struct stokes *sys,
		      int *flag_Q,
		      double *Ui, double *Oi, double *Ei,
		      double *F, double *T, double *E,
		      double *uf, double *of, double *ef,
		      double *xf,
		      double *lat,
		      int *shear_mode, double *shear_rate,
		      int *flag_rng);


#endif /* !_STOKES_NC_H_ */
