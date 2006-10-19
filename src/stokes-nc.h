/* header file for stokes-nc.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.h,v 5.3 2006/10/19 02:44:30 ichiki Exp $
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

  int xf_id;
  int uf_id;
  int of_id;
  int ef_id;
  int ff_id;
  int tf_id;
  int sf_id;

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

  int flag_xf;
  int flag_uf;
  int flag_of;
  int flag_ef;
  int flag_ff;
  int flag_tf;
  int flag_sf;
};


void
stokes_nc_error (int status, const char * message, const char * varname);

void
stokes_nc_print_actives (struct stokes_nc * nc,
			 FILE * out);


/* initialize NetCDF file for libstokes for mob_F problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, x, u.
 */
struct stokes_nc *
stokes_nc_mob_f_init (const char * filename, int np);
/* initialize NetCDF file for libstokes for mob_FT problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, t0, x, u, o.
 */
struct stokes_nc *
stokes_nc_mob_ft_init (const char * filename, int np);
/* initialize NetCDF file for libstokes for mob_FTS problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, t0, e0, x, u, o, s.
 */
struct stokes_nc *
stokes_nc_mob_fts_init (const char * filename, int np);

/* initialize NetCDF file for libstokes for mob_fix_F problem
 * with a constant imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui0, xf0, f0, uf0, x, u, ff.
 */
struct stokes_nc *
stokes_nc_mob_fix_f_i0_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FT problem
 * with a constant imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui0, oi0, xf0, f0, t0, uf0, of0, x, u, o, ff, tf.
 */
struct stokes_nc *
stokes_nc_mob_fix_ft_i0_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FTS problem
 * with a constant imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui0, oi0, ei0,
 *                         xf0, f0, t0, e0, uf0, of0, ef0,
 *                         x, u, o, s, ff, tf, sf.
 */
struct stokes_nc *
stokes_nc_mob_fix_fts_i0_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_F problem
 * with a time-changing imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui, xf0, f0, uf0, x, u, ff.
 */
struct stokes_nc *
stokes_nc_mob_fix_f_it_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FT problem
 * with a time-changing imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui, oi, xf0, f0, t0, uf0, of0, x, u, o, ff, tf.
 */
struct stokes_nc *
stokes_nc_mob_fix_ft_it_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FTS problem
 * with a time-changing imposed flow
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, ui, oi, ei,
 *                         xf0, f0, t0, e0, uf0, of0, ef0,
 *                         x, u, o, s, ff, tf, sf.
 */
struct stokes_nc *
stokes_nc_mob_fix_fts_it_init (const char * filename, int nm, int nf);


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
/* set ui
 */
void
stokes_nc_set_ui (struct stokes_nc * nc,
		  const double * ui);
/* set oi
 */
void
stokes_nc_set_oi (struct stokes_nc * nc,
		  const double * oi);
/* set ei
 */
void
stokes_nc_set_ei (struct stokes_nc * nc,
		  const double * ei);
/* set x0
 */
void
stokes_nc_set_x0 (struct stokes_nc * nc,
		  const double * x0);
/* set u0
 */
void
stokes_nc_set_u0 (struct stokes_nc * nc,
		  const double * u0);
/* set o0 data
 */
void
stokes_nc_set_o0 (struct stokes_nc * nc,
		  const double * o0);
/* set e0 data
 */
void
stokes_nc_set_e0 (struct stokes_nc * nc,
		  const double * e0);
/* set f0 data
 */
void
stokes_nc_set_f0 (struct stokes_nc * nc,
		  const double * f0);
/* set t0 data
 */
void
stokes_nc_set_t0 (struct stokes_nc * nc,
		  const double * t0);
/* set s0 data
 */
void
stokes_nc_set_s0 (struct stokes_nc * nc,
		  const double * s0);
/* set xf0
 */
void
stokes_nc_set_xf0 (struct stokes_nc * nc,
		   const double * xf0);
/* set uf0
 */
void
stokes_nc_set_uf0 (struct stokes_nc * nc,
		   const double * uf0);
/* set of0 data
 */
void
stokes_nc_set_of0 (struct stokes_nc * nc,
		   const double * of0);
/* set ef0 data
 */
void
stokes_nc_set_ef0 (struct stokes_nc * nc,
		   const double * ef0);
/* set ff0 data
 */
void
stokes_nc_set_ff0 (struct stokes_nc * nc,
		   const double * ff0);
/* set tf0 data
 */
void
stokes_nc_set_tf0 (struct stokes_nc * nc,
		   const double * tf0);
/* set sf0 data
 */
void
stokes_nc_set_sf0 (struct stokes_nc * nc,
		   const double * sf0);

/* set time (step)
 */
void
stokes_nc_set_time (struct stokes_nc * nc,
		    int step, double time);
/* set x at time (step)
 */
void
stokes_nc_set_x (struct stokes_nc * nc,
		 int step, double time,
		 const double * x);
/* set u at time (step)
 */
void
stokes_nc_set_u (struct stokes_nc * nc,
		 int step, double time,
		 const double * u);
/* set o at time (step)
 */
void
stokes_nc_set_o (struct stokes_nc * nc,
		 int step, double time,
		 const double * o);
/* set e at time (step)
 */
void
stokes_nc_set_e (struct stokes_nc * nc,
		 int step, double time,
		 const double * e);
/* set f at time (step)
 */
void
stokes_nc_set_f (struct stokes_nc * nc,
		 int step, double time,
		 const double * f);
/* set t at time (step)
 */
void
stokes_nc_set_t (struct stokes_nc * nc,
		 int step, double time,
		 const double * t);
/* set s at time (step)
 */
void
stokes_nc_set_s (struct stokes_nc * nc,
		 int step, double time,
		 const double * s);;
/* set xf at time (step)
 */
void
stokes_nc_set_xf (struct stokes_nc * nc,
		  int step, double time,
		  const double * xf);
/* set uf at time (step)
 */
void
stokes_nc_set_uf (struct stokes_nc * nc,
		  int step, double time,
		  const double * uf);
/* set of at time (step)
 */
void
stokes_nc_set_of (struct stokes_nc * nc,
		  int step, double time,
		  const double * of);
/* set ef at time (step)
 */
void
stokes_nc_set_ef (struct stokes_nc * nc,
		  int step, double time,
		  const double * ef);
/* set ff at time (step)
 */
void
stokes_nc_set_ff (struct stokes_nc * nc,
		  int step, double time,
		  const double * ff);
/* set tf at time (step)
 */
void
stokes_nc_set_tf (struct stokes_nc * nc,
		  int step, double time,
		  const double * tf);
/* set sf at time (step)
 */
void
stokes_nc_set_sf (struct stokes_nc * nc,
		  int step, double time,
		  const double * sf);


#endif /* !_STOKES_NC_H_ */
