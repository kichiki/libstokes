/* header file for stokes-nc.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.h,v 5.1 2006/10/11 20:28:20 ichiki Exp $
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

/* initialize NetCDF file for libstokes for mob_F problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_f_init (const char * filename, int np);
/* initialize NetCDF file for libstokes for mob_FT problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_ft_init (const char * filename, int np);
/* initialize NetCDF file for libstokes for mob_FTS problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fts_init (const char * filename, int np);
/* initialize NetCDF file for libstokes for mob_fix_F problem
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_f_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FT problem
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_ft_init (const char * filename, int nm, int nf);
/* initialize NetCDF file for libstokes for mob_fix_FTS problem
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_fts_init (const char * filename, int nm, int nf);


/* close (and write if necessary) NetCDF file for libstokes
 */
void
stokes_nc_free (struct stokes_nc * nc);


void
stokes_nc_append (struct stokes_nc * nc,
		  int step, double time,
		  const double * x,
		  const double * u, const double * o, const double * e,
		  const double * f, const double * t, const double * s);



#endif /* !_STOKES_NC_H_ */
