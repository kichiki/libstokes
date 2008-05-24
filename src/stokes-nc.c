/* NetCDF interface for libstokes
 * Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.c,v 5.16 2008/05/24 06:00:28 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // fabs()
#include "memory-check.h" // macro CHECK_MALLOC

#include <netcdf.h>

#include "stokes-nc-read.h"
#include "stokes-nc.h"


void
stokes_nc_error (int status, const char *message, const char *varname)
{
  if (varname != NULL)
    {
      fprintf (stderr, "NetCDF error %d %s: %s for %s\n",
	       status, nc_strerror (status),
	       message, varname);
    }
  else
    {
      fprintf (stderr, "NetCDF error %d %s: %s\n",
	       status, nc_strerror (status),
	       message);
    }
  exit (0);
}


void
stokes_nc_print_actives (struct stokes_nc *nc,
			 FILE *out)
{
  fprintf (out, "    | x | u o e | f t s | ui oi ei\n");
  fprintf (out, "----+---+-------+-------+----------\n");
  fprintf (out, "x0  | %d | %d %d %d | %d %d %d | %d  %d  %d\n",
	   nc->flag_x0,
	   nc->flag_u0,
	   nc->flag_o0,
	   nc->flag_e0,
	   nc->flag_f0,
	   nc->flag_t0,
	   nc->flag_s0,
	   nc->flag_ui0,
	   nc->flag_oi0,
	   nc->flag_ei0);
  fprintf (out, "xf0 | %d | %d %d %d | %d %d %d |\n",
	   nc->flag_xf0,
	   nc->flag_uf0,
	   nc->flag_of0,
	   nc->flag_ef0,
	   nc->flag_ff0,
	   nc->flag_tf0,
	   nc->flag_sf0);
  fprintf (out, "----+---+-------+-------+----------\n");
  fprintf (out, "x   | %d | %d %d %d | %d %d %d | %d  %d  %d\n",
	   nc->flag_x,
	   nc->flag_u,
	   nc->flag_o,
	   nc->flag_e,
	   nc->flag_f,
	   nc->flag_t,
	   nc->flag_s,
	   nc->flag_ui,
	   nc->flag_oi,
	   nc->flag_ei);
  fprintf (out, "xf  | %d | %d %d %d | %d %d %d |\n",
	   nc->flag_xf,
	   nc->flag_uf,
	   nc->flag_of,
	   nc->flag_ef,
	   nc->flag_ff,
	   nc->flag_tf,
	   nc->flag_sf);
  fprintf (out, "----+---+-------+-------+----------\n");
  fprintf (out, "q   | %d |\n",
	   nc->flag_q);
  fprintf (out, "a   | %d |\n",
	   nc->flag_a);
  fprintf (out, "af  | %d |\n",
	   nc->flag_af);

  // (int)shear_mode
  int shear_mode;
  int status;
  int index = 0;
  status = nc_get_var1_int (nc->id, nc->shear_mode_id, &index,
			    &shear_mode);
  if (status != NC_NOERR)
    {
      fprintf (stderr,
	       "# stokes_nc_print_actives() : nc_get_var1_int() failed"
	       " for shear_mode\n");
    }
  fprintf (out, "----+---+-------+-------+----------\n");
  fprintf (out, "shear_mode      | %d \n", shear_mode);
}


/* init a variable x[time][p][vec]
 */
static void
stokes_nc_init_t_p_vec (struct stokes_nc *nc,
			const char *name,
			const char *longname,
			int *var_id)
{
  int status;
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->p_dim;
  dimids[2] = nc->vec_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_p_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_p_vec", name);
    }
}
/* init a variable x[time][p][stt]
 * (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_t_p_stt (struct stokes_nc *nc,
			const char *name,
			const char *longname,
			int *var_id)
{
  int status;
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->p_dim;
  dimids[2] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_p_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_p_stt", name);
    }
}

/* init a variable x[time][p][quat]
 */
static void
stokes_nc_init_t_p_quat (struct stokes_nc *nc,
			 const char *name,
			 const char *longname,
			 int *var_id)
{
  int status;
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->p_dim;
  dimids[2] = nc->quat_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_p_quat", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_p_quat", name);
    }
}

/* init a variable x[time][pf][vec]
 */
static void
stokes_nc_init_t_pf_vec (struct stokes_nc *nc,
			 const char *name,
			 const char *longname,
			 int *var_id)
{
  int status;
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->pf_dim;
  dimids[2] = nc->vec_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_pf_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_pf_vec", name);
    }
}
/* init a variable x[time][pf][stt]
 * (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_t_pf_stt (struct stokes_nc *nc,
			 const char *name,
			 const char *longname,
			 int *var_id)
{
  int status;
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->pf_dim;
  dimids[2] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_pf_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_pf_stt", name);
    }
}

/* init a variable x[p][vec]
 */
static void
stokes_nc_init_p_vec (struct stokes_nc *nc,
		      const char *name,
		      const char *longname,
		      int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->p_dim;
  dimids[1] = nc->vec_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_p_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_p_vec", name);
    }
}
/* init a variable x[p][stt]
 * (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_p_stt (struct stokes_nc *nc,
		      const char *name,
		      const char *longname,
		      int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->p_dim;
  dimids[1] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_p_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_p_stt", name);
    }
}

/* init a variable x[pf][vec]
 */
static void
stokes_nc_init_pf_vec (struct stokes_nc *nc,
		       const char *name,
		       const char *longname,
		       int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->pf_dim;
  dimids[1] = nc->vec_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_pf_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_pf_vec", name);
    }
}
/* init a variable x[pf][stt]
 * (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_pf_stt (struct stokes_nc *nc,
		       const char *name,
		       const char *longname,
		       int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->pf_dim;
  dimids[1] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_pf_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_pf_stt", name);
    }
}

/* init a variable x[vec]
 */
static void
stokes_nc_init_vec (struct stokes_nc *nc,
		    const char *name,
		    const char *longname,
		    int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 1, &nc->vec_dim, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_vec", name);
    }
}
/* init a variable x[stt]
 * (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_stt (struct stokes_nc *nc,
		    const char *name,
		    const char *longname,
		    int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 1, &nc->stt_dim, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_stt", name);
    }
}
/* init a variable x[p]
 */
static void
stokes_nc_init_p (struct stokes_nc *nc,
		  const char *name,
		  const char *longname,
		  int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 1, &nc->p_dim, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_p", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_p", name);
    }
}
/* init a variable x[pf]
 */
static void
stokes_nc_init_pf (struct stokes_nc *nc,
		   const char *name,
		   const char *longname,
		   int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 1, &nc->pf_dim, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_pf", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_pf", name);
    }
}

/* init a variable x[time][vec]
 */
static void
stokes_nc_init_t_vec (struct stokes_nc *nc,
		      const char *name,
		      const char *longname,
		      int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->vec_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_vec", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_vec", name);
    }
}
/* init a variable x[time][stt] 
* (stt = symmetric traceless tensor)
 */
static void
stokes_nc_init_t_stt (struct stokes_nc *nc,
		      const char *name,
		      const char *longname,
		      int *var_id)
{
  int status;
  int dimids[2];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 2, dimids, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t_stt", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t_stt", name);
    }
}
/* init a variable x[time]
 */
static void
stokes_nc_init_t (struct stokes_nc *nc,
		  const char *name,
		  const char *longname,
		  int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 1, &nc->time_dim, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_t", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_t", name);
    }
}
/* init a variable x (just a scalar)
 */
static void
stokes_nc_init_scalar (struct stokes_nc *nc,
		       const char *name,
		       const char *longname,
		       int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 0, NULL, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_scalar", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_scalar", name);
    }
}
/* init a variable x (just a scalar) with NC_INT
 */
static void
stokes_nc_init_scalar_int (struct stokes_nc *nc,
			   const char *name,
			   const char *longname,
			   int *var_id)
{
  int status;

  status = nc_def_var
    (nc->id, name, NC_INT, 0, NULL, var_id);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_def_var() in stokes_nc_init_scalar_int", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_scalar_int", name);
    }
}

/* initialize NetCDF file for libstokes (internal use only)
 * INPUT
 *  nm : number of MOBILE particles
 *       (note that np in NetCDF is only MOBILE particles)
 *  nf : number of fixed particles
 *  flag_* : active/inactive flags for each entry
 *  shear_mode : 0 == imposed flow is given by Ui,Oi,Ei.
 *               1 == x = flow dir, y = grad dir
 *               2 == x = flow dir, z = grad dir
 *               NOTE: shear_rate and shear_shift[t] are defined only for
 *               shear_mode = 1 or 2.
 * OUTPUT
 *  (returned value) : ncid
 */
static struct stokes_nc *
stokes_nc_init_ (const char *filename, int nm, int nf,
		 int flag_ui0,
		 int flag_oi0,
		 int flag_ei0,
		 int flag_ui,
		 int flag_oi,
		 int flag_ei,
		 int flag_x0,
		 int flag_u0,
		 int flag_o0,
		 int flag_e0,
		 int flag_f0,
		 int flag_t0,
		 int flag_s0,
		 int flag_xf0,
		 int flag_uf0,
		 int flag_of0,
		 int flag_ef0,
		 int flag_ff0,
		 int flag_tf0,
		 int flag_sf0,
		 int flag_x,
		 int flag_u,
		 int flag_o,
		 int flag_e,
		 int flag_f,
		 int flag_t,
		 int flag_s,
		 int flag_q,
		 int flag_xf,
		 int flag_uf,
		 int flag_of,
		 int flag_ef,
		 int flag_ff,
		 int flag_tf,
		 int flag_sf,
		 int flag_a,
		 int flag_af,
		 int shear_mode)
{
  struct stokes_nc *nc
    = (struct stokes_nc *)malloc (sizeof (struct stokes_nc));
  CHECK_MALLOC (nc, "stokes_nc_init");

  nc->flag_ui0 = flag_ui0;
  nc->flag_oi0 = flag_oi0;
  nc->flag_ei0 = flag_ei0;
  nc->flag_ui = flag_ui;
  nc->flag_oi = flag_oi;
  nc->flag_ei = flag_ei;
  nc->flag_x0 = flag_x0;
  nc->flag_u0 = flag_u0;
  nc->flag_o0 = flag_o0;
  nc->flag_e0 = flag_e0;
  nc->flag_f0 = flag_f0;
  nc->flag_t0 = flag_t0;
  nc->flag_s0 = flag_s0;
  nc->flag_xf0 = flag_xf0;
  nc->flag_uf0 = flag_uf0;
  nc->flag_of0 = flag_of0;
  nc->flag_ef0 = flag_ef0;
  nc->flag_ff0 = flag_ff0;
  nc->flag_tf0 = flag_tf0;
  nc->flag_sf0 = flag_sf0;
  nc->flag_x = flag_x;
  nc->flag_u = flag_u;
  nc->flag_o = flag_o;
  nc->flag_e = flag_e;
  nc->flag_f = flag_f;
  nc->flag_t = flag_t;
  nc->flag_s = flag_s;
  nc->flag_q = flag_q;
  nc->flag_xf = flag_xf;
  nc->flag_uf = flag_uf;
  nc->flag_of = flag_of;
  nc->flag_ef = flag_ef;
  nc->flag_ff = flag_ff;
  nc->flag_tf = flag_tf;
  nc->flag_sf = flag_sf;
  nc->flag_a = flag_a;
  nc->flag_af = flag_af;

  /* create netCDF dataset: enter define mode */
  int status;
  status = nc_create (filename,
		      NC_NOCLOBBER|NC_64BIT_OFFSET,
		      &(nc->id));
  // NC_NOCLOBBER : not clobber (overwrite) an existing dataset
  // NC_64BIT_OFFSET : for very large (i.e. over 2 GB) data files
  if (status != NC_NOERR)
    {
      stokes_nc_error (status, "at nc_create() in stokes_nc_init", NULL);
    }

  /* define dimensions: from name and length */
  nc->np  = nm;
  nc->npf = nf;
  nc->nvec = 3; // number of vector components
  nc->nstt = 5; // number of symmetric traceless tensor components
  nc->nquat= 4; // number of components for quaternion

  /* particle index */
  status = nc_def_dim (nc->id, "p", nc->np, &(nc->p_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for p in stokes_nc_init", NULL);
    }
  int dimids[3];
  dimids[0] = nc->p_dim;
  status = nc_def_var (nc->id, "p", NC_INT, 1, dimids, &(nc->p_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_var() for p in stokes_nc_init", NULL);
    }

  /* fixed particle index */
  if (nf > 0)
    {
      status = nc_def_dim (nc->id, "pf", nc->npf, &(nc->pf_dim));
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
 			   "at nc_def_dim() for pf in stokes_nc_init", NULL);
	}
      dimids[0] = nc->pf_dim;
      status = nc_def_var (nc->id, "pf", NC_INT, 1, dimids, &(nc->pf_id));
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_def_dim() for pf in stokes_nc_init", NULL);
	}
    }

  /* vector index */
  status = nc_def_dim (nc->id, "vec", nc->nvec, &(nc->vec_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for vec in stokes_nc_init", NULL);
    }
  dimids[0] = nc->vec_dim;
  status = nc_def_var (nc->id, "vec", NC_INT, 1, dimids, &(nc->vec_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for vec in stokes_nc_init", NULL);
    }

  /* symmetric traceless tensor index */
  status = nc_def_dim (nc->id, "stt", nc->nstt, &(nc->stt_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for stt in stokes_nc_init", NULL);
    }
  dimids[0] = nc->stt_dim;
  status = nc_def_var (nc->id, "stt", NC_INT, 1, dimids, &(nc->stt_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for stt in stokes_nc_init", NULL);
    }

  /* quaternion index */
  status = nc_def_dim (nc->id, "quat", nc->nquat, &(nc->quat_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for quat in stokes_nc_init", NULL);
    }
  dimids[0] = nc->quat_dim;
  status = nc_def_var (nc->id, "quat", NC_INT, 1, dimids, &(nc->quat_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for quat in stokes_nc_init", NULL);
    }

  /* time index */
  status = nc_def_dim(nc->id, "time", NC_UNLIMITED, &(nc->time_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for time in stokes_nc_init", NULL);
    }
  dimids[0] = nc->time_dim;
  status = nc_def_var (nc->id, "time", NC_DOUBLE, 1, dimids, &(nc->time_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for time in stokes_nc_init", NULL);
    }
 
  /* lattice vector */
  stokes_nc_init_vec (nc, "l",
		      "lattice vector",
		      &nc->l_id);


  /* data of imposed flow */
  /* shear_mode */
  stokes_nc_init_scalar_int (nc, "shear_mode",
			     "shear mode",
			     &nc->shear_mode_id);
  if (shear_mode != 0)
    {
      /* shear_rate */
      stokes_nc_init_scalar (nc, "shear_rate",
			     "shear rate",
			     &nc->shear_rate_id);
      /* shift[t] */
      stokes_nc_init_t (nc, "shear_shift",
			"cell shift for shear boundary condition",
			&nc->shear_shift_id);
    }
  /* Ui0 */
  if (flag_ui0 != 0)
    {
      stokes_nc_init_vec (nc, "Ui0",
			  "constant imposed translational velocity",
			  &nc->ui0_id);
    }
  /* Oi0 */
  if (flag_oi0 != 0)
    {
      stokes_nc_init_vec (nc, "Oi0",
			  "constant imposed angular velocity",
			  &nc->oi0_id);
    }
  /* Ei0 */
  if (flag_ei0 != 0)
    {
      stokes_nc_init_stt (nc, "Ei0",
			  "constant imposed strain",
			  &nc->ei0_id);
    }
  /* Ui */
  if (flag_ui != 0)
    {
      stokes_nc_init_t_vec (nc, "Ui",
			    "(changing) imposed translational velocity",
			    &nc->ui_id);
    }
  /* Oi */
  if (flag_oi != 0)
    {
      stokes_nc_init_t_vec (nc, "Oi",
			    "(changing) imposed angular velocity",
			    &nc->oi_id);
    }
  /* Ei */
  if (flag_ei != 0)
    {
      stokes_nc_init_t_stt (nc, "Ei",
			    "(changing) imposed strain",
			    &nc->ei_id);
    }

  /* data for time-independent particles */
  /* x0 */
  if (flag_x0 != 0)
    {
      stokes_nc_init_p_vec (nc, "x0",
			    "positions of particle center at t0",
			    &nc->x0_id);
    }
  /* U0 */
  if (flag_u0 != 0)
    {
      stokes_nc_init_p_vec (nc, "U0",
			    "translational velocities at t0",
			    &nc->u0_id);
    }
  /* O0 */
  if (flag_o0 != 0)
    {
      stokes_nc_init_p_vec (nc, "O0",
			    "angular velocities at t0",
			    &nc->o0_id);
    }
  /* E0 */
  if (flag_e0 != 0)
    {
      stokes_nc_init_p_stt (nc, "E0",
			    "rate of strain at t0",
			    &nc->e0_id);
    }
  /* F0 */
  if (flag_f0 != 0)
    {
      stokes_nc_init_p_vec (nc, "F0",
			    "forces at t0",
			    &nc->f0_id);
    }
  /* T0 */
  if (flag_t0 != 0)
    {
      stokes_nc_init_p_vec (nc, "T0",
			    "torques at t0",
			    &nc->t0_id);
    }
  /* S0 */
  if (flag_s0 != 0)
    {
      stokes_nc_init_p_stt (nc, "S0",
			    "stresslets at t0",
			    &nc->s0_id);
    }

  /* data for (mobile) particles */
  /* x */
  if (flag_x != 0)
    {
      stokes_nc_init_t_p_vec (nc, "x",
			      "positions of particle center",
			      &nc->x_id);
    }
  /* q */
  if (flag_q != 0)
    {
      stokes_nc_init_t_p_quat (nc, "q",
			       "quaternion of particle direction",
			       &nc->q_id);
    }
  /* a */
  if (flag_a != 0)
    {
      stokes_nc_init_p (nc, "a",
			"radii of particle",
			&nc->a_id);
    }
  /* U */
  if (flag_u != 0)
    {
      stokes_nc_init_t_p_vec (nc, "U",
			      "translational velocities",
			      &nc->u_id);
    }
  /* O */
  if (flag_o != 0)
    {
      stokes_nc_init_t_p_vec (nc, "O",
			      "angular velocities",
			      &nc->o_id);
    }
  /* E */
  if (flag_e != 0)
    {
      stokes_nc_init_t_p_stt (nc, "E",
			      "rate of strain",
			      &nc->e_id);
    }
  /* F */
  if (flag_f != 0)
    {
      stokes_nc_init_t_p_vec (nc, "F",
			      "forces",
			      &nc->f_id);
    }
  /* T */
  if (flag_t != 0)
    {
      stokes_nc_init_t_p_vec (nc, "T",
			      "torques",
			      &nc->t_id);
    }
  /* S */
  if (flag_s != 0)
    {
      stokes_nc_init_t_p_stt (nc, "S",
			      "stresslets",
			      &nc->s_id);
    }

  /* data for fixed particles */
  if (nf > 0)
    {
      /* data for time-independent fixed particles */
      /* xf0 */
      if (flag_xf0 != 0)
	{
	  stokes_nc_init_pf_vec
	    (nc, "xf0",
	     "positions of fixed particle center at t0",
	     &nc->xf0_id);
	}
      /* af */
      if (flag_af != 0)
	{
	  stokes_nc_init_pf (nc, "af",
			     "radii of particle",
			     &nc->af_id);
	}
      /* Uf0 */
      if (flag_uf0 != 0)
	{
	  stokes_nc_init_pf_vec
	    (nc, "Uf0",
	     "translational velocities for fixed particles at t0",
	     &nc->uf0_id);
	}
      /* Of0 */
      if (flag_of0 != 0)
	{
	  stokes_nc_init_pf_vec
	    (nc, "Of0",
	     "angular velocities for fixed particles at t0",
	     &nc->of0_id);
	}
      /* Ef0 */
      if (flag_ef0 != 0)
	{
	  stokes_nc_init_pf_stt
	    (nc, "Ef0",
	     "rate of strain for fixed particles at t0",
	     &nc->ef0_id);
	}
      /* Ff0 */
      if (flag_ff0 != 0)
	{
	  stokes_nc_init_pf_vec
	    (nc, "Ff0",
	     "forces for fixed particles at t0",
	     &nc->ff0_id);
	}
      /* Tf0 */
      if (flag_tf0 != 0)
	{
	  stokes_nc_init_pf_vec
	    (nc, "Tf0",
	     "torques for fixed particles at t0",
	     &nc->tf0_id);
	}
      /* Sf0 */
      if (flag_sf0 != 0)
	{
	  stokes_nc_init_pf_stt
	    (nc, "Sf0",
	     "stresslets for fixed particles at t0",
	     &nc->sf0_id);
	}
      /* end of data for time-independent fixed particles */

      /* xf */
      if (flag_xf != 0)
	{
	  stokes_nc_init_t_pf_vec (nc, "xf",
			      "positions of fixed particle center",
			      &nc->xf_id);
	}
      /* Uf */
      if (flag_uf != 0)
	{
	  stokes_nc_init_t_pf_vec (nc, "Uf",
			      "translational velocities for fixed particles",
			      &nc->uf_id);
	}
      /* Of */
      if (flag_of != 0)
	{
	  stokes_nc_init_t_pf_vec (nc, "Of",
			      "angular velocities for fixed particles",
			      &nc->of_id);
	}
      /* Ef */
      if (flag_ef != 0)
	{
	  stokes_nc_init_t_pf_stt (nc, "Ef",
			      "rate of strain for fixed particles",
			      &nc->ef_id);
	}
      /* Ff */
      if (flag_ff != 0)
	{
	  stokes_nc_init_t_pf_vec (nc, "Ff",
			      "forces for fixed particles",
			      &nc->ff_id);
	}
      /* Tf */
      if (flag_tf != 0)
	{
	  stokes_nc_init_t_pf_vec (nc, "Tf",
			      "torques for fixed particles",
			      &nc->tf_id);
	}
      /* S */
      if (flag_sf != 0)
	{
	  stokes_nc_init_t_pf_stt (nc, "Sf",
			      "stresslets for fixed particles",
			      &nc->sf_id);
	}
    }

  /* end definitions: leave define mode */
  status = nc_enddef (nc->id);  /*leave define mode*/
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_enddef() in stokes_nc_init", NULL);
    }


  // set (0, 0, 0) initially. this means non-periodic system
  double lattice [3] = {0.0, 0.0, 0.0};
  stokes_nc_set_l (nc, lattice);


  /* set dim data */
  int i;
  size_t idx[1];
  for (i = 0; i < nc->np; i ++)
    {
      idx[0] = i;
      status = nc_put_var1_int (nc->id, nc->p_id, idx, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1_int() for p in stokes_nc_init", NULL);
	}
    }
  for (i = 0; i < nc->nvec; i ++)
    {
      idx[0] = i;
      status = nc_put_var1_int (nc->id, nc->vec_id, idx, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1() for vec in stokes_nc_init", NULL);
	}
    }
  for (i = 0; i < nc->nstt; i ++)
    {
      idx[0] = i;
      status = nc_put_var1_int (nc->id, nc->stt_id, idx, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1() for stt in stokes_nc_init", NULL);
	}
    }
  for (i = 0; i < nc->nquat; i ++)
    {
      idx[0] = i;
      status = nc_put_var1_int (nc->id, nc->quat_id, idx, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1() for quat in stokes_nc_init", NULL);
	}
    }

  if (nf > 0)
    {
      for (i = 0; i < nc->npf; i ++)
	{
          idx[0] = i;
	  status = nc_put_var1_int (nc->id, nc->pf_id, idx, &i);
	  if (status != NC_NOERR)
	    {
	      stokes_nc_error
		(status,
		 "at nc_put_var1() for pf in stokes_nc_init", NULL);
	    }
	}
    }

  /* shear_mode */
  status = nc_put_var1_int (nc->id, nc->shear_mode_id, NULL, &shear_mode);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_var1() for shear_mode in stokes_nc_init", NULL);
    }

  return (nc);
}

/* initialize NetCDF file for libstokes for just x
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  the only activated entry is x.
 */
struct stokes_nc *
stokes_nc_x_init (const char * filename, int np)
{
  return stokes_nc_init_ (filename, np, 0,
			  0, // ui0
			  0, // oi0
			  0, // ei0
			  0, // ui
			  0, // oi
			  0, // ei
			  0, // x0
			  0, // u0
			  0, // o0
			  0, // e0
			  0, // f0
			  0, // t0
			  0, // s0
			  0, // xf0
			  0, // uf0
			  0, // of0
			  0, // ef0
			  0, // ff0
			  0, // tf0
			  0, // sf0
			  1, // x
			  0, // u
			  0, // o
			  0, // e
			  0, // f
			  0, // t
			  0, // s
			  0, // q
			  0, // xf
			  0, // uf
			  0, // of
			  0, // ef
			  0, // ff
			  0, // tf
			  0, // sf
			  0, // a
			  0, // af
			  0);// shear_mode
}

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
		int shear_mode)
{
  int flag_ui0 = 0;
  int flag_oi0 = 0;
  int flag_ei0 = 0;
  int flag_ui = 0;
  int flag_oi = 0;
  int flag_ei = 0;
  int flag_x0 = 0;
  int flag_u0 = 0;
  int flag_o0 = 0;
  int flag_e0 = 0;
  int flag_f0 = 0;
  int flag_t0 = 0;
  int flag_s0 = 0;
  int flag_xf0 = 0;
  int flag_uf0 = 0;
  int flag_of0 = 0;
  int flag_ef0 = 0;
  int flag_ff0 = 0;
  int flag_tf0 = 0;
  int flag_sf0 = 0;
  int flag_x = 0;
  int flag_u = 0;
  int flag_o = 0;
  int flag_e = 0;
  int flag_f = 0;
  int flag_t = 0;
  int flag_s = 0;
  int flag_xf = 0;
  int flag_uf = 0;
  int flag_of = 0;
  int flag_ef = 0;
  int flag_ff = 0;
  int flag_tf = 0;
  int flag_sf = 0;
  int flag_a = 0;
  int flag_af = 0;


  if (flag_it == 0)
    {
      flag_ui0 = 1;
      flag_oi0 = 1;
      flag_ei0 = 1;
    }
  else
    {
      flag_ui = 1;
      flag_oi = 1;
      flag_ei = 1;
    }

  if (nf == 0)
    {
      flag_a = flag_poly;
      if (version == 0) // F version
	{
	  flag_f0 = 1;
	  flag_x = 1;
	  flag_u = 1;
	}
      else if (version == 1) // FT version
	{
	  flag_f0 = 1;
	  flag_t0 = 1;
	  flag_x = 1;
	  flag_u = 1;
	  flag_o = 1;
	}
      else // FTS version
	{
	  flag_f0 = 1;
	  flag_t0 = 1;
	  flag_e0 = 1;
	  flag_x = 1;
	  flag_u = 1;
	  flag_o = 1;
	  flag_s = 1;
	}
    }
  else if (nf > 0)
    {
      flag_a = flag_poly;
      flag_af = flag_poly;
      if (version == 0) // F version
	{
	  flag_f0 = 1;
	  flag_x = 1;
	  flag_u = 1;

	  flag_uf0 = 1;
	  flag_xf0 = 1;
	  flag_ff = 1;
	}
      else if (version == 1) // FT version
	{
	  flag_f0 = 1;
	  flag_t0 = 1;
	  flag_x = 1;
	  flag_u = 1;
	  flag_o = 1;

	  flag_uf0 = 1;
	  flag_of0 = 1;
	  flag_xf0 = 1;
	  flag_ff = 1;
	  flag_tf = 1;
	}
      else // FTS version
	{
	  flag_f0 = 1;
	  flag_t0 = 1;
	  flag_e0 = 1;
	  flag_x = 1;
	  flag_u = 1;
	  flag_o = 1;
	  flag_s = 1;

	  flag_uf0 = 1;
	  flag_of0 = 1;
	  flag_ef0 = 1;
	  flag_xf0 = 1;
	  flag_ff = 1;
	  flag_tf = 1;
	  flag_sf = 1;
	}
    }

  return stokes_nc_init_ (filename, np, nf,
			  flag_ui0, flag_oi0, flag_ei0,
			  flag_ui,  flag_oi,  flag_ei,
			  flag_x0,
			  flag_u0,  flag_o0,  flag_e0,
			  flag_f0,  flag_t0,  flag_s0,
			  flag_xf0,
			  flag_uf0, flag_of0, flag_ef0,
			  flag_ff0, flag_tf0, flag_sf0,
			  flag_x,
			  flag_u,   flag_o,   flag_e,
			  flag_f,   flag_t,   flag_s,
			  flag_Q,
			  flag_xf,
			  flag_uf,  flag_of,  flag_ef,
			  flag_ff,  flag_tf,  flag_sf,
			  flag_a,   flag_af,
			  shear_mode);
}


/* close (and write if necessary) NetCDF file for libstokes
 */
void
stokes_nc_free (struct stokes_nc *nc)
{
  if (nc == NULL) return;

  /* close: save new netCDF dataset */
  int status;
  status = nc_close (nc->id);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_close() in stokes_nc_free", NULL);
    }

  free (nc);
}

/** set nc data **/
/* set shear_rate (just a scalar)
 */
void
stokes_nc_set_shear_rate (struct stokes_nc *nc,
			  double shear_rate)
{
  int status;
  status = nc_put_var1_double (nc->id, nc->shear_rate_id, NULL, &shear_rate);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_shear_rate", NULL);
    }
}
/* set l = (lx, ly, lz)
 */
void
stokes_nc_set_l (struct stokes_nc *nc,
		 const double *l)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->l_id, start, count, l);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_l", NULL);
    }
}
/* set ui0[vec]
 */
void
stokes_nc_set_ui0 (struct stokes_nc *nc,
		   const double * ui0)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->ui0_id, start, count, ui0);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_ui0", NULL);
    }
}
/* set oi[vec]
 */
void
stokes_nc_set_oi0 (struct stokes_nc *nc,
		   const double * oi0)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->oi0_id, start, count, oi0);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_oi0", NULL);
    }
}
/* set ei0[stt]
 */
void
stokes_nc_set_ei0 (struct stokes_nc *nc,
		   const double * ei0)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->ei0_id, start, count, ei0);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_ei", NULL);
    }
}

/* set x0[np][vec]
 */
void
stokes_nc_set_x0 (struct stokes_nc *nc,
		  const double * x0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->x0_id, start, count, x0);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_x0", NULL);
    }
}
/* set u0[np][vec]
 */
void
stokes_nc_set_u0 (struct stokes_nc *nc,
		  const double * u0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->u0_id, start, count, u0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_u0", NULL);
    }
}
/* set o0[np][vec]
 */
void
stokes_nc_set_o0 (struct stokes_nc *nc,
		  const double * o0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->o0_id, start, count, o0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_o0", NULL);
    }
}
/* set e0[np][stt]
 */
void
stokes_nc_set_e0 (struct stokes_nc *nc,
		  const double * e0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->e0_id, start, count, e0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_e0", NULL);
    }
}
/* set f0[np][vec]
 */
void
stokes_nc_set_f0 (struct stokes_nc *nc,
		  const double * f0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->f0_id, start, count, f0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_f0", NULL);
    }
}
/* set t0[np][vec]
 */
void
stokes_nc_set_t0 (struct stokes_nc *nc,
		  const double * t0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->t0_id, start, count, t0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_t0", NULL);
    }
}
/* set s0[np][stt]
 */
void
stokes_nc_set_s0 (struct stokes_nc *nc,
		  const double * s0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->np;
  count[1] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->s0_id, start, count, s0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_s0", NULL);
    }
}
/* set xf0[npf][vec]
 */
void
stokes_nc_set_xf0 (struct stokes_nc *nc,
		   const double * xf0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->xf0_id, start, count, xf0);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_xf0", NULL);
    }
}
/* set uf0[npf][vec]
 */
void
stokes_nc_set_uf0 (struct stokes_nc *nc,
		   const double * uf0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->uf0_id, start, count, uf0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_uf0", NULL);
    }
}
/* set of0[npf][vec]
 */
void
stokes_nc_set_of0 (struct stokes_nc *nc,
		   const double * of0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->of0_id, start, count, of0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_of0", NULL);
    }
}
/* set ef0[npf][stt]
 */
void
stokes_nc_set_ef0 (struct stokes_nc *nc,
		   const double * ef0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->ef0_id, start, count, ef0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_ef0", NULL);
    }
}
/* set ff0[npf][vec]
 */
void
stokes_nc_set_ff0 (struct stokes_nc *nc,
		   const double * ff0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->ff0_id, start, count, ff0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_ff0", NULL);
    }
}
/* set tf0[npf][vec]
 */
void
stokes_nc_set_tf0 (struct stokes_nc *nc,
		   const double * tf0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->tf0_id, start, count, tf0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_tf0", NULL);
    }
}
/* set sf0[npf][stt]
 */
void
stokes_nc_set_sf0 (struct stokes_nc *nc,
		   const double * sf0)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = 0;
  start[1] = 0;

  count[0] = nc->npf;
  count[1] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->sf0_id, start, count, sf0);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_sf0", NULL);
    }
}
/* set a[np]
 */
void
stokes_nc_set_a (struct stokes_nc *nc,
		 const double * a)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->np;

  status = nc_put_vara_double(nc->id, nc->a_id, start, count, a);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_a", NULL);
    }
}
/* set af[npf]
 */
void
stokes_nc_set_af (struct stokes_nc *nc,
		  const double * af)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->npf;

  status = nc_put_vara_double(nc->id, nc->af_id, start, count, af);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_af", NULL);
    }
}

/* set time (step)
 */
void
stokes_nc_set_time (struct stokes_nc *nc,
		    int step, double time)
{
  int status;
  size_t time_index[1]; // where to put the value
  time_index[0] = step;

  //status = nc_put_var1_double (nc->id, nc->time_id, (size_t *)&step, &time);
  status = nc_put_var1_double (nc->id, nc->time_id, time_index, &time);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for time in stokes_nc_set_time", NULL);
    }
}
/* set shear_shift[time]
 */
void
stokes_nc_set_shear_shift (struct stokes_nc *nc,
			   int step,
			   double shear_shift)
{
  size_t index[1]; // where to put the value
  index[0] = step;

  int status;
  status = nc_put_var1_double (nc->id, nc->shear_shift_id, index, &shear_shift);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_shear_shift", NULL);
    }
}

/* set ui[time][vec] at time (step)
 */
void
stokes_nc_set_ui (struct stokes_nc *nc,
		  int step,
		  const double * ui)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = step;
  start[1] = 0;

  count[0] = 1;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->ui_id, start, count, ui);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_ui", NULL);
    }
}
/* set oi[time][vec] at time (step)
 */
void
stokes_nc_set_oi (struct stokes_nc *nc,
		  int step,
		  const double * oi)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = step;
  start[1] = 0;

  count[0] = 1;
  count[1] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->oi_id, start, count, oi);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_oi", NULL);
    }
}
/* set ei[time][stt] at time (step)
 */
void
stokes_nc_set_ei (struct stokes_nc *nc,
		  int step,
		  const double * ei)
{
  size_t start[2];
  size_t count[2];

  int status;

  start[0] = step;
  start[1] = 0;

  count[0] = 1;
  count[1] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->ei_id, start, count, ei);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_ei", NULL);
    }
}
/* set x at time (step)
 */
void
stokes_nc_set_x (struct stokes_nc *nc,
		 int step,
		 const double * x)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->x_id, start, count, x);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_x", NULL);
    }
}
/* set q at time (step)
 */
void
stokes_nc_set_q (struct stokes_nc *nc,
		 int step,
		 const double * q)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nquat;

  status = nc_put_vara_double(nc->id, nc->q_id, start, count, q);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() in stokes_nc_set_q", NULL);
    }
}
/* set u at time (step)
 */
void
stokes_nc_set_u (struct stokes_nc *nc,
		 int step,
		 const double * u)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->u_id, start, count, u);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_u", NULL);
    }
}
/* set o at time (step)
 */
void
stokes_nc_set_o (struct stokes_nc *nc,
		 int step,
		 const double * o)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->o_id, start, count, o);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_o", NULL);
    }
}
/* set e at time (step)
 */
void
stokes_nc_set_e (struct stokes_nc *nc,
		 int step,
		 const double * e)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->np;
  count[2] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->e_id, start, count, e);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_e", NULL);
    }
}
/* set f at time (step)
 */
void
stokes_nc_set_f (struct stokes_nc *nc,
		 int step,
		 const double * f)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->f_id, start, count, f);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_f", NULL);
    }
}
/* set t at time (step)
 */
void
stokes_nc_set_t (struct stokes_nc *nc,
		 int step,
		 const double * t)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->t_id, start, count, t);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_t", NULL);
    }
}
/* set s at time (step)
 */
void
stokes_nc_set_s (struct stokes_nc *nc,
		 int step,
		 const double * s)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->np;
  count[2] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->s_id, start, count, s);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_s", NULL);
    }
}

/* set xf at time (step)
 */
void
stokes_nc_set_xf (struct stokes_nc *nc,
		  int step,
		  const double * xf)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->npf;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->xf_id, start, count, xf);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for xf in stokes_nc_set_xf", NULL);
    }
}
/* set uf at time (step)
 */
void
stokes_nc_set_uf (struct stokes_nc *nc,
		  int step,
		  const double * uf)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->npf;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->uf_id, start, count, uf);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_uf", NULL);
    }
}
/* set of at time (step)
 */
void
stokes_nc_set_of (struct stokes_nc *nc,
		  int step,
		  const double * of)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->npf;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->of_id, start, count, of);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_of", NULL);
    }
}
/* set ef at time (step)
 */
void
stokes_nc_set_ef (struct stokes_nc *nc,
		  int step,
		  const double * ef)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->npf;
  count[2] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->ef_id, start, count, ef);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_ef", NULL);
    }
}
/* set ff at time (step)
 */
void
stokes_nc_set_ff (struct stokes_nc *nc,
		  int step,
		  const double * ff)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->npf;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->ff_id, start, count, ff);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_ff", NULL);
    }
}
/* set tf at time (step)
 */
void
stokes_nc_set_tf (struct stokes_nc *nc,
		  int step,
		  const double * tf)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = nc->npf;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->tf_id, start, count, tf);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_tf", NULL);
    }
}
/* set sf at time (step)
 */
void
stokes_nc_set_sf (struct stokes_nc *nc,
		  int step,
		  const double * sf)
{
  size_t start[3];
  size_t count[3];

  int status;

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->npf;
  count[2] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->sf_id, start, count, sf);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_put_vara_double() in stokes_nc_set_sf", NULL);
    }
}


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
			 int shear_mode, double shear_rate)
{
  int nm = sys->nm;
  int nf = sys->np - nm;

  int flag_poly = 0;
  if (sys->a != NULL) flag_poly = 1;

  // create new stokes_nc
  struct stokes_nc *nc
    = stokes_nc_init (out_file, nm, nf,
		      sys->version, flag_poly, flag_Q, 0,
		      shear_mode);
  CHECK_MALLOC (nc, "stokes_nc_set_by_params");

  stokes_nc_set_ui0 (nc, Ui);
  stokes_nc_set_oi0 (nc, Oi);
  stokes_nc_set_ei0 (nc, Ei);
  // non-periodic system
  if (sys->periodic == 1)
    {
      stokes_nc_set_l (nc, lat);
    }

  // first for mobile particle informations (for both mob and mix)
  if (flag_poly != 0)
    {
      stokes_nc_set_a (nc, sys->a);
    }
  if (sys->version == 0)
    {
      // F version
      stokes_nc_set_f0 (nc, F);
    }
  else if (sys->version == 1)
    {
      // FT version
      stokes_nc_set_f0 (nc, F);
      stokes_nc_set_t0 (nc, T);
    }
  else
    {
      // FTS version
      stokes_nc_set_f0 (nc, F);
      stokes_nc_set_t0 (nc, T);
      stokes_nc_set_e0 (nc, E);
    }

  // mix problem (with fixed particles)
  if (nf > 0)
    {
      if (flag_poly != 0)
	{
	  stokes_nc_set_af (nc, sys->a + sys->nm);
	}

      if (sys->version == 0)
	{
	  // F version
	  stokes_nc_set_uf0 (nc, uf);
	}
      else if (sys->version == 1)
	{
	  // FT version
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	}
      else
	{
	  // FTS version
	  stokes_nc_set_uf0 (nc, uf);
	  stokes_nc_set_of0 (nc, of);
	  stokes_nc_set_ef0 (nc, ef);
	}
      stokes_nc_set_xf0 (nc, xf);
    }

  if (shear_mode != 0)
    {
      stokes_nc_set_shear_rate (nc, shear_rate);
    }

  return (nc);
}

/*
 * INPUT
 *  x[n] :
 *  y[n] :
 *  tiny : if (fabs (x[i] - y[i]) < tiny), x[i] == y[i].
 * OUTPUT
 *  0 : x[] == y[]
 *  1 : x[] != y[]
 */
static int
comp_array (const double *x, const double *y, int n, double tiny)
{
  int check = 0;
  int i;
  for (i = 0; i < n; i ++)
    {
      double d = fabs (x[i] - y[i]);
      if (d >= tiny)
	{
	  check = 1;
	  break;
	}
    }
  return (check);
}

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
			double tiny)
{
  int nm = sys->nm;
  int nf = sys->np - nm;

  int nm3 = nm * 3;
  int nm5 = nm * 5;
  int nf3 = nf * 3;
  int nf5 = nf * 5;

  // check the parameters
  int check = 0;
  if (nm != nc->np)
    {
      fprintf (stderr, "mismatch np %d != %d\n", nm, nc->np);
      check ++;
    }
  if (nf != nc->npf)
    {
      fprintf (stderr, "mismatch npf %d != %d\n", nf, nc->npf);
      check ++;
    }
  if (nf != nc->npf)
    {
      fprintf (stderr, "mismatch npf %d != %d\n", nf, nc->npf);
      check ++;
    }
  // because flag_it is always 1 now
  if (nc->flag_ui0 != 1 ||
      nc->flag_oi0 != 1 ||
      nc->flag_ei0 != 1)
    {
      fprintf (stderr, "imposed flows are not set (%d,%d,%d)\n",
	       nc->flag_ui0,
	       nc->flag_oi0,
	       nc->flag_ei0);
      check ++;
    }
  if (nf == 0)
    {
      if ((sys->a == NULL && nc->flag_a != 0) ||
	  (sys->a != NULL && nc->flag_a == 0))
	{
	  fprintf (stderr, "mismatch flag_a\n");
	  check ++;
	}
      if (sys->version == 0) // F version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1)
	    {
	      fprintf (stderr, "F params are not set properly.\n");
	      check ++;
	    }
	}
      else if (sys->version == 1) // FT version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_t0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1 ||
	      nc->flag_o != 1)
	    {
	      fprintf (stderr, "FT params are not set properly.\n");
	      check ++;
	    }
	}
      else // FTS version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_t0 != 1 ||
	      nc->flag_e0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1 ||
	      nc->flag_o != 1 ||
	      nc->flag_s != 1)
	    {
	      fprintf (stderr, "FTS params are not set properly.\n");
	      check ++;
	    }
	}
    }
  else if (nf > 0)
    {
      if ((sys->a == NULL && nc->flag_af != 0) ||
	  (sys->a != NULL && nc->flag_af == 0))
	{
	  fprintf (stderr, "mismatch flag_af\n");
	  check ++;
	}
      if (sys->version == 0) // F version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1 ||
	      nc->flag_uf0 != 1 ||
	      nc->flag_xf0 != 1 ||
	      nc->flag_ff != 1)
	    {
	      fprintf (stderr, "F-fix params are not set properly.\n");
	      check ++;
	    }
	}
      else if (sys->version == 1) // FT version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_t0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1 ||
	      nc->flag_o != 1 ||
	      nc->flag_uf0 != 1 ||
	      nc->flag_of0 != 1 ||
	      nc->flag_xf0 != 1 ||
	      nc->flag_ff != 1 ||
	      nc->flag_tf != 1)
	    {
	      fprintf (stderr, "FT-fix params are not set properly.\n");
	      check ++;
	    }
	}
      else // FTS version
	{
	  if (nc->flag_f0 != 1 ||
	      nc->flag_t0 != 1 ||
	      nc->flag_e0 != 1 ||
	      nc->flag_x != 1 ||
	      nc->flag_u != 1 ||
	      nc->flag_o != 1 ||
	      nc->flag_s != 1 ||
	      nc->flag_uf0 != 1 ||
	      nc->flag_of0 != 1 ||
	      nc->flag_ef0 != 1 ||
	      nc->flag_xf0 != 1 ||
	      nc->flag_ff != 1 ||
	      nc->flag_tf != 1 ||
	      nc->flag_sf != 1)
	    {
	      fprintf (stderr, "FTS-fix params are not set properly.\n");
	      check ++;
	    }
	}
    }
  double check_array [5];
  stokes_nc_get_array1d (nc, "Ui0", check_array);
  if (comp_array (Ui, check_array, 3, tiny) != 0)
    {
      fprintf (stderr, "Ui mismatch\n");
      check ++;
    }
  stokes_nc_get_array1d (nc, "Oi0", check_array);
  if (comp_array (Oi, check_array, 3, tiny) != 0)
    {
      fprintf (stderr, "Oi mismatch\n");
      check ++;
    }
  stokes_nc_get_array1d (nc, "Ei0", check_array);
  if (comp_array (Ei, check_array, 5, tiny) != 0)
    {
      fprintf (stderr, "Ei mismatch\n");
      check ++;
    }

  // non-periodic system
  if (sys->periodic == 1)
    {
      stokes_nc_get_array1d (nc, "l", check_array);
      if (comp_array (lat, check_array, 3, tiny) != 0)
	{
	  fprintf (stderr, "l mismatch\n");
	  check ++;
	}
    }

  // auxiliary imposed-flow parameters for simple shear
  if (sys->shear_mode != 0)
    {
      // (int)shear_mode
      int check_int;
      int status;
      status = nc_get_var1_int (nc->id, nc->shear_mode_id, NULL,
				&check_int);
      if (status != NC_NOERR)
	{
	  fprintf (stderr,
		   "at nc_get_var1_int() for shear_mode"
		   " in stokes_nc_check_params\n");
	}
      if (sys->shear_mode != check_int)
	{
	  fprintf (stderr, "shear_mode mismatch\n");
	  check ++;
	}

      // shear_rate
      status = nc_get_var1_double (nc->id, nc->shear_rate_id, NULL,
				   check_array);
      if (status != NC_NOERR)
	{
	  fprintf (stderr,
		   "at nc_get_var1_double() for shear_rate"
		   " in stokes_nc_check_params\n");
	}
      if (comp_array (&(sys->shear_rate), check_array, 1, tiny) != 0)
	{
	  fprintf (stderr, "shear_rate mismatch\n");
	  check ++;
	}
    }

  // first for mobile particle informations (for both mob and mix)
  double *check_parray = (double *)malloc (sizeof (double) * sys->np * 3);
  CHECK_MALLOC (check_parray, "stokes_nc_check_params");
  if (sys->a != NULL)
    {
      stokes_nc_get_array1d (nc, "a", check_parray);
      if (comp_array (sys->a, check_parray, nm, tiny) != 0)
	{
	  fprintf (stderr, "a mismatch\n");
	  check ++;
	}
    }

  if (sys->version == 0)
    {
      // F version
      stokes_nc_get_data0 (nc, "F0", check_parray);
      if (comp_array (F, check_parray, nm3, tiny) != 0)
	{
	  fprintf (stderr, "F0 mismatch\n");
	  check ++;
	}
    }
  else if (sys->version == 1)
    {
      // FT version
      stokes_nc_get_data0 (nc, "F0", check_parray);
      if (comp_array (F, check_parray, nm3, tiny) != 0)
	{
	  fprintf (stderr, "F0 mismatch\n");
	  check ++;
	}
      stokes_nc_get_data0 (nc, "T0", check_parray);
      if (comp_array (T, check_parray, nm3, tiny) != 0)
	{
	  fprintf (stderr, "T0 mismatch\n");
	  check ++;
	}
    }
  else
    {
      // FTS version
      stokes_nc_get_data0 (nc, "F0", check_parray);
      if (comp_array (F, check_parray, nm3, tiny) != 0)
	{
	  fprintf (stderr, "F0 mismatch\n");
	  check ++;
	}
      stokes_nc_get_data0 (nc, "T0", check_parray);
      if (comp_array (T, check_parray, nm3, tiny) != 0)
	{
	  fprintf (stderr, "T0 mismatch\n");
	  check ++;
	}
      stokes_nc_get_data0 (nc, "E0", check_parray);
      if (comp_array (E, check_parray, nm5, tiny) != 0)
	{
	  fprintf (stderr, "E0 mismatch\n");
	  check ++;
	}
    }
  if (nf > 0)
    {
      // mix problem (with fixed particles)
      if (sys->a != NULL)
	{
	  stokes_nc_get_array1d (nc, "af", check_parray);
	  if (comp_array (sys->a + sys->nm,
			  check_parray, nf, tiny) != 0)
	    {
	      fprintf (stderr, "af mismatch\n");
	      check ++;
	    }
	}

      if (sys->version == 0)
	{
	  // F version
	  stokes_nc_get_data0 (nc, "Uf0", check_parray);
	  if (comp_array (uf, check_parray, nf3, tiny) != 0)
	    {
	      fprintf (stderr, "Uf0 mismatch\n");
	      check ++;
	    }
	}
      else if (sys->version == 1)
	{
	  // FT version
	  stokes_nc_get_data0 (nc, "Uf0", check_parray);
	  if (comp_array (uf, check_parray, nf3, tiny) != 0)
	    {
	      fprintf (stderr, "Uf0 mismatch\n");
	      check ++;
	    }
	  stokes_nc_get_data0 (nc, "Of0", check_parray);
	  if (comp_array (of, check_parray, nf3, tiny) != 0)
	    {
	      fprintf (stderr, "Of0 mismatch\n");
	      check ++;
	    }
	}
      else
	{
	  // FTS version
	  stokes_nc_get_data0 (nc, "Uf0", check_parray);
	  if (comp_array (uf, check_parray, nf3, tiny) != 0)
	    {
	      fprintf (stderr, "Uf0 mismatch\n");
	      check ++;
	    }
	  stokes_nc_get_data0 (nc, "Of0", check_parray);
	  if (comp_array (of, check_parray, nf3, tiny) != 0)
	    {
	      fprintf (stderr, "Of0 mismatch\n");
	      check ++;
	    }
	  stokes_nc_get_data0 (nc, "Ef0", check_parray);
	  if (comp_array (ef, check_parray, nf5, tiny) != 0)
	    {
	      fprintf (stderr, "Ef0 mismatch\n");
	      check ++;
	    }
	}
      stokes_nc_get_data0 (nc, "xf0", check_parray);
      if (comp_array (xf, check_parray, nf3, tiny) != 0)
	{
	  fprintf (stderr, "xf0 mismatch\n");
	  check ++;
	}
    }
  free (check_parray);

  return (check);
}
