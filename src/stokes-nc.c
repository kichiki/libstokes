/* NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.c,v 5.4 2006/10/22 22:35:42 kichiki Exp $
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

#include <netcdf.h>

#include "stokes-nc.h"


void
stokes_nc_error (int status, const char * message, const char * varname)
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
stokes_nc_print_actives (struct stokes_nc * nc,
			 FILE * out)
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
}


static void
stokes_nc_init_t_p_vec (struct stokes_nc * nc,
			const char * name,
			const char * longname,
			int * var_id)
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
static void
stokes_nc_init_t_p_stt (struct stokes_nc * nc,
			const char * name,
			const char * longname,
			int * var_id)
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

static void
stokes_nc_init_t_pf_vec (struct stokes_nc * nc,
			 const char * name,
			 const char * longname,
			 int * var_id)
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
static void
stokes_nc_init_t_pf_stt (struct stokes_nc * nc,
			 const char * name,
			 const char * longname,
			 int * var_id)
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

static void
stokes_nc_init_p_vec (struct stokes_nc * nc,
		      const char * name,
		      const char * longname,
		      int * var_id)
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
static void
stokes_nc_init_p_stt (struct stokes_nc * nc,
		      const char * name,
		      const char * longname,
		      int * var_id)
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

static void
stokes_nc_init_pf_vec (struct stokes_nc * nc,
		       const char * name,
		       const char * longname,
		       int * var_id)
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
static void
stokes_nc_init_pf_stt (struct stokes_nc * nc,
		       const char * name,
		       const char * longname,
		       int * var_id)
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

static void
stokes_nc_init_vec (struct stokes_nc * nc,
		    const char * name,
		    const char * longname,
		    int * var_id)
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
static void
stokes_nc_init_stt (struct stokes_nc * nc,
		    const char * name,
		    const char * longname,
		    int * var_id)
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

static void
stokes_nc_init_t_vec (struct stokes_nc * nc,
		      const char * name,
		      const char * longname,
		      int * var_id)
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
static void
stokes_nc_init_t_stt (struct stokes_nc * nc,
		      const char * name,
		      const char * longname,
		      int * var_id)
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

/* initialize NetCDF file for libstokes
 * INPUT
 *  nm : number of MOBILE particles
 *       (note that np in NetCDF is only MOBILE particles)
 *  nf : number of fixed particles
 *  flag_* : active/inactive flags for each entry
 * OUTPUT
 *  (returned value) : ncid
 */
static struct stokes_nc *
stokes_nc_init (const char * filename, int nm, int nf,
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
		int flag_xf,
		int flag_uf,
		int flag_of,
		int flag_ef,
		int flag_ff,
		int flag_tf,
		int flag_sf)
{
  struct stokes_nc * nc = NULL;
  int status;
  int dimids[3];

  nc = (struct stokes_nc *) malloc (sizeof (struct stokes_nc));
  if (nc == NULL)
    {
      fprintf (stderr, "stokes_nc_init: allocation error on nc\n");
      exit (0);
    }

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
  nc->flag_xf = flag_xf;
  nc->flag_uf = flag_uf;
  nc->flag_of = flag_of;
  nc->flag_ef = flag_ef;
  nc->flag_ff = flag_ff;
  nc->flag_tf = flag_tf;
  nc->flag_sf = flag_sf;

  /* create netCDF dataset: enter define mode */
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

  /* particle index */
  status = nc_def_dim (nc->id, "p", nc->np, &(nc->p_dim));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for p in stokes_nc_init", NULL);
    }
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


  /* set dim data */
  int i;
  for (i = 0; i < nc->np; i ++)
    {
      status = nc_put_var1_int (nc->id, nc->p_id, &i, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1_int() for p in stokes_nc_init", NULL);
	}
    }
  for (i = 0; i < nc->nvec; i ++)
    {
      status = nc_put_var1_int (nc->id, nc->vec_id, &i, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1() for vec in stokes_nc_init", NULL);
	}
    }
  for (i = 0; i < nc->nstt; i ++)
    {
      status = nc_put_var1_int (nc->id, nc->stt_id, &i, &i);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_put_var1() for stt in stokes_nc_init", NULL);
	}
    }

  if (nf > 0)
    {
      for (i = 0; i < nc->npf; i ++)
	{
	  status = nc_put_var1_int (nc->id, nc->pf_id, &i, &i);
	  if (status != NC_NOERR)
	    {
	      stokes_nc_error
		(status,
		 "at nc_put_var1() for pf in stokes_nc_init", NULL);
	    }
	}
    }

  return (nc);
}

/* initialize NetCDF file for libstokes for mob_F problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, x, u.
 */
struct stokes_nc *
stokes_nc_mob_f_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
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
			 1, // f0
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
			 1, // u
			 0, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 0, // ff
			 0, // tf
			 0);// sf
}
/* initialize NetCDF file for libstokes for mob_FT problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, t0, x, u, o.
 */
struct stokes_nc *
stokes_nc_mob_ft_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
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
			 1, // f0
			 1, // t0
			 0, // s0
			 0, // xf0
			 0, // uf0
			 0, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 0, // ff
			 0, // tf
			 0);// sf
}
/* initialize NetCDF file for libstokes for mob_FTS problem
 * INPUT
 *  np : number of MOBILE particles
 * OUTPUT
 *  (returned value) : ncid
 *  activated entries are, f0, t0, e0, x, u, o, s.
 */
struct stokes_nc *
stokes_nc_mob_fts_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
			 0, // ui0
			 0, // oi0
			 0, // ei0
			 0, // ui
			 0, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 1, // e0
			 1, // f0
			 1, // t0
			 0, // s0
			 0, // xf0
			 0, // uf0
			 0, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 1, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 0, // ff
			 0, // tf
			 0);// sf
}
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
stokes_nc_mob_fix_f_i0_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 1, // ui0
			 0, // oi0
			 0, // ei0
			 0, // ui
			 0, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 0, // e0
			 1, // f0
			 0, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 0, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 0, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 0, // tf
			 0);// sf
}
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
stokes_nc_mob_fix_ft_i0_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 1, // ui0
			 1, // oi0
			 0, // ei0
			 0, // ui
			 0, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 0, // e0
			 1, // f0
			 1, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 1, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 1, // tf
			 0);// sf
}
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
stokes_nc_mob_fix_fts_i0_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 1, // ui0
			 1, // oi0
			 1, // ei0
			 0, // ui
			 0, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 1, // e0
			 1, // f0
			 1, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 1, // of0
			 1, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 1, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 1, // tf
			 1);// sf
}
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
stokes_nc_mob_fix_f_it_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 0, // ui0
			 0, // oi0
			 0, // ei0
			 1, // ui
			 0, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 0, // e0
			 1, // f0
			 0, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 0, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 0, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 0, // tf
			 0);// sf
}
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
stokes_nc_mob_fix_ft_it_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 0, // ui0
			 0, // oi0
			 0, // ei0
			 1, // ui
			 1, // oi
			 0, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 0, // e0
			 1, // f0
			 1, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 1, // of0
			 0, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 0, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 1, // tf
			 0);// sf
}
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
stokes_nc_mob_fix_fts_it_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
			 0, // ui0
			 0, // oi0
			 0, // ei0
			 1, // ui
			 1, // oi
			 1, // ei
			 0, // x0
			 0, // u0
			 0, // o0
			 1, // e0
			 1, // f0
			 1, // t0
			 0, // s0
			 1, // xf0
			 1, // uf0
			 1, // of0
			 1, // ef0
			 0, // ff0
			 0, // tf0
			 0, // sf0
			 1, // x
			 1, // u
			 1, // o
			 0, // e
			 0, // f
			 0, // t
			 1, // s
			 0, // xf
			 0, // uf
			 0, // of
			 0, // ef
			 1, // ff
			 1, // tf
			 1);// sf
}


/* close (and write if necessary) NetCDF file for libstokes
 */
void
stokes_nc_free (struct stokes_nc * nc)
{
  if (nc == NULL) return;

  int status;
  /* close: save new netCDF dataset */
  status = nc_close (nc->id);
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_close() in stokes_nc_free", NULL);
    }

  free (nc);
}

/** set nc data **/
/* set l = (lx, ly, lz)
 */
void
stokes_nc_set_l (struct stokes_nc * nc,
		 const double * l)
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
	 "at nc_put_vara_double() for l in stokes_nc_append", NULL);
    }
}
/* set ui
 */
void
stokes_nc_set_ui (struct stokes_nc * nc,
		  const double * ui)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->ui_id, start, count, ui);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for ui in stokes_nc_append", NULL);
    }
}
/* set oi
 */
void
stokes_nc_set_oi (struct stokes_nc * nc,
		  const double * oi)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->oi_id, start, count, oi);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for oi in stokes_nc_append", NULL);
    }
}
/* set ei
 */
void
stokes_nc_set_ei (struct stokes_nc * nc,
		  const double * ei)
{
  size_t start[1];
  size_t count[1];

  int status;

  start[0] = 0;

  count[0] = nc->nstt;

  status = nc_put_vara_double(nc->id, nc->ei_id, start, count, ei);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for ei in stokes_nc_append", NULL);
    }
}

/* set x0
 */
void
stokes_nc_set_x0 (struct stokes_nc * nc,
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
	 "at nc_put_vara_double() for x0 in stokes_nc_append", NULL);
    }
}
/* set u0
 */
void
stokes_nc_set_u0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for u0"
		       " in stokes_nc_append", NULL);
    }
}
/* set o0 data
 */
void
stokes_nc_set_o0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for o0"
		       " in stokes_nc_append", NULL);
    }
}
/* set e0 data
 */
void
stokes_nc_set_e0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for e0"
		       " in stokes_nc_append", NULL);
    }
}
/* set f0 data
 */
void
stokes_nc_set_f0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for f0"
		       " in stokes_nc_append", NULL);
    }
}
/* set t0 data
 */
void
stokes_nc_set_t0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for t0"
		       " in stokes_nc_append", NULL);
    }
}
/* set s0 data
 */
void
stokes_nc_set_s0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for s0"
		       " in stokes_nc_append", NULL);
    }
}
/* set xf0
 */
void
stokes_nc_set_xf0 (struct stokes_nc * nc,
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
	 "at nc_put_vara_double() for xf0 in stokes_nc_append", NULL);
    }
}
/* set uf0
 */
void
stokes_nc_set_uf0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for uf0"
		       " in stokes_nc_append", NULL);
    }
}
/* set of0 data
 */
void
stokes_nc_set_of0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for of0"
		       " in stokes_nc_append", NULL);
    }
}
/* set ef0 data
 */
void
stokes_nc_set_ef0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for ef0"
		       " in stokes_nc_append", NULL);
    }
}
/* set ff0 data
 */
void
stokes_nc_set_ff0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for ff0"
		       " in stokes_nc_append", NULL);
    }
}
/* set tf0 data
 */
void
stokes_nc_set_tf0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for tf0"
		       " in stokes_nc_append", NULL);
    }
}
/* set sf0 data
 */
void
stokes_nc_set_sf0 (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for sf0"
		       " in stokes_nc_append", NULL);
    }
}

/* set time (step)
 */
void
stokes_nc_set_time (struct stokes_nc * nc,
		    int step, double time)
{
  int status;

  /* write values into netCDF variable */
  status = nc_put_var1_double (nc->id, nc->time_id, &step, &time);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for time in stokes_nc_append", NULL);
    }
}
/* set x at time (step)
 */
void
stokes_nc_set_x (struct stokes_nc * nc,
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
	 "at nc_put_vara_double() for x in stokes_nc_append", NULL);
    }
}
/* set u at time (step)
 */
void
stokes_nc_set_u (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for u"
		       " in stokes_nc_append", NULL);
    }
}
/* set o at time (step)
 */
void
stokes_nc_set_o (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for o"
		       " in stokes_nc_append", NULL);
    }
}
/* set e at time (step)
 */
void
stokes_nc_set_e (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for e"
		       " in stokes_nc_append", NULL);
    }
}
/* set f at time (step)
 */
void
stokes_nc_set_f (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for f"
		       " in stokes_nc_append", NULL);
    }
}
/* set t at time (step)
 */
void
stokes_nc_set_t (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for t"
		       " in stokes_nc_append", NULL);
    }
}
/* set s at time (step)
 */
void
stokes_nc_set_s (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for s"
		       " in stokes_nc_append", NULL);
    }
}

/* set xf at time (step)
 */
void
stokes_nc_set_xf (struct stokes_nc * nc,
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
	 "at nc_put_vara_double() for xf in stokes_nc_append", NULL);
    }
}
/* set uf at time (step)
 */
void
stokes_nc_set_uf (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for uf"
		       " in stokes_nc_append", NULL);
    }
}
/* set of at time (step)
 */
void
stokes_nc_set_of (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for of"
		       " in stokes_nc_append", NULL);
    }
}
/* set ef at time (step)
 */
void
stokes_nc_set_ef (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for ef"
		       " in stokes_nc_append", NULL);
    }
}
/* set ff at time (step)
 */
void
stokes_nc_set_ff (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for ff"
		       " in stokes_nc_append", NULL);
    }
}
/* set tf at time (step)
 */
void
stokes_nc_set_tf (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for tf"
		       " in stokes_nc_append", NULL);
    }
}
/* set sf at time (step)
 */
void
stokes_nc_set_sf (struct stokes_nc * nc,
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
		       "at nc_put_vara_double() for sf"
		       " in stokes_nc_append", NULL);
    }
}
