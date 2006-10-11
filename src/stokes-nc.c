/* NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc.c,v 5.1 2006/10/11 20:27:57 ichiki Exp $
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


static void
stokes_nc_error (int status, const char * message, const char * varname)
{
  if (varname != NULL)
    {
      fprintf (stderr, "NetCDF error %d : %s for %s\n",
	       status, message, varname);
    }
  else
    {
      fprintf (stderr, "NetCDF error %d : %s\n", status, message);
    }
  exit (0);
}

static void
stokes_nc_init_vec (struct stokes_nc * nc,
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
  int dimids[3];

  dimids[0] = nc->time_dim;
  dimids[1] = nc->pf_dim;
  dimids[2] = nc->stt_dim;
  status = nc_def_var
    (nc->id, name, NC_DOUBLE, 3, dimids, var_id);
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
stokes_nc_init_vec0 (struct stokes_nc * nc,
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
	(status, "at nc_def_var() in stokes_nc_init_vec0", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_vec0", name);
    }
}
static void
stokes_nc_init_stt0 (struct stokes_nc * nc,
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
	(status, "at nc_def_var() in stokes_nc_init_stt0", name);
    }
  status = nc_put_att_text
    (nc->id, *var_id, "long_name",
     strlen(longname), longname);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_put_att() in stokes_nc_init_stt0", name);
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
  status = nc_def_var (nc->id, "time", NC_INT, 1, dimids, &(nc->time_id));
  if (status != NC_NOERR)
    {
      stokes_nc_error (status,
		       "at nc_def_dim() for time in stokes_nc_init", NULL);
    }
 
  /* data for time-independent particles */
  /* x0 */
  if (flag_x0 != 0)
    {
      stokes_nc_init_vec0 (nc, "x0",
			   "positions of particle center at t0",
			   &nc->x0_id);
    }
  /* U0 */
  if (flag_u0 != 0)
    {
      stokes_nc_init_vec0 (nc, "U0",
			   "translational velocities at t0",
			   &nc->u0_id);
    }
  /* O0 */
  if (flag_o0 != 0)
    {
      stokes_nc_init_vec0 (nc, "O0",
			   "angular velocities at t0",
			   &nc->o0_id);
    }
  /* E0 */
  if (flag_e0 != 0)
    {
      stokes_nc_init_stt0 (nc, "E0",
			   "rate of strain at t0",
			   &nc->e0_id);
    }
  /* F0 */
  if (flag_f0 != 0)
    {
      stokes_nc_init_vec0 (nc, "F0",
			   "forces at t0",
			   &nc->f0_id);
    }
  /* T0 */
  if (flag_t0 != 0)
    {
      stokes_nc_init_vec0 (nc, "T0",
			   "torques at t0",
			   &nc->t0_id);
    }
  /* S0 */
  if (flag_s0 != 0)
    {
      stokes_nc_init_stt0 (nc, "S0",
			   "stresslets at t0",
			   &nc->s0_id);
    }

  /* data for (mobile) particles */
  /* x */
  if (flag_x != 0)
    {
      stokes_nc_init_vec (nc, "x",
			  "positions of particle center",
			  &nc->x_id);
    }
  /* U */
  if (flag_u != 0)
    {
      stokes_nc_init_vec (nc, "U",
			  "translational velocities",
			  &nc->u_id);
    }
  /* O */
  if (flag_o != 0)
    {
      stokes_nc_init_vec (nc, "O",
			  "angular velocities",
			  &nc->o_id);
    }
  /* E */
  if (flag_e != 0)
    {
      stokes_nc_init_stt (nc, "E",
			  "rate of strain",
			  &nc->e_id);
    }
  /* F */
  if (flag_f != 0)
    {
      stokes_nc_init_vec (nc, "F",
			  "forces",
			  &nc->f_id);
    }
  /* T */
  if (flag_t != 0)
    {
      stokes_nc_init_vec (nc, "T",
			  "torques",
			  &nc->t_id);
    }
  /* S */
  if (flag_s != 0)
    {
      stokes_nc_init_stt (nc, "S",
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
	  stokes_nc_init_vec0
	    (nc, "xf0",
	     "positions of fixed particle center at t0",
	     &nc->xf0_id);
	}
      /* Uf0 */
      if (flag_uf0 != 0)
	{
	  stokes_nc_init_vec0
	    (nc, "Uf0",
	     "translational velocities for fixed particles at t0",
	     &nc->uf0_id);
	}
      /* Of0 */
      if (flag_of0 != 0)
	{
	  stokes_nc_init_vec0
	    (nc, "Of0",
	     "angular velocities for fixed particles at t0",
	     &nc->of0_id);
	}
      /* Ef0 */
      if (flag_ef0 != 0)
	{
	  stokes_nc_init_stt0
	    (nc, "Ef0",
	     "rate of strain for fixed particles at t0",
	     &nc->ef0_id);
	}
      /* Ff0 */
      if (flag_ff0 != 0)
	{
	  stokes_nc_init_vec0
	    (nc, "Ff0",
	     "forces for fixed particles at t0",
	     &nc->ff0_id);
	}
      /* Tf0 */
      if (flag_tf0 != 0)
	{
	  stokes_nc_init_vec0
	    (nc, "Tf0",
	     "torques for fixed particles at t0",
	     &nc->tf0_id);
	}
      /* Sf0 */
      if (flag_sf0 != 0)
	{
	  stokes_nc_init_stt0
	    (nc, "Sf0",
	     "stresslets for fixed particles at t0",
	     &nc->sf0_id);
	}
      /* end of data for time-independent fixed particles */

      /* xf */
      if (flag_xf != 0)
	{
	  stokes_nc_init_vec (nc, "xf",
			      "positions of fixed particle center",
			      &nc->xf_id);
	}
      /* Uf */
      if (flag_uf != 0)
	{
	  stokes_nc_init_vec (nc, "Uf",
			      "translational velocities for fixed particles",
			      &nc->uf_id);
	}
      /* Of */
      if (flag_of != 0)
	{
	  stokes_nc_init_vec (nc, "Of",
			      "angular velocities for fixed particles",
			      &nc->of_id);
	}
      /* Ef */
      if (flag_ef != 0)
	{
	  stokes_nc_init_stt (nc, "Ef",
			      "rate of strain for fixed particles",
			      &nc->ef_id);
	}
      /* Ff */
      if (flag_ff != 0)
	{
	  stokes_nc_init_vec (nc, "Ff",
			      "forces for fixed particles",
			      &nc->ff_id);
	}
      /* Tf */
      if (flag_tf != 0)
	{
	  stokes_nc_init_vec (nc, "Tf",
			      "torques for fixed particles",
			      &nc->tf_id);
	}
      /* S */
      if (flag_sf != 0)
	{
	  stokes_nc_init_stt (nc, "Sf",
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
 */
struct stokes_nc *
stokes_nc_mob_f_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
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
 */
struct stokes_nc *
stokes_nc_mob_ft_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
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
 */
struct stokes_nc *
stokes_nc_mob_fts_init (const char * filename, int np)
{
  return stokes_nc_init (filename, np, 0,
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
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_f_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
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
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_ft_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
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
 * INPUT
 *  nm : number of MOBILE particles
 *  nf : number of fixed particles
 * OUTPUT
 *  (returned value) : ncid
 */
struct stokes_nc *
stokes_nc_mob_fix_fts_init (const char * filename, int nm, int nf)
{
  return stokes_nc_init (filename, nm, nf,
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


void
stokes_nc_append (struct stokes_nc * nc,
		  int step, double time,
		  const double * x,
		  const double * u, const double * o, const double * e,
		  const double * f, const double * t, const double * s)
{
  size_t start[3];
  size_t count[3];

  int status;
  /* write values into netCDF variable */
  status = nc_put_var1_double (nc->id, nc->time_id, &step, &time);

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->np;
  count[2] = nc->nvec;

  status = nc_put_vara_double(nc->id, nc->x_id, start, count, x);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_put_vara_double() for x in stokes_nc_append", NULL);
    }

  if (u != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->u_id, start, count, u);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for u"
			   " in stokes_nc_append", NULL);
	}
    }
  if (f != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->f_id, start, count, f);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for f"
			   " in stokes_nc_append", NULL);
	}
    }
  if (o != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->o_id, start, count, o);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for o"
			   " in stokes_nc_append", NULL);
	}
    }
  if (t != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->t_id, start, count, t);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for t"
			   " in stokes_nc_append", NULL);
	}
    }

  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 0;
  count[1] = nc->np;
  count[2] = nc->nstt;

  if (e != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->e_id, start, count, e);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for e"
			   " in stokes_nc_append", NULL);
	}
    }
  if (s != NULL)
    {
      status = nc_put_vara_double(nc->id, nc->s_id, start, count, s);
      if (status != NC_NOERR)
	{
	  stokes_nc_error (status,
			   "at nc_put_vara_double() for s"
			   " in stokes_nc_append", NULL);
	}
    }
}


