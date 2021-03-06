/* NetCDF interface for libstokes
 * Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc-read.c,v 5.12 2008/06/03 02:33:51 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include <netcdf.h>

#include "stokes-nc.h"
#include "stokes-nc-read.h"


static void
stokes_nc_check_1 (int id,
		   nc_type xtype,
		   int ndims, const int *dimids,
		   int dim0_id,
		   int * var_id,
		   int * flag)
{
  if (xtype != NC_DOUBLE ||
      ndims != 1 ||
      dimids[0] != dim0_id)
    {
      fprintf (stderr, "invalid data for time-dependent variable"
	       " in stokes_nc_check_1()\n");
    }
  else
    {
      *var_id = id;
      *flag   = 1;
    }
}
static void
stokes_nc_check_1i (int id,
		    nc_type xtype,
		    int ndims, const int *dimids,
		    int dim0_id,
		    int * var_id,
		    int * flag)
{
  if (xtype != NC_INT ||
      ndims != 1 ||
      dimids[0] != dim0_id)
    {
      fprintf (stderr, "invalid data for time-dependent variable"
	       " in stokes_nc_check_1i()\n");
    }
  else
    {
      *var_id = id;
      *flag   = 1;
    }
}
static void
stokes_nc_check_2 (int id,
		   nc_type xtype,
		   int ndims, const int *dimids,
		   int dim0_id, int dim1_id,
		   int * var_id,
		   int * flag)
{
  if (xtype != NC_DOUBLE ||
      ndims != 2 ||
      dimids[0] != dim0_id ||
      dimids[1] != dim1_id)
    {
      fprintf (stderr, "invalid data for time-dependent variable"
	       " in stokes_nc_check_2()\n");
    }
  else
    {
      *var_id = id;
      *flag   = 1;
    }
}
static void
stokes_nc_check_2i (int id,
		    nc_type xtype,
		    int ndims, const int *dimids,
		    int dim0_id, int dim1_id,
		    int * var_id,
		    int * flag)
{
  if (xtype != NC_INT ||
      ndims != 2 ||
      dimids[0] != dim0_id ||
      dimids[1] != dim1_id)
    {
      fprintf (stderr, "invalid data for time-dependent variable"
	       " in stokes_nc_check_2i()\n");
    }
  else
    {
      *var_id = id;
      *flag   = 1;
    }
}
static void
stokes_nc_check_3 (int id,
		   nc_type xtype,
		   int ndims, const int *dimids,
		   int dim0_id, int dim1_id, int dim2_id,
		   int * var_id,
		   int * flag)
{
  if (xtype != NC_DOUBLE ||
      ndims != 3 ||
      dimids[0] != dim0_id ||
      dimids[1] != dim1_id ||
      dimids[2] != dim2_id)
    {
      fprintf (stderr, "invalid data for time-dependent variable"
	       " in stokes_nc_check_3()\n");
    }
  else
    {
      *var_id = id;
      *flag   = 1;
    }
}


static struct stokes_nc *
stokes_nc_open_ (const char * filename, int omode)
{
  struct stokes_nc *nc
    = (struct stokes_nc *)malloc (sizeof (struct stokes_nc));
  CHECK_MALLOC (nc, "stokes_nc_open_");

  int status = nc_open (filename, omode, &(nc->id));
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_open() in stokes_nc_open", NULL);
    }

  /* ndimsp
   *   Pointer to location for returned number of dimensions
   *   defined for this netCDF dataset.
   * nvarsp
   *   Pointer to location for returned number of variables
   *   defined for this netCDF dataset.
   * ngattsp
   *   Pointer to location for returned number of global
   *   attributes defined for this netCDF dataset.
   * unlimdimidp
   *   Pointer to location for returned ID of the unlimited
   *   dimension, if there is one for this netCDF dataset.
   *   If no unlimited length dimension has been defined, -1 is returned.
   * formatp
   *   Pointer to location for returned format version,
   *   one of NC_FORMAT_CLASSIC, NC_FORMAT_64BIT,
   *   NC_FORMAT_NETCDF4, NC_FORMAT_NETCDF4_CLASSIC. 
   */
  int ndims;
  int nvars;
  int ngatts;
  int unlimdimid;
  status = nc_inq (nc->id,
		   &ndims,
		   &nvars,
		   &ngatts,
		   &unlimdimid);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status, "at nc_inq() in stokes_nc_open", NULL);
    }
  /*
  fprintf (stdout, "ndims      = %d\n", ndims);
  fprintf (stdout, "nvars      = %d\n", nvars);
  fprintf (stdout, "ngatts     = %d\n", ngatts);
  fprintf (stdout, "unlimdimid = %d\n", unlimdimid);
  */

  int i;
  char name[NC_MAX_NAME+1];
  size_t len[1];
  nc->npf = 0; // to detect the case of nf == 0
  for (i = 0; i < ndims; i ++)
    {
      status = nc_inq_dim (nc->id, i, name, len);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status, "at nc_inq_dim() in stokes_nc_open", NULL);
	}
      if (strcmp ("p", name) == 0)
	{
	  nc->p_dim = i;
	  nc->np = len[0];
	}
      else if (strcmp ("pf", name) == 0)
	{
	  nc->pf_dim = i;
	  nc->npf = len[0];
	}
      else if (strcmp ("vec", name) == 0)
	{
	  nc->vec_dim = i;
	  nc->nvec = len[0];
	}
      else if (strcmp ("stt", name) == 0)
	{
	  nc->stt_dim = i;
	  nc->nstt = len[0];
	}
      else if (strcmp ("quat", name) == 0)
	{
	  nc->quat_dim = i;
	  nc->nquat = len[0];
	}
      else if (strcmp ("rng", name) == 0)
	{
	  nc->rng_dim = i;
	  nc->nrng = len[0];
	}
      else if (strcmp ("time", name) == 0)
	{
	  nc->time_dim = i;
	  nc->ntime = len[0];
	}
    }
  /*
  fprintf (stdout, "np   = %d\n", nc->np);
  fprintf (stdout, "npf  = %d\n", nc->npf);
  fprintf (stdout, "vec  = %d\n", nc->nvec);
  fprintf (stdout, "stt  = %d\n", nc->nstt);
  fprintf (stdout, "time = %d\n", nc->ntime);
  */

  /* check active vars */
  nc->flag_ui0 = 0;
  nc->flag_oi0 = 0;
  nc->flag_ei0 = 0;
  nc->flag_ui = 0;
  nc->flag_oi = 0;
  nc->flag_ei = 0;

  nc->flag_x0 = 0;
  nc->flag_u0 = 0;
  nc->flag_o0 = 0;
  nc->flag_e0 = 0;
  nc->flag_f0 = 0;
  nc->flag_t0 = 0;
  nc->flag_s0 = 0;

  nc->flag_xf0 = 0;
  nc->flag_uf0 = 0;
  nc->flag_of0 = 0;
  nc->flag_ef0 = 0;
  nc->flag_ff0 = 0;
  nc->flag_tf0 = 0;
  nc->flag_sf0 = 0;

  nc->flag_x = 0;
  nc->flag_u = 0;
  nc->flag_o = 0;
  nc->flag_e = 0;
  nc->flag_f = 0;
  nc->flag_t = 0;
  nc->flag_s = 0;
  nc->flag_q = 0;

  nc->flag_xf = 0;
  nc->flag_uf = 0;
  nc->flag_of = 0;
  nc->flag_ef = 0;
  nc->flag_ff = 0;
  nc->flag_tf = 0;
  nc->flag_sf = 0;

  nc->flag_a = 0;
  nc->flag_af = 0;

  nc_type xtype;
  int nd;
  int d_ids[NC_MAX_VAR_DIMS];
  int natts;
  int flag_dummy;
  for (i = 0; i < nvars; i ++)
    {
      status = nc_inq_var (nc->id, i, name, &xtype, &nd, d_ids, &natts);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status, "at nc_inq_var() in stokes_nc_open", NULL);
	}
      //fprintf (stdout, "# name = %s\n", name);
      /* dimension data */
      if (strcmp ("p", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable p"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->p_id = d_ids[0];
	    }
	}
      else if (strcmp ("pf", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable pf"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->pf_id = d_ids[0];
	    }
	}
      else if (strcmp ("vec", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable vec"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->vec_id = d_ids[0];
	    }
	}
      else if (strcmp ("stt", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable stt"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->stt_id = d_ids[0];
	    }
	}
      else if (strcmp ("quat", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable quat"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->quat_id = d_ids[0];
	    }
	}
      else if (strcmp ("rng", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 1)
	    {
	      fprintf (stderr, "invalid data for variable rng"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->rng_id = d_ids[0];
	    }
	}
      /* system data */
      else if (strcmp ("shear_mode", name) == 0)
	{
	  if (xtype != NC_INT ||
	      nd != 0)
	    {
	      fprintf (stderr, "invalid data for variable shear_mode"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->shear_mode_id = d_ids[0];
	    }
	}
      else if (strcmp ("shear_rate", name) == 0)
	{
	  if (xtype != NC_DOUBLE ||
	      nd != 0)
	    {
	      fprintf (stderr, "invalid data for variable shear_rate"
		       " in stokes_nc_open()\n");
	    }
	  else
	    {
	      nc->shear_rate_id = d_ids[0];
	    }
	}
      else if (strcmp ("shear_shift", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->time_dim,
			     &(nc->shear_shift_id),
			     &flag_dummy);
	}
      else if (strcmp ("l", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->vec_dim,
			     &(nc->l_id),
			     &flag_dummy);
	}
      else if (strcmp ("time", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->time_dim,
			     &(nc->time_id),
			     &flag_dummy);
	}
      /* data of imposed flow */
      else if (strcmp ("Ui0", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->vec_dim,
			     &(nc->ui0_id),
			     &(nc->flag_ui0));
	}
      else if (strcmp ("Oi0", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->vec_dim,
			     &(nc->oi0_id),
			     &(nc->flag_oi0));
	}
      else if (strcmp ("Ei0", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->stt_dim,
			     &(nc->ei0_id),
			     &(nc->flag_ei0));
	}
      else if (strcmp ("Ui", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->vec_dim,
			     &(nc->ui_id),
			     &(nc->flag_ui));
	}
      else if (strcmp ("Oi", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->vec_dim,
			     &(nc->oi_id),
			     &(nc->flag_oi));
	}
      else if (strcmp ("Ei", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->stt_dim,
			     &(nc->ei_id),
			     &(nc->flag_ei));
	}
      /* data for time-independent particles */
      else if (strcmp ("x0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->vec_dim,
			     &(nc->x0_id),
			     &(nc->flag_x0));
	}
      else if (strcmp ("U0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->vec_dim,
			     &(nc->u0_id),
			     &(nc->flag_u0));
	}
      else if (strcmp ("O0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->vec_dim,
			     &(nc->o0_id),
			     &(nc->flag_o0));
	}
      else if (strcmp ("E0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->stt_dim,
			     &(nc->e0_id),
			     &(nc->flag_e0));
	}
      else if (strcmp ("F0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->vec_dim,
			     &(nc->f0_id),
			     &(nc->flag_f0));
	}
      else if (strcmp ("T0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->vec_dim,
			     &(nc->t0_id),
			     &(nc->flag_t0));
	}
      else if (strcmp ("S0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->p_dim, nc->stt_dim,
			     &(nc->s0_id),
			     &(nc->flag_s0));
	}
      /* data for (mobile) particles */
      else if (strcmp ("x", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->vec_dim,
			     &(nc->x_id),
			     &(nc->flag_x));
	}
      else if (strcmp ("q", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->quat_dim,
			     &(nc->q_id),
			     &(nc->flag_q));
	}
      else if (strcmp ("a", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->p_dim,
			     &(nc->a_id),
			     &(nc->flag_a));
	}
      else if (strcmp ("U", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->vec_dim,
			     &(nc->u_id),
			     &(nc->flag_u));
	}
      else if (strcmp ("O", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->vec_dim,
			     &(nc->o_id),
			     &(nc->flag_o));
	}
      else if (strcmp ("E", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->stt_dim,
			     &(nc->e_id),
			     &(nc->flag_e));
	}
      else if (strcmp ("F", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->vec_dim,
			     &(nc->f_id),
			     &(nc->flag_f));
	}
      else if (strcmp ("T", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->vec_dim,
			     &(nc->t_id),
			     &(nc->flag_t));
	}
      else if (strcmp ("S", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->p_dim, nc->stt_dim,
			     &(nc->s_id),
			     &(nc->flag_s));
	}
      /* data for time-independent fixed particles */
      else if (strcmp ("xf0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->vec_dim,
			     &(nc->xf0_id),
			     &(nc->flag_xf0));
	}
      else if (strcmp ("af", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->pf_dim,
			     &(nc->af_id),
			     &(nc->flag_af));
	}
      else if (strcmp ("Uf0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->vec_dim,
			     &(nc->uf0_id),
			     &(nc->flag_uf0));
	}
      else if (strcmp ("Of0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->vec_dim,
			     &(nc->of0_id),
			     &(nc->flag_of0));
	}
      else if (strcmp ("Ef0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->stt_dim,
			     &(nc->ef0_id),
			     &(nc->flag_ef0));
	}
      else if (strcmp ("Ff0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->vec_dim,
			     &(nc->ff0_id),
			     &(nc->flag_ff0));
	}
      else if (strcmp ("Tf0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->vec_dim,
			     &(nc->tf0_id),
			     &(nc->flag_tf0));
	}
      else if (strcmp ("Sf0", name) == 0)
	{
	  stokes_nc_check_2 (i, xtype, nd, d_ids,
			     nc->pf_dim, nc->stt_dim,
			     &(nc->sf0_id),
			     &(nc->flag_sf0));
	}
      /* data for fixed particles */
      else if (strcmp ("xf", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->vec_dim,
			     &(nc->xf_id),
			     &(nc->flag_xf));
	}
      else if (strcmp ("Uf", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->vec_dim,
			     &(nc->uf_id),
			     &(nc->flag_uf));
	}
      else if (strcmp ("Of", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->vec_dim,
			     &(nc->of_id),
			     &(nc->flag_of));
	}
      else if (strcmp ("Ef", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->stt_dim,
			     &(nc->ef_id),
			     &(nc->flag_ef));
	}
      else if (strcmp ("Ff", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->vec_dim,
			     &(nc->ff_id),
			     &(nc->flag_ff));
	}
      else if (strcmp ("Tf", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->vec_dim,
			     &(nc->tf_id),
			     &(nc->flag_tf));
	}
      else if (strcmp ("Sf", name) == 0)
	{
	  stokes_nc_check_3 (i, xtype, nd, d_ids,
			     nc->time_dim, nc->pf_dim, nc->stt_dim,
			     &(nc->sf_id),
			     &(nc->flag_sf));
	}
      /* data for random number generator */
      else if (strcmp ("mt", name) == 0)
	{
	  stokes_nc_check_2i (i, xtype, nd, d_ids,
			      nc->time_dim, nc->rng_dim,
			      &(nc->mt_id),
			      &(nc->flag_rng));
	}
      else if (strcmp ("mti", name) == 0)
	{
	  stokes_nc_check_1i (i, xtype, nd, d_ids,
			      nc->time_dim,
			      &(nc->mti_id),
			      &(nc->flag_rng));
	}
      else if (strcmp ("mt_Ghs", name) == 0)
	{
	  stokes_nc_check_1i (i, xtype, nd, d_ids,
			      nc->time_dim,
			      &(nc->mt_Ghs_id),
			      &(nc->flag_rng));
	}
      else if (strcmp ("mt_Gs", name) == 0)
	{
	  stokes_nc_check_1 (i, xtype, nd, d_ids,
			     nc->time_dim,
			     &(nc->mt_Gs_id),
			     &(nc->flag_rng));
	}
      else
	{
	  fprintf (stderr, "variable %s is not a member of stokes-nc\n",
		   name);
	}
    }

  return (nc);
}

/* open stokes_nc file in NC_NOWRITE mode
 * this is for usual analysis
 */
struct stokes_nc *
stokes_nc_open (const char * filename)
{
  return (stokes_nc_open_ (filename, NC_NOWRITE));
}

/* open stokes_nc file in NC_WRITE mode
 * this is for continuation (appending the results)
 */
struct stokes_nc *
stokes_nc_reopen (const char * filename)
{
  return (stokes_nc_open_ (filename, NC_WRITE));
}


/* read 1d array [vec/stt/np/npf]
 * INPUT
 *  name : either one of them, Ui0, Oi0, Ei0, Ui, Oi, Ei, a, af, l
 * OUTPUT
 *  x[]
 */
void
stokes_nc_get_array1d (const struct stokes_nc *nc,
		       const char * name,
		       double * x)
{
  size_t start[1];
  size_t count[1];

  int status;


  start[0] = 0;

  if (strcmp ("Ui0", name) == 0)
    {
      count[0] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->ui0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("Oi0", name) == 0)
    {
      count[0] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->oi0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("Ei0", name) == 0)
    {
      count[0] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->ei0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("Ui", name) == 0)
    {
      count[0] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->ui_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("Oi", name) == 0)
    {
      count[0] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->oi_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("Ei", name) == 0)
    {
      count[0] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->ei0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("a", name) == 0)
    {
      count[0] = nc->np;

      status = nc_get_vara_double (nc->id, nc->a_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("af", name) == 0)
    {
      count[0] = nc->npf;

      status = nc_get_vara_double (nc->id, nc->af_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else if (strcmp ("l", name) == 0)
    {
      count[0] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->l_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() for in stokes_nc_get_array1d", name);
	}
    }
  else
    {
      fprintf (stderr, "invalid name for stokes_nc_get_array1d()\n"
	       "which is for Ui0, Oi0, Ei0, Ui, Oi, Ei, a, af, or l.\n");
    }
}

/* read constant data for particles in 2d array [np/npf][vec/stt]
 */
void
stokes_nc_get_data0 (const struct stokes_nc *nc,
		     const char * name,
		     double * x)
{
  size_t start[2];
  size_t count[2];

  int status;


  start[0] = 0;
  start[1] = 0;

  if (strcmp ("x0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->x0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("U0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->u0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("O0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->o0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("E0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->e0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("F0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->f0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("T0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->t0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("S0", name) == 0)
    {
      count[0] = nc->np;
      count[1] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->s0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("xf0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->xf0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Uf0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->uf0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Of0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->of0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Ef0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->ef0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Ff0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->ff0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Tf0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->tf0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else if (strcmp ("Sf0", name) == 0)
    {
      count[0] = nc->npf;
      count[1] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->sf0_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data0", name);
	}
    }
  else
    {
      fprintf (stderr, "invalid name for stokes_nc_get_data0()\n"
	       "which is for x0, U0, O0, E0, F0, T0, S0,\n"
	       "xf0, Uf0, Of0, Ef0, Ff0, Tf0, or Sf0.\n");
    }
}

/* read time-dep. particle data at step in 3d array [step][np/npf][vec/stt]
 */
void
stokes_nc_get_data (const struct stokes_nc *nc,
		    const char * name,
		    int step,
		    double * x)
{
  size_t start[3];
  size_t count[3];

  int status;


  start[0] = step;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;

  if (strcmp ("x", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->x_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("q", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nquat;

      status = nc_get_vara_double (nc->id, nc->q_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("U", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->u_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("O", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->o_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("E", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->e_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("F", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->f_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("T", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->t_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("S", name) == 0)
    {
      count[1] = nc->np;
      count[2] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->s_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("xf", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->xf_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Uf", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->uf_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Of", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->of_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Ef", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->ef_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Ff", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->ff_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Tf", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nvec;

      status = nc_get_vara_double (nc->id, nc->tf_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else if (strcmp ("Sf", name) == 0)
    {
      count[1] = nc->npf;
      count[2] = nc->nstt;

      status = nc_get_vara_double (nc->id, nc->sf_id,
				   start, count, x);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status,
	     "at nc_get_vara_double() in stokes_nc_get_data", name);
	}
    }
  else
    {
      fprintf (stderr, "invalid name for stokes_nc_get_data()\n"
	       "which is for x, U, O, E, F, T, S,\n"
	       "xf, Uf, Of, Ef, Ff, Tf, or Sf.\n");
    }
}

/* read (the whole) time vector
 * INPUT
 *  time[nc->ntime]
 * OUTPUT
 *  time[nc->ntime]
 */
void
stokes_nc_get_time (const struct stokes_nc *nc,
		    double * time)
{
  size_t start;
  size_t count;
  int status;


  start = 0;
  count = nc->ntime;

  status = nc_get_vara_double (nc->id, nc->time_id,
			       &start, &count,
			       time);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_double() in stokes_nc_get_time", "time");
    }
}

/* read time at a step
 * INPUT
 *  step
 * OUTPUT
 *  returned value : time[step]
 */
double
stokes_nc_get_time_step (const struct stokes_nc *nc,
			 int step)
{
  double time;
  size_t index;
  int status;


  index = step;
  status = nc_get_var1_double (nc->id, nc->time_id,
			       &index,
			       &time);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_double() in stokes_nc_get_time", "time");
    }
  return (time);
}

/* read rng data at time (step)
 * INPUT
 *  step
 */
void
stokes_nc_get_rng (struct stokes_nc *nc,
		   int step,
		   struct KIrand *rng)
{
  size_t start[2];
  size_t count[2];

  start[0] = step;
  start[1] = 0;

  count[0] = 1;
  count[1] = nc->nrng;

  size_t index[1];
  index[0] = step;

  int status;


  // mt[nrng]
  int *mt = (int *)malloc (sizeof(int) * MTRNG_N);
  CHECK_MALLOC (mt, "stokes_nc_set_rng");

  status = nc_get_vara_int(nc->id, nc->mt_id, start, count, mt);
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_int() in stokes_nc_get_rng", "mt");
    }

  // copy (unsigned long) array to (int) array
  bcopy (mt, rng->mt, sizeof(int) * MTRNG_N);

  // mti
  status = nc_get_var1_int (nc->id, nc->mti_id, index, &(rng->mti));
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_int() in stokes_nc_get_rng", "mti");
    }

  // Gaussian_has_saved
  status = nc_get_var1_int (nc->id, nc->mt_Ghs_id, index,
			    &(rng->Gaussian_has_saved));
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_int() in stokes_nc_get_rng", "mt_Ghs");
    }

  // Gaussian_saved
  status = nc_get_var1_double (nc->id, nc->mt_Gs_id, index,
			       &(rng->Gaussian_saved));
  if (status != NC_NOERR)
    {
      stokes_nc_error
	(status,
	 "at nc_get_vara_double() in stokes_nc_get_rng", "mt_Gs");
    }
}
