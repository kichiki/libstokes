/* NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc-read.c,v 5.2 2006/10/22 22:34:05 kichiki Exp $
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

struct stokes_nc *
stokes_nc_open (const char * filename)
{
  struct stokes_nc * nc;
  int status;


  nc = (struct stokes_nc *) malloc (sizeof (struct stokes_nc));
  if (nc == NULL)
    {
      fprintf (stderr, "allocation error in stokes_nc_open()\n");
      exit (1);
    }

  status = nc_open (filename, NC_NOWRITE, &(nc->id));
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
  fprintf (stdout, "ndims      = %d\n", ndims);
  fprintf (stdout, "nvars      = %d\n", nvars);
  fprintf (stdout, "ngatts     = %d\n", ngatts);
  fprintf (stdout, "unlimdimid = %d\n", unlimdimid);

  int i;
  char name[NC_MAX_NAME+1];
  int len;
    nc->npf = 0; // to detect the case of nf == 0
  for (i = 0; i < ndims; i ++)
    {
      status = nc_inq_dim (nc->id, i, name, &len);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status, "at nc_inq_dim() in stokes_nc_open", NULL);
	}
      if (strcmp ("p", name) == 0)
	{
	  nc->p_dim = i;
	  nc->np = len;
	}
      else if (strcmp ("pf", name) == 0)
	{
	  nc->pf_dim = i;
	  nc->npf = len;
	}
      else if (strcmp ("vec", name) == 0)
	{
	  nc->vec_dim = i;
	  nc->nvec = len;
	}
      else if (strcmp ("stt", name) == 0)
	{
	  nc->stt_dim = i;
	  nc->nstt = len;
	}
      else if (strcmp ("time", name) == 0)
	{
	  nc->time_dim = i;
	  nc->ntime = len;
	}
    }

  fprintf (stdout, "np   = %d\n", nc->np);
  fprintf (stdout, "npf  = %d\n", nc->npf);
  fprintf (stdout, "vec  = %d\n", nc->nvec);
  fprintf (stdout, "stt  = %d\n", nc->nstt);
  fprintf (stdout, "time = %d\n", nc->ntime);

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

  nc->flag_xf = 0;
  nc->flag_uf = 0;
  nc->flag_of = 0;
  nc->flag_ef = 0;
  nc->flag_ff = 0;
  nc->flag_tf = 0;
  nc->flag_sf = 0;

  nc_type xtype;
  int nd;
  int d_ids[NC_MAX_VAR_DIMS];
  int natts;
  for (i = 0; i < nvars; i ++)
    {
      status = nc_inq_var (nc->id, i, name, &xtype, &nd, d_ids, &natts);
      if (status != NC_NOERR)
	{
	  stokes_nc_error
	    (status, "at nc_inq_var() in stokes_nc_open", NULL);
	}
      /* data of imposed flow */
      if (strcmp ("Ui0", name) == 0)
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
    }

  return (nc);
}


void
stokes_nc_get_data0 (struct stokes_nc * nc,
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
}

void
stokes_nc_get_data (struct stokes_nc * nc,
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
}
