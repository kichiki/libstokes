/* Brownian dynamics code
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: brownian.c,v 1.4 2007/11/04 03:21:32 kichiki Exp $
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
#include <math.h>

#include <stokes.h>
#include <ewald-3f.h>   // solve_res_3f()
#include <ewald-3ft.h>  // solve_res_3ft()
#include <ewald-3fts.h> // solve_res_3fts()

#include <ewald-3f-matrix.h>   // solve_res_3f_matrix()
#include <ewald-3ft-matrix.h>  // solve_res_3ft_matrix()
#include <ewald-3fts-matrix.h> // solve_res_3fts_matrix()

#include <ewald.h> // atimes_3all()
#include <ode.h> // solve_mix_3all()
#include <ode-quaternion.h> // quaternion_dQdt()
#include "memory-check.h"

#include <KIrand.h> // struct KIrand

#include <chebyshev.h>
#include <dnaupd_c.h>
#include <dsaupd_c.h>
#include <dgetri_c.h> // lapack_inv_()
#include <lub-matrix.h> // make_matrix_lub_3f()
#include <dpotrf_c.h> // dpotrf_wrap()
#include <matrix.h> // mul_matrices()
// check
#include <dgeev_c.h>
#include <libiter.h>
#include <lub.h>

#include <bonds.h> // bonds_calc_force ()
#include <excluded-volume.h> // EV_calc_force ()

#include "brownian.h"


/* set the parameters to struct BD_params
 * INPUT
 *  ** NOTE ** the following pointers are just pointers.
 *             you have to take care of them! (free, for example.)
 *  (struct stokes *)sys -- initialize before calling!
 *  (struct KIrand *)rng -- initialize before calling!
 *  (double *)pos_fix (the position of fixed particles)
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_mat
 *  (int) flag_lub
 *  (double) stokes -- currently this is just a place holder
 *  (struct bonds *)bonds
 *  (double) gamma
 *  (struct EV *)ev
 *  (int) flag_Q
 *  (double) peclet
 *  (double) eps
 *  (int) n_minv
 *  (int) n_lub
 *  (int) scheme
 *  (double) BB_n
 * OUTPUT :
 *  (struct ode_params) params
 */
struct BD_params *
BD_params_init (struct stokes *sys,
		struct KIrand *rng,
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
		double gamma,
		struct EV *ev,
		int flag_Q,
		double peclet,
		double eps,
		int    n_minv,
		int    n_lub,
		int    scheme,
		double BB_n)
{
  struct BD_params *BD
    = (struct BD_params *)malloc (sizeof (struct BD_params));
  CHECK_MALLOC (BD, "BD_params_init");

  BD->sys = sys;
  BD->rng = rng;
  BD->pos_fixed = pos_fixed;
  BD->F = F;
  BD->T = T;
  BD->E = E;
  BD->uf = uf;
  BD->of = of;
  BD->ef = ef;

  BD->flag_lub = flag_lub;
  BD->flag_mat = flag_mat;
  BD->st       = st;
  BD->bonds    = bonds;
  BD->gamma    = gamma;
  BD->ev       = ev;

  BD->flag_Q   = flag_Q;

  BD->peclet = peclet;
  BD->eps = eps;

  BD->n_minv = n_minv;
  BD->eig_minv[0] = 0.0;
  BD->eig_minv[1] = 0.0;
  BD->a_minv = (double *)malloc (sizeof (double) * n_minv);
  CHECK_MALLOC (BD->a_minv, "BD_params_init");
  int i;
  for (i = 0; i < n_minv; i ++)
    {
      BD->a_minv [i] = 0.0;
    }

  BD->n_lub = n_lub;
  BD->eig_lub[0] = 0.0;
  BD->eig_lub[1] = 0.0;
  BD->a_lub = (double *)malloc (sizeof (double) * n_lub);
  CHECK_MALLOC (BD->a_lub, "BD_params_init");
  for (i = 0; i < n_lub; i ++)
    {
      BD->a_lub [i] = 0.0;
    }


  BD->scheme = scheme;
  BD->BB_n = BB_n;

  return (BD);
}

void
BD_params_free (struct BD_params *BD)
{
  if (BD == NULL) return;

  if (BD->a_minv != NULL) free (BD->a_minv);
  if (BD->a_lub  != NULL) free (BD->a_lub);
  free (BD);
}


static double
BD_inv_sqrt (double x)
{
  return (1.0 / sqrt (x));
}

/* similar function (for mono-disperse) is in ewald-2.c */
/* return numbers of overlapping pairs
 */
struct overlap {
  double r2; // square of distance for the overlapping pair
  double a2; // (a_i + a_j)^2 for the overlapping pair
  int i; // particle index
  int j; // particle index
  int k; // lattice index
};

static int
check_overlap (struct stokes *sys, const double *pos,
	       struct overlap *ol)
{
  int flag = 0;
  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      int j;
      for (j = i; j < sys->nm; j ++)
	{
	  double x = pos[i*3+0] - pos[j*3+0];
	  double y = pos[i*3+1] - pos[j*3+1];
	  double z = pos[i*3+2] - pos[j*3+2];
	  double a2;
	  if (sys->a == NULL) a2 = 4.0;
	  else a2 = (sys->a[i] + sys->a[j])*(sys->a[i] + sys->a[j]);
	  double r2;

	  if (i != j)
	    {
	      r2 = x*x + y*y + z*z;
	      if (r2 <= a2)
		{
		  fprintf (stderr, "# overlap %d - %d (0) %e\n",
			   i, j, r2);
		  if (flag == 0)
		    {
		      ol->r2 = r2;
		      ol->a2 = a2;
		      ol->i  = i;
		      ol->j  = j;
		      ol->k  = 0;
		    }
		  else if (r2 < ol->r2)
		    {
		      ol->r2 = r2;
		      ol->a2 = a2;
		      ol->i  = i;
		      ol->j  = j;
		      ol->k  = 0;
		    }
		  flag ++;
		}
	    }
	  if (sys->periodic != 0)
	    {
	      int k;
	      for (k = 1; k < 27; k ++) // excluding the primary cell (k=0)
		{
		  double xx = x - sys->llx[k];
		  double yy = y - sys->lly[k];
		  double zz = z - sys->llz[k];

		  r2 = xx*xx + yy*yy + zz*zz;
		  if (r2 <= a2)
		    {
		      fprintf (stderr, "# overlap %d - %d (%d) %e\n",
			       i, j, k, r2);
		      if (flag == 0)
			{
			  ol->r2 = r2;
			  ol->a2 = a2;
			  ol->i  = i;
			  ol->j  = j;
			  ol->k  = 0;
			}
		      else if (r2 < ol->r2)
			{
			  ol->r2 = r2;
			  ol->a2 = a2;
			  ol->i  = i;
			  ol->j  = j;
			  ol->k  = 0;
			}
		      flag ++;
		    }
		}
	    }
	}
    }
  return (flag);
}



/* calculate y = A.x for Brownian force, where A = M^inf
 * so that give 1/sqrt(x) for chebyshev.
 * Note that for FTS version, A = M_{UF} - M_{US}.(M_{ES})^{-1}.M_{EF}.
 * INPUT
 *  n    : dimension (n = 3 for F, n = 6 for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : U (and O)
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : F (and T)
 */
void
BD_atimes_mob_FU (int n, const double *x, double *y, void *user_data)
{
  struct BD_params *BD = (struct BD_params *)user_data;

  int nm  = BD->sys->nm;
  int nm3 = nm * 3;
  int nm6 = nm * 6;

  int np3 = BD->sys->np * 3;
  int np6 = BD->sys->np * 6;

  int i;

  if (BD->sys->version == 0)
    {
      // F version
      if (n != nm3)
	{
	  fprintf (stderr, "invalid n = %d != 3 * %d\n",
		   n, nm);
	  exit (1);
	}
      // resistance problem in F version
      if (BD->sys->np == nm)
	{
	  atimes_3all (n, x, y, (void *)BD->sys);
	}
      else
	{
	  // with fixed particles
	  double *f = (double *)malloc (sizeof (double) * np3);
	  double *u = (double *)malloc (sizeof (double) * np3);
	  CHECK_MALLOC (f, "BD_atimes_mob_FU");
	  CHECK_MALLOC (u, "BD_atimes_mob_FU");
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = x[i];
	    }
	  for (; i < np3; i ++)
	    {
	      f[i] = 0.0;
	    }

	  atimes_3all (n, f, u, (void *)BD->sys);

	  for (i = 0; i < nm3; i ++)
	    {
	      y[i] = u[i];
	    }
	  free (f);
	  free (u);
	}
    }
  else if (BD->sys->version == 1)
    {
      // FT version
      if (n != 2 * nm3)
	{
	  fprintf (stderr, "invalid n = %d != 6 * %d\n",
		   n, nm);
	  exit (1);
	}
      // resistance problem in FT version
      if (BD->sys->np == nm)
	{
	  atimes_3all (n, x, y, (void *)BD->sys);
	}
      else
	{
	  // with fixed particles
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (ft, "BD_atimes_mob_FU");
	  CHECK_MALLOC (uo, "BD_atimes_mob_FU");
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      ft[i]       = x[i];
	      ft[i + np3] = x[i + nm3];
	    }
	  for (; i < np3; i ++)
	    {
	      ft[i]       = 0.0;
	      ft[i + np3] = 0.0;
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      ft[i] = x[i];
	    }
	  for (; i < np6; i ++)
	    {
	      ft[i] = 0.0;
	    }

	  atimes_3all (n, ft, uo, (void *)BD->sys);

	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      y[i]       = uo[i];
	      y[i + nm3] = uo[i + np3];
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      y[i] = uo[i];
	    }
	  free (ft);
	  free (uo);
	}
    }
  else if (BD->sys->version == 2)
    {
      static struct iter *iter = NULL;
      if (iter == NULL)
	{
	  iter = iter_init ("gmres", 2000, 20, 1.0e-6,
			    BD->sys->np * 11, NULL, 1, // guess
			    0, stderr);
	}

      int ii;
      int nm5 = nm * 5;
      // FTS version
      if (BD->sys->np == nm)
	{
	  double *f = (double *)malloc (sizeof (double) * nm3);
	  double *t = (double *)malloc (sizeof (double) * nm3);
	  double *e = (double *)malloc (sizeof (double) * nm5);
	  double *u = (double *)malloc (sizeof (double) * nm3);
	  double *o = (double *)malloc (sizeof (double) * nm3);
	  double *s = (double *)malloc (sizeof (double) * nm5);
	  CHECK_MALLOC (f, "BD_atimes_mob_FU");
	  CHECK_MALLOC (t, "BD_atimes_mob_FU");
	  CHECK_MALLOC (e, "BD_atimes_mob_FU");
	  CHECK_MALLOC (u, "BD_atimes_mob_FU");
	  CHECK_MALLOC (o, "BD_atimes_mob_FU");
	  CHECK_MALLOC (s, "BD_atimes_mob_FU");
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = x[i];
	      t[i] = x[i + nm3];
	    }
	  for (i = 0; i < nm5; i ++)
	    {
	      e[i] = 0.0;
	    }
	  */
	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  f[i*3+ii] = x[i*6+  ii];
		  t[i*3+ii] = x[i*6+3+ii];
		}
	      for (ii = 0; ii < 5; ii ++)
		{
		  e[i*5+ii] = 0.0;
		}
	    }

	  solve_mob_3fts_0 (BD->sys, iter, f, t, e,
			    u, o, s);
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      y[i      ] = u[i];
	      y[i + nm3] = o[i];
	    }
	  */
	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  y[i*6+  ii] = u[i*3+ii];
		  y[i*6+3+ii] = o[i*3+ii];
		}
	    }
	  free (f);
	  free (t);
	  free (e);
	  free (u);
	  free (o);
	  free (s);
	}
      else
	{
	  int nf = BD->sys->np - nm;
	  int nf3 = nf * 3;
	  int nf5 = nf * 5;
	  double *f = (double *)malloc (sizeof (double) * nm3);
	  double *t = (double *)malloc (sizeof (double) * nm3);
	  double *e = (double *)malloc (sizeof (double) * nm5);
	  double *u = (double *)malloc (sizeof (double) * nm3);
	  double *o = (double *)malloc (sizeof (double) * nm3);
	  double *s = (double *)malloc (sizeof (double) * nm5);
	  double *uf = (double *)malloc (sizeof (double) * nf3);
	  double *of = (double *)malloc (sizeof (double) * nf3);
	  double *ef = (double *)malloc (sizeof (double) * nf5);
	  double *ff = (double *)malloc (sizeof (double) * nf3);
	  double *tf = (double *)malloc (sizeof (double) * nf3);
	  double *sf = (double *)malloc (sizeof (double) * nf5);
	  CHECK_MALLOC (f, "BD_atimes_mob_FU");
	  CHECK_MALLOC (t, "BD_atimes_mob_FU");
	  CHECK_MALLOC (e, "BD_atimes_mob_FU");
	  CHECK_MALLOC (u, "BD_atimes_mob_FU");
	  CHECK_MALLOC (o, "BD_atimes_mob_FU");
	  CHECK_MALLOC (s, "BD_atimes_mob_FU");
	  CHECK_MALLOC (uf, "BD_atimes_mob_FU");
	  CHECK_MALLOC (of, "BD_atimes_mob_FU");
	  CHECK_MALLOC (ef, "BD_atimes_mob_FU");
	  CHECK_MALLOC (ff, "BD_atimes_mob_FU");
	  CHECK_MALLOC (ff, "BD_atimes_mob_FU");
	  CHECK_MALLOC (ff, "BD_atimes_mob_FU");
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      f[i] = x[i];
	      t[i] = x[i + nm3];
	    }
	  */
	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  f[i*3+ii] = x[i*6+  ii];
		  t[i*3+ii] = x[i*6+3+ii];
		}
	    }
	  for (i = 0; i < nm5; i ++)
	    {
	      e[i] = 0.0;
	    }
	  for (i = 0; i < nf3; i ++)
	    {
	      uf[i] = 0.0;
	      of[i] = 0.0;
	    }
	  for (i = 0; i < nf5; i ++)
	    {
	      ef[i] = 0.0;
	    }

	  solve_mix_3fts_0 (BD->sys, iter, f, t, e, uf, of, ef,
			    u, o, s, ff, tf, sf);
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      y[i      ] = u[i];
	      y[i + nm3] = o[i];
	    }
	  */
	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  y[i*6+  ii] = u[i*3+ii];
		  y[i*6+3+ii] = o[i*3+ii];
		}
	    }
	  free (f);
	  free (t);
	  free (e);
	  free (u);
	  free (o);
	  free (s);
	  free (uf);
	  free (of);
	  free (ef);
	  free (ff);
	  free (tf);
	  free (sf);
	}
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", BD->sys->version);
      exit (1);
    }
}

/* calculate y = A.x for Brownian force, where A = L (lubrication)
 * so that give sqrt(x) for chebyshev.
 * note that FU part (for mobile particles) are just extracted
 * (in other words, F for the fixed particles is set by zero,
 *  and E is also set by zero for FTS case).
 * INPUT
 *  n    : dimension (n = 3 for F, n = 6 for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : U (and O)
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : F (and T)
 */
void
BD_atimes_lub_FU (int n, const double *x, double *y, void *user_data)
{
  struct BD_params *BD = (struct BD_params *)user_data;

  int nm3 = BD->sys->nm * 3;
  int nm6 = BD->sys->nm * 6;
  int np3 = BD->sys->np * 3;
  int np6 = BD->sys->np * 6;

  int i;

  if (BD->sys->version == 0)
    {
      // F version
      if (n != nm3)
	{
	  fprintf (stderr, "invalid n = %d != 3 * %d\n",
		   n, BD->sys->nm);
	  exit (1);
	}
      // resistance problem in F version
      if (BD->sys->np == BD->sys->nm)
	{
	  calc_lub_3f (BD->sys, x, y);
	}
      else
	{
	  // with fixed particles
	  double *u = (double *)malloc (sizeof (double) * np3);
	  double *f = (double *)malloc (sizeof (double) * np3);
	  CHECK_MALLOC (u, "BD_atimes_lub_FU");
	  CHECK_MALLOC (f, "BD_atimes_lub_FU");
	  for (i = 0; i < nm3; i ++)
	    {
	      u[i] = x[i];
	    }
	  for (; i < np3; i ++)
	    {
	      u[i] = 0.0;
	    }

	  calc_lub_3f (BD->sys, u, f);

	  for (i = 0; i < nm3; i ++)
	    {
	      y[i] = f[i];
	    }
	  free (u);
	  free (f);
	}
    }
  else if (BD->sys->version == 1)
    {
      // FT version
      if (n != 2 * nm3)
	{
	  fprintf (stderr, "invalid n = %d != 6 * %d\n",
		   n, BD->sys->nm);
	  exit (1);
	}
      // resistance problem in FT version
      if (BD->sys->np == BD->sys->nm)
	{
	  calc_lub_3ft (BD->sys, x, y);
	}
      else
	{
	  // with fixed particles
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (ft, "BD_atimes_lub_FU");
	  CHECK_MALLOC (uo, "BD_atimes_lub_FU");
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      uo[i]       = x[i];
	      uo[i + np3] = x[i + nm3];
	    }
	  for (; i < np3; i ++)
	    {
	      uo[i]       = 0.0;
	      uo[i + np3] = 0.0;
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      uo[i] = x[i];
	    }
	  for (; i < np6; i ++)
	    {
	      uo[i] = 0.0;
	    }

	  calc_lub_3ft (BD->sys, uo, ft);

	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      y[i]       = ft[i];
	      y[i + nm3] = ft[i + np3];
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      y[i] = ft[i];
	    }
	  free (uo);
	  free (ft);
	}
    }
  else if (BD->sys->version == 2)
    {
      // FTS version
      // make struct stokes in FT version
      struct stokes *sys_ft = stokes_copy(BD->sys);
      CHECK_MALLOC (sys_ft, "BD_atimes_mob_FU");
      sys_ft->version = 1; // FT

      if (n != 2 * nm3)
	{
	  fprintf (stderr, "invalid n = %d != 6 * %d\n",
		   n, BD->sys->nm);
	  exit (1);
	}
      // resistance problem in FTS version with E = 0
      // and extract F and T
      if (BD->sys->np == BD->sys->nm)
	{
	  calc_lub_3ft (sys_ft, x, y);
	}
      else
	{
	  // with fixed particles
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (uo, "BD_atimes_lub_FU");
	  CHECK_MALLOC (ft, "BD_atimes_lub_FU");
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      uo[i]       = x[i];
	      uo[i + np3] = x[i + nm3];
	    }
	  for (; i < np3; i ++)
	    {
	      uo[i]       = 0.0;
	      uo[i + np3] = 0.0;
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      uo[i]       = x[i];
	    }
	  for (; i < np6; i ++)
	    {
	      uo[i]       = 0.0;
	    }


	  calc_lub_3ft (sys_ft, uo, ft);
	  /*
	  for (i = 0; i < nm3; i ++)
	    {
	      y[i]       = ft[i];
	      y[i + nm3] = ft[i + np3];
	    }
	  */
	  for (i = 0; i < nm6; i ++)
	    {
	      y[i]       = ft[i];
	    }
	  free (uo);
	  free (ft);
	}

      stokes_free (sys_ft);
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", BD->sys->version);
      exit (1);
    }
}

void
check_symmetric (int n, const double *m, double tiny)
{
  int i, j;
  for (i = 0; i < n; i ++)
    {
      for (j = i+1; i < n; i ++)
	{
	  double d = fabs (m[i*n+j] - m[j*n+i]);
	  if (d > tiny)
	    {
	      fprintf (stderr, "# NO, non-symmetric...\n");
	      return;
	    }
	}
    }
  fprintf (stderr, "# YES, symmetric!\n");
}

/* calc (M^{-1})_{FU} in FTS version
 * where M=(a b) and M^{-1}=(A B), 
 *         (c d)            (C D)
 * A = (a - b.d^{-1}.c)^{-1}.
 * INPUT
 *  np : number of particles
 *  m[np11 * np11] : mobility matrix in FTS
 * OUTPUT
 *  minv_FU[np6 * np6]
 */
void
BD_minv_FU_in_FTS (int np, const double *m, double *minv_FU)
{
  int np5  = np * 5;
  int np6  = np * 6;
  double *a = (double *)malloc (sizeof (double) * np6 * np6);
  double *b = (double *)malloc (sizeof (double) * np6 * np5);
  double *c = (double *)malloc (sizeof (double) * np5 * np6);
  double *d = (double *)malloc (sizeof (double) * np6 * np6);
  double *x = (double *)malloc (sizeof (double) * np5 * np6);
  CHECK_MALLOC (a, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (b, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (c, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (d, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (x, "BD_minv_FU_in_FTS");

  // decompose m into a,b,c,d
  split_matrix_3fts (np, m, a, b, c, d);

  // d = d^{-1}
  lapack_inv_ (np5, d);

  // x[np5,np6] = d[np5,np5] * c[np5,np6]
  mul_matrices (d, np5, np5,
		c, np5, np6,
		x);
  // minv_FU[] = a[] - b[] * x[]
  // D[] = a * A[] + b * B[] . C[]
  add_and_mul (a, np6, np6,
	       b, np6, np5,
	       x, np5, np6,
	       1.0, -1.0,
	       minv_FU);
  /*
  double *y = (double *)malloc (sizeof (double) * np6 * np6);
  CHECK_MALLOC (y, "BD_minv_FU_in_FTS");
  // y[np6,np6] = b[np6,np5] * x[np5,np6]
  mul_matrices (b, np6, np5,
		x, np5, np6,
		y);
  // a = a - b.d^{-1}.c
  int i;
  for (i = 0; i < np6 * np6; i ++)
    {
      minv_FU [i] = a[i] - y[i];
    }
  free (y);
  */
  lapack_inv_ (np6, minv_FU);

  free (a);
  free (b);
  free (c);
  free (d);
  free (x);
}

/* make mobility matrix (M^inf)^{-1} in UF part (for mobile particles)
 * (in other words, F for the fixed particles is set by zero,
 *  and S is also set by zero for FTS case).
 * INPUT
 *  BD   : struct BD_params
 * OUTPUT
 *  minv : (M^inf)^{-1} in UF part
 */
void
BD_matrix_minv_FU (struct BD_params *BD, double *minv)
{
  int nm = BD->sys->nm;
  int nm3 = nm * 3;

  int np3 = BD->sys->np * 3;
  int np6 = BD->sys->np * 6;
  int np11 = BD->sys->np * 11;

  int n;
  int i, j;
  int ii, jj;
  if (BD->sys->version == 0)
    {
      // F version
      n = nm3;

      // resistance problem in F version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_mob_3all (BD->sys, minv);
	  lapack_inv_ (n, minv);
	}
      else
	{
	  // with fixed particles
	  double *m = (double *)malloc (sizeof (double) * np3 * np3);
	  CHECK_MALLOC (m, "BD_matrix_minv_FU");

	  make_matrix_mob_3all (BD->sys, m);
	  lapack_inv_ (np3, m);

	  // extract mobile part
	  for (i = 0; i < nm; i ++)
	    {
	      for (j = 0; j < nm; j ++)
		{
		  for (ii = 0; ii < 3; ii ++)
		    {
		      for (jj = 0; jj < 3; jj ++)
			{
			  minv [(i*3+ii)*n+(j*3+jj)]
			    = m[(i*3+ii)*np3+(j*3+jj)];
			}
		    }
		}
	    }
	  free (m);
	}
    }
  else if (BD->sys->version == 1)
    {
      // FT version
      n = 2 * nm3;

      // resistance problem in FT version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_mob_3all (BD->sys, minv);
	  lapack_inv_ (n, minv);
	}
      else
	{
	  double *m = (double *)malloc (sizeof (double) * np6 * np6);
	  CHECK_MALLOC (m, "BD_matrix_minv_FU");

	  make_matrix_mob_3all (BD->sys, m);
	  lapack_inv_ (np6, m);

	  // extract mobile part
	  for (i = 0; i < nm; i ++)
	    {
	      for (j = 0; j < nm; j ++)
		{
		  for (ii = 0; ii < 6; ii ++)
		    {
		      for (jj = 0; jj < 6; jj ++)
			{
			  minv [(i*6+ii)*n+(j*6+jj)]
			    = m[(i*6+ii)*np6+(j*6+jj)];
			}
		    }
		}
	    }
	  free (m);
	}
    }
  else if (BD->sys->version == 2)
    {
      // FTS version
      n = 2 * nm3;

      // resistance problem in FTS version with E = 0
      double *m  = (double *)malloc (sizeof (double) * np11 * np11);
      double *mi = (double *)malloc (sizeof (double) * np6 * np6);
      CHECK_MALLOC (m,  "BD_matrix_minv_FU");
      CHECK_MALLOC (mi, "BD_matrix_minv_FU");

      make_matrix_mob_3all (BD->sys, m);
      BD_minv_FU_in_FTS (BD->sys->np, m, mi);

      // extract FT mobile part
      for (i = 0; i < nm; i ++)
	{
	  for (j = 0; j < nm; j ++)
	    {
	      for (ii = 0; ii < 6; ii ++)
		{
		  for (jj = 0; jj < 6; jj ++)
		    {
		      minv [(i*6+ii)*n+(j*6+jj)]
			= mi[(i*6+ii)*np6+(j*6+jj)];
		    }
		}
	    }
	}
      free (m);
      free (mi);
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", BD->sys->version);
      exit (1);
    }
}

/* make lubrication matrix L in UF part (for mobile particles)
 * (in other words, F for the fixed particles is set by zero,
 *  and S is also set by zero for FTS case).
 * INPUT
 *  n    : dimension (n = 3 for F, n = 6 for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : U (and O)
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : F (and T)
 */
void
BD_matrix_lub_FU (struct BD_params *BD, double *lub)
{
  int nm = BD->sys->nm;
  int nm3 = nm * 3;

  int np3 = BD->sys->np * 3;
  int np6 = BD->sys->np * 6;
  int np9 = np3 * np3;
  int np36 = np6 * np6;

  int n;
  int i, j;
  int ii, jj;
  if (BD->sys->version == 0)
    {
      // F version
      n = nm3;

      // resistance problem in F version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3f (BD->sys, lub);
	}
      else
	{
	  // with fixed particles
	  double *m = (double *)malloc (sizeof (double) * np9);
	  CHECK_MALLOC (m, "BD_matrix_lub_FU");

	  make_matrix_lub_3f (BD->sys, m);

	  // extract mobile part
	  for (i = 0; i < nm; i ++)
	    {
	      for (j = 0; j < nm; j ++)
		{
		  for (ii = 0; ii < 3; ii ++)
		    {
		      for (jj = 0; jj < 3; jj ++)
			{
			  lub [(i*3+ii)*n+(j*3+jj)]
			    = m[(i*3+ii)*np3+(j*3+jj)];
			}
		    }
		}
	    }
	  free (m);
	}
    }
  else if (BD->sys->version == 1)
    {
      // FT version
      n = 2 * nm3;

      // resistance problem in FT version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3ft (BD->sys, lub);
	}
      else
	{
	  double *m = (double *)malloc (sizeof (double) * np36);
	  CHECK_MALLOC (m, "BD_matrix_lub_FU");

	  make_matrix_lub_3ft (BD->sys, m);

	  // extract mobile part
	  for (i = 0; i < nm; i ++)
	    {
	      for (j = 0; j < nm; j ++)
		{
		  for (ii = 0; ii < 6; ii ++)
		    {
		      for (jj = 0; jj < 6; jj ++)
			{
			  lub [(i*6+ii)*n+(j*6+jj)]
			    = m[(i*6+ii)*np6+(j*6+jj)];
			}
		    }
		}
	    }
	  free (m);
	}
    }
  else if (BD->sys->version == 2)
    {
      // FTS version
      n = 2 * nm3;

      // make struct stokes in FT version
      struct stokes *sys_ft = stokes_copy(BD->sys);
      CHECK_MALLOC (sys_ft, "BD_matrix_lub_FU");
      sys_ft->version = 1; // FT

      // resistance problem in FTS version with E = 0
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3ft (sys_ft, lub);
	}
      else
	{
	  double *m = (double *)malloc (sizeof (double) * np6 * np6);
	  CHECK_MALLOC (m, "BD_matrix_lub_FU");

	  make_matrix_lub_3ft (sys_ft, m);

	  // extract mobile part
	  for (i = 0; i < nm; i ++)
	    {
	      for (j = 0; j < nm; j ++)
		{
		  for (ii = 0; ii < 6; ii ++)
		    {
		      for (jj = 0; jj < 6; jj ++)
			{
			  lub [(i*6+ii)*n+(j*6+jj)]
			    = m[(i*6+ii)*np6+(j*6+jj)];
			}
		    }
		}
	    }
	  free (m);
	}
      stokes_free (sys_ft);
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", BD->sys->version);
      exit (1);
    }
}

static void
atimes_lub_fts (int n, const double *x, double *b, void *user_data)
{
  struct stokes *sys = (struct stokes *)user_data;
  calc_lub_3fts (sys, x, b);
}

static void
set_matrix_by_atimes (int n,
		      void (*atimes)(int, const double *, double *, void *),
		      void *user_data,
		      double *a)
{
  double *x = (double *)malloc (sizeof (double) * n);
  double *b = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "set_matrix_by_atimes");
  CHECK_MALLOC (b, "set_matrix_by_atimes");

  int i, j;
  for (j = 0; j < n; j ++)
    {
      for (i = 0; i < n; i ++)
	{
	  x[i] = 0.0;
	}
      x[j] = 1.0;
      // x_k = delta_{kj}

      // b_i = A_{ik} x_k = A_{ij}
      atimes (n, x, b, user_data);

      // set A_{ij}
      for (i = 0; i < n; i ++)
	{
	  a[i*n+j] = b[i];
	}
    }

  free (x);
  free (b);
}

static void
dgeev_min_max (int n,
	       void (*atimes)(int, const double *, double *, void *),
	       void *user_data)
{
  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *wr = (double *)malloc (sizeof (double) * n);
  double *wi = (double *)malloc (sizeof (double) * n);
  double *v  = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "dgeev_min_max");
  CHECK_MALLOC (wr, "dgeev_min_max");
  CHECK_MALLOC (wi, "dgeev_min_max");
  CHECK_MALLOC (v,  "dgeev_min_max");

  set_matrix_by_atimes (n, atimes, user_data, a);

  dgeev_wrap (n, a, wr, wi, v);

  double lmin = wr[0];
  double lmax = wr[0];
  int i;
  for (i = 1; i < n; i ++)
    {
      if (lmin > wr[i]) lmin = wr[i];
      if (lmax < wr[i]) lmax = wr[i];
    }

  fprintf (stderr, "# dgeev eig = %e, %e\n", lmin, lmax);

  free (a);
  free (wr);
  free (wi);
  free (v);
}


static void
transpose (int n, const double *a, double *at)
{
  int i, j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  at[i*n+j] = a[j*n+i];
	}
    }
}

/* calculate sqrt of the matrix a[n*n] by dgeev()
 */
int
BD_sqrt_by_dgeev (int n, const double *a, double *s)
{
  double tiny = 1.0e-15;

  double *wr = (double *)malloc (sizeof (double) * n);
  double *wi = (double *)malloc (sizeof (double) * n);
  double *v  = (double *)malloc (sizeof (double) * n * n);
  double *vi = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (wr, "BD_sqrt_by_dgeev");
  CHECK_MALLOC (wi, "BD_sqrt_by_dgeev");
  CHECK_MALLOC (v,  "BD_sqrt_by_dgeev");
  CHECK_MALLOC (vi, "BD_sqrt_by_dgeev");

  dgeev_wrap (n, a, wr, wi, v);


  // vi = (v^t)^{-1}
  transpose (n, v, vi);
  lapack_inv_ (n, vi);

  int i, j, k;
  // s = A^(1/2)
  for (i = 0; i < n; i ++)
    {
      if (fabs (wi[i]) > tiny || wr[i] < 0.0)
	{
	  free (wr);
	  free (wi);
	  free (v);
	  free (vi);

	  return (-1);
	} 
      for (j = 0; j < n; j ++)
	{
	  s[i*n+j] = 0.0;
	  for (k = 0; k < n; k ++)
	    {
	      s[i*n+j] += v[k*n+i] * vi[k*n+j] * sqrt (wr[k]);
	    }
	}
    }

  free (wr);
  free (wi);
  free (v);
  free (vi);

  return (0);
}


/*
 * INPUT
 *  BD   : struct BD_params
 *         (sys, rng, eig, n_minv, a_minv, n_lub, a_lub, eps are used)
 * OUTPUT
 *  z[n] : random vector, with which F^B = z * sqrt(2/(peclet * dt))
 *         in FT and FTS, first nm3 are the force, the next nm3 are the torque
 *         (different strage from FTS where f,t,s are ordered particle-wise).
 */
void
calc_brownian_force (struct BD_params *BD,
		     double *z)
{
  int n = BD_get_n (BD->sys);

  // this is ordered particle-wise
  double *zp = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (zp, "calc_brownian_force");

  /**
   * M^{-1} part
   */
  double *y_minv = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y_minv, "calc_brownian_force");
  int i;
  //int j;
  for (i = 0; i < n; i ++)
    {
      y_minv[i] = KIrand_Gaussian (BD->rng);
    }

  int status;
  if (BD->n_minv <= 0)
    {
      // matrix version
      double *minv = (double *)malloc (sizeof (double) * n * n);
      double *l    = (double *)malloc (sizeof (double) * n * n);
      BD_matrix_minv_FU (BD, minv);
      // check
      //check_symmetric (n, minv, 1.0e-12);
      status = dpotrf_wrap (n, minv, l);
      if (status != 0)
	{
	  fprintf (stderr, "minv : failed in dpotrf() : %d\n", status);
	  status = dpotf2_wrap (n, minv, l);
	  fprintf (stderr, "minv : how about in dpotf2() : %d\n", status);
	  if (status != 0)
	    {
	      status = BD_sqrt_by_dgeev (n, minv, l);
	      fprintf (stderr,
		       "minv : how about in sqrt_by_dgeev() : %d\n",
		       status);
	    }
	}

      /*
      for (i = 0; i < n; i ++)
	{
	  zp[i] = 0.0;
	  for (j = 0; j < n; j ++)
	    {
	      zp[i] += l[i*n+j] * y_minv[j];
	    }
	}
      */
      dot_prod_matrix (l, n, n, y_minv, zp);

      free (minv);
      free (l);
    }
  else
    {
      // atimes version
      // check a_minv[] is defined.
      if (BD->eig_minv[0] == BD->eig_minv[1])
	{
	  // a_minv[] is not defined.
	  dnaupd_wrap_min_max (n, BD->eig_minv,
			       BD_atimes_mob_FU, (void *)BD, BD->eps);
	  fprintf (stderr, "# eig_minv = %e, %e\n",
		   BD->eig_minv[0], BD->eig_minv[1]);
	  chebyshev_coef (BD->n_minv, BD_inv_sqrt,
			  BD->eig_minv[0], BD->eig_minv[1], BD->a_minv);
	}

      double err_minv;
      // first try (without updating a_minv[])
      chebyshev_eval_atimes (BD->n_minv, BD->a_minv,
			     n, y_minv, zp,
			     BD->eig_minv[0], BD->eig_minv[1],
			     BD_atimes_mob_FU, (void *)BD);
      err_minv = chebyshev_error_minvsqrt (n, y_minv, zp,
					   BD_atimes_mob_FU, (void *)BD);
      // check the accuracy
      if (err_minv > BD->eps)
	{
	  fprintf (stderr, "# re-calculate a_minv (err_cheb = %e)", err_minv);
	  // re-estimate a_minv[]
	  dnaupd_wrap_min_max (n, BD->eig_minv,
			       BD_atimes_mob_FU, (void *)BD, BD->eps);
	  chebyshev_coef (BD->n_minv, BD_inv_sqrt,
			  BD->eig_minv[0], BD->eig_minv[1], BD->a_minv);

	  // then, re-evaluate the random vector zp[]
	  chebyshev_eval_atimes (BD->n_minv, BD->a_minv,
				 n, y_minv, zp, 
				 BD->eig_minv[0], BD->eig_minv[1],
				 BD_atimes_mob_FU, (void *)BD);
	  err_minv = chebyshev_error_minvsqrt (n, y_minv, zp,
					       BD_atimes_mob_FU, (void *)BD);
	  fprintf (stderr, " => err = %e\n", err_minv);
	}
      free (y_minv);
    }

  /**
   * lubrication part
   */
  if (BD->flag_lub != 0)
    {
      double *y_lub = (double *)malloc (sizeof (double) * n);
      double *z_lub = (double *)malloc (sizeof (double) * n);
      CHECK_MALLOC (y_lub, "calc_brownian_force");
      CHECK_MALLOC (z_lub, "calc_brownian_force");
      int i;
      for (i = 0; i < n; i ++)
	{
	  y_lub[i] = KIrand_Gaussian (BD->rng);
	}

      if (BD->n_lub <= 0)
	{
	  // matrix version
	  double *lub = (double *)malloc (sizeof (double) * n * n);
	  double *l   = (double *)malloc (sizeof (double) * n * n);
	  BD_matrix_lub_FU (BD, lub);
	  // check
	  //check_symmetric (n, lub, 1.0e-12);
	  status = dpotrf_wrap (n, lub, l);
	  if (status != 0)
	    {
	      fprintf (stderr, "lub : failed in dpotrf() : %d\n", status);
	      status = dpotf2_wrap (n, lub, l);
	      fprintf (stderr, "lub : how about in dpotf2() : %d\n", status);
	      if (status != 0)
		{
		  status = BD_sqrt_by_dgeev (n, lub, l);
		  fprintf (stderr,
			   "lub : how about in sqrt_by_dgeev() : %d\n",
			   status);
		}
	    }

	  /*
	  for (i = 0; i < n; i ++)
	    {
	      z_lub[i] = 0.0;
	      for (j = 0; j < n; j ++)
		{
		  z_lub[i] += l[i*n+j] * y_lub[j];
		}
	    }
	  */
	  dot_prod_matrix (l, n, n, y_lub, z_lub);

	  free (lub);
	  free (l);
	}
      else
	{
	  // atimes version
	  // check a_lub[] is defined.
	  if (BD->eig_lub[0] == BD->eig_lub[1])
	    {
	      dnaupd_wrap_min_max (n, BD->eig_lub,
				   BD_atimes_lub_FU, (void *)BD, BD->eps);
	      fprintf (stderr, "# eig_lub = %e, %e\n",
		       BD->eig_lub[0], BD->eig_lub[1]);
	      chebyshev_coef (BD->n_lub, sqrt,
			      BD->eig_lub[0], BD->eig_lub[1], BD->a_lub);
	    }

	  double err_lub;
	  // first try (without updating a_lub[])
	  chebyshev_eval_atimes (BD->n_lub, BD->a_lub,
				 n, y_lub, z_lub, 
				 BD->eig_lub[0], BD->eig_lub[1],
				 BD_atimes_lub_FU, (void *)BD);
	  err_lub = chebyshev_error_Rsqrt (n, y_lub, z_lub,
					   BD_atimes_lub_FU, (void *)BD);
	  // check the accuracy
	  if (err_lub > BD->eps)
	    {
	      fprintf (stderr, "# re-calculate a_lub (err_cheb = %e)", err_lub);

	      // re-estimate a_lub[]
	      dnaupd_wrap_min_max (n, BD->eig_lub,
				   BD_atimes_lub_FU, (void *)BD, BD->eps);
	      chebyshev_coef (BD->n_lub, sqrt,
			      BD->eig_lub[0], BD->eig_lub[1], BD->a_lub);

	      // then, re-evaluate the random vector z[]
	      chebyshev_eval_atimes (BD->n_lub, BD->a_lub,
				     n, y_lub, z_lub, 
				     BD->eig_lub[0], BD->eig_lub[1],
				     BD_atimes_lub_FU, (void *)BD);
	      err_lub = chebyshev_error_Rsqrt (n, y_lub, z_lub,
					       BD_atimes_lub_FU, (void *)BD);
	      fprintf (stderr, " => err = %e\n", err_lub);
	    }
	}
      // add z_lub[] into z[], total random force
      for (i = 0; i < n; i ++)
	{
	  zp[i] += z_lub[i];
	}
      free (y_lub);
      free (z_lub);
    }

  // re-order zp[] to z[]
  if (BD->sys->version == 0)
    {
      for (i = 0; i < n; i ++)
	{
	  z[i] = zp[i];
	}
    }
  else
    {
      int nm3 = BD->sys->nm * 3;
      for (i = 0; i < BD->sys->nm; i ++)
	{
	  int ii;
	  for (ii = 0; ii < 3; ii ++)
	    {
	      z[      i*3 + ii] = zp[i*6 +     ii];
	      z[nm3 + i*3 + ii] = zp[i*6 + 3 + ii];
	    }
	}
    }
  free (zp);
  // now, z[i] * sqrt(2/(peclet * dt)) is F^B_n
}


/*
 * INPUT
 *  sys : struct stokes
 *  u[nm*3] : velocity ( = dx/dt)
 *  o[nm*3] : angluar velocity, which is converted to dq/dt.
 *  x[nm*3] : present position
 *  q[nm*4] : present quaternion
 *  dt      : time step
 * OUTPUT
 *  x[nm*3] : updated position (x += u * dt)
 *  q[nm*4] : present quaternion (q += dq/dt * dt)
 */
void
evolve_Euler_3all (struct stokes *sys,
		   const double *u, const double *o,
		   double dt,
		   double *x, double *q)
{
  int i;
  int nm3 = sys->nm * 3;
  int nm4 = nm3 + sys->nm;
  double *dQdt = NULL;

  if (q != NULL)
    {
      dQdt = (double *)malloc (sizeof (double) * nm4);
      CHECK_MALLOC (dQdt, "evolve_Euler_3all");
    }

  if (sys->version == 0)
    {
      // F version
      for (i = 0; i < nm3; i ++)
	{
	  x[i] += dt * u[i];
	}
    }
  else if (sys->version == 1)
    {
      // FT version
      for (i = 0; i < nm3; i ++)
	{
	  x[i] += dt * u[i];
	}
      // calculate dQdt
      if (q != NULL)
	{
	  for (i = 0; i < sys->nm; i ++)
	    {
	      int i3 = i * 3;
	      int i4 = i3 + i;
	      quaternion_dQdt (q + i4, o + i3, dQdt + i4);
	    }
	  for (i = 0; i < nm4; i ++)
	    {
	      q[i] += dt * dQdt[i];
	    }
	}
    }
  else
    {
      for (i = 0; i < nm3; i ++)
	{
	  x[i] += dt * u[i];
	}
      // calculate dQdt
      if (q != NULL)
	{
	  for (i = 0; i < sys->nm; i ++)
	    {
	      int i3 = i * 3;
	      int i4 = i3 + i;
	      quaternion_dQdt (q + i4, o + i3, dQdt + i4);
	    }
	  for (i = 0; i < nm4; i ++)
	    {
	      q[i] += dt * dQdt[i];
	    }
	}
    }

  if (q != NULL) free (dQdt);
}


/*
 * INPUT
 *  sys : struct stokes
 *        version, np, nm are used.
 */
struct FTS *
FTS_init (struct stokes *sys)
{
  int nm = sys->nm;
  int nf = sys->np - nm;

  struct FTS *FTS = (struct FTS *)malloc (sizeof (struct FTS));
  CHECK_MALLOC (FTS, "FTS_init");
  FTS->f  = NULL;
  FTS->t  = NULL;
  FTS->s  = NULL;
  FTS->ff = NULL;
  FTS->tf = NULL;
  FTS->sf = NULL;
  FTS->u  = NULL;
  FTS->o  = NULL;
  FTS->e  = NULL;
  FTS->uf = NULL;
  FTS->of = NULL;
  FTS->ef = NULL;

  int nm3 = nm * 3;
  int nm4 = nm3 + nm;
  int nm5 = nm4 + nm;
  FTS->f = (double *)malloc (sizeof (double) * nm3);
  FTS->u = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (FTS->f, "FTS_init");
  CHECK_MALLOC (FTS->u, "FTS_init");

  int nf3 = nf * 3;
  int nf4 = nf3 + nf;
  int nf5 = nf4 + nf;
  if (nf > 0)
    {
      FTS->ff = (double *)malloc (sizeof (double) * nf3);
      CHECK_MALLOC (FTS->ff, "FTS_init");
    }
  if (sys->version > 0)
    {
      FTS->t = (double *)malloc (sizeof (double) * nm3);
      FTS->o = (double *)malloc (sizeof (double) * nm3);
      CHECK_MALLOC (FTS->t, "FTS_init");
      CHECK_MALLOC (FTS->o, "FTS_init");
      if (nf > 0)
	{
	  FTS->tf = (double *)malloc (sizeof (double) * nf3);
	  CHECK_MALLOC (FTS->tf, "FTS_init");
	}
    }
  if (sys->version > 1)
    {
      FTS->e = (double *)malloc (sizeof (double) * nm5);
      FTS->s = (double *)malloc (sizeof (double) * nm5);
      CHECK_MALLOC (FTS->e, "FTS_init");
      CHECK_MALLOC (FTS->s, "FTS_init");
      if (nf > 0)
	{
	  FTS->sf = (double *)malloc (sizeof (double) * nf5);
	  CHECK_MALLOC (FTS->sf, "FTS_init");
	}
    }
  return (FTS);
}

void
FTS_free (struct FTS *FTS)
{
  if (FTS == NULL) return;
  if (FTS->f  != NULL) free (FTS->f);
  if (FTS->t  != NULL) free (FTS->t);
  if (FTS->s  != NULL) free (FTS->s);
  if (FTS->ff != NULL) free (FTS->ff);
  if (FTS->tf != NULL) free (FTS->tf);
  if (FTS->sf != NULL) free (FTS->sf);
  if (FTS->u  != NULL) free (FTS->u);
  if (FTS->o  != NULL) free (FTS->o);
  if (FTS->e  != NULL) free (FTS->e);
  if (FTS->uf != NULL) free (FTS->uf);
  if (FTS->of != NULL) free (FTS->of);
  if (FTS->ef != NULL) free (FTS->ef);
  free (FTS);
}


/* dimension for BD scheme:
 * n = nm * 3 for F version -- velocity of mobile particles
 *   = nm * 6 for FT and FTS versions
 *                -- translational and angular velocities of mobile particles
 */
int
BD_get_n (struct stokes *sys)
{
  int n = 0;
  if (sys->version == 0)
    {
      // F version
      n = sys->nm * 3;
    }
  else if (sys->version == 1)
    {
      // FT version
      n = sys->nm * 6;
    }
  else if (sys->version == 2)
    {
      // FTS version
      n = sys->nm * 6;
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", sys->version);
      exit (1);
    }
  return (n);
}



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
	       double *y)
{
  int nm3 = BD->sys->nm * 3;
  // asign local pointers x[] and q[] by y[].
  double *x = y;
  double *q = NULL;
  if (BD->flag_Q != 0)
    {
      q = y + nm3;
    }

  // local time step
  double dt_local = *dt;

  // loop until *t == t_out
  do
    {
      if ((*t) + dt_local > t_out)
	{
	  dt_local = t_out - (*t);
	}

      switch (BD->scheme)
	{
	case 1: // Banchio-Brady (2003)
	  BD_evolve_BB03 (BD, x, q, dt_local);
	  break;

	case 2: // Ball-Melrose (1997)
	  BD_evolve_BM97 (BD, x, q, dt_local);
	  break;

	case 0: // the mid-point algorithm
	default:
	  BD_evolve_mid (BD, x, q, dt_local);
	  break;
	}

      (*t) += dt_local;
    }
  while ((*t) < t_out);
}

/*
 * INPUT
 *  sys : struct stokes
 *        (llx, lly, llz) are used.
 *  x0[] : initial configuration
 */
static double
reset_dt_by_ol (struct stokes *sys,
		double dt0, const double *x0, const struct overlap *ol)
{
  int ix = ol->i * 3;
  int iy = ix + 1;
  int iz = iy + 1;
  int jx = ol->j * 3;
  int jy = jx + 1;
  int jz = jy + 1;
  double dx0 = x0[ix] - x0[jx];
  double dy0 = x0[iy] - x0[jy];
  double dz0 = x0[iz] - x0[jz];
  double r02 = dx0*dx0 + dy0*dy0 + dz0*dz0;
  double dr2 = r02 - ol->r2;
  /* for dt0, r02 => ol->r2 (reduces by dr2).
   * therefore the velocity = dr2 / dt0.
   * if the pair reaches a2 for dt, vel = (r02 - a2) / dt
   * => dt = (r02 - a2) / vel = dt0 * (r02 - a2) / dr2.
   */
  double dt = dt0; // for compiler warning...
  if (r02 < ol->a2) // initial configuration is already overlapping
    {
      fprintf (stderr, "initial configuration is already overlapping\n");
      dt *= 0.5;
    }
  else
    {
      dt = 0.5 * dt0 * (r02 - ol->a2) / dr2;
      // 0.9 is the safety factor
    }
  return (dt);
}

/* evolve position of particles -- the mid-point scheme
 * INPUT
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 */
void
BD_evolve_mid (struct BD_params *BD,
	       double *x, double *q,
	       double dt)
{
  struct overlap ol;
  if (check_overlap (BD->sys, x, &ol) > 0)
    {
      fprintf (stderr, "# mid-point: overlap in the initial step\n");
    }

  /**
   * memory allocation for working area
   */
  int n = BD_get_n (BD->sys);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "BD_evolve_mid");

  // fts vectors for local use
  // uf, of, ef are not used (those in BD are used)
  struct FTS *FTS = FTS_init (BD->sys);
  CHECK_MALLOC (FTS, "BD_evolve_mid");

  int nm3 = BD->sys->nm * 3;
  int nm4 = nm3 + BD->sys->nm;
  int nm5 = nm4 + BD->sys->nm;

  // to keep the initial configuration for re-do.
  double *xmid = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (xmid, "BD_evolve_mid");
  double *qmid = NULL;
  if (BD->sys->version > 0 && q != NULL)
    {
      qmid = (double *)malloc (sizeof (double) * nm4);
      CHECK_MALLOC (qmid, "BD_evolve_mid");
    }

  double fact = sqrt(2.0 / (BD->peclet * dt));

  int i;

  //BD_evolve_mid_REDO:
  // set configuration for calc_brownian_force()
  stokes_set_pos (BD->sys, x);
  calc_brownian_force (BD, z);
  // now, F^B_n = fact * z[i]

 BD_evolve_mid_REDO_scale:
  // (re-)set the configuration (for re-do).
  for (i = 0; i < nm3; i ++)
    {
      xmid[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qmid[i] = q[i];
	}
    }


  /*******************************************
   * 1 : solve for the initial configuration *
   *******************************************/
  // set force (and torque)
  if (BD->sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
      for (i = 0; i < nm3; i ++)
	{
	  FTS->t[i] = BD->T[i] + fact * z[nm3 + i];
	}
      if (BD->sys->version == 2)
	{
	  // FTS version
	  for (i = 0; i < nm5; i ++)
	    {
	      FTS->e[i] = BD->E[i];
	    }
	}
    }
  if (BD->bonds->n > 0)
    {
      // calc bond (spring) force
      bonds_calc_force (BD->bonds, BD->sys,
			FTS->f, 1/* add */);
    }
  if (BD->ev != NULL)
    {
      // calc EV force
      EV_calc_force (BD->ev, BD->sys,
		     FTS->f, 1/* add */);
    }

  // calc dydt
  // BD->sys->pos is already definted above
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f, FTS->t, FTS->e,
		  BD->uf, BD->of, BD->ef,
		  FTS->u, FTS->o, FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  // evolve by Euler for dt/2 from the initial configuration
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt * 0.5,
		     xmid, qmid);
  if (check_overlap (BD->sys, xmid, &ol) > 0)
    {
      //fprintf (stderr, "# mid-point: overlap in 1st step\n");
      // just reject the random force z[]
      //goto BD_evolve_mid_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# mid-point: overlap in 1st step. dt = %e\n", dt);
      goto BD_evolve_mid_REDO_scale;
    }


  /************************************************************
   * 2 : solve for the intermediate (mid-point) configuration *
   ************************************************************/
  // set the intermediate configuration
  stokes_set_pos (BD->sys, xmid);
  // set force (and torque) for xmid[]
  if (BD->sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
      for (i = 0; i < nm3; i ++)
	{
	  FTS->t[i] = BD->T[i] + fact * z[nm3 + i];
	}
      if (BD->sys->version == 2)
	{
	  // FTS version
	  for (i = 0; i < nm5; i ++)
	    {
	      FTS->e[i] = BD->E[i];
	    }
	}
    }
  if (BD->bonds->n > 0)
    {
      // calc force on the mobile particles
      bonds_calc_force (BD->bonds, BD->sys,
			FTS->f, 1/* add */);
    }
  if (BD->ev != NULL)
    {
      // calc EV force
      EV_calc_force (BD->ev, BD->sys,
		     FTS->f, 1/* add */);
    }

  // calc dydt
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f, FTS->t, FTS->e,
		  BD->uf, BD->of, BD->ef,
		  FTS->u, FTS->o, FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /*************************************************************
   * 3 : evolve by Euler for dt with the velocity from x and q *
   *************************************************************/
  // (re-)set the configuration by the initial one
  for (i = 0; i < nm3; i ++)
    {
      xmid[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qmid[i] = q[i];
	}
    }
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt,
		     xmid, qmid);
  if (check_overlap (BD->sys, xmid, &ol) > 0)
    {
      //fprintf (stderr, "# mid-point: overlap in 2nd step\n");
      // just reject the random force z[]
      //goto BD_evolve_mid_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# mid-point: overlap in 2nd step. dt = %e\n", dt);
      goto BD_evolve_mid_REDO_scale;
    }

  // final accepted configuration
  for (i = 0; i < nm3; i ++)
    {
      x[i] = xmid[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  q[i] = qmid[i];
	}
    }

  free (z);
  FTS_free (FTS);
  free (xmid);
  if (qmid != NULL) free (qmid);
}

/* evolve position of particles -- Banchio-Brady scheme
 * reference : Banchio and Brady (2003) Phys. Fluids
 * INPUT
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 */
void
BD_evolve_BB03 (struct BD_params *BD,
		double *x, double *q,
		double dt)
{
  struct overlap ol;

  /**
   * memory allocation for working area
   */
  int n = BD_get_n (BD->sys);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "BD_evolve_BB03");

  // fts vectors for local use
  struct FTS *FTS = FTS_init (BD->sys);

  int nm3 = BD->sys->nm * 3;
  int nm4 = nm3 + BD->sys->nm;
  int nm5 = nm4 + BD->sys->nm;
  int nf = BD->sys->np - BD->sys->nm;
  int nf3 = nf * 3;
  int nf5 = nf * 5;

  // for the intermediate configurations
  double *xBB = (double *)malloc (sizeof (double) * nm3);
  double *uBB = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (xBB, "BD_evolve_BB03");
  CHECK_MALLOC (uBB, "BD_evolve_BB03");
  double *qBB = NULL;
  double *oBB = NULL;
  if (BD->sys->version > 0)
    {
      oBB = (double *)malloc (sizeof (double) * nm4);
      CHECK_MALLOC (oBB, "BD_evolve_BB03");
      if (q != NULL)
	{
	  qBB = (double *)malloc (sizeof (double) * nm4);
	  CHECK_MALLOC (qBB, "BD_evolve_BB03");
	}
    }

  double fact = sqrt(2.0 / (BD->peclet * dt));


  int i;

  //BD_evolve_BB03_REDO:
  // (re-)set brownian force
  stokes_set_pos (BD->sys, x);
  calc_brownian_force (BD, z);
  // now, F^B_n = fact * z[i]

 BD_evolve_BB03_REDO_scale:
  // (re-)set the configuration (for re-do).
  for (i = 0; i < nm3; i ++)
    {
      xBB[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qBB[i] = q[i];
	}
    }


  /*****************************
   * 1 : solve stochastic part *
   *****************************/
  // set force (and torque)
  if (BD->sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = fact * z[i];
	}
      for (i = 0; i < nf3; i ++)
	{
	  FTS->uf[i] = 0.0;
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = fact * z[i];
	  FTS->t[i] = fact * z[nm3 + i];
	}
      for (i = 0; i < nf3; i ++)
	{
	  FTS->uf[i] = 0.0;
	  FTS->of[i] = 0.0;
	}
      if (BD->sys->version == 2)
	{
	  // FTS version
	  for (i = 0; i < nm5; i ++)
	    {
	      FTS->e[i] = 0.0;
	    }
	  for (i = 0; i < nf5; i ++)
	    {
	      FTS->ef[i] = 0.0;
	    }
	}
    }

  /* 1-1 : solve the Brownian velocity for the initial configuration
   */
  // BD->sys->pos is already definted above
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  FTS->uf, FTS->of, FTS->ef,
		  uBB,     oBB,     FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /* 1-2 : evolve by Euler for dt/n with the Brownian velocity
   *       to obtain the intermediate configuration
   */
  evolve_Euler_3all (BD->sys, uBB, oBB, dt / BD->BB_n,
		     xBB, qBB);
  if (check_overlap (BD->sys, xBB, &ol) > 0)
    {
      //fprintf (stderr, "# BB03 : overlap in the intermediate step\n");
      // just reject the random force z[]
      //goto BD_evolve_BB03_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# BB03 : overlap in the intermediate step."
	       " dt = %e\n", dt);
      goto BD_evolve_BB03_REDO_scale;
    }

  /* 1-3 : solve the Brownian velocity for the intermediate configuration
   */
  // set the intermediate configuration
  stokes_set_pos (BD->sys, xBB);
  // calc dydt with the same Brownian force
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  FTS->uf, FTS->of, FTS->ef,
		  FTS->u,  FTS->o,  FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /* 1-4 : evaluate Brownian velocity uBB := u_drift + uBB_0,
   *       where u_drift = (n/2)(uBB_int - uBB_0).
   */
  double BBfact = 0.5 * BD->BB_n;
  for (i = 0; i < nm3; i ++)
    {
      uBB[i] = uBB[i] + BBfact * (FTS->u[i] - uBB[i]);
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  oBB[i] = oBB[i] + BBfact * (FTS->o[i] - oBB[i]);
	}
    }

  /********************************
   * 2 : solve deterministic part *
   ********************************/
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  BD->F,   BD->T,   BD->E,
		  BD->uf,  BD->of,  BD->ef,
		  FTS->u,  FTS->o,  FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /* evaluate the total velocity 
   *   FTS->(u,o) := uDet + (n/2)(uBB_int - uBB_0),
   * where (uDet) and (uBB) are in the labo frame,
   * that is, the results of solve_mix_3all()
   * and it is equal to
   *   uinf + vDet + (n/2)(vBB_int - vBB_0),
   * where (vDet = uDet - uinf) and (vBB = uBB - uinf).
   */
  for (i = 0; i < nm3; i ++)
    {
      FTS->u[i] += uBB[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  FTS->o[i] += oBB[i];
	}
    }

  /******************************************************
   * 3 : evolve by Euler for dt with the total velocity *
   *     to obtain the final configuration              *
   ******************************************************/
  /* reset the configuration by x[] and q[]
   * (because this step is from the x[] and q[])
   */
  for (i = 0; i < nm3; i ++)
    {
      xBB[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qBB[i] = q[i];
	}
    }
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt,
		     xBB, qBB);
  if (check_overlap (BD->sys, xBB, &ol) > 0)
    {
      //fprintf (stderr, "# BB03 : overlap in the final step\n");
      // just reject the random force z[]
      //goto BD_evolve_BB03_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# BB03 : overlap in the final step."
	       " dt = %e\n", dt);
      goto BD_evolve_BB03_REDO_scale;
    }

  // final accepted configuration
  for (i = 0; i < nm3; i ++)
    {
      x[i] = xBB[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  q[i] = qBB[i];
	}
    }

  free (z);
  free (xBB);
  free (uBB);
  if (qBB != NULL) free (qBB);
  if (oBB != NULL) free (oBB);

  FTS_free (FTS);
}

/* evolve position of particles -- Ball-Melrose scheme
 * reference : Ball and Melrose (1997)
 * INPUT
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 */
void
BD_evolve_BM97 (struct BD_params *BD,
		double *x, double *q,
		double dt)
{
  struct overlap ol;

  /**
   * memory allocation for working area
   */
  int n = BD_get_n (BD->sys);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "BD_evolve_BM97");

  // fts vectors for local use
  struct FTS *FTS = FTS_init (BD->sys);

  int nm3 = BD->sys->nm * 3;
  int nm4 = nm3 + BD->sys->nm;
  int nm5 = nm4 + BD->sys->nm;

  // for the intermediate configurations
  double *xBM = (double *)malloc (sizeof (double) * nm3);
  double *uBM = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (xBM, "BD_evolve_BM97");
  CHECK_MALLOC (uBM, "BD_evolve_BM97");
  double *qBM = NULL;
  double *oBM = NULL;
  if (BD->sys->version > 0)
    {
      oBM = (double *)malloc (sizeof (double) * nm4);
      CHECK_MALLOC (oBM, "BD_evolve_BM97");
      if (q != NULL)
	{
	  qBM = (double *)malloc (sizeof (double) * nm4);
	  CHECK_MALLOC (qBM, "BD_evolve_BM97");
	}
    }

  double fact = sqrt(2.0 / (BD->peclet * dt));

  int i;


  //BD_evolve_BM97_REDO:
  // (re-)set brownian force
  stokes_set_pos (BD->sys, x);
  calc_brownian_force (BD, z);
  // now, F^B_n = fact * z[i]


 BD_evolve_BM97_REDO_scale:
  // (re-)set the configuration (for re-do).
  for (i = 0; i < nm3; i ++)
    {
      xBM[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qBM[i] = q[i];
	}
    }


  /*******************************************
   * 1 : solve for the initial configuration *
   *******************************************/
  // set force (and torque)
  if (BD->sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i] + fact * z[i];
	}
      for (i = 0; i < nm3; i ++)
	{
	  FTS->t[i] = BD->T[i] + fact * z[nm3 + i];
	}
      if (BD->sys->version == 2)
	{
	  // FTS version
	  for (i = 0; i < nm5; i ++)
	    {
	      FTS->e[i] = BD->E[i];
	    }
	}
    }

  // calc dydt
  // BD->sys->pos is already definted above
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  BD->uf,  BD->of,  BD->ef,
		  uBM,     oBM,     FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  // evolve by Euler for dt from the initial configuration
  evolve_Euler_3all (BD->sys, uBM, oBM, dt,
		     xBM, qBM);
  if (check_overlap (BD->sys, xBM, &ol) > 0)
    {
      //fprintf (stderr, "# BM97: overlap in 1st step\n");
      // just reject the random force z[]
      //goto BD_evolve_BM97_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# BM97: overlap in 1st step. dt = %e\n", dt);
      goto BD_evolve_BM97_REDO_scale;
    }

  /************************************************************
   * 2 : solve for the intermediate (mid-point) configuration *
   ************************************************************/
  // set the intermediate configuration
  stokes_set_pos (BD->sys, xBM);
  /* now BD->F[i] etc. are independent of the configuration,
   * we just skip the update the force f[].
  for (i = 0; i < n; i ++)
    {
      FTS->f[i] = BD->F[i] + fact * z[i];
    }
  */
  solve_mix_3all (BD->sys,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  BD->uf,  BD->of,  BD->ef,
		  FTS->u,  FTS->o,  FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /* evaluate the total velocity 
   *   FTS->(u,o) := (1/2)(uBM_int + uBM_0),
   * where (uBM)'s are in the labo frame.
   */
  for (i = 0; i < nm3; i ++)
    {
      FTS->u[i] = 0.5 * (uBM[i] + FTS->u[i]);
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  FTS->o[i] = 0.5 * (oBM[i] + FTS->o[i]);
	}
    }

  // evolve by Euler for dt from the initial configuration
  /* reset the configuration by x[] and q[]
   * (because this step is for the x[] and q[]
   */
  for (i = 0; i < nm3; i ++)
    {
      xBM[i] = x[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  qBM[i] = q[i];
	}
    }
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt,
		     xBM, qBM);
  if (check_overlap (BD->sys, xBM, &ol) > 0)
    {
      //fprintf (stderr, "# BM97: overlap in 2nd step\n");
      // just reject the random force z[]
      //goto BD_evolve_BM97_REDO;
      dt = reset_dt_by_ol (BD->sys, dt, x, &ol);
      fact = sqrt(2.0 / (BD->peclet * dt));
      fprintf (stderr, "# BM97: overlap in 2nd step. dt = %e\n", dt);
      goto BD_evolve_BM97_REDO_scale;
    }

  // final accepted configuration
  for (i = 0; i < nm3; i ++)
    {
      x[i] = xBM[i];
    }
  if (BD->sys->version > 0 && q != NULL)
    {
      for (i = 0; i < nm4; i ++)
	{
	  q[i] = qBM[i];
	}
    }

  free (z);
  free (xBM);
  free (uBM);
  if (qBM != NULL) free (qBM);
  if (oBM != NULL) free (oBM);

  FTS_free (FTS);
}
