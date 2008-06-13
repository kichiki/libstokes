/* Brownian dynamics code
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: brownian.c,v 1.26 2008/06/13 05:07:52 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

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

#include <bead-rod.h> // struct BeadRod
#include <bonds.h> // bonds_calc_force ()
#include <excluded-volume.h> // EV_calc_force ()
#include <angles.h> // angles_calc_force ()
#include <ev-dh.h> // EV_DH_calc_force ()
#include <ev-LJ.h> // EV_LJ_calc_force ()
#include <confinement.h> // CF_calc_force ()

#include "brownian.h"


/* set the parameters to struct BD_params
 * INPUT
 *  ** NOTE ** the following pointers are just pointers.
 *             you have to take care of them! (free, for example.)
 *  (struct stokes *)sys -- initialize before calling!
 *  seed : for random number generator
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_noHI
 *  (int) flag_lub
 *  (int) flag_mat
 *        NOTE, flag_lub_B is used for BD_calc_FB() where
 *        lub among ONLY mobile particles are taken.
 *        therefore, check for the existance of mobile pair(s)
 *        not excluded by sys->ex_lub for flag_lub_B.
 *  (double) stokes -- currently this is just a place holder
 *  (struct BeadRod *)br
 *  (struct bonds *)bonds
 *  (double) gamma
 *  (struct EV *)ev
 *  (struct angles *)ang
 *  (struct EV_DH *)ev_dh
 *  (struct EV_LJ *)ev_LJ
 *  (struct confinement *)cf
 *  (int) flag_Q
 *  (double) peclet
 *  (double) eps
 *  (int) n_minv
 *  (int) n_lub
 *  (int) scheme
 *  (double) BB_n
 *  (double) rmin
 *  (double) dt_lim
 * OUTPUT :
 *  (struct ode_params) params
 */
struct BD_params *
BD_params_init (struct stokes *sys,
		unsigned long seed,
		double *F,
		double *T,
		double *E,
		double *uf,
		double *of,
		double *ef,
		int flag_noHI,
		int flag_lub,
		int flag_mat,
		double st,
		struct BeadRod *br,
		struct bonds *bonds,
		double gamma,
		struct EV *ev,
		struct angles *ang,
		struct EV_DH *ev_dh,
		struct EV_LJ *ev_LJ,
		struct confinement *cf,
		int flag_Q,
		double peclet,
		double eps,
		int    n_minv,
		int    n_lub,
		int    scheme,
		double BB_n,
		double rmin,
		double dt_lim)
{
  struct BD_params *BD
    = (struct BD_params *)malloc (sizeof (struct BD_params));
  CHECK_MALLOC (BD, "BD_params_init");

  BD->sys = sys;

  BD->rng = KIrand_init ();
  KIrand_init_genrand (BD->rng, seed);

  BD->F = F;
  BD->T = T;
  BD->E = E;
  BD->uf = uf;
  BD->of = of;
  BD->ef = ef;

  BD->flag_noHI = flag_noHI;
  BD->flag_mat = flag_mat;

  BD->flag_lub = flag_lub;
  /* BD->flag_lub_B is used for BD_calc_FB() where, 
   * lubrication among mobile particles are taken into account.
   * So, we need to check there is "at least" one mobile particle
   * not excluded for lub calculation to set BD->flag_lub_B = 1. */
  BD->flag_lub_B = 0;
  int i;
  if (flag_lub != 0 && bonds->n != 0)
  //if (flag_lub != 0 && bonds != NULL)
    {
      for (i = 0; i < sys->nm; i ++)
	{
	  int j;
	  for (j = i+1; j < sys->nm; j ++)
	    {
	      if (list_ex_check (sys->ex_lub, i, j) == 0)
		{
		  // not excluded!
		  BD->flag_lub_B = 1;
		  break;
		}
	    }
	  if (BD->flag_lub_B == 1)
	    {
	      // we've already found the non-excluded pair in mobile particles
	      break;
	    }
	}
    }

  BD->st       = st;

  BD->br       = br;
  BD->bonds    = bonds;
  BD->gamma    = gamma;
  BD->ev       = ev;
  BD->ang      = ang;
  BD->ev_dh    = ev_dh;
  BD->ev_LJ    = ev_LJ;
  BD->cf       = cf;

  BD->flag_Q   = flag_Q;

  BD->peclet = peclet;
  BD->eps = eps;

  BD->n_minv = n_minv;
  BD->eig_minv[0] = 0.0;
  BD->eig_minv[1] = 0.0;
  if (n_minv == 0)
    {
      BD->a_minv = NULL;
    }
  else
    {
      BD->a_minv = (double *)malloc (sizeof (double) * n_minv);
      CHECK_MALLOC (BD->a_minv, "BD_params_init");
      for (i = 0; i < n_minv; i ++)
	{
	  BD->a_minv [i] = 0.0;
	}
    }

  BD->n_lub = n_lub;
  BD->eig_lub[0] = 0.0;
  BD->eig_lub[1] = 0.0;
  if (n_lub == 0)
    {
      BD->a_lub = NULL;
    }
  else
    {
      BD->a_lub = (double *)malloc (sizeof (double) * n_lub);
      CHECK_MALLOC (BD->a_lub, "BD_params_init");
      for (i = 0; i < n_lub; i ++)
	{
	  BD->a_lub [i] = 0.0;
	}
    }


  BD->scheme = scheme;
  BD->BB_n = BB_n;

  BD->rmin   = rmin;
  BD->dt_lim = dt_lim;

  return (BD);
}

void
BD_params_free (struct BD_params *BD)
{
  if (BD == NULL) return;

  KIrand_free (BD->rng);
  if (BD->a_minv != NULL) free (BD->a_minv);
  if (BD->a_lub  != NULL) free (BD->a_lub);
  free (BD);
}

/* set the reference for cell-shift (shear_mode = 1 and 2)
 * INTPUT
 *  t0 : reference time for s0
 *  s0 : cell shift at time t0 (for shear_mode = 1 or 2)
 * OUTPUT
 *  BD->t0, BD->s0 :
 */
void
BD_set_shear_shift_ref (struct BD_params *BD,
			double t0, double s0)
{
  BD->t0 = t0;
  BD->s0 = s0;
}

static double
BD_inv_sqrt (double x)
{
  return (1.0 / sqrt (x));
}


/* return numbers of overlapping pairs
 * INPUT
 *  rmin : factor to judge "overlap"
 *         where (r2 <= rmin * a2) are overlapping
 *         rmin == 1 => contact point is the condition
 *         rmin == 0 => all pairs are NOT overlapping
 * OUTPUT
 *  returned value : number of "overlapping" pairs (r2 <= rmin * a2)
 *                   (0 == no "overlapping" pairs)
 */
int
check_overlap (struct stokes *sys, const double *pos, double rmin,
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
	      if (r2 < rmin * a2)
		{
		  fprintf (stderr, "# overlap %d - %d (0) %e\n",
			   i, j, r2);
		  if (flag == 0 || r2 < ol->r2)
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
		  double xx = x + (double)sys->ilx[k] * sys->lx;
		  double yy = y + (double)sys->ily[k] * sys->ly;
		  double zz = z + (double)sys->ilz[k] * sys->lz;

		  // shift for shear
		  if (sys->shear_mode == 1)
		    {
		      xx += (double)sys->ily[k] * sys->shear_shift;
		    }
		  else if (sys->shear_mode == 2)
		    {
		      xx += (double)sys->ilz[k] * sys->shear_shift;
		    }

		  r2 = xx*xx + yy*yy + zz*zz;
		  if (r2 < rmin * a2)
		    {
		      fprintf (stderr, "# overlap %d - %d (%d) %e\n",
			       i, j, k, r2);
		      if (flag == 0 || r2 < ol->r2)
			{
			  ol->r2 = r2;
			  ol->a2 = a2;
			  ol->i  = i;
			  ol->j  = j;
			  ol->k  = k;
			}
		      flag ++;
		    }
		}
	    }
	}
    }
  return (flag);
}

static void
BD_scale_G_to_P_uo (struct BD_params *BD,
		    const double *uo_G,
		    double *uo_P)
{
  int nm  = BD->sys->nm;
  if (BD->sys->version == 0)
    {
      // F version
      int i;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];

	  int ix = i * 3;
	  int iy = ix + 1;
	  int iz = ix + 2;
	  uo_P[ix] = uo_G[ix] * a;
	  uo_P[iy] = uo_G[iy] * a;
	  uo_P[iz] = uo_G[iz] * a;
	}
    }
  else
    {
      // FT or FTS version
      int i;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];
	  double a3 = a * a * a;

	  int iux = i * 6;
	  int iuy = iux + 1;
	  int iuz = iux + 2;
	  int iox = iux + 3;
	  int ioy = iux + 4;
	  int ioz = iux + 5;

	  uo_P[iux] = uo_G[iux] * a;
	  uo_P[iuy] = uo_G[iuy] * a;
	  uo_P[iuz] = uo_G[iuz] * a;

	  uo_P[iox] = uo_G[iox] * a3;
	  uo_P[ioy] = uo_G[ioy] * a3;
	  uo_P[ioz] = uo_G[ioz] * a3;
	}
    }
}
static void
BD_scale_P_to_G_uo (struct BD_params *BD,
		    const double *uo_P,
		    double *uo_G)
{
  int nm  = BD->sys->nm;
  if (BD->sys->version == 0)
    {
      // F version
      int i;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];

	  int ix = i * 3;
	  int iy = ix + 1;
	  int iz = ix + 2;
	  uo_G[ix] = uo_P[ix] / a;
	  uo_G[iy] = uo_P[iy] / a;
	  uo_G[iz] = uo_P[iz] / a;
	}
    }
  else
    {
      // FT or FTS version
      int i;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];
	  double a3 = a * a * a;

	  int iux = i * 6;
	  int iuy = iux + 1;
	  int iuz = iux + 2;
	  int iox = iux + 3;
	  int ioy = iux + 4;
	  int ioz = iux + 5;

	  uo_G[iux] = uo_P[iux] / a;
	  uo_G[iuy] = uo_P[iuy] / a;
	  uo_G[iuz] = uo_P[iuz] / a;

	  uo_G[iox] = uo_P[iox] / a3;
	  uo_G[ioy] = uo_P[ioy] / a3;
	  uo_G[ioz] = uo_P[ioz] / a3;
	}
    }
}


/* calculate y = A.x for Brownian force, where A = M^inf
 * so that give 1/sqrt(x) for chebyshev.
 * Note that for FTS version, A = M_{UF} - M_{US}.(M_{ES})^{-1}.M_{EF}.
 * INPUT
 *  n    : dimension (n = 3*nm for F, n = 6*nm for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : F (and T) in Global scaling, rather than Particle-wise scaling
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : U (and O) in Global scaling, rather than Particle-wise scaling
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

  // convert from Global scaling to Particle-wise scaling
  double *y_P = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y_P, "BD_atimes_mob_FU");

  if (BD->sys->version == 0)
    {
      // F version
      if (n != nm3)
	{
	  fprintf (stderr, "invalid n = %d != 3 * %d\n",
		   n, nm);
	  exit (1);
	}
      // mobility problem in F version
      if (BD->sys->np == nm)
	{
	  atimes_3all (n, x, y_P, (void *)BD->sys);
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
	      //f[i] = x_P[i];
	    }
	  for (; i < np3; i ++)
	    {
	      f[i] = 0.0;
	    }

	  atimes_3all (n, f, u, (void *)BD->sys);

	  for (i = 0; i < nm3; i ++)
	    {
	      y_P[i] = u[i];
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
      // mobility problem in FT version
      if (BD->sys->np == nm)
	{
	  atimes_3all (n, x, y_P, (void *)BD->sys);
	}
      else
	{
	  // with fixed particles
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (ft, "BD_atimes_mob_FU");
	  CHECK_MALLOC (uo, "BD_atimes_mob_FU");

	  for (i = 0; i < nm6; i ++)
	    {
	      ft[i] = x[i];
	    }
	  for (; i < np6; i ++)
	    {
	      ft[i] = 0.0;
	    }

	  atimes_3all (n, ft, uo, (void *)BD->sys);

	  for (i = 0; i < nm6; i ++)
	    {
	      y_P[i] = uo[i];
	    }
	  free (ft);
	  free (uo);
	}
    }
  else if (BD->sys->version == 2)
    {
      struct iter *iter = iter_init ("gmres", 2000, 20, 1.0e-6,
				     BD->sys->np * 11, NULL, 1, // guess
				     0, stderr);
      CHECK_MALLOC (iter, "BD_atimes_mob_FU");

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

	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  y_P[i*6+  ii] = u[i*3+ii];
		  y_P[i*6+3+ii] = o[i*3+ii];
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

	  for (i = 0; i < nm; i ++)
	    {
	      for (ii = 0; ii < 3; ii ++)
		{
		  y_P[i*6+  ii] = u[i*3+ii];
		  y_P[i*6+3+ii] = o[i*3+ii];
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
      iter_free (iter);
    }
  else
    {
      fprintf (stderr, "invalid version %d\n", BD->sys->version);
      exit (1);
    }

  // convert from Particle-wise scaling to Global scaling
  BD_scale_P_to_G_uo (BD, y_P, y);
  //free (x_P);
  free (y_P);
}

/* calculate y = A.x for Brownian force, where A = L (lubrication)
 * so that give sqrt(x) for chebyshev.
 * note that FU part (for mobile particles) are just extracted
 * (in other words, F for the fixed particles is set by zero,
 *  and E is also set by zero for FTS case).
 * INPUT
 *  n    : dimension (n = 3*nm for F, n = 6*nm for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : U (and O) in Global scaling, rather than Particle-wise scaling
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : F (and T) in Global scaling, rather than Particle-wise scaling
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

  // convert from Global scaling to Particle-wise scaling
  double *x_P = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x_P, "BD_atimes_mob_FU");
  BD_scale_G_to_P_uo (BD, x, x_P);

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
	  calc_lub_3f (BD->sys, x_P, y);
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
	      u[i] = x_P[i];
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
	  calc_lub_3ft (BD->sys, x_P, y);
	}
      else
	{
	  // with fiexd particles
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (ft, "BD_atimes_lub_FU");
	  CHECK_MALLOC (uo, "BD_atimes_lub_FU");

	  for (i = 0; i < nm6; i ++)
	    {
	      uo[i] = x_P[i];
	    }
	  for (; i < np6; i ++)
	    {
	      uo[i] = 0.0;
	    }

	  calc_lub_3ft (BD->sys, uo, ft);

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
	  calc_lub_3ft (sys_ft, x_P, y);
	}
      else
	{
	  // with fixed particles
	  double *uo = (double *)malloc (sizeof (double) * np6);
	  double *ft = (double *)malloc (sizeof (double) * np6);
	  CHECK_MALLOC (uo, "BD_atimes_lub_FU");
	  CHECK_MALLOC (ft, "BD_atimes_lub_FU");

	  for (i = 0; i < nm6; i ++)
	    {
	      uo[i] = x_P[i];
	    }
	  for (; i < np6; i ++)
	    {
	      uo[i] = 0.0;
	    }


	  calc_lub_3ft (sys_ft, uo, ft);

	  for (i = 0; i < nm6; i ++)
	    {
	      y[i] = ft[i];
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

  free (x_P);
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

static void
BD_scale_P_to_G_res (struct BD_params *BD,
		     const double *r_P,
		     double *r_G)
{
  int nm  = BD->sys->nm;
  if (BD->sys->version == 0)
    {
      // F version
      int n = nm * 3;
      int i, j;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];

	  int ifx = i * 3;
	  int ify = ifx + 1;
	  int ifz = ifx + 2;

	  for (j = 0; j < n; j ++)
	    {
	      // A part
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;
	    }
	}
    }
  else
    {
      int n = nm * 6;
      // FT or FTS version
      int i, j;
      for (i = 0; i < nm; i ++)
	{
	  double a;
	  if (BD->sys->a == NULL) a = 1.0;
	  else                    a = BD->sys->a[i];
	  double a3 = a * a * a;

	  int ifx = i * 6;
	  int ify = ifx + 1;
	  int ifz = ifx + 2;
	  int itx = ifx + 3;
	  int ity = ifx + 4;
	  int itz = ifx + 5;

	  for (j = 0; j < n; j ++)
	    {
	      // A part
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ify * n + j] = r_P[ify * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a;

	      // tilde(B) part
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a3;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a3;
	      r_G[ifx * n + j] = r_P[ifx * n + j] * a3;
	      r_G[ify * n + j] = r_P[ify * n + j] * a3;
	      r_G[ify * n + j] = r_P[ify * n + j] * a3;
	      r_G[ify * n + j] = r_P[ify * n + j] * a3;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a3;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a3;
	      r_G[ifz * n + j] = r_P[ifz * n + j] * a3;

	      // B part
	      r_G[itx * n + j] = r_P[itx * n + j] * a;
	      r_G[itx * n + j] = r_P[itx * n + j] * a;
	      r_G[itx * n + j] = r_P[itx * n + j] * a;
	      r_G[ity * n + j] = r_P[ity * n + j] * a;
	      r_G[ity * n + j] = r_P[ity * n + j] * a;
	      r_G[ity * n + j] = r_P[ity * n + j] * a;
	      r_G[itz * n + j] = r_P[itz * n + j] * a;
	      r_G[itz * n + j] = r_P[itz * n + j] * a;
	      r_G[itz * n + j] = r_P[itz * n + j] * a;

	      // C part
	      r_G[itx * n + j] = r_P[itx * n + j] * a3;
	      r_G[itx * n + j] = r_P[itx * n + j] * a3;
	      r_G[itx * n + j] = r_P[itx * n + j] * a3;
	      r_G[ity * n + j] = r_P[ity * n + j] * a3;
	      r_G[ity * n + j] = r_P[ity * n + j] * a3;
	      r_G[ity * n + j] = r_P[ity * n + j] * a3;
	      r_G[itz * n + j] = r_P[itz * n + j] * a3;
	      r_G[itz * n + j] = r_P[itz * n + j] * a3;
	      r_G[itz * n + j] = r_P[itz * n + j] * a3;
	    }
	}
    }
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
  int np6  = np * 6;
  int np5  = np * 5;
  double *a = (double *)malloc (sizeof (double) * np6 * np6);
  double *b = (double *)malloc (sizeof (double) * np6 * np5);
  double *c = (double *)malloc (sizeof (double) * np5 * np6);
  double *d = (double *)malloc (sizeof (double) * np5 * np5);
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
 *  minv : UF part of (M^inf)^{-1} in Global scaling not Particle-wise.
 */
void
BD_matrix_minv_FU (struct BD_params *BD, double *minv)
{
  int nm = BD->sys->nm;
  int nm3 = nm * 3;

  int np3 = BD->sys->np * 3;
  int np6 = BD->sys->np * 6;
  int np11 = BD->sys->np * 11;

  double *minv_P = NULL;

  int n;
  int i, j;
  int ii, jj;
  if (BD->sys->version == 0)
    {
      // F version
      n = nm3;
      minv_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (minv_P, "BD_matrix_minv_FU");

      // resistance problem in F version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_mob_3all (BD->sys, minv_P);
	  lapack_inv_ (n, minv_P);
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
			  minv_P [(i*3+ii)*n+(j*3+jj)]
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
      minv_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (minv_P, "BD_matrix_minv_FU");

      // resistance problem in FT version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_mob_3all (BD->sys, minv_P);
	  lapack_inv_ (n, minv_P);
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
			  minv_P [(i*6+ii)*n+(j*6+jj)]
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
      minv_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (minv_P, "BD_matrix_minv_FU");

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
		      minv_P [(i*6+ii)*n+(j*6+jj)]
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

  // convert from Particle-wise scaling to Global scaling
  BD_scale_P_to_G_res (BD, minv_P, minv);
  free (minv_P);
}

/* make lubrication matrix L in UF part (for mobile particles)
 * (in other words, F for the fixed particles is set by zero,
 *  and S is also set by zero for FTS case).
 * INPUT
 *  BD   : struct BD_params
 * OUTPUT
 *  lub  : UF part of L in Global scaling not Particle-wise.
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

  double *lub_P = NULL;

  int n;
  int i, j;
  int ii, jj;
  if (BD->sys->version == 0)
    {
      // F version
      n = nm3;
      lub_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (lub_P, "BD_matrix_lub_FU");

      // resistance problem in F version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3f (BD->sys, lub_P);
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
			  lub_P [(i*3+ii)*n+(j*3+jj)]
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
      lub_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (lub_P, "BD_matrix_lub_FU");

      // resistance problem in FT version
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3ft (BD->sys, lub_P);
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
			  lub_P [(i*6+ii)*n+(j*6+jj)]
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
      lub_P = (double *)malloc (sizeof (double) * n * n);
      CHECK_MALLOC (lub_P, "BD_matrix_lub_FU");

      // make struct stokes in FT version
      struct stokes *sys_ft = stokes_copy(BD->sys);
      CHECK_MALLOC (sys_ft, "BD_matrix_lub_FU");
      sys_ft->version = 1; // FT

      // resistance problem in FTS version with E = 0
      if (BD->sys->np == BD->sys->nm)
	{
	  make_matrix_lub_3ft (sys_ft, lub_P);
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
			  lub_P [(i*6+ii)*n+(j*6+jj)]
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

  // convert from Particle-wise scaling to Global scaling
  BD_scale_P_to_G_res (BD, lub_P, lub);
  free (lub_P);
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

	  fprintf (stderr, "BD_sqrt_by_dgeev : %e %e\n",
		   fabs (wi[i]), wr[i]);
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
BD_calc_FB (struct BD_params *BD,
	    double *z)
{
  int n = BD_get_n (BD->sys);
  int i;

  if (BD->flag_noHI == 1)
    {
      for (i = 0; i < BD->sys->nm; i ++)
	{
	  double as;
	  if (BD->sys->a == NULL) as = 1.0;
	  else                    as = sqrt (BD->sys->a[i]);

	  if (BD->sys->version == 0)
	    {
	      int ix = i * 3;
	      int iy = ix + 1;
	      int iz = ix + 2;

	      z[ix] = KIrand_Gaussian (BD->rng) * as;
	      z[iy] = KIrand_Gaussian (BD->rng) * as;
	      z[iz] = KIrand_Gaussian (BD->rng) * as;
	    }
	  else
	    {
	      double as3 = as * as * as;
	      // sqrt(3/4) = 0.8660254037844386
	      as3 /= 0.8660254037844386;
	      /* now as3 = sqrt (a^3 / (3/4)) */

	      int ifx = i * 6;
	      int ify = ifx + 1;
	      int ifz = ifx + 2;
	      int itx = ifx + 3;
	      int ity = ifx + 4;
	      int itz = ifx + 5;

	      z[ifx] = KIrand_Gaussian (BD->rng) * as;
	      z[ify] = KIrand_Gaussian (BD->rng) * as;
	      z[ifz] = KIrand_Gaussian (BD->rng) * as;

	      z[itx] = KIrand_Gaussian (BD->rng) * as3;
	      z[ity] = KIrand_Gaussian (BD->rng) * as3;
	      z[itz] = KIrand_Gaussian (BD->rng) * as3;
	    }
	}
      // done
      return;
    }

  // this is ordered particle-wise
  double *zp = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (zp, "BD_calc_FB");

  /**
   * M^{-1} part
   */
  double *y_minv = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y_minv, "BD_calc_FB");
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
      CHECK_MALLOC (minv, "BD_calc_FB");
      CHECK_MALLOC (l,    "BD_calc_FB");
      BD_matrix_minv_FU (BD, minv);
      // check
      status = dpotrf_wrap (n, minv, l);
      if (status != 0)
	{
	  check_symmetric (n, minv, 1.0e-12);
	  struct overlap ol;
	  if (check_overlap (BD->sys, BD->sys->pos, 1.0, &ol) > 0)
	    {
	      fprintf (stderr, "overlap : [%d %d] %e < %e\n",
		       ol.i, ol.j,
		       ol.r2, ol.a2);
	    }
	  fprintf (stderr, "minv : failed in dpotrf() : %d\n", status);
	  status = dpotf2_wrap (n, minv, l);
	  fprintf (stderr, "minv : how about in dpotf2() : %d\n", status);
	  if (status != 0)
	    {
	      // CHECK estimate all eigenvalues by dgeev
	      double *wr = (double *)malloc (sizeof (double) * n);
	      double *wi = (double *)malloc (sizeof (double) * n);
	      double *ev = (double *)malloc (sizeof (double) * n * n);
	      dgeev_wrap (n, minv, wr, wi, ev);
	      for (i = 0; i < n; i ++)
		{
		  fprintf (stderr, "# dgeev : w[%d] = %f + i %f\n",
			   i, wr[i], wi[i]);
		}
	      free (wr);
	      free (wi);
	      free (ev);

	      status = BD_sqrt_by_dgeev (n, minv, l);
	      fprintf (stderr,
		       "minv : how about in sqrt_by_dgeev() : %d\n",
		       status);
	      if (status != 0)
		{
		  /* CHECK */
		  for (i = 0; i < BD->sys->np; i ++)
		    {
		      fprintf (stderr, "pos[%d] = %e %e %e\n",
			       i,
			       BD->sys->pos[i*3+0],
			       BD->sys->pos[i*3+1],
			       BD->sys->pos[i*3+2]);
		    }
		  for (i = 0; i < n; i ++)
		    {
		      int j;
		      for (j = 0; j < n; j ++)
			{
			  fprintf (stderr, "minv[%d,%d] = %e\n",
				   i, j,
				   minv [i*n+j]);
			}
		    }
		  exit (1);
		}
	    }
	}

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
    }
  free (y_minv);

  /**
   * lubrication part
   */
  if (BD->flag_lub_B != 0)
    {
      double *y_lub = (double *)malloc (sizeof (double) * n);
      double *z_lub = (double *)malloc (sizeof (double) * n);
      CHECK_MALLOC (y_lub, "BD_calc_FB");
      CHECK_MALLOC (z_lub, "BD_calc_FB");
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
	  CHECK_MALLOC (lub, "BD_calc_FB");
	  CHECK_MALLOC (l,   "BD_calc_FB");
	  BD_matrix_lub_FU (BD, lub);
	  // check
	  status = dpotrf_wrap (n, lub, l);
	  if (status != 0)
	    {
	      check_symmetric (n, lub, 1.0e-12);
	      fprintf (stderr, "lub : failed in dpotrf() : %d\n", status);
	      status = dpotf2_wrap (n, lub, l);
	      fprintf (stderr, "lub : how about in dpotf2() : %d\n", status);
	      if (status != 0)
		{
		  status = BD_sqrt_by_dgeev (n, lub, l);
		  fprintf (stderr,
			   "lub : how about in sqrt_by_dgeev() : %d\n",
			   status);
		  if (status != 0)
		    {
		      /* CHECK */
		      for (i = 0; i < BD->sys->np; i ++)
			{
			  fprintf (stderr, "pos[%d] = %e %e %e\n",
				   i,
				   BD->sys->pos[i*3+0],
				   BD->sys->pos[i*3+1],
				   BD->sys->pos[i*3+2]);
			}
		      for (i = 0; i < n; i ++)
			{
			  int j;
			  for (j = 0; j < n; j ++)
			    {
			      fprintf (stderr, "lub[%d,%d] = %e\n",
				       i, j,
				       lub [i*n+j]);
			    }
			}
		      exit (1);
		    }
		}
	    }

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
      // add z_lub[] into zp[], total random force in the particle-wise order
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
  double dQdt[4];
  int i;
  int nm3 = sys->nm * 3;

  if (sys->version == 0)
    {
      // F version
      for (i = 0; i < nm3; i ++)
	{
	  x[i] += dt * u[i];
	}
    }
  else
    {
      // FT and FTS version
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
	      quaternion_dQdt (q + i4, o + i3, dQdt);

	      q[i4]   += dt * dQdt[0];
	      q[i4+1] += dt * dQdt[1];
	      q[i4+2] += dt * dQdt[2];
	      q[i4+3] += dt * dQdt[3];
	    }
	}
    }
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
      FTS->uf = (double *)malloc (sizeof (double) * nf3);
      CHECK_MALLOC (FTS->ff, "FTS_init");
      CHECK_MALLOC (FTS->uf, "FTS_init");
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
	  FTS->of = (double *)malloc (sizeof (double) * nf3);
	  CHECK_MALLOC (FTS->tf, "FTS_init");
	  CHECK_MALLOC (FTS->of, "FTS_init");
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
	  FTS->ef = (double *)malloc (sizeof (double) * nf5);
	  CHECK_MALLOC (FTS->sf, "FTS_init");
	  CHECK_MALLOC (FTS->ef, "FTS_init");
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
 *                n = nm*3 + nm*4 with quaternion (by BD->flag_Q == 1)
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

      double dt_done;
      switch (BD->scheme)
	{
	case 1: // Banchio-Brady (2003)
	  dt_done = BD_evolve_BB03 (*t, BD, x, q, dt_local);
	  break;

	case 2: // Ball-Melrose (1997)
	  dt_done = BD_evolve_BM97 (*t, BD, x, q, dt_local);
	  break;

	case 0: // the mid-point algorithm
	default:
	  dt_done = BD_evolve_mid (*t, BD, x, q, dt_local);
	  break;
	}

      (*t) += dt_done;
    }
  while ((*t) < t_out);
}

/*
 * INPUT
 *  sys  : struct stokes
 *  dt0  : the initial time step
 *  x0[] : initial configuration
 *  ol   : struct overlap
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

/* return factor for BD->z[]
 */
double
BD_params_get_fact (struct BD_params *BD,
		    double dt)
{
  double fact = 0.0;
  if (BD->peclet > 0.0)
    {
      fact = sqrt(2.0 / (BD->peclet * dt));
    }
  return (fact);
}


/* add interaction forces on each particle including
 *   bond force F^bond(sys->pos[]) by bonds_calc_force()
 *   EV force F^EV(sys->pos[]) by EV_calc_force()
 * this is the position where new forces are added for Brownian Dynamics.
 * INPUT
 *  BD : struct BD_params
 *  pos[] : position for F^bond and F^EV
 * OUTPUT
 *  FTS : struct FTS
 */
void
BD_add_FP (struct BD_params *BD,
	   const double *pos,
	   struct FTS *FTS)
{
  stokes_set_pos_mobile (BD->sys, pos);

  //if (BD->bonds != NULL)
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
  if (BD->ang != NULL)
    {
      // calc angle (bending) force
      angles_calc_force (BD->ang, BD->sys,
			 FTS->f, 1/* add */);
    }
  if (BD->ev_dh != NULL)
    {
      // calc EV_DH force
      EV_DH_calc_force (BD->ev_dh, BD->sys,
			FTS->f, 1/* add */);
    }
  if (BD->ev_LJ != NULL)
    {
      // calc EV_LJ force
      EV_LJ_calc_force (BD->ev_LJ, BD->sys,
			FTS->f, 1/* add */);
    }
  if (BD->cf != NULL)
    {
      // calc confinement force
      CF_calc_force (BD->cf, BD->sys,
		     FTS->f, 1/* add */);
    }
}


/* calculate forces on each particle including
 *   constant force in BDimp->BD->F[],T[]
 *   Brownian force BDimp->fact * BDimp->z[],
 *   interaction forces by BD_add_FP().
 * INPUT
 *  BD : struct BD_params
 *  z[] : random vector obtained by BD_calc_FB()
 *  fact : FB factor for z[] evaluated by BD_params_get_fact()
 * OUTPUT
 *  FTS : struct FTS
 */
static void
BD_calc_forces (struct BD_params *BD,
		double *z, double fact,
		struct FTS *FTS)
{
  int nm3 = BD->sys->nm * 3;
  int nm5 = BD->sys->nm * 5;

  int i;
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

  // interaction forces depending on the configuration by sys->pos[]
  BD_add_FP (BD, BD->sys->pos, FTS);
}


/* evolve position of particles -- the mid-point scheme
 * INPUT
 *  t       : current time
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
 *  returned value : the integrated time duration
 */
double
BD_evolve_mid (double t,
	       struct BD_params *BD,
	       double *x, double *q,
	       double dt)
{
  double dt_local = dt;
  struct overlap ol;

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

  // to keep the initial configuration for re-do.
  double *xmid = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (xmid, "BD_evolve_mid");
  double *qmid = NULL;
  if (BD->sys->version > 0 && q != NULL)
    {
      qmid = (double *)malloc (sizeof (double) * nm4);
      CHECK_MALLOC (qmid, "BD_evolve_mid");
    }

  double fact = BD_params_get_fact (BD, dt_local);

  int i;

 BD_evolve_mid_REDO:
  // set configuration for BD_calc_FB()
  stokes_set_pos_mobile (BD->sys, x);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t,
			  BD->t0, BD->s0);
 BD_evolve_mid_REDO_FB:
  BD_calc_FB (BD, z);
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
  BD_calc_forces (BD, z, fact, FTS);

  // calc dydt
  // BD->sys->pos is already definted above
  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f, FTS->t, FTS->e,
		  BD->uf, BD->of, BD->ef,
		  FTS->u, FTS->o, FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  // evolve by Euler for dt/2 from the initial configuration
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt_local * 0.5,
		     xmid, qmid);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xmid, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# mid-point: overlap in 1st step. "
		   "dt = %e < %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_mid_REDO_FB;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# mid-point: overlap in 1st step. "
		   "dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_mid_REDO_scale;
	}
    }


  /************************************************************
   * 2 : solve for the intermediate (mid-point) configuration *
   ************************************************************/
  // set the intermediate configuration
  stokes_set_pos_mobile (BD->sys, xmid);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t + dt_local * 0.5,
			  BD->t0, BD->s0);

  // set force (and torque) for xmid[]
  BD_calc_forces (BD, z, fact, FTS);

  // calc dydt
  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
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
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt_local,
		     xmid, qmid);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xmid, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos should be reset by the initial config x[] */
	  fprintf (stderr, "# mid-point: overlap in 2nd step. "
		   "dt = %e < %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_mid_REDO;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# mid-point: overlap in 2nd step. "
		   "dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_mid_REDO_scale;
	}
    }

  // constraints
  if (BD->br != NULL)
    {
      /* at this point, 
       * x[]    is the initial configuration (at t) and
       * xmid[] is the final configuration (at t+dt) without constraints
       */
      BeadRod_set_coefs (BD->br, dt_local, 1.0);
      BeadRod_adjust_by_constraints (BD->br,
				     BD->sys->nm,
				     xmid,
				     x);
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
  free (xmid);
  if (qmid != NULL) free (qmid);
  FTS_free (FTS);

  return (dt_local);
}

/* evolve position of particles -- Banchio-Brady scheme
 * reference : Banchio and Brady (2003) Phys. Fluids
 * INPUT
 *  t       : current time
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
 *  returned value : the integrated time duration
 */
double
BD_evolve_BB03 (double t,
		struct BD_params *BD,
		double *x, double *q,
		double dt)
{
  double dt_local = dt;
  struct overlap ol;

  /**
   * memory allocation for working area
   */
  int n = BD_get_n (BD->sys);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "BD_evolve_BB03");

  // fts vectors for local use
  struct FTS *FTS = FTS_init (BD->sys);
  CHECK_MALLOC (FTS, "BD_evolve_BB03");

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
      oBB = (double *)malloc (sizeof (double) * nm3);
      CHECK_MALLOC (oBB, "BD_evolve_BB03");
      if (q != NULL)
	{
	  qBB = (double *)malloc (sizeof (double) * nm4);
	  CHECK_MALLOC (qBB, "BD_evolve_BB03");
	}
    }

  double fact = BD_params_get_fact (BD, dt_local);

  int i;

 BD_evolve_BB03_REDO:
  // (re-)set brownian force
  stokes_set_pos_mobile (BD->sys, x);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t,
			  BD->t0, BD->s0);
 BD_evolve_BB03_REDO_FB:
  BD_calc_FB (BD, z);
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
  // set force only by the Brownian force
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
		  BD->flag_noHI,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  FTS->uf, FTS->of, FTS->ef,
		  uBB,     oBB,     FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  /* 1-2 : evolve by Euler for dt/n with the Brownian velocity
   *       to obtain the intermediate configuration
   */
  evolve_Euler_3all (BD->sys, uBB, oBB, dt_local / BD->BB_n,
		     xBB, qBB);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xBB, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# BB03 : overlap in the intermediate step."
		   " dt = %e < %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_BB03_REDO_FB;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# BB03 : overlap in the intermediate step."
		   " dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_BB03_REDO_scale;
	}
    }

  /* 1-3 : solve the Brownian velocity for the intermediate configuration
   */
  // set the intermediate configuration
  stokes_set_pos_mobile (BD->sys, xBB);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t + dt_local / BD->BB_n,
			  BD->t0, BD->s0);

  // calc dydt with the same Brownian force
  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
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
      for (i = 0; i < nm3; i ++)
	{
	  oBB[i] = oBB[i] + BBfact * (FTS->o[i] - oBB[i]);
	}
    }

  /********************************
   * 2 : solve deterministic part *
   ********************************/
  // set force except for the Brownian force
  if (BD->sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i];
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  FTS->f[i] = BD->F[i];
	}
      for (i = 0; i < nm3; i ++)
	{
	  FTS->t[i] = BD->T[i];
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

  // (re-)set brownian force
  //stokes_set_pos_mobile (BD->sys, x);
  // set in the BD_add_FP() below.

  // (re-)set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t,
			  BD->t0, BD->s0);
  // interaction forces
  BD_add_FP (BD, x, FTS);

  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
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
      for (i = 0; i < nm3; i ++)
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
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt_local,
		     xBB, qBB);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xBB, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos should be reset by the initial config x[] */
	  fprintf (stderr, "# BB03 : overlap in the final step."
		   " dt = %e > %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_BB03_REDO;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# BB03 : overlap in the final step."
		   " dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_BB03_REDO_scale;
	}
    }

  // constraints
  if (BD->br != NULL)
    {
      /* at this point, 
       * x[]    is the initial configuration (at t) and
       * xBB[]  is the final configuration (at t+dt) without constraints
       */
      BeadRod_set_coefs (BD->br, dt_local, 1.0);
      BeadRod_adjust_by_constraints (BD->br,
				     BD->sys->nm,
				     xBB,
				     x);
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

  return (dt_local);
}

/* evolve position of particles -- Ball-Melrose scheme
 * reference : Ball and Melrose (1997)
 * INPUT
 *  t       : current time
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
 *  returned value : the integrated time duration
 */
double
BD_evolve_BM97 (double t,
		struct BD_params *BD,
		double *x, double *q,
		double dt)
{
  double dt_local = dt;
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

  // for the intermediate configurations
  double *xBM = (double *)malloc (sizeof (double) * nm3);
  double *uBM = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (xBM, "BD_evolve_BM97");
  CHECK_MALLOC (uBM, "BD_evolve_BM97");
  double *qBM = NULL;
  double *oBM = NULL;
  if (BD->sys->version > 0)
    {
      oBM = (double *)malloc (sizeof (double) * nm3);
      CHECK_MALLOC (oBM, "BD_evolve_BM97");
      if (q != NULL)
	{
	  qBM = (double *)malloc (sizeof (double) * nm4);
	  CHECK_MALLOC (qBM, "BD_evolve_BM97");
	}
    }

  double fact = BD_params_get_fact (BD, dt_local);

  int i;


 BD_evolve_BM97_REDO:
  // (re-)set brownian force
  stokes_set_pos_mobile (BD->sys, x);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t,
			  BD->t0, BD->s0);
 BD_evolve_BM97_REDO_FB:
  BD_calc_FB (BD, z);
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
  BD_calc_forces (BD, z, fact, FTS);

  // calc dydt
  // BD->sys->pos is already definted above
  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
		  BD->flag_lub, BD->flag_mat,
		  FTS->f,  FTS->t,  FTS->e,
		  BD->uf,  BD->of,  BD->ef,
		  uBM,     oBM,     FTS->s,
		  FTS->ff, FTS->tf, FTS->sf);

  // evolve by Euler for dt from the initial configuration
  evolve_Euler_3all (BD->sys, uBM, oBM, dt_local,
		     xBM, qBM);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xBM, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# BM97: overlap in 1st step. "
		   "dt = %e > %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_BM97_REDO_FB;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# BM97: overlap in 1st step. "
		   "dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_BM97_REDO_scale;
	}
    }

  /************************************************************
   * 2 : solve for the intermediate (mid-point) configuration *
   ************************************************************/
  // set the intermediate configuration
  stokes_set_pos_mobile (BD->sys, xBM);
  // set sys->shear_shift if necessary
  stokes_set_shear_shift (BD->sys, t + dt_local,
			  BD->t0, BD->s0);

  // set force (and torque) for new position xBM[]
  BD_calc_forces (BD, z, fact, FTS);

  solve_mix_3all (BD->sys,
		  BD->flag_noHI,
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
      for (i = 0; i < nm3; i ++)
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
  evolve_Euler_3all (BD->sys, FTS->u, FTS->o, dt_local,
		     xBM, qBM);

  /* dt-ajustment process -- overlap check */
  if (BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BD->sys, xBM, BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos should be reset by the initial config x[] */
	  fprintf (stderr, "# BM97: overlap in 2nd step. "
		   "dt = %e < %e => reconstruct FB\n", dt_local, BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_BM97_REDO;
	}
      else
	{
	  fact = BD_params_get_fact (BD, dt_local);

	  fprintf (stderr, "# BM97: overlap in 2nd step. "
		   "dt = %e > %e\n", dt_local, BD->dt_lim);
	  goto BD_evolve_BM97_REDO_scale;
	}
    }

  // constraints
  if (BD->br != NULL)
    {
      /* at this point, 
       * x[]    is the initial configuration (at t) and
       * xBM[]  is the final configuration (at t+dt) without constraints
       */
      BeadRod_set_coefs (BD->br, dt_local, 1.0);
      BeadRod_adjust_by_constraints (BD->br,
				     BD->sys->nm,
				     xBM,
				     x);
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

  return (dt_local);
}
