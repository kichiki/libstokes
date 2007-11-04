/* test code for brownian.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-brownian.c,v 1.2 2007/11/04 00:17:08 kichiki Exp $
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
#include "memory-check.h"
#include "check.h" // compare()

#include <dnaupd_c.h> // dnaupd_wrap_min_max()
#include <chebyshev.h>
#include <dgetri_c.h> // lapack_inv_()


#include <dpotrf_c.h> // dpotrf_wrap()


#include <stokes.h>
#include <brownian.h>
#include <ewald-3fts-matrix.h>
#include <ewald.h> // make_matrix_mob_3all()
#include <matrix.h> // multiply_extmat_with_extvec_3fts()

#include <bench.h> // ptime_ms_d()


static void
atimes_by_matrix (int n, const double *x, double *b, void *user_data)
{
  double *a = (double *)user_data;

  int i, j;
  for (i = 0; i < n; i ++)
    {
      b[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  b[i] += a[i*n+j] * x[j];
	}
    }
}

static double
inv_sqrt (double x)
{
  return (1.0 / sqrt (x));
}

/* compare two process to obtain vector z[i]
 * which satisfies z.z = y.(A^{-1}).y
 */
int
check_cheb_minv (int n,
		 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_cheb_minv : start\n");
      fprintf (stdout, "# z[i] satisfy z.z = y.(A^{-1}).y\n");
    }

  int check = 0;


  double err_minv;
  /**
   * initialization
   */
  double *y = (double *)malloc (sizeof (double) * n);
  double *z_mat = (double *)malloc (sizeof (double) * n);
  double *a = (double *)malloc (sizeof (double) * n * n);
  double *ainv = (double *)malloc (sizeof (double) * n * n);
  double *l = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (y, "check_cheb_minv");
  CHECK_MALLOC (z_mat, "check_cheb_minv");
  CHECK_MALLOC (a, "check_cheb_minv");
  CHECK_MALLOC (ainv, "check_cheb_minv");
  CHECK_MALLOC (l, "check_cheb_minv");

  // set the matrix 'a'
  int i, j;
  for (i = 0; i < n; i ++)
    {
      a[i*n+i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  a[i*n+j] = 1.0 / (double)(i+j+1);
	  a[j*n+i] = a[i*n+j];
	}
    }
  // set the vector 'y'
  for (i = 0; i < n; i ++)
    {
      y[i] = 1.0;
    }

  /**
   * matrix version
   */
  lapack_inv (n, a, ainv);
  dpotrf_wrap (n, ainv, l);
  for (i = 0; i < n; i ++)
    {
      z_mat[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  z_mat[i] += l[j*n+i] * y[j];
	}
    }
  /* note that here M^-1 = L . L^T, so that 
   * z_i = y_j L_{ji} has the property that
   * z . z = y_j L_{ji} y_k L_{ki} = M^-1_{jk} y_j y_k 
   * therefore, use not chebyshev_error_minvsqrt() but chebyshev_error_Rsqrt() 
   */
  err_minv = chebyshev_error_Rsqrt (n, y, z_mat,
				    atimes_by_matrix, (void *)ainv);
  if (verbose != 0)
    {
      fprintf (stderr, "# matrix : err_minv = %e\n", err_minv);
    }

  /**
   * atimes version
   */
  double *z_ax = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z_ax, "check_cheb_minv");

  int n_minv = 50;
  double *a_minv = (double *)malloc (sizeof (double) * n_minv);
  CHECK_MALLOC (a_minv, "check_cheb_minv");

  double eig[2];
  dnaupd_wrap_min_max (n, eig,
		       atimes_by_matrix, (void *)a, 1.0e-6);
  chebyshev_coef (n_minv, inv_sqrt,
		  eig[0], eig[1], a_minv);
  chebyshev_eval_atimes (n_minv, a_minv,
			 n, y, z_ax,
			 eig[0], eig[1],
			 atimes_by_matrix, (void *)a);
  err_minv = chebyshev_error_minvsqrt (n, y, z_ax,
				       atimes_by_matrix, (void *)a);
  if (verbose != 0)
    {
      fprintf (stderr, "# atimes : err_minv = %e\n", err_minv);
    }

  double zz_mat = 0.0;
  double zz_ax  = 0.0;
  for (i = 0; i < n; i ++)
    {
      zz_mat += z_mat[i] * z_mat[i];
      zz_ax  += z_ax[i]  * z_ax[i];
    }
  check += compare (zz_mat, zz_ax, "# zz(mat) vs zz(ax)", verbose, tiny);


  free (y);
  free (z_mat);
  free (a);
  free (ainv);
  free (l);

  free (z_ax);
  free (a_minv);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* compare two process to obtain vector z[i]
 * which satisfies z.z = y.A.y
 */
int
check_cheb_lub (int n,
		int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_cheb_lub : start\n");
      fprintf (stdout, "# z[i] satisfy z.z = y.A.y\n");
    }

  int check = 0;


  double err_lub;
  /**
   * initialization
   */
  double *y = (double *)malloc (sizeof (double) * n);
  double *z_mat = (double *)malloc (sizeof (double) * n);
  double *a = (double *)malloc (sizeof (double) * n * n);
  double *l = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (y, "check_cheb_lub");
  CHECK_MALLOC (z_mat, "check_cheb_lub");
  CHECK_MALLOC (a, "check_cheb_lub");
  CHECK_MALLOC (l, "check_cheb_lub");

  // set the matrix 'a'
  int i, j;
  for (i = 0; i < n; i ++)
    {
      a[i*n+i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  a[i*n+j] = 1.0 / (double)(i+j+1);
	  a[j*n+i] = a[i*n+j];
	}
    }
  // set the vector 'y'
  for (i = 0; i < n; i ++)
    {
      y[i] = 1.0;
    }

  /**
   * matrix version
   */
  dpotrf_wrap (n, a, l);
  for (i = 0; i < n; i ++)
    {
      z_mat[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  z_mat[i] += l[j*n+i] * y[j];
	}
    }
  err_lub = chebyshev_error_Rsqrt (n, y, z_mat,
				   atimes_by_matrix, (void *)a);
  if (verbose != 0)
    {
      fprintf (stderr, "# matrix : err_lub = %e\n", err_lub);
    }

  /**
   * atimes version
   */
  double *z_ax = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z_ax, "check_cheb_lub");

  int n_lub = 50;
  double *a_lub = (double *)malloc (sizeof (double) * n_lub);
  CHECK_MALLOC (a_lub, "check_cheb_lub");

  double eig[2];
  dnaupd_wrap_min_max (n, eig,
		       atimes_by_matrix, (void *)a, 1.0e-6);
  chebyshev_coef (n_lub, sqrt,
		  eig[0], eig[1], a_lub);
  chebyshev_eval_atimes (n_lub, a_lub,
			 n, y, z_ax,
			 eig[0], eig[1],
			 atimes_by_matrix, (void *)a);
  err_lub = chebyshev_error_Rsqrt (n, y, z_ax,
				   atimes_by_matrix, (void *)a);
  if (verbose != 0)
    {
      fprintf (stderr, "# atimes : err_lub = %e\n", err_lub);
    }

  double zz_mat = 0.0;
  double zz_ax  = 0.0;
  for (i = 0; i < n; i ++)
    {
      zz_mat += z_mat[i] * z_mat[i];
      zz_ax  += z_ax[i]  * z_ax[i];
    }
  check += compare (zz_mat, zz_ax, "# zz(mat) vs zz(ax)", verbose, tiny);


  free (y);
  free (z_mat);
  free (a);
  free (l);

  free (z_ax);
  free (a_lub);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* compare minv routines in the matrix and atimes versions
 *  BD_matrix_minv_FU() and BT_atimes_mob_FU()
 */
int
check_minv_FU (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_minv_FU : start\n");
      fprintf (stdout, "# M^{-1}_{FU} by matrix and atimes routines\n");
    }

  int check = 0;


  /**
   * initialization
   */
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_minv_FU");

  sys->version = 2; // FTS
  int np = 9;
  stokes_set_np (sys, np, np);

  double *pos = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (pos, "check_minv_FU");
  pos[ 0] = 4.841648; pos[ 1] = 4.345871; pos[ 2] = 2.411424;
  pos[ 3] = 1.484547; pos[ 4] = 2.931487; pos[ 5] = 2.321764;
  pos[ 6] = 4.494053; pos[ 7] = 1.739099; pos[ 8] = 2.732656;
  pos[ 9] = 1.436324; pos[10] = 5.216197; pos[11] = 5.634981;
  pos[12] = 2.707117; pos[13] = 0.577389; pos[14] = 2.489404;
  pos[15] = 3.502632; pos[16] = 5.234916; pos[17] = 5.542215;
  pos[18] = 4.615775; pos[19] = 3.146068; pos[20] = 0.080511;
  pos[21] = 5.510249; pos[22] = 1.060571; pos[23] = 0.710376;
  pos[24] = 1.143043; pos[25] = 2.055843; pos[26] = 4.194753;
  stokes_set_pos (sys, pos);
  free (pos);
  
  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  sys->periodic = 1;
  stokes_set_l (sys, 5.733683, 5.733683, 5.733683);
  double xi = xi_by_tratio (sys, 1.0);
  stokes_set_xi (sys, xi, 1.0e-6);

  /**
   * 1) atimes calculation of M_{FU}
   *    x[n] = M_{FU} . y[n],
   *    where y[n] is random vector.
   */
  int n = np * 6;

  struct BD_params *BD = BD_params_init
    (sys,
     NULL, // struct KIrand *rng,
     NULL, // double *pos_fixed,
     NULL, // double *F,
     NULL, // double *T,
     NULL, // double *E,
     NULL, // double *uf,
     NULL, // double *of,
     NULL, // double *ef,
     0,    // int flag_lub,
     1,    // int flag_mat,
     0.0,  // double st,
     NULL, // struct bonds *bonds,
     0.0,  // double gamma,
     NULL, // struct EV *ev,
     0,    // int flag_Q,
     -1.0, // double peclet,
     0.0,  // double eps,
     0,    // int    n_minv,
     0,    // int    n_lub,
     0,    // int    scheme,
     0.0   // double BB_n
     );
  CHECK_MALLOC (BD, "check_minv_FU");

  double *x = (double *)malloc (sizeof (double) * n);
  double *y = (double *)malloc (sizeof (double) * n);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_minv_FU");
  CHECK_MALLOC (y, "check_minv_FU");
  CHECK_MALLOC (z, "check_minv_FU");

  int i;
  srand48(0);
  for (i = 0; i < n; i ++)
    {
      y[i] = drand48();
    }
  BD_atimes_mob_FU (n, y, x, (void *)BD);


  /**
   * 2) matrix version of M_{FU} by solve_mob_3fts_matrix()
   *    z[n] = M_{FU} . y[n]
   */
  // compare with solve_mob_3fts_matrix() in ewlad-3fts-matrix.c
  int np3 = np * 3;
  int np5 = np * 5;
  double *f = (double *)malloc (sizeof (double) * np3);
  double *t = (double *)malloc (sizeof (double) * np3);
  double *s = (double *)malloc (sizeof (double) * np5);
  double *u = (double *)malloc (sizeof (double) * np3);
  double *o = (double *)malloc (sizeof (double) * np3);
  double *e = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (f, "check_minv_FU");
  CHECK_MALLOC (t, "check_minv_FU");
  CHECK_MALLOC (s, "check_minv_FU");
  CHECK_MALLOC (u, "check_minv_FU");
  CHECK_MALLOC (o, "check_minv_FU");
  CHECK_MALLOC (e, "check_minv_FU");
  int ii, jj;
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  f[i*3+ii] = y[i*6+  ii];
	  t[i*3+ii] = y[i*6+3+ii];
	}
      for (ii = 0; ii < 5; ii ++)
	{
	  // e = 0
	  e[i*5+ii] = 0.0;
	}
    }

  solve_mob_3fts_matrix (sys, f, t, e, u, o, s);
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  z[i*6+  ii] = u[i*3+ii];
	  z[i*6+3+ii] = o[i*3+ii];
	}
    }

  fprintf (stdout, "# 1) & 2) BD_atimes_mob_FU() vs solve_mob_3fts_matrix()\n");
  char label[80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# %d", i);
      check += compare (z[i], x[i], label, verbose, tiny);
    }

  /**
   * 3) matrix version of (M_{FU})^{-1} by solve_res_3fts_matrix_0()
   * z0[n] = R_{FU} . x[n]
   * because x[n] = M_{FU} . y[n], z0[n] == y[n].
   */
  // compare with solve_res_3fts_matrix() in ewlad-3fts-matrix.c
  // z0 = R . x  == R . M . y = y
  double *z0 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z0, "check_minv_FU");
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  u[i*3+ii] = x[i*6+ii];
	  o[i*3+ii] = x[i*6+3+ii];
	}
      for (ii = 0; ii < 5; ii ++)
	{
	  // e = 0
	  e[i*5+ii] = 0.0;
	}
    }

  solve_res_3fts_matrix_0 (sys, u, o, e, f, t, s);

  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  z0[i*6+  ii] = f[i*3+ii];
	  z0[i*6+3+ii] = t[i*3+ii];
	}
    }

  fprintf (stdout, "# 3) compare for (M_FU)^{-1}"
	   " by solve_res_3fts_matrix_0()\n");
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# solve_res_3fts_matrix_0() : %d", i);
      check += compare (z0[i], y[i], label, verbose, tiny);
    }

  free (f);
  free (t);
  free (s);
  free (u);
  free (o);
  free (e);

  /**
   * 4) matrix make_matrix_mob_3all() -> lapack_inv_()
   * 
   */
  int np11 = np * 11;
  double *m = (double *)malloc (sizeof (double) * np11 * np11);
  CHECK_MALLOC (m, "check_minv_FU");
  make_matrix_mob_3all (BD->sys, m);
  lapack_inv_ (np11, m);
  // m is the resistance matrix in the INVERSE form
  trans_ext (np, m);
  // now m is the resistance matrix in the EXTRACTED form

  double *uoe = (double *)malloc (sizeof (double) * np11);
  double *fts = (double *)malloc (sizeof (double) * np11);
  CHECK_MALLOC (uoe, "check_minv_FU");
  CHECK_MALLOC (fts, "check_minv_FU");
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 6; ii ++)
	{
	  uoe[i*11+ii] = x[i*6+ii];
	}
      for (ii = 6; ii < 11; ii ++)
	{
	  // e = 0
	  uoe[i*11+ii] = 0.0;
	}
    }
  multiply_extmat_with_extvec_3fts (np, m, uoe, fts);
  // y = R(ext) . x in FTS version

  double *z2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z2, "check_minv_FU");
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 6; ii ++)
	{
	  z2[i*6+ii] = fts[i*11+ii];
	}
    }
  fprintf (stdout, "# 4) compare with direct calc by"
	   " make_matrix_mob_3all() -> lapack_inv_()\n");
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# z0 and z2 : %d", i);
      check += compare (z0[i], z2[i], label, verbose, tiny);
    }

  free (z2);
  free (uoe);
  free (fts);

  // extract FU part from m[] to minv0[]
  double *minv0 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (minv0, "check_minv_FU");

  int j;
  for (i = 0; i < np; i ++)
    {
      for (j = 0; j < np; j ++)
	{
	  for (ii = 0; ii < 6; ii ++)
	    {
	      for (jj = 0; jj < 6; jj ++)
		{
		  minv0 [(i*6+ii)*n+(j*6+jj)]
		    = m[(i*11+ii)*np11+(j*11+jj)];
		}
	    }
	}
    }

  double *z1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z1, "check_minv_FU");
  for (i = 0; i < n; i ++)
    {
      z1[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  z1[i] += minv0[i*n+j] * x[j];
	}
    }

  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# z0 and z1 : %d", i);
      check += compare (z0[i], z1[i], label, verbose, tiny);
    }

  /**
   * 5) matrix BD_matrix_minv_FU()
   */
  double *minv = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (minv, "check_minv_FU");

  BD_matrix_minv_FU (BD, minv);

  // first compare the matrices
  fprintf (stdout, "# 5) matrix BD_matrix_minv_FU()\n");
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# minv [%d, %d]", i, j);
	  check += compare (minv0[i*n+j], minv[i*n+j], label, verbose, tiny);
	}
    }


  for (i = 0; i < n; i ++)
    {
      z[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  z[i] += minv[i*n+j] * x[j];
	}
    }
  // z = M^{-1} . x should be equal to y

  fprintf (stdout, "# 6) Compare z and y\n");
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# %d", i);
      check += compare (z[i], y[i], label, verbose, tiny);
    }

  fprintf (stdout, "# 7) Compare z0 and z\n");
  // comparison between solve_... and BD_...
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# solve_... vs BD_... : %d", i);
      check += compare (z0[i], z[i], label, verbose, tiny);
    }

  free (m);
  free (minv);
  free (minv0);
  free (z0);
  free (z1);
  free (x);
  free (y);
  free (z);
  BD_params_free (BD);
  stokes_free (sys);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

/* compare lub routines in the matrix and atimes versions
 *  BD_matrix_lub_FU() and BT_atimes_lub_FU()
 */
int
check_lub_FU (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_lub_FU : start\n");
      fprintf (stdout, "# L_{FU} by matrix and atimes routines\n");
    }

  int check = 0;


  /**
   * initialization
   */
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_minv_FU");

  sys->version = 2; // FTS
  int np = 9;
  stokes_set_np (sys, np, np);

  double *pos = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (pos, "check_minv_FU");
  pos[ 0] = 4.841648; pos[ 1] = 4.345871; pos[ 2] = 2.411424;
  pos[ 3] = 1.484547; pos[ 4] = 2.931487; pos[ 5] = 2.321764;
  pos[ 6] = 4.494053; pos[ 7] = 1.739099; pos[ 8] = 2.732656;
  pos[ 9] = 1.436324; pos[10] = 5.216197; pos[11] = 5.634981;
  pos[12] = 2.707117; pos[13] = 0.577389; pos[14] = 2.489404;
  pos[15] = 3.502632; pos[16] = 5.234916; pos[17] = 5.542215;
  pos[18] = 4.615775; pos[19] = 3.146068; pos[20] = 0.080511;
  pos[21] = 5.510249; pos[22] = 1.060571; pos[23] = 0.710376;
  pos[24] = 1.143043; pos[25] = 2.055843; pos[26] = 4.194753;
  stokes_set_pos (sys, pos);
  free (pos);
  
  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  sys->periodic = 1;
  stokes_set_l (sys, 5.733683, 5.733683, 5.733683);
  double xi = xi_by_tratio (sys, 1.0);
  stokes_set_xi (sys, xi, 1.0e-6);

  struct BD_params *BD = BD_params_init
    (sys,
     NULL, // struct KIrand *rng,
     NULL, // double *pos_fixed,
     NULL, // double *F,
     NULL, // double *T,
     NULL, // double *E,
     NULL, // double *uf,
     NULL, // double *of,
     NULL, // double *ef,
     1,    // int flag_lub,
     1,    // int flag_mat,
     0.0,  // double st,
     NULL, // struct bonds *bonds,
     0.0,  // double gamma,
     NULL, // struct EV *ev,
     0,    // int flag_Q,
     -1.0, // double peclet,
     0.0,  // double eps,
     0,    // int    n_minv,
     0,    // int    n_lub,
     0,    // int    scheme,
     0.0   // double BB_n
     );
  CHECK_MALLOC (BD, "check_lub_FU");

  int n = np * 6;
  double *x = (double *)malloc (sizeof (double) * n);
  double *y = (double *)malloc (sizeof (double) * n);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_lub_FU");
  CHECK_MALLOC (y, "check_lub_FU");
  CHECK_MALLOC (z, "check_lub_FU");

  int i;
  srand48(0);
  for (i = 0; i < n; i ++)
    {
      y[i] = drand48();
    }

  /**
   * 1) atimes calculation of L_{FU}
   *    x[n] = L_{FU} . y[n],
   *    where y[n] is random vector.
   */
  BD_atimes_lub_FU (n, y, x, (void *)BD);


  /**
   * 2) matrix version of L_{FU}
   *    z[n] = L_{FU} . y[n]
   */
  double *lub = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (lub, "check_lub_FU");

  BD_matrix_lub_FU (BD, lub);

  int j;
  for (i = 0; i < n; i ++)
    {
      z[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  z[i] += lub[i*n+j] * y[j];
	}
    }


  char label[80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "# %d", i);
      check += compare (z[i], x[i], label, verbose, tiny);
    }


  free (lub);
  free (x);
  free (y);
  free (z);
  BD_params_free (BD);
  stokes_free (sys);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* speed test on BD_minv_FU_in_FTS()
 */
void
BD_minv_FU_in_FTS_ (int np, const double *m, double *minv_FU)
{
  int np5  = np * 5;
  int np6  = np * 6;
  double *a = (double *)malloc (sizeof (double) * np6 * np6);
  double *b = (double *)malloc (sizeof (double) * np6 * np5);
  double *c = (double *)malloc (sizeof (double) * np5 * np6);
  double *d = (double *)malloc (sizeof (double) * np6 * np6);
  double *x = (double *)malloc (sizeof (double) * np5 * np6);
  double *y = (double *)malloc (sizeof (double) * np6 * np6);
  CHECK_MALLOC (a, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (b, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (c, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (d, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (x, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (y, "BD_minv_FU_in_FTS");

  // decompose m into a,b,c,d
  split_matrix_3fts (np, m, a, b, c, d);

  // d = d^{-1}
  lapack_inv_ (np5, d);

  // x[np5,np6] = d[np5,np5] * c[np5,np6]
  /*
  mul_matrices (d, np5, np5,
		c, np5, np6,
		x);
  */
  int i, j;
  int k;
  for (i = 0; i < np5; i ++)
    {
      for (j = 0; j < np6; j ++)
	{
	  x[i*np6+j] = 0.0;
	  for (k = 0; k < np5; k ++)
	    {
	      x[i*np6+j] += d[i*np5+k] * c[k*np6+j];
	    }
	}
    }

  // y[np6,np6] = b[np6,np5] * x[np5,np6]
  /*
  mul_matrices (b, np6, np5,
		x, np5, np6,
		y);
  */
  for (i = 0; i < np6; i ++)
    {
      for (j = 0; j < np6; j ++)
	{
	  y[i*np6+j] = 0.0;
	  for (k = 0; k < np5; k ++)
	    {
	      y[i*np6+j] += b[i*np5+k] * x[k*np6+j];
	    }
	}
    }
  // a = a - b.d^{-1}.c
  for (i = 0; i < np6 * np6; i ++)
    {
      minv_FU [i] = a[i] - y[i];
    }
  lapack_inv_ (np6, minv_FU);

  free (a);
  free (b);
  free (c);
  free (d);
  free (x);
  free (y);
}

int
benchmark_BD_minv_FU_in_FTS (int np, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "benchmark_BD_minv_FU_in_FTS : start\n");
      fprintf (stdout, "# benchmark test on benchmark_BD_minv_FU_in_FTS()\n");
      fprintf (stdout, "# np, t1, t2, where\n"
	       "# t1 : with mul_matrices in libstokes\n"
	       "# t2 : with double for loop\n");
    }

  int check = 0;


  /**
   * initialization
   */
  int n11 = np * 11;
  int n6  = np * 6;

  double *a = (double *)malloc (sizeof (double) * n11 * n11);
  double *ainv1 = (double *)malloc (sizeof (double) * n6 * n6);
  double *ainv2 = (double *)malloc (sizeof (double) * n6 * n6);
  CHECK_MALLOC (a, "benchmark_BD_minv_FU_in_FTS");
  CHECK_MALLOC (ainv1, "benchmark_BD_minv_FU_in_FTS");
  CHECK_MALLOC (ainv2, "benchmark_BD_minv_FU_in_FTS");

  // set the matrix 'a'
  int i, j;
  srand48(0);
  for (i = 0; i < n11; i ++)
    {
      a[i*n11+i] = 1.0;
      for (j = i+1; j < n11; j ++)
	{
	  a[i*n11+j] = drand48();
	  a[j*n11+i] = a[i*n11+j];
	}
    }

  // take an average for 10 calculations
  double t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      BD_minv_FU_in_FTS (np, a, ainv1);
    }
  double t1 = ptime_ms_d ();
  t1 -= t0;
  t1 *= 0.1;


  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      BD_minv_FU_in_FTS_ (np, a, ainv2);
    }
  double t2 = ptime_ms_d ();
  t2 -= t0;
  t2 *= 0.1;

  fprintf (stdout, "%d %f %f\n", np, t1, t2);

  char label[80];
  for (i = 0; i < n6; i ++)
    {
      for (j = 0; j < n6; j ++)
	{
	  sprintf (label, "# ainv[%d %d]", i, j);
	  check += compare (ainv1[i*n6+j], ainv2[i*n6+j], label, verbose, tiny);
	}
    }

  free (a);
  free (ainv1);
  free (ainv2);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}



/* check for inversion by sub-matrices
 */
void
solve_ainv_1 (int n1, int n2,
	      const double *a, const double *b,
	      const double *c, const double *d,
	      double *ainv)
{
  int np6  = n1;
  int np5  = n2;
  double *dinv = (double *)malloc (sizeof (double) * np5 * np5);
  double *x = (double *)malloc (sizeof (double) * np5 * np6);
  double *y = (double *)malloc (sizeof (double) * np6 * np6);
  CHECK_MALLOC (dinv, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (x, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (y, "BD_minv_FU_in_FTS");

  // d = d^{-1}
  lapack_inv (np5, d, dinv);

  // x[np5,np6] = d[np5,np5] * c[np5,np6]
  mul_matrices (d, np5, np5,
		c, np5, np6,
		x);

  // y[np6,np6] = b[np6,np5] * x[np5,np6]
  mul_matrices (b, np6, np5,
		x, np5, np6,
		y);

  // a = a - b.d^{-1}.c
  int i;
  for (i = 0; i < np6 * np6; i ++)
    {
      ainv [i] = a[i] - y[i];
    }
  lapack_inv_ (np6, ainv);

  free (dinv);
  free (x);
  free (y);
}
void
solve_ainv_2 (int n1, int n2,
	      const double *a, const double *b,
	      const double *c, const double *d,
	      double *ainv)
{
  int np6  = n1;
  int np5  = n2;
  double *dinv = (double *)malloc (sizeof (double) * np5 * np5);
  double *x = (double *)malloc (sizeof (double) * np5 * np6);
  double *y = (double *)malloc (sizeof (double) * np6 * np6);
  CHECK_MALLOC (dinv, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (x, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (y, "BD_minv_FU_in_FTS");

  // d = d^{-1}
  lapack_inv (np5, d, dinv);

  // x[np5,np6] = d[np5,np5] * c[np5,np6]
  int i, j, k;
  for (i = 0; i < np5; i ++)
    {
      for (j = 0; j < np6; j ++)
	{
	  x[i*np6+j] = 0.0;
	  for (k = 0; k < np5; k ++)
	    {
	      x[i*np6+j] += d[i*np5+k] * c[k*np6+j];
	    }
	}
    }

  // y[np6,np6] = b[np6,np5] * x[np5,np6]
  for (i = 0; i < np6; i ++)
    {
      for (j = 0; j < np6; j ++)
	{
	  y[i*np6+j] = 0.0;
	  for (k = 0; k < np5; k ++)
	    {
	      y[i*np6+j] += b[i*np5+k] * x[k*np6+j];
	    }
	}
    }
  // a = a - b.d^{-1}.c
  for (i = 0; i < np6 * np6; i ++)
    {
      ainv [i] = a[i] - y[i];
    }
  lapack_inv_ (np6, ainv);

  free (dinv);
  free (x);
  free (y);
}
int
check_inv_by_submatrices (int n1, int n2, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_inv_by_submatrices : start\n");
    }

  int check = 0;

  /**
   * initialization
   */
  int n = n1 + n2;
  double *m = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (m, "check_inv_by_submatrices");

  // set the matrix 'a'
  int i, j;
  srand48(0);
  for (i = 0; i < n; i ++)
    {
      m[i*n+i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  m[i*n+j] = drand48();
	  m[j*n+i] = m[i*n+j];
	}
    }

  // simply split into 4 parts
  double *a = (double *)malloc (sizeof (double) * n1 * n1);
  double *b = (double *)malloc (sizeof (double) * n1 * n2);
  double *c = (double *)malloc (sizeof (double) * n2 * n1);
  double *d = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (a, "check_inv_by_submatrices");
  CHECK_MALLOC (b, "check_inv_by_submatrices");
  CHECK_MALLOC (c, "check_inv_by_submatrices");
  CHECK_MALLOC (d, "check_inv_by_submatrices");
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  a[i*n1+j] = m[i*n+j];
	}
      for (j = n1; j < n; j ++)
	{
	  b[i*n2+(j-n1)] = m[i*n+j];
	}
    }
  for (i = n1; i < n; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  c[(i-n1)*n1+j] = m[i*n+j];
	}
      for (j = n1; j < n; j ++)
	{
	  d[(i-n1)*n2+(j-n1)] = m[i*n+j];
	}
    }
  free (m);


  /**
   * locally compare step by step
   */
  double *dinv = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (dinv, "check_inv_by_submatrices");
  lapack_inv (n2, d, dinv);

  double *x1 = (double *)malloc (sizeof (double) * n2 * n1);
  double *x2 = (double *)malloc (sizeof (double) * n2 * n1);
  CHECK_MALLOC (x1, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (x2, "BD_minv_FU_in_FTS");

  // x[np5,np6] = d[np5,np5] * c[np5,np6]
  mul_matrices (d, n2, n2,
		c, n2, n1,
		x1);
  int k;
  for (i = 0; i < n2; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  x2[i*n1+j] = 0.0;
	  for (k = 0; k < n2; k ++)
	    {
	      x2[i*n1+j] += d[i*n2+k] * c[k*n1+j];
	    }
	}
    }
  char label[80];
  for (i = 0; i < n2; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "# x1,x2[%d %d]", i, j);
	  check += compare (x1[i*n1+j], x2[i*n1+j], label, verbose, tiny);
	}
    }

  double *y1 = (double *)malloc (sizeof (double) * n1 * n1);
  double *y2 = (double *)malloc (sizeof (double) * n1 * n1);
  CHECK_MALLOC (y1, "BD_minv_FU_in_FTS");
  CHECK_MALLOC (y2, "BD_minv_FU_in_FTS");
  // y[np6,np6] = b[np6,np5] * x[np5,np6]
  mul_matrices (b, n2, n2,
		x1, n2, n1,
		y1);
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  y2[i*n1+j] = 0.0;
	  for (k = 0; k < n2; k ++)
	    {
	      y2[i*n1+j] += b[i*n2+k] * x2[k*n1+j];
	    }
	}
    }
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "# y1,y2[%d %d]", i, j);
	  check += compare (y1[i*n1+j], y2[i*n1+j], label, verbose, tiny);
	}
    }

  double *ainv1 = (double *)malloc (sizeof (double) * n1 * n1);
  double *ainv2 = (double *)malloc (sizeof (double) * n1 * n1);
  CHECK_MALLOC (ainv1, "check_inv_by_submatrices");
  CHECK_MALLOC (ainv2, "check_inv_by_submatrices");
  for (i = 0; i < n1 * n1; i ++)
    {
      ainv1 [i] = a[i] - y1[i];
      ainv2 [i] = a[i] - y2[i];
    }

  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "# ainv1,ainv2[%d %d]", i, j);
	  check += compare (ainv1[i*n1+j], ainv2[i*n1+j], label, verbose, tiny);
	}
    }

  lapack_inv_ (n1, ainv1);
  lapack_inv_ (n1, ainv2);

  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "# after inv - ainv1,ainv2[%d %d]", i, j);
	  check += compare (ainv1[i*n1+j], ainv2[i*n1+j], label, verbose, tiny);
	}
    }

  free (x1);
  free (x2);
  free (y1);
  free (y2);
  free (dinv);

  /**
   * top-level comparison
   */
  solve_ainv_1 (n1, n2, a, b, c, d, ainv1);
  solve_ainv_2 (n1, n2, a, b, c, d, ainv2);

  free (a);
  free (b);
  free (c);
  free (d);

  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "# ainv[%d %d]", i, j);
	  check += compare (ainv1[i*n1+j], ainv2[i*n1+j], label, verbose, tiny);
	}
    }
  free (ainv1);
  free (ainv2);

  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
