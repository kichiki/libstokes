/* test code for brownian.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-brownian.c,v 1.1 2007/10/27 23:09:31 kichiki Exp $
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
