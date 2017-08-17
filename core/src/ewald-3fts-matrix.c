/* Solvers for 3 dimensional FTS version problems by MATRIX procedure
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts-matrix.c,v 2.18 2007/12/22 17:36:56 kichiki Exp $
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "memory-check.h" // CHECK_MALLOC

#include "dgetri_c.h" /* lapack_inv_() */

#include "stokes.h"
#include "bench.h"
#include "fts.h"
#include "ft.h"
#include "f.h"

#include "ewald.h" // make_matrix_mob_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3fts-matrix.h"


/** natural resistance problem **/
/* this is just a test routine */
static void
test_symmetric (int n, const double * mat, double tiny)
{
  int i, j;
  double d;


  for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
	{
	  d = fabs (mat [i * n + j] - mat [j * n + i]);
	  if (d > tiny)
	    fprintf (stderr, "mat [%d, %d] != mat [%d, %d], "
		     "|%f - %f| = %e\n",
		     i, j, j, i,
		     mat [i * n + j], mat [j * n + i],
		     d);
		     
	}
    }
}
/** copy from NR/src/FINITE/stokes-fts.c Rev 1.4 **/
/*
 * INPUT
 *  r [np * 11 * np * 11] : this is INVERSED form
 * OUTPUT
 *  r [np * 11 * np * 11] : this is EXTRACTED form
 */
void
trans_ext (int np, double *r)
{
  static double
    tinv [121] = {
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  2.0/3.0, 0.0, 0.0, 0.0, -1.0/3.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.5, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.5, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.5,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/3.0, 0.0, 0.0, 0.0, 2.0/3.0};

  int i, ii, i11, i0;
  int j, jj, j11, j0;
  int k;

  int n = np * 11;
  double *tmp = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (tmp, "trans_ext");

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  tmp [i0 * n + j0] = 0.0;
		  for (k=0; k<11; k++)
		    {
		      tmp [i0 * n + j0] +=
			r [i0 * n + j11 + k]
			* tinv [k * 11 + jj];
		    }
		}
	    }
	}
    }

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  r [i0 * n + j0] = 0.0;
		  for (k = 0; k < 11; ++k)
		    {
		      r [i0 * n + j0] +=
			tinv [ii * 11 + k]
			* tmp [(i11 + k) * n + j0];
		    }
		}
	    }
	}
    }

  free (tmp);
}
/* solve natural resistance problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts_matrix_0 (struct stokes * sys,
			 const double *u, const double *o, const double *e,
			 double *f, double *t, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_res_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int n11 = np * 11;
  double *mat = (double *)malloc (sizeof (double) * n11 * n11);
  double *b = (double *)malloc (sizeof (double) * n11);
  double *x = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_res_3fts_matrix");
  CHECK_MALLOC (b,   "solve_res_3fts_matrix");
  CHECK_MALLOC (x,   "solve_res_3fts_matrix");

  /* b := (UOE) */
  set_fts_by_FTS (np, b, u, o, e);

  // mobility matrix in EXTRACTED form
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  /* for test
  test_symmetric (n11, mat, 1.0e-12);
  */

  /* resistance matrix in INVERSED form */
  lapack_inv_ (n11, mat);
  trans_ext (np, mat); // resistance matrix in EXTRACTED form

  multiply_extmat_with_extvec_3fts (np, mat, b, x);

  set_FTS_by_fts (np, f, t, s, x);

  free (mat);
  free (b);
  free (x);
}
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts_matrix (struct stokes * sys,
		       const double *u, const double *o, const double *e,
		       double *f, double *t, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_res_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *e0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0,  "solve_res_3fts_matrix");
  CHECK_MALLOC (o0,  "solve_res_3fts_matrix");
  CHECK_MALLOC (e0,  "solve_res_3fts_matrix");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_3fts_matrix_0 (sys,
			   u0, o0, e0,
			   f, t, s);

  free (u0);
  free (o0);
  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}



/* solve natural resistance problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_lub_3fts_matrix_0 (struct stokes * sys,
			     const double *u, const double *o,
			     const double *e,
			     double *f, double *t, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_res_lub_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int n11 = np * 11;

  double *mat = (double *)malloc(sizeof (double) * n11 * n11);
  double *b = (double *)malloc(sizeof (double) * n11);
  double *x = (double *)malloc(sizeof (double) * n11);
  double *y = (double *)malloc(sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_res_lub_3fts_matrix");
  CHECK_MALLOC (b, "solve_res_lub_3fts_matrix");
  CHECK_MALLOC (x, "solve_res_lub_3fts_matrix");
  CHECK_MALLOC (y, "solve_res_lub_3fts_matrix");

  /* b := (UOE) */
  set_fts_by_FTS (np, b, u, o, e);

  make_matrix_lub_3fts (sys, mat); // lub matrix in EXTRACTED form
  multiply_extmat_with_extvec_3fts (np, mat, b, x); // x := L.(UOE)

  // mobility matrix in EXTRACTED form
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  multiply_extmat_with_extvec_3fts (np, mat, x, y); // y := M.L.(UOE)

  /* y := (I + M.L).(UOE) */
  int i;
  for (i = 0; i < n11; ++i)
    {
      b [i] += y [i];
    }

  /* resistance matrix in INVERSED form */
  lapack_inv_ (n11, mat);
  trans_ext (np, mat); // resistance matrix in EXTRACTED form

  /* x := (M^-1).(I + M.L).(UOE) */
  multiply_extmat_with_extvec_3fts (np, mat, b, x);

  set_FTS_by_fts (np, f, t, s, x);

  free (mat);
  free (b);
  free (x);
  free (y);
}
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_lub_3fts_matrix (struct stokes * sys,
			   const double *u, const double *o,
			   const double *e,
			   double *f, double *t, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_res_lub_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *e0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0, "solve_res_lub_3fts_matrix");
  CHECK_MALLOC (o0, "solve_res_lub_3fts_matrix");
  CHECK_MALLOC (e0, "solve_res_lub_3fts_matrix");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_lub_3fts_matrix_0 (sys,
			       u0, o0, e0,
			       f, t, s);

  free (u0);
  free (o0);
  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}



/** natural mobility problem **/
/* multiply transformation matrix from right-hand-side,
 * so that this could be multiplied by extracted vector
 * and return the extracted vector
 * INPUT
 *  mat [np * 11 * np * 11] :
 * OUTPUT
 *  mat [np * 11 * np * 11] := mat . t, where t.E = E~
 */
static void
trans_mat_ext2ext (int np, double * mat)
{
  static double
    t [121] = {
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  2.0,  0.0, 0.0, 0.0,  1.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  2.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 2.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 2.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  1.0,  0.0, 0.0, 0.0,  2.0};

  int i, ii, i11, i0;
  int j, jj, j11, j0;
  int k;

  int n = np * 11;
  double *tmp = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (tmp, "trans_mat_ext2ext");

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  tmp [i0 * n + j0] = 0.0;
		  for (k=0; k<11; k++)
		    {
		      tmp [i0 * n + j0] +=
			mat [i0 * n + j11 + k]
			* t [k * 11 + jj];
		    }
		}
	    }
	}
    }

  for (i = 0; i < n * n; ++i)
    {
      mat [i] = tmp [i];
    }

  free (tmp);
}
void
split_matrix_3fts (int np, const double *mat,
		   double * mat_ll, double * mat_lh,
		   double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11, n6, n5;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_ll [(i6 + ii) * n6 + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_lh [(i6 + ii) * n5 + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(i5 + ii) * n6 + j6 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(i5 + ii) * n5 + j5 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }
}
static void
merge_matrix_3fts (int np,
		   const double * mat_ll, const double * mat_lh,
		   const double * mat_hl, const double * mat_hh,
		   double *mat)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11, n6, n5;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_ll [(i6 + ii) * n6 + j6 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_lh [(i6 + ii) * n5 + j5 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hl [(i5 + ii) * n6 + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(i5 + ii) * n5 + j5 + jj];
		}
	    }
	}
    }
}
/* solve natural mobility problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_3fts_matrix_0 (struct stokes * sys,
			 const double *f, const double *t, const double *e,
			 double *u, double *o, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mob_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int n11 = np * 11;
  int n6 = np * 6;
  int n5 = np * 5;


  double *mat = (double *) malloc (sizeof (double) * n11 * n11);
  double *mat_ll = (double *) malloc (sizeof (double) * n6 * n6);
  double *mat_lh = (double *) malloc (sizeof (double) * n6 * n5);
  double *mat_hl = (double *) malloc (sizeof (double) * n5 * n6);
  double *mat_hh = (double *) malloc (sizeof (double) * n5 * n5);
  double *mob_ll = (double *) malloc (sizeof (double) * n6 * n6);
  double *mob_lh = (double *) malloc (sizeof (double) * n6 * n5);
  double *mob_hl = (double *) malloc (sizeof (double) * n5 * n6);
  double *mob_hh = (double *) malloc (sizeof (double) * n5 * n5);
  double *b = (double *) malloc (sizeof (double) * n11);
  double *x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mat_ll, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mat_lh, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mat_hl, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mat_hh, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mob_ll, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mob_lh, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mob_hl, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (mob_hh, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (b, "solve_mob_3fts_matrix_0");
  CHECK_MALLOC (x, "solve_mob_3fts_matrix_0");

  /* b := (FTE) */
  set_fts_by_FTS (np, b, f, t, e);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_3fts (np, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (n5, n6,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  merge_matrix_3fts (np, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (np, u, o, s, x);

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);
}
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  s [np * 5] :
 */
void
solve_mob_3fts_matrix (struct stokes * sys,
		       const double *f, const double *t, const double *e,
		       double *u, double *o, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mob_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  double *e0 = (double *) malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (e0, "solve_mob_3fts_matrix");

  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_mob_3fts_matrix_0 (sys,
			   f, t, e0,
			   u, o, s);

  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}


/* return A_ij = A_ik . B_kj
 * INPUT
 *  a [n * n] :
 *  b [n * n] :
 * OUTPUT
 *  a [n * n] :
 */
static void
multiply_matrices (int n, double *a, const double *b)
{
  double *tmp = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (tmp, "multiply_matrices");

  int i, j, k;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
	{
	  tmp [i * n + j] = 0.0;
	  for (k = 0; k < n; ++k)
	    {
	      tmp [i * n + j] += a [i * n + k] * b [k * n + j];
	    }
	}
    }

  for (i = 0; i < n * n; ++i)
    a [i] = tmp [i];

  free (tmp);
}
/* solve natural mobility problem in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_lub_3fts_matrix_0 (struct stokes * sys,
			     const double *f, const double *t,
			     const double *e,
			     double *u, double *o, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int n11 = np * 11;
  int n6 = np * 6;
  int n5 = np * 5;

  double *mat = (double *) malloc (sizeof (double) * n11 * n11);
  double *lub = (double *) malloc (sizeof (double) * n11 * n11);
  double *mat_ll = (double *) malloc (sizeof (double) * n6 * n6);
  double *mat_lh = (double *) malloc (sizeof (double) * n6 * n5);
  double *mat_hl = (double *) malloc (sizeof (double) * n5 * n6);
  double *mat_hh = (double *) malloc (sizeof (double) * n5 * n5);
  double *mob_ll = (double *) malloc (sizeof (double) * n6 * n6);
  double *mob_lh = (double *) malloc (sizeof (double) * n6 * n5);
  double *mob_hl = (double *) malloc (sizeof (double) * n5 * n6);
  double *mob_hh = (double *) malloc (sizeof (double) * n5 * n5);
  double *b = (double *) malloc (sizeof (double) * n11);
  double *x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (lub, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mat_ll, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mat_lh, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mat_hl, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mat_hh, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mob_ll, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mob_lh, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mob_hl, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (mob_hh, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (b, "solve_mob_lub_3fts_matrix_0");
  CHECK_MALLOC (x, "solve_mob_lub_3fts_matrix_0");

  /* b := (FTE) */
  set_fts_by_FTS (np, b, f, t, e); // (f, t, e) => b

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_3fts (np, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub matrix in EXTRACTED form */
  make_matrix_lub_3fts (sys, lub);
  /* lub := L.T, where T.(UOE) = (UOE~) */
  trans_mat_ext2ext (np, lub);
  /* mat := (M.T).(L.T) */
  multiply_matrices (n11, mat, lub);
  /* mat := I + (M.T).(L.T) */
  int i;
  for (i = 0; i < n11; ++i)
    {
      mat [i * n11 + i] += 1.0;
    }
  /* note: at this point, lub[] is free to use.
   * so that I_?? use lub [] */
  double *I_ll = lub;
  double *I_lh = I_ll + n6 * n6;
  double *I_hl = I_lh + n6 * n5;
  double *I_hh = I_hl + n5 * n6;
  split_matrix_3fts (np, mat, I_ll, I_lh, I_hl, I_hh);

  solve_gen_linear (n6, n5,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);
  /* note: at this point, I_??, therefore lub[], is free to use. */
  free (lub);

  merge_matrix_3fts (np, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (np, u, o, s, x); // x => (u, o, s)

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);
}
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
solve_mob_lub_3fts_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   const double *e,
			   double *u, double *o, double *s)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;

  double *e0 = (double *) malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (e0, "solve_mob_lub_3fts_matrix");

  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_mob_lub_3fts_matrix_0 (sys, f, t, e0,
			       u, o, s);

  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}

/** natural mobility problem with fixed particles **/
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat [np * 11 * np * 11] : matrix to split
 * OUTPUT
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 6
 */
static void
split_matrix_fix_3fts (int np, int nm,
		       const double * mat,
		       double * mat_ll, double * mat_lh,
		       double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11;
  int nl, nh;
  int nm5;


  n11 = np * 11;
  nl = nm * 6;
  nh = n11 - nl;
  nm5 = nm * 5;

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_ll [(i6 + ii) * nl + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_lh [(i6 + ii) * nh + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(i5 + ii) * nl + j6 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(i5 + ii) * nh + j5 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* lh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_lh [(i6 + ii) * nh + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_hh [(i5 + ii) * nh + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(nm5 + (i - nm) * 11 + ii) * nl + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(nm5 + (i - nm) * 11 + ii) * nh + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_hh [(nm5 + (i - nm) * 11 + ii) * nh
			 + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	    }
	}
    }
}
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 6
 * OUTPUT
 *  mat [np * 11 * np * 11] : matrix to split
 */
static void
merge_matrix_fix_3fts (int np, int nm,
		       const double * mat_ll, const double * mat_lh,
		       const double * mat_hl, const double * mat_hh,
		       double * mat)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11;
  int nl, nh;
  int nm5;


  n11 = np * 11;
  nl = nm * 6;
  nh = n11 - nl;
  nm5 = nm * 5;

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_ll [(i6 + ii) * nl + j6 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_lh [(i6 + ii) * nh + j5 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hl [(i5 + ii) * nl + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(i5 + ii) * nh + j5 + jj];
		}
	    }
	}
    }

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* lh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_lh [(i6 + ii) * nh + nm5 + (j - nm) * 11 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hh [(i5 + ii) * nh + nm5 + (j - nm) * 11 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_hl [(nm5 + (i - nm) * 11 + ii) * nl + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(nm5 + (i - nm) * 11 + ii) * nh + j5 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_hh [(nm5 + (i - nm) * 11 + ii) * nh
			     + nm5 + (j - nm) * 11 + jj];
		}
	    }
	}
    }
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : in the fluid-rest frame
 *  of [nf * 3] : in the fluid-rest frame
 *  ef [nf * 5] : in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts_matrix_0 (struct stokes * sys,
			 const double *f, const double *t,
			 const double *e,
			 const double *uf, const double *of,
			 const double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mix_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3fts_matrix_0 (sys,
			       f, t, e,
			       u, o, s);
      return;
    }

  int n11 = np * 11;
  int nf = np - nm;
  int nm11 = nm * 11;
  int nl = nm * 6;
  int nh = n11 - nl;


  double *mat = (double *) malloc (sizeof (double) * n11 * n11);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n11);
  double *x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_3fts_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_3fts_matrix");
  CHECK_MALLOC (b, "solve_mix_3fts_matrix");
  CHECK_MALLOC (x, "solve_mix_3fts_matrix");

  /* b := (FTE) */
  set_fts_by_FTS (nm, b, f, t, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_fix_3fts (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (nh, nl,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  merge_matrix_fix_3fts (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_3fts_matrix (struct stokes * sys,
		       const double *f, const double *t,
		       const double *e,
		       const double *uf, const double *of,
		       const double *ef,
		       double *u, double *o, double *s,
		       double *ff, double *tf, double *sf)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mix_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3fts_matrix (sys,
			     f, t, e,
			     u, o, s);
      return;
    }

  int nf = np - nm;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *ef0 = (double *) malloc (sizeof (double) * nf * 5);
  double *e0  = (double *) malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (uf0, "solve_mix_3fts_matrix");
  CHECK_MALLOC (of0, "solve_mix_3fts_matrix");
  CHECK_MALLOC (ef0, "solve_mix_3fts_matrix");
  CHECK_MALLOC (e0, "solve_mix_3fts_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  shift_labo_to_rest_E (sys, nf, ef, ef0);
  shift_labo_to_rest_E (sys, nm, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_mix_3fts_matrix_0 (sys, f, t, e0,
			   uf0, of0, ef0,
			   u, o, s, ff, tf, sf);

  free (uf0);
  free (of0);
  free (ef0);
  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : in the fluid-rest frame
 *  of [nf * 3] : in the fluid-rest frame
 *  ef [nf * 5] : in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_lub_3fts_matrix_0 (struct stokes * sys,
			     const double *f, const double *t,
			     const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3fts_matrix_0 :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3fts_matrix_0 (sys,
				   f, t, e,
				   u, o, s);
      return;
    }

  int n11 = np * 11;
  int nf = np - nm;
  int nm11 = nm * 11;
  int nl = nm * 6;
  int nh = n11 - nl;

  double *mat = (double *) malloc (sizeof (double) * n11 * n11);
  double *lub = (double *) malloc (sizeof (double) * n11 * n11);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n11);
  double *x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (mat, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (lub, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (b, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (x, "solve_mix_lub_3fts_matrix");

  /* b := (FTE) */
  set_fts_by_FTS (nm, b, f, t, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_3all (sys, mat); // sys->version is 2 (FTS)
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_fix_3fts (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub matrix in EXTRACTED form */
  make_matrix_lub_3fts (sys, lub);
  /* lub := L.T, where T.(UOE) = (UOE~) */
  trans_mat_ext2ext (np, lub);
  /* mat := (M.T).(L.T) */
  multiply_matrices (n11, mat, lub);
  /* mat := I + (M.T).(L.T) */
  int i;
  for (i = 0; i < n11; ++i)
    {
      mat [i * n11 + i] += 1.0;
    }
  /* note: at this point, lub[] is free to use.
   * so that I_?? use lub [] */
  double *I_ll = lub;
  double *I_lh = I_ll + nl * nl;
  double *I_hl = I_lh + nl * nh;
  double *I_hh = I_hl + nh * nl;
  split_matrix_fix_3fts (np, nm, mat, I_ll, I_lh, I_hl, I_hh);

  solve_gen_linear (nl, nh,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);
  /* note: at this point, I_??, therefore lub[], is free to use. */
  free (lub);

  merge_matrix_fix_3fts (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_lub_3fts_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   const double *e,
			   const double *uf, const double *of,
			   const double *ef,
			   double *u, double *o, double *s,
			   double *ff, double *tf, double *sf)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3fts_matrix :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3fts_matrix (sys,
				 f, t, e,
				 u, o, s);
      return;
    }

  int nf = np - nm;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *ef0 = (double *) malloc (sizeof (double) * nf * 5);
  double *e0  = (double *) malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (uf0, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (of0, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (ef0, "solve_mix_lub_3fts_matrix");
  CHECK_MALLOC (e0, "solve_mix_lub_3fts_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  shift_labo_to_rest_E (sys, nf, ef, ef0);
  shift_labo_to_rest_E (sys, nm, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_mix_lub_3fts_matrix_0 (sys, f, t, e0,
			       uf0, of0, ef0,
			       u, o, s, ff, tf, sf);

  free (uf0);
  free (of0);
  free (ef0);
  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
