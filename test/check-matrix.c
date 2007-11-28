/* test code for matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-matrix.c,v 1.2 2007/11/28 03:42:14 kichiki Exp $
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
#include "../config.h"

#include <stdio.h>
#include <stdlib.h>
#include "memory-check.h"
#include "check.h" // compare()

#include <matrix.h> // multiply_extmat_with_extvec_3fts()

#include <bench.h> // ptime_ms_d()

#ifdef HAVE_CBLAS_H
#include <cblas.h>
#endif

#include "dgemm_c.h"


void
mul_matrices_ (const double *a, int na1, int na2,
	       const double *b, int nb1, int nb2,
	       double *c)
{
  if (na2 != nb1)
    {
      fprintf (stderr, "mul_matrices_ : invalid dimension"
	       " (na2)%d != (nb1)%d\n",
	       na2, nb1);
      exit (1);
    }
  if (a == c || b == c)
    {
      fprintf (stderr, "illeagal output pointer in mul_matrices().\n");
      exit (1);
    }

  int i, j, k;
  for (i = 0; i < na1; i ++)
    {
      for (j = 0; j < nb2; j ++)
	{
	  c[i*nb2+j] = 0.0;
	  for (k = 0; k < na2; k ++)
	    {
	      c[i*nb2+j] += a[i*na2+k] * b[k*nb2+j];
	    }
	}
    }
}
void
mul_matrices_blas (const double *A, int na1, int na2,
		   const double *B, int nb1, int nb2,
		   double *C)
{
  if (na2 != nb1)
    {
      fprintf (stderr, "illeagal multiplication in mul_matrices().\n");
      exit (1);
    }
  if (A == C || B == C)
    {
      fprintf (stderr, "illeagal output pointer in mul_matrices().\n");
      exit (1);
    }

  /* use Fortran BLAS routines */
  dgemm_wrap (na1, nb2, na2,
	      1.0, A,
	      B,
	      0.0, C);
}

#ifdef HAVE_CBLAS_H
void
mul_matrices_cblas (const double *A, int na1, int na2,
		    const double *B, int nb1, int nb2,
		    double *C)
{
  if (na2 != nb1)
    {
      fprintf (stderr, "illeagal multiplication in mul_matrices().\n");
      exit (1);
    }
  if (A == C || B == C)
    {
      fprintf (stderr, "illeagal output pointer in mul_matrices().\n");
      exit (1);
    }

  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       na1, nb2, na2,
	       1.0, A, na2,
	       B, nb2,
	       0.0, C, nb2);
}
#endif // HAVE_CBLAS_H

	       
/* C = A.B + beta C
 * INPUT
 *  c[m,n]
 *  a[m,k]
 *  b[k,n]
 */
static void
my_dgemm (int m, int n, int k,
	  double alpha, const double *a, int lda,
	  const double *b, int ldb,
	  double beta, double *c, int ldc)
{
  int i, j, l;
  for (i = 0; i < m; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  c[i*ldc + j] = beta * c[i*ldc + j];
	  for (l = 0; l < k; l ++)
	    {
	      c[i*ldc + j] += a[i*lda + l] * b[l*ldb + j];
	    }
	}
    }
}
int
benchmark_dgemm (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "benchmark_dgemm : start\n");
#ifndef HAVE_CBLAS_H
      fprintf (stdout, "# n, t1, t2, t3, where\n"
	       "# t1 : by mul_matrices_ -- double for loop\n"
	       "# t2 : by my_dgemm -- another double for loop\n"
	       "# t3 : by dgemm in Fortran BLAS\n");
#else // HAVE_CBLAS_H
      fprintf (stdout, "# n, t1, t2, t3, t4, where\n"
	       "# t1 : by mul_matrices_ -- double for loop\n"
	       "# t2 : by my_dgemm -- another double for loop\n"
	       "# t3 : by dgemm in Fortran BLAS\n"
	       "# t4 : by cblas_dgemm in ALTLAS\n");
#endif // HAVE_CBLAS_H
    }

  int check = 0;


  /**
   * initialization
   */
  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  double *c1 = (double *)malloc (sizeof (double) * n * n);
  double *c2 = (double *)malloc (sizeof (double) * n * n);
  double *c3 = (double *)malloc (sizeof (double) * n * n);
  double *c4 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "benchmark_mul_matrices");
  CHECK_MALLOC (b,  "benchmark_mul_matrices");
  CHECK_MALLOC (c1, "benchmark_mul_matrices");
  CHECK_MALLOC (c2, "benchmark_mul_matrices");
  CHECK_MALLOC (c3, "benchmark_mul_matrices");
  CHECK_MALLOC (c4, "benchmark_mul_matrices");

  // set the matrix 'a'
  srand48(0);
  int i;
  for (i = 0; i < n * n; i ++)
    {
      a[i] = drand48();
      b[i] = drand48();
    }

  /**
   * mul_matrices_
   */
  double t0;
  // warming up
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_ (a, n, n,
		     b, n, n,
		     c1);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_ (a, n, n,
		     b, n, n,
		     c1);
    }
  double t1 = ptime_ms_d ();
  t1 -= t0;
  t1 *= 0.1;


  /**
   * my_dgemm
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      my_dgemm (n, n, n,
		1.0, a, n,
		b, n,
		0.0, c2, n);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      my_dgemm (n, n, n,
		1.0, a, n,
		b, n,
		0.0, c2, n);
    }
  double t2 = ptime_ms_d ();
  t2 -= t0;
  t2 *= 0.1;


  /**
   * dgemm_wrap (Fortran BLAS)
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      dgemm_wrap (n, n, n,
		  1.0, a, b,
		  0.0, c3);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      dgemm_wrap (n, n, n,
		  1.0, a, b,
		  0.0, c3);
    }
  double t3 = ptime_ms_d ();
  t3 -= t0;
  t3 *= 0.1;


#ifndef HAVE_CBLAS_H
  fprintf (stdout, "%d %f %f %f\n", n, t1, t2, t3);
#else // HAVE_CBLAS_H
  /**
   * cblas_dgemm (ATLAS dgemm)
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
		   n, n, n,
		   1.0, a, n,
		   b, n,
		   0.0, c4, n);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
		   n, n, n,
		   1.0, a, n,
		   b, n,
		   0.0, c4, n);
    }
  double t4 = ptime_ms_d ();
  t4 -= t0;
  t4 *= 0.1;

  fprintf (stdout, "%d %f %f %f %f\n", n, t1, t2, t3, t4);
#endif // HAVE_CBLAS_H

  char label[80];
  int j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# c1,c2 [%d %d]", i, j);
	  check += compare (c1[i*n+j], c2[i*n+j], label, verbose, tiny);
	}
    }
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# c1,c3 [%d %d]", i, j);
	  check += compare (c1[i*n+j], c3[i*n+j], label, verbose, tiny);
	}
    }
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# c1,c3 [%d %d]", i, j);
	  check += compare (c1[i*n+j], c3[i*n+j], label, verbose, tiny);
	}
    }
#ifdef HAVE_CBLAS_H
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# c1,c4 [%d %d]", i, j);
	  check += compare (c1[i*n+j], c4[i*n+j], label, verbose, tiny);
	}
    }
#endif // HAVE_CBLAS_H

  free (a);
  free (b);
  free (c1);
  free (c2);
  free (c3);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}




int
benchmark_mul_matrices (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "benchmark_mul_matrices : start\n");
#ifndef HAVE_CBLAS_H
      fprintf (stdout, "# n, t1, t2, t3, where\n"
	       "# t1 : by mul_matrices in libstokes\n"
	       "# t2 : by mul_matrices_ -- double for loop\n"
	       "# t3 : by mul_matrices_blas by dgemm in Fortran BLAS\n");
#else // HAVE_CBLAS_H
      fprintf (stdout, "# n, t1, t2, t3, t4, where\n"
	       "# t1 : by mul_matrices in libstokes\n"
	       "# t2 : by mul_matrices_ -- double for loop\n"
	       "# t3 : by mul_matrices_blas by dgemm in Fortran BLAS\n"
	       "# t4 : by mul_matrices_cblas by cblas_dgemm in ALTLAS\n");
#endif // HAVE_CBLAS_H
    }

  int check = 0;


  /**
   * initialization
   */
  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  double *c1 = (double *)malloc (sizeof (double) * n * n);
  double *c2 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "benchmark_mul_matrices");
  CHECK_MALLOC (b,  "benchmark_mul_matrices");
  CHECK_MALLOC (c1, "benchmark_mul_matrices");
  CHECK_MALLOC (c2, "benchmark_mul_matrices");

  // set the matrix 'a'
  srand48(0);
  int i;
  for (i = 0; i < n * n; i ++)
    {
      a[i] = drand48();
      b[i] = drand48();
    }

  double t0;
  /**
   * by mul_matrices (with some BLAS routine)
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      mul_matrices (a, n, n, b, n, n, c1);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      mul_matrices (a, n, n, b, n, n, c1);
    }
  double t1 = ptime_ms_d ();
  t1 -= t0;
  t1 *= 0.1;

  /**
   * by direct for-loops
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_ (a, n, n, b, n, n, c2);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_ (a, n, n, b, n, n, c2);
    }
  double t2 = ptime_ms_d ();
  t2 -= t0;
  t2 *= 0.1;

  /**
   * by blas version
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_blas (a, n, n, b, n, n, c2);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_blas (a, n, n, b, n, n, c2);
    }
  double t3 = ptime_ms_d ();
  t3 -= t0;
  t3 *= 0.1;

#ifndef HAVE_CBLAS_H
  fprintf (stdout, "%d %f %f %f\n", n, t1, t2, t3);
#else // HAVE_CBLAS_H

  /**
   * by cblas version
   */
  // warming up
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_cblas (a, n, n, b, n, n, c2);
    }
  // take an average for 10 calculations
  t0 = ptime_ms_d ();
  for (i = 0; i < 10; i ++)
    {
      mul_matrices_cblas (a, n, n, b, n, n, c2);
    }
  double t4 = ptime_ms_d ();
  t4 -= t0;
  t4 *= 0.1;

  fprintf (stdout, "%d %f %f %f %f\n", n, t1, t2, t3, t4);
#endif // HAVE_CBLAS_H

  char label[80];
  int j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "# c[%d %d]", i, j);
	  check += compare (c1[i*n+j], c2[i*n+j], label, verbose, tiny);
	}
    }

  free (a);
  free (b);
  free (c1);
  free (c2);


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
