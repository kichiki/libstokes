/* matrix-manipulating routines
 * Copyright (C) 2001-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: matrix.c,v 1.9 2007/11/04 01:12:37 kichiki Exp $
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
#include "../config.h" // for HAVE_*_H

#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "dgetri_c.h" /* lapack_inv_() */

#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */

#include <cblas.h>

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
/* use Fortran BLAS routines */

#include "dgemm_c.h"
#include "dgemv_c.h"

# else // !HAVE_BLAS_H
/* use local BLAS routines */

/* y = alpha * A.x + beta y
 * INPUT
 *  a[m,n] : (i,j) element is allocated at a[i*lda+j]
 *  y[m]   : (i)   element is allocated at y[i*incy]
 *  x[n]   : (j)   element is allocated at x[i*incx]
 */
static void
my_dgemv (int m, int n,
	  double alpha, const double *a, int lda,
	  const double *x, int incx, double beta,
	  double *y, int incy)
{
  int i, j;

  for (i = 0; i < m; i ++)
    {
      y[i] = beta * y[i];
      for (j = 0; j < n; j ++)
	{
	  y[i] += a[i*lda+j] * x[j];
	}
    }
}

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


# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


#include "memory-check.h" // macro CHECK_MALLOC
#include "matrix.h"


/* return B . A in A, where B is square matrix
 * Remark!! this is not A.B!!
 * INPUT
 *  A [na1, na2]
 *  B [nb, nb] : square matrix!
 *  where (nb == na1)
 * OUTPUT
 *  A [na1, na2] = B [nb, nb] . A [na1, na2]
 */
static void
mul_left_sq (double * A, int na1, int na2,
	     const double * B)
{
  int i;
  int nb = na1;

  double *tmp = (double *) malloc (sizeof (double) * na1 * na2);
  CHECK_MALLOC (tmp, "mul_left_sq");

  mul_matrices (B, nb, nb,
		A, na1, na2,
		tmp);

  for (i = 0; i < na1 * na2; ++i)
    {
      A [i] = tmp [i];
    }

  free (tmp);
}

/*
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  C [nc1, nc2]
 *  a
 *  b
 *  where (nb2 == nc1) && (na1 = nb1) && (na2 == nc2)
 * OUTPUT
 *  D [na1, na2] = a * A [na1, na2] + b * B [nb1, nb2] . C [nc1, nc2]
 */
static void
add_and_mul (const double *A, int na1, int na2,
	     const double *B, int nb1, int nb2,
	     const double *C, int nc1, int nc2,
	     double a, double b,
	     double *D)
{
  if (nb2 != nc1
      || na1 != nb1
      || na2 != nc2)
    {
      fprintf (stderr, "illeagal multiplication in add_and_mul().\n");
      exit (1);
    }
  if (A == D)
    {
      fprintf (stderr, "illeagal pointer D in add_and_mul().\n");
      exit (1);
    }

  mul_matrices (B, nb1, nb2, C, nc1, nc2, D);

  int i;
  for (i = 0; i < na1*na2; i ++)
    {
      D[i] = a * A[i] + b * D[i];
    }
}

/* solve generalized linear set of equations using LU-decomposition
 * INPUT
 *  n1, n2 : dimension
 *  A [n1 * n1] :
 *  B [n1 * n2] :
 *  C [n2 * n1] :
 *  D [n2 * n2] :
 *
 *  E [n1 * n1] :
 *  F [n1 * n2] :
 *  G [n2 * n1] :
 *  H [n2 * n2] :
 *  where the generalized linear set of equations is
 *   [A B](x) = [E F](b)
 *   [C D](y)   [G H](c)
 * OUTPUT
 *  I [n1 * n1] :
 *  J [n1 * n2] :
 *  K [n2 * n1] :
 *  L [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (x) = [I J](b)
 *   (c)   [K L](y)
 *  note that A-D, E-H are destroyed!
 */
void
solve_gen_linear (int n1, int n2,
		  double * A, double * B, double * C, double * D,
		  double * E, double * F, double * G, double * H,
		  double * I, double * J, double * K, double * L)
{
  /* H := H^-1 */
  lapack_inv_ (n2, H);

  /* C := (H^-1) . C */
  mul_left_sq (C, n2, n1, H);
  /* G := (H^-1) . G */
  mul_left_sq (G, n2, n1, H);
  /* D := (H^-1) . D */
  mul_left_sq (D, n2, n2, H);


  double *a = (double *)malloc (sizeof (double) * n1 * n1);
  double *e = (double *)malloc (sizeof (double) * n1 * n1);
  double *b = (double *)malloc (sizeof (double) * n1 * n2);
  CHECK_MALLOC (a, "solve_gen_linear");
  CHECK_MALLOC (e, "solve_gen_linear");
  CHECK_MALLOC (b, "solve_gen_linear");

  /* a [n1, n1] := A - F.(H^-1).C */
  add_and_mul (A, n1, n1, F, n1, n2, C, n2, n1, 1.0, -1.0, a);
  // e [n1, n1] := E - F.(H^-1).G
  add_and_mul (E, n1, n1, F, n1, n2, G, n2, n1, 1.0, -1.0, e);
  // b [n1, n2] := - B + F.(H^-1).D
  add_and_mul (B, n1, n2, F, n1, n2, D, n2, n2, -1.0, 1.0, b);

  /* a := (A-F.(H^-1).C)^-1 */
  lapack_inv_ (n1, a);

  /* I := a.e */
  mul_matrices (a, n1, n1, e, n1, n1, I);
  /* J := a.b */
  mul_matrices (a, n1, n1, b, n1, n2, J);

  free (a);
  free (e);
  free (b);


  // K := - G + C.I
  add_and_mul (G, n2, n1, C, n2, n1, I, n1, n1, -1.0, 1.0, K);
  // L := D + C.J
  add_and_mul (D, n2, n2, C, n2, n1, J, n1, n2, 1.0, 1.0, L);
}

/* solve linear set of equations using LU-decomposition
 * INPUT
 *  n1, n2 : dimension
 *  A [n1 * n1] :
 *  B [n1 * n2] :
 *  C [n2 * n1] :
 *  D [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (b) = [A B](x)
 *   (c)   [C D](y)
 * OUTPUT
 *  I [n1 * n1] :
 *  J [n1 * n2] :
 *  K [n2 * n1] :
 *  L [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (x) = [I J](b)
 *   (c)   [K L](y)
 *  note that A-D are destroyed!
 */
void
solve_linear (int n1, int n2,
	      double * A, double * B, double * C, double * D,
	      double * I, double * J, double * K, double * L)
{
  int i;
  /* type 1 */

  /* A := A^-1 */
  lapack_inv_ (n1, A);

  /* B := (A^-1).B */
  double *tmp = (double *) malloc (sizeof (double) * n1 * n2);
  CHECK_MALLOC (tmp, "solve_linear");

  mul_matrices (A, n1, n1,
		B, n1, n2,
		tmp);

  for (i = 0; i < n1 * n2; i ++)
    {
      B[i] = tmp[i];
    }
  free (tmp);

  // L [n2, n2] := D - C.(A^-1).B
  add_and_mul (D, n2, n2, C, n2, n1, B, n1, n2, 1.0, -1.0, L);
  // K := C.A
  mul_matrices (C, n2, n1, A, n1, n1, K);

  for (i = 0; i < n1 * n1; ++i)
    {
      I [i] = A [i];
    }
  for (i = 0; i < n1 * n2; ++i)
    {
      J [i] = - B [i];
    }
}

/* multiply two matrices (a wrapper to BLAS routine)
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  where (na2 == nb1)
 * OUTPUT
 *  C [na1, nb2] = A [na1, na2] . B [nb1, nb2]
 */
void
mul_matrices (const double *A, int na1, int na2,
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

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       na1, nb2, na2,
	       1.0, A, na2,
	       B, nb2,
	       0.0, C, nb2);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemm_wrap (na1, nb2, na2,
	      1.0, A,
	      B,
	      0.0, C);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemm (na1, nb2, na2,
	    1.0, A, na2,
	    B, nb2,
	    0.0, C, nb2);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
}


/*
 * INPUT
 *  mat [n1, n2]
 *  x [n2]
 * OUTPUT
 *  y [n1] = mat [n1, n2] . x [n2]
 */
void
dot_prod_matrix (const double * mat, int n1, int n2,
		 const double * x,
		 double * y)
{
#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemv (CblasRowMajor, CblasNoTrans,
	       n1, n2,
	       1.0, mat, n2,
	       x, 1,
	       0.0, y, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemv_wrap (n1, n2,
	      1.0, mat,
	      x,
	      0.0, y);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemv (n1, n2,
	    1.0, mat, n2,
	    x, 1,
	    0.0, y, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
}


/* utility routine for matrix in the extracted form
 * INPUT
 *  np : # particles (not # elements!)
 *  m [np *11 * np *11] : matrix in the extracted form
 *  x [np *11] : vector in the extracted form
 * INPUT
 *  y [np *11] : output vector in the extracted form (:= m.x)
 */
void
multiply_extmat_with_extvec_3fts (int np, const double * m, const double * x,
				  double * y)
{
  int n11;
  int i;
  int j;
  int j11;
  int jj;


  n11 = np * 11;

  for (i = 0; i < n11; ++i)
    {
      y [i] = 0.0;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (jj = 0; jj < 6; ++jj)
	    {
	      y [i] += m [i * n11 + j11 + jj] * x [j11 + jj];
	    }
	  y [i] += m [i * n11 + j11 + 6]
	    * (2.0 * x [j11 + 6] + x [j11 + 10]);
	  y [i] += m [i * n11 + j11 + 7] * 2.0 * x [j11 + 7];
	  y [i] += m [i * n11 + j11 + 8] * 2.0 * x [j11 + 8];
	  y [i] += m [i * n11 + j11 + 9] * 2.0 * x [j11 + 9];
	  y [i] += m [i * n11 + j11 + 10]
	    * (2.0 * x [j11 + 10] + x [j11 + 6]);
	}
    }
}
