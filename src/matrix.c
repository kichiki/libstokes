/* matrix-manipulating routines
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: matrix.c,v 1.5 2006/10/22 22:47:10 kichiki Exp $
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


#include "matrix.h"


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
  /* type 2 */

  /* H := H^-1 */
  lapack_inv_ (n2, H);

  /* C := (H^-1) . C */
  mul_left_sq (C, n2, n1, H, n2);
  /* G := (H^-1) . G */
  mul_left_sq (G, n2, n1, H, n2);
  /* D := (H^-1) . D */
  mul_left_sq (D, n2, n2, H, n2);

  /* A [n1, n1] := A - F.(H^-1).C */
  add_and_mul (A, n1, n1, F, n1, n2, C, n2, n1, 1.0, -1.0, A);
  /* E [n1, n1] := E - F.(H^-1).G */
  add_and_mul (E, n1, n1, F, n1, n2, G, n2, n1, 1.0, -1.0, E);
  /* B [n1, n2] := - B + F.(H^-1).D */
  add_and_mul (B, n1, n2, F, n1, n2, D, n2, n2, -1.0, 1.0, B);

  /* A := (A-F.(H^-1).C)^-1 */
  //lu_inv (A, n1);
  lapack_inv_ (n1, A);

  /* I := A.E */
  mul_matrices (A, n1, n1, E, n1, n1, I);
  /* J := A.B */
  mul_matrices (A, n1, n1, B, n1, n2, J);

  /* K := - G + C.I */
  add_and_mul (G, n2, n1, C, n2, n1, I, n1, n1, -1.0, 1.0, K);
  /* L := D + C.J */
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
  double *tmp;
  tmp = (double *) malloc (sizeof (double) * n1 * n2);

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       n1, n2, n1,
	       1.0, A, n1,
	       B, n2,
	       0.0, tmp, n2);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemm_wrap (n1, n2, n1,
	      1.0, A,
	      b,
	      0.0, tmp);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemm (n1, n2, n1,
	    1.0, A, n1,
	    B, n2,
	    0.0, tmp, n2);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  for (i = 0; i < n1 * n2; i ++)
    {
      B[i] = tmp[i];
    }
  free (tmp);

  /* L [n2, n2] := D - C.(A^-1).B */
  for (i = 0; i < n2 * n2; i ++)
    {
      L[i] = D[i];
    }
#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       n2, n2, n1,
	       -1.0, C, n1,
	       B, n2,
	       1.0, L, n2);

  /* K := C.A */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       n2, n1, n1,
	       1.0, C, n1,
	       A, n1,
	       0.0, K, n1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemm_wrap (n2, n2, n1,
	      -1.0, C,
	      B,
	      1.0, L);

  /* K := C.A */
  dgemm_wrap (n2, n1, n1,
	      1.0, C,
	      A,
	      0.0, K);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemm (n2, n2, n1,
	    -1.0, C, n1,
	    B, n2,
	    1.0, L, n2);

  my_dgemm (n2, n1, n1,
	    1.0, C, n1,
	    A, n1,
	    0.0, K, n1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  for (i = 0; i < n1 * n1; ++i)
    {
      I [i] = A [i];
    }
  for (i = 0; i < n1 * n2; ++i)
    {
      J [i] = - B [i];
    }
}

/*
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  where (na2 == nb1)
 * OUTPUT
 *  C [na1, nb2] = A [na1, na2] . B [nb1, nb2]
 */
void
mul_matrices (const double * A, int na1, int na2,
	      const double * B, int nb1, int nb2,
	      double * C)
{
  if (na2 != nb1)
    {
      fprintf (stderr, "illeagal multiplication in mul_left().\n");
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
	      1.0, C);

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
 *  A [na1, na2]
 *  B [nb, nb] : square matrix!
 *  where (nb == na1)
 * OUTPUT
 *  A [na1, na2] = B [nb, nb] . A [na1, na2]
 */
void
mul_left_sq (double * A, int na1, int na2,
	     const double * B, int nb)
{
  int i;
  double * tmp;


  if (nb != na1)
    {
      fprintf (stderr, "illeagal multiplication in mul_left_sq().\n");
      exit (1);
    }

  tmp = (double *) malloc (sizeof (double) * na1 * na2);
  if (tmp == NULL)
    {
      fprintf (stderr, "allocation error in mul_left_sq().\n");
      exit (1);
    }
#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       na1, na2, nb,
	       1.0, B, nb,
	       A, na2,
	       0.0, tmp, na2);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemm_wrap (na1, na2, nb,
	      1.0, B,
	      A,
	      0.0, tmp);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemm (na1, na2, nb,
	    1.0, B, nb,
	    A, na2,
	    0.0, tmp, na2);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

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
 *  D could be same to A itself!
 */
void
add_and_mul (const double * A, int na1, int na2,
	     const double * B, int nb1, int nb2,
	     const double * C, int nc1, int nc2,
	     double a, double b,
	     double * D)
{
  int i;

  if (nb2 != nc1
      || na1 != nb1
      || na2 != nc2)
    {
      fprintf (stderr, "illeagal multiplication in add_and_mul().\n");
      exit (1);
    }

  for (i = 0; i < na1 * na2; i ++)
    {
      D[i] = A[i];
    }

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
	       nb1, nc2, nb2,
	       b, B, nb2,
	       C, nc2,
	       a, D, na2);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */
  dgemm_wrap (nb1, nc2, nb2,
	      b, B,
	      C,
	      a, D);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */
  my_dgemm (nb1, nc2, nb2,
	    b, B, nb2,
	    C, nc2,
	    a, D, na2);

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
