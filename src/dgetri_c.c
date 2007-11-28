/* C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgetri_c.c,v 1.7 2007/11/28 03:32:43 kichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h>
#include "memory-check.h" // CHECK_MALLOC


void dcopy_ (int *, double *, int *, double *, int *);

/*
DGETRI(l)                              )                             DGETRI(l)



NAME
       DGETRI  -  compute  the  inverse of a matrix using the LU factorization
       computed by DGETRF

SYNOPSIS
       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

           INTEGER        INFO, LDA, LWORK, N

           INTEGER        IPIV( * )

           DOUBLE         PRECISION A( LDA, * ), WORK( * )

PURPOSE
       DGETRI computes the inverse of a matrix using the LU factorization com-
       puted  by  DGETRF.   This  method inverts U and then computes inv(A) by
       solving the system inv(A)*L = inv(U) for inv(A).


ARGUMENTS
       N       (input) INTEGER
               The order of the matrix A.  N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On entry, the factors L and U from the factorization A =  P*L*U
               as  computed  by  DGETRF.  On exit, if INFO = 0, the inverse of
               the original matrix A.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,N).

       IPIV    (input) INTEGER array, dimension (N)
               The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix
               was interchanged with row IPIV(i).

       WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
               On exit, if INFO=0, then WORK(1) returns the optimal LWORK.

       LWORK   (input) INTEGER
               The dimension of the array WORK.  LWORK >= max(1,N).  For opti-
               mal performance LWORK >= N*NB, where NB is the  optimal  block-
               size returned by ILAENV.

               If  LWORK  = -1, then a workspace query is assumed; the routine
               only calculates the optimal size of  the  WORK  array,  returns
               this  value  as the first entry of the WORK array, and no error
               message related to LWORK is issued by XERBLA.

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is singu-
               lar and its inverse could not be computed.



LAPACK version 3.0               15 June 2000                        DGETRI(l)
*/

void dgetri_ (int *n, double *a, int *lda, int *ipiv,
	      double *work, int *lwork,
	      int *info);


/*
DGETRF(l)                              )                             DGETRF(l)



NAME
       DGETRF - compute an LU factorization of a general M-by-N matrix A using
       partial pivoting with row interchanges

SYNOPSIS
       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

           INTEGER        INFO, LDA, M, N

           INTEGER        IPIV( * )

           DOUBLE         PRECISION A( LDA, * )

PURPOSE
       DGETRF computes an LU factorization of a general M-by-N matrix A  using
       partial pivoting with row interchanges.  The factorization has the form
          A = P * L * U
       where P is a permutation matrix, L is lower triangular with unit diago-
       nal  elements  (lower  trapezoidal if m > n), and U is upper triangular
       (upper trapezoidal if m < n).

       This is the right-looking Level 3 BLAS version of the algorithm.


ARGUMENTS
       M       (input) INTEGER
               The number of rows of the matrix A.  M >= 0.

       N       (input) INTEGER
               The number of columns of the matrix A.  N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On entry, the M-by-N matrix to be factored.  On exit, the  fac-
               tors  L and U from the factorization A = P*L*U; the unit diago-
               nal elements of L are not stored.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,M).

       IPIV    (output) INTEGER array, dimension (min(M,N))
               The pivot indices; for 1 <= i <= min(M,N), row i of the  matrix
               was interchanged with row IPIV(i).

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               >  0:   if  INFO = i, U(i,i) is exactly zero. The factorization
               has been completed, but the factor U is exactly  singular,  and
               division  by zero will occur if it is used to solve a system of
               equations.



LAPACK version 3.0               15 June 2000                        DGETRF(l)
*/
void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv,
	      int *info);



void lapack_inv (int n, const double *a,
		 double *ai)
{
  int *ipvt = (int *)malloc (sizeof (int) * n);
  double *work = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ipvt, "lapack_inv");
  CHECK_MALLOC (work, "lapack_inv");

  int i_1 = 1;
  int nn = n * n;
  dcopy_ (&nn, a, &i_1, ai, &i_1);

  int info;
  dgetrf_ (&n, &n, ai, &n, ipvt, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetrf (info = %d)\n", info);
      exit (1);
    }

  dgetri_ (&n, ai, &n, ipvt, work, &n, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetri (info = %d)\n", info);
      exit (1);
    }

  free (ipvt);
  free (work);
}

/* the version that a[n*n] is input AND output
 */
void lapack_inv_ (int n, double *a)
{
  int info;
  int *ipvt = (int *)malloc (sizeof (int) * n);
  double *work = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ipvt, "lapack_inv_");
  CHECK_MALLOC (work, "lapack_inv_");

  dgetrf_ (&n, &n, a, &n, ipvt, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetrf (info = %d)\n", info);
      exit (1);
    }

  dgetri_ (&n, a, &n, ipvt, work, &n, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetri (info = %d)\n", info);
      exit (1);
    }

  free (ipvt);
  free (work);
}
