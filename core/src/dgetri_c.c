/* C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgetri_c.c,v 1.8 2007/12/01 18:42:05 kichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h>
#include "memory-check.h" // CHECK_MALLOC


void dcopy_ (int *, double *, int *, double *, int *);

/*
DGETRI(1)                LAPACK routine (version 3.1)                DGETRI(1)



NAME
       DGETRI - the inverse of a matrix using the LU factorization computed by
       DGETRF

SYNOPSIS
       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

           INTEGER        INFO, LDA, LWORK, N

           INTEGER        IPIV( * )

           DOUBLE         PRECISION A( LDA, * ), WORK( * )

PURPOSE
       DGETRI computes the inverse of a matrix using the LU factorization com-
       puted by DGETRF.

       This  method  inverts  U and then computes inv(A) by solving the system
       inv(A)*L = inv(U) for inv(A).


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

       WORK       (workspace/output)   DOUBLE   PRECISION   array,   dimension
       (MAX(1,LWORK))
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



 LAPACK routine (version 3.1)    November 2006                       DGETRI(1)
*/

void dgetri_ (int *n, double *a, int *lda, int *ipiv,
	      double *work, int *lwork,
	      int *info);


/*
DGETRF(1)                LAPACK routine (version 3.1)                DGETRF(1)



NAME
       DGETRF - an LU factorization of a general M-by-N matrix A using partial
       pivoting with row interchanges

SYNOPSIS
       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

           INTEGER        INFO, LDA, M, N

           INTEGER        IPIV( * )

           DOUBLE         PRECISION A( LDA, * )

PURPOSE
       DGETRF computes an LU factorization of a general M-by-N matrix A  using
       partial pivoting with row interchanges.

       The factorization has the form
          A = P * L * U
       where P is a permutation matrix, L is lower triangular with unit diago-
       nal elements (lower trapezoidal if m > n), and U  is  upper  triangular
       (upper trapezoidal if m < n).

       This is the right-looking Level 3 BLAS version of the algorithm.


ARGUMENTS
       M       (input) INTEGER
               The number of rows of the matrix A.  M >= 0.

       N       (input) INTEGER
               The number of columns of the matrix A.  N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On  entry, the M-by-N matrix to be factored.  On exit, the fac-
               tors L and U from the factorization A = P*L*U; the unit  diago-
               nal elements of L are not stored.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,M).

       IPIV    (output) INTEGER array, dimension (min(M,N))
               The  pivot indices; for 1 <= i <= min(M,N), row i of the matrix
               was interchanged with row IPIV(i).

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, U(i,i) is exactly  zero.  The  factorization
               has  been  completed, but the factor U is exactly singular, and
               division by zero will occur if it is used to solve a system  of
               equations.



 LAPACK routine (version 3.1)    November 2006                       DGETRF(1)
*/
void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv,
	      int *info);


/*
DGETRS(1)                LAPACK routine (version 3.1)                DGETRS(1)



NAME
       DGETRS  -  a system of linear equations  A * X = B or A' * X = B with a
       general N-by-N matrix A using the LU factorization computed by DGETRF

SYNOPSIS
       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

           CHARACTER      TRANS

           INTEGER        INFO, LDA, LDB, N, NRHS

           INTEGER        IPIV( * )

           DOUBLE         PRECISION A( LDA, * ), B( LDB, * )

PURPOSE
       DGETRS solves a system of linear equations
          A * X = B  or  A' * X = B with a general N-by-N matrix A  using  the
       LU factorization computed by DGETRF.


ARGUMENTS
       TRANS   (input) CHARACTER*1
               Specifies the form of the system of equations:
               = 'N':  A * X = B  (No transpose)
               = 'T':  A'* X = B  (Transpose)
               = 'C':  A'* X = B  (Conjugate transpose = Transpose)

       N       (input) INTEGER
               The order of the matrix A.  N >= 0.

       NRHS    (input) INTEGER
               The  number of right hand sides, i.e., the number of columns of
               the matrix B.  NRHS >= 0.

       A       (input) DOUBLE PRECISION array, dimension (LDA,N)
               The factors L and U from the factorization A =  P*L*U  as  com-
               puted by DGETRF.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,N).

       IPIV    (input) INTEGER array, dimension (N)
               The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix
               was interchanged with row IPIV(i).

       B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
               On entry, the right hand side matrix B.  On exit, the  solution
               matrix X.

       LDB     (input) INTEGER
               The leading dimension of the array B.  LDB >= max(1,N).

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value



 LAPACK routine (version 3.1)    November 2006                       DGETRS(1)
*/
void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda,
	      int *ipiv, double *b, int *ldb, int *info);



void lapack_inv (int n, const double *a,
		 double *ai)
{
  int *ipiv = (int *)malloc (sizeof (int) * n);
  double *work = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ipiv, "lapack_inv");
  CHECK_MALLOC (work, "lapack_inv");

  int i_1 = 1;
  int nn = n * n;
  dcopy_ (&nn, a, &i_1, ai, &i_1);

  int info;
  dgetrf_ (&n, &n, ai, &n, ipiv, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetrf (info = %d)\n", info);
      exit (1);
    }

  dgetri_ (&n, ai, &n, ipiv, work, &n, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetri (info = %d)\n", info);
      exit (1);
    }

  free (ipiv);
  free (work);
}

/* the version that a[n*n] is input AND output
 */
void lapack_inv_ (int n, double *a)
{
  int info;
  int *ipiv = (int *)malloc (sizeof (int) * n);
  double *work = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ipiv, "lapack_inv_");
  CHECK_MALLOC (work, "lapack_inv_");

  dgetrf_ (&n, &n, a, &n, ipiv, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetrf (info = %d)\n", info);
      exit (1);
    }

  dgetri_ (&n, a, &n, ipiv, work, &n, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetri (info = %d)\n", info);
      exit (1);
    }

  free (ipiv);
  free (work);
}

/* just solve one problem A.x = b
 * INPUT
 *  n : the order of the matrix A
 *  a[n*n] : coefficient matrix
 *  b[n]   : given vector
 * OUTPUT
 *  x[n]   : the solution
 */
void lapack_solve_lin (int n, const double *a, const double *b,
		       double *x)
{
  int *ipiv = (int *)malloc (sizeof (int) * n);
  double *work = (double *)malloc (sizeof (double) * n);
  double *lu   = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (ipiv, "lapack_solve_lin");
  CHECK_MALLOC (work, "lapack_solve_lin");
  CHECK_MALLOC (lu,   "lapack_solve_lin");

  int i_1 = 1;
  int nn = n * n;
  dcopy_ (&nn, a, &i_1, lu, &i_1);

  int info;
  dgetrf_ (&n, &n, lu, &n, ipiv, &info);
  if (info != 0)
    {
      fprintf (stderr, "singular matrix met at dgetrf (info = %d)\n", info);
      exit (1);
    }

  dcopy_ (&n, b, &i_1, x, &i_1);

  char trans = 'T'; // fortran's memory allocation is transposed
  dgetrs_ (&trans, &n, &i_1,
	   lu, &n,
	   ipiv,
	   x, &n, &info);
  if (info != 0)
    {
      fprintf (stderr, "failed in dgetrs (info = %d)\n", info);
    }

  free (ipiv);
  free (work);
  free (lu);
}

