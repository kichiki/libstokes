/* C wrappers for LAPACK's dpotf2_()
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dpotrf_c.c,v 1.4 2007/11/04 02:57:09 kichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h>

/*
DPOTRF(1)                LAPACK routine (version 3.1)                DPOTRF(1)



NAME
       DPOTRF  - the Cholesky factorization of a real symmetric positive defi-
       nite matrix A

SYNOPSIS
       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )

           CHARACTER      UPLO

           INTEGER        INFO, LDA, N

           DOUBLE         PRECISION A( LDA, * )

PURPOSE
       DPOTRF computes the Cholesky factorization of a real symmetric positive
       definite matrix A.

       The factorization has the form
          A = U**T * U,  if UPLO = 'U', or
          A = L  * L**T,  if UPLO = 'L',
       where U is an upper triangular matrix and L is lower triangular.

       This is the block version of the algorithm, calling Level 3 BLAS.


ARGUMENTS
       UPLO    (input) CHARACTER*1
               = 'U':  Upper triangle of A is stored;
               = 'L':  Lower triangle of A is stored.

       N       (input) INTEGER
               The order of the matrix A.  N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On  entry,  the symmetric matrix A.  If UPLO = 'U', the leading
               N-by-N upper triangular part of A contains the upper triangular
               part of the matrix A, and the strictly lower triangular part of
               A is not referenced.  If UPLO = 'L', the leading  N-by-N  lower
               triangular  part of A contains the lower triangular part of the
               matrix A, and the strictly upper triangular part of  A  is  not
               referenced.

               On  exit, if INFO = 0, the factor U or L from the Cholesky fac-
               torization A = U**T*U or A = L*L**T.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,N).

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value
               > 0:  if INFO = i, the leading minor of order i is not positive
               definite, and the factorization could not be completed.



 LAPACK routine (version 3.1)    November 2006                       DPOTRF(1)
*/
void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);


/*
DPOTF2(1)                LAPACK routine (version 3.1)                DPOTF2(1)



NAME
       DPOTF2  - the Cholesky factorization of a real symmetric positive defi-
       nite matrix A

SYNOPSIS
       SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )

           CHARACTER      UPLO

           INTEGER        INFO, LDA, N

           DOUBLE         PRECISION A( LDA, * )

PURPOSE
       DPOTF2 computes the Cholesky factorization of a real symmetric positive
       definite matrix A.

       The factorization has the form
          A = U' * U ,  if UPLO = 'U', or
          A = L  * L',  if UPLO = 'L',
       where U is an upper triangular matrix and L is lower triangular.

       This is the unblocked version of the algorithm, calling Level 2 BLAS.


ARGUMENTS
       UPLO    (input) CHARACTER*1
               Specifies  whether  the  upper  or lower triangular part of the
               symmetric matrix A is stored.  = 'U':  Upper triangular
               = 'L':  Lower triangular

       N       (input) INTEGER
               The order of the matrix A.  N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On entry, the symmetric matrix A.  If UPLO = 'U', the leading n
               by  n  upper triangular part of A contains the upper triangular
               part of the matrix A, and the strictly lower triangular part of
               A  is  not referenced.  If UPLO = 'L', the leading n by n lower
               triangular part of A contains the lower triangular part of  the
               matrix  A,  and  the strictly upper triangular part of A is not
               referenced.

               On exit, if INFO = 0, the factor U or L from the Cholesky  fac-
               torization A = U'*U  or A = L*L'.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,N).

       INFO    (output) INTEGER
               = 0: successful exit
               < 0: if INFO = -k, the k-th argument had an illegal value
               >  0: if INFO = k, the leading minor of order k is not positive
               definite, and the factorization could not be completed.



 LAPACK routine (version 3.1)    November 2006                       DPOTF2(1)
 */
void dpotf2_(char *uplo, int *n, double *a, int *lda, int *info);


/*
 * INPUT
 *  a[i*n+j] := A_{ij}
 * OUTPUT
 *  wr[n],wi[n] : real and imaginary part of the eigenvalues
 *  v[n*n] : eigenvectors
 *           k-th eigen vector (corresponding to l[k])
 *           is v[k*n+i] with i = 0 ~ n-1.
 *  returned value : info
 */
int
dpotrf_wrap (int n, const double *a, double *l)
{
  char uplo = 'L';

  int i;
  for (i = 0; i < n * n; i ++)
    {
      l[i] = a[i];
    }

  int info;
  dpotrf_(&uplo, &n, l, &n, &info);

  /*
  if (info != 0)
    {
      fprintf (stderr, "dpotrf : info = %d\n", info);
      exit (1);
    }
  */

  for (i = 0; i < n; i ++)
    {
      int j;
      for (j = i+1; j < n; j ++)
	{
	  /* convert from FORTRAN from to C form
	   */
	  l[j*n+i] = l[i*n+j];
	  l[i*n+j] = 0.0;
	}
    }

  return (info);
}

/*
 * INPUT
 *  a[i*n+j] := A_{ij}
 * OUTPUT
 *  wr[n],wi[n] : real and imaginary part of the eigenvalues
 *  v[n*n] : eigenvectors
 *           k-th eigen vector (corresponding to l[k])
 *           is v[k*n+i] with i = 0 ~ n-1.
 *  returned value : info
 */
int
dpotf2_wrap (int n, const double *a, double *l)
{
  char uplo = 'L';

  int i;
  for (i = 0; i < n * n; i ++)
    {
      l[i] = a[i];
    }

  int info;
  dpotf2_(&uplo, &n, l, &n, &info);

  /*
  if (info != 0)
    {
      fprintf (stderr, "dpotf2 : info = %d\n", info);
      exit (1);
    }
  */

  for (i = 0; i < n; i ++)
    {
      int j;
      for (j = i+1; j < n; j ++)
	{
	  /* convert from FORTRAN from to C form
	   */
	  l[j*n+i] = l[i*n+j];
	  l[i*n+j] = 0.0;
	}
    }

  return (info);
}
