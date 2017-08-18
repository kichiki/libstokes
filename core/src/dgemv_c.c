/* C wrappers for BLAS' dgemv_()
 * Copyright (C) 2005-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgemv_c.c,v 1.1 2006/10/22 22:43:44 kichiki Exp $
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

/*
DGEMV(l)                         BLAS routine                         DGEMV(l)



NAME
       DGEMV  - perform one of the matrix-vector operations   y := alpha*A*x +
       beta*y, or y := alpha*A'*x + beta*y,

SYNOPSIS
       SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

           DOUBLE       PRECISION ALPHA, BETA

           INTEGER      INCX, INCY, LDA, M, N

           CHARACTER*1  TRANS

           DOUBLE       PRECISION A( LDA, * ), X( * ), Y( * )

PURPOSE
       DGEMV  performs one of the matrix-vector operations

       where  alpha and beta are scalars, x and y are vectors and A is an m by
       n matrix.


PARAMETERS
       TRANS  - CHARACTER*1.
              On entry, TRANS specifies the operation to be performed as  fol-
              lows:

              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.

              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.

              Unchanged on exit.

       M      - INTEGER.
              On  entry,  M  specifies  the number of rows of the matrix A.  M
              must be at least zero.  Unchanged on exit.

       N      - INTEGER.
              On entry, N specifies the number of columns of the matrix A.   N
              must be at least zero.  Unchanged on exit.

       ALPHA  - DOUBLE PRECISION.
              On  entry, ALPHA specifies the scalar alpha.  Unchanged on exit.

       A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
              Before entry, the leading m by n part of the array A  must  con-
              tain the matrix of coefficients.  Unchanged on exit.

       LDA    - INTEGER.
*/
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
	    int *lda, double *x, int *incx,
	    double *beta, double *y, int *incy);

/* y := alpha*A*x + beta*y.
 * INPUT
 *  m : the number of rows   of the matrix A, that is, y[m]
 *  n : the number of colums of the matrix A, that is, x[n]
 *      A is 'm by n' matrix.
 *  dgemv_() is implemented in FORTRAN, so that,
 *  for 'N' case,
 *     y[i] = sum_j A[i,j] x[j]
 *          = sum_J a[I+m*J] * x[J], where I:=i-1, etc.
 *          = sum_J a[J*m+I] * x[J]
 *          = A_C[J,I] * x[J]
 *  for 'T' case,
 *     y[i] = sum_j A[j,i] x[j]
 *          = sum_J a[J+n*I] * x[J], where I:=i-1, etc.
 *          = sum_J a[I*n+J] * x[J]
 *          = A_C[I,J] * x[J]
 *  NOTE, in this case, m and n also should be exchanged...
 */
void dgemv_wrap (int m, int n, double alpha, double *a,
		 double *x, double beta,
		 double *y)
{
  char trans = 'T'; /* fortran's memory allocation is transposed */
  int one = 1;

  dgemv_ (&trans, &n, &m, &alpha, a, &n, x, &one, &beta, y, &one);
}
