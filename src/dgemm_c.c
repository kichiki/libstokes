/* C wrappers for BLAS' dgemm_()
 * Copyright (C) 2005-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgemm_c.c,v 1.1 2006/10/22 22:43:44 kichiki Exp $
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
#include <stdlib.h>

/*
DGEMM(l)                         BLAS routine                         DGEMM(l)



NAME
       DGEMM  - perform one of the matrix-matrix operations   C := alpha*op( A
       )*op( B ) + beta*C,

SYNOPSIS
       SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K,  ALPHA,  A,  LDA,  B,  LDB,
                        BETA, C, LDC )

           CHARACTER*1  TRANSA, TRANSB

           INTEGER      M, N, K, LDA, LDB, LDC

           DOUBLE       PRECISION ALPHA, BETA

           DOUBLE       PRECISION A( LDA, * ), B( LDB, * ), C( LDC, * )

PURPOSE
       DGEMM  performs one of the matrix-matrix operations

       where  op( X ) is one of

          op( X ) = X   or   op( X ) = X',

       alpha  and  beta are scalars, and A, B and C are matrices, with op( A )
       an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.


PARAMETERS
       TRANSA - CHARACTER*1.  On entry, TRANSA specifies the form of op(  A  )
       to be used in the matrix multiplication as follows:

       TRANSA = 'N' or 'n',  op( A ) = A.

       TRANSA = 'T' or 't',  op( A ) = A'.

       TRANSA = 'C' or 'c',  op( A ) = A'.

       Unchanged on exit.

       TRANSB  -  CHARACTER*1.  On entry, TRANSB specifies the form of op( B )
       to be used in the matrix multiplication as follows:

       TRANSB = 'N' or 'n',  op( B ) = B.

       TRANSB = 'T' or 't',  op( B ) = B'.

       TRANSB = 'C' or 'c',  op( B ) = B'.

       Unchanged on exit.

       M      - INTEGER.
              On entry,  M  specifies  the number  of rows  of the  matrix op(
              A  )   and  of  the   matrix   C.   M   must  be at least  zero.
              Unchanged on exit.

       N      - INTEGER.
              On entry,  N  specifies the number  of columns of the matrix op(
              B  )  and  the  number  of columns of the matrix C. N must be at
              least zero.  Unchanged on exit.

       K      - INTEGER.
              On entry,  K  specifies  the number of columns of the matrix op(
              A  )  and the number of rows of the matrix op( B ). K must be at
              least  zero.  Unchanged on exit.

       ALPHA  - DOUBLE PRECISION.
              On entry, ALPHA specifies the scalar alpha.  Unchanged on  exit.

       A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
              k   when   TRANSA  =  'N' or 'n',  and is  m  otherwise.  Before
              entry with  TRANSA = 'N' or 'n',  the leading  m by  k  part  of
              the array  A  must contain the matrix  A,  otherwise the leading
              k by m  part of the  array   A   must  contain   the  matrix  A.
              Unchanged on exit.

       LDA    - INTEGER.
              On  entry, LDA specifies the first dimension of A as declared in
              the calling (sub) program. When  TRANSA = 'N' or  'n'  then  LDA
              must  be  at least  max( 1, m ), otherwise  LDA must be at least
              max( 1, k ).  Unchanged on exit.

       B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
              n  when  TRANSB = 'N' or 'n',  and  is   k   otherwise.   Before
              entry  with   TRANSB  = 'N' or 'n',  the leading  k by n part of
              the array  B  must contain the matrix  B,  otherwise the leading
              n  by  k   part  of  the  array   B  must contain  the matrix B.
              Unchanged on exit.

       LDB    - INTEGER.
              On entry, LDB specifies the first dimension of B as declared  in
              the  calling  (sub)  program. When  TRANSB = 'N' or 'n' then LDB
              must be at least  max( 1, k ), otherwise  LDB must be  at  least
              max( 1, n ).  Unchanged on exit.

       BETA   - DOUBLE PRECISION.
              On  entry,   BETA   specifies  the scalar  beta.  When  BETA  is
              supplied as zero then C need not be set on input.  Unchanged  on
              exit.

       C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
              Before  entry,  the  leading   m by n  part of the array  C must
              contain the matrix  C,  except when  beta   is  zero,  in  which
              case  C  need  not  be  set on entry.  On exit, the array  C  is
              overwritten by the  m by n  matrix ( alpha*op( A  )*op(  B  )  +
              beta*C ).

       LDC    - INTEGER.
              On  entry, LDC specifies the first dimension of C as declared in
              the  calling  (sub)  program.   LDC  must  be  at  least max( 1,
              m ).  Unchanged on exit.

              Level 3 Blas routine.

              --  Written on 8-February-1989.  Jack Dongarra, Argonne National
              Laboratory.  Iain Duff, AERE Harwell.  Jeremy Du Croz, Numerical
              Algorithms  Group  Ltd.   Sven  Hammarling, Numerical Algorithms
              Group Ltd.







BLAS routine                    16 October 1992                       DGEMM(l)
*/
void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	    double *alpha, double *a, int *lda,
	    double *b, int *ldb,
	    double *beta, double *c, int *ldc);

/* C := alpha*op( A )*op( B ) + beta*C
 * INPUT
 *  m : the number of rows   of the matrix op(A),C => op(A)[ij],i=0~m.
 *  n : the number of colums of the matrix op(B),C => op(B)[ij],j=0~n.
 *  k : the number of colums of the matrix op(A) and
 *      the number of rows   of the matrix op(B)   => op(A)[ij],j=0~k.
 *      that is, C[m,n] = alpha A[m,k] . B[k,n] + beta * C[m,n]
 *  op(A) is 'm by k' matrix => op(A)[m,k]
 *  op(B)    'k by n' matrix => op(B)[k,n]
 *     C     'm by n' matrix =>    C [m,n],
 *            where n is the number of columns (colum means 'tate-block')
 *            and   m is the number of rows.   (row   means 'line')
 *  dgemm_() is implemented in FORTRAN, so that,
 *  for 'N' case,
 *     C[i,j] = sum_l A[i,l] B[l,j]
 *            = sum_L a[I+m*L] * b[L+k*J], where I:=i-1, etc.
 *            = sum_L a[L*m+I] * b[J*k+L]
 *            = A_C[L,I] * B_C[J,L]
 *  for 'T' case,
 *     C[i,j] = sum_l A[l,i] B[j,l]
 *            = sum_L a[L+k*I] * b[J+n*L] where I:=i-1, etc.
 *            = sum_L a[I*k+L] * b[L*n+J]
 *            = A_C[I,L] * B_C[L,J]
 *  for both cases,
 *     C[i,j] = c[I+m*J] = C[J*m+I] = C_C[J,I]
 * SUMMARY in C:
 *   for 'N' : C[I,J] = A[L,J] * B[I,L]
 *   for 'T' : C[I,J] = A[J,L] * B[L,I]
 * here, call dgemm_() in 'T' mode, and then, transpose the result.
 */
void dgemm_wrap (int m, int n, int k,
		 double alpha, double *a, double *b,
		 double beta, double *c)
{
  char trans = 'T'; /* fortran's memory allocation is transposed */
  double *d = NULL;
  int i,j;

  d = (double *) malloc (sizeof (double) * m*n);
  for (i = 0; i < m; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  d[j*m+i] = c[i*n+j];
	}
    }

  dgemm_ (&trans, &trans, &m, &n, &k,
	  &alpha, a, &k, b, &n,
	  &beta, d, &m);

  for (i = 0; i < m; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  c[i*n+j] = d[j*m+i];
	}
    }
  free (d);
}
