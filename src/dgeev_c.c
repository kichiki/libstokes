/* C wrappers for LAPACK's dgeev_()
 * Copyright (C) 2005 Kengo Ichiki <kichiki@uwo.ca>
 * $Id: dgeev_c.c,v 1.1 2007/09/29 20:38:32 kichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h>

/*
DGEEV(l)                               )                              DGEEV(l)



NAME
       DGEEV - compute for an N-by-N real nonsymmetric matrix A, the eigenval-
       ues and, optionally, the left and/or right eigenvectors

SYNOPSIS
       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,  LDVR,
                         WORK, LWORK, INFO )

           CHARACTER     JOBVL, JOBVR

           INTEGER       INFO, LDA, LDVL, LDVR, LWORK, N

           DOUBLE        PRECISION  A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
                         WI( * ), WORK( * ), WR( * )

PURPOSE
       DGEEV computes for an N-by-N real nonsymmetric matrix A, the  eigenval-
       ues  and,  optionally,  the  left and/or right eigenvectors.  The right
       eigenvector v(j) of A satisfies
                        A * v(j) = lambda(j) * v(j)
       where lambda(j) is its eigenvalue.
       The left eigenvector u(j) of A satisfies
                     u(j)**H * A = lambda(j) * u(j)**H
       where u(j)**H denotes the conjugate transpose of u(j).

       The computed eigenvectors are normalized to have Euclidean  norm  equal
       to 1 and largest component real.


ARGUMENTS
       JOBVL   (input) CHARACTER*1
               = 'N': left eigenvectors of A are not computed;
               = 'V': left eigenvectors of A are computed.

       JOBVR   (input) CHARACTER*1
               = 'N': right eigenvectors of A are not computed;
               = 'V': right eigenvectors of A are computed.

       N       (input) INTEGER
               The order of the matrix A. N >= 0.

       A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
               On  entry,  the N-by-N matrix A.  On exit, A has been overwrit-
               ten.

       LDA     (input) INTEGER
               The leading dimension of the array A.  LDA >= max(1,N).

       WR      (output) DOUBLE PRECISION array, dimension (N)
       WI      (output) DOUBLE PRECISION array, dimension (N)
               WR  and
               WI  contain  the real and imaginary parts, respectively, of the
               computed eigenvalues.  Complex conjugate pairs  of  eigenvalues
               appear  consecutively  with  the eigenvalue having the positive
               imaginary part first.

       VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
               If JOBVL = 'V', the left eigenvectors u(j) are stored one after
               another in the columns of VL, in the same order as their eigen-
               values.  If JOBVL = 'N', VL is not referenced.  If the j-th ei-
               genvalue  is  real, then u(j) = VL(:,j), the j-th column of VL.
               If the j-th and (j+1)-st eigenvalues form a  complex  conjugate
               pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
               u(j+1) = VL(:,j) - i*VL(:,j+1).

       LDVL    (input) INTEGER
               The  leading  dimension of the array VL.  LDVL >= 1; if JOBVL =
               'V', LDVL >= N.

       VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
               If JOBVR = 'V', the right  eigenvectors  v(j)  are  stored  one
               after  another in the columns of VR, in the same order as their
               eigenvalues.  If JOBVR = 'N', VR is not referenced.  If the  j-
               th  eigenvalue is real, then v(j) = VR(:,j), the j-th column of
               VR.  If the j-th and (j+1)-st eigenvalues form a complex conju-
               gate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
               v(j+1) = VR(:,j) - i*VR(:,j+1).

       LDVR    (input) INTEGER
               The  leading  dimension of the array VR.  LDVR >= 1; if JOBVR =
               'V', LDVR >= N.

       WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
               On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

       LWORK   (input) INTEGER
               The dimension of the array WORK.  LWORK >= max(1,3*N),  and  if
               JOBVL  =  'V'  or  JOBVR = 'V', LWORK >= 4*N.  For good perfor-
               mance, LWORK must generally be larger.

               If LWORK = -1, then a workspace query is assumed;  the  routine
               only  calculates  the  optimal  size of the WORK array, returns
               this value as the first entry of the WORK array, and  no  error
               message related to LWORK is issued by XERBLA.

       INFO    (output) INTEGER
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an illegal value.
               >  0:   if INFO = i, the QR algorithm failed to compute all the
               eigenvalues, and no eigenvectors have been  computed;  elements
               i+1:N of WR and WI contain eigenvalues which have converged.



LAPACK version 3.0               15 June 2000                         DGEEV(l)
 */
void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
	    int *lda, double *wr, double *wi, double *vl,
	    int *ldvl, double *vr, int *ldvr, double *work,
	    int *lwork, int *info);

/*
 * INPUT
 *  a[i*n+j] := A_{ij}
 * OUTPUT
 *  wr[n],wi[n] : real and imaginary part of the eigenvalues
 *  v[n*n] : eigenvectors
 *           k-th eigen vector (corresponding to l[k])
 *           is v[k*n+i] with i = 0 ~ n-1.
 */
void dgeev_wrap (int n, const double *a,
		 double *wr, double *wi, double *v)
{
  char jobvl = 'N';
  char jobvr = 'V';
  int lda;
  int ldvl;
  int ldvr;
  int lwork;    /* >= 4*n for 'V' case */
  int info;

  double *vl;   /* vl(ldvl,n) */
  double *work; /* work(lwork) */

  double *aa;
  int i, j;

  lda = n;
  ldvl = 1;
  ldvr = n;
  lwork = 4 * n;

  aa = malloc (sizeof (double) * n*n);
  vl = malloc (sizeof (double) * ldvl * n);
  work = malloc (sizeof (double) * lwork);

  /* AA = A^t (for fortran handling) */
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  aa[i*n+j] = a[j*n+i];
	}
    }

  dgeev_(&jobvl, &jobvr,
	 &n, aa, &lda,
	 wr, wi,
	 vl, &ldvl,
	 v, &ldvr,
	 work, &lwork,
	 &info);

  if (info != 0)
    {
      fprintf (stderr, "info = %d\n", info);
    }

  free (aa);
  free (vl);
  free (work);
}

/* obtian both left and right eigen vectors
 * INPUT
 *  a[i*n+j] := A_{ij}
 * OUTPUT
 *  wr[n],wi[n] : real and imaginary part of the eigenvalues
 *  vl[n*n],vr[n*n] : left- and right-eigenvectors
 *                    k-th eigen vector (corresponding to wr[k],wi[k])
 *                    is vl[k*n+i] and vr[k*n+i] with i = 0 ~ n-1.
 */
void dgeev_wrap_ (int n, const double *a,
		  double *wr, double *wi, double *vl, double *vr)
{
  char jobvl = 'V';
  char jobvr = 'V';
  int lda;
  int ldvl;
  int ldvr;
  int lwork;    /* >= 4*n for 'V' case */
  int info;

  double *work; /* work(lwork) */

  double *aa;
  int i, j;

  lda = n;
  ldvl = n;
  ldvr = n;
  lwork = 4 * n;

  aa = malloc (sizeof (double) * n*n);
  work = malloc (sizeof (double) * lwork);

  /* AA = A^t (for fortran handling) */
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  aa[i*n+j] = a[j*n+i];
	}
    }

  dgeev_(&jobvl, &jobvr,
	 &n, aa, &lda,
	 wr, wi,
	 vl, &ldvl,
	 vr, &ldvr,
	 work, &lwork,
	 &info);

  if (info != 0)
    {
      fprintf (stderr, "info = %d\n", info);
    }

  free (aa);
  free (work);
}
