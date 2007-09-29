/* header file for dgeev.c --
 * C wrappers for LAPACK's dgeev_()
 * Copyright (C) 2005 Kengo Ichiki <kichiki@uwo.ca>
 * $Id: dgeev_c.h,v 1.1 2007/09/29 20:38:32 kichiki Exp $
 */
#ifndef	_DGEEV_C_H_
#define	_DGEEV_C_H_

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
		 double *wr, double *wi, double *v);

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
		  double *wr, double *wi, double *vl, double *vr);

#endif /* !_DGEEV_C_H_ */
