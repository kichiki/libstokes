/* header file for dpotrf_c.c --
 * C wrappers for LAPACK's dpotrf_()
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dpotrf_c.h,v 1.2 2007/10/31 06:05:18 kichiki Exp $
 */
#ifndef	_DPOTRF_C_H_
#define	_DPOTRF_C_H_


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
dpotrf_wrap (int n, const double *a, double *l);

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
dpotf2_wrap (int n, const double *a, double *l);


#endif /* !_DPOTRF_C_H_ */
