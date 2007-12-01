/* header file for dgetri.c --
 * C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgetri_c.h,v 1.4 2007/12/01 18:42:23 kichiki Exp $
 */
#ifndef	_DGETRI_C_H_
#define	_DGETRI_C_H_

void lapack_inv (int n, const double *a,
		 double *ai);

/* the version that a[n*n] is input AND output
 */
void lapack_inv_ (int n, double *a);

/* just solve one problem A.x = b
 * INPUT
 *  n : the order of the matrix A
 *  a[n*n] : coefficient matrix
 *  b[n]   : given vector
 * OUTPUT
 *  x[n]   : the solution
 */
void lapack_solve_lin (int n, const double *a, const double *b,
		       double *x);


#endif /* !_DGETRI_C_H_ */
