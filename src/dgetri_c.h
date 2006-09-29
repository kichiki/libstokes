/* header file for dgetri.c --
 * C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgetri_c.h,v 1.3 2006/09/29 03:33:18 ichiki Exp $
 */
#ifndef	_DGETRI_C_H_
#define	_DGETRI_C_H_

void lapack_inv (int n, const double *a,
		 double *ai);

/* the version that a[n*n] is input AND output
 */
void lapack_inv_ (int n, double *a);

#endif /* !_DGETRI_C_H_ */
