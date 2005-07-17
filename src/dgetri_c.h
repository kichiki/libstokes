/* header file for dgetri.c --
 * C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005 Kengo Ichiki <kichiki@uwo.ca>
 * $Id: dgetri_c.h,v 1.2 2005/07/17 18:58:46 ichiki Exp $
 */
#ifndef	_DGETRI_C_H_
#define	_DGETRI_C_H_

void lapack_inv (int n, const double *a,
		 double *ai);

#endif /* !_DGETRI_C_H_ */
