/* header file for dgetri.c --
 * C wrappers for LAPACK's dgetri_() and dgetrt_()
 * Copyright (C) 2005 Kengo Ichiki <kichiki@uwo.ca>
 * $Id: dgetri_c.h,v 1.1 2005/03/29 18:48:17 ichiki Exp $
 */

void lapack_inv (int n, const double *a,
		 double *ai);
