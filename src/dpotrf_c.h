/* header file for dpotrf_c.c --
 * C wrappers for LAPACK's dpotrf_()
 * Copyright (C) 2007 Kengo Ichiki <kichiki@uwo.ca>
 * $Id: dpotrf_c.h,v 1.1 2007/10/27 22:55:08 kichiki Exp $
 */
#ifndef	_DPOTRF_C_H_
#define	_DPOTRF_C_H_


void dpotrf_wrap (int n, const double *a, double *l);
void dpotf2_wrap (int n, const double *a, double *l);


#endif /* !_DPOTRF_C_H_ */
