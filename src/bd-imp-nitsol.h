/* header file for bd-imp-nitsol.c --
 * implicit Brownian dynamics algorithms by NITSOL
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp-nitsol.h,v 1.2 2008/06/06 03:46:21 kichiki Exp $
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
#ifndef	_BD_IMP_NITSOL_H_
#define	_BD_IMP_NITSOL_H_

#include <bd-imp.h> // struct BD_imp


/* jacobian and preconditioner
 * at this time, nothing will be done for both
 */
void BD_imp_NITSOL_jacv (int *n, double *xcur, double *fcur,
			 int *ijob, double *v, double *z,
			 double *rpar, int *ipar, int *itrmjv);

/* a wrapper of the nonlinear equations for the semi-implicit algorithms
 * (both JGdP and siPC) for NITSOL
 * INTPUT
 *  xcur[n] : = (x[nm*3])          for F version
 *            = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  rpar    : (struct BD_imp *)BDimp
 *  ipar    : not used.
 * OUTPUT
 *  itrmf   : 0 means success
 *  fcur[n] := x - x0
 *           - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for JGdP, or
 *  fcur[n] := x - x0
 *           - (dt/2) * (U^pr
 *                       + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for siPC, where U^pr is the predictor obtained by Euler scheme as 
 *  U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
void
BD_imp_NITSOL_f (int *n, double *xcur, double *fcur,
		 double *rpar, int *ipar, int *itrmf);

void
BD_imp_NITSOL_wrap (struct BD_imp *BDimp,
		    double *x, double *q);


#endif /* !_BD_IMP_NITSOL_H_ */
