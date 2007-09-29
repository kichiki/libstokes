/* header file for chebyshev.c --
 * Chebyshev polynomial routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: chebyshev.h,v 1.1 2007/09/29 20:20:47 kichiki Exp $
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
#ifndef	_CHEBYSHEV_H_
#define	_CHEBYSHEV_H_


/* calc Chebyshev coefficient
 * INPUT
 *  ncheb  : number of coefficients
 *  f(x)   : function
 *  x0, x1 : the range of the function
 * OUTPUT
 *  a[ncheb] :
 */
void
chebyshev_coef (int ncheb, double (*func)(double), double x0, double x1,
		double *a);


/* eval Chebyshev polynomials by Clenshaw's algorithm
 * INPUT
 *  ncheb    : number of coefficients
 *  a[ncheb] : Chebyshev coefficients
 *  x0, x1   : the range of the function
 *  x        : 
 * OUTPUT
 *  (returned value) : = sum a_i C_i(x),
 *           where C_i is the Chebyshev polynomial.
 */
double
chebyshev_eval (int ncheb, const double *a, double x0, double x1, double x);

/* eval Chebyshev polynomials
 * INPUT
 *  ncheb    : number of coefficients
 *  a[ncheb] : Chebyshev coefficients
 *  x0, x1   : the range of the function
 *  x        : 
 * OUTPUT
 *  (returned value) : = sum a_i C_i(x),
 *           where C_i is the Chebyshev polynomial.
 */
double
chebyshev_eval_ (int ncheb, const double *a, double x0, double x1, double x);


/* eval vector by Chebyshev polynomials with Clenshaw's algorithm
 * INPUT
 *  ncheb         : number of coefficients
 *  a[ncheb]      : Chebyshev coefficients
 *  n             : dimension of the vector
 *  y[n]          : vector multiplied to the matrix M
 *  eign0, eign1  : the min and max of the eigenvalues of M
 *  atimes(n,b,x) : routine to calculate b = M.x
 *  user_data     : pointer for atimes()
 * OUTPUT
 *  z[n]          : vector approximation of (M)^{-1/2}.y
 */
void
chebyshev_eval_atimes (int ncheb, const double *a,
		       int n, const double *y,
		       double *z,
		       double eign0, double eign1,
		       void (*atimes)(int, const double *, double *, void *),
		       void *user_data);

/* estimate the error of Chebyshev approximation
 * INPUT
 *  n             : dimension of the vector
 *  y[n]          : vector given to chebyshev_eval_atimes()
 *  z[n]          : vector obtained by chebyshev_eval_atimes()
 *  atimes(n,b,x) : routine to calculate b = M.x
 *  user_data     : pointer for atimes()
 * OUTPUT
 *  (returned value) : error^2 := |zz - yy| / zz,
 *                     where zz = z . M . z, and yy = y . y.
 */
double
chebyshev_error_minvsqrt (int n, const double *y, const double *z,
			  void (*atimes)(int, const double *, double *, void *),
			  void *user_data);


#endif /* !_CHEBYSHEV_H_ */
