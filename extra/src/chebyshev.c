/* Chebyshev polynomial routines
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: chebyshev.c,v 1.3 2008/05/24 05:58:31 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory-check.h" // macro CHECK_MALLOC


/* calc Chebyshev coefficient
 * INPUT
 *  ncheb  : number of coefficients
 *  func(x)   : function
 *  x0, x1 : the range of the function
 * OUTPUT
 *  a[ncheb] :
 */
void
chebyshev_coef (int ncheb, double (*func)(double), double x0, double x1,
		double *a)
{
  double y0 = 0.5 * (x1 - x0);
  double y1 = 0.5 * (x1 + x0);
  double twon = 2.0 / (double)ncheb;

  int j;
  for (j = 0; j < ncheb; j ++)
    {
      double x = 0.0;
      int k;
      for (k = 1; k <= ncheb; k ++)
	{
	  double pikn = M_PI * ((double)k - 0.5) / (double)ncheb;
	  /* y is in [-1,1] */
	  double y = cos (pikn);
	  // shift y to yy which is in [x0,x1]
	  double yy = y * y0 + y1;
	  double z = cos ((double)j * pikn);
	  x += func(yy) * z;
	}
      a[j] = twon * x;
    }

  // adjust a[0]
  a[0] *= 0.5;
}


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
chebyshev_eval (int ncheb, const double *a, double x0, double x1, double x)
{
  double y0 = 0.5 * (x1 - x0);
  double y1 = 0.5 * (x1 + x0);
  // y is in [-1,1]
  double y = (x - y1) / y0;

  double b[3]; // we only need three coefficients
  // zero clear (for sure)
  b [0] = 0.0;
  b [1] = 0.0;
  b [2] = 0.0;
  int i, j;
  int i2, i1, i0;
  for (i = 0; i < ncheb - 1; i ++) // i runs from 0 to ncheb-2.
    {
      i2 = (i+2) % 3;
      i1 = (i+1) % 3;
      i0 = (i  ) % 3;
      j = ncheb - 1 - i; // runs from (ncheb-1) to 1.
      b[i2] = 2 * y * b[i1] - b[i0] + a[j];
    }
  i = ncheb - 1;
  //i2 = (i+2) % 3;
  i1 = (i+1) % 3;
  i0 = (i  ) % 3;
  j = 0;
  double f = y * b[i1] - b[i0] + a[j];

  return (f);
}

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
chebyshev_eval_ (int ncheb, const double *a, double x0, double x1, double x)
{
  double y0 = 0.5 * (x1 - x0);
  double y1 = 0.5 * (x1 + x0);
  // y is in [-1,1]
  double y = (x - y1) / y0;

  double f = 0.0;
  double c; // chebyshev polynomials
  double c0, c1;

  // k = 0
  c = 1.0;
  f += a[0];

  // k = 1;
  c0 = c;
  c = y;
  f += a[1] * c;

  int k;
  for (k = 2; k < ncheb; k ++)
    {
      c1 = c0;
      c0 = c;
      c  = 2.0 * y * c0 - c1;
      f += a[k] * c;
    }

  return (f);
}


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
		       void *user_data)
{
  // x := M . y
  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "chebyshev_eval_atimes");

  double da = 2.0 / (eign1 - eign0);
  double db = - (eign1 + eign0) / (eign1 - eign0);


  double *b = (double *)malloc (sizeof (double) * 3 * n);
  CHECK_MALLOC (b, "chebyshev_eval_atimes");
  int k;
  for (k = 0; k < n; k ++)
    {
      b[0*n+k] = 0.0;
      b[1*n+k] = 0.0;
    }
  int i, j;
  int i2, i1, i0;
  for (i = 0; i < ncheb - 1; i ++) // i runs from 0 to ncheb-2.
    {
      i2 = ((i+2) % 3) * n;
      i1 = ((i+1) % 3) * n;
      i0 = ((i  ) % 3) * n;
      j = ncheb - 1 - i; // runs from (ncheb-1) to 1.

      // x = M . b[i1]
      atimes (n, b + i1, x, user_data);
      for (k = 0; k < n; k ++)
	{
	  b[i2 + k] = 2 * da * x[k] // M. b[i1]
	    + 2 * db * b[i1 + k]
	    - b[i0 + k]
	    + a[j] * y[k];
	}
    }

  // for the last term
  i = ncheb - 1;
  //i2 = ((i+2) % 3) * n;
  i1 = ((i+1) % 3) * n;
  i0 = ((i  ) % 3) * n;
  j = 0;

  // x = M . b[i1]
  atimes (n, b + i1, x, user_data);
  // z[] is the final result
  for (k = 0; k < n; k ++)
    {
      z[k] = da * x[k] // M. b[i1]
	+ db * b[i1 + k]
	- b[i0 + k]
	+ a[j] * y[k];
    }

  // house-keeping
  free (b);
  free (x);
}

/* estimate the error of Chebyshev approximation of sqrt of Minv
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
			  void *user_data)
{
  // x := M . z
  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "chebyshev_error_minvsqrt");
  atimes (n, z, x, user_data); // x = M.z

  int k;
  double yy = 0.0;
  double zz = 0.0;
  for (k = 0; k < n; k ++)
    {
      yy += y[k] * y[k];
      zz += z[k] * x[k];
    }
  // house-keeping
  free (x);

  if (zz <= 0.0)
    {
      fprintf (stderr, "negative norm! %e\n", zz);
      exit (1);
    }

  return (fabs (zz - yy) / zz);
}

/* estimate the error of Chebyshev approximation of sqrt() of R
 * INPUT
 *  n             : dimension of the vector
 *  y[n]          : vector given to chebyshev_eval_atimes()
 *  z[n]          : vector obtained by chebyshev_eval_atimes()
 *  atimes(n,b,x) : routine to calculate b = R.x
 *  user_data     : pointer for atimes()
 * OUTPUT
 *  (returned value) : error^2 := |zz - yy| / zz,
 *                     where zz = z . M . z, and yy = y . y.
 */
double
chebyshev_error_Rsqrt (int n, const double *y, const double *z,
		       void (*atimes)(int, const double *, double *, void *),
		       void *user_data)
{
  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "chebyshev_error_Rsqrt");
  atimes (n, y, x, user_data); // x = R.y

  int k;
  double yy = 0.0;
  double zz = 0.0;
  for (k = 0; k < n; k ++)
    {
      yy += y[k] * x[k]; // yy = y . R . y
      zz += z[k] * z[k]; // zz = z . z = (R^{1/2}.y).(R^{1/2}.y)
    }
  // house-keeping
  free (x);
  //fprintf (stderr, "# yy = %e, zz = %e\n", yy, zz);

  if (yy <= 0.0)
    {
      fprintf (stderr, "negative norm! %e\n", yy);
      exit (1);
    }

  return (fabs (yy - zz) / yy);
}
