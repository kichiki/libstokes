/* test code for solve_cubic() with GSL routine poly_solve_cubic()
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include <gsl/gsl_poly.h> // gsl_poly_solve_cubic()

#include <libstokes-core.h> // ptime_ms_d()

#include "memory-check.h" // macro CHECK_MALLOC

#include "check.h" // compare_max()


double
calc_cubic (double a, double b, double c, double d, double x)
{
  double y = d + x * (c + x * (b + x * a));

  return (y);
}


void
sort3 (double *a, double *b, double *c)
{
  if ((*a) > (*b))
    {
      double x = (*a);
      (*a) = (*b);
      (*b) = x;
    }
  // now (*a) <= (*b)
  if ((*b) > (*c))
    {
      if ((*a) > (*c))
	{
	  // now (*a) <= (*b), but (*a) > (*c)
	  double x = (*a);
	  (*a) = (*c);
	  (*c) = x;
	}
      else
	{
	  // now (*a) <= (*b), but (*b) > (*c) > (*a)
	  double x = (*b);
	  (*b) = (*c);
	  (*c) = x;
	}
    }
}


/* my own implementation for solving cubic equation
 * INPUT
 * OUTPUT
 *  returned value : number of real solutions
 */
int
solve_cubic (double a, double b, double c, double d,
	     double *x1, double *x2, double *x3)
{
  if (a == 0.0) return (0); // fail

  double A = 3.0 * a;
  double A2 = A * A;
  double A3 = A2 * A;

  double b2 = b * b;

  double P = (-b2 + A * c) / A2;
  double Q = 0.5 * (b * (2.0 * b2 - 9.0 * a * c) + 3.0 * A2 * d) / A3;

  double P3 = P * P * P;
  double Q2 = Q * Q;

  double D = Q2 + P3;

  if (D <= 0.0)
    {
      double r = -Q;
      double y = sqrt (-D);

      //double R3 = r*r + y*y;
      double R3 = Q2 - D;
      R3 = pow (R3, 1.0 / 6.0);

      double t3 = atan2 (y, r);
      //t3 += 2.0 * M_PI;
      t3 /= 3.0;

      double X = R3 * cos (t3);
      double Y = R3 * sin (t3);

      *x1 = 2.0 * X - b / A;

      double sq3 = sqrt (3.0);
      *x2 = - X - b / A - sq3 * Y;
      *x3 = - X - b / A + sq3 * Y;

      return (3); // three real solutions
    }

  if (P == 0.0)
    {
      // sqD = sqrt(Q^2) = Q
      double um = -pow (2.0 * Q, 1.0 / 3.0);
      *x1 = um - b / A;
    }
  else
    {
      double sqD = sqrt (D);
      double up =  pow (-Q + sqD, 1.0 / 3.0);
      *x1 = up - P / up - b / A;

      /*
      double um = -pow (+Q + sqD, 1.0 / 3.0);
      double xm = um - P / um - b / A;
      if (fabs (xp - xm) > 1.0e-5)
	{
	  fprintf (stderr, "xp, xm = %f, %f, %e\n",
		   xp, xm, fabs (xp - xm));
	  exit (1);
	}
      */
    }

  return (1); // one real solution
}

int
check_solve_cubic (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_cubic\n");
    }

  int check = 0;
  double max = 0.0;

  int num;
  double x1, x2, x3;

  double a;
  double b;
  double c;
  double d;

  /*
// single real root
(%i8) expand((x-1)*(x^2+1));
                                 3    2
(%o8)                           x  - x  + x - 1
  */
  a = 1.0;
  b = -1.0;
  c = 1.0;
  d = -1.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 1.0, " num for x^3-x^2+x-1",
			verbose, tiny, &max);
  check += compare_max (x1, 1.0, " x1 for x^3-x^2+x-1",
			verbose, tiny, &max);

  /*
// single real root
(%i1) expand((x+1)*(x^2+2));
                                3    2
(%o1)                          x  + x  + 2 x + 2
  */
  a = 1.0;
  b = 1.0;
  c = 2.0;
  d = 2.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 1.0, " num for x^3+x^2+2x+2",
			verbose, tiny, &max);
  check += compare_max (x1, -1.0, " x1 for x^3+x^2+2x+2",
			verbose, tiny, &max);

  /* trivial case : (x-1)^3 = x^3-3x^2+3x+1
  */
  a = 1.0;
  b = -3.0;
  c = +3.0;
  d = -1.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  //sort3 (&x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 3.0, " num for (x-1)^3",
			verbose, tiny, &max);
  check += compare_max (x1, 1.0, " x1 for (x-1)^3",
			verbose, tiny, &max);
  check += compare_max (x2, 1.0, " x2 for (x-1)^3",
			verbose, tiny, &max);
  check += compare_max (x3, 1.0, " x3 for (x-1)^3",
			verbose, tiny, &max);

  /* trivial case : (x-1)^2(x+1) = (x^2-1)*(x-1) = x^3-x^2-x+1
  */
  a = 1.0;
  b = -1.0;
  c = -1.0;
  d = 1.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  //sort3 (&x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 3.0, " num for (x+1)(x-1)^2",
			verbose, tiny, &max);
  check += compare_max (x1, -1.0, " x1 for (x+1)(x-1)^2",
			verbose, tiny, &max);
  check += compare_max (x2, 1.0, " x2 for (x+1)(x-1)^2",
			verbose, tiny, &max);
  check += compare_max (x3, 1.0, " x3 for (x+1)(x-1)^2",
			verbose, tiny, &max);

  /*
(%i5) expand((x-1)*(x-2)*(x-3));
                              3      2
(%o5)                        x  - 6 x  + 11 x - 6
  */
  a = 1.0;
  b = -6.0;
  c = 11.0;
  d = -6.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  //sort3 (&x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 3.0, " num for (x-1)(x-2)(x-3)",
			verbose, tiny, &max);
  check += compare_max (x1, 1.0, " x1 for (x-1)(x-2)(x-3)",
			verbose, tiny, &max);
  check += compare_max (x2, 2.0, " x2 for (x-1)(x-2)(x-3)",
			verbose, tiny, &max);
  check += compare_max (x3, 3.0, " x3 for (x-1)(x-2)(x-3)",
			verbose, tiny, &max);

  /*
(%i6) expand((x+1)*(x+2)*(x+3));
                              3      2
(%o6)                        x  + 6 x  + 11 x + 6

(%i7) expand((x-4)*(x-5)*(x-6));
                             3       2
(%o7)                       x  - 15 x  + 74 x - 120
  */
  a = 1.0;
  b = -15.0;
  c = 74.0;
  d = -120.0;
  //num = solve_cubic (a, b, c, d, &x1, &x2, &x3);
  //sort3 (&x1, &x2, &x3);
  num = gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);

  check += compare_max ((double)num, 3.0, " num for (x-4)(x-5)(x-6)",
			verbose, tiny, &max);
  check += compare_max (x1, 4.0, " x1 for (x-4)(x-5)(x-6)",
			verbose, tiny, &max);
  check += compare_max (x2, 5.0, " x2 for (x-4)(x-5)(x-6)",
			verbose, tiny, &max);
  check += compare_max (x3, 6.0, " x3 for (x-4)(x-5)(x-6)",
			verbose, tiny, &max);


  // benchmark test
  int i;
  double t0 = ptime_ms_d ();
  for (i = 0; i < 10000; i ++)
    {
      solve_cubic (a, b, c, d, &x1, &x2, &x3);
    }
  double t1 = ptime_ms_d ();

  double t2 = ptime_ms_d ();
  for (i = 0; i < 10000; i ++)
    {
      gsl_poly_solve_cubic (b/a, c/a, d/a, &x1, &x2, &x3);
    }
  double t3 = ptime_ms_d ();
  if (verbose != 0)
    {
      fprintf (stdout, " CPU times: %f (local) vs %f (GSL)\n",
	       t1 - t0,
	       t3 - t2);
    }

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

