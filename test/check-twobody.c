/* test code for twobody.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-twobody.c,v 1.1 2007/04/25 05:56:11 kichiki Exp $
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
#include "memory-check.h"

#include <twobody.h>
#include <two-body-res.h>
#include <bench.h> // ptime_ms_d()

#include "check.h" // compare()


/* calc scalar functions of resistance problem by lub form
 * all-in-one form (to reduce calculating the same parameters)
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  far [22] : scalar functions
 *      0, 1 : (XA11, XA12)
 *      2, 3 : (YA11, YA12)
 *      4, 5 : (YB11, YB12)
 *      6, 7 : (XC11, XC12)
 *      8  9 : (YC11, YC12)
 *     10,11 : (XG11, XG12)
 *     12,13 : (YG11, YG12)
 *     14,15 : (YH11, YH12)
 *     16,17 : (XM11, XM12)
 *     18,19 : (YM11, YM12)
 *     20,21 : (ZM11, ZM12)
 */
void twobody_lub_ (int version,
		   int n, double l, double s,
		   double *lub)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XA_lub");

  double l1 = 1.0 + l;
  double l12 = l1 * l1;
  double l13 = l12 * l1;
  double l14 = l12 * l12;

  double s1, s2, s3, s4, s5;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;

  double ls1, ls2, ls5;
  ls1 = log (s1);
  ls2 = log (s2);

  double g1, g2, g3, g4, g5, g6;

  double s1l;
  double s1lm;
  double s21;
  double s2m;

  double a;
  int m;

  // XA
  twobody_XA (n, l, f);

  lub [0] = 0.0;
  lub [1] = 0.0;

  g1 = 2.0 * l * l / l13;
  g2 = 0.2 * l * (1.0 + l * (7.0 + l)) / l13;
  g3 = (1.0 + l * (18.0 + l * (- 29.0 + l * (18.0 + l)))) / l13 / 42.0;

  lub [0] = g1 / s1
    - g2 * ls1
    - g3 * s1 * ls1
    + f[0]
    - g1;
  // next is m = 2

  lub [1] = 2.0 / s * g1 / s1
    + g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = -g1 - 2.0 * g2 / (double)m;
	  if (m == 2)
	    {
	      a += 4.0 * g3 / (double)(-2*m);
	    }
	  else
	    {
	      a += 4.0 * g3 / (double)(m*(m-2));
	    }
	  lub [0] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = -g1 - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  lub [1] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub [1] *= -2.0 / l1;


  // YA
  twobody_YA (n, l, f);

  lub [2] = 0.0;
  lub [3] = 0.0;

  g2 = (4.0 / 15.0) * l * (2.0 + l * (1.0 + l * 2.0)) / l13;
  g3 = (2.0 / 375.0)
    * (16.0 + l * (-45.0 + l * (58.0 + l * (-45.0 + l * 16.0)))) / l13;

  lub [2] = - g2 * ls1
    - g3 * s1 * ls1
    + f[0];
  // next is m = 2

  lub [3] = g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - 2.0 * g2 / (double)m;
	  if (m == 2)
	    {
	      a += 4.0 * g3 / (double)(-2*m);
	    }
	  else
	    {
	      a += 4.0 * g3 / (double)(m*(m-2));
	    }
	  lub [2] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  lub [3] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub [3] *= -2.0 / l1;

  // F version
  if (version == 0) return;


  // YB
  twobody_YB (n, l, f);

  lub [4] = 0.0;
  lub [5] = 0.0;

  g2 = -0.2 * l * (4.0 + l) / l12;
  g3 = -0.004 * (32.0 + l * (-33.0 + l * (83.0 + l * 43.0))) / l12;

  lub [4] = g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  lub [5] = - g2 * ls1
    - g3 * s1 * ls1;
  // next is m = 2

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - 2.0 * g2 / (double)m;
	  if (m == 2)
	    {
	      a += 4.0 * g3 / (double)(-2*m);
	    }
	  else
	    {
	      a += 4.0 * g3 / (double)(m*(m-2));
	    }
	  lub [5] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  lub [4] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub [5] *= -4.0 / l12;

  // XC
  twobody_XC (n, l, f);

  lub [6] = 0.0;
  lub [7] = 0.0;

  lub [6] = l * l / l1 * (0.5 * ls1 + ls2 / s)
    + 1.0;
  // next is m = 2

  lub [7] = - l * l / l1 * (0.5 * ls2 + (ls1 - 2.0) / s);
  // next is m = 3

  s1l = s * (1.0 + l);
  s1lm = s1l * s1l; // starting from m=2
  s21  = 2.0 / s;
  s2m  = s21 * s21; // starting from m=2
  for (m = 2; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - l * l / l1 / (double)(m * (m-1));
	  lub [6] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - l * l / l1 / (double)(m * (m-1));
	  lub [7] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  //lub [7] *= -8.0 / l1 / l1 / l1;
  lub [7] *= -8.0 / l13;

  // YC
  twobody_YC (n, l, f);

  lub [8] = 0.0;
  lub [9] = 0.0;

  g2 = 0.4 * l / l1;
  g3 = 0.008 * (8.0 + l * (6.0 + l * 33.0)) / l1;
  g4 = 0.8 * l * l / l14;
  g5 = 0.016 * l * (43.0 + l * (-24.0 + l * 43.0)) / l14;


  lub [8] = - g2 * ls1
    - g3 * s1 * ls1
    + f[0];
  // next is m = 2

  lub [9] = g4 * ls2
    + g5 * s1 * ls2
    + 4.0 * g5 / s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - 2.0 * g2 / (double)m;
	  if (m == 2)
	    {
	      a += 4.0 * g3 / (double)(-2*m);
	    }
	  else
	    {
	      a += 4.0 * g3 / (double)(m*(m-2));
	    }
	  lub [8] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g4 / (double)m
	    + 4.0 * g5 / (double)(m*(m-2));
	  lub [9] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  //lub [9] *= 8.0 / l1 / l1 / l1;
  lub [9] *= 8.0 / l13;

  // FT version
  if (version == 1) return;


  // redefine some parameters for G, H, M parts
  s1 = 2.0 * s / (s * s - 4.0);
  //s2 = (s + 2.0) / (s - 2.0); // this is the same before
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  //ls2 = log (s2); // this is the same before
  ls5 = log (s5);
  

  // XG
  twobody_XG (n, l, f);

  lub[10] = 0.0;
  lub[11] = 0.0;

  g1 = 3.0 * l * l / l13;
  g2 = 0.3 * l * (1.0 + l * (12.0 - 4.0 * l)) / l13;
  g3 = (5.0 + l*(181.0 + l*(-453.0 + l*(566.0 - 65.0 *l)))) / 140.0 / l13;

  lub[10] = g1*s1 + (g2 + g3*s3)*ls2 - g3*s;
  // next is m = 1

  lub[11] = -g1*s4 - (g2 + g3*s3)*ls5 + g3;
  // next is m = 2

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[11] -= f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[10] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  //lub[11] *= 4.0 / l1 / l1;
  lub[11] *= 4.0 / l12;

  // YG
  twobody_YG (n, l, f);

  lub[12] = 0.0;
  lub[13] = 0.0;

  g2 = 0.1 * l * (4.0 + l * (-1.0 + 7.0 * l)) / l13;
  g3 = (32.0 + l*(-179.0 + l*(532.0 + l*(-356.0 + 221.0 *l)))) / 500.0 / l13;

  lub[12] = (g2 + g3*s3)*ls2 - g3*s;
  // next is m = 1

  lub[13] = - (g2 + g3*s3)*ls5 + g3;
  // next is m = 2

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[13] -= f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[12] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  //lub[13] *= 4.0 / l1 / l1;
  lub[13] *= 4.0 / l12;

  // YH
  twobody_YH (n, l, f);

  lub[14] = 0.0;
  lub[15] = 0.0;

  g2 = 0.1 * l * (2.0 - l) / l12;
  g3 = (16.0 + l*(-61.0 + l*(180.0 + 2.0 *l))) / 500.0 / l12;
  g5 = 0.05 * l * l * (1.0 + 7.0 * l) / l12;
  g6 = l * (43.0 + l*(147.0 + l*(-185.0 + 221.0 *l))) / 1000.0 / l12;

  lub[14] = (g2 + g3*s3)*ls5 - g3;
  // next is m = 2

  lub[15] = (g5 + g6*s3)*ls2 - g6*s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[14] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  lub[15] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  //lub[15] *= 8.0 / l12 / l1;
  lub[15] *= 8.0 / l13;

  // XM
  twobody_XM (n, l, f);

  lub[16] = 0.0;
  lub[17] = 0.0;

  g1 = 1.2 * l * l / l13;
  g2 = 0.12 * l * (1.0 + l*(17.0 -9.0 * l)) / l13;
  g3 = (5.0 + l*(272.0 + l*(-831.0 + l*(1322.0 - 415.0 *l)))) / 350.0 / l13;
  g4 = g1 * l;
  g5 = 0.12 * l * l * (-4.0 + l*(17.0 -4.0 * l)) / l13;
  g6 = l*(-65.0 + l*(832.0 + l*(-1041.0 + l*(832.0 - 65.0 *l)))) / 350.0 / l13;

  lub[16] = g1*s4 + (g2 + g3*s3)*ls5 - g3 + 1.0;
  // next is m = 2

  lub[17] = g4*s1 + (g5 + g6*s3)*ls2 - g6*s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[16] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - g4
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  lub[17] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub[17] *= 8.0 / l13;

  // YM
  twobody_YM (n, l, f);

  lub[18] = 0.0;
  lub[19] = 0.0;

  g2 = 0.24 * l * (1.0 + l*(-1.0 +4.0 * l)) / l13;
  g3 = (24.0 + l*(-201.0 + l*(882.0 + l*(-1182 + 591.0 *l)))) / 625.0 / l13;
  g5 = 0.06 * l * l * (7.0 + l*(-10.0 +7.0 * l)) / l13;
  g6 = 3.0*l*(221.0 + l*(-728.0 + l*(1902.0 + l*(-728 + 221.0 *l))))
    / 2500.0 / l13;

  lub[18] = (g2 + g3*s3)*ls5 - g3 + 1.0;
  // next is m = 2

  lub[19] = (g5 + g6*s3)*ls2 - g6*s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[18] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  lub[19] += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub[19] *= 8.0 / l13;

  // ZM
  twobody_ZM (n, l, f);

  lub[20] = 0.0;
  lub[21] = 0.0;

  g3 = -0.3*l*l*(1.0 + l*l) / l13;

  lub[20] = g3*s3*ls5 - g3 + 1.0;
  // next is m = 2

  lub[21] = - g3*s3*ls2 + g3*s;
  // next is m = 1

  s1l = s * (1.0 + l);
  s1lm = s1l;
  s21  = 2.0 / s;
  s2m  = s21;
  for (m = 1; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a =
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[20] += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    + 4.0 * g3 / (double)(m*(m+2));
	  lub[21] -= f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  lub[21] *= 8.0 / l13;


  free (f);
}

/* check twobody_lub()
 */
int
check_twobody_lub (int version, double r, double a1, double a2, int nmax,
		   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_lub : start\n");
    }

  int check = 0;

  double l = a2 / a1;
  double s = 2.0 * r / (a1 + a2);

  double lub1[22];
  double t0, t;
  t0 = ptime_ms_d ();
  twobody_lub (nmax, l, s, lub1);
  t = ptime_ms_d ();
  double ptime_lub1 = t - t0;

  double lub2[22];
  struct twobody_f *f = twobody_f_init (nmax, l);
  t0 = ptime_ms_d ();
  /*
  twobody_lub_ (2, // FTS version
		nmax, l, s, lub2);
  */
  twobody_lub_with_f (f, // struct twobody_f
		      2, // FTS version
		      nmax, l, s, lub2);
  t = ptime_ms_d ();
  double ptime_lub2 = t - t0;
  twobody_f_free (f);
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare (lub1 [0], lub2 [0], "check_twobody_lub: XA11", verbose, tiny);
  check += compare (lub1 [1], lub2 [1], "check_twobody_lub: XA12", verbose, tiny);
  check += compare (lub1 [2], lub2 [2], "check_twobody_lub: YA11", verbose, tiny);
  check += compare (lub1 [3], lub2 [3], "check_twobody_lub: YA12", verbose, tiny);
  check += compare (lub1 [4], lub2 [4], "check_twobody_lub: YB11", verbose, tiny);
  check += compare (lub1 [5], lub2 [5], "check_twobody_lub: YB12", verbose, tiny);
  check += compare (lub1 [6], lub2 [6], "check_twobody_lub: XC11", verbose, tiny);
  check += compare (lub1 [7], lub2 [7], "check_twobody_lub: XC12", verbose, tiny);
  check += compare (lub1 [8], lub2 [8], "check_twobody_lub: YC11", verbose, tiny);
  check += compare (lub1 [9], lub2 [9], "check_twobody_lub: YC12", verbose, tiny);
  check += compare (lub1[10], lub2[10], "check_twobody_lub: XG11", verbose, tiny);
  check += compare (lub1[11], lub2[11], "check_twobody_lub: XG12", verbose, tiny);
  check += compare (lub1[12], lub2[12], "check_twobody_lub: YG11", verbose, tiny);
  check += compare (lub1[13], lub2[13], "check_twobody_lub: YG12", verbose, tiny);
  check += compare (lub1[14], lub2[14], "check_twobody_lub: YH11", verbose, tiny);
  check += compare (lub1[15], lub2[15], "check_twobody_lub: YH12", verbose, tiny);
  check += compare (lub1[16], lub2[16], "check_twobody_lub: XM11", verbose, tiny);
  check += compare (lub1[17], lub2[17], "check_twobody_lub: XM12", verbose, tiny);
  check += compare (lub1[18], lub2[18], "check_twobody_lub: YM11", verbose, tiny);
  check += compare (lub1[19], lub2[19], "check_twobody_lub: YM12", verbose, tiny);
  check += compare (lub1[20], lub2[20], "check_twobody_lub: ZM11", verbose, tiny);
  check += compare (lub1[21], lub2[21], "check_twobody_lub: ZM12", verbose, tiny);

  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_lub :"
	       " ptime lub1, lub2 = %.3f %.3f, lub1/lub2 = %f\n",
	       ptime_lub1, ptime_lub2, ptime_lub1 / ptime_lub2);

      if (check == 0)
	fprintf (stdout, "check_twobody_lub :"
		 " PASSED\n\n");
      else
	fprintf (stdout, "check_twobody_lub :"
		 " FAILED\n\n");
    }

  return (check);
}



/* check scalars_minv_f_poly() with scalar_minv_f() for equal sphere
 */
int
check_twobody_scalars_res_with_equal (double r, int nmax,
				      int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_scalars_res_with_equal : start\n");
    }

  int check = 0;

  double mono[22];
  double t0, t;
  t0 = ptime_ms_d ();
  scalar_two_body_res (r, mono);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;

  double poly[22];
  t0 = ptime_ms_d ();
  twobody_scalars_res (r, 1.0, 1.0, nmax,
		       1, /* lub form */
		       1, /* dimensional form */
		       poly);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare (mono [0], poly [0], "check_twobody_scalars_res_with_equal: XA11", verbose, tiny);
  check += compare (mono [1], poly [1], "check_twobody_scalars_res_with_equal: XA12", verbose, tiny);
  check += compare (mono [2], poly [2], "check_twobody_scalars_res_with_equal: YA11", verbose, tiny);
  check += compare (mono [3], poly [3], "check_twobody_scalars_res_with_equal: YA12", verbose, tiny);
  check += compare (mono [4], poly [4], "check_twobody_scalars_res_with_equal: YB11", verbose, tiny);
  check += compare (mono [5], poly [5], "check_twobody_scalars_res_with_equal: YB12", verbose, tiny);
  check += compare (mono [6], poly [6], "check_twobody_scalars_res_with_equal: XC11", verbose, tiny);
  check += compare (mono [7], poly [7], "check_twobody_scalars_res_with_equal: XC12", verbose, tiny);
  check += compare (mono [8], poly [8], "check_twobody_scalars_res_with_equal: YC11", verbose, tiny);
  check += compare (mono [9], poly [9], "check_twobody_scalars_res_with_equal: YC12", verbose, tiny);
  check += compare (mono[10], poly[10], "check_twobody_scalars_res_with_equal: XG11", verbose, tiny);
  check += compare (mono[11], poly[11], "check_twobody_scalars_res_with_equal: XG12", verbose, tiny);
  check += compare (mono[12], poly[12], "check_twobody_scalars_res_with_equal: YG11", verbose, tiny);
  check += compare (mono[13], poly[13], "check_twobody_scalars_res_with_equal: YG12", verbose, tiny);
  check += compare (mono[14], poly[14], "check_twobody_scalars_res_with_equal: YH11", verbose, tiny);
  check += compare (mono[15], poly[15], "check_twobody_scalars_res_with_equal: YH12", verbose, tiny);
  check += compare (mono[16], poly[16], "check_twobody_scalars_res_with_equal: XM11", verbose, tiny);
  check += compare (mono[17], poly[17], "check_twobody_scalars_res_with_equal: XM12", verbose, tiny);
  check += compare (mono[18], poly[18], "check_twobody_scalars_res_with_equal: YM11", verbose, tiny);
  check += compare (mono[19], poly[19], "check_twobody_scalars_res_with_equal: YM12", verbose, tiny);
  check += compare (mono[20], poly[20], "check_twobody_scalars_res_with_equal: ZM11", verbose, tiny);
  check += compare (mono[21], poly[21], "check_twobody_scalars_res_with_equal: ZM12", verbose, tiny);

  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_scalars_res_with_equal :"
	       " ptime mono, poly = %.3f %.3f, poly/mono = %f\n",
	       ptime_mono, ptime_poly, ptime_poly / ptime_mono);

      if (check == 0)
	fprintf (stdout, "check_twobody_scalars_res_with_equal :"
		 " PASSED\n\n");
      else
	fprintf (stdout, "check_twobody_scalars_res_with_equal :"
		 " FAILED\n\n");
    }

  return (check);
}
