/* RYUON-twobody : exact 2-body resistance scalar functions
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: twobody.c,v 1.2 2007/03/27 07:02:29 kichiki Exp $
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
#include <math.h> // log()
#include "memory-check.h" // macro CHECK_MALLOC

#include "twobody.h"


/** far form **/

/* calc XA11 and XA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_far (int n, double l, double s,
		     double *XA11, double *XA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XA_far");

  twobody_XA (n, l, f);

  *XA11 = 0.0;
  *XA12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *XA11 += f[i*2  ] / s11;
      *XA12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *XA12 *= -2.0 / (1.0 + l);

  free (f);
}

/* calc YA11 and YA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_far (int n, double l, double s,
		     double *YA11, double *YA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YA_far");

  twobody_YA (n, l, f);

  *YA11 = 0.0;
  *YA12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YA11 += f[i*2  ] / s11;
      *YA12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YA12 *= -2.0 / (1.0 + l);

  free (f);
}

/* calc YB11 and YB12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_far (int n, double l, double s,
		     double *YB11, double *YB12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YB_far");

  twobody_YB (n, l, f);

  *YB11 = 0.0;
  *YB12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = s1l; // odd
  s12 = 1.0; // even
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YB11 += f[i*2+1] / s11;
      *YB12 += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YB12 *= -4.0 / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc XC11 and XC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_far (int n, double l, double s,
		     double *XC11, double *XC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XC_far");

  twobody_XC (n, l, f);

  *XC11 = 0.0;
  *XC12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *XC11 += f[i*2  ] / s11;
      *XC12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *XC12 *= -8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc YC11 and YC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_far (int n, double l, double s,
		     double *YC11, double *YC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YC_far");

  twobody_YC (n, l, f);

  *YC11 = 0.0;
  *YC12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YC11 += f[i*2  ] / s11;
      *YC12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YC12 *= 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc XG11 and XG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XG11
 *  XG12
 */
void twobody_XG_far (int n, double l, double s,
		     double *XG11, double *XG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XG_far");

  twobody_XG (n, l, f);

  *XG11 = 0.0;
  *XG12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = s1l; // odd
  s12 = 1.0; // even
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *XG11 += f[i*2+1] / s11;
      *XG12 += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *XG12 *= -4.0 / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc YG11 and YG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YG11
 *  YG12
 */
void twobody_YG_far (int n, double l, double s,
		     double *YG11, double *YG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YG_far");

  twobody_YG (n, l, f);

  *YG11 = 0.0;
  *YG12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = s1l; // odd
  s12 = 1.0; // even
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YG11 += f[i*2+1] / s11;
      *YG12 += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YG12 *= - 4.0 / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc YH11 and YH12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YH11
 *  YH12
 */
void twobody_YH_far (int n, double l, double s,
		     double *YH11, double *YH12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YH_far");

  twobody_YH (n, l, f);

  *YH11 = 0.0;
  *YH12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YH11 += f[i*2  ] / s11;
      *YH12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YH12 *= 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc XM11 and XM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XM11
 *  XM12
 */
void twobody_XM_far (int n, double l, double s,
		     double *XM11, double *XM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XM_far");

  twobody_XM (n, l, f);

  *XM11 = 0.0;
  *XM12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *XM11 += f[i*2  ] / s11;
      *XM12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *XM12 *= 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc YM11 and YM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YM11
 *  YM12
 */
void twobody_YM_far (int n, double l, double s,
		     double *YM11, double *YM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YM_far");

  twobody_YM (n, l, f);

  *YM11 = 0.0;
  *YM12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *YM11 += f[i*2  ] / s11;
      *YM12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *YM12 *= 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc ZM11 and ZM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  ZM11
 *  ZM12
 */
void twobody_ZM_far (int n, double l, double s,
		     double *ZM11, double *ZM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_ZM_far");

  twobody_ZM (n, l, f);

  *ZM11 = 0.0;
  *ZM12 = 0.0;

  double s1l = s * (1.0 + l);
  double s11, s12;
  s11 = 1.0; // even
  s12 = s1l; // odd
  double s2 = s1l * s1l;

  int i;
  for (i = 0; i < n/2; i ++)
    {
      *ZM11 += f[i*2  ] / s11;
      *ZM12 += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }

  *ZM12 *= - 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc scalar functions of resistance problem by 1/s expansion
 * INPUT
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
void twobody_far (int n, double l, double s,
		  double *far)
{
  twobody_XA_far (n, l, s, far +  0, far +  1);
  twobody_YA_far (n, l, s, far +  2, far +  3);
  twobody_YB_far (n, l, s, far +  4, far +  5);
  twobody_XC_far (n, l, s, far +  6, far +  7);
  twobody_YC_far (n, l, s, far +  8, far +  9);
  twobody_XG_far (n, l, s, far + 10, far + 11);
  twobody_YG_far (n, l, s, far + 12, far + 13);
  twobody_YH_far (n, l, s, far + 14, far + 15);
  twobody_XM_far (n, l, s, far + 16, far + 17);
  twobody_YM_far (n, l, s, far + 18, far + 19);
  twobody_ZM_far (n, l, s, far + 20, far + 21);
}


/** lubrication form **/

/* calc XA11 and XA12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_lub (int n, double l, double s,
		     double *XA11, double *XA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XA_lub");

  twobody_XA (n, l, f);

  *XA11 = 0.0;
  *XA12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g1, g2, g3;
  g1 = 2.0 * l * l / l13;
  g2 = 0.2 * l * (1.0 + l * (7.0 + l)) / l13;
  g3 = (1.0 + l * (18.0 + l * (- 29.0 + l * (18.0 + l)))) / l13 / 42.0;

  double s1, s2, s3, s4, s5, ls1, ls2;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  *XA11 = g1 / s1
    - g2 * ls1
    - g3 * s1 * ls1
    + f[0]
    - g1;
  // next is m = 2

  *XA12 = 2.0 / s * g1 / s1
    + g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  double s1l = s * (1.0 + l);
  double s1lm = s1l;

  double s21  = 2.0 / s;
  double s2m  = s21;

  double a;
  int m;
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
	  *XA11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = -g1 - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  *XA12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *XA12 *= -2.0 / (1.0 + l);

  free (f);
}

/* calc YA11 and YA12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_lub (int n, double l, double s,
		     double *YA11, double *YA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YA_lub");

  twobody_YA (n, l, f);

  *YA11 = 0.0;
  *YA12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g2, g3;
  g2 = (4.0 / 15.0) * l * (2.0 + l * (1.0 + l * 2.0)) / l13;
  g3 = (2.0 / 375.0)
    * (16.0 + l * (-45.0 + l * (58.0 + l * (-45.0 + l * 16.0)))) / l13;

  double s1, s2, s3, s4, s5, ls1, ls2;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  *YA11 = - g2 * ls1
    - g3 * s1 * ls1
    + f[0];
  // next is m = 2

  *YA12 = g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  double s1l = s * (1.0 + l);
  double s1lm = s1l;

  double s21  = 2.0 / s;
  double s2m  = s21;

  double a;
  int m;
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
	  *YA11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  *YA12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YA12 *= -2.0 / (1.0 + l);

  free (f);
}

/* calc YB11 and YB12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_lub (int n, double l, double s,
		     double *YB11, double *YB12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YB_lub");

  twobody_YB (n, l, f);

  *YB11 = 0.0;
  *YB12 = 0.0;

  double l1 = 1.0 + l;
  double l12 = l1 * l1;
  double g2, g3;
  g2 = -0.2 * l * (4.0 + l) / l12;
  g3 = -0.004 * (32.0 + l * (-33.0 + l * (83.0 + l * 43.0))) / l12;

  double s1, s2, s3, s4, s5, ls1, ls2;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  *YB11 = g2 * ls2
    + g3 * s1 * ls2
    + 4.0 * g3 / s;
  // next is m = 1

  *YB12 = - g2 * ls1
    - g3 * s1 * ls1;
  // next is m = 2

  double s1l = s * (1.0 + l);
  double s1lm = s1l;

  double s21  = 2.0 / s;
  double s2m  = s21;

  double a;
  int m;
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
	  *YB12 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m-2));
	  *YB11 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YB12 *= -4.0 / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc XC11 and XC12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_lub (int n, double l, double s,
		     double *XC11, double *XC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XC_lub");

  twobody_XC (n, l, f);

  *XC11 = 0.0;
  *XC12 = 0.0;

  double l1 = 1.0 + l;

  double s1, s2, s3, s4, s5, ls1, ls2;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  *XC11 = l * l / l1 * (0.5 * ls1 + ls2 / s)
    + 1.0;
  // next is m = 2

  *XC12 = - l * l / l1 * (0.5 * ls2 + (ls1 - 2.0) / s);
  // next is m = 3

  double s1l = s * (1.0 + l);
  double s1lm = s1l * s1l; // starting from m=2

  double s21  = 2.0 / s;
  double s2m  = s21 * s21; // starting from m=2

  double a;
  int m;
  for (m = 2; m < n; m ++)
    {
      if (m%2 == 0)
	{
	  a = - l * l / l1 / (double)(m * (m-1));
	  *XC11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - l * l / l1 / (double)(m * (m-1));
	  *XC12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *XC12 *= -8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}

/* calc YC11 and YC12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_lub (int n, double l, double s,
		     double *YC11, double *YC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YC_lub");

  twobody_YC (n, l, f);

  *YC11 = 0.0;
  *YC12 = 0.0;

  double l1 = 1.0 + l;
  double l14 = l1 * l1 * l1 * l1;
  double g2, g3, g4, g5;
  g2 = 0.4 * l / l1;
  g3 = 0.008 * (8.0 + l * (6.0 + l * 33.0)) / l1;
  g4 = 0.8 * l * l / l14;
  g5 = 0.016 * l * (43.0 + l * (-24.0 + l * 43.0)) / l14;


  double s1, s2, s3, s4, s5, ls1, ls2;
  s1 = 1.0 - 4.0 / s / s;
  s2 = (s + 2.0) / (s - 2.0);
  s3 = s * s;
  s4 =  1.0 / s3;
  s5 = s * s / 4.0 - 1.0;
  ls1 = log (s1);
  ls2 = log (s2);
  
  *YC11 = - g2 * ls1
    - g3 * s1 * ls1
    + f[0];
  // next is m = 2

  *YC12 = g4 * ls2
    + g5 * s1 * ls2
    + 4.0 * g5 / s;
  // next is m = 1

  double s1l = s * (1.0 + l);
  double s1lm = s1l;

  double s21  = 2.0 / s;
  double s2m  = s21;

  double a;
  int m;
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
	  *YC11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - 2.0 * g4 / (double)m
	    + 4.0 * g5 / (double)(m*(m-2));
	  *YC12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YC12 *= 8.0 / (1.0 + l) / (1.0 + l) / (1.0 + l);

  free (f);
}




/* calc scalar functions of resistance problem by lub form
 * INPUT
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
 *     10,11 : (XG11, XG12) not implemented
 *     12,13 : (YG11, YG12) not implemented
 *     14,15 : (YH11, YH12) not implemented
 *     16,17 : (XM11, XM12) not implemented
 *     18,19 : (YM11, YM12) not implemented
 *     20,21 : (ZM11, ZM12) not implemented
 */
void twobody_lub (int n, double l, double s,
		  double *far)
{
  twobody_XA_lub (n, l, s, far +  0, far +  1);
  twobody_YA_lub (n, l, s, far +  2, far +  3);
  twobody_YB_lub (n, l, s, far +  4, far +  5);
  twobody_XC_lub (n, l, s, far +  6, far +  7);
  twobody_YC_lub (n, l, s, far +  8, far +  9);
  /*
  twobody_XG_lub (n, l, s, far + 10, far + 11);
  twobody_YG_lub (n, l, s, far + 12, far + 13);
  twobody_YH_lub (n, l, s, far + 14, far + 15);
  twobody_XM_lub (n, l, s, far + 16, far + 17);
  twobody_YM_lub (n, l, s, far + 18, far + 19);
  twobody_ZM_lub (n, l, s, far + 20, far + 21);
  */
}
