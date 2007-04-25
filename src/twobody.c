/* RYUON-twobody : exact 2-body resistance scalar functions
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: twobody.c,v 1.6 2007/04/25 05:34:59 kichiki Exp $
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


/** utility routines for struct twobody_f and twobody_f_list **/

struct twobody_f *
twobody_f_init (int nmax, double lambda)
{
  struct twobody_f *f
    = (struct twobody_f *)malloc (sizeof (struct twobody_f));
  CHECK_MALLOC (f, "twobody_f_init");

  f->nmax = nmax;
  f->lambda = lambda;

  f->XA = (double *)malloc (sizeof (double) * nmax);
  f->YA = (double *)malloc (sizeof (double) * nmax);
  f->YB = (double *)malloc (sizeof (double) * nmax);
  f->XC = (double *)malloc (sizeof (double) * nmax);
  f->YC = (double *)malloc (sizeof (double) * nmax);
  f->XG = (double *)malloc (sizeof (double) * nmax);
  f->YG = (double *)malloc (sizeof (double) * nmax);
  f->YH = (double *)malloc (sizeof (double) * nmax);
  f->XM = (double *)malloc (sizeof (double) * nmax);
  f->YM = (double *)malloc (sizeof (double) * nmax);
  f->ZM = (double *)malloc (sizeof (double) * nmax);

  CHECK_MALLOC (f->XA, "twobody_f_init");
  CHECK_MALLOC (f->YA, "twobody_f_init");
  CHECK_MALLOC (f->YB, "twobody_f_init");
  CHECK_MALLOC (f->XC, "twobody_f_init");
  CHECK_MALLOC (f->YC, "twobody_f_init");
  CHECK_MALLOC (f->XG, "twobody_f_init");
  CHECK_MALLOC (f->YG, "twobody_f_init");
  CHECK_MALLOC (f->YH, "twobody_f_init");
  CHECK_MALLOC (f->XM, "twobody_f_init");
  CHECK_MALLOC (f->YM, "twobody_f_init");
  CHECK_MALLOC (f->ZM, "twobody_f_init");

  twobody_XA (nmax, lambda, f->XA);
  twobody_YA (nmax, lambda, f->YA);
  twobody_YB (nmax, lambda, f->YB);
  twobody_XC (nmax, lambda, f->XC);
  twobody_YC (nmax, lambda, f->YC);
  twobody_XG (nmax, lambda, f->XG);
  twobody_YG (nmax, lambda, f->YG);
  twobody_YH (nmax, lambda, f->YH);
  twobody_XM (nmax, lambda, f->XM);
  twobody_YM (nmax, lambda, f->YM);
  twobody_ZM (nmax, lambda, f->ZM);

  return (f);
}

void
twobody_f_free (struct twobody_f *f)
{
  if (f == NULL) return;

  if (f->XA != NULL) free (f->XA);
  if (f->YA != NULL) free (f->YA);
  if (f->YB != NULL) free (f->YB);
  if (f->XC != NULL) free (f->XC);
  if (f->YC != NULL) free (f->YC);
  if (f->XG != NULL) free (f->XG);
  if (f->YG != NULL) free (f->YG);
  if (f->YH != NULL) free (f->YH);
  if (f->XM != NULL) free (f->XM);
  if (f->YM != NULL) free (f->YM);
  if (f->ZM != NULL) free (f->ZM);
  free (f);
}


struct twobody_f_list *
twobody_f_list_init (void)
{
  struct twobody_f_list *list
    = (struct twobody_f_list *)malloc (sizeof (struct twobody_f_list));
  CHECK_MALLOC (list, "twobody_f_list_init");

  list->n = 0;
  list->l = NULL;
  list->f = NULL;

  return (list);
}

void
twobody_f_list_append (struct twobody_f_list *list,
		       int nmax, double lambda)
{
  list->n++;

  list->l = (double *)realloc (list->l, sizeof (double) * (list->n));
  list->f = (struct twobody_f **)realloc (list->f,
					 sizeof (struct twobody_f *)
					 * (list->n));
  CHECK_MALLOC (list->l, "twobody_f_list_append");
  CHECK_MALLOC (list->f, "twobody_f_list_append");

  // i is the last element in the list
  int i = list->n - 1;
  list->l [i] = lambda;
  list->f [i] = twobody_f_init (nmax, lambda);
}

void
twobody_f_list_free (struct twobody_f_list *list)
{
  if (list == NULL) return;
  if (list->l != NULL) free (list->l);
  if (list->f != NULL)
    {
      int i;
      for (i = 0; i < list->n; i ++)
	{
	  if (list->f[i] != NULL) twobody_f_free (list->f[i]);
	}
      free (list->f);
    }
  free (list);
}


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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  double *f = (double *)malloc (sizeof (double) * n);
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

  *XA12 *= -2.0 / l1;

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

  double *f = (double *)malloc (sizeof (double) * n);
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

  *YA12 *= -2.0 / l1;

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

  double *f = (double *)malloc (sizeof (double) * n);
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

  *YB12 *= -4.0 / l12;

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

  double *f = (double *)malloc (sizeof (double) * n);
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

  *XC12 *= -8.0 / l1 / l1 / l1;

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

  double *f = (double *)malloc (sizeof (double) * n);
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

  *YC12 *= 8.0 / l1 / l1 / l1;

  free (f);
}

/* calc XG11 and XG12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_XG_lub (int n, double l, double s,
		     double *XG11, double *XG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XG_lub");

  twobody_XG (n, l, f);

  *XG11 = 0.0;
  *XG12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g1, g2, g3;
  g1 = 3.0 * l * l / l13;
  g2 = 0.3 * l * (1.0 + l * (12.0 - 4.0 * l)) / l13;
  g3 = (5.0 + l*(181.0 + l*(-453.0 + l*(566.0 - 65.0 *l)))) / 140.0 / l13;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *XG11 = g1*s1 + (g2 + g3*s3)*ls2 - g3*s;
  // next is m = 1

  *XG12 = -g1*s4 - (g2 + g3*s3)*ls5 + g3;
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
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *XG12 -= f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *XG11 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *XG12 *= 4.0 / l1 / l1;

  free (f);
}

/* calc YG11 and YG12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YG_lub (int n, double l, double s,
		     double *YG11, double *YG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YG_lub");

  twobody_YG (n, l, f);

  *YG11 = 0.0;
  *YG12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g2, g3;
  g2 = 0.1 * l * (4.0 + l * (-1.0 + 7.0 * l)) / l13;
  g3 = (32.0 + l*(-179.0 + l*(532.0 + l*(-356.0 + 221.0 *l)))) / 500.0 / l13;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *YG11 = (g2 + g3*s3)*ls2 - g3*s;
  // next is m = 1

  *YG12 = - (g2 + g3*s3)*ls5 + g3;
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
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *YG12 -= f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *YG11 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YG12 *= 4.0 / l1 / l1;

  free (f);
}

/* calc YH11 and YH12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YH_lub (int n, double l, double s,
		     double *YH11, double *YH12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YH_lub");

  twobody_YH (n, l, f);

  *YH11 = 0.0;
  *YH12 = 0.0;

  double l1 = 1.0 + l;
  double l12 = l1 * l1;
  double g2, g3, g5, g6;
  g2 = 0.1 * l * (2.0 - l) / l12;
  g3 = (16.0 + l*(-61.0 + l*(180.0 + 2.0 *l))) / 500.0 / l12;
  g5 = 0.05 * l * l * (1.0 + 7.0 * l) / l12;
  g6 = l * (43.0 + l*(147.0 + l*(-185.0 + 221.0 *l))) / 1000.0 / l12;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *YH11 = (g2 + g3*s3)*ls5 - g3;
  // next is m = 2

  *YH12 = (g5 + g6*s3)*ls2 - g6*s;
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
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *YH11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  *YH12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YH12 *= 8.0 / l12 / l1;

  free (f);
}

/* calc XM11 and XM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_XM_lub (int n, double l, double s,
		     double *XM11, double *XM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_XM_lub");

  twobody_XM (n, l, f);

  *XM11 = 0.0;
  *XM12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g1, g2, g3, g4, g5, g6;
  g1 = 1.2 * l * l / l13;
  g2 = 0.12 * l * (1.0 + l*(17.0 -9.0 * l)) / l13;
  g3 = (5.0 + l*(272.0 + l*(-831.0 + l*(1322.0 - 415.0 *l)))) / 350.0 / l13;
  g4 = g1 * l;
  g5 = 0.12 * l * l * (-4.0 + l*(17.0 -4.0 * l)) / l13;
  g6 = l*(-65.0 + l*(832.0 + l*(-1041.0 + l*(832.0 - 65.0 *l)))) / 350.0 / l13;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *XM11 = g1*s4 + (g2 + g3*s3)*ls5 - g3 + 1.0;
  // next is m = 2

  *XM12 = g4*s1 + (g5 + g6*s3)*ls2 - g6*s;
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
	  a = - g1
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *XM11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a = - g4
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  *XM12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *XM12 *= 8.0 / l13;

  free (f);
}

/* calc YM11 and YM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YM_lub (int n, double l, double s,
		     double *YM11, double *YM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_YM_lub");

  twobody_YM (n, l, f);

  *YM11 = 0.0;
  *YM12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g2, g3, g5, g6;
  g2 = 0.24 * l * (1.0 + l*(-1.0 +4.0 * l)) / l13;
  g3 = (24.0 + l*(-201.0 + l*(882.0 + l*(-1182 + 591.0 *l)))) / 625.0 / l13;
  g5 = 0.06 * l * l * (7.0 + l*(-10.0 +7.0 * l)) / l13;
  g6 = 3.0*l*(221.0 + l*(-728.0 + l*(1902.0 + l*(-728 + 221.0 *l))))
    / 2500.0 / l13;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *YM11 = (g2 + g3*s3)*ls5 - g3 + 1.0;
  // next is m = 2

  *YM12 = (g5 + g6*s3)*ls2 - g6*s;
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
	  a =
	    - 2.0 * g2 / (double)m
	    + 4.0 * g3 / (double)(m*(m+2));
	  *YM11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    - 2.0 * g5 / (double)m
	    + 4.0 * g6 / (double)(m*(m+2));
	  *YM12 += f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *YM12 *= 8.0 / l13;

  free (f);
}

/* calc ZM11 and ZM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_ZM_lub (int n, double l, double s,
		     double *ZM11, double *ZM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (f, "twobody_ZM_lub");

  twobody_ZM (n, l, f);

  *ZM11 = 0.0;
  *ZM12 = 0.0;

  double l1 = 1.0 + l;
  double l13 = l1 * l1 * l1;
  double g3;
  g3 = -0.3*l*l*(1.0 + l*l) / l13;

  double s1, s2, s3, s4, s5;
  double ls2, ls5;
  s1 = 2.0 * s / (s * s - 4.0);
  s2 = (s + 2.0) / (s - 2.0);
  s3 = (s * s - 4.0) / 4.0;
  s4 = 4.0 / (s * s - 4.0);
  s5 = s * s / (s * s - 4.0);
  ls2 = log (s2);
  ls5 = log (s5);
  
  *ZM11 = g3*s3*ls5 - g3 + 1.0;
  // next is m = 2

  *ZM12 = - g3*s3*ls2 + g3*s;
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
	  a =
	    + 4.0 * g3 / (double)(m*(m+2));
	  *ZM11 += f[m] / s1lm
	    + a * s2m;
	}
      else
	{
	  a =
	    + 4.0 * g3 / (double)(m*(m+2));
	  *ZM12 -= f[m] / s1lm
	    + a * s2m;
	}

      s1lm *= s1l;
      s2m  *= s21;
    }

  *ZM12 *= 8.0 / l13;

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
 *     10,11 : (XG11, XG12)
 *     12,13 : (YG11, YG12)
 *     14,15 : (YH11, YH12)
 *     16,17 : (XM11, XM12)
 *     18,19 : (YM11, YM12)
 *     20,21 : (ZM11, ZM12)
 */
void twobody_lub (int n, double l, double s,
		  double *lub)
{
  twobody_XA_lub (n, l, s, lub +  0, lub +  1);
  twobody_YA_lub (n, l, s, lub +  2, lub +  3);
  twobody_YB_lub (n, l, s, lub +  4, lub +  5);
  twobody_XC_lub (n, l, s, lub +  6, lub +  7);
  twobody_YC_lub (n, l, s, lub +  8, lub +  9);
  twobody_XG_lub (n, l, s, lub + 10, lub + 11);
  twobody_YG_lub (n, l, s, lub + 12, lub + 13);
  twobody_YH_lub (n, l, s, lub + 14, lub + 15);
  twobody_XM_lub (n, l, s, lub + 16, lub + 17);
  twobody_YM_lub (n, l, s, lub + 18, lub + 19);
  twobody_ZM_lub (n, l, s, lub + 20, lub + 21);
}

/* calc scalar functions of resistance problem by lub form
 * all-in-one form (to reduce calculating the same parameters)
 * INPUT
 *  f2b     : struct twobody_f for the pair
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
void twobody_lub_with_f (struct twobody_f *f2b,
			 int version,
			 int n, double l, double s,
			 double *lub)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = NULL;
  if (f2b == NULL)
    {
      f = (double *)malloc (sizeof (double) * n);
      CHECK_MALLOC (f, "twobody_XA_lub");
    }

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
  if (f2b == NULL)
    {
      twobody_XA (n, l, f);
    }
  else
    {
      f = f2b->XA;
    }

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
  if (f2b == NULL)
    {
      twobody_YA (n, l, f);
    }
  else
    {
      f = f2b->YA;
    }

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
  if (f2b == NULL)
    {
      twobody_YB (n, l, f);
    }
  else
    {
      f = f2b->YB;
    }

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
  if (f2b == NULL)
    {
      twobody_XC (n, l, f);
    }
  else
    {
      f = f2b->XC;
    }

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
  if (f2b == NULL)
    {
      twobody_YC (n, l, f);
    }
  else
    {
      f = f2b->YC;
    }

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
  if (f2b == NULL)
    {
      twobody_XG (n, l, f);
    }
  else
    {
      f = f2b->XG;
    }

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
  if (f2b == NULL)
    {
      twobody_YG (n, l, f);
    }
  else
    {
      f = f2b->YG;
    }

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
  if (f2b == NULL)
    {
      twobody_YH (n, l, f);
    }
  else
    {
      f = f2b->YH;
    }

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
  if (f2b == NULL)
    {
      twobody_XM (n, l, f);
    }
  else
    {
      f = f2b->XM;
    }

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
  if (f2b == NULL)
    {
      twobody_YM (n, l, f);
    }
  else
    {
      f = f2b->YM;
    }

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
  if (f2b == NULL)
    {
      twobody_ZM (n, l, f);
    }
  else
    {
      f = f2b->ZM;
    }

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


  if (f2b == NULL)
    {
      free (f);
    }
}


/* scale the scalar functions from Jeffrey-Onishi to Stokesian dynamics
 * INPUT
 *  two [22] : scalar functions in Jeffrey form
 *  l        : lambda = ab / aa,
 *             where aa and ab are radii for particles a(alpha) and b(beta)
 *             Note that the scalar functions are for "a-b" interaction.
 * OUTPUT
 *  two [22] : scalar functions in SD form
 */
void
twobody_scale_SD (double *two, double l)
{
  double l1 = 1.0 + l;
  double l2 = l1 * l1;
  double l3 = l2 * l1;

  l1 /= 2.0; // l1 = (1 + lambda)   / 2 (= 1   for lambda == 1)
  l2 /= 6.0; // l2 = (1 + lambda)^2 / 6 (= 2/3 for lambda == 1)
  l3 /= 6.0; // l3 = (1 + lambda)^3 / 6 (= 4/3 for lambda == 1)
  double lm;
  lm = 5.0 / 6.0 * l3; // l3 = (5/36)(1+lambda)^3 (= 10/9 for lambda == 1)

  // A part
  two [1] *= l1;// XA12
  two [3] *= l1;// YA12

  // B part
  two [4] *= 2.0 / 3.0; // YB11
  two [5] *= l2;        // YB12

  // C part
  two [6] *= 4.0 / 3.0; // XC11
  two [7] *= l3;        // XC12
  two [8] *= 4.0 / 3.0; // YC11
  two [9] *= l3;        // YC12

  // G part
  two [10] *= 2.0 / 3.0; // XG11
  two [11] *= l2;        // XG12
  two [12] *= 2.0 / 3.0; // YG11
  two [13] *= l2;        // YG12

  // H part
  two [14] *= 4.0 / 3.0; // YH11
  two [15] *= l3;        // YH12

  // M part
  two [16] *= 10.0 / 9.0; // XM11
  two [17] *= lm;         // XM12
  two [18] *= 10.0 / 9.0; // YM11
  two [19] *= lm;         // YM12
  two [20] *= 10.0 / 9.0; // ZM11
  two [21] *= lm;         // ZM12
}

/* scale the scalar functions from Jeffrey-Onishi to the dimensional form
 * INPUT
 *  two [22] : scalar functions in Jeffrey form
 *  l        : lambda = ab / aa,
 *             where aa and ab are radii for particles a(alpha) and b(beta)
 *             Note that the scalar functions are for "a-b" interaction.
 * OUTPUT
 *  two [22] : scalar functions in SD form
 */
void
twobody_scale (double *two, double a1, double l)
{
  double a12 = a1  * a1;
  double a13 = a12 * a1;

  double l1 = 1.0 + l;
  double l2 = l1 * l1;
  double l3 = l2 * l1;

  l1 /= 2.0; // l1 = (1 + lambda)   / 2 (= 1   for lambda == 1)
  l2 /= 6.0; // l2 = (1 + lambda)^2 / 6 (= 2/3 for lambda == 1)
  l3 /= 6.0; // l3 = (1 + lambda)^3 / 6 (= 4/3 for lambda == 1)
  double lm;
  lm = 5.0 / 6.0 * l3; // l3 = (5/36)(1+lambda)^3 (= 10/9 for lambda == 1)

  // A part
  two [0] *= a1;     // XA11
  two [1] *= a1 * l1;// XA12
  two [2] *= a1;     // YA11
  two [3] *= a1 * l1;// YA12

  // B part
  two [4] *= a12 * 2.0 / 3.0; // YB11
  two [5] *= a12 * l2;        // YB12

  // C part
  two [6] *= a13 * 4.0 / 3.0; // XC11
  two [7] *= a13 * l3;        // XC12
  two [8] *= a13 * 4.0 / 3.0; // YC11
  two [9] *= a13 * l3;        // YC12

  // G part
  two [10] *= a12 * 2.0 / 3.0; // XG11
  two [11] *= a12 * l2;        // XG12
  two [12] *= a12 * 2.0 / 3.0; // YG11
  two [13] *= a12 * l2;        // YG12

  // H part
  two [14] *= a13 * 4.0 / 3.0; // YH11
  two [15] *= a13 * l3;        // YH12

  // M part
  two [16] *= a13 * 10.0 / 9.0; // XM11
  two [17] *= a13 * lm;         // XM12
  two [18] *= a13 * 10.0 / 9.0; // YM11
  two [19] *= a13 * lm;         // YM12
  two [20] *= a13 * 10.0 / 9.0; // ZM11
  two [21] *= a13 * lm;         // ZM12
}


/* calc scalar functions of two-body exact solution in resistance problem
 * INPUT
 *  r          : distance between the two := x_b - x_a
 *  aa, ab     : radii for particles a(alpha) and b(beta)
 *  n          : max order for the coefficients
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 *  flag_scale : 0 no scaling, that is, in Jeffrey form
 *               1 for the dimensional form
 *               2 for the Stokesian dynamics form
 *  res [22]   :
 * OUTPUT
 *  res [22]   : scalar functions. the scaling is given by flag_scale.
 */
void
twobody_scalars_res (double r,
		     double aa, double ab,
		     int n, int flag_lub, int flag_scale,
		     double *res)
{
  double l = ab / aa;
  double s = 2.0 * r / (aa + ab);

  if (flag_lub == 0)
    {
      twobody_far (n, l, s, res);
    }
  else
    {
      twobody_lub (n, l, s, res);
    }

  if (flag_scale == 1)
    {
      twobody_scale (res, aa, l);
    }
  else if (flag_scale == 2)
    {
      twobody_scale_SD (res, l);      
    }
}
