/* twobody solutions for slip particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: twobody-slip.c,v 1.1 2007/08/17 04:31:19 kichiki Exp $
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
#include <math.h> // sqrt()
#include <stokes.h> // struct stokes
#include <twobody.h> // twobody_scale(), twobody_scale_SD()
#include <minv-poly.h> // scalars_minv_f_poly()
#include "memory-check.h" // macro CHECK_MALLOC

#include "twobody-slip.h"


/** global variable **/
double SL_hat_g1;
double SL_hat_g2;

/* 
 * INPUT
 *  (extern) lambda : slip length scaled by the radius
 *                    if negative, the slip length is infinity (perfect slip)
 */
double SL_G1 (int m, int n)
{
  extern double SL_hat_g1;

  double L;
  if (SL_hat_g1 >= 0.0)
    {
      L = (1.0 + (double)m * SL_hat_g1) / (1.0 + (double)n * SL_hat_g1);
    }
  else
    {
      L = (double)m / (double)n;
    }

  return (L);
}
double SL_G2 (int m, int n)
{
  extern double SL_hat_g2;

  double L;
  if (SL_hat_g2 >= 0.0)
    {
      L = (1.0 + (double)m * SL_hat_g2) / (1.0 + (double)n * SL_hat_g2);
    }
  else
    {
      L = (double)m / (double)n;
    }

  return (L);
}


/** utility routines for struct twobody_f and twobody_f_list **/

struct twobody_slip_f *
twobody_slip_f_init (int nmax,
		     double lambda,
		     double hat_g1, double hat_g2,
		     double slip_a1, double slip_a2)
{
  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;

  struct twobody_slip_f *f
    = (struct twobody_slip_f *)malloc (sizeof (struct twobody_slip_f));
  CHECK_MALLOC (f, "twobody_slip_f_init");

  f->nmax = nmax;
  f->lambda = lambda;
  f->hat_g1 = hat_g1;
  f->hat_g2 = hat_g2;
  f->slip_a1 = slip_a1;
  f->slip_a2 = slip_a2;

  f->XA = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YA = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YB = (double *)malloc (sizeof (double) * (nmax + 1));
  f->XC = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YC = (double *)malloc (sizeof (double) * (nmax + 1));
  f->XG = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YG = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YH = (double *)malloc (sizeof (double) * (nmax + 1));
  f->XM = (double *)malloc (sizeof (double) * (nmax + 1));
  f->YM = (double *)malloc (sizeof (double) * (nmax + 1));
  f->ZM = (double *)malloc (sizeof (double) * (nmax + 1));

  CHECK_MALLOC (f->XA, "twobody_slip_f_init");
  CHECK_MALLOC (f->YA, "twobody_slip_f_init");
  CHECK_MALLOC (f->YB, "twobody_slip_f_init");
  CHECK_MALLOC (f->XC, "twobody_slip_f_init");
  CHECK_MALLOC (f->YC, "twobody_slip_f_init");
  CHECK_MALLOC (f->XG, "twobody_slip_f_init");
  CHECK_MALLOC (f->YG, "twobody_slip_f_init");
  CHECK_MALLOC (f->YH, "twobody_slip_f_init");
  CHECK_MALLOC (f->XM, "twobody_slip_f_init");
  CHECK_MALLOC (f->YM, "twobody_slip_f_init");
  CHECK_MALLOC (f->ZM, "twobody_slip_f_init");

  twobody_XA_slip (nmax, lambda, f->XA);
  twobody_YA_slip (nmax, lambda, f->YA);
  twobody_YB_slip (nmax, lambda, f->YB);
  twobody_XC_slip (nmax, lambda, f->XC);
  twobody_YC_slip (nmax, lambda, f->YC);
  twobody_XG_slip (nmax, lambda, f->XG);
  twobody_YG_slip (nmax, lambda, f->YG);
  twobody_YH_slip (nmax, lambda, f->YH);
  twobody_XM_slip (nmax, lambda, f->XM);
  twobody_YM_slip (nmax, lambda, f->YM);
  twobody_ZM_slip (nmax, lambda, f->ZM);

  return (f);
}

void
twobody_slip_f_free (struct twobody_slip_f *f)
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


struct twobody_slip_f_list *
twobody_slip_f_list_init (void)
{
  struct twobody_slip_f_list *list
    = (struct twobody_slip_f_list *)malloc
    (sizeof (struct twobody_slip_f_list));
  CHECK_MALLOC (list, "twobody_slip_f_list_init");

  list->n = 0;
  list->f = NULL;

  return (list);
}

void
twobody_slip_f_list_append (struct twobody_slip_f_list *list,
			    int nmax,
			    double lambda,
			    double hat_g1, double hat_g2,
			    double slip_a1, double slip_a2)
{
  list->n++;
  list->f = (struct twobody_slip_f **)realloc
    (list->f,
     sizeof (struct twobody_slip_f *) * (list->n));
  CHECK_MALLOC (list->f, "twobody_slip_f_list_append");

  // i is the last element in the list
  int i = list->n - 1;
  list->f [i] = twobody_slip_f_init (nmax,
				     lambda, hat_g1, hat_g2, slip_a1, slip_a2);
}

void
twobody_slip_f_list_free (struct twobody_slip_f_list *list)
{
  if (list == NULL) return;
  if (list->f != NULL)
    {
      int i;
      for (i = 0; i < list->n; i ++)
	{
	  if (list->f[i] != NULL) twobody_slip_f_free (list->f[i]);
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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XA11, double *XA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_XA_slip_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_XA_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YA11, double *YA12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YA_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YA_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YB11, double *YB12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YB_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YB_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XC11, double *XC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_XC_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_XC_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YC11, double *YC12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YC_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YC_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XG11
 *  XG12
 */
void twobody_XG_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XG11, double *XG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_XG_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_XG_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YG11
 *  YG12
 */
void twobody_YG_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YG11, double *YG12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YG_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YG_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YH11
 *  YH12
 */
void twobody_YH_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YH11, double *YH12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YH_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YH_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XM11
 *  XM12
 */
void twobody_XM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XM11, double *XM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_XM_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_XM_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YM11
 *  YM12
 */
void twobody_YM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YM11, double *YM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_YM_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_YM_slip (n, l, f);

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
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  ZM11
 *  ZM12
 */
void twobody_ZM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *ZM11, double *ZM12)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  double *f = (double *)malloc (sizeof (double) * (n + 1));
  CHECK_MALLOC (f, "twobody_ZM_far");

  extern double SL_hat_g1;
  extern double SL_hat_g2;
  SL_hat_g1 = hat_g1;
  SL_hat_g2 = hat_g2;
  twobody_ZM_slip (n, l, f);

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
 *  version : 0=F, 1=FT, 2=FTS.
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
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
void twobody_slip_far (int version, int n, double l, double s,
		       double hat_g1, double hat_g2,
		       double *far)
{
  twobody_XA_slip_far (n, l, s, hat_g1, hat_g2, far +  0, far +  1);
  twobody_YA_slip_far (n, l, s, hat_g1, hat_g2, far +  2, far +  3);
  if (version == 0) return; // F version

  twobody_YB_slip_far (n, l, s, hat_g1, hat_g2, far +  4, far +  5);
  twobody_XC_slip_far (n, l, s, hat_g1, hat_g2, far +  6, far +  7);
  twobody_YC_slip_far (n, l, s, hat_g1, hat_g2, far +  8, far +  9);
  if (version == 1) return; // FT version

  twobody_XG_slip_far (n, l, s, hat_g1, hat_g2, far + 10, far + 11);
  twobody_YG_slip_far (n, l, s, hat_g1, hat_g2, far + 12, far + 13);
  twobody_YH_slip_far (n, l, s, hat_g1, hat_g2, far + 14, far + 15);
  twobody_XM_slip_far (n, l, s, hat_g1, hat_g2, far + 16, far + 17);
  twobody_YM_slip_far (n, l, s, hat_g1, hat_g2, far + 18, far + 19);
  twobody_ZM_slip_far (n, l, s, hat_g1, hat_g2, far + 20, far + 21);
}

/* calc scalar functions of resistance problem by 1/s expansion
 * all-in-one form (to reduce calculating the same parameters)
 * and with struct twobody_f *f12 table (to avoid recalculating them)
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  f12     : struct twobody_f for the pair
 *            you can give NULL for them.
 *            then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
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
void twobody_slip_far_with_f (int version,
			      struct twobody_slip_f *f12,
			      int n, double l, double s,
			      double hat_g1, double hat_g2,
			      double *far)
{
  if (n%2 != 0) n--; // make n even for fail-safe

  extern double SL_hat_g1;
  extern double SL_hat_g2;

  double *f = NULL;
  if (f12 == NULL)
    {
      SL_hat_g1 = hat_g1;
      SL_hat_g2 = hat_g2;

      f = (double *)malloc (sizeof (double) * (n + 1));
      CHECK_MALLOC (f, "twobody_slip_far_with_f");
    }

  double l1 = 1.0 + l;
  double l12 = l1 * l1;
  double l13 = l12 * l1;

  double s1l = s * l1;
  double s11, s12;
  double s2;

  int i;

  // XA
  if (f12 == NULL)
    {
      twobody_XA_slip (n, l, f);
    }
  else
    {
      f = f12->XA;
    }
  far [0] = 0.0;
  far [1] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far [0] += f[i*2  ] / s11;
      far [1] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far [1] *= -2.0 / l1;

  // YA
  if (f12 == NULL)
    {
      twobody_YA_slip (n, l, f);
    }
  else
    {
      f = f12->YA;
    }
  far [2] = 0.0;
  far [3] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far [2] += f[i*2  ] / s11;
      far [3] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far [3] *= -2.0 / l1;

  // F version
  if (version == 0) return;

  // YB
  if (f12 == NULL)
    {
      twobody_YB_slip (n, l, f);
    }
  else
    {
      f = f12->YB;
    }
  far [4] = 0.0;
  far [5] = 0.0;
  s11 = s1l; // odd
  s12 = 1.0; // even
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far [4] += f[i*2+1] / s11;
      far [5] += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far [5] *= -4.0 / l12;

  // XC
  if (f12 == NULL)
    {
      twobody_XC_slip (n, l, f);
    }
  else
    {
      f = f12->XC;
    }
  far[6] = 0.0;
  far[7] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[6] += f[i*2  ] / s11;
      far[7] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[7] *= -8.0 / l13;

  // YC
  if (f12 == NULL)
    {
      twobody_YC_slip (n, l, f);
    }
  else
    {
      f = f12->YC;
    }
  far[8] = 0.0;
  far[9] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[8] += f[i*2  ] / s11;
      far[9] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[9] *= 8.0 / l13;

  // FT version
  if (version == 1) return;

  // XG
  if (f12 == NULL)
    {
      twobody_XG_slip (n, l, f);
    }
  else
    {
      f = f12->XG;
    }
  far[10] = 0.0;
  far[11] = 0.0;
  s11 = s1l; // odd
  s12 = 1.0; // even
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[10] += f[i*2+1] / s11;
      far[11] += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[11] *= -4.0 / l12;

  // YG
  if (f12 == NULL)
    {
      twobody_YG_slip (n, l, f);
    }
  else
    {
      f = f12->YG;
    }
  far[12] = 0.0;
  far[13] = 0.0;
  s11 = s1l; // odd
  s12 = 1.0; // even
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[12] += f[i*2+1] / s11;
      far[13] += f[i*2  ] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[13] *= - 4.0 / l12;

  // YH
  if (f12 == NULL)
    {
      twobody_YH_slip (n, l, f);
    }
  else
    {
      f = f12->YH;
    }
  far[14] = 0.0;
  far[15] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[14] += f[i*2  ] / s11;
      far[15] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[15] *= 8.0 / l13;

  // XM
  if (f12 == NULL)
    {
      twobody_XM_slip (n, l, f);
    }
  else
    {
      f = f12->XM;
    }
  far[16] = 0.0;
  far[17] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[16] += f[i*2  ] / s11;
      far[17] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[17] *= 8.0 / l13;

  // YM
  if (f12 == NULL)
    {
      twobody_YM_slip (n, l, f);
    }
  else
    {
      f = f12->YM;
    }
  far[18] = 0.0;
  far[19] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[18] += f[i*2  ] / s11;
      far[19] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[19] *= 8.0 / l13;

  // ZM
  if (f12 == NULL)
    {
      twobody_ZM_slip (n, l, f);
    }
  else
    {
      f = f12->ZM;
    }
  far[20] = 0.0;
  far[21] = 0.0;
  s11 = 1.0; // even
  s12 = s1l; // odd
  s2 = s1l * s1l;
  for (i = 0; i < n/2; i ++)
    {
      far[20] += f[i*2  ] / s11;
      far[21] += f[i*2+1] / s12;

      s11 *= s2;
      s12 *= s2;
    }
  far[21] *= - 8.0 / l13;

  if (f12 == NULL)
    {
      free (f);
    }
}


/* calc scalar functions of two-body exact solution in resistance problem
 * INPUT
 *  version    : 0=F, 1=FT, 2=FTS.
 *  r          : distance between the two := x_b - x_a
 *  aa, ab     : radii for particles a(alpha) and b(beta)
 *  f12        : (struct twobody_slip_f *).
 *               you can give NULL for them.
 *               then, the coefs are calculated on-the-fly (terribly slow).
 *  n          : max order for the coefficients
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 *               (*** currently, lub form is not implemented ***)
 *  flag_scale : 0 no scaling, that is, in Jeffrey form
 *               1 for the dimensional form
 *               2 for the Stokesian dynamics form
 *  res [22]   :
 * OUTPUT
 *  res [22]   : scalar functions. the scaling is given by flag_scale.
 */
void
twobody_slip_scalars_res (int version,
			  double r,
			  double aa, double ab,
			  struct twobody_slip_f *f12,
			  int n, int flag_lub, int flag_scale,
			  double *res)
{
  double l = ab / aa;
  double s = 2.0 * r / (aa + ab);

  /* currently we have only far form
  if (flag_lub == 0)
    {
  */
      twobody_slip_far_with_f (version, f12, n, l, s,
			       f12->hat_g1, f12->hat_g2,
			       res);
  /*
    }
  else
    {
      twobody_slip_lub_with_f (version, f12, n, l, s,
			       f12->hat_g1, f12->hat_g2,
			       res);
    }
  */

  if (flag_scale == 1)
    {
      twobody_scale (version, res, aa, l);
    }
  else if (flag_scale == 2)
    {
      twobody_scale_SD (version, res, l);      
    }
}

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12,f21: (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  flag_lub   : 0 to use twobody_slip_far()
 *               1 to use twobody_slip_lub()
 *               (*** currently, lub form is not implemented ***)
 * OUTPUT
 *  lub [44] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *    8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *   12,13,14,15 : (XC11, XC12, XC21, XC22)
 *   16,17,18,19 : (YC11, YC12, YC21, YC22)
 *   20,21,22,23 : (XG11, XG12, XG21, XG22)
 *   24,25,26,27 : (YG11, YG12, YG21, YG22)
 *   28,29,30,31 : (YH11, YH12, YH21, YH22)
 *   32,33,34,35 : (XM11, XM12, XM21, XM22)
 *   36,37,38,39 : (YM11, YM12, YM21, YM22)
 *   40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 */
void
scalars_lub_slip_full (int version,
		       double r, double a1, double a2,
		       struct twobody_slip_f *f12,
		       struct twobody_slip_f *f21,
		       int n, int flag_lub,
		       double *lub)
{
  if (f12 == NULL || f21 == NULL)
    {
      fprintf (stderr, "minv needs twobody_slip_f table\n");
      exit (1);
    }

  /* some trivial check for sure */
  if (f12->hat_g1 != f21->hat_g2 ||
      f12->hat_g2 != f21->hat_g1 ||
      f12->slip_a1 != f21->slip_a2 ||
      f12->slip_a2 != f21->slip_a1)
    {
      fprintf (stderr, "inconsistent table f12 and f21\n");
      exit (1);
    }

  // zero clear
  int i;
  for (i = 0; i < 44; i ++)
    {
      lub[i] = 0.0;
    }

  double res12 [22];
  twobody_slip_scalars_res (version,
			    r, a1, a2, f12,
			    n, flag_lub, 1, // dimensional
			    res12);

  double res21 [22];
  twobody_slip_scalars_res (version,
			    r, a2, a1, f21,
			    n, flag_lub, 1, // dimensional
			    res21);

  double minv [44];
  if (version == 0)
    {
      // F version
      scalars_minv_f_poly (r, f12->slip_a1, f12->slip_a2, minv);
    }
  else if (version == 1)
    {
      // FT version
      scalars_minv_ft_poly (r, f12->slip_a1, f12->slip_a2, minv);
    }
  else
    {
      // FTS version
      scalars_minv_fts_poly (r, f12->slip_a1, f12->slip_a2, minv);
    }

  lub [ 0] = res12 [ 0] - minv [0]; // XA11
  lub [ 1] = res12 [ 1] - minv [1]; // XA12
  lub [ 2] = res21 [ 1] - minv [2]; // XA21
  lub [ 3] = res21 [ 0] - minv [3]; // XA22

  lub [ 4] = res12 [ 2] - minv [4]; // YA11
  lub [ 5] = res12 [ 3] - minv [5]; // YA12
  lub [ 6] = res21 [ 3] - minv [6]; // YA21
  lub [ 7] = res21 [ 2] - minv [7]; // YA22
  if (version == 0) return;

  lub  [8] = res12 [ 4] - minv [8]; // YB11
  lub  [9] = res12 [ 5] - minv [9]; // YB12
  lub [10] = res21 [ 5] - minv[10]; // YB21
  lub [11] = res21 [ 4] - minv[11]; // YB22

  lub [12] = res12 [ 6] - minv[12]; // XC11
  lub [13] = res12 [ 7] - minv[13]; // XC12
  lub [14] = res21 [ 7] - minv[14]; // XC21
  lub [15] = res21 [ 6] - minv[15]; // XC22

  lub [16] = res12 [ 8] - minv[16]; // YC11
  lub [17] = res12 [ 9] - minv[17]; // YC12
  lub [18] = res21 [ 9] - minv[18]; // YC21
  lub [19] = res21 [ 8] - minv[19]; // YC22
  if (version == 1) return;

  lub [20] = res12 [10] - minv[20]; // XG11
  lub [21] = res12 [11] - minv[21]; // XG12
  lub [22] = res21 [11] - minv[22]; // XG21
  lub [23] = res21 [10] - minv[23]; // XG22

  lub [24] = res12 [12] - minv[24]; // YG11
  lub [25] = res12 [13] - minv[25]; // YG12
  lub [26] = res21 [13] - minv[26]; // YG21
  lub [27] = res21 [12] - minv[27]; // YG22

  lub [28] = res12 [14] - minv[28]; // YH11
  lub [29] = res12 [15] - minv[29]; // YH12
  lub [30] = res21 [15] - minv[30]; // YH21
  lub [31] = res21 [14] - minv[31]; // YH22

  lub [32] = res12 [16] - minv[32]; // XM11
  lub [33] = res12 [17] - minv[33]; // XM12
  lub [34] = res21 [17] - minv[34]; // XM21
  lub [35] = res21 [16] - minv[35]; // XM22

  lub [36] = res12 [18] - minv[36]; // YM11
  lub [37] = res12 [19] - minv[37]; // YM12
  lub [38] = res21 [19] - minv[38]; // YM21
  lub [39] = res21 [18] - minv[39]; // YM22

  lub [40] = res12 [20] - minv[40]; // ZM11
  lub [41] = res12 [21] - minv[41]; // ZM12
  lub [42] = res21 [21] - minv[42]; // ZM21
  lub [43] = res21 [20] - minv[43]; // ZM22
}


/** F version **/

/* calculate f by u for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   u1 [3] : velocity of particle 1
 *   u2 [3] : velocity of particle 2
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   f1 [3] : force of particle 1
 *   f2 [3] : force of particle 2
 */
void
calc_lub_f_2b_slip (struct stokes *sys,
		    const double *u1, const double *u2,
		    const double *x1, const double *x2,
		    int i1, int i2,
		    double *f1, double *f2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (0, // F version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (0, // F version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];

  matrix_f_atimes (u1, f1,
		   ex, ey, ez,
		   xa11, ya11);
  matrix_f_atimes (u2, f1,
		   ex, ey, ez,
		   xa12, ya12);

  matrix_f_atimes (u2, f2,
		   -ex, -ey, -ez,
		   xa22, ya22);
  matrix_f_atimes (u1, f2,
		   -ex, -ey, -ez,
		   xa21, ya21);
}

/* calculate lub-matrix in F version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b_slip (struct stokes *sys,
		      int i, int j,
		      const double *x1, const double *x2,
		      int i1, int i2,
		      int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (0, // F version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (0, // F version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];

  matrix_f_ij (i, i,
	       ex, ey, ez,
	       xa11, ya11,
	       n, mat);
  matrix_f_ij (i, j,
	       ex, ey, ez,
	       xa12, ya12,
	       n, mat);

  matrix_f_ij (j, j,
	       -ex, -ey, -ez,
	       xa22, ya22,
	       n, mat);
  matrix_f_ij (j, i,
	       -ex, -ey, -ez,
	       xa21, ya21,
	       n, mat);
}


/** FT version **/

/* calculate ft by uo for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uo1 [6] : velocity, angular velocity
 *   uo2 [6] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   ft1 [6] : force, torque
 *   ft2 [6] :
 */
void
calc_lub_ft_2b_slip (struct stokes *sys,
		     const double *uo1, const double *uo2,
		     const double *x1, const double *x2,
		     int i1, int i2,
		     double *ft1, double *ft2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (1, // FT version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (1, // FT version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];
  yb11 = lub [8];
  yb12 = lub [9];
  yb21 = lub[10];
  yb22 = lub[11];
  xc11 = lub[12];
  xc12 = lub[13];
  xc21 = lub[14];
  xc22 = lub[15];
  yc11 = lub[16];
  yc12 = lub[17];
  yc21 = lub[18];
  yc22 = lub[19];

  matrix_ft_self_atimes (uo1, ft1,
			 ex, ey, ez,
			 xa11, ya11,
			 yb11,
			 xc11, yc11);
  matrix_ft_atimes (uo2, ft1,
		    ex, ey, ez,
		    xa12, ya12,
		    yb12,
		    xc12, yc12);

  matrix_ft_self_atimes (uo2, ft2,
			 -ex, -ey, -ez,
			 xa22, ya22,
			 yb22,
			 xc22, yc22);
  matrix_ft_atimes (uo1, ft2,
		    -ex, -ey, -ez,
		    xa21, ya21,
		    yb21,
		    xc21, yc21);
}

/* calculate lub-matrix in FT version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*6)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_ft_2b_slip (struct stokes *sys,
		       int i, int j,
		       const double *x1, const double *x2,
		       int i1, int i2,
		       int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (1, // FT version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (1, // FT version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];
  yb11 = lub [8];
  yb12 = lub [9];
  yb21 = lub[10];
  yb22 = lub[11];
  xc11 = lub[12];
  xc12 = lub[13];
  xc21 = lub[14];
  xc22 = lub[15];
  yc11 = lub[16];
  yc12 = lub[17];
  yc21 = lub[18];
  yc22 = lub[19];

  matrix_ft_ij (i, i,
		ex, ey, ez,
		xa11, ya11,
		yb11,
		xc11, yc11,
		n, mat);
  matrix_ft_ij (i, j,
		ex, ey, ez,
		xa12, ya12,
		yb12,
		xc12, yc12,
		n, mat);

  matrix_ft_ij (j, j,
		-ex, -ey, -ez,
		xa22, ya22,
		yb22,
		xc22, yc22,
		n, mat);
  matrix_ft_ij (j, i,
		-ex, -ey, -ez,
		xa21, ya21,
		yb21,
		xc21, yc21,
		n, mat);
}


/** FTS version **/

/* calculate fts by uoe for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b_slip (struct stokes *sys,
		      const double *uoe1, const double *uoe2,
		      const double *x1, const double *x2,
		      int i1, int i2,
		      double *fts1, double *fts2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (2, // FTS version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (2, // FTS version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  double xg11, xg12, xg21, xg22;
  double yg11, yg12, yg21, yg22;
  double yh11, yh12, yh21, yh22;
  double xm11, xm12, xm21, xm22;
  double ym11, ym12, ym21, ym22;
  double zm11, zm12, zm21, zm22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];
  yb11 = lub [8];
  yb12 = lub [9];
  yb21 = lub[10];
  yb22 = lub[11];
  xc11 = lub[12];
  xc12 = lub[13];
  xc21 = lub[14];
  xc22 = lub[15];
  yc11 = lub[16];
  yc12 = lub[17];
  yc21 = lub[18];
  yc22 = lub[19];
  xg11 = lub[20];
  xg12 = lub[21];
  xg21 = lub[22];
  xg22 = lub[23];
  yg11 = lub[24];
  yg12 = lub[25];
  yg21 = lub[26];
  yg22 = lub[27];
  yh11 = lub[28];
  yh12 = lub[29];
  yh21 = lub[30];
  yh22 = lub[31];
  xm11 = lub[32];
  xm12 = lub[33];
  xm21 = lub[34];
  xm22 = lub[35];
  ym11 = lub[36];
  ym12 = lub[37];
  ym21 = lub[38];
  ym22 = lub[39];
  zm11 = lub[40];
  zm12 = lub[41];
  zm21 = lub[42];
  zm22 = lub[43];

  matrix_fts_self_atimes (uoe1, fts1,
			  ex, ey, ez,
			  xa11, ya11,
			  yb11,
			  xc11, yc11,
			  xg11, yg11,
			  yh11,
			  xm11, ym11, zm11);
  matrix_fts_atimes (uoe2, fts1,
		     ex, ey, ez,
		     xa12, ya12,
		     yb12,
		     xc12, yc12,
		     xg12, yg12,
		     yh12,
		     xm12, ym12, zm12);

  matrix_fts_self_atimes (uoe2, fts2,
			  -ex, -ey, -ez,
			  xa22, ya22,
			  yb22,
			  xc22, yc22,
			  xg22, yg22,
			  yh22,
			  xm22, ym22, zm22);
  matrix_fts_atimes (uoe1, fts2,
		     -ex, -ey, -ez,
		     xa21, ya21,
		     yb21,
		     xc21, yc21,
		     xg21, yg21,
		     yh21,
		     xm21, ym21, zm21);
}

/* calculate lub-matrix in FTS version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*11)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_fts_2b_slip (struct stokes *sys,
			int i, int j,
			const double *x1, const double *x2,
			int i1, int i2,
			int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (2, // FTS version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (2, // FTS version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  double xg11, xg12, xg21, xg22;
  double yg11, yg12, yg21, yg22;
  double yh11, yh12, yh21, yh22;
  double xm11, xm12, xm21, xm22;
  double ym11, ym12, ym21, ym22;
  double zm11, zm12, zm21, zm22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];
  yb11 = lub [8];
  yb12 = lub [9];
  yb21 = lub[10];
  yb22 = lub[11];
  xc11 = lub[12];
  xc12 = lub[13];
  xc21 = lub[14];
  xc22 = lub[15];
  yc11 = lub[16];
  yc12 = lub[17];
  yc21 = lub[18];
  yc22 = lub[19];
  xg11 = lub[20];
  xg12 = lub[21];
  xg21 = lub[22];
  xg22 = lub[23];
  yg11 = lub[24];
  yg12 = lub[25];
  yg21 = lub[26];
  yg22 = lub[27];
  yh11 = lub[28];
  yh12 = lub[29];
  yh21 = lub[30];
  yh22 = lub[31];
  xm11 = lub[32];
  xm12 = lub[33];
  xm21 = lub[34];
  xm22 = lub[35];
  ym11 = lub[36];
  ym12 = lub[37];
  ym21 = lub[38];
  ym22 = lub[39];
  zm11 = lub[40];
  zm12 = lub[41];
  zm21 = lub[42];
  zm22 = lub[43];

  matrix_fts_ij (i, i,
		 ex, ey, ez,
		 xa11, ya11,
		 yb11,
		 xc11, yc11,
		 xg11, yg11,
		 yh11,
		 xm11, ym11, zm11,
		 n, mat);
  matrix_fts_ij (i, j,
		 ex, ey, ez,
		 xa12, ya12,
		 yb12,
		 xc12, yc12,
		 xg12, yg12,
		 yh12,
		 xm12, ym12, zm12,
		 n, mat);

  matrix_fts_ij (j, j,
		 -ex, -ey, -ez,
		 xa22, ya22,
		 yb22,
		 xc22, yc22,
		 xg22, yg22,
		 yh22,
		 xm22, ym22, zm22,
		 n, mat);
  matrix_fts_ij (j, i,
		 -ex, -ey, -ez,
		 xa21, ya21,
		 yb21,
		 xc21, yc21,
		 xg21, yg21,
		 yh21,
		 xm21, ym21, zm21,
		 n, mat);
}

