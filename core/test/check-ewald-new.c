/* test code for polydisperse bugs in ewald-new.c
 * Copyright (C) 2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include <stokes.h> // struct stokes
#include <ewald.h> // atimes_ewald_3all()

#include "check.h" // compare()


/** check routines **/

// compare mono (without a) with polydisperse systems with a = 1.0
int
check_atimes_ewald_3all_new_1b
(int version,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_ewald_3all_new_1b\n"
	       "(%d)"
	       ": start\n",
	       version);
    }


  double phi = 0.3;
  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  // mono
  struct stokes *sys0 = stokes_init ();
  CHECK_MALLOC (sys0, "check_atimes_ewald_3all_new_1b");
  sys0->version = version;

  // N=8
  int np = 8;
  int n = n0 * np;
  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->pos [3] = l;
  sys0->pos [4] = 0.0;
  sys0->pos [5] = 0.0;

  sys0->pos [6] = 0.0;
  sys0->pos [7] = l;
  sys0->pos [8] = 0.0;

  sys0->pos [9] = 0.0;
  sys0->pos[10] = 0.0;
  sys0->pos[11] = l;

  sys0->pos[12] = l;
  sys0->pos[13] = l;
  sys0->pos[14] = 0.0;

  sys0->pos[15] = l;
  sys0->pos[16] = 0.0;
  sys0->pos[17] = l;

  sys0->pos[18] = 0.0;
  sys0->pos[19] = l;
  sys0->pos[20] = l;

  sys0->pos[21] = l;
  sys0->pos[22] = l;
  sys0->pos[23] = l;


  sys0->periodic = 1;
  stokes_set_l (sys0, 2.0*l, 2.0*l, 2.0*l);


  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;

  double xi = xi_by_tratio (sys0, ewald_tr);
  stokes_set_xi (sys0, xi, ewald_eps);


  // poly (with a=1)
  struct stokes *sys1 = stokes_init ();
  CHECK_MALLOC (sys1, "check_atimes_ewald_3all_new_1b");
  sys1->version = version;

  stokes_set_np (sys1, np, np);

  for (int i = 0; i < np * 3; i ++)
    {
      sys1->pos [i] = sys0->pos [i];
    }

  double *a = (double *)malloc (sizeof (double) * np);
  for (int i = 0; i < np; i ++)
    {
      a[i] = 1.0;
    }
  stokes_set_radius (sys1, a);

  sys1->periodic = 1;
  stokes_set_l (sys1, 2.0*l, 2.0*l, 2.0*l);

  stokes_set_xi (sys1, xi, ewald_eps);


  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_ewald_3all_new_1b");
  for (int i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }


  // mono
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_ewald_3all_new_1b");

  atimes_ewald_3all (n, x, y1, (void *)sys0);

  // poly (a=1)
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_ewald_3all_new_1b");

  atimes_ewald_3all (n, x, y2, (void *)sys1);


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];
  // between y1 and y2
  for (int i = 0; i < n; i ++)
    {
      sprintf (label, " (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }


  free (y1);
  free (y2);
  free (x);
  free (a);
  stokes_free (sys0);
  stokes_free (sys1);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


// compare with polydisperse systems
// scaled by the factor "scale"
// F, T, S and U, O, E are properly scaled for comparison
// compare with polydisperse systems with a = 1.0 (monodisperse)
int
check_atimes_ewald_3all_new_2
(int version,
 double scale,
 double phi,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_ewald_3all_new_2\n"
	       "(%d)"
	       ": start\n",
	       version);
    }


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_ewald_3all_new_2");
  sys->version = version;


  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;


  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  // N=8
  int np = 8;
  int n = n0 * np;
  stokes_set_np (sys, np, np);


  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_nonewald_3all_new_2");

  double *a = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (a, "check_atimes_nonewald_3all_new_2");


  // system 1
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_nonewald_3all_new_2");

  {
    sys->pos [0] = 0.0;
    sys->pos [1] = 0.0;
    sys->pos [2] = 0.0;

    sys->pos [3] = l;
    sys->pos [4] = 0.0;
    sys->pos [5] = 0.0;

    sys->pos [6] = 0.0;
    sys->pos [7] = l;
    sys->pos [8] = 0.0;

    sys->pos [9] = 0.0;
    sys->pos[10] = 0.0;
    sys->pos[11] = l;

    sys->pos[12] = l;
    sys->pos[13] = l;
    sys->pos[14] = 0.0;

    sys->pos[15] = l;
    sys->pos[16] = 0.0;
    sys->pos[17] = l;

    sys->pos[18] = 0.0;
    sys->pos[19] = l;
    sys->pos[20] = l;

    sys->pos[21] = l;
    sys->pos[22] = l;
    sys->pos[23] = l;


    for (int i = 0; i < np; i ++)
      {
	a[i] = 1.0;
      }
    stokes_set_radius (sys, a);


    sys->periodic = 1;
    stokes_set_l (sys, 2.0*l, 2.0*l, 2.0*l);


    double xi = xi_by_tratio (sys, ewald_tr);
    stokes_set_xi (sys, xi, ewald_eps);


    for (int i = 0; i < n; i ++)
      {
	x[i] = 1.0;
      }


    atimes_ewald_3all (n, x, y1, (void *)sys);
  }


  // system 2
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_nonewald_3all_new_2");

  {
    for (int i = 0; i < np; i ++) {
      a[i] *= scale;

      int i3 = i * 3;
      sys->pos [i3 + 0] *= scale;
      sys->pos [i3 + 1] *= scale;
      sys->pos [i3 + 2] *= scale;
    }

    stokes_set_radius (sys, a);


    l *= scale;
    stokes_set_l (sys, 2.0*l, 2.0*l, 2.0*l);

    double xi = xi_by_tratio (sys, ewald_tr);
    stokes_set_xi (sys, xi, ewald_eps);


    if (version == 0)
      {
	// F version
	/*
	for (int i = 0; i < np; i ++)
	  {
	    int i3 = i * 3;
	    // F
	    x[i3    ] *= 1.0;
	    x[i3 + 1] *= 1.0;
	    x[i3 + 2] *= 1.0;
	  }
	*/
      }
    else if (version == 1)
      {
	// FT version
	for (int i = 0; i < np; i ++)
	  {
	    int i6 = i * 6;
	    // F
	    /*
	    x[i3    ] *= 1.0;
	    x[i3 + 1] *= 1.0;
	    x[i3 + 2] *= 1.0;
	    */
	    // T
	    x[i6 + 3] *= scale;
	    x[i6 + 4] *= scale;
	    x[i6 + 5] *= scale;
	  }
      }
    else
      {
	// FTS version
	for (int i = 0; i < np; i ++)
	  {
	    int i11 = i * 11;
	    // F
	    /*
	    x[i3    ] *= 1.0;
	    x[i3 + 1] *= 1.0;
	    x[i3 + 2] *= 1.0;
	    */
	    // T
	    x[i11 + 3] *= scale;
	    x[i11 + 4] *= scale;
	    x[i11 + 5] *= scale;
	    // S
	    x[i11 + 6] *= scale;
	    x[i11 + 7] *= scale;
	    x[i11 + 8] *= scale;
	    x[i11 + 9] *= scale;
	    x[i11 + 10] *= scale;
	  }
      }


    atimes_ewald_3all (n, x, y2, (void *)sys);


    // scale U, O, E for comparison
    double scale2 = scale * scale;
    if (version == 0)
      {
	// F version
	for (int i = 0; i < np; i ++)
	  {
	    int i3 = i * 3;
	    // U
	    y2[i3    ] *= scale;
	    y2[i3 + 1] *= scale;
	    y2[i3 + 2] *= scale;
	  }
      }
    else if (version == 1)
      {
	// FT version
	for (int i = 0; i < np; i ++)
	  {
	    int i6 = i * 6;
	    // U
	    y2[i6    ] *= scale;
	    y2[i6 + 1] *= scale;
	    y2[i6 + 2] *= scale;
	    // O
	    y2[i6 + 3] *= scale2;
	    y2[i6 + 4] *= scale2;
	    y2[i6 + 5] *= scale2;
	  }
      }
    else
      {
	// FTS version
	for (int i = 0; i < np; i ++)
	  {
	    int i11 = i * 11;
	    // U
	    y2[i11    ] *= scale;
	    y2[i11 + 1] *= scale;
	    y2[i11 + 2] *= scale;
	    // O
	    y2[i11 + 3] *= scale2;
	    y2[i11 + 4] *= scale2;
	    y2[i11 + 5] *= scale2;
	    // E
	    y2[i11 + 6] *= scale2;
	    y2[i11 + 7] *= scale2;
	    y2[i11 + 8] *= scale2;
	    y2[i11 + 9] *= scale2;
	    y2[i11 + 10] *= scale2;
	  }
      }
  }


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];
  // between y1 and y2
  for (int i = 0; i < n; i ++)
    {
      sprintf (label, " (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }


  free (y1);
  free (y2);
  free (a);
  free (x);
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
