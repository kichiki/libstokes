/* test code for polydisperse bugs in non-ewald-new.c
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
#include <non-ewald.h> // atimes_nonewald_3all()
#include <non-ewald-new.h> // atimes_nonewald_3all_new()

#include "check.h" // compare()


/** check routines **/

// compare with monodisperse systems
int
check_atimes_nonewald_3all_new_0
(int version,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_nonewald_3all_new_0\n"
	       "(%d)"
	       ": start\n",
	       version);
    }


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_nonewald_3all_new_0");
  sys->version = version;


  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  // N = 3
  int np = 3;
  int n = n0 * np;

  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;

  sys->pos [3] = 10.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = 0.0;

  sys->pos [6] = -10.0;
  sys->pos [7] = 0.0;
  sys->pos [8] = 0.0;


  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_nonewald_3all_new_0");
  for (int i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }


  // bug fixed version
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_nonewald_3all_new_0");

  atimes_nonewald_3all_new (n, x, y1, (void *)sys);

  // the old version
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_nonewald_3all_new_0");

  atimes_nonewald_3all (n, x, y2, (void *)sys);


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
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

// compare with polydisperse systems with a = 1.0 (monodisperse)
int
check_atimes_nonewald_3all_new_1
(int version,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_nonewald_3all_new_1\n"
	       "(%d)"
	       ": start\n",
	       version);
    }


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_nonewald_3all_new_1");
  sys->version = version;


  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  // N = 3
  int np = 3;
  int n = n0 * np;

  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;

  sys->pos [3] = 10.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = 0.0;

  sys->pos [6] = -10.0;
  sys->pos [7] = 0.0;
  sys->pos [8] = 0.0;

  double a [3] = {1.0, 1.0, 1.0};
  stokes_set_radius (sys, a);


  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_nonewald_3all_new_1");
  for (int i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }


  // bug fixed version
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_nonewald_3all_new_1");

  atimes_nonewald_3all_new (n, x, y1, (void *)sys);

  // the old version
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_nonewald_3all_new_1");

  atimes_nonewald_3all (n, x, y2, (void *)sys);


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
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

// compare with polydisperse systems
// system 1) a1 = 1  at (0, 0, 0), a2 = 10  at (20, 0, 0)
// system 2) scaled by the factor "scale"
// F, T, S and U, O, E are properly scaled for comparison
int
check_atimes_nonewald_3all_new_2
(int version,
 double scale,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_nonewald_3all_new_2\n"
	       "(%d)"
	       ": start\n",
	       version);
    }


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_nonewald_3all_new_2");
  sys->version = version;


  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  // N = 2
  //int np = 2;
  int np = 3;
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
    a[0] = 1.0;
    sys->pos [0] = 0.0;
    sys->pos [1] = 0.0;
    sys->pos [2] = 0.0;

    a[1] = 10.0;
    sys->pos [3] = 20.0;
    sys->pos [4] = 10.0;
    sys->pos [5] = 5.0;

    a[2] = 0.1;
    sys->pos [6] = -2.0;
    sys->pos [7] = 2.0;
    sys->pos [8] = 2.0;

    stokes_set_radius (sys, a);


    for (int i = 0; i < n; i ++)
      {
	//x[i] = 5.0;
	x[i] = 5.0 + 2.0 * (drand48 () - 0.5);
      }


    atimes_nonewald_3all_new (n, x, y1, (void *)sys);
    //atimes_nonewald_3all (n, x, y1, (void *)sys);
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


    atimes_nonewald_3all_new (n, x, y2, (void *)sys);
    //atimes_nonewald_3all (n, x, y2, (void *)sys);

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


