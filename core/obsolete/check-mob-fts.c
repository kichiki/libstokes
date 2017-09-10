/* test code for mob-fts problem
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-mob-fts.c,v 1.1 2007/12/22 18:23:35 kichiki Exp $
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
#include "check.h" // compare()

#include <stokes.h>
#include <ewald-3fts.h>
#include <ewald-3fts-matrix.h>


/* compare solve_mob_3fts_matrix() with solve_mob_3fts_matrix_0()
 * and solve_mob_3fts_0().
 */
int
check_mob_fts (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_mob_fts : start\n");
    }

  int check = 0;
  double max = 0.0;


  /**
   * initialization
   */
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_minv_FU");

  sys->version = 2; // FTS
  int np = 9;
  stokes_set_np (sys, np, np);

  double *pos = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (pos, "check_minv_FU");
  pos[ 0] = 4.841648; pos[ 1] = 4.345871; pos[ 2] = 2.411424;
  pos[ 3] = 1.484547; pos[ 4] = 2.931487; pos[ 5] = 2.321764;
  pos[ 6] = 4.494053; pos[ 7] = 1.739099; pos[ 8] = 2.732656;
  pos[ 9] = 1.436324; pos[10] = 5.216197; pos[11] = 5.634981;
  pos[12] = 2.707117; pos[13] = 0.577389; pos[14] = 2.489404;
  pos[15] = 3.502632; pos[16] = 5.234916; pos[17] = 5.542215;
  pos[18] = 4.615775; pos[19] = 3.146068; pos[20] = 0.080511;
  pos[21] = 5.510249; pos[22] = 1.060571; pos[23] = 0.710376;
  pos[24] = 1.143043; pos[25] = 2.055843; pos[26] = 4.194753;
  stokes_set_pos (sys, pos);
  free (pos);
  
  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-12, 1, stdout);

  sys->periodic = 1;
  stokes_set_l (sys, 5.733683, 5.733683, 5.733683);
  double xi = xi_by_tratio (sys, 1.0);
  stokes_set_xi (sys, xi, 1.0e-12);


  int np3 = np * 3;
  int np5 = np * 5;
  double *f = (double *)malloc (sizeof (double) * np3);
  double *t = (double *)malloc (sizeof (double) * np3);
  double *s = (double *)malloc (sizeof (double) * np5);
  double *u = (double *)malloc (sizeof (double) * np3);
  double *o = (double *)malloc (sizeof (double) * np3);
  double *e = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (f, "check_mob_fts");
  CHECK_MALLOC (t, "check_mob_fts");
  CHECK_MALLOC (s, "check_mob_fts");
  CHECK_MALLOC (u, "check_mob_fts");
  CHECK_MALLOC (o, "check_mob_fts");
  CHECK_MALLOC (e, "check_mob_fts");

  int i;
  int ii;
  srand48(0);
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  f[i*3+ii] = drand48();
	  t[i*3+ii] = drand48();
	}
      for (ii = 0; ii < 5; ii ++)
	{
	  // e = 0
	  e[i*5+ii] = 0.0;
	}
    }

  // part 0 -- by solve_mob_3fts_matrix()
  solve_mob_3fts_matrix (sys, f, t, e, u, o, s);

  // part 1 -- by solve_mob_3fts_matrix_0()
  double *u1 = (double *)malloc (sizeof (double) * np3);
  double *o1 = (double *)malloc (sizeof (double) * np3);
  double *s1 = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (u1, "check_mob_fts");
  CHECK_MALLOC (o1, "check_mob_fts");
  CHECK_MALLOC (s1, "check_mob_fts");

  solve_mob_3fts_matrix_0 (sys, f, t, e, u1, o1, s1);

  char label_[80];
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "u1[%d]", i);
      check += compare_max (u[i], u1[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "o1[%d]", i);
      check += compare_max (o[i], o1[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*5; i ++)
    {
      sprintf (label_, "s1[%d]", i);
      check += compare_max (s[i], s1[i], label_, verbose, tiny, &max);
    }
  free (u1);
  free (o1);
  free (s1);

  // part 2 -- by solve_mob_3fts_0()
  double *u2 = (double *)malloc (sizeof (double) * np3);
  double *o2 = (double *)malloc (sizeof (double) * np3);
  double *s2 = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (u2, "check_mob_fts");
  CHECK_MALLOC (o2, "check_mob_fts");
  CHECK_MALLOC (s2, "check_mob_fts");

  solve_mob_3fts_0 (sys, NULL, f, t, e, u2, o2, s2);

  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "u2[%d]", i);
      check += compare_max (u[i], u2[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "o2[%d]", i);
      check += compare_max (o[i], o2[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*5; i ++)
    {
      sprintf (label_, "s2[%d]", i);
      check += compare_max (s[i], s2[i], label_, verbose, tiny, &max);
    }
  free (u2);
  free (o2);
  free (s2);

  free (f);
  free (t);
  free (s);
  free (u);
  free (o);
  free (e);

  stokes_free (sys);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

/* compare solve_mob_lub_3fts_matrix()
 * with solve_mob_lub_3fts_matrix_0()
 * and solve_mob_lub_3fts_0().
 */
int
check_mob_lub_fts (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_mob_lub_fts : start\n");
    }

  int check = 0;
  double max = 0.0;


  /**
   * initialization
   */
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_minv_FU");

  sys->version = 2; // FTS
  int np = 9;
  stokes_set_np (sys, np, np);

  double *pos = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (pos, "check_minv_FU");
  pos[ 0] = 4.841648; pos[ 1] = 4.345871; pos[ 2] = 2.411424;
  pos[ 3] = 1.484547; pos[ 4] = 2.931487; pos[ 5] = 2.321764;
  pos[ 6] = 4.494053; pos[ 7] = 1.739099; pos[ 8] = 2.732656;
  pos[ 9] = 1.436324; pos[10] = 5.216197; pos[11] = 5.634981;
  pos[12] = 2.707117; pos[13] = 0.577389; pos[14] = 2.489404;
  pos[15] = 3.502632; pos[16] = 5.234916; pos[17] = 5.542215;
  pos[18] = 4.615775; pos[19] = 3.146068; pos[20] = 0.080511;
  pos[21] = 5.510249; pos[22] = 1.060571; pos[23] = 0.710376;
  pos[24] = 1.143043; pos[25] = 2.055843; pos[26] = 4.194753;
  stokes_set_pos (sys, pos);
  free (pos);
  
  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-12, 1, stdout);

  sys->periodic = 1;
  stokes_set_l (sys, 5.733683, 5.733683, 5.733683);
  double xi = xi_by_tratio (sys, 1.0);
  stokes_set_xi (sys, xi, 1.0e-12);


  int np3 = np * 3;
  int np5 = np * 5;
  double *f = (double *)malloc (sizeof (double) * np3);
  double *t = (double *)malloc (sizeof (double) * np3);
  double *s = (double *)malloc (sizeof (double) * np5);
  double *u = (double *)malloc (sizeof (double) * np3);
  double *o = (double *)malloc (sizeof (double) * np3);
  double *e = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (f, "check_mob_lub_fts");
  CHECK_MALLOC (t, "check_mob_lub_fts");
  CHECK_MALLOC (s, "check_mob_lub_fts");
  CHECK_MALLOC (u, "check_mob_lub_fts");
  CHECK_MALLOC (o, "check_mob_lub_fts");
  CHECK_MALLOC (e, "check_mob_lub_fts");

  int i;
  int ii;
  srand48(0);
  for (i = 0; i < np; i ++)
    {
      for (ii = 0; ii < 3; ii ++)
	{
	  f[i*3+ii] = drand48();
	  t[i*3+ii] = drand48();
	}
      for (ii = 0; ii < 5; ii ++)
	{
	  // e = 0
	  e[i*5+ii] = 0.0;
	}
    }

  // part 0 -- by solve_mob_lub_3fts_matrix()
  solve_mob_lub_3fts_matrix (sys, f, t, e, u, o, s);

  // part 1 -- by solve_mob_3fts_matrix_0()
  double *u1 = (double *)malloc (sizeof (double) * np3);
  double *o1 = (double *)malloc (sizeof (double) * np3);
  double *s1 = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (u1, "check_mob_lub_fts");
  CHECK_MALLOC (o1, "check_mob_lub_fts");
  CHECK_MALLOC (s1, "check_mob_lub_fts");

  solve_mob_lub_3fts_matrix_0 (sys, f, t, e, u1, o1, s1);

  char label_[80];
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "u1[%d]", i);
      check += compare_max (u[i], u1[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "o1[%d]", i);
      check += compare_max (o[i], o1[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*5; i ++)
    {
      sprintf (label_, "s1[%d]", i);
      check += compare_max (s[i], s1[i], label_, verbose, tiny, &max);
    }
  free (u1);
  free (o1);
  free (s1);

  // part 2 -- by solve_mob_3fts_0()
  double *u2 = (double *)malloc (sizeof (double) * np3);
  double *o2 = (double *)malloc (sizeof (double) * np3);
  double *s2 = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (u2, "check_mob_lub_fts");
  CHECK_MALLOC (o2, "check_mob_lub_fts");
  CHECK_MALLOC (s2, "check_mob_lub_fts");

  solve_mob_lub_3fts_0 (sys, f, t, e, u2, o2, s2);

  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "u2[%d]", i);
      check += compare_max (u[i], u2[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*3; i ++)
    {
      sprintf (label_, "o2[%d]", i);
      check += compare_max (o[i], o2[i], label_, verbose, tiny, &max);
    }
  for (i = 0; i < np*5; i ++)
    {
      sprintf (label_, "s2[%d]", i);
      check += compare_max (s[i], s2[i], label_, verbose, tiny, &max);
    }
  //iter_free (iter);
  free (u2);
  free (o2);
  free (s2);

  free (f);
  free (t);
  free (s);
  free (u);
  free (o);
  free (e);

  stokes_free (sys);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
