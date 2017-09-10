/* test code for ewald-3fts-new.c
 * natural mobility problem with lubrication
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

#include <ewald-3fts.h> // solve_res_3fts_0()

#include "check.h" // compare()


/** check routines **/

// compare with new (mono) and new (poly a=1)
int
check_solve_mob_lub_3fts_b
(int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_mob_lub_3fts_b\n"
	       "start\n");
    }


  srand48(0);


  // monodisperse
  struct stokes *sys0 = stokes_init ();
  CHECK_MALLOC (sys0, "check_solve_mob_lub_3fts_b");
  sys0->version = 2; // FTS

  // N = 3
  int np = 3;

  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->pos [3] = 3.0 + drand48();
  sys0->pos [4] = 3.0 + drand48();
  sys0->pos [5] = 3.0 + drand48();

  sys0->pos [6] = -3.0 + drand48();
  sys0->pos [7] = -3.0 + drand48();
  sys0->pos [8] = -3.0 + drand48();

  stokes_set_iter (sys0, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  double Ui[3] = {0.1, 0.2, 0.3};
  double Oi[3] = {-0.1, -0.2, -0.3};
  double Ei[5] = {0.2, 0.3, -0.1, -0.2, -0.3};
  stokes_set_Ui (sys0, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys0, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys0, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  // polydisperse
  struct stokes *sys1 = stokes_init ();
  CHECK_MALLOC (sys1, "check_solve_mob_lub_3fts_b");
  sys1->version = 2; // FTS

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

  stokes_set_iter (sys1, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  stokes_set_Ui (sys1, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys1, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys1, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  double *f = (double *)malloc (sizeof (double) * np * 3);
  double *t = (double *)malloc (sizeof (double) * np * 3);
  double *e = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (t, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (e, "check_solve_mob_lub_3fts_b");
  for (int i = 0; i < np * 3; i ++)
    {
      f[i] = drand48();
      t[i] = drand48();
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e[i] = drand48();
    }

  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *s0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (o0, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (s0, "check_solve_mob_lub_3fts_b");

  double *u1 = (double *)malloc (sizeof (double) * np * 3);
  double *o1 = (double *)malloc (sizeof (double) * np * 3);
  double *s1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u1, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (o1, "check_solve_mob_lub_3fts_b");
  CHECK_MALLOC (s1, "check_solve_mob_lub_3fts_b");


  solve_mob_lub_3fts (sys0,
		      f, t, e,
		      u0, o0, s0);

  solve_mob_lub_3fts (sys1,
		      f, t, e,
		      u1, o1, s1);


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];

  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " U(%d) ", i);
      check += compare_max (u0[i], u1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " O(%d) ", i);
      check += compare_max (o0[i], o1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 5; i ++)
    {
      sprintf (label, " S(%d) ", i);
      check += compare_max (s0[i], s1[i], label, verbose, tiny, &max);
    }


  free (a);
  free (f);
  free (t);
  free (e);
  free (u0);
  free (o0);
  free (s0);
  free (u1);
  free (o1);
  free (s1);
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


// compare with mono and scaled poly (a = "scale")
int
check_solve_mob_lub_3fts_b2
(double scale,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_mob_lub_3fts_b2\n"
	       "(scale=%f)"
	       ": start\n",
	       scale);
    }


  srand48(0);


  // monodisplers
  struct stokes *sys0 = stokes_init ();
  CHECK_MALLOC (sys0, "check_solve_mob_lub_3fts_b2");
  sys0->version = 2; // FTS

  // N = 3
  int np = 3;

  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->pos [3] = 3.0 + drand48();
  sys0->pos [4] = 3.0 + drand48();
  sys0->pos [5] = 3.0 + drand48();

  sys0->pos [6] = -3.0 + drand48();
  sys0->pos [7] = -3.0 + drand48();
  sys0->pos [8] = -3.0 + drand48();

  stokes_set_iter (sys0, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  double Ui[3] = {0.1, 0.2, 0.3};
  double Oi[3] = {-0.1, -0.2, -0.3};
  double Ei[5] = {0.2, 0.3, -0.1, -0.2, -0.3};
  stokes_set_Ui (sys0, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys0, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys0, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);

  double *f0 = (double *)malloc (sizeof (double) * np * 3);
  double *t0 = (double *)malloc (sizeof (double) * np * 3);
  double *e0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f0, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (t0, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (e0, "check_solve_mob_lub_3fts_b2");
  for (int i = 0; i < np * 3; i ++)
    {
      f0[i] = drand48();
      t0[i] = drand48();
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e0[i] = drand48();
    }


  // polydisperse (a = "scale")
  struct stokes *sys1 = stokes_init ();
  CHECK_MALLOC (sys1, "check_solve_mob_lub_3fts_b2");
  sys1->version = 2; // FTS

  stokes_set_np (sys1, np, np);

  for (int i = 0; i < np * 3; i ++)
    {
      sys1->pos [i] = scale * sys0->pos [i];
    }

  double *a = (double *)malloc (sizeof (double) * np);
  for (int i = 0; i < np; i ++)
    {
      a[i] = scale;
    }
  stokes_set_radius (sys1, a);

  stokes_set_iter (sys1, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  double Ui1[3];
  double Oi1[3];
  double Ei1[5];
  // scale
  for (int i = 0; i < 3; i ++)
    {
      Ui1[i] = Ui[i] / scale;
      Oi1[i] = Oi[i] / (scale * scale);
    }
  for (int i = 0; i < 5; i ++)
    {
      Ei1[i] = Ei[i] / (scale * scale);
    }

  stokes_set_Ui (sys1, Ui1[0], Ui1[1], Ui1[2]);
  stokes_set_Oi (sys1, Oi1[0], Oi1[1], Oi1[2]);
  stokes_set_Ei (sys1, Ei1[0], Ei1[1], Ei1[2], Ei1[3], Ei1[4]);

  double *f1 = (double *)malloc (sizeof (double) * np * 3);
  double *t1 = (double *)malloc (sizeof (double) * np * 3);
  double *e1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f1, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (t1, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (e1, "check_solve_mob_lub_3fts_b2");
  // scale
  for (int i = 0; i < np * 3; i ++)
    {
      // F
      f1[i] = f0[i];
      // T
      t1[i] = t0[i] * scale;
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e1[i] = e0[i] / (scale * scale);
    }

  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *s0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (o0, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (s0, "check_solve_mob_lub_3fts_b2");

  double *u1 = (double *)malloc (sizeof (double) * np * 3);
  double *o1 = (double *)malloc (sizeof (double) * np * 3);
  double *s1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u1, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (o1, "check_solve_mob_lub_3fts_b2");
  CHECK_MALLOC (s1, "check_solve_mob_lub_3fts_b2");


  solve_mob_lub_3fts (sys0,
		      f0, t0, e0,
		      u0, o0, s0);

  solve_mob_lub_3fts (sys1,
		      f1, t1, e1,
		      u1, o1, s1);


  // scale
  for (int i = 0; i < np * 3; i ++)
    {
      // U
      u1[i] *= scale;
      // O
      o1[i] *= scale * scale;
    }
  for (int i = 0; i < np * 5; i ++)
    {
      s1[i] /= scale;
    }


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];

  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " U(%d) ", i);
      check += compare_max (u0[i], u1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " O(%d) ", i);
      check += compare_max (o0[i], o1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 5; i ++)
    {
      sprintf (label, " S(%d) ", i);
      check += compare_max (s0[i], s1[i], label, verbose, tiny, &max);
    }


  free (a);
  free (f0);
  free (t0);
  free (e0);
  free (f1);
  free (t1);
  free (e1);
  free (u0);
  free (o0);
  free (s0);
  free (u1);
  free (o1);
  free (s1);
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


// ** peroidic systems **

// compare with new (mono) and new (poly a=1)
int
check_solve_mob_lub_3fts_d
(int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_mob_lub_3fts_d\n"
	       "start\n");
    }


  srand48(0);


  // mono
  struct stokes *sys0 = stokes_init ();
  CHECK_MALLOC (sys0, "check_solve_mob_lub_3fts_d");
  sys0->version = 2; // FTS

  double phi = 0.3;
  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  /*
  // N=1
  int np = 1;
  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->periodic = 1;
  stokes_set_l (sys0, l, l, l);
  */
  /*
  // N=2
  int np = 2;
  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->pos [3] = l;
  sys0->pos [4] = 0.0;
  sys0->pos [5] = 0.0;

  sys0->periodic = 1;
  stokes_set_l (sys0, 2.0 * l, l, l);
  */
  /*
  // N=4
  int np = 4;
  stokes_set_np (sys0, np, np);

  sys0->pos [0] = 0.0;
  sys0->pos [1] = 0.0;
  sys0->pos [2] = 0.0;

  sys0->pos [3] = l;
  sys0->pos [4] = 0.0;
  sys0->pos [5] = 0.0;

  sys0->pos [6] = 2.0 * l;
  sys0->pos [7] = 0.0;
  sys0->pos [8] = 0.0;

  sys0->pos [9] = 3.0 * l;
  sys0->pos[10] = 0.0;
  sys0->pos[11] = 0.0;

  sys0->periodic = 1;
  stokes_set_l (sys0, 4.0 * l, l, l);
  */
  /*
  // N=4
  int np = 4;
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

  sys0->pos [9] = l;
  sys0->pos[10] = l;
  sys0->pos[11] = 0.0;

  sys0->periodic = 1;
  stokes_set_l (sys0, 2.0 * l, 2.0 * l, l);
  */
  // N=8
  int np = 8;
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

  stokes_set_iter (sys0, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  double Ui[3] = {0.1, 0.2, 0.3};
  double Oi[3] = {-0.1, -0.2, -0.3};
  double Ei[5] = {0.2, 0.3, -0.1, -0.2, -0.3};
  stokes_set_Ui (sys0, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys0, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys0, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  // poly (a=1)
  struct stokes *sys1 = stokes_init ();
  CHECK_MALLOC (sys1, "check_solve_mob_lub_3fts_d");
  sys1->version = 2; // FTS

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
  stokes_set_l (sys1, 2.0*l, 2.0*l, 2.0*l); // N=8
  //stokes_set_l (sys1, 2.0 * l, l, l); // N=2
  //stokes_set_l (sys1, 4.0 * l, l, l); // N=4
  //stokes_set_l (sys1, 2.0 * l, 2.0 * l, l); // N=4
  //stokes_set_l (sys1, l, l, l); // N=1

  stokes_set_xi (sys1, xi, ewald_eps);

  stokes_set_iter (sys1, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  stokes_set_Ui (sys1, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys1, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys1, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  double *f = (double *)malloc (sizeof (double) * np * 3);
  double *t = (double *)malloc (sizeof (double) * np * 3);
  double *e = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (t, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (e, "check_solve_mob_lub_3fts_d");
  for (int i = 0; i < np * 3; i ++)
    {
      f[i] = drand48();
      t[i] = drand48();
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e[i] = drand48();
    }

  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *s0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (o0, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (s0, "check_solve_mob_lub_3fts_d");

  double *u1 = (double *)malloc (sizeof (double) * np * 3);
  double *o1 = (double *)malloc (sizeof (double) * np * 3);
  double *s1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u1, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (o1, "check_solve_mob_lub_3fts_d");
  CHECK_MALLOC (s1, "check_solve_mob_lub_3fts_d");


  // mono
  solve_mob_lub_3fts (sys0,
		      f, t, e,
		      u0, o0, s0);

  // poly(a=1)
  solve_mob_lub_3fts (sys1,
		      f, t, e,
		      u1, o1, s1);


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];

  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " U(%d) ", i);
      check += compare_max (u0[i], u1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " O(%d) ", i);
      check += compare_max (o0[i], o1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 5; i ++)
    {
      sprintf (label, " S(%d) ", i);
      check += compare_max (s0[i], s1[i], label, verbose, tiny, &max);
    }


  free (a);
  free (f);
  free (t);
  free (e);
  free (u0);
  free (o0);
  free (s0);
  free (u1);
  free (o1);
  free (s1);
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


// compare with mono and scaled poly (a = "scale")
int
check_solve_mob_lub_3fts_d2
(double scale,
 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_mob_lub_3fts_d2\n"
	       "(scale=%f)"
	       ": start\n",
	       scale);
    }


  srand48(0);


  // mono
  struct stokes *sys0 = stokes_init ();
  CHECK_MALLOC (sys0, "check_solve_mob_lub_3fts_d2");
  sys0->version = 2; // FTS

  double phi = 0.3;
  double l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  // N=8
  int np = 8;
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


  double *a0 = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (a0, "check_solve_mob_lub_3fts_d2");
  for (int i = 0; i < np; i ++)
    {
      a0[i] = 1.0;
    }
  stokes_set_radius (sys0, a0);


  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;

  double xi = xi_by_tratio (sys0, ewald_tr);
  stokes_set_xi (sys0, xi, ewald_eps);

  stokes_set_iter (sys0, "gmres", 2000, 20, 1.0e-6, 1, stderr);


  double Ui[3] = {0.1, 0.2, 0.3};
  double Oi[3] = {-0.1, -0.2, -0.3};
  double Ei[5] = {0.2, 0.3, -0.1, -0.2, -0.3};
  stokes_set_Ui (sys0, Ui[0], Ui[1], Ui[2]);
  stokes_set_Oi (sys0, Oi[0], Oi[1], Oi[2]);
  stokes_set_Ei (sys0, Ei[0], Ei[1], Ei[2], Ei[3], Ei[4]);


  double *f0 = (double *)malloc (sizeof (double) * np * 3);
  double *t0 = (double *)malloc (sizeof (double) * np * 3);
  double *e0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f0, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (t0, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (e0, "check_solve_mob_lub_3fts_d2");
  for (int i = 0; i < np * 3; i ++)
    {
      f0[i] = drand48();
      t0[i] = drand48();
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e0[i] = drand48();
    }


  // polydisperse (a = "scale")
  struct stokes *sys1 = stokes_init ();
  CHECK_MALLOC (sys1, "check_solve_mob_lub_3fts_d2");
  sys1->version = 2; // FTS

  stokes_set_np (sys1, np, np);

  for (int i = 0; i < np * 3; i ++)
    {
      sys1->pos [i] = scale * sys0->pos [i];
    }

  double *a1 = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (a1, "check_solve_mob_lub_3fts_d2");
  for (int i = 0; i < np; i ++)
    {
      a1[i] = a0[i] * scale;
    }
  stokes_set_radius (sys1, a1);

  sys1->periodic = 1;
  double l1 = l * scale;
  stokes_set_l (sys1, 2.0*l1, 2.0*l1, 2.0*l1);

  //double xi1 = xi_by_tratio (sys1, ewald_tr);
  double xi1 = xi / scale;
  stokes_set_xi (sys1, xi1, ewald_eps);

  stokes_set_iter (sys1, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  double Ui1[3];
  double Oi1[3];
  double Ei1[5];
  // scale
  for (int i = 0; i < 3; i ++)
    {
      Ui1[i] = Ui[i] / scale;
      Oi1[i] = Oi[i] / (scale * scale);
    }
  for (int i = 0; i < 5; i ++)
    {
      Ei1[i] = Ei[i] / (scale * scale);
    }

  stokes_set_Ui (sys1, Ui1[0], Ui1[1], Ui1[2]);
  stokes_set_Oi (sys1, Oi1[0], Oi1[1], Oi1[2]);
  stokes_set_Ei (sys1, Ei1[0], Ei1[1], Ei1[2], Ei1[3], Ei1[4]);


  double *f1 = (double *)malloc (sizeof (double) * np * 3);
  double *t1 = (double *)malloc (sizeof (double) * np * 3);
  double *e1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (f1, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (t1, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (e1, "check_solve_mob_lub_3fts_d2");
  // scale
  for (int i = 0; i < np * 3; i ++)
    {
      // F
      f1[i] = f0[i];
      // T
      t1[i] = t0[i] * scale;
    }
  for (int i = 0; i < np * 5; i ++)
    {
      e1[i] = e0[i] / (scale * scale);
    }


  double *u0 = (double *)malloc (sizeof (double) * np * 3);
  double *o0 = (double *)malloc (sizeof (double) * np * 3);
  double *s0 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u0, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (o0, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (s0, "check_solve_mob_lub_3fts_d2");

  double *u1 = (double *)malloc (sizeof (double) * np * 3);
  double *o1 = (double *)malloc (sizeof (double) * np * 3);
  double *s1 = (double *)malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (u1, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (o1, "check_solve_mob_lub_3fts_d2");
  CHECK_MALLOC (s1, "check_solve_mob_lub_3fts_d2");


  solve_mob_lub_3fts (sys0,
		      f0, t0, e0,
		      u0, o0, s0);

  solve_mob_lub_3fts (sys1,
		      f1, t1, e1,
		      u1, o1, s1);


  // scale
  for (int i = 0; i < np * 3; i ++)
    {
      // U
      u1[i] *= scale;
      // O
      o1[i] *= scale * scale;
    }
  for (int i = 0; i < np * 5; i ++)
    {
      s1[i] /= scale;
    }


  // compare
  int check = 0;
  double max = 0.0;

  char label [80];

  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " U(%d) ", i);
      check += compare_max (u0[i], u1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 3; i ++)
    {
      sprintf (label, " O(%d) ", i);
      check += compare_max (o0[i], o1[i], label, verbose, tiny, &max);
    }
  for (int i = 0; i < np * 5; i ++)
    {
      sprintf (label, " S(%d) ", i);
      check += compare_max (s0[i], s1[i], label, verbose, tiny, &max);
    }


  free (a0);
  free (a1);

  free (f0);
  free (t0);
  free (e0);
  free (f1);
  free (t1);
  free (e1);
  free (u0);
  free (o0);
  free (s0);
  free (u1);
  free (o1);
  free (s1);

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
