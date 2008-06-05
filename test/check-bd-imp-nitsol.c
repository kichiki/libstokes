/* test code for bd-imp-nitsol.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd-imp-nitsol.c,v 1.1 2008/06/05 03:23:54 kichiki Exp $
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
#include "check.h" // compare_max()

#include <stokes.h>
#include <bench.h>
#include <brownian.h>
#include <bd-imp.h>
#include "bd-imp-nitsol.h"
#include <bonds.h> // bonds_init(), bonds_free()
#include <angles.h> // angles_init(), angles_free()



int
check_BD_imp_NITSOL (int version, int np,
		     int flag_lub, int flag_mat, int flag_Q, double dt,
		     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BD_imp_NITSOL\n"
	       "(ver=%d,N=%d,lub=%d,mat=%d,Q=%d,dt=%f)",
	       version, np, flag_lub, flag_mat, flag_Q, dt);
    }

  int check = 0;
  double max = 0.0;


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_BD_imp_NITSOL");

  sys->version = version;

  //int np = 3;
  int nm = np;
  stokes_set_np (sys, np, nm);
  sys->periodic = 0; // non-periodic

  int np3 = np * 3;
  double *pos = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (pos, "check_BD_imp_NITSOL");

  // straight configuration
  int i;
  srand48 (0);
  for (i = 0; i < nm; i ++)
    {
      int i3 = i * 3;
      pos[i3 + 0] = (double)i * 2.5;
      pos[i3 + 1] = drand48() - 0.5;
      pos[i3 + 2] = drand48() - 0.5;
    }

  stokes_set_pos (sys, pos);

  int nm3 = nm * 3;
  int nm5 = nm * 5;
  double *F = (double *)malloc (sizeof (double) * nm3);
  double *T = (double *)malloc (sizeof (double) * nm3);
  double *E = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (F, "check_BD_imp_NITSOL");
  CHECK_MALLOC (T, "check_BD_imp_NITSOL");
  CHECK_MALLOC (E, "check_BD_imp_NITSOL");
  for (i = 0; i < nm; i ++)
    {
      int i3 = i * 3;
      int i5 = i * 5;

      F [i3  ] = 0.0;
      F [i3+1] = 0.0;
      F [i3+2] =-1.0;

      T [i3  ] = 0.0;
      T [i3+1] = 0.0;
      T [i3+2] = 0.0;

      E [i5  ] = 0.0;
      E [i5+1] = 0.0;
      E [i5+2] = 0.0;
      E [i5+3] = 0.0;
      E [i5+4] = 0.0;
    }

  if (flag_mat == 0)
    {
      stokes_set_iter (sys, "otmk", 2000, 20, 1.0e-6, 0, NULL);
    }

  struct bonds *bonds = bonds_init ();
  struct BD_params *BD
    = BD_params_init (sys,
		      0, // BD_seed,
		      F, T, E,
		      NULL, NULL, NULL, // uf, of, ef,
		      0, // int flag_noHI,
		      flag_lub,
		      flag_mat,
		      0.0, // st,
		      bonds,
		      1.0, // gamma,
		      NULL, // ev,
		      NULL, // ang
		      NULL, // ev_dh,
		      NULL, // ev_LJ,
		      NULL, // cf,
		      flag_Q,
		      -1.0, // peclet,
		      1.0e-6, // ode_eps,
		      0, // n_minv, // n of chebyshev for minv
		      0, // n_lub,  // n of chebyshev for lub
		      0, // BD_scheme, (0=mid, 1=BB, 2=BM, 3=JGdP, 4=siPC)
		      100, // BB_n,
		      1.0,  // double rmin
		      1.0e-12 //dt_lim
		      );
  CHECK_MALLOC (BD, "check_BD_imp_NITSOL");

  double *x_JGdP = (double *)malloc (sizeof (double) * np3);
  double *x_siPC  = (double *)malloc (sizeof (double) * np3);
  double *x_NJGdP   = (double *)malloc (sizeof (double) * np3);
  double *x_NsiPC   = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (x_JGdP, "check_BD_imp_NITSOL");
  CHECK_MALLOC (x_siPC,  "check_BD_imp_NITSOL");
  CHECK_MALLOC (x_NJGdP,   "check_BD_imp_NITSOL");
  CHECK_MALLOC (x_NsiPC,   "check_BD_imp_NITSOL");
  for (i = 0; i < np3; i ++)
    {
      x_JGdP[i]  = pos[i];
      x_siPC[i]  = pos[i];
      x_NJGdP[i] = pos[i];
      x_NsiPC[i] = pos[i];
    }

  int np4 = np * 4;
  double *q_JGdP = NULL;
  double *q_siPC  = NULL;
  double *q_NJGdP   = NULL;
  double *q_NsiPC   = NULL;
  if (sys->version > 0 && flag_Q != 0)
    {
      q_JGdP = (double *)malloc (sizeof (double) * np4);
      q_siPC  = (double *)malloc (sizeof (double) * np4);
      q_NJGdP   = (double *)malloc (sizeof (double) * np4);
      q_NsiPC   = (double *)malloc (sizeof (double) * np4);
      CHECK_MALLOC (q_JGdP, "check_BD_imp_NITSOL");
      CHECK_MALLOC (q_siPC,  "check_BD_imp_NITSOL");
      CHECK_MALLOC (q_NJGdP,   "check_BD_imp_NITSOL");
      CHECK_MALLOC (q_NsiPC,   "check_BD_imp_NITSOL");
      for (i = 0; i < np4; i ++)
	{
	  q_JGdP[i]  = 0.0;
	  q_siPC[i]  = 0.0;
	  q_NJGdP[i] = 0.0;
	  q_NsiPC[i] = 0.0;
	}
      // set Q4 = 1 as an initial condition
      for (i = 0; i < np; i ++)
	{
	  q_JGdP[i*4+3]  = 1.0;
	  q_siPC[i*4+3]  = 1.0;
	  q_NJGdP[i*4+3] = 1.0;
	  q_NsiPC[i*4+3] = 1.0;
	}
    }


  /**
   * BD implicit scheme
   */
  struct BD_imp *BDimp;

  BD->scheme = 3; // JGdP
  BDimp = BD_imp_init (BD,
		       1000,  // itmax
		       1.0e-6 // eps
		       );
  CHECK_MALLOC (BDimp, "check_BD_imp_NITSOL");

  double t = 0.0;
  double t_1 = ptime_ms_d();
  double dt_JGdP  = BD_evolve_JGdP00        (t, BDimp, x_JGdP,  q_JGdP,  dt);
  double t_2 = ptime_ms_d();
  double dt_NJGdP = BD_evolve_NITSOL_JGdP00 (t, BDimp, x_NJGdP, q_NJGdP, dt);
  double t_3 = ptime_ms_d();

  BD_imp_free (BDimp);
  BD->scheme = 4; // siPC
  BDimp = BD_imp_init (BD,
		       1000,  // itmax
		       1.0e-6 // eps
		       );
  CHECK_MALLOC (BDimp, "check_BD_imp_NITSOL");


  double t_4 = ptime_ms_d();
  double dt_siPC  = BD_evolve_imp_PC        (t, BDimp, x_siPC,  q_siPC,  dt);
  double t_5 = ptime_ms_d();
  double dt_NsiPC = BD_evolve_NITSOL_imp_PC (t, BDimp, x_NsiPC, q_NsiPC, dt);
  double t_6 = ptime_ms_d();

  fprintf (stdout, "BENCHMARK:\n");
  fprintf (stdout, "  JGdP : %e %e, ratio = %f\n",
	   (t_2 - t_1),
	   (t_3 - t_2),
	   (t_2 - t_1) / (t_3 - t_2));
  fprintf (stdout, "  siPC : %e %e, ratio = %f\n",
	   (t_5 - t_4),
	   (t_6 - t_5),
	   (t_5 - t_4) / (t_6 - t_5));

  check += compare_max (dt, dt_JGdP,  "dt_JGdP",  verbose, tiny, &max);
  check += compare_max (dt, dt_siPC,  "dt_siPC",  verbose, tiny, &max);
  check += compare_max (dt, dt_NJGdP, "dt_NJGdP", verbose, tiny, &max);
  check += compare_max (dt, dt_NsiPC, "dt_NsiPC", verbose, tiny, &max);

  // the error is estimated relative to the radius of particle
  for (i = 0; i < np; i ++)
    {
      int i3 = i * 3;
      char label[80];

      sprintf (label, "siPC-JGdP [%d] x", i);
      check += compare_max (x_siPC[i3]-x_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "siPC-JGdP [%d] y", i);
      check += compare_max (x_siPC[i3+1]-x_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "siPC-JGdP [%d] z", i);
      check += compare_max (x_siPC[i3+2]-x_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "NsiPC-NJGdP [%d] x", i);
      check += compare_max (x_NsiPC[i3]-x_NJGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NsiPC-NJGdP [%d] y", i);
      check += compare_max (x_NsiPC[i3+1]-x_NJGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NsiPC-NJGdP [%d] z", i);
      check += compare_max (x_NsiPC[i3+2]-x_NJGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "NJGdP-JGdP [%d] x", i);
      check += compare_max (x_NJGdP[i3]-x_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NJGdP-JGdP [%d] y", i);
      check += compare_max (x_NJGdP[i3+1]-x_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NJGdP-JGdP [%d] z", i);
      check += compare_max (x_NJGdP[i3+2]-x_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "NsiPC-siPC [%d] x", i);
      check += compare_max (x_NsiPC[i3]-x_siPC[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NsiPC-siPC [%d] y", i);
      check += compare_max (x_NsiPC[i3+1]-x_siPC[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "NsiPC-siPC [%d] z", i);
      check += compare_max (x_NsiPC[i3+2]-x_siPC[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

    }

  if (version > 0 && flag_Q != 0)
    {
      for (i = 0; i < np; i ++)
	{
	  int i4 = i * 4;
	  char label[80];

	  sprintf (label, "siPC-JGdP [%d] Q1", i);
	  check += compare_max (q_siPC[i4]-q_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "siPC-JGdP [%d] Q2", i);
	  check += compare_max (q_siPC[i4+1]-q_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "siPC-JGdP [%d] Q3", i);
	  check += compare_max (q_siPC[i4+2]-q_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "siPC-JGdP [%d] Q4", i);
	  check += compare_max (q_siPC[i4+3]-q_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "NsiPC-NJGdP [%d] Q1", i);
	  check += compare_max (q_NsiPC[i4]-q_NJGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-NJGdP [%d] Q2", i);
	  check += compare_max (q_NsiPC[i4+1]-q_NJGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-NJGdP [%d] Q3", i);
	  check += compare_max (q_NsiPC[i4+2]-q_NJGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-NJGdP [%d] Q4", i);
	  check += compare_max (q_NsiPC[i4+3]-q_NJGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "NJGdP-JGdP [%d] Q1", i);
	  check += compare_max (q_NJGdP[i4]-q_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NJGdP-JGdP [%d] Q2", i);
	  check += compare_max (q_NJGdP[i4+1]-q_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NJGdP-JGdP [%d] Q3", i);
	  check += compare_max (q_NJGdP[i4+2]-q_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NJGdP-JGdP [%d] Q4", i);
	  check += compare_max (q_NJGdP[i4+3]-q_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "NsiPC-siPC [%d] Q1", i);
	  check += compare_max (q_NsiPC[i4]-q_siPC[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-siPC [%d] Q2", i);
	  check += compare_max (q_NsiPC[i4+1]-q_siPC[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-siPC [%d] Q3", i);
	  check += compare_max (q_NsiPC[i4+2]-q_siPC[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "NsiPC-siPC [%d] Q4", i);
	  check += compare_max (q_NsiPC[i4+3]-q_siPC[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);
	}
    }


  stokes_free (sys);
  free (pos);
  free (F);
  free (T);
  free (E);
  bonds_free (bonds);
  BD_params_free (BD);
  free (x_siPC);
  free (x_NJGdP);
  free (x_NsiPC);
  free (x_JGdP);
  if (q_siPC != NULL)
    {
      free (q_siPC);
      free (q_NJGdP);
      free (q_NsiPC);
      free (q_JGdP);
    }
  BD_imp_free (BDimp);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
