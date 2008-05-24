/* test code for bd-imp.c
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd-imp.c,v 1.6 2008/05/24 06:09:35 kichiki Exp $
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
#include <brownian.h>
#include <bd-imp.h>
#include <bonds.h> // bonds_init(), bonds_free()
#include <angles.h> // angles_init(), angles_free()



int
check_BD_evolve_JGdP00 (int version, int flag_lub, int flag_mat,
			int flag_Q, double dt,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BD_evolve_JGdP00\n"
	       "(ver=%d,lub=%d,mat=%d,Q=%d,dt=%f)"
	       " : start\n",
	       version, flag_lub, flag_mat, flag_Q, dt);
    }

  int check = 0;
  double max = 0.0;


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_BD_evolve_JGdP00");

  sys->version = version;

  int np = 3;
  int nm = np;
  stokes_set_np (sys, np, nm);
  sys->periodic = 0; // non-periodic

  int np3 = np * 3;
  double *pos = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (pos, "check_BD_evolve_JGdP00");

  pos[0] = -5.0;
  pos[1] = 0.0;
  pos[2] = 0.0;

  pos[3] = 0.0;
  pos[4] = 0.0;
  pos[5] = 0.0;

  pos[6] = 7.0;
  pos[7] = 0.0;
  pos[8] = 0.0;

  stokes_set_pos (sys, pos);

  int nm3 = nm * 3;
  int nm5 = nm * 5;
  double *F = (double *)malloc (sizeof (double) * nm3);
  double *T = (double *)malloc (sizeof (double) * nm3);
  double *E = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (F, "check_BD_evolve_JGdP00");
  CHECK_MALLOC (T, "check_BD_evolve_JGdP00");
  CHECK_MALLOC (E, "check_BD_evolve_JGdP00");
  int i;
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
  CHECK_MALLOC (BD, "check_BD_evolve_JGdP00");

  double *x_mid  = (double *)malloc (sizeof (double) * np3);
  double *x_BB   = (double *)malloc (sizeof (double) * np3);
  double *x_BM   = (double *)malloc (sizeof (double) * np3);
  double *x_JGdP = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (x_mid,  "check_BD_evolve_JGdP00");
  CHECK_MALLOC (x_BB,   "check_BD_evolve_JGdP00");
  CHECK_MALLOC (x_BM,   "check_BD_evolve_JGdP00");
  CHECK_MALLOC (x_JGdP, "check_BD_evolve_JGdP00");
  for (i = 0; i < np3; i ++)
    {
      x_mid[i]  = pos[i];
      x_BB[i]   = pos[i];
      x_BM[i]   = pos[i];
      x_JGdP[i] = pos[i];
    }

  int np4 = np * 4;
  double *q_mid  = NULL;
  double *q_BB   = NULL;
  double *q_BM   = NULL;
  double *q_JGdP = NULL;
  if (sys->version > 0 && flag_Q != 0)
    {
      q_mid  = (double *)malloc (sizeof (double) * np4);
      q_BB   = (double *)malloc (sizeof (double) * np4);
      q_BM   = (double *)malloc (sizeof (double) * np4);
      q_JGdP = (double *)malloc (sizeof (double) * np4);
      CHECK_MALLOC (q_mid,  "check_BD_evolve_JGdP00");
      CHECK_MALLOC (q_BB,   "check_BD_evolve_JGdP00");
      CHECK_MALLOC (q_BM,   "check_BD_evolve_JGdP00");
      CHECK_MALLOC (q_JGdP, "check_BD_evolve_JGdP00");
      for (i = 0; i < np4; i ++)
	{
	  q_mid[i]  = 0.0;
	  q_BB[i]   = 0.0;
	  q_BM[i]   = 0.0;
	  q_JGdP[i] = 0.0;
	}
      // set Q4 = 1 as an initial condition
      for (i = 0; i < np; i ++)
	{
	  q_mid[i*4+3]  = 1.0;
	  q_BB[i*4+3]   = 1.0;
	  q_BM[i*4+3]   = 1.0;
	  q_JGdP[i*4+3] = 1.0;
	}
    }


  /**
   * BD explicit schemes
   */
  double t = 0.0;
  double dt_mid  = BD_evolve_mid    (t, BD, x_mid,  q_mid,  dt);
  double dt_BB   = BD_evolve_BB03   (t, BD, x_BB,   q_BB,   dt);
  double dt_BM   = BD_evolve_BM97   (t, BD, x_BM,   q_BM,   dt);

  /**
   * BD implicit scheme
   */
  BD->scheme = 3; // JGdP
  struct BD_imp *BDimp = BD_imp_init (BD,
				      1000,  // itmax
				      1.0e-6 // eps
				      );
  CHECK_MALLOC (BDimp, "check_BD_evolve_JGdP00");

  double dt_JGdP = BD_evolve_JGdP00 (t, BDimp, x_JGdP, q_JGdP, dt);

  check += compare_max (dt, dt_mid,  "dt_mid",  verbose, tiny, &max);
  check += compare_max (dt, dt_BB,   "dt_BB",   verbose, tiny, &max);
  check += compare_max (dt, dt_BM,   "dt_BM",   verbose, tiny, &max);
  check += compare_max (dt, dt_JGdP, "dt_JGdP", verbose, tiny, &max);

  // the error is estimated relative to the radius of particle
  for (i = 0; i < np; i ++)
    {
      int i3 = i * 3;
      char label[80];

      sprintf (label, "mid [%d] x", i);
      check += compare_max (x_mid[i3]-x_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid [%d] y", i);
      check += compare_max (x_mid[i3+1]-x_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid [%d] z", i);
      check += compare_max (x_mid[i3+2]-x_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BB [%d] x", i);
      check += compare_max (x_BB[i3]-x_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "BB [%d] y", i);
      check += compare_max (x_BB[i3+1]-x_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "BB [%d] z", i);
      check += compare_max (x_BB[i3+2]-x_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BM [%d] x", i);
      check += compare_max (x_BM[i3]-x_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "BM [%d] y", i);
      check += compare_max (x_BM[i3+1]-x_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "BM [%d] z", i);
      check += compare_max (x_BM[i3+2]-x_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      /* check */
      sprintf (label, "mid-BB [%d] x", i);
      check += compare_max (x_mid[i3]-x_BB[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid-BB [%d] y", i);
      check += compare_max (x_mid[i3+1]-x_BB[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid-BB [%d] z", i);
      check += compare_max (x_mid[i3+2]-x_BB[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "mid-BM [%d] x", i);
      check += compare_max (x_mid[i3]-x_BM[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid-BM [%d] y", i);
      check += compare_max (x_mid[i3+1]-x_BM[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);
      sprintf (label, "mid-BM [%d] z", i);
      check += compare_max (x_mid[i3+2]-x_BM[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);
    }

  if (version > 0 && flag_Q != 0)
    {
      for (i = 0; i < np; i ++)
	{
	  int i4 = i * 4;
	  char label[80];

	  sprintf (label, "mid [%d] Q1", i);
	  check += compare_max (q_mid[i4]-q_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid [%d] Q2", i);
	  check += compare_max (q_mid[i4+1]-q_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid [%d] Q3", i);
	  check += compare_max (q_mid[i4+2]-q_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid [%d] Q4", i);
	  check += compare_max (q_mid[i4+3]-q_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BB [%d] Q1", i);
	  check += compare_max (q_BB[i4]-q_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BB [%d] Q2", i);
	  check += compare_max (q_BB[i4+1]-q_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BB [%d] Q3", i);
	  check += compare_max (q_BB[i4+2]-q_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BB [%d] Q4", i);
	  check += compare_max (q_BB[i4+3]-q_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BM [%d] Q1", i);
	  check += compare_max (q_BM[i4]-q_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BM [%d] Q2", i);
	  check += compare_max (q_BM[i4+1]-q_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BM [%d] Q3", i);
	  check += compare_max (q_BM[i4+2]-q_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "BM [%d] Q4", i);
	  check += compare_max (q_BM[i4+3]-q_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  /* check */
	  sprintf (label, "mid-BB [%d] Q1", i);
	  check += compare_max (q_mid[i4]-q_BB[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BB [%d] Q2", i);
	  check += compare_max (q_mid[i4+1]-q_BB[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BB [%d] Q3", i);
	  check += compare_max (q_mid[i4+2]-q_BB[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BB [%d] Q4", i);
	  check += compare_max (q_mid[i4+3]-q_BB[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "mid-BM [%d] Q1", i);
	  check += compare_max (q_mid[i4]-q_BM[i4]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BM [%d] Q2", i);
	  check += compare_max (q_mid[i4+1]-q_BM[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BM [%d] Q3", i);
	  check += compare_max (q_mid[i4+2]-q_BM[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);
	  sprintf (label, "mid-BM [%d] Q4", i);
	  check += compare_max (q_mid[i4+3]-q_BM[i4+3]+1.0, 1.0,
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
  free (x_mid);
  free (x_BB);
  free (x_BM);
  free (x_JGdP);
  if (q_mid != NULL)
    {
      free (q_mid);
      free (q_BB);
      free (q_BM);
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


int
check_BD_imp_ode_evolve (int version, int flag_lub, int flag_mat,
			 int flag_Q, double h, double t_out,
			 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BD_imp_ode_evolve\n"
	       "(ver=%d,lub=%d,mat=%d,Q=%d,h=%f,t_out=%f)"
	       " : start\n",
	       version, flag_lub, flag_mat, flag_Q, h, t_out);
    }

  int check = 0;
  double max = 0.0;


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_BD_imp_ode_evolve");

  sys->version = version;

  int np = 3;
  int nm = np;
  stokes_set_np (sys, np, nm);
  sys->periodic = 0; // non-periodic

  int np3 = np * 3;
  double *pos = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (pos, "check_BD_imp_ode_evolve");

  pos[0] = -5.0;
  pos[1] = 0.0;
  pos[2] = 0.0;

  pos[3] = 0.0;
  pos[4] = 0.0;
  pos[5] = 0.0;

  pos[6] = 7.0;
  pos[7] = 0.0;
  pos[8] = 0.0;

  stokes_set_pos (sys, pos);

  int nm3 = nm * 3;
  int nm5 = nm * 5;
  double *F = (double *)malloc (sizeof (double) * nm3);
  double *T = (double *)malloc (sizeof (double) * nm3);
  double *E = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (F, "check_BD_imp_ode_evolve");
  CHECK_MALLOC (T, "check_BD_imp_ode_evolve");
  CHECK_MALLOC (E, "check_BD_imp_ode_evolve");
  int i;
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
  CHECK_MALLOC (BD, "check_BD_imp_ode_evolve");

  int n;
  if (version == 0)
    {
      n = np * 3;
    }
  else if (flag_Q == 0)
    {
      n = np * 3;
    }
  else
    {
      n = np * 7;
    }
  double *y_mid  = (double *)malloc (sizeof (double) * n);
  double *y_BB   = (double *)malloc (sizeof (double) * n);
  double *y_BM   = (double *)malloc (sizeof (double) * n);
  double *y_JGdP = (double *)malloc (sizeof (double) * n);
  double *y_siPC = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y_mid,  "check_BD_imp_ode_evolve");
  CHECK_MALLOC (y_BB,   "check_BD_imp_ode_evolve");
  CHECK_MALLOC (y_BM,   "check_BD_imp_ode_evolve");
  CHECK_MALLOC (y_JGdP, "check_BD_imp_ode_evolve");
  CHECK_MALLOC (y_siPC, "check_BD_imp_ode_evolve");
  for (i = 0; i < np3; i ++)
    {
      y_mid [i] = pos[i];
      y_BB  [i] = pos[i];
      y_BM  [i] = pos[i];
      y_JGdP[i] = pos[i];
      y_siPC[i] = pos[i];
    }
  if (version > 0 && flag_Q != 0)
    {
      for (i = 0; i < np * 4; i ++)
	{
	  y_mid [np3 + i] = 0.0;
	  y_BB  [np3 + i] = 0.0;
	  y_BM  [np3 + i] = 0.0;
	  y_JGdP[np3 + i] = 0.0;
	  y_siPC[np3 + i] = 0.0;
	}
      // set Q4 = 1 as an initial condition
      for (i = 0; i < np; i ++)
	{
	  y_mid [np3 + i*4+3] = 1.0;
	  y_BB  [np3 + i*4+3] = 1.0;
	  y_BM  [np3 + i*4+3] = 1.0;
	  y_JGdP[np3 + i*4+3] = 1.0;
	  y_siPC[np3 + i*4+3] = 1.0;
	}
    }



  /**
   * BD explicit schemes
   */
  double t_mid = 0.0;
  double dt_mid = h;
  BD->scheme = 0; // mid point
  BD_ode_evolve (BD, &t_mid, t_out, &dt_mid, y_mid);

  double t_BB = 0.0;
  double dt_BB = h;
  BD->scheme = 1; // Banchio-Brady
  BD_ode_evolve (BD, &t_BB,  t_out, &dt_BB,  y_BB);

  double t_BM = 0.0;
  double dt_BM = h;
  BD->scheme = 2; // Ball-Melrose
  BD_ode_evolve (BD, &t_BM,  t_out, &dt_BM,  y_BM);

  /**
   * BD implicit scheme
   */
  struct BD_imp *BDimp = NULL;
  BD->scheme = 3; // JGdP
  BDimp = BD_imp_init (BD,
		       1000,  // itmax
		       1.0e-6 // eps
		       );
  CHECK_MALLOC (BDimp, "check_BD_imp_ode_evolve");

  double t_JGdP = 0.0;
  double dt_JGdP = h;
  BD_imp_ode_evolve (BDimp, &t_JGdP, t_out, &dt_JGdP, y_JGdP);
  BD_imp_free (BDimp);


  BD->scheme = 4; // siPC
  BDimp = BD_imp_init (BD,
		       1000,  // itmax
		       1.0e-6 // eps
		       );
  CHECK_MALLOC (BDimp, "check_BD_imp_ode_evolve");

  double t_siPC = 0.0;
  double dt_siPC = h;
  BD_imp_ode_evolve (BDimp, &t_siPC, t_out, &dt_siPC, y_siPC);
  BD_imp_free (BDimp);


  check += compare_max (dt_JGdP, dt_mid, "dt_mid",  verbose, tiny, &max);
  check += compare_max (dt_JGdP, dt_BB,  "dt_BB",   verbose, tiny, &max);
  check += compare_max (dt_JGdP, dt_BM,  "dt_BM",   verbose, tiny, &max);
  check += compare_max (dt_JGdP, dt_siPC,"dt_siPC", verbose, tiny, &max);

  check += compare_max (t_out, t_mid,  "t_mid",  verbose, tiny, &max);
  check += compare_max (t_out, t_BB,   "t_BB",   verbose, tiny, &max);
  check += compare_max (t_out, t_BM,   "t_BM",   verbose, tiny, &max);
  check += compare_max (t_out, t_JGdP, "t_JGdP", verbose, tiny, &max);
  check += compare_max (t_out, t_siPC, "t_siPC", verbose, tiny, &max);

  // the error is estimated relative to the radius of particle
  for (i = 0; i < np; i ++)
    {
      int i3 = i * 3;
      char label[80];

      sprintf (label, "mid [%d] x", i);
      check += compare_max (y_mid[i3]-y_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "mid [%d] y", i);
      check += compare_max (y_mid[i3+1]-y_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "mid [%d] z", i);
      check += compare_max (y_mid[i3+2]-y_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BB [%d] x", i);
      check += compare_max (y_BB[i3]-y_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BB [%d] y", i);
      check += compare_max (y_BB[i3+1]-y_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BB [%d] z", i);
      check += compare_max (y_BB[i3+2]-y_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BM [%d] x", i);
      check += compare_max (y_BM[i3]-y_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BM [%d] y", i);
      check += compare_max (y_BM[i3+1]-y_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "BM [%d] z", i);
      check += compare_max (y_BM[i3+2]-y_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "siPC [%d] x", i);
      check += compare_max (y_siPC[i3]-y_JGdP[i3]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "siPC [%d] y", i);
      check += compare_max (y_siPC[i3+1]-y_JGdP[i3+1]+1.0, 1.0,
			    label, verbose, tiny, &max);

      sprintf (label, "siPC [%d] z", i);
      check += compare_max (y_siPC[i3+2]-y_JGdP[i3+2]+1.0, 1.0,
			    label, verbose, tiny, &max);
    }
  if (version > 0 && flag_Q != 0)
    {
      for (i = 0; i < np; i ++)
	{
	  int i4 = i * 4;
	  char label[80];

	  sprintf (label, "mid [%d] Q1", i);
	  check += compare_max (y_mid[i4]-y_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "mid [%d] Q2", i);
	  check += compare_max (y_mid[i4+1]-y_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "mid [%d] Q3", i);
	  check += compare_max (y_mid[i4+2]-y_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "mid [%d] Q4", i);
	  check += compare_max (y_mid[i4+3]-y_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BB [%d] Q1", i);
	  check += compare_max (y_BB[i4]-y_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BB [%d] Q2", i);
	  check += compare_max (y_BB[i4+1]-y_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BB [%d] Q3", i);
	  check += compare_max (y_BB[i4+2]-y_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BB [%d] Q4", i);
	  check += compare_max (y_BB[i4+3]-y_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BM [%d] Q1", i);
	  check += compare_max (y_BM[i4]-y_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BM [%d] Q2", i);
	  check += compare_max (y_BM[i4+1]-y_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BM [%d] Q3", i);
	  check += compare_max (y_BM[i4+2]-y_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "BM [%d] Q4", i);
	  check += compare_max (y_BM[i4+3]-y_JGdP[i4+3]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "siPC [%d] Q1", i);
	  check += compare_max (y_siPC[i4]-y_JGdP[i4]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "siPC [%d] Q2", i);
	  check += compare_max (y_siPC[i4+1]-y_JGdP[i4+1]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "siPC [%d] Q3", i);
	  check += compare_max (y_siPC[i4+2]-y_JGdP[i4+2]+1.0, 1.0,
				label, verbose, tiny, &max);

	  sprintf (label, "siPC [%d] Q4", i);
	  check += compare_max (y_siPC[i4+3]-y_JGdP[i4+3]+1.0, 1.0,
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
  free (y_mid);
  free (y_BB);
  free (y_BM);
  free (y_JGdP);
  free (y_siPC);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
