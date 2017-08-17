/* test code for bd-imp-fast.c
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include <gsl/gsl_poly.h> // gsl_poly_solve_cubic()

#include "memory-check.h" // macro CHECK_MALLOC

#include <libstokes-core.h> // ptime_ms_d()

#include <bonds.h>    // struct BONDS
#include <brownian.h> // struct BD_params
#include <bd-imp.h>   // struct BD_imp

#include "check.h" // compare_max()

#include <bd-imp-fast.h>


struct BD_imp *
BD_imp_init_for_test (int np,
		      int type,
		      int flag_noHI,
		      double radius,
		      double peclet,
		      double dt)
{
  // struct stokes
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "BD_imp_init_for_test");
  sys->version = 0; // F version
  stokes_set_np (sys, np, np);
  sys->periodic = 0; // non-periodic
  //sys->rmin = radius;

  double *a = (double *)calloc (np, sizeof (double));
  CHECK_MALLOC (a, "BD_imp_init_for_test");
  int i;
  for (i = 0; i < np; i ++)
    {
      a[i] = radius;
    }
  stokes_set_radius (sys, a);
  free (a);

  // for flag_mat == 0
  stokes_set_iter (sys, "otmk", 2000, 20, 1.0e-6, 0, NULL);


  // struct BONDS
  struct BONDS *bonds = BONDS_init ();
  CHECK_MALLOC (bonds, "BD_imp_init_for_test");

  // straight configuration
  double A  = 1.0;
  double Q0 = 1.0;
  double s  = 1.0e-2;
  for (i = 0; i < np - 1; i ++)
    {
      BONDS_append (bonds,
		    type,
		    0, // fene, in this case, (p1,p2) = (A,Q0)
		    A,
		    Q0,
		    s,
		    i, i + 1);
    }


  // struct BD_params
  int np3 = np * 3;
  int np5 = np * 5;
  double *F = (double *)malloc (sizeof (double) * np3);
  double *T = (double *)malloc (sizeof (double) * np3);
  double *E = (double *)malloc (sizeof (double) * np5);
  CHECK_MALLOC (F, "BD_imp_init_for_test");
  CHECK_MALLOC (T, "BD_imp_init_for_test");
  CHECK_MALLOC (E, "BD_imp_init_for_test");
  for (i = 0; i < np; i ++)
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

  struct BD_params *BD
    = BD_params_init (sys,
		      0, // unsigned long seed
		      F, T, E,
		      NULL, NULL, NULL, // uf, of, ef,
		      flag_noHI,
		      0, // int flag_lub,
		      0, // int flag_mat,
		      0.0, // st,
		      NULL, // struct BeadRod *br,
		      bonds,
		      1.0, // gamma,
		      NULL, // ev,
		      NULL, // ang
		      NULL, // ev_dh,
		      NULL, // ev_LJ,
		      NULL, // cf,
		      0, // int flag_Q,
		      peclet,
		      1.0e-6, // ode_eps,
		      0, // n_minv, // n of chebyshev for minv
		      0, // n_lub,  // n of chebyshev for lub
		      0, // BD_scheme, (0=mid, 1=BB, 2=BM, 3=JGdP, 4=siPC)
		      100, // BB_n,
		      1.0,  // double rmin
		      1.0e-12 //dt_lim
		      );
  CHECK_MALLOC (BD, "BD_imp_init_for_test");
  BD->scheme = 3; // JGdP

  // struct BD_imp
  struct BD_imp *b = BD_imp_init (BD,
				  2,       // fastSI,
				  1000,    // itmax
				  1.0e-6); // eps
  CHECK_MALLOC (b, "BD_imp_init_for_test");
  for (i = 0; i < np3; i ++)
    {
      b->z[i] = 0.0;
    }
  BD_imp_set_dt (b, dt);

  return (b);
}

void
BD_imp_free_for_test (struct BD_imp *b)
{
  if (b == NULL) return;

  if (b->BD->sys != NULL) stokes_free (b->BD->sys);
  if (b->BD->bonds != NULL) BONDS_free (b->BD->bonds);
  if (b->BD->F != NULL) free (b->BD->F);
  if (b->BD->T != NULL) free (b->BD->T);
  if (b->BD->E != NULL) free (b->BD->E);
  if (b->BD != NULL) BD_params_free (b->BD);
  BD_imp_free (b);
}


/*
 * INPUT
 *  solver : 0 == GSL multiroot
 *           1 == NITSOL
 *           2 == fastSI
 */
struct BD_imp *
BD_imp_reset_solver (struct BD_imp *b,
		     int solver,
		     int itmax, double eps)
{
  if (b == NULL)
    {
      fprintf (stderr, "# BD_imp_reset_solver:"
	       " b should be prepared before calling\n");
      exit (1);
    }

  // keep the necessary informations
  struct BD_params *BD = b->BD;
  double dt = b->dt;

  BD_imp_free (b);

  // re-init
  b = BD_imp_init (BD,
		   solver,
		   itmax,
		   eps);
  int i;
  for (i = 0; i < BD->sys->np * 3; i ++)
    {
      b->z[i] = 0.0;
    }
  // set dt
  BD_imp_set_dt (b, dt);

  return (b);
  // b is allocated so that it should be returned
}


void
set_q (struct BD_imp *b,
       double *q,
       double dq,
       long seed)
{
  srand48(seed);

  // loop for groups
  int ig0 = 0; // head pointer of the group "ig"
  int ng = b->BD->groups->n;
  int ig;
  for (ig = 0; ig < ng; ig ++)
    {
      struct BONDS_GROUP *g = b->BD->groups->group[ig];

      if (g->np == 1)
	{
	  // isolated particle
	  int iq = ig0;
	  int iq3 = iq * 3;
	  q[iq3  ] = drand48();
	  q[iq3+1] = drand48();
	  q[iq3+2] = drand48();
	}
      else // if (g->np > 1)
	{
	  // loop for bonds in group "ig"
	  int i;
	  for (i = 0; i < g->np - 1; i ++)
	    {
	      int ib = g->bonds[i];
	      int iq = ig0 + i;
	      int iq3 = iq * 3;

	      double Q;
	      int type  = b->BD->bonds->type[ib];
	      double Q0 = b->BD->bonds->p2[ib];
	      if (type == 0 ||
		  type == 6)
		{
		  Q = Q0 * (1.0 + dq * 2.0 * (drand48() - 0.5));
		}
	      else if (type == 7)
		{
		  double s = b->BD->bonds->p3[ib];
		  Q = Q0 * (1.0 + s * dq * 2.0 * (drand48() - 0.5));
		}
	      else
		{
		  Q = Q0 * dq * drand48();
		}

	      if (i % 3 == 0)
		{
		  q[iq3  ] = Q;
		  q[iq3+1] = 0.0;
		  q[iq3+2] = 0.0;
		}
	      else if (i % 3 == 1)
		{
		  q[iq3  ] = 0.0;
		  q[iq3+1] = Q;
		  q[iq3+2] = 0.0;
		}
	      else
		{
		  q[iq3  ] = 0.0;
		  q[iq3+1] = 0.0;
		  q[iq3+2] = Q;
		}
	    }
	  // the COM component
	  int iCOM  = ig0 + i;
	  int iCOM3 = iCOM * 3;
	  q[iCOM3  ] = drand48();
	  q[iCOM3+1] = drand48();
	  q[iCOM3+2] = drand48();
	}
      ig0 += g->np;
    }
}



/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_rhs (int type, int np, int flag_noHI,
		  int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_fastSI_rhs: type = %d, np = %d, noHI = %d\n",
	       type, np, flag_noHI);
    }

  int check = 0;
  double max = 0.0;

  double radius = 0.1;
  double peclet = 1.0;
  double dt = 1.0e-2;
  struct BD_imp *b
    = BD_imp_init_for_test (np, type, flag_noHI,
			    radius, peclet, dt);
  CHECK_MALLOC (b, "check_fastSI_rhs");

  int np3 = np * 3;
  double *q0   = (double *)calloc (np3, sizeof (double));
  double *rhs  = (double *)calloc (np3, sizeof (double));
  double *f    = (double *)calloc (np3, sizeof (double));
  double *lhs  = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (q0,   "check_fastSI_rhs");
  CHECK_MALLOC (rhs,  "check_fastSI_rhs");
  CHECK_MALLOC (f,    "check_fastSI_rhs");
  CHECK_MALLOC (lhs,  "check_fastSI_rhs");


  set_q (b, q0, 0.8, 0);
  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     q0,
		     b->BD->sys->pos);

  int ns = np - 1;
  int i;
  for (i = 0; i < ns; i ++)
    {
      int ix = i * 3;
      fastSI_rhs (b, q0, q0, q0, q0,
		  i, // ib
		  i, // index for connector (divided by 3)
		  rhs  + ix);
      // here use q0[] for all connector vectors

      // calc lhs
      double fs[3];
      BONDS_calc_force_spring_i (b->BD->bonds, i, q0 + ix, fs);

      double B;
      if (b->BD->sys->a == NULL)
	{
	  B = 2.0;
	}
      else
	{
	  int ia = b->BD->bonds->ia[i];
	  int ib = b->BD->bonds->ib[i];
	  B = 1.0 / b->BD->sys->a[ia]
	    + 1.0 / b->BD->sys->a[ib];
	}
      lhs[ix  ] = q0[ix  ] + dt * B * fs[0];
      lhs[ix+1] = q0[ix+1] + dt * B * fs[1];
      lhs[ix+2] = q0[ix+2] + dt * B * fs[2];
    }

  fastSI_set_Q0 (b, q0); // for fastSI_f()
  fastSI_f (b, q0, f);


  char label[80];
  for (i = 0; i < ns * 3; i ++)
    {
      sprintf (label, " f vs rhs - lhs [%d]", i);
      check += compare_max (1.0 + f[i],
			    1.0 + rhs[i] - lhs[i],
			    label,  verbose, tiny, &max);
    }

  BD_imp_free_for_test (b);
  free (q0);
  free (rhs);
  free (f);
  free (lhs);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_solve_cubic (int type, int np, int flag_noHI,
			  int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_fastSI_solve_cubic: type = %d, np = %d, noHI = %d\n",
	       type, np, flag_noHI);
    }

  int check = 0;
  double max = 0.0;

  double radius = 0.1;
  double peclet = 1.0;
  double dt = 0.01;
  struct BD_imp *b = BD_imp_init_for_test (np, type, flag_noHI,
					   radius, peclet, dt);
  CHECK_MALLOC (b, "check_fastSI_solve");

  double *q0 = (double *)calloc (np*3, sizeof (double));
  CHECK_MALLOC (q0, "check_fastSI_solve");

  set_q (b, q0, 0.8, 0);
  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     q0,
		     b->BD->sys->pos);

  char label[80];
  int i;
  for (i = 0; i < np - 1; i ++)
    {
      double rhs[3];
      fastSI_rhs (b, q0, q0, q0, q0,
		  i, // ib
		  i, // index for connector (divided by 3)
		  rhs);
      double r_norm = sqrt (rhs[0] * rhs[0]
			    + rhs[1] * rhs[1]
			    + rhs[2] * rhs[2]);
      double q_norm = fastSI_solve_cubic (b, i, r_norm);
      double q[3];
      q[0] = q_norm * rhs[0] / r_norm;
      q[1] = q_norm * rhs[1] / r_norm;
      q[2] = q_norm * rhs[2] / r_norm;

      double f[3];
      BONDS_calc_force_spring_i (b->BD->bonds, i, q, f);

      double B;
      if (b->BD->sys->a == NULL)
	{
	  B = 2.0;
	}
      else
	{
	  int ia = b->BD->bonds->ia[i];
	  int ib = b->BD->bonds->ib[i];
	  B = 1.0 / b->BD->sys->a[ia]
	    + 1.0 / b->BD->sys->a[ib];
	}
      double lhs[3];
      lhs[0] = q[0] + b->dt * B * f[0];
      lhs[1] = q[1] + b->dt * B * f[1];
      lhs[2] = q[2] + b->dt * B * f[2];

      sprintf (label, " type=%d [%d] x", type, i);
      check += compare_max (lhs[0], rhs[0], label,  verbose, tiny, &max);

      sprintf (label, " type=%d [%d] y", type, i);
      check += compare_max (lhs[1], rhs[1], label,  verbose, tiny, &max);

      sprintf (label, " type=%d [%d] z", type, i);
      check += compare_max (lhs[2], rhs[2], label,  verbose, tiny, &max);
    }

  BD_imp_free_for_test (b);
  free (q0);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

/* 
 * INPUT
 *  type : the following type (except for ILC)
 *         0 == Hookean (Fraenkel)
 *         1 == WLC
 *         2 == ILC
 *         3 == Cohen
 *         4 == Werner
 *         5 == Hookean
 *         6 == another Hookean (Fraenkel)
 *         7 == FENE-Fraenkel
 *  ns   : number of spring
 */
int
check_fastSI_solve (int type, int np, int flag_noHI, double dt,
		    int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_fastSI_solve: type = %d, np = %d, noHI = %d, dt = %e\n",
	       type, np, flag_noHI, dt);
    }

  int check = 0;
  double max = 0.0;

  double radius = 0.1;
  double peclet = 1.0;
  struct BD_imp *b
    = BD_imp_init_for_test (np, type, flag_noHI,
			    radius, peclet, dt);
  CHECK_MALLOC (b, "check_fastSI_solve");

  int np3 = np * 3;
  double *q0 = (double *)calloc (np3, sizeof (double));
  double *q1 = (double *)calloc (np3, sizeof (double));
  double *q2 = (double *)calloc (np3, sizeof (double));
  double *q3 = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (q0, "check_fastSI_solve");
  CHECK_MALLOC (q1, "check_fastSI_solve");
  CHECK_MALLOC (q2, "check_fastSI_solve");
  CHECK_MALLOC (q3, "check_fastSI_solve");

  set_q (b, q0, 0.8, 0);
  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     q0,
		     b->BD->sys->pos);

  // fastSI
  b = BD_imp_reset_solver (b,
			   2,       // fastSI
			   1000,    // itmax
			   1.0e-8); // eps
  double t_fastSI_0 = ptime_ms_d ();
  fastSI_solve (b, q0, q1);
  double t_fastSI_1 = ptime_ms_d ();

  // NITSOL
  b = BD_imp_reset_solver (b,
			   1,       // NITSOL
			   1000,    // itmax
			   1.0e-8); // eps
  NITSOL_set_iplvl (b->nit,
		    1,  // iplvl
		    6); // stdout
  b->nit->f = fastSI_NITSOL_f;
  double t_NITSOL_0 = ptime_ms_d ();
  fastSI_NITSOL_wrap (b, q0, q2);
  double t_NITSOL_1 = ptime_ms_d ();

  // GSL
  b = BD_imp_reset_solver (b,
			   0,       // GSL multiroot
			   1000,    // itmax
			   1.0e-8); // eps
  b->F->f = fastSI_GSL_MULTIROOT_func;
  double t_GSL_0 = ptime_ms_d ();
  fastSI_GSL_MULTIROOT_wrap (b, q0, q3);
  double t_GSL_1 = ptime_ms_d ();

  double *f1 = (double *)calloc (np3, sizeof (double));
  double *f2 = (double *)calloc (np3, sizeof (double));
  double *f3 = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (f1, "check_fastSI_solve");
  CHECK_MALLOC (f2, "check_fastSI_solve");
  CHECK_MALLOC (f3, "check_fastSI_solve");

  fastSI_set_Q0 (b, q0);
  fastSI_f (b, q1, f1);
  fastSI_f (b, q2, f2);
  fastSI_f (b, q3, f3);

  char label[80];
  int i;
  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " NITSOL vs GSL q[%d]", i);
      check += compare_max (1.0+q2[i], 1.0+q3[i], label,  verbose, tiny, &max);
    }
  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " fastSI vs NITSOL q[%d]", i);
      check += compare_max (1.0+q1[i], 1.0+q2[i], label,  verbose, tiny, &max);
    }
  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " fastSI vs GSL q[%d]", i);
      check += compare_max (1.0+q1[i], 1.0+q3[i], label,  verbose, tiny, &max);
    }

  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " f for fastSI [%d]", i);
      check += compare_max (f1[i]+1.0, 1.0, label,  verbose, tiny, &max);
    }
  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " f for NITSOL [%d]", i);
      check += compare_max (f2[i]+1.0, 1.0, label,  verbose, tiny, &max);
    }
  for (i = 0; i < np3; i ++)
    {
      sprintf (label, " f for GSL [%d]", i);
      check += compare_max (f3[i]+1.0, 1.0, label,  verbose, tiny, &max);
    }

  BD_imp_free_for_test (b);
  free (q0);
  free (q1);
  free (q2);
  free (q3);
  free (f1);
  free (f2);
  free (f3);

  if (verbose != 0)
    {
      // benchmark
      fprintf (stdout, " CPU times (fastSI) %f\n", t_fastSI_1 - t_fastSI_0);
      fprintf (stdout, " CPU times (NITSOL) %f\n", t_NITSOL_1 - t_NITSOL_0);
      fprintf (stdout, " CPU times (GSL)    %f\n", t_GSL_1    - t_GSL_0);
      fprintf (stdout, " speed of fastSI: %f (NITSOL)  %f (GSL)\n",
	       (t_NITSOL_1 - t_NITSOL_0) / (t_fastSI_1 - t_fastSI_0),
	       (t_GSL_1    - t_GSL_0   ) / (t_fastSI_1 - t_fastSI_0));
      fprintf (stdout, "\n");

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
