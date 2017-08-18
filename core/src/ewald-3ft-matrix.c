/* Solvers for 3 dimensional FT version problems by MATRIX procedure
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft-matrix.c,v 2.16 2007/12/01 18:43:28 kichiki Exp $
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
#include "memory-check.h" // CHECK_MALLOC

#include "dgetri_c.h" /* lapack_inv_() */

#include "stokes.h"
#include "bench.h"
#include "ft.h"
#include "f.h"

#include "ewald.h" // make_matrix_mob_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3ft-matrix.h"


/** natural resistance problem **/
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : in the labo frame.
 *  o [np * 3] : in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft_matrix (struct stokes * sys,
		      const double *u, const double *o,
		      double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  double *u0 = (double *) malloc (sizeof (double) * np * 3);
  double *o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_3ft_matrix");
  CHECK_MALLOC (o0, "solve_res_3ft_matrix");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_3ft_matrix_0 (sys,
			  u0, o0, 
			  f, t);

  free (u0);
  free (o0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem in FT version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft_matrix_0 (struct stokes * sys,
			const double *u, const double *o,
			double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_3ft_matrix_0 :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *mat = (double *) malloc (sizeof (double) * n6 * n6);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (mat, "solve_res_3ft_matrix");
  CHECK_MALLOC (b, "solve_res_3ft_matrix");
  CHECK_MALLOC (x, "solve_res_3ft_matrix");

  /* b := (UO) */
  set_ft_by_FT (np, b, u, o);

  make_matrix_mob_3all (sys, mat); // sys->version is 1 (FT)
  // x := M^-1.b
  /*
  lapack_inv_ (n6, mat);
  dot_prod_matrix (mat, n6, n6, b, x);
  */
  lapack_solve_lin (n6, mat, b, x);

  /* x := (FT) */
  set_FT_by_ft (np, f, t, x);

  free (mat);
  free (b);
  free (x);
}


/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : in the labo frame.
 *  o [np * 3] : in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_lub_3ft_matrix (struct stokes * sys,
			  const double *u, const double *o,
			  double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_lub_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  double *u0 = (double *) malloc (sizeof (double) * np * 3);
  double *o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_lub_3ft_matrix");
  CHECK_MALLOC (o0, "solve_res_lub_3ft_matrix");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_lub_3ft_matrix_0 (sys, u0, o0,
			      f, t);

  free (u0);
  free (o0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem in FT version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_lub_3ft_matrix_0 (struct stokes * sys,
			    const double *u, const double *o,
			    double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_lub_3ft_matrix_0 :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *mob = (double *) malloc (sizeof (double) * n6 * n6);
  double *lub = (double *) malloc (sizeof (double) * n6 * n6);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (mob, "solve_res_lub_3ft_matrix");
  CHECK_MALLOC (lub, "solve_res_lub_3ft_matrix");
  CHECK_MALLOC (b, "solve_res_lub_3ft_matrix");
  CHECK_MALLOC (x, "solve_res_lub_3ft_matrix");

  // M matrix
  make_matrix_mob_3all (sys, mob); // sys->version is 1 (FT)
  // M^-1
  lapack_inv_ (n6, mob);

  // L matrix
  make_matrix_lub_3ft (sys, lub);

  // M^-1 + L
  int i;
  for (i = 0; i < n6 * n6; i ++)
    {
      lub [i] += mob [i];
    }
  free (mob);

  /* b := (UO) */
  set_ft_by_FT (np, b, u, o);

  // x := (M^-1 + L).(UO)
  dot_prod_matrix (lub, n6, n6, b, x);

  // (FT) = x
  set_FT_by_ft (np, f, t, x);

  free (lub);
  free (b);
  free (x);
}


/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft_matrix (struct stokes * sys,
		      const double *f, const double *t,
		      double *u, double *o)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mob_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *mat = (double *)malloc (sizeof (double) * n6 * n6);
  double *b = (double *)malloc (sizeof (double) * n6);
  double *x = (double *)malloc (sizeof (double) * n6);
  CHECK_MALLOC (mat, "solve_mob_3ft_matrix");
  CHECK_MALLOC (b, "solve_mob_3ft_matrix");
  CHECK_MALLOC (x, "solve_mob_3ft_matrix");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  /* b := (FTE) */
  set_ft_by_FT (np, b, f, t);

  // M
  make_matrix_mob_3all (sys, mat); // sys->version is 1 (FT)
  // x = M.b
  dot_prod_matrix (mat, n6, n6, b, x);

  set_FT_by_ft (np, u, o, x);

  free (mat);
  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}

/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft_matrix (struct stokes * sys,
			  const double *f, const double *t,
			  double *u, double *o)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *mat = (double *) malloc (sizeof (double) * n6 * n6);
  double *lub = (double *) malloc (sizeof (double) * n6 * n6);
  double *iml = (double *) malloc (sizeof (double) * n6 * n6);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (mat, "solve_mob_lub_3ft_matrix");
  CHECK_MALLOC (lub, "solve_mob_lub_3ft_matrix");
  CHECK_MALLOC (iml, "solve_mob_lub_3ft_matrix");
  CHECK_MALLOC (b, "solve_mob_lub_3ft_matrix");
  CHECK_MALLOC (x, "solve_mob_lub_3ft_matrix");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  // M
  make_matrix_mob_3all (sys, mat); // sys->version is 1 (FT)
  // L
  make_matrix_lub_3ft (sys, lub);
  // IML := M.L
  mul_matrices (mat, n6, n6,
		lub, n6, n6,
		iml);
  free (lub);

  /* b := (FT) */
  set_ft_by_FT (np, b, f, t);
  // x := M.(FT)
  dot_prod_matrix (mat, n6, n6, b, x);

  // IML = I + M.L
  int i;
  for (i = 0; i < n6; ++i)
    {
      iml [i * n6 + i] += 1.0;
    }
  // b := (I+M.L)^-1.M.(FT)
  /*
  // IML^-1
  lapack_inv_ (n6, iml);
  dot_prod_matrix (iml, n6, n6, x, b);
  */
  lapack_solve_lin (n6, iml, x, b);

  set_FT_by_ft (np, u, o, b);

  free (mat);
  free (iml);
  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}

/** natural mobility problem with fixed particles **/
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat [np * 6 * np * 6] : matrix to split
 * OUTPUT
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 6 - nm * 6
 */
static void
split_matrix_fix_3ft (int np, int nm,
		      const double * mat,
		      double * mat_ll, double * mat_lh,
		      double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i6, j6;
  int n6;
  int nl, nh;


  n6 = np * 6;
  nl = nm * 6;
  nh = n6 - nl;

  /* ll -- mobile,mobile */
  for (i = 0; i < nm; ++i)
    {
      i6 = i * 6;
      for (j = 0; j < nm; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_ll [(i6 + ii) * nl + j6 + jj]
		    = mat [(i6 + ii) * n6 + j6 + jj];
		}
	    }
	}
    }

  /* lh -- mobile,fixed */
  for (i = 0; i < nm; ++i)
    {
      i6 = i * 6;
      for (j = nm; j < np; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_lh [(i6 + ii) * nh + (j - nm) *6 + jj]
		    = mat [(i6 + ii) * n6 + j6 + jj];
		}
	    }
	}
    }

  /* hl -- fixed,mobile */
  for (i = nm; i < np; ++i)
    {
      i6 = i * 6;
      for (j = 0; j < nm; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [((i - nm) * 6 + ii) * nl + j6 + jj]
		    = mat [(i6 + ii) * n6 + j6 + jj];
		}
	    }
	}
    }

  /* hh -- fixed,fixed */
  for (i = nm; i < np; ++i)
    {
      i6 = i * 6;
      for (j = nm; j < np; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hh [((i - nm) * 6 + ii) * nh + (j - nm) * 6 + jj]
		    = mat [(i6 + ii) * n6 + j6 + jj];
		}
	    }
	}
    }
}
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 6
 * OUTPUT
 *  mat [np * 11 * np * 11] : matrix to split
 */
static void
merge_matrix_fix_3ft (int np, int nm,
		      const double * mat_ll, const double * mat_lh,
		      const double * mat_hl, const double * mat_hh,
		      double * mat)
{
  int i, j;
  int ii, jj;
  int i6, j6;
  int n6;
  int nl, nh;


  n6 = np * 6;
  nl = nm * 6;
  nh = n6 - nl;

  for (i = 0; i < nm; ++i)
    {
      i6 = i * 6;
      for (j = 0; j < nm; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i6 + ii) * n6 + j6 + jj]
		    = mat_ll [(i6 + ii) * nl + j6 + jj];
		}
	    }
	}
    }

  /* lh -- mobile,fixed */
  for (i = 0; i < nm; ++i)
    {
      i6 = i * 6;
      for (j = nm; j < np; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i6 + ii) * n6 + j6 + jj]
		    = mat_lh [(i6 + ii) * nh + (j - nm) *6 + jj];
		}
	    }
	}
    }

  /* hl -- fixed,mobile */
  for (i = nm; i < np; ++i)
    {
      i6 = i * 6;
      for (j = 0; j < nm; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i6 + ii) * n6 + j6 + jj]
		    = mat_hl [((i - nm) * 6 + ii) * nl + j6 + jj];
		}
	    }
	}
    }

  /* hh -- fixed,fixed */
  for (i = nm; i < np; ++i)
    {
      i6 = i * 6;
      for (j = nm; j < np; ++j)
	{
	  j6 = j * 6;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i6 + ii) * n6 + j6 + jj]
		    = mat_hh [((i - nm) * 6 + ii) * nh + (j - nm) * 6 + jj];
		}
	    }
	}
    }
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_3ft_matrix (struct stokes * sys,
		      const double *f, const double *t,
		      const double *uf, const double *of,
		      double *u, double *o,
		      double *ff, double *tf)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mix_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3ft_matrix (sys, f, t,
			    u, o);
      return;
    }

  int n6 = np * 6;
  int nf = np - nm;
  int nm6 = nm * 6;
  int nl = nm * 6;
  int nh = n6 - nl;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *mat    = (double *) malloc (sizeof (double) * n6 * n6);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (uf0, "solve_mix_3ft_matrix");
  CHECK_MALLOC (of0, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mat, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_3ft_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_3ft_matrix");
  CHECK_MALLOC (b, "solve_mix_3ft_matrix");
  CHECK_MALLOC (x, "solve_mix_3ft_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  /* b := (FT,UfOf) */
  set_ft_by_FT (nm, b, f, t);
  //set_ft_by_FT (nf, b + nm6, uf, of);
  set_ft_by_FT (nf, b + nm6, uf0, of0);
  free (uf0);
  free (of0);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_3all (sys, mat); // sys->version is 1 (FT)
  split_matrix_fix_3ft (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (nh, nl,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  merge_matrix_fix_3ft (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n6, n6,
		   b, x);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_lub_3ft_matrix (struct stokes * sys,
			  const double *f, const double *t,
			  const double *uf, const double *of,
			  double *u, double *o,
			  double *ff, double *tf)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3ft_matrix :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3ft_matrix (sys, f, t,
				u, o);
      return;
    }

  int n6 = np * 6;
  int nf = np - nm;
  int nm6 = nm * 6;
  int nl = nm * 6;
  int nh = n6 - nl;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *mat = (double *) malloc (sizeof (double) * n6 * n6);
  double *lub = (double *) malloc (sizeof (double) * n6 * n6);
  double *tmp = (double *) malloc (sizeof (double) * n6 * n6);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (uf0, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (of0, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mat, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (lub, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (tmp, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (b, "solve_mix_lub_3ft_matrix");
  CHECK_MALLOC (x, "solve_mix_lub_3ft_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  /* b := (FT,UfOf) */
  set_ft_by_FT (nm, b, f, t);
  set_ft_by_FT (nf, b + nm6, uf0, of0);
  free (uf0);
  free (of0);

  /* mob */
  make_matrix_mob_3all (sys, mat); // sys->version is 1 (FT)
  split_matrix_fix_3ft (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub */
  make_matrix_lub_3ft (sys, lub);
  /* tmp := M.L */
  mul_matrices (mat, n6, n6, lub, n6, n6, tmp);
  /* note: at this point, lub[] is free to use.
   * so that I_?? use lub [] */
  double *I_ll = lub;
  double *I_lh = I_ll + nl * nl;
  double *I_hl = I_lh + nl * nh;
  double *I_hh = I_hl + nh * nl;

  /* tmp := I + (M.T).(L.T) */
  int i;
  for (i = 0; i < n6; ++i)
    {
      tmp [i * n6 + i] += 1.0;
    }
  split_matrix_fix_3ft (np, nm, tmp, I_ll, I_lh, I_hl, I_hh);
  free (tmp);

  solve_gen_linear (nl, nh,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);
  /* note: at this point, I_??, therefore lub[], is free to use. */
  free (lub);

  merge_matrix_fix_3ft (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n6, n6,
		   b, x);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (mat);
  free (mat_ll);
  free (mat_lh);
  free (mat_hl);
  free (mat_hh);
  free (mob_ll);
  free (mob_lh);
  free (mob_hl);
  free (mob_hh);
  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
