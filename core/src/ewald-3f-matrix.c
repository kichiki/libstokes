/* Solvers for 3 dimensional F version problems by MATRIX procedure
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f-matrix.c,v 2.17 2007/12/01 18:42:45 kichiki Exp $
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
#include "f.h"

#include "ewald.h" // make_matrix_mob_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3f-matrix.h"


/** natural resistance problem **/
/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] : in the labo frame.
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_3f_matrix (struct stokes * sys,
		     const double *u,
		     double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int n3 = sys->np * 3;
  double *u0 = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (u0, "solve_res_3f_matrix");

  shift_labo_to_rest_U (sys, sys->np, u, u0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_3f_matrix_0 (sys, u0, f);

  free (u0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem in F version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f_matrix_0 (struct stokes * sys,
		       const double *u,
		       double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_3f_matrix_0 :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int n3 = sys->np * 3;

  double *mat = (double *) malloc (sizeof (double) * n3 * n3);
  CHECK_MALLOC (mat, "solve_res_3f_matrix");

  make_matrix_mob_3all (sys, mat); // sys->version is 0 (F)
  // f := M^-1.u
  /*
  lapack_inv_ (n3, mat);
  dot_prod_matrix (mat, n3, n3, u, f);
  */
  lapack_solve_lin (n3, mat, u, f);

  free (mat);
}


/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_lub_3f_matrix (struct stokes * sys,
			 const double *u,
			 double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_lub_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  double *u0 = (double *) malloc (sizeof (double) * sys->np * 3);
  CHECK_MALLOC (u0, "solve_res_lub_3f_matrix");

  shift_labo_to_rest_U (sys, sys->np, u, u0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_lub_3f_matrix_0 (sys, u0, f);

  free (u0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem in F version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_lub_3f_matrix_0 (struct stokes * sys,
			   const double *u,
			   double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_lub_3f_matrix_0 :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int n3 = np * 3;

  double *mob = (double *) malloc (sizeof (double) * n3 * n3);
  double *lub = (double *) malloc (sizeof (double) * n3 * n3);
  CHECK_MALLOC (mob, "solve_res_lub_3f_matrix");
  CHECK_MALLOC (lub, "solve_res_lub_3f_matrix");

  // M matrix
  make_matrix_mob_3all (sys, mob); // sys->version is 0 (F)
  // M^-1
  lapack_inv_ (n3, mob);

  // L matrix
  make_matrix_lub_3f (sys, lub);

  // M^-1 + L
  int i;
  for (i = 0; i < n3 * n3; i ++)
    {
      lub [i] += mob [i];
    }
  free (mob);

  // x := (M^-1 + L).(UO)
  dot_prod_matrix (lub, n3, n3, u, f);

  free (lub);
}


/** natural mobility problem **/
/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_3f_matrix (struct stokes * sys,
		     const double *f,
		     double *u)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mob_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  atimes_3all_matrix (sys->np * 3, f, u, (void *) sys);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
}

/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_lub_3f_matrix (struct stokes * sys,
			 const double *f,
			 double *u)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;

  int n3 = np * 3;

  double *mat = (double *) malloc (sizeof (double) * n3 * n3);
  double *lub = (double *) malloc (sizeof (double) * n3 * n3);
  double *iml = (double *) malloc (sizeof (double) * n3 * n3);
  double *x   = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (mat, "solve_mob_lub_3f_matrix");
  CHECK_MALLOC (lub, "solve_mob_lub_3f_matrix");
  CHECK_MALLOC (iml, "solve_mob_lub_3f_matrix");
  CHECK_MALLOC (x, "solve_mob_lub_3f_matrix");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  // M
  make_matrix_mob_3all (sys, mat); // sys->version is 0 (F)
  // L
  make_matrix_lub_3f (sys, lub);
  // IML := M.L
  mul_matrices (mat, n3, n3,
		lub, n3, n3,
		iml);
  free (lub);

  // IML = I + M.L
  int i;
  for (i = 0; i < n3; ++i)
    {
      iml [i * n3 + i] += 1.0;
    }

  // x := M.F
  dot_prod_matrix (mat, n3, n3, f, x);

  // u := (I+M.L)^-1.M.F
  /*
  // IML^-1
  lapack_inv_ (n3, iml);
  dot_prod_matrix (iml, n3, n3, x, u);
  */
  lapack_solve_lin (n3, iml, x, u);

  free (mat);
  free (iml);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
}

/** natural mobility problem with fixed particles **/
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat [np * 3 * np * 3] : matrix to split
 * OUTPUT
 *  mat_ll [nm * 3 * nm * 3] :
 *  mat_lh [nm * 3 * n'    ] :
 *  mat_hl [n'     * nm * 3] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 3 - nm * 3
 */
static void
split_matrix_fix_3f (int np, int nm,
		     const double * mat,
		     double * mat_ll, double * mat_lh,
		     double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i3, j3;
  int n3;
  int nl, nh;


  n3 = np * 3;
  nl = nm * 3;
  nh = n3 - nl;

  /* ll -- mobile,mobile */
  for (i = 0; i < nm; ++i)
    {
      i3 = i * 3;
      for (j = 0; j < nm; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat_ll [(i3 + ii) * nl + j3 + jj]
		    = mat [(i3 + ii) * n3 + j3 + jj];
		}
	    }
	}
    }

  /* lh -- mobile,fixed */
  for (i = 0; i < nm; ++i)
    {
      i3 = i * 3;
      for (j = nm; j < np; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat_lh [(i3 + ii) * nh + (j - nm) *3 + jj]
		    = mat [(i3 + ii) * n3 + j3 + jj];
		}
	    }
	}
    }

  /* hl -- fixed,mobile */
  for (i = nm; i < np; ++i)
    {
      i3 = i * 3;
      for (j = 0; j < nm; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat_hl [((i - nm) * 3 + ii) * nl + j3 + jj]
		    = mat [(i3 + ii) * n3 + j3 + jj];
		}
	    }
	}
    }

  /* hh -- fixed,fixed */
  for (i = nm; i < np; ++i)
    {
      i3 = i * 3;
      for (j = nm; j < np; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat_hh [((i - nm) * 3 + ii) * nh + (j - nm) * 3 + jj]
		    = mat [(i3 + ii) * n3 + j3 + jj];
		}
	    }
	}
    }
}
/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat_ll [nm * 3 * nm * 3] :
 *  mat_lh [nm * 3 * n'    ] :
 *  mat_hl [n'     * nm * 3] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 3
 * OUTPUT
 *  mat [np * 11 * np * 11] : matrix to split
 */
static void
merge_matrix_fix_3f (int np, int nm,
		     const double * mat_ll, const double * mat_lh,
		     const double * mat_hl, const double * mat_hh,
		     double * mat)
{
  int i, j;
  int ii, jj;
  int i3, j3;
  int n3;
  int nl, nh;


  n3 = np * 3;
  nl = nm * 3;
  nh = n3 - nl;

  for (i = 0; i < nm; ++i)
    {
      i3 = i * 3;
      for (j = 0; j < nm; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat [(i3 + ii) * n3 + j3 + jj]
		    = mat_ll [(i3 + ii) * nl + j3 + jj];
		}
	    }
	}
    }

  /* lh -- mobile,fixed */
  for (i = 0; i < nm; ++i)
    {
      i3 = i * 3;
      for (j = nm; j < np; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat [(i3 + ii) * n3 + j3 + jj]
		    = mat_lh [(i3 + ii) * nh + (j - nm) *3 + jj];
		}
	    }
	}
    }

  /* hl -- fixed,mobile */
  for (i = nm; i < np; ++i)
    {
      i3 = i * 3;
      for (j = 0; j < nm; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat [(i3 + ii) * n3 + j3 + jj]
		    = mat_hl [((i - nm) * 3 + ii) * nl + j3 + jj];
		}
	    }
	}
    }

  /* hh -- fixed,fixed */
  for (i = nm; i < np; ++i)
    {
      i3 = i * 3;
      for (j = nm; j < np; ++j)
	{
	  j3 = j * 3;
	  for (ii = 0; ii < 3; ++ii)
	    {
	      for (jj = 0; jj < 3; ++jj)
		{
		  mat [(i3 + ii) * n3 + j3 + jj]
		    = mat_hh [((i - nm) * 3 + ii) * nh + (j - nm) * 3 + jj];
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
solve_mix_3f_matrix (struct stokes * sys,
		     const double *f, const double *uf,
		     double *u, double *ff)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mix_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int nm = sys->nm;

  if (np == nm)
    {
      solve_mob_3f_matrix (sys, f, u);
      return;
    }

  int n3 = np * 3;
  int nf = np - nm;
  int nm3 = nm * 3;
  int nl = nm * 3;
  int nh = n3 - nl;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *mat    = (double *) malloc (sizeof (double) * n3 * n3);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n3);
  double *x = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (uf0, "solve_mix_3f_matrix");
  CHECK_MALLOC (mat, "solve_mix_3f_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_3f_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_3f_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_3f_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_3f_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_3f_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_3f_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_3f_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_3f_matrix");
  CHECK_MALLOC (b, "solve_mix_3f_matrix");
  CHECK_MALLOC (x, "solve_mix_3f_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  /* b := (F,Uf) */
  // note: set_F_by_f() == set_f_by_F() so not implemented explicitely
  set_F_by_f (nm, b, f);
  set_F_by_f (nf, b + nm3, uf0);
  free (uf0);

  /* mob */
  make_matrix_mob_3all (sys, mat); // sys->version is 0 (F)
  split_matrix_fix_3f (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (nh, nl,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  merge_matrix_fix_3f (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n3, n3,
		   b, x);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

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
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_lub_3f_matrix (struct stokes * sys,
			 const double *f, const double *uf,
			 double *u, double *ff)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3f_matrix :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int nm = sys->nm;

  if (np == nm)
    {
      solve_mob_lub_3f_matrix (sys, f, u);
      return;
    }

  int n3 = np * 3;
  int nf = np - nm;
  int nm3 = nm * 3;
  int nl = nm * 3;
  int nh = n3 - nl;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *mat = (double *) malloc (sizeof (double) * n3 * n3);
  double *lub = (double *) malloc (sizeof (double) * n3 * n3);
  double *tmp = (double *) malloc (sizeof (double) * n3 * n3);
  double *mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  double *mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  double *mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  double *mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  double *b = (double *) malloc (sizeof (double) * n3);
  double *x = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (uf0, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mat, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (lub, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (tmp, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mat_ll, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mat_lh, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mat_hl, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mat_hh, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mob_ll, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mob_lh, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mob_hl, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (mob_hh, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (b, "solve_mix_lub_3f_matrix");
  CHECK_MALLOC (x, "solve_mix_lub_3f_matrix");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  /* b := (F,Uf) */
  // note: set_F_by_f() == set_f_by_F() so not implemented explicitely
  set_F_by_f (nm, b, f);
  set_F_by_f (nf, b + nm3, uf0);
  free (uf0);

  /* mob */
  make_matrix_mob_3all (sys, mat); // sys->version is 0 (F)
  split_matrix_fix_3f (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub */
  make_matrix_lub_3f (sys, lub);
  /* tmp := M.L */
  mul_matrices (mat, n3, n3, lub, n3, n3, tmp);
  /* note: at this point, lub[] is free to use.
   * so that I_?? use lub [] */
  double *I_ll = lub;
  double *I_lh = I_ll + nl * nl;
  double *I_hl = I_lh + nl * nh;
  double *I_hh = I_hl + nh * nl;

  /* tmp := I + (M.T).(L.T) */
  int i;
  for (i = 0; i < n3; ++i)
    {
      tmp [i * n3 + i] += 1.0;
    }
  split_matrix_fix_3f (np, nm, tmp, I_ll, I_lh, I_hl, I_hh);
  /* note: at this point, tmp[] is free to use. */
  free (tmp);

  solve_gen_linear (nl, nh,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);
  /* note: at this point, I_??, therefore lub[], is free to use. */
  free (lub);

  merge_matrix_fix_3f (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n3, n3,
		   b, x);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

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
}
