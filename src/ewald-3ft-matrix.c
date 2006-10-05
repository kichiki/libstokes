/* Ewald summation technique with FT version -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft-matrix.c,v 2.4 2006/10/05 04:53:14 ichiki Exp $
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

#include "dgetri_c.h" /* lapack_inv_() */

#include <libstokes.h> /* struct stokeks */
#include "/home/ichiki/WORK/SF/ryuon/libstokes/bench.h"
#include "/home/ichiki/WORK/SF/ryuon/libstokes/ft.h"

#include "ewald.h" // make_matrix_mob_ewald_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3ft-matrix.h"


/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_ewald_3ft_matrix (struct stokes * sys,
			   const double *u, const double *o,
			   double *f, double *t)
{
  int np;
  int n6;

  double * mat;
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;

  n6 = np * 6;
  mat = (double *) malloc (sizeof (double) * n6 * n6);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  /* b := (UO) */
  set_ft_by_FT (np, b, u, o);

  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 1 (FT)
  lapack_inv_ (n6, mat);

  // x := M^-1.b
  dot_prod_matrix (mat, n6, n6, b, x);

  /* x := (FT) */
  set_FT_by_ft (np, f, t, x);

  free (mat);
  free (b);
  free (x);
}


/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_lub_ewald_3ft_matrix (struct stokes * sys,
			       const double *u, const double *o,
			       double *f, double *t)
{
  int np;
  int i;
  int n6;

  double * mob;
  double * lub;
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;
  n6 = np * 6;
  mob = (double *) malloc (sizeof (double) * n6 * n6);
  lub = (double *) malloc (sizeof (double) * n6 * n6);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  // M matrix
  make_matrix_mob_ewald_3all (sys, mob); // sys->version is 1 (FT)
  // M^-1
  lapack_inv_ (n6, mob);

  // L matrix
  make_matrix_lub_ewald_3ft (sys, lub);

  // M^-1 + L
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
/* solve natural mobility problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
calc_mob_ewald_3ft_matrix (struct stokes * sys,
			   const double *f, const double *t,
			   double *u, double *o)
{
  int np;
  int n6;

  double * mat;
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;
  n6 = np * 6;
  mat = malloc (sizeof (double) * n6 * n6);
  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);

  /* b := (FTE) */
  set_ft_by_FT (np, b, f, t);

  // M
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 1 (FT)
  // x = M.b
  dot_prod_matrix (mat, n6, n6, b, x);

  set_FT_by_ft (np, u, o, x);

  free (mat);
  free (b);
  free (x);
}

/* solve natural mobility problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
calc_mob_lub_ewald_3ft_matrix (struct stokes * sys,
			       const double *f, const double *t,
			       double *u, double *o)
{
  int np;
  int i;
  int n6;

  double * mat;
  double * lub;
  double * iml;
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;

  n6 = np * 6;
  mat = (double *) malloc (sizeof (double) * n6 * n6);
  lub = (double *) malloc (sizeof (double) * n6 * n6);
  iml = (double *) malloc (sizeof (double) * n6 * n6);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  // M
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 1 (FT)
  // L
  make_matrix_lub_ewald_3ft (sys, lub);
  // IML := M.L
  mul_matrices (mat, n6, n6,
		lub, n6, n6,
		iml);
  free (lub);
  // IML = I + M.L
  for (i = 0; i < n6; ++i)
    {
      iml [i * n6 + i] += 1.0;
    }
  // IML^-1
  lapack_inv_ (n6, iml);

  /* b := (FT) */
  set_ft_by_FT (np, b, f, t);
  // x := M.(FT)
  dot_prod_matrix (mat, n6, n6, b, x);

  // b := (I+M.L)^-1.M.(FT)
  dot_prod_matrix (iml, n6, n6, x, b);

  set_FT_by_ft (np, u, o, b);

  free (mat);
  free (iml);
  free (b);
  free (x);
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
		  mat_hh [((i - nm) * 6 + ii) * nl + (j - nm) * 6 + jj]
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
		    = mat_hh [((i - nm) * 6 + ii) * nl + (j - nm) * 6 + jj];
		}
	    }
	}
    }
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
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
calc_mob_fix_ewald_3ft_matrix (struct stokes * sys,
			       const double *f, const double *t,
			       const double *uf, const double *of,
			       double *u, double *o,
			       double *ff, double *tf)
{
  int np, nm;
  int n6;
  int nf, nm6;
  int nl, nh;

  double * mat;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;
  nm = sys->nm;

  n6 = np * 6;
  nf = np - nm;
  nm6 = nm * 6;
  nl = nm * 6;
  nh = n6 - nl;

  mat    = (double *) malloc (sizeof (double) * n6 * n6);
  mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  /* b := (FT,UfOf) */
  set_ft_by_FT (nm, b, f, t);
  set_ft_by_FT (nf, b + nm6, uf, of);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 1 (FT)
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
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
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
calc_mob_lub_fix_ewald_3ft_matrix (struct stokes * sys,
				   const double *f, const double *t,
				   const double *uf, const double *of,
				   double *u, double *o,
				   double *ff, double *tf)
{
  int np, nm;

  int i;
  int n6;
  int nf, nm6;
  int nl, nh;

  double * mat;
  double * lub;
  double * tmp;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * I_ll, * I_lh, * I_hl, * I_hh; /* used at lub [] */
  double * b;
  double * x;


  sys->version = 1; // FT version
  np = sys->np;
  nm = sys->nm;

  n6 = np * 6;
  nf = np - nm;
  nm6 = nm * 6;
  nl = nm * 6;
  nh = n6 - nl;

  mat = (double *) malloc (sizeof (double) * n6 * n6);
  lub = (double *) malloc (sizeof (double) * n6 * n6);
  tmp = (double *) malloc (sizeof (double) * n6 * n6);
  mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  /* used at lub [] */
  I_ll = lub;
  I_lh = I_ll + nl * nl;
  I_hl = I_lh + nl * nh;
  I_hh = I_hl + nh * nl;

  /* b := (FT,UfOf) */
  set_ft_by_FT (nm, b, f, t);
  set_ft_by_FT (nf, b + nm6, uf, of);

  /* mob */
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 1 (FT)
  split_matrix_fix_3ft (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub */
  make_matrix_lub_ewald_3ft (sys, lub);
  /* tmp := M.L */
  mul_matrices (mat, n6, n6, lub, n6, n6, tmp);
  /* note: at this point, lub[] is free to use. */
  free (lub);
  /* lub := I + (M.T).(L.T) */
  for (i = 0; i < n6; ++i)
    {
      tmp [i * n6 + i] += 1.0;
    }
  split_matrix_fix_3ft (np, nm, tmp, I_ll, I_lh, I_hl, I_hh);
  /* note: at this point, tmp[] is free to use. */
  free (tmp);

  solve_gen_linear (nl, nh,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);

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
}
