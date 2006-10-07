/* Ewald summation technique with F version -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f-matrix.c,v 2.7 2006/10/07 01:04:41 kichiki Exp $
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

#include "stokes.h"
#include "bench.h"
#include "f.h"

#include "ewald.h" // make_matrix_mob_ewald_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3f-matrix.h"


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_ewald_3f_matrix (struct stokes * sys,
			  const double *u,
			  double *f)
{
  int n3;
  double * mat;


  sys->version = 0; // F version
  n3 = sys->np * 3;
  mat = malloc (sizeof (double) * n3 * n3);

  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 0 (F)
  lapack_inv_ (n3, mat);

  // f := M^-1.u
  dot_prod_matrix (mat, n3, n3, u, f);

  free (mat);
}


/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *u,
			      double *f)
{
  int np;
  int i;
  int n3;

  double * mob;
  double * lub;


  sys->version = 0; // F version
  np = sys->np;
  n3 = np * 3;
  mob = malloc (sizeof (double) * n3 * n3);
  lub = malloc (sizeof (double) * n3 * n3);

  // M matrix
  make_matrix_mob_ewald_3all (sys, mob); // sys->version is 0 (F)
  // M^-1
  lapack_inv_ (n3, mob);

  // L matrix
  make_matrix_lub_ewald_3f (sys, lub);

  // M^-1 + L
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
/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_ewald_3f_matrix (struct stokes * sys,
			  const double *f,
			  double *u)
{
  sys->version = 0; // F version
  atimes_ewald_3all_matrix (sys->np * 3, f, u, (void *) sys);
}

/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
calc_mob_lub_ewald_3f_matrix (struct stokes * sys,
			      const double *f,
			      double *u)
{
  int np;
  int i;
  int n3;

  double * mat;
  double * lub;
  double * iml;
  double * x;


  sys->version = 0; // F version
  np = sys->np;

  n3 = np * 6;
  mat = malloc (sizeof (double) * n3 * n3);
  lub = malloc (sizeof (double) * n3 * n3);
  iml = malloc (sizeof (double) * n3 * n3);
  x = malloc (sizeof (double) * n3);

  // M
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 0 (F)
  // L
  make_matrix_lub_ewald_3f (sys, lub);
  // IML := M.L
  mul_matrices (mat, n3, n3,
		lub, n3, n3,
		iml);
  free (lub);
  // IML = I + M.L
  for (i = 0; i < n3; ++i)
    {
      iml [i * n3 + i] += 1.0;
    }
  // IML^-1
  lapack_inv_ (n3, iml);

  // x := M.(FT)
  dot_prod_matrix (mat, n3, n3, f, x);

  // u := (I+M.L)^-1.M.(FT)
  dot_prod_matrix (iml, n3, n3, x, u);

  free (mat);
  free (iml);
  free (x);
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
calc_mob_fix_ewald_3f_matrix (struct stokes * sys,
			      const double *f, const double *uf,
			      double *u, double *ff)
{
  int np, nm;
  int n3;
  int nf, nm3;
  int nl, nh;

  double * mat;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * b;
  double * x;


  sys->version = 0; // F version
  np = sys->np;
  nm = sys->nm;

  if (np == nm)
    {
      calc_mob_ewald_3f_matrix (sys, f, u);
      return;
    }

  n3 = np * 3;
  nf = np - nm;
  nm3 = nm * 3;
  nl = nm * 3;
  nh = n3 - nl;

  mat    = (double *) malloc (sizeof (double) * n3 * n3);
  mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  b = (double *) malloc (sizeof (double) * n3);
  x = (double *) malloc (sizeof (double) * n3);

  /* b := (F,Uf) */
  // note: set_F_by_f() == set_f_by_F() so not implemented explicitely
  set_F_by_f (nm, b, f);
  set_F_by_f (nf, b + nm3, uf);

  /* mob */
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 0 (F)
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
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_3f_matrix (struct stokes * sys,
				  const double *f, const double *uf,
				  double *u, double *ff)
{
  int np, nm;

  int i;
  int n3;
  int nf, nm3;
  int nl, nh;

  double * mat;
  double * lub;
  double * tmp;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * I_ll, * I_lh, * I_hl, * I_hh; /* used at lub [] */
  double * b;
  double * x;


  sys->version = 0; // F version
  np = sys->np;
  nm = sys->nm;

  n3 = np * 3;
  nf = np - nm;
  nm3 = nm * 3;
  nl = nm * 3;
  nh = n3 - nl;

  mat = (double *) malloc (sizeof (double) * n3 * n3);
  lub = (double *) malloc (sizeof (double) * n3 * n3);
  tmp = (double *) malloc (sizeof (double) * n3 * n3);
  mat_ll = (double *) malloc (sizeof (double) * nl * nl);
  mat_lh = (double *) malloc (sizeof (double) * nl * nh);
  mat_hl = (double *) malloc (sizeof (double) * nh * nl);
  mat_hh = (double *) malloc (sizeof (double) * nh * nh);
  mob_ll = (double *) malloc (sizeof (double) * nl * nl);
  mob_lh = (double *) malloc (sizeof (double) * nl * nh);
  mob_hl = (double *) malloc (sizeof (double) * nh * nl);
  mob_hh = (double *) malloc (sizeof (double) * nh * nh);
  b = (double *) malloc (sizeof (double) * n3);
  x = (double *) malloc (sizeof (double) * n3);

  /* used at lub [] */
  I_ll = lub;
  I_lh = I_ll + nl * nl;
  I_hl = I_lh + nl * nh;
  I_hh = I_hl + nh * nl;

  /* b := (F,Uf) */
  // note: set_F_by_f() == set_f_by_F() so not implemented explicitely
  set_F_by_f (nm, b, f);
  set_F_by_f (nf, b + nm3, uf);

  /* mob */
  make_matrix_mob_ewald_3all (sys, mat); // sys->version is 0 (F)
  split_matrix_fix_3f (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub */
  make_matrix_lub_ewald_3f (sys, lub);
  /* tmp := M.L */
  mul_matrices (mat, n3, n3, lub, n3, n3, tmp);
  /* note: at this point, lub[] is free to use. */
  free (lub);
  /* lub := I + (M.T).(L.T) */
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
}
