/* Ewald summation technique under FTS version -- MATRIX procedure
 * Copyright (C) 1993-1996,1999-2001 Kengo Ichiki
 *   <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3fts-matrix.c,v 1.5 2001/02/15 04:57:13 ichiki Exp $
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */

#ifdef ZETA
#include "../bench.h" // ptime_ms_d()
#endif /* ZETA */

#include "../cholesky.h" /* cholesky() */
#include "../ludcmp.h" /* ludcmp() */

#include "fts.h"
#include "matrix.h"
#include "ewald-3fts-matrix.h"


/* (local) global variable */
/*int NUMB_mobile_particles;*/ /* this is dirty, though ... */


/** function prototypes for local routines **/
static void
make_matrix_mob_ewald_3fts (int np, double * mat);
static void
make_matrix_lub_ewald_3fts (int np, double * mat);
static int
cond_lub (double * x1, double * x2);
static void
trans_ext (int np, double *r);
static void
trans_mat_ext2ext (int np, double * mat);
static void
multiply_extmat_with_extvec_3fts (int np, double * m, double * x,
				  double * y);
static void
split_matrix_3fts (int np, double *mat,
		   double * mat_ll, double * mat_lh,
		   double * mat_hl, double * mat_hh);
static void
merge_matrix_3fts (int np,
		   double * mat_ll, double * mat_lh,
		   double * mat_hl, double * mat_hh,
		   double *mat);
static void
split_matrix_fix_3fts (int np, int nm,
		       double * mat,
		       double * mat_ll, double * mat_lh,
		       double * mat_hl, double * mat_hh);
static void
merge_matrix_fix_3fts (int np, int nm,
		       double * mat_ll, double * mat_lh,
		       double * mat_hl, double * mat_hh,
		       double * mat);

static void
test_symmetric (int n, double * mat, double tiny);
static void
multiply_matrices (int n, double *a, double *b);


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FTS version
 * INPUT
 *  (global) pos [] : position of particles
 *  n := np * 11
 *  x [n * 11] : FTS
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3fts_matrix (int n, double *x, double *y)
{
  int np;
  double * mat;


  np = n / 11;
  mat = malloc (sizeof (double) * n * n);
  if (mat == NULL)
    {
      fprintf (stderr, "allocation error in atimes_ewald_3fts_matrix ().\n");
      exit (1);
    }

  make_matrix_mob_ewald_3fts (np, mat);
  multiply_extmat_with_extvec_3fts (np, mat, x, y);

  free (mat);
}

/** natural resistance problem **/
/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_ewald_3fts_matrix (int np,
			    double *u, double *o, double *e,
			    double *f, double *t, double *s)
{
  int n11;

  double * mat;
  double * b;
  double * x;


  n11 = np * 11;
  mat = malloc (sizeof (double) * n11 * n11);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (mat == NULL
      || b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* b := (UOE) */
  set_fts_by_FTS (np, b, u, o, e);

  make_matrix_mob_ewald_3fts (np, mat); // mobility matrix in EXTRACTED form
  /* for test */
  test_symmetric (n11, mat, 1.0e-12);

  /* resistance matrix in INVERSED form */
  /*cholesky (mat, n11);*/
  lu_inv (mat, n11);
  trans_ext (np, mat); // resistance matrix in EXTRACTED form

  multiply_extmat_with_extvec_3fts (np, mat, b, x);

  set_FTS_by_fts (np, f, t, s, x);

  free (mat);
  free (b);
  free (x);
}

/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_lub_ewald_3fts_matrix (int np,
				double *u, double *o, double *e,
				double *f, double *t, double *s)
{
  int i;
  int n11;

  double * mat;
  double * b;
  double * x;
  double * y;


  n11 = np * 11;
  mat = malloc (sizeof (double) * n11 * n11);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  y = malloc (sizeof (double) * n11);
  if (mat == NULL
      || b == NULL
      || x == NULL
      || y == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* b := (UOE) */
  set_fts_by_FTS (np, b, u, o, e);

  make_matrix_lub_ewald_3fts (np, mat); // lub matrix in EXTRACTED form
  multiply_extmat_with_extvec_3fts (np, mat, b, x); // x := L.(UOE)

  make_matrix_mob_ewald_3fts (np, mat); // mobility matrix in EXTRACTED form
  multiply_extmat_with_extvec_3fts (np, mat, x, y); // y := M.L.(UOE)

  /* y := (I + M.L).(UOE) */
  for (i = 0; i < n11; ++i)
    b [i] += y [i];

  /* resistance matrix in INVERSED form */
  /*cholesky (mat, n11);*/
  lu_inv (mat, n11);
  trans_ext (np, mat); // resistance matrix in EXTRACTED form

  /* x := (M^-1).(I + M.L).(UOE) */
  multiply_extmat_with_extvec_3fts (np, mat, b, x);

  set_FTS_by_fts (np, f, t, s, x);

  free (mat);
  free (b);
  free (x);
  free (y);
}

/** natural mobility problem **/
/* solve natural mobility problem in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_3fts_matrix (int np,
			    double *f, double *t, double *e,
			    double *u, double *o, double *s)
{
  int n11, n6, n5;

  double * mat;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * b;
  double * x;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;
  mat = malloc (sizeof (double) * n11 * n11);
  mat_ll = malloc (sizeof (double) * n6 * n6);
  mat_lh = malloc (sizeof (double) * n6 * n5);
  mat_hl = malloc (sizeof (double) * n5 * n6);
  mat_hh = malloc (sizeof (double) * n5 * n5);
  mob_ll = malloc (sizeof (double) * n6 * n6);
  mob_lh = malloc (sizeof (double) * n6 * n5);
  mob_hl = malloc (sizeof (double) * n5 * n6);
  mob_hh = malloc (sizeof (double) * n5 * n5);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (mat == NULL
      || mat_ll == NULL
      || mat_lh == NULL
      || mat_hl == NULL
      || mat_hh == NULL
      || mob_ll == NULL
      || mob_lh == NULL
      || mob_hl == NULL
      || mob_hh == NULL
      || b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* b := (FTE) */
  set_fts_by_FTS (np, b, f, t, e);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_ewald_3fts (np, mat);
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_3fts (np, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (n5, n6,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  /* STEP 6 */
  merge_matrix_3fts (np, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (np, u, o, s, x);

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

/* solve natural mobility problem in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_lub_ewald_3fts_matrix (int np,
				double *f, double *t, double *e,
				double *u, double *o, double *s)
{
  int i;
  int n11, n6, n5;

  double * mat;
  double * lub;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * I_ll, * I_lh, * I_hl, * I_hh; /* used at lub [] */
  double * b;
  double * x;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;
  mat = malloc (sizeof (double) * n11 * n11);
  lub = malloc (sizeof (double) * n11 * n11);
  mat_ll = malloc (sizeof (double) * n6 * n6);
  mat_lh = malloc (sizeof (double) * n6 * n5);
  mat_hl = malloc (sizeof (double) * n5 * n6);
  mat_hh = malloc (sizeof (double) * n5 * n5);
  mob_ll = malloc (sizeof (double) * n6 * n6);
  mob_lh = malloc (sizeof (double) * n6 * n5);
  mob_hl = malloc (sizeof (double) * n5 * n6);
  mob_hh = malloc (sizeof (double) * n5 * n5);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (mat == NULL
      || lub == NULL
      || mat_ll == NULL
      || mat_lh == NULL
      || mat_hl == NULL
      || mat_hh == NULL
      || mob_ll == NULL
      || mob_lh == NULL
      || mob_hl == NULL
      || mob_hh == NULL
      || b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* used at lub [] */
  I_ll = lub;
  I_lh = I_ll + n6 * n6;
  I_hl = I_lh + n6 * n5;
  I_hh = I_hl + n5 * n6;

  /* b := (FTE) */
  set_fts_by_FTS (np, b, f, t, e);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_ewald_3fts (np, mat);
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_3fts (np, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub matrix in EXTRACTED form */
  make_matrix_lub_ewald_3fts (np, lub);
  /* lub := L.T, where T.(UOE) = (UOE~) */
  trans_mat_ext2ext (np, lub);
  /* lub := (M.T).(L.T) */
  multiply_matrices (n11, mat, lub);
  /* lub := I + (M.T).(L.T) */
  for (i = 0; i < n11; ++i)
    mat [i * n11 + i] += 1.0;
  /* note: at this point, lub[] is free to use. */
  split_matrix_3fts (np, mat, I_ll, I_lh, I_hl, I_hh);

  solve_gen_linear (n6, n5,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);

  /* STEP 6 */
  merge_matrix_3fts (np, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (np, u, o, s, x);

  free (mat);
  free (lub);
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

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_fix_ewald_3fts_matrix (int np, int nm,
				double *f, double *t, double *e,
				double *uf, double *of, double *ef,
				double *u, double *o, double *s,
				double *ff, double *tf, double *sf)
{
  int n11;
  int nf, nm11;
  int nl, nh;

  double * mat;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * b;
  double * x;


  n11 = np * 11;
  nf = np - nm;
  nm11 = nm * 11;
  nl = nm * 6;
  nh = n11 - nl;

  mat = malloc (sizeof (double) * n11 * n11);
  mat_ll = malloc (sizeof (double) * nl * nl);
  mat_lh = malloc (sizeof (double) * nl * nh);
  mat_hl = malloc (sizeof (double) * nh * nl);
  mat_hh = malloc (sizeof (double) * nh * nh);
  mob_ll = malloc (sizeof (double) * nl * nl);
  mob_lh = malloc (sizeof (double) * nl * nh);
  mob_hl = malloc (sizeof (double) * nh * nl);
  mob_hh = malloc (sizeof (double) * nh * nh);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (mat == NULL
      || mat_ll == NULL
      || mat_lh == NULL
      || mat_hl == NULL
      || mat_hh == NULL
      || mob_ll == NULL
      || mob_lh == NULL
      || mob_hl == NULL
      || mob_hh == NULL
      || b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* b := (FTE) */
  set_fts_by_FTS (nm, b, f, t, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_ewald_3fts (np, mat);
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_fix_3fts (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  solve_linear (nh, nl,
		mat_hh, mat_hl, mat_lh, mat_ll,
		mob_hh, mob_hl, mob_lh, mob_ll);

  /* STEP 6 */
  merge_matrix_fix_3fts (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

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
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_lub_fix_ewald_3fts_matrix (int np, int nm,
				    double *f, double *t, double *e,
				    double *uf, double *of, double *ef,
				    double *u, double *o, double *s,
				    double *ff, double *tf, double *sf)
{
  int i;
  int n11;
  int nf, nm11;
  int nl, nh;

  double * mat;
  double * lub;
  double * mat_ll, * mat_lh, * mat_hl, * mat_hh;
  double * mob_ll, * mob_lh, * mob_hl, * mob_hh;
  double * I_ll, * I_lh, * I_hl, * I_hh; /* used at lub [] */
  double * b;
  double * x;


  n11 = np * 11;
  nf = np - nm;
  nm11 = nm * 11;
  nl = nm * 6;
  nh = n11 - nl;

  mat = malloc (sizeof (double) * n11 * n11);
  lub = malloc (sizeof (double) * n11 * n11);
  mat_ll = malloc (sizeof (double) * nl * nl);
  mat_lh = malloc (sizeof (double) * nl * nh);
  mat_hl = malloc (sizeof (double) * nh * nl);
  mat_hh = malloc (sizeof (double) * nh * nh);
  mob_ll = malloc (sizeof (double) * nl * nl);
  mob_lh = malloc (sizeof (double) * nl * nh);
  mob_hl = malloc (sizeof (double) * nh * nl);
  mob_hh = malloc (sizeof (double) * nh * nh);
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (mat == NULL
      || lub == NULL
      || mat_ll == NULL
      || mat_lh == NULL
      || mat_hl == NULL
      || mat_hh == NULL
      || mob_ll == NULL
      || mob_lh == NULL
      || mob_hl == NULL
      || mob_hh == NULL
      || b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts_matrix().\n");
      exit (1);
    }

  /* used at lub [] */
  I_ll = lub;
  I_lh = I_ll + nl * nl;
  I_hl = I_lh + nl * nh;
  I_hh = I_hl + nh * nl;

  /* b := (FTE) */
  set_fts_by_FTS (nm, b, f, t, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* mobility matrix in EXTRACTED form */
  make_matrix_mob_ewald_3fts (np, mat);
  /* mat := M.T, where T.(FTS) = (FTS~) */
  trans_mat_ext2ext (np, mat);
  split_matrix_fix_3fts (np, nm, mat, mat_ll, mat_lh, mat_hl, mat_hh);

  /* lub matrix in EXTRACTED form */
  make_matrix_lub_ewald_3fts (np, lub);
  /* lub := L.T, where T.(UOE) = (UOE~) */
  trans_mat_ext2ext (np, lub);
  /* lub := (M.T).(L.T) */
  multiply_matrices (n11, mat, lub);
  /* lub := I + (M.T).(L.T) */
  for (i = 0; i < n11; ++i)
    mat [i * n11 + i] += 1.0;
  /* note: at this point, lub[] is free to use. */
  split_matrix_fix_3fts (np, nm, mat, I_ll, I_lh, I_hl, I_hh);

  solve_gen_linear (nl, nh,
		    I_ll, I_lh, I_hl, I_hh,
		    mat_ll, mat_lh, mat_hl, mat_hh,
		    mob_ll, mob_lh, mob_hl, mob_hh);

  /* STEP 6 */
  merge_matrix_fix_3fts (np, nm, mob_ll, mob_lh, mob_hl, mob_hh, mat);
  dot_prod_matrix (mat, n11, n11,
		   b, x);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (mat);
  free (lub);
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


/* make ewald-summed mobility matrix for FTS version
 * INPUT
 *  (global) pos [] : position of particles
 *  np : # particles
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
static void
make_matrix_mob_ewald_3fts (int np, double * mat)
{
  extern int pcellx, pcelly, pcellz;
  extern int kmaxx, kmaxy, kmaxz;

  extern double zeta, zeta2, zaspi, za2;
  extern double pi2;
  extern double pivol;
  extern double lx, ly, lz; /* cell size */

  extern double *pos;

#ifdef ZETA
  extern double cpu1, cpu2, cpu3;
#endif /* ZETA */

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double kexp;

  double erfczr;
  double expzr2;

  double a2, c2;


  n = np * 11;

  /* clear result */
  for (i = 0; i < n * n; ++i)
    mat [i] = 0.0;

  /* diagonal part ( self part ) */
  xa = ya = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
  xc = yc = 0.75 - zaspi * za2 * 10.0;
  xm = ym = zm = 0.9 - zaspi * za2 * (12.0 - 30.24 * za2);

  for (i = 0; i < np; i++)
    {
      matrix_fts_ij (i, i,
		     0.0, 0.0, 0.0,
		     xa, ya,
		     0.0,
		     xc, yc,
		     0.0, 0.0,
		     0.0,
		     xm, ym, zm,
		     n, mat);
    }

#ifdef ZETA
  ptime_ms_d ();
#endif /* ZETA */

  /* first Ewald part ( real space ) */
  for (i = 0; i < np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < np; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  for (m1 = - pcellx; m1 <= pcellx; m1++)
	    {
	      rlx = lx * (double) m1;
	      for (m2 = - pcelly; m2 <= pcelly; m2++)
		{
		  rly = ly * (double) m2;
		  for (m3 = - pcellz; m3 <= pcellz; m3++)
		    {
		      rlz = lz * (double) m3;
  
		      xx = pos [jx] - pos [ix] + rlx;
		      yy = pos [jy] - pos [iy] + rly;
		      zz = pos [jz] - pos [iz] + rlz;
		      rr = sqrt (xx * xx + yy * yy + zz * zz);

		      if (rr > 0.0)
			{
			  zr = zeta * rr;
			  zr2 = zr * zr;
			  s  = rr;
			  s2 = s * s;

			  erfczr = erfc (zr);
			  expzr2 = zaspi * exp (- zr2);

			  ex = xx / rr;
			  ey = yy / rr;
			  ez = zz / rr;

			  ya = (0.75 + 0.5 / s2) / s * erfczr
			    + ((1.0 + zr2 *
				(14.0 + 4.0 * zr2 *
				 (- 5.0 + zr2))) / s2
			       - 4.5 + 3.0 * zr2)
			    * expzr2;
			  a2 = (0.75 - 1.5 / s2) / s * erfczr
			    + ((- 3.0 + zr2 *
				(- 2.0 + 4.0 * zr2 *
				 (4.0 - zr2))) / s2
			       + 1.5 - 3.0 * zr2)
			    * expzr2;
			  xa = a2 + ya;
	      
			  yb = - 0.75 / s2 * erfczr
			    - 1.5 * (+ 1.0 + zr2 *
				     (- 6.0 + zr2 *
				      (+ 2.0)))
			    / s * expzr2;

			  yc = - 3.0 / 8.0 / s2 / s * erfczr
			    - 0.75 * (+ 1.0 + zr2 *
				      (+ 14.0 + zr2 *
				       (-20.0 + zr2 *
					( + 4.0))))
			    / s2 * expzr2;
			  c2 = 9.0 / 8.0 / s2 / s * erfczr
			    - 0.75 * (- 3.0 + zr2 *
				      (- 2.0 + zr2 *
				       (+ 16.0 + zr2 *
					(- 4.0))))
			    / s2 * expzr2;
			  xc = c2 + yc;
	      
			  xg = (2.25 - 3.6 / s2) / s2 * erfczr
			    + (- 1.5 * (- 3.0 + zr2 *
					(+ 6.0))
			       - 0.8 * (+ 9.0 + zr2 *
					(+ 6.0 + zr2 *
					 (- 48.0 + zr2 *
					  (+ 12.0)))) / s2)
			    / s * expzr2;
			  yg = 1.2 / s2 / s2 * erfczr
			    + (- 3.0 * ( zr2 *
					 (2.0 + zr2 *
					  (- 1.0)))
			       - 0.8 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s * expzr2;

			  yh = - 9.0 / 8.0 / s2 / s * erfczr
			    + 1.5 * (- 1.5 + zr2 *
				     (- 1.0 + zr2 *
				      (+ 8.0 + zr2 *
				       (- 2.0))))
			    / s2 * expzr2;

			  xm = (- 4.5 + 10.8 / s2) / s / s2 * erfczr
			    + (+ 1.5 * (- 6.0 +  zr2 *
					(- 12.0 + zr2 *
					 (+ 12.0)))
			       + 1.2 * (+ 18.0 + zr2 *
					(+ 12.0 + zr2 *
					 (+ 30.0 + zr2 *
					  (- 66.0 + zr2 *
					   (+ 12.0))))) / s2)
			    / s2 * expzr2;
			  ym = (+ 2.25 - 7.2 / s2) / s / s2 * erfczr
			    + (- 1.5 * (- 3.0 +  zr2 *
					(+ 6.0 + zr2 *
					 (- 12.0 + zr2 *
					  (+ 4.0))))
			       - 1.2 * (+ 12.0 + zr2 *
					(+ 8.0 + zr2 *
					 (- 22.0 + zr2 *
					  (+ 58.0 + zr2 *
					   (- 34.0 + zr2 *
					    (+ 4.0)))))) / s2)
			    / s2 * expzr2;
			  zm = + 1.8 / s2 / s / s2 * erfczr
			    + (- 1.5 * (+ 0.0 +  zr2 *
					(+ 8.0 + zr2 *
					 (- 4.0)))
			       - 1.2 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s2 * expzr2;
	      
			  matrix_fts_ij (j, i,
					 ex, ey, ez,
					 xa, ya,
					 yb,
					 xc, yc,
					 xg, yg,
					 yh,
					 xm, ym, zm,
					 n, mat);
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  cpu2 = ptime_ms_d ();
#endif /* ZETA */

  /* Second Ewald part ( reciprocal space ) */
  for (m1 = - kmaxx; m1 <= kmaxx; m1++)
    {
      k1 = pi2 * (double) m1 / lx;
      for (m2 = - kmaxy; m2 <= kmaxy; m2++)
	{
	  k2 = pi2 * (double) m2 / ly;
	  for (m3 = - kmaxz; m3 <= kmaxz; m3++)
	    {
	      k3 = pi2 * (double) m3 / lz;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  kk = k1 * k1 + k2 * k2 + k3 * k3;
		  k = sqrt (kk);
		  k4z = kk / 4.0 / zeta2;
		  kexp = pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ex = k1 / k;
		  ey = k2 / k;
		  ez = k3 / k;

		  ya = 6.0 * (1.0 - kk / 3.0) * kexp;
		  yb = 3.0 * k * kexp;
		  yc = 3.0 / 2.0 * kk * kexp;
		  yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
		  yh = 3.0 / 2.0 * kk * kexp;
		  ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
      
		  for (i = 0; i < np; i++)
		    {
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      for (j = 0; j < np; j++)
			{
			  jx = j * 3;
			  jy = jx + 1;
			  jz = jx + 2;

			  xx = pos [jx] - pos [ix];
			  yy = pos [jy] - pos [iy];
			  zz = pos [jz] - pos [iz];

			  cf = cos (+ k1 * xx
				    + k2 * yy
				    + k3 * zz);

			  sf = - sin (+ k1 * xx
				      + k2 * yy
				      + k3 * zz);

			  matrix_fts_ij (j, i,
					 ex, ey, ez,
					 0.0, cf * ya,
					 sf * yb,
					 0.0, cf * yc,
					 0.0, sf * yg,
					 cf * yh,
					 0.0, cf * ym, 0.0,
					 n, mat);
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  cpu3 = ptime_ms_d ();
  cpu1 = cpu2 + cpu3;
#endif /* ZETA */
}


/* make lubrication matrix for FTS version for all particles
 * under the periodic boundary condition
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
static void
make_matrix_lub_ewald_3fts (int np, double * mat)
{
  extern double * pos;
  extern double llx [27], lly [27], llz [27];

  int i, j, k;
  int i3, i11;
  int j3, j11;
  int n;

  double * tmp_pos;


  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3fts().\n");
      exit (1);
    }

  n = np * 11;

  /* clear result */
  for (i = 0; i < n * n; ++i)
    mat [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j11 = j * 11;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = pos [j3 + 0] + llx [k];
	      tmp_pos [1] = pos [j3 + 1] + lly [k];
	      tmp_pos [2] = pos [j3 + 2] + llz [k];
	      if (cond_lub (pos + i3, tmp_pos) == 0)
		{
		  matrix_lub_fts_2b (i, j,
				     pos + i3, tmp_pos,
				     n, mat);
		}
	    }
	}
    }

  free (tmp_pos);
}

/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub (double * x1, double * x2)
{
  double x, y, z;
  double r2;


  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  r2 = x * x
    + y * y
    + z * z;

  if (r2 != 0.0
      && r2 < 9.0) // r = 3.0 is the critical separation for lubrication now.
    return 0;
  else
    return 1;
}


/** copy from NR/src/FINITE/stokes-fts.c Rev 1.4 **/
/*
 * INPUT
 *  r [np * 11 * np * 11] : this is INVERSED form
 * OUTPUT
 *  r [np * 11 * np * 11] : this is EXTRACTED form
 */
static void
trans_ext (int np, double *r)
{
  static double
    tinv [121] = {
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  0.0,     0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  2.0/3.0, 0.0, 0.0, 0.0, -1.0/3.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.5, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.5, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,     0.0, 0.0, 0.5,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/3.0, 0.0, 0.0, 0.0, 2.0/3.0};

  int i, ii, i11, i0;
  int j, jj, j11, j0;
  int k;
  int n;

  double * tmp;

  n = np * 11;

  tmp = malloc (sizeof (double) * n * n);
  if (tmp == NULL)
    {
      fprintf (stderr, "allocation error in trans_ext().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  tmp [i0 * n + j0] = 0.0;
		  for (k=0; k<11; k++)
		    {
		      tmp [i0 * n + j0] +=
			r [i0 * n + j11 + k]
			* tinv [k * 11 + jj];
		    }
		}
	    }
	}
    }

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  r [i0 * n + j0] = 0.0;
		  for (k = 0; k < 11; ++k)
		    {
		      r [i0 * n + j0] +=
			tinv [ii * 11 + k]
			* tmp [(i11 + k) * n + j0];
		    }
		}
	    }
	}
    }

  free (tmp);
}

/* multiply transformation matrix from right-hand-side,
 * so that this could be multiplied by extracted vector
 * and return the extracted vector
 * INPUT
 *  mat [np * 11 * np * 11] :
 * OUTPUT
 *  mat [np * 11 * np * 11] := mat . t, where t.E = E~
 */
static void
trans_mat_ext2ext (int np, double * mat)
{
  static double
    t [121] = {
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  0.0,  0.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  2.0,  0.0, 0.0, 0.0,  1.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  2.0, 0.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 2.0, 0.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 2.0,  0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  1.0,  0.0, 0.0, 0.0,  2.0};

  int i, ii, i11, i0;
  int j, jj, j11, j0;
  int k;
  int n;

  double * tmp;

  n = np * 11;

  tmp = malloc (sizeof (double) * n * n);
  if (tmp == NULL)
    {
      fprintf (stderr, "allocation error in trans_ext().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      i0 = i11 + ii;
	      for (jj = 0; jj < 11; ++jj)
		{
		  j0 = j11 + jj;
		  tmp [i0 * n + j0] = 0.0;
		  for (k=0; k<11; k++)
		    {
		      tmp [i0 * n + j0] +=
			mat [i0 * n + j11 + k]
			* t [k * 11 + jj];
		    }
		}
	    }
	}
    }

  for (i = 0; i < n * n; ++i)
    {
      mat [i] = tmp [i];
    }

  free (tmp);
}

/** copy from test-fts-atimes.c Rev 1.5 **/
/* utility routine for matrix in the extracted form
 * INPUT
 *  np : # particles (not # elements!)
 *  m [np *11 * np *11] : matrix in the extracted form
 *  x [np *11] : vector in the extracted form
 * INPUT
 *  y [np *11] : output vector in the extracted form (:= m.x)
 */
static void
multiply_extmat_with_extvec_3fts (int np, double * m, double * x,
				  double * y)
{
  int n11;
  int i;
  int j;
  int j11;
  int jj;


  n11 = np * 11;

  for (i = 0; i < n11; ++i)
    {
      y [i] = 0.0;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (jj = 0; jj < 6; ++jj)
	    {
	      y [i] += m [i * n11 + j11 + jj] * x [j11 + jj];
	    }
	  y [i] += m [i * n11 + j11 + 6]
	    * (2.0 * x [j11 + 6] + x [j11 + 10]);
	  y [i] += m [i * n11 + j11 + 7] * 2.0 * x [j11 + 7];
	  y [i] += m [i * n11 + j11 + 8] * 2.0 * x [j11 + 8];
	  y [i] += m [i * n11 + j11 + 9] * 2.0 * x [j11 + 9];
	  y [i] += m [i * n11 + j11 + 10]
	    * (2.0 * x [j11 + 10] + x [j11 + 6]);
	}
    }
}

static void
split_matrix_3fts (int np, double *mat,
		   double * mat_ll, double * mat_lh,
		   double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11, n6, n5;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_ll [(i6 + ii) * n6 + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_lh [(i6 + ii) * n5 + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(i5 + ii) * n6 + j6 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(i5 + ii) * n5 + j5 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }
}

static void
merge_matrix_3fts (int np,
		   double * mat_ll, double * mat_lh,
		   double * mat_hl, double * mat_hh,
		   double *mat)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11, n6, n5;


  n11 = np * 11;
  n6 = np * 6;
  n5 = np * 5;

  for (i = 0; i < np; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_ll [(i6 + ii) * n6 + j6 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_lh [(i6 + ii) * n5 + j5 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hl [(i5 + ii) * n6 + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(i5 + ii) * n5 + j5 + jj];
		}
	    }
	}
    }
}

/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat [np * 11 * np * 11] : matrix to split
 * OUTPUT
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 6
 */
static void
split_matrix_fix_3fts (int np, int nm,
		       double * mat,
		       double * mat_ll, double * mat_lh,
		       double * mat_hl, double * mat_hh)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11;
  int nl, nh;
  int nm5;


  n11 = np * 11;
  nl = nm * 6;
  nh = n11 - nl;
  nm5 = nm * 5;

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_ll [(i6 + ii) * nl + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_lh [(i6 + ii) * nh + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(i5 + ii) * nl + j6 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(i5 + ii) * nh + j5 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* lh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_lh [(i6 + ii) * nh + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_hh [(i5 + ii) * nh + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + 6 + ii) * n11 + j11 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat_hl [(nm5 + (i - nm) * 11 + ii) * nl + j6 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat_hh [(nm5 + (i - nm) * 11 + ii) * nh + j5 + jj]
		    = mat [(i11 + ii) * n11 + j11 + 6 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat_hh [(nm5 + (i - nm) * 11 + ii) * nh
			 + nm5 + (j - nm) * 11 + jj]
		    = mat [(i11 + ii) * n11 + j11 + jj];
		}
	    }
	}
    }
}

/*
 * INPUT
 *  np : # ALL particles
 *  nm : # MOBILE particles
 *  mat [np * 11 * np * 11] : matrix to split
 * OUTPUT
 *  mat_ll [nm * 6 * nm * 6] :
 *  mat_lh [nm * 6 * n'    ] :
 *  mat_hl [n'     * nm * 6] :
 *  mat_hh [n'     * n'    ] :
 *  where n' = np * 11 - nm * 6
 */
static void
merge_matrix_fix_3fts (int np, int nm,
		       double * mat_ll, double * mat_lh,
		       double * mat_hl, double * mat_hh,
		       double * mat)
{
  int i, j;
  int ii, jj;
  int i11, i6, i5;
  int j11, j6, j5;
  int n11;
  int nl, nh;
  int nm5;


  n11 = np * 11;
  nl = nm * 6;
  nh = n11 - nl;
  nm5 = nm * 5;

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* ll */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_ll [(i6 + ii) * nl + j6 + jj];
		}
	      /* lh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_lh [(i6 + ii) * nh + j5 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hl [(i5 + ii) * nl + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(i5 + ii) * nh + j5 + jj];
		}
	    }
	}
    }

  for (i = 0; i < nm; ++i)
    {
      i11 = i * 11;
      i6 = i * 6;
      i5 = i * 5;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 6; ++ii)
	    {
	      /* lh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_lh [(i6 + ii) * nh + nm5 + (j - nm) * 11 + jj];
		}
	    }
	  for (ii = 0; ii < 5; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + 6 + ii) * n11 + j11 + jj]
		    = mat_hh [(i5 + ii) * nh + nm5 + (j - nm) * 11 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = 0; j < nm; ++j)
	{
	  j11 = j * 11;
	  j6 = j * 6;
	  j5 = j * 5;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hl */
	      for (jj = 0; jj < 6; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_hl [(nm5 + (i - nm) * 11 + ii) * nl + j6 + jj];
		}
	      /* hh */
	      for (jj = 0; jj < 5; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + 6 + jj]
		    = mat_hh [(nm5 + (i - nm) * 11 + ii) * nh + j5 + jj];
		}
	    }
	}
    }

  for (i = nm; i < np; ++i)
    {
      i11 = i * 11;
      for (j = nm; j < np; ++j)
	{
	  j11 = j * 11;
	  for (ii = 0; ii < 11; ++ii)
	    {
	      /* hh */
	      for (jj = 0; jj < 11; ++jj)
		{
		  mat [(i11 + ii) * n11 + j11 + jj]
		    = mat_hh [(nm5 + (i - nm) * 11 + ii) * nh
			     + nm5 + (j - nm) * 11 + jj];
		}
	    }
	}
    }
}

/* this is just a test routine */
static void
test_symmetric (int n, double * mat, double tiny)
{
  int i, j;
  double d;


  for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
	{
	  d = fabs (mat [i * n + j] - mat [j * n + i]);
	  if (d > tiny)
	    fprintf (stderr, "mat [%d, %d] != mat [%d, %d], "
		     "|%f - %f| = %e\n",
		     i, j, j, i,
		     mat [i * n + j], mat [j * n + i],
		     d);
		     
	}
    }
}

/* return A_ij = A_ik . B_kj
 * INPUT
 *  a [n * n] :
 *  b [n * n] :
 * OUTPUT
 *  a [n * n] :
 */
static void
multiply_matrices (int n, double *a, double *b)
{
  int i, j, k;
  double * tmp;

  tmp = malloc (sizeof (double) * n * n);

  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
	{
	  tmp [i * n + j] = 0.0;
	  for (k = 0; k < n; ++k)
	    {
	      tmp [i * n + j] += a [i * n + k] * b [k * n + j];
	    }
	}
    }

  for (i = 0; i < n * n; ++i)
    a [i] = tmp [i];

  free (tmp);
}
