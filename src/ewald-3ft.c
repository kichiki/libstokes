/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2001 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3ft.c,v 3.2 2001/02/03 08:57:59 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 3D configuration
 * periodic boundary condition in 3 direction,
 * FT version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */

#ifdef ZETA
#include <time.h> /* clock() */
#endif /* ZETA */

#include <libiter.h> /* solve_iter_stab (), gpb () */

#include "ft.h"
#include "lub.h"
#include "ewald-3ft.h"


/* (local) global variable */
int NUMB_mobile_particles; /* this is dirty, though ... */

/** function prototypes for local routines **/
/* utility routines for calc_mob_fix_ewald_3ft () */
static void
calc_b_mob_fix_ewald_3ft (int np, int nm,
			  double *f, double *t,
			  double *uf, double *of,
			  double *b);
static void
atimes_mob_fix_ewald_3ft (int n, double *x, double *y);


/* utility routines for calc_mob_lub_fix_ewald_3ft () */
static void
calc_b_mob_lub_fix_ewald_3ft (int np, int nm,
			      double *f, double *t,
			      double *uf, double *of,
			      double *b);
static void
atimes_mob_lub_fix_ewald_3ft (int n, double *x, double *y);
static void
calc_lub_ewald_3ft (int np, double * uo, double * ft);
static int
cond_lub (double * x1, double * x2);


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FT version
 * INPUT
 *  (global) pos [] : position of particles
 *  n := np * 11
 *  x [n * 6] : FT
 * OUTPUT
 *  y [n * 6] : UO
 */
void
atimes_ewald_3ft (int n, double *x, double *y)
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
  clock_t ctmp1, ctmp2, ctmp3;
#endif /* ZETA */

  double xa, ya; 
  double yb;
  double xc, yc;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int np;
  int i, j;
  int i6, j6;
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


  np = n / 6;

  /* clear result */
  for (i = 0; i < n; i ++)
    y [i] = 0.0;

  /* diagonal part ( self part ) */
  xa = ya = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
  xc = yc = 0.75 - zaspi * za2 * 10.0;

  for (i = 0; i < np; i++)
    {
      i6 = i * 6;
      matrix_ft_atimes (x + i6, y + i6,
			 0.0, 0.0, 0.0,
			 xa, ya,
			 0.0,
			 xc, yc);
    }

#ifdef ZETA
  ctmp1 = clock ();
#endif /* ZETA */

  /* first Ewald part ( real space ) */
  for (i = 0; i < np; i++)
    {
      i6 = i * 6;
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < np; j++)
	{
	  j6 = j * 6;
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
	      
			  matrix_ft_atimes (x + i6, y + j6,
					     ex, ey, ez,
					     xa, ya,
					     yb,
					     xc, yc);
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  ctmp2 = clock ();
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
      
		  for (i = 0; i < np; i++)
		    {
		      i6 = i * 6;
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      for (j = 0; j < np; j++)
			{
			  j6 = j * 6;
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

			  matrix_ft_atimes (x + i6, y + j6,
					     ex, ey, ez,
					     0.0, cf * ya,
					     sf * yb,
					     0.0, cf * yc);
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  ctmp3 = clock ();

  cpu1 = (double) (ctmp3 - ctmp1);
  cpu2 = (double) (ctmp2 - ctmp1);
  cpu3 = (double) (ctmp3 - ctmp2);
#endif /* ZETA */
}


/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_ewald_3ft (int np,
		     double *u, double *o,
		     double *f, double *t)
{
  int i;
  int n6;

  double *b;
  double *x;


  n6 = np * 6;
  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3ft ().\n");
      exit (1);
    }

  set_ft_by_FT (np, b, u, o);

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  solve_iter_stab (n6, b, x, atimes_ewald_3ft,
		   gpb);

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
}

/** natural mobility problem **/
/* solve natural mobility problem in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
calc_mob_ewald_3ft (int np,
		    double *f, double *t,
		    double *u, double *o)
{
  int n6;

  double *b;
  double *x;


  n6 = np * 6;
  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3ft ().\n");
      exit (1);
    }

  set_FT_by_ft (np, f, t, x);
  atimes_ewald_3ft (n6, x, b);
  set_ft_by_FT (np, b, u, o);

  free (b);
  free (x);
}


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
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
calc_mob_fix_ewald_3ft (int np, int nm,
			 double *f, double *t,
			 double *uf, double *of,
			 double *u, double *o,
			 double *ff, double *tf)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int n6;
  int nf;
  int nm6;

  double *b;
  double *x;


  nf = np - nm;
  n6 = np * 6;
  nm6 = nm * 6;

  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3ft ().\n");
      exit (1);
    }

  calc_b_mob_fix_ewald_3ft (np, nm, f, t, uf, of, b);

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  NUMB_mobile_particles = nm;
  solve_iter_stab (n6, b, x, atimes_mob_fix_ewald_3ft,
		   gpb);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
 * where b := - [(0,0)_m,(u,o)_f] + M.[(f,t)_m,(0,0)_f].
 * INPUT
 *  np : # all particles (not # elements in b[]!)
 *  nm : # mobile particles (not # elements in b[]!)
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  b [np * 6] : constant vector
 */
static void
calc_b_mob_fix_ewald_3ft (int np, int nm,
			  double *f, double *t,
			  double *uf, double *of,
			  double *b)
{
  int i;
  int i3, i6;
  int j;
  int nf;
  int n3, n6;
  int nm6;

  double *x;
  double *v3_0;


  nf = np - nm;
  n3 = np * 3;
  n6 = np * 6;
  nm6 = nm * 6;

  x = malloc (sizeof (double) * n6);
  v3_0 = malloc (sizeof (double) * n3);
  if (x == NULL
      || v3_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_ewald_3ft ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);
  atimes_ewald_3ft (n6, x, b);

  /* set b := M.x - [(0,0)_m,(U,O)_f] */
  for (i = 0; i < nf; ++i)
    {
      i3 = i * 3;
      i6 = (i + nm) * 6;
      for (j = 0; j < 3; ++j)
	{
	  b [i6 + j] -= uf [i3 + j];
	  b [i6 + 3 + j] -= of [i3 + j];
	}
    }

  free (x);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  (global) : NUMB_mobile_particles -- this is dirty, though ...
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_fix_ewald_3ft (int n, double *x, double *y)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3, nm6;

  double *z;
  double *v3_0;
  double *u;
  double *o;
  double *ff;
  double *tf;


  np = n / 6;
  np3 = np * 3;
  nm = NUMB_mobile_particles;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;
  nm6 = nm * 6;

  z = malloc (sizeof (double) * n);
  v3_0 = malloc (sizeof (double) * np3);
  u = malloc (sizeof (double) * nm3);
  o = malloc (sizeof (double) * nm3);
  ff = malloc (sizeof (double) * nf3);
  tf = malloc (sizeof (double) * nf3);
  if (z == NULL
      || v3_0 == NULL
      || u == NULL
      || o == NULL
      || ff == NULL
      || tf == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_ewald_3ft ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  /* set y := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, y, v3_0, v3_0);
  set_ft_by_FT (nf, y + nm6, ff, tf);
  atimes_ewald_3ft (n, y, z);

  /* set y := [(U,O,0)_mobile,(0,0,0)_fixed] */
  set_ft_by_FT (nm, y, u, o);
  set_ft_by_FT (nf, y + nm6, v3_0, v3_0);

  /* set y := [(U,O)_m,(0,0)_f] - M.[(0,0)_m,(F,T)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v3_0);
  free (u);
  free (o);
  free (ff);
  free (tf);
}


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
calc_res_lub_ewald_3ft (int np,
			 double *u, double *o,
			 double *f, double *t)
{
  int i;
  int n6;

  double *b;
  double *x;
  double *lub;


  n6 = np * 6;
  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  lub = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL
      || lub == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_lub_ewald_3ft ().\n");
      exit (1);
    }

  set_ft_by_FT (np, b, u, o);
  calc_lub_ewald_3ft (np, b, lub);
  atimes_ewald_3ft (n6, lub, x); // x[] is used temporaly
  for (i = 0; i < n6; ++i)
    b [i] += x [i];

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  solve_iter_stab (n6, b, x, atimes_ewald_3ft,
		   gpb);

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
  free (lub);
}


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
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
calc_mob_lub_fix_ewald_3ft (int np, int nm,
			     double *f, double *t, double *e,
			     double *uf, double *of, double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int n6;
  int nf;
  int nm6;

  double *b;
  double *x;


  nf = np - nm;
  n6 = np * 6;
  nm6 = nm * 6;

  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3ft ().\n");
      exit (1);
    }

  calc_b_mob_lub_fix_ewald_3ft (np, nm, f, t, uf, of, b);

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  NUMB_mobile_particles = nm;
  solve_iter_stab (n6, b, x, atimes_mob_lub_fix_ewald_3ft,
		   gpb/*sta*//*gpb_chk*/);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  np : # all particles (not # elements in b[]!)
 *  nm : # mobile particles (not # elements in b[]!)
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  b [np * 6] : constant vector
 */
static void
calc_b_mob_lub_fix_ewald_3ft (int np, int nm,
			      double *f, double *t,
			      double *uf, double *of,
			      double *b)
{
  int i;
  int nf;
  int n3, n6;
  int nm6;

  double *x;
  double *y;
  double *v3_0;


  nf = np - nm;
  n3 = np * 3;
  n6 = np * 6;
  nm6 = nm * 6;

  x = malloc (sizeof (double) * n6);
  y = malloc (sizeof (double) * n6);
  v3_0 = malloc (sizeof (double) * n3);
  if (x == NULL
      || y == NULL
      || v3_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_3ft ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b :=  - (I + M.L).[(0,0)_m,(U,O)_f] */
  set_ft_by_FT (nm, b, v3_0, v3_0);
  set_ft_by_FT (nf, b + nm6, uf, of);
  calc_lub_ewald_3ft (np, b, y);
  atimes_ewald_3ft (n6, y, x); // x[] is used temporaly
  for (i = 0; i < n6; ++i)
    {
      b [i] = - (b [i] + x [i]);
    }

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);
  atimes_ewald_3ft (n6, x, y);

  /* set b := - (I + M.L).[(0,0,E)_m,(U,O,E)_f] + [(F,T,0)_m,(0,0,0)_f] */
  for (i = 0; i < n6; ++i)
    {
      b [i] += y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  (global) : NUMB_mobile_particles -- this is dirty, though ...
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_3ft (int n, double *x, double *y)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3, nm6;

  double *w;
  double *z;
  double *v3_0;
  double *u;
  double *o;
  double *ff;
  double *tf;


  np = n / 6;
  np3 = np * 3;
  nm = NUMB_mobile_particles;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;
  nm6 = nm * 6;

  w = malloc (sizeof (double) * n);
  z = malloc (sizeof (double) * n);
  v3_0 = malloc (sizeof (double) * np3);
  u = malloc (sizeof (double) * nm3);
  o = malloc (sizeof (double) * nm3);
  ff = malloc (sizeof (double) * nf3);
  tf = malloc (sizeof (double) * nf3);
  if (z == NULL
      || v3_0 == NULL
      || u == NULL
      || o == NULL
      || ff == NULL
      || tf == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_3ft ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  /* set y := (I + M.L).[(U,O)_mobile,(0,0)_fixed] */
  set_ft_by_FT (nm, y, u, o);
  set_ft_by_FT (nf, y + nm6, v3_0, v3_0);
  calc_lub_ewald_3ft (np, y, w);
  atimes_ewald_3ft (n, w, z);
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  /* set z := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, w, v3_0, v3_0);
  set_ft_by_FT (nf, w + nm6, ff, tf);
  atimes_ewald_3ft (n, w, z);

  /* set y := (I + M.L).[(U,O)_m,(0,0)_f] - M.[(0,0)_m,(F,T)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (w);
  free (z);
  free (v3_0);
  free (u);
  free (o);
  free (ff);
  free (tf);
}


/* calculate lubrication ft by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
static void
calc_lub_ewald_3ft (int np, double * uo, double * ft)
{
  extern double * pos;
  extern double llx [27], lly [27], llz [27];

  int i, j, k;
  int i3, i6;
  int j3, j6;

  double * tmp_pos;


  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3ft().\n");
      exit (1);
    }

  /* clear ft [np * 6] */
  for (i = 0; i < np * 6; ++i)
    ft [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j6 = j * 6;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = pos [j3 + 0] + llx [k];
	      tmp_pos [1] = pos [j3 + 1] + lly [k];
	      tmp_pos [2] = pos [j3 + 2] + llz [k];
	      if (cond_lub (pos + i3, tmp_pos) == 0)
		{
		  calc_lub_2b (uo + i6, uo + j6,
			       pos + i3, tmp_pos,
			       ft + i6, ft + j6);
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
