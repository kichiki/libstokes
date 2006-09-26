/* header file for f.c --
 * subroutine for the procedure of F version
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: f.h,v 1.1 2006/09/26 01:09:51 ichiki Exp $
 */
#ifndef	_F_H_
#define	_F_H_


/* store matrix in F format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n3 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n3 * n3] :
 */
void
matrix_f_ij (int i, int j,
	     double ex, double ey, double ez,
	     double xa, double ya,
	     int n3, double *mat);

/* store matrix in A part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_A (int n, double *mat,
	     double ex, double ey, double ez,
	     double xa, double ya);

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in F format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [3] : F of particle 'i'
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [3] : U of particle 'j'
 */
void
matrix_f_atimes (double *x, double *y,
		 double ex, double ey, double ez,
		 double xa, double ya);

/* convert f[] to f[] (this is applicable for U)
 * INPUT
 *  n : # particles
 *  f2 [n * 3] :
 * OUTPUT
 *  f1 [n * 3] := f2 [n * 3]
 */
void
set_F_by_f (int n,
	    double *f1,
	    double *f2);

/* calc scalar functions of (M^inf)^-1 in FT
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *  scalar_f [4] :
 */
void
scalar_minv_f (double s, double * scalar_f);

/* calculate lubrication f by u for all particles
 * INPUT
 *  sys : system parameters
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
void
calc_lub_3f (struct stokes * sys,
	     double * u, double * f);

/* calculate f by u for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   u1 [3] : velocity
 *   u2 [3] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   f1 [3] : force
 *   f2 [3] :
 */
void
calc_lub_f_2b (double * u1, double * u2,
	       double * x1, double * x2,
	       double * f1, double * f2);

#endif /* !_F_H_ */
