/* header file for fts.c --
 * subroutine for the procedure of FTS version
 * Copyright (C) 2000-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: fts.h,v 1.1 2006/09/26 01:12:11 ichiki Exp $
 */
#ifndef	_FTS_H_
#define	_FTS_H_


/* store matrix in FTS format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n11 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n11 * n11] :
 */
void
matrix_fts_ij (int i, int j,
	       double ex, double ey, double ez,
	       double xa, double ya,
	       double yb,
	       double xc, double yc,
	       double xg, double yg,
	       double yh,
	       double xm, double ym, double zm,
	       int n11, double *mat);

/* store matrix in G part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xg, yg : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2,3,4) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_G (int n, double *mat,
	     double ex, double ey, double ez,
	     double xg, double yg);

/* store matrix in G-tilde part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xg, yg : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2,3,4)] : added (not cleared!)
 */
void
matrix_ij_Gt (int n, double *mat,
	      double ex, double ey, double ez,
	      double xg, double yg);

/* store matrix in H part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   yh : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2,3,4) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_H (int n, double *mat,
	     double ex, double ey, double ez,
	     double yh);

/* store matrix in H-tilde part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   yh : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2,3,4)] : added (not cleared!)
 */
void
matrix_ij_Ht (int n, double *mat,
	      double ex, double ey, double ez,
	      double yh);

/* store matrix in M part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xm, ym, zm : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2,3,4) * n + (0,1,2,3,4)] : added (not cleared!)
 */
void
matrix_ij_M (int n, double *mat,
	     double ex, double ey, double ez,
	     double xm, double ym, double zm);


/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FTS format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [11] : FTS of particle 'i' (extracted form)
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [11] : UOE of particle 'j' (extracted form)
 */
void
matrix_fts_atimes (double *x, double *y,
		   double ex, double ey, double ez,
		   double xa, double ya,
		   double yb,
		   double xc, double yc,
		   double xg, double yg,
		   double yh,
		   double xm, double ym, double zm);

/* convert fts[] to f[], t[], s[] (this is applicable for UOE)
 * INPUT
 *  n : # particles
 *  fts [n * 11] :
 * OUTPUT
 *  f[n * 3] :
 *  t[n * 3] :
 *  s[n * 5] :
 */
void
set_FTS_by_fts (int n,
		double *f, double *t, double *s,
		double *fts);

/* convert fts[] to f[], t[], s[] (this is applicable for UOE)
 * INPUT
 *  n : # particles
 *  f[n * 3] :
 *  t[n * 3] :
 *  s[n * 5] :
 * OUTPUT
 *  fts [n * 11] :
 */
void
set_fts_by_FTS (int n,
		double *fts,
		double *f, double *t, double *s);

/* calc scalar functions of (M^inf)^-1 in FTS
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   lub [22] : scalar functions
 */
void
scalar_minv_fts (double s,  double * scalar_fts);

/* calculate lubrication fts by uoe for all particles
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_3fts (struct stokes * sys,
	       double * uoe, double * fts);

/* calculate fts by uoe for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b (double *uoe1, double *uoe2,
		 double *x1, double *x2,
		 double *fts1, double *fts2);

/* calculate fts by uoe for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   i : particle index for '1'
 *   j : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   n : dimension of matrix 'mat'
 * OUTPUT
 *   mat [n * n] : add for (i,j)-pair
 */
void
matrix_lub_fts_2b (int i, int j,
		   double *x1, double *x2,
		   int n, double * mat);

#endif /* !_FTS_H_ */
