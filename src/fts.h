/* header file for fts.c --
 * subroutine for the procedure of FTS version
 * Copyright (C) 2000-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: fts.h,v 2.8 2007/04/26 05:13:34 kichiki Exp $
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

/* store matrix in G-transpose part with scalar functions
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
matrix_ij_GT (int n, double *mat,
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

/* store matrix in H-transpose part with scalar functions
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
matrix_ij_HT (int n, double *mat,
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
matrix_fts_atimes (const double *x,
		   double *y,
		   double ex, double ey, double ez,
		   double xa, double ya,
		   double yb,
		   double xc, double yc,
		   double xg, double yg,
		   double yh,
		   double xm, double ym, double zm);

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FTS format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- alpha(i)' interaction is stored.
 * INPUT
 *   x [11] : FTS of particle 'i' (extracted form)
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [11] : UOE of particle 'j' (extracted form)
 */
void
matrix_fts_self_atimes (const double *x,
			double *y,
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
		const double *fts);

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
		const double *f, const double *t, const double *s);

/* calc scalar functions of (M^inf)^-1 in FTS
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   lub [22] : scalar functions
 */
void
scalar_minv_fts (double s,  double * scalar_fts);

/* calculate fts by uoe for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin is used.
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b (struct stokes * sys,
		 const double *uoe1, const double *uoe2,
		 const double *x1, const double *x2,
		 double *fts1, double *fts2);

/* calculate lub-matrix in FTS version for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin is used.
 *   i : particle index for '1'
 *   j : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   n : dimension of matrix 'mat' (must be np*11)
 * OUTPUT
 *   mat [n * n] : add for (i,j)-pair
 */
void
matrix_lub_fts_2b (struct stokes * sys,
		   int i, int j,
		   const double *x1, const double *x2,
		   int n, double * mat);

/* calculate fts by uoe for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin       : min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b_poly (struct stokes *sys,
		      const double *uoe1, const double *uoe2,
		      const double *x1, const double *x2,
		      int i1, int i2,
		      double *fts1, double *fts2);

/* calculate lub-matrix in FTS version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin       : min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*11)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_fts_2b_poly (struct stokes *sys,
			int i, int j,
			const double *x1, const double *x2,
			int i1, int i2,
			int n, double *mat);

/* pre-process for imposed flow shifting, that is, converting E
 * from the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * to the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  e[np*5] : strain in the fluid-rest frame
 *            (data is preserved)
 * OUTPUT
 *  e0[np*5] : strain in the fluid-rest frame
 */
void
shift_labo_to_rest_E (struct stokes * sys,
		      int np, const double *e,
		      double *e0);


#endif /* !_FTS_H_ */
