/* header file for brownian.c --
 * Brownian dynamics code
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: brownian.h,v 1.16 2008/05/24 05:49:11 kichiki Exp $
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
#ifndef	_BROWNIAN_H_
#define	_BROWNIAN_H_


#include <KIrand.h>

struct BD_params
{
  /* note that the following pointers are just pointers, therefore, 
   * you have to take care of them (to free, for example).
   */
  struct stokes *sys;
  struct KIrand *rng;
  double *F;
  double *T;
  double *E;
  double *uf;
  double *of;
  double *ef;

  int flag_noHI;
  int flag_mat;
  int flag_lub;
  int flag_lub_B; // for BD_calc_FB(), lub among mobile particles

  // auxiliary imposed-flow parameters for simple shear
  double t0; // reference time for s0
  double s0; // cell shift at time t0 (for shear_mode = 1 or 2)

  double st; // currently this is just place holders

  struct bonds *bonds;
  double gamma;
  struct EV *ev;
  struct angles *ang;
  struct EV_DH *ev_dh;
  struct EV_LJ *ev_LJ;
  struct confinement *cf;

  int flag_Q;

  // parameters for Brownian dynamics
  double peclet;
  double eps;

  int n_minv;
  double eig_minv[2];
  double *a_minv;
  int n_lub;
  double eig_lub[2];
  double *a_lub;

  int scheme;  /* 0 : the mid-point algorithm
		* 1 : Banchio-Brady (2003)
		* 2 : Ball-Melrose (1997)
		* 3 : Jendrejack et al (2000)
		* 4 : semi-implicit predictor-corrector
		* note that, for 3 and 4, BD_imp_ode_evolve() should be used.
		*/
  double BB_n; // step parameter for BB03 algorithm

  double rmin;   /* minimum distance for dt-adjustment
		  * if (r < rmin*(ai+aj)), dt is adjusted
		  * rmin == 0 corresponds to "no dt-adjustment"
		  * rmin == 1 corresponds to dt-adjustment for overlap
		  * NOTE: if sys->rmin is defined, dt-adjustment is ignored.
		  */
  double dt_lim; /* lower bound to shrink dt to prevent overlaps
		  * set "dt" if you don't want to adjust dt but just reject
		  */
};


/* set the parameters to struct BD_params
 * INPUT
 *  ** NOTE ** the following pointers are just pointers.
 *             you have to take care of them! (free, for example.)
 *  (struct stokes *)sys -- initialize before calling!
 *  seed : for random number generator
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_noHI
 *  (int) flag_lub
 *  (int) flag_mat
 *        NOTE, flag_lub_B is used for BD_calc_FB() where
 *        lub among ONLY mobile particles are taken.
 *        therefore, check for the existance of mobile pair(s)
 *        not excluded by sys->ex_lub for flag_lub_B.
 *  (double) stokes -- currently this is just a place holder
 *  (struct bonds *)bonds
 *  (double) gamma
 *  (struct EV *)ev
 *  (struct angles *)ang
 *  (struct EV_DH *)ev_dh
 *  (struct EV_LJ *)ev_LJ
 *  (struct confinement *)cf
 *  (int) flag_Q
 *  (double) peclet
 *  (double) eps
 *  (int) n_minv
 *  (int) n_lub
 *  (int) scheme
 *  (double) BB_n
 *  (double) rmin
 *  (double) dt_lim
 * OUTPUT :
 *  (struct ode_params) params
 */
struct BD_params *
BD_params_init (struct stokes *sys,
		unsigned long seed,
		double *F,
		double *T,
		double *E,
		double *uf,
		double *of,
		double *ef,
		int flag_noHI,
		int flag_lub,
		int flag_mat,
		double st,
		struct bonds *bonds,
		double gamma,
		struct EV *ev,
		struct angles *ang,
		struct EV_DH *ev_dh,
		struct EV_LJ *ev_LJ,
		struct confinement *cf,
		int flag_Q,
		double peclet,
		double eps,
		int    n_minv,
		int    n_lub,
		int    scheme,
		double BB_n,
		double rmin,
		double dt_lim);

void
BD_params_free (struct BD_params *BD);

/* set the reference for cell-shift (shear_mode = 1 and 2)
 * INTPUT
 *  t0 : reference time for s0
 *  s0 : cell shift at time t0 (for shear_mode = 1 or 2)
 * OUTPUT
 *  BD->t0, BD->s0 :
 */
void
BD_set_shear_shift_ref (struct BD_params *BD,
			double t0, double s0);

/* for check_overlap()
 */
struct overlap {
  double r2; // square of distance for the overlapping pair
  double a2; // (a_i + a_j)^2 for the overlapping pair
  int i; // particle index
  int j; // particle index
  int k; // lattice index
};

/* return numbers of overlapping pairs
 * INPUT
 *  rmin : factor to judge "overlap"
 *         where (r2 <= rmin * a2) are overlapping
 *         rmin == 1 => contact point is the condition
 *         rmin == 0 => all pairs are NOT overlapping
 * OUTPUT
 *  returned value : number of "overlapping" pairs (r2 <= rmin * a2)
 *                   (0 == no "overlapping" pairs)
 */
int
check_overlap (struct stokes *sys, const double *pos, double rmin,
	       struct overlap *ol);


/* calculate y = A.x for Brownian force, where A = M^inf
 * so that give 1/sqrt(x) for chebyshev.
 * Note that for FTS version, A = M_{UF} - M_{US}.(M_{ES})^{-1}.M_{EF}.
 * INPUT
 *  n    : dimension (n = 3*nm for F, n = 6*nm for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : F (and T) in Global scaling, rather than Particle-wise scaling
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : U (and O) in Global scaling, rather than Particle-wise scaling
 */
void
BD_atimes_mob_FU (int n, const double *x, double *y, void *user_data);


/* calculate y = A.x for Brownian force, where A = L (lubrication)
 * so that give sqrt(x) for chebyshev.
 * note that FU part (for mobile particles) are just extracted
 * (in other words, F for the fixed particles is set by zero,
 *  and E is also set by zero for FTS case).
 * INPUT
 *  n    : dimension (n = 3*nm for F, n = 6*nm for FT and FTS)
 *         S component is not included (because it is not random variable.)
 *  x[n] : U (and O) in Global scaling, rather than Particle-wise scaling
 *  user_data : (struct BD_params) *BD
 * OUTPUT
 *  y[n] : F (and T) in Global scaling, rather than Particle-wise scaling
 */
void
BD_atimes_lub_FU (int n, const double *x, double *y, void *user_data);


/* calc (M^{-1})_{FU} in FTS version
 * where M=(a b) and M^{-1}=(A B), 
 *         (c d)            (C D)
 * A = (a - b.d^{-1}.c)^{-1}.
 * INPUT
 *  np : number of particles
 *  m[np11 * np11] : mobility matrix in FTS
 * OUTPUT
 *  minv_FU[np6 * np6]
 */
void
BD_minv_FU_in_FTS (int np, const double *m, double *minv_FU);


/* make mobility matrix (M^inf)^{-1} in UF part (for mobile particles)
 * (in other words, F for the fixed particles is set by zero,
 *  and S is also set by zero for FTS case).
 * INPUT
 *  BD   : struct BD_params
 * OUTPUT
 *  minv : UF part of (M^inf)^{-1} in Global scaling not Particle-wise.
 */
void
BD_matrix_minv_FU (struct BD_params *BD, double *minv);


/* make lubrication matrix L in UF part (for mobile particles)
 * (in other words, F for the fixed particles is set by zero,
 *  and S is also set by zero for FTS case).
 * INPUT
 *  BD   : struct BD_params
 * OUTPUT
 *  lub  : UF part of L in Global scaling not Particle-wise.
 */
void
BD_matrix_lub_FU (struct BD_params *BD, double *lub);


/* calculate sqrt of the matrix a[n*n] by dgeev()
 */
int
BD_sqrt_by_dgeev (int n, const double *a, double *s);


/*
 * INPUT
 *  BD   : struct BD_params
 *         (sys, rng, eig, n_minv, a_minv, n_lub, a_lub, eps are used)
 * OUTPUT
 *  z[n] : random vector, with which F^B = z * sqrt(2/(peclet * dt))
 *         in FT and FTS, first nm3 are the force, the next nm3 are the torque
 *         (different strage from FTS where f,t,s are ordered particle-wise).
 */
void
BD_calc_FB (struct BD_params *BD,
	    double *z);


/*
 * INPUT
 *  sys : struct stokes
 *  u[nm*3] : velocity ( = dx/dt)
 *  o[nm*3] : angluar velocity, which is converted to dq/dt.
 *  x[nm*3] : present position
 *  q[nm*4] : present quaternion
 *  dt      : time step
 * OUTPUT
 *  x[nm*3] : updated position (x += u * dt)
 *  q[nm*4] : present quaternion (q += dq/dt * dt)
 */
void
evolve_Euler_3all (struct stokes *sys,
		   const double *u, const double *o,
		   double dt,
		   double *x, double *q);


struct FTS
{
  int nm;
  int nf;
  double *f;
  double *t;
  double *s;
  double *ff;
  double *tf;
  double *sf;
  double *u;
  double *o;
  double *e;
  double *uf;
  double *of;
  double *ef;
};

/*
 * INPUT
 *  sys : struct stokes
 *        version, np, nm are used.
 */
struct FTS *
FTS_init (struct stokes *sys);

void
FTS_free (struct FTS *FTS);


/* dimension for BD scheme:
 * n = nm * 3 for F version -- velocity of mobile particles
 *   = nm * 6 for FT and FTS versions
 *                -- translational and angular velocities of mobile particles
 */
int
BD_get_n (struct stokes *sys);


/* wrapper for BD_evolve()
 * INPUT
 *  *t    : (input) current time
 *  t_out : output time
 *  *dt   : (input) current inner time-step
 *  y[n]  : (input) current configuration at (*t),
 *          where n = nm*3 for F version
 *                    (y[] in the first (nm*3) elements = velocity),
 *                n = nm*3 + nm*4 with quaternion
 *                    (y[] in the first (nm*3) elements = velocity,
 *                     y[] in the next (nm*4) elements = quaternion).
 * OUTPUT
 *  *t    : (output) output time (= t_out)
 *  *dt   : (output) current (updated, if necessary) inner time-step
 *  y[n]  : (output) updated configuration at t_out
 */
void
BD_ode_evolve (struct BD_params *BD,
	       double *t, double t_out, double *dt,
	       double *y);

/* return factor for BD->z[]
 */
double
BD_params_get_fact (struct BD_params *BD,
		    double dt);

/* add interaction forces on each particle including
 *   bond force F^bond(sys->pos[]) by bonds_calc_force()
 *   EV force F^EV(sys->pos[]) by EV_calc_force()
 * INPUT
 *  BD : struct BD_params
 *  pos[] : position for F^bond and F^EV
 * OUTPUT
 *  FTS : struct FTS
 */
void
BD_add_FP (struct BD_params *BD,
	   const double *pos,
	   struct FTS *FTS);

/* evolve position of particles -- the mid-point scheme
 * INPUT
 *  t       : current time
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 *  returned value : the integrated time duration
 */
double
BD_evolve_mid (double t,
	       struct BD_params *BD,
	       double *x, double *q,
	       double dt);

/* evolve position of particles -- Banchio-Brady scheme
 * reference : Banchio and Brady (2003) Phys. Fluids
 * INPUT
 *  t       : current time
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 *  returned value : the integrated time duration
 */
double
BD_evolve_BB03 (double t,
		struct BD_params *BD,
		double *x, double *q,
		double dt);

/* evolve position of particles -- Ball-Melrose scheme
 * reference : Ball and Melrose (1997)
 * INPUT
 *  t       : current time
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 *  returned value : the integrated time duration
 */
double
BD_evolve_BM97 (double t,
		struct BD_params *BD,
		double *x, double *q,
		double dt);


#endif /* !_BROWNIAN_H_ */
