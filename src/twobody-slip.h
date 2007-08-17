/* header file for twobody-slip.c --
 * twobody solutions for slip particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: twobody-slip.h,v 1.1 2007/08/17 04:31:19 kichiki Exp $
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
#ifndef	_TWOBODY_SLIP_H_
#define	_TWOBODY_SLIP_H_


struct twobody_slip_f {
  int nmax;
  double lambda;  // = a2 / a1
  double hat_g1;  // = gamma1 / a1
  double hat_g2;  // = gamma2 / a2
  double slip_a1; // = a1 * sqrt(Gamma^(1)(0,2)
  double slip_a2; // = a2 * sqrt(Gamma^(2)(0,2)

  double *XA;
  double *YA;
  double *YB;
  double *XC;
  double *YC;
  double *XG;
  double *YG;
  double *YH;
  double *XM;
  double *YM;
  double *ZM;
};

struct twobody_slip_f_list {
  int n;
  struct twobody_slip_f **f;
};



void twobody_XA_slip (int n, double l, double *f);
void twobody_YA_slip (int n, double l, double *f);
void twobody_YB_slip (int n, double l, double *f);
void twobody_XC_slip (int n, double l, double *f);
void twobody_YC_slip (int n, double l, double *f);
void twobody_XG_slip (int n, double l, double *f);
void twobody_YG_slip (int n, double l, double *f);
void twobody_YH_slip (int n, double l, double *f);
void twobody_XM_slip (int n, double l, double *f);
void twobody_YM_slip (int n, double l, double *f);
void twobody_ZM_slip (int n, double l, double *f);


/* 
 * INPUT
 *  (extern) lambda : slip length scaled by the radius
 *                    if negative, the slip length is infinity (perfect slip)
 */
double SL_G1 (int m, int n);
double SL_G2 (int m, int n);


/** utility routines for struct twobody_f and twobody_f_list **/
struct twobody_slip_f *
twobody_slip_f_init (int nmax,
		     double lambda,
		     double hat_g1, double hat_g2,
		     double slip_a1, double slip_a2);

void
twobody_slip_f_free (struct twobody_slip_f *f);

struct twobody_slip_f_list *
twobody_slip_f_list_init (void);

void
twobody_slip_f_list_append (struct twobody_slip_f_list *list,
			    int nmax,
			    double lambda,
			    double hat_g1, double hat_g2,
			    double slip_a1, double slip_a2);

void
twobody_slip_f_list_free (struct twobody_slip_f_list *list);


/** far form **/

/* calc XA11 and XA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XA11, double *XA12);

/* calc YA11 and YA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YA11, double *YA12);

/* calc YB11 and YB12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YB11, double *YB12);

/* calc XC11 and XC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XC11, double *XC12);

/* calc YC11 and YC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YC11, double *YC12);

/* calc XG11 and XG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XG11
 *  XG12
 */
void twobody_XG_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XG11, double *XG12);

/* calc YG11 and YG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YG11
 *  YG12
 */
void twobody_YG_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YG11, double *YG12);

/* calc YH11 and YH12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YH11
 *  YH12
 */
void twobody_YH_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YH11, double *YH12);

/* calc XM11 and XM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  XM11
 *  XM12
 */
void twobody_XM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *XM11, double *XM12);

/* calc YM11 and YM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  YM11
 *  YM12
 */
void twobody_YM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *YM11, double *YM12);

/* calc ZM11 and ZM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  ZM11
 *  ZM12
 */
void twobody_ZM_slip_far (int n, double l, double s,
			  double hat_g1, double hat_g2,
			  double *ZM11, double *ZM12);

/* calc scalar functions of resistance problem by 1/s expansion
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  far [22] : scalar functions
 *      0, 1 : (XA11, XA12)
 *      2, 3 : (YA11, YA12)
 *      4, 5 : (YB11, YB12)
 *      6, 7 : (XC11, XC12)
 *      8  9 : (YC11, YC12)
 *     10,11 : (XG11, XG12)
 *     12,13 : (YG11, YG12)
 *     14,15 : (YH11, YH12)
 *     16,17 : (XM11, XM12)
 *     18,19 : (YM11, YM12)
 *     20,21 : (ZM11, ZM12)
 */
void twobody_slip_far (int version, int n, double l, double s,
		       double hat_g1, double hat_g2,
		       double *far);

/* calc scalar functions of resistance problem by 1/s expansion
 * all-in-one form (to reduce calculating the same parameters)
 * and with struct twobody_f *f12 table (to avoid recalculating them)
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  f12     : struct twobody_f for the pair
 *            you can give NULL for them.
 *            then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 *  hat_g1 := gamma1 / a1
 *  hat_g2 := gamma2 / a2
 * OUTPUT
 *  far [22] : scalar functions
 *      0, 1 : (XA11, XA12)
 *      2, 3 : (YA11, YA12)
 *      4, 5 : (YB11, YB12)
 *      6, 7 : (XC11, XC12)
 *      8  9 : (YC11, YC12)
 *     10,11 : (XG11, XG12)
 *     12,13 : (YG11, YG12)
 *     14,15 : (YH11, YH12)
 *     16,17 : (XM11, XM12)
 *     18,19 : (YM11, YM12)
 *     20,21 : (ZM11, ZM12)
 */
void twobody_slip_far_with_f (int version,
			      struct twobody_slip_f *f12,
			      int n, double l, double s,
			      double hat_g1, double hat_g2,
			      double *far);


/* calc scalar functions of two-body exact solution in resistance problem
 * INPUT
 *  version    : 0=F, 1=FT, 2=FTS.
 *  r          : distance between the two := x_b - x_a
 *  aa, ab     : radii for particles a(alpha) and b(beta)
 *  f12        : (struct twobody_slip_f *).
 *               you can give NULL for them.
 *               then, the coefs are calculated on-the-fly (terribly slow).
 *  n          : max order for the coefficients
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 *               (*** currently, lub form is not implemented ***)
 *  flag_scale : 0 no scaling, that is, in Jeffrey form
 *               1 for the dimensional form
 *               2 for the Stokesian dynamics form
 *  res [22]   :
 * OUTPUT
 *  res [22]   : scalar functions. the scaling is given by flag_scale.
 */
void
twobody_slip_scalars_res (int version,
			  double r,
			  double aa, double ab,
			  struct twobody_slip_f *f12,
			  int n, int flag_lub, int flag_scale,
			  double *res);

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12,f21: (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  flag_lub   : 0 to use twobody_slip_far()
 *               1 to use twobody_slip_lub()
 *               (*** currently, lub form is not implemented ***)
 * OUTPUT
 *  lub [44] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *    8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *   12,13,14,15 : (XC11, XC12, XC21, XC22)
 *   16,17,18,19 : (YC11, YC12, YC21, YC22)
 *   20,21,22,23 : (XG11, XG12, XG21, XG22)
 *   24,25,26,27 : (YG11, YG12, YG21, YG22)
 *   28,29,30,31 : (YH11, YH12, YH21, YH22)
 *   32,33,34,35 : (XM11, XM12, XM21, XM22)
 *   36,37,38,39 : (YM11, YM12, YM21, YM22)
 *   40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 */
void
scalars_lub_slip_full (int version,
		       double r, double a1, double a2,
		       struct twobody_slip_f *f12,
		       struct twobody_slip_f *f21,
		       int n, int flag_lub,
		       double *lub);

#include "stokes.h" // struct stokes
/** F version **/

/* calculate f by u for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   u1 [3] : velocity of particle 1
 *   u2 [3] : velocity of particle 2
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   f1 [3] : force of particle 1
 *   f2 [3] : force of particle 2
 */
void
calc_lub_f_2b_slip (struct stokes *sys,
		    const double *u1, const double *u2,
		    const double *x1, const double *x2,
		    int i1, int i2,
		    double *f1, double *f2);

/* calculate lub-matrix in F version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b_slip (struct stokes *sys,
		      int i, int j,
		      const double *x1, const double *x2,
		      int i1, int i2,
		      int n, double *mat);

/** FT version **/

/* calculate ft by uo for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uo1 [6] : velocity, angular velocity
 *   uo2 [6] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   ft1 [6] : force, torque
 *   ft2 [6] :
 */
void
calc_lub_ft_2b_slip (struct stokes *sys,
		     const double *uo1, const double *uo2,
		     const double *x1, const double *x2,
		     int i1, int i2,
		     double *ft1, double *ft2);

/* calculate lub-matrix in FT version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*6)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_ft_2b_slip (struct stokes *sys,
		       int i, int j,
		       const double *x1, const double *x2,
		       int i1, int i2,
		       int n, double *mat);

/** FTS version **/

/* calculate fts by uoe for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
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
calc_lub_fts_2b_slip (struct stokes *sys,
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
 *            sys->lubmin2      : square of min distance for lub calculation.
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
matrix_lub_fts_2b_slip (struct stokes *sys,
			int i, int j,
			const double *x1, const double *x2,
			int i1, int i2,
			int n, double *mat);


#endif /* !_TWOBODY_SLIP_H_ */
