/* header file for library 'libstokes-core'
 * Copyright (C) 1993-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_LIBSTOKES_CORE_H_
#define	_LIBSTOKES_CORE_H_

#include <stdio.h> // FILE


/***********************************
 ** system parameters             **
 ** from stokes.h                 **
 ***********************************/
struct stokes {
  int version; /* 0 = F, 1 = FT, 2 = FTS  */

  int np;      /* number of all particles  */
  int nm;      /* number of mobile particles  */
  double *pos; /* position of particles  */

  /**
   * parameters for the polydisperse system
   */
  double *a;   /* radius of particles
		* Note : NULL (default) is for monodisperse system  */
  double rmin; /* minimum distance for HI calculation 
		* as in Rotne-Prager for r < 2a.
		* if (r < rmin*(ai+aj)), r = rmin*(ai+aj)
		* (default is 0.0) */
  int twobody_nmax;// max order for the coefficient for twobody_scalars_res()
  int twobody_lub; // 0 (far form) or 1 (lub form) for twobody_scalars_res()
  int *poly_table; // for (i,j), [i*np+j] gives the index of "twobody_f_list"
  struct twobody_f_list *twobody_f_list;

  /**
   * slip parameters
   */
  /* Note: if slip == NULL, the system is treated as the no-slip */
  double *slip;     // slip length
  double *slip_a;   // effective radius for the laplacian terms
  double *slip_G32; // Lambda(3,2) = 1/Lambda(2,3) for a-self.
  double *slip_G30; // Lambda(3,0) = 1/Lambda(0,3) for c-self.
  double *slip_G52; // Lambda(5,2) = 1/Lambda(2,5) for m-self.
  // slip table -- twobody_nmax and twobody_lub below are used, too
  int *slip_table; /* for (i,j), [i*np+j] gives the index
		    * of "twobody_slip_f_list" */
  struct twobody_slip_f_list *twobody_slip_f_list;

  /**
   * imposed flow
   */
  double Ui[3];
  double Oi[3];
  double Ei[5];

  /* auxiliary imposed-flow parameters for simple shear */
  int shear_mode; /* 0 : imposed flow is given by Ui,Oi,Ei.
		   * 1 : x = flow dir, y = grad dir
		   * 2 : x = flow dir, z = grad dir
		   */
  double shear_rate;
  double shear_shift; /* shift of the periodic cell (H.rate.t) at the time t
		       * which is in [0, L)
		       * H : cell size in the grad dir
		       * L : cell size in the flow dir
		       */

  /**
   * periodic parameters
   */
  int periodic; // 0 = non periodic, 1 = periodic

  /* for ewald codes */
  double ewald_eps;
  double rmax2;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double xi, xi2, xiaspi, xia2;
  double pivol;

  double lx, ly, lz;
  int ilx[27], ily[27], ilz[27];
  // note: ll[xyz][i] = l[xyz] * (double)il[xyz][i]

  // self part
  double self_a;
  double self_c;
  double self_m;

  // table for lattice summation
  int flag_table; // 0 = inactive, 1 = active
  // real space
  int nr; // number of lattice points
  double * rlx;
  double * rly;
  double * rlz;
  int * rmx;
  int * rmy;
  int * rmz;
  // reciprocal space
  int nk; // number of lattice points
  double * ex;
  double * ey;
  double * ez;
  double * k;
  double * k1;
  double * k2;
  double * k3;
  double * ya;
  double * yb;
  double * yc;
  double * yg;
  double * yh;
  double * ym;

  /**
   * for lubrication
   */
  double lubmin2;/* square of min distance for lub
		  * (distance is replaced by its root-square)
		  */
  double lubmax; /* max distance for lub;
		  * for the pair beyond this is just ignored.
		  * 0 means no limit for open systems and 
		  * all particles within +/-1 cells in x,y,z for the periodic
		  */
  /* exclusion list for lubrication due to the bonding */
  //struct list_ex *ex_lub;

  /**
   * for zeta program
   */
  double cpu1, cpu2, cpu3;

  /**
   * for iterative solvers
   */
  struct iter * it;
};


/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void);

void
stokes_free (struct stokes * sys);

/* set np and nm and allocate the memory for pos[np*3]
 * also struct list_ex *ex_lub is allocated here
 * (becuase np is necessary for ex_lub).
 */
void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_Ui (struct stokes * sys,
	       double uix, double uiy, double uiz);
void
stokes_set_Oi (struct stokes * sys,
	       double oix, double oiy, double oiz);
void
stokes_set_Ei (struct stokes * sys,
	       double eixx, double eixy, double eixz,
	       double eiyz, double eiyy);

/* set auxiliary imposed flow parameters for simple shear
 * and overwrite the imposed parameters Ui, Oi, Ei.
 * INPUT
 *  shear_mode : 1 (x = flow dir, y = grad dir)
 *               2 (x = flow dir, z = grad dir)
 *  shear_rate :
 * OUTPUT
 *  sys : struct stokes
 */
void
stokes_set_shear (struct stokes *sys,
		  int shear_mode,
		  double shear_rate);
/* get shear_shift at the given time t
 * INPUT
 *  sys : struct stokes
 *  t   : current time 
 *  t0  : time of the reference
 *  s0  : shift at t0
 * OUTPUT
 *  returned value : shift in the range [-lx/2, lx/2)
 */
double
stokes_get_shear_shift (struct stokes *sys,
			double t,
			double t0, double s0);
/* set shear_shift at the given time t
 * INPUT
 *  sys : struct stokes
 *  t   : current time 
 *  t0  : time of the reference
 *  s0  : shift at t0
 * OUTPUT
 *  sys->shear_shift : set and place in the range [-lx/2, lx/2)
 */
void
stokes_set_shear_shift (struct stokes *sys,
			double t,
			double t0, double s0);

void
stokes_set_l (struct stokes * sys,
	      double lx, double ly, double lz);

void
stokes_set_xi (struct stokes * sys,
	       double xi, double ewald_eps);

double
xi_by_tratio (struct stokes * sys,
	      double tratio);


/* set iter param
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *   out   : FILE * to output debug info.
 */
void
stokes_set_iter (struct stokes *sys,
		 const char *solver,
		 int max,
		 int restart,
		 double eps,
		 int debug,
		 FILE *out);

/* set pos for all particles safely by another array
 * INPUT
 *  pos[np*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < np*3)
 */
void
stokes_set_pos (struct stokes *sys,
		const double *pos);

/* set pos for mobile particles safely by another array
 * INPUT
 *  pos[nm*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < nm*3)
 */
void
stokes_set_pos_mobile (struct stokes *sys,
		       const double *pos);

/* set pos for fixed particles safely by another array
 * INPUT
 *  pos[nf*3] : only fixed particles are set, where nf = np - nm
 * OUTPUT
 *  sys->pos[i] for (i = nm*3; i < np*3)
 */
void
stokes_set_pos_fixed (struct stokes *sys,
		      const double *pos);

/* set radius (sys->a[], sys->twobody_f_list, and sys->poly_table).
 * Note that the default setting (sys->a == NULL) is for monodisperse system
 * where a=1 for all particles
 * INPUT
 *  a[np] :
 *  sys->twobody_nmax : define sys->twobody_nmax before calling.
 * OUTPUT
 *  sys->a[np]             :
 *  sys->poly_table[np*np] :
 *  sys->twobody_f_list[]  :
 */
void
stokes_set_radius (struct stokes *sys,
		   const double *a);
/* unset radius (sys->a[], sys->twobody_f_list, and sys->poly_table).
 * that is, the system is treated as a monodisperse system
 * where a=1 for all particles as in the default setting.
 * INPUT
 *  sys                    : struct stokes
 * OUTPUT
 *  sys->a[np]             : freed and set NULL
 *  sys->poly_table[np*np] : freed and set NULL
 *  sys->twobody_f_list[]  : freed and set NULL
 */
void
stokes_unset_radius (struct stokes *sys);


/** slip parameters **/
/* set slip parameters (slip[], slip_a[], slip_G32[], slip_G30[], slip_G52[])
 * Note that the default setting (sys->slip == NULL) is for no-slip system
 * where gamma=0 for all particles
 * INPUT
 *  gamma[np] : slip length
 * OUTPUT
 *  sys->slip[np]     : slip length
 *  sys->slip_a[np]   : effective radius for the laplacian terms
 *  sys->slip_G32[np] : Lambda(3,2) = 1/Lambda(2,3) for a-self.
 *  sys->slip_G30[np] : Lambda(3,0) = 1/Lambda(0,3) for c-self.
 *  sys->slip_G52[np] : Lambda(5,2) = 1/Lambda(2,5) for m-self.
 */
void
stokes_set_slip (struct stokes *sys,
		 const double *gamma);

/* unset slip params (slip[], slip_a[], slip_G32[], slip_G30[], slip_G52[])
 * that is, the system is treated as a no-slip system
 * where gamma=0 for all particles as in the default setting.
 * INPUT
 *  sys                    : struct stokes
 * OUTPUT
 *  sys->slip[np]     : freed and set NULL
 *  sys->slip_a[np]   : freed and set NULL
 *  sys->slip_G32[np] : freed and set NULL
 *  sys->slip_G30[np] : freed and set NULL
 *  sys->slip_G52[np] : freed and set NULL
 */
void
stokes_unset_slip (struct stokes *sys);


/* make a copy of struct stokes s0
 */
struct stokes *
stokes_copy (struct stokes *s0);



/******************/
/* from twobody.h */
/******************/
struct twobody_f {
  int nmax;
  double lambda;

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

struct twobody_f_list {
  int n;
  double *l; // lambda
  struct twobody_f **f;
};


void twobody_XA (int n, double l, double * f);
void twobody_YA (int n, double l, double * f);

void twobody_YB (int n, double l, double * f);

void twobody_XC (int n, double l, double * f);
void twobody_YC (int n, double l, double * f);

void twobody_XG (int n, double l, double * f);
void twobody_YG (int n, double l, double * f);

void twobody_YH (int n, double l, double * f);

void twobody_XM (int n, double l, double * f);
void twobody_YM (int n, double l, double * f);
void twobody_ZM (int n, double l, double * f);

void twobody_XP (int n, double l, double * f);
void twobody_XQ (int n, double l, double * f);


/** utility routines for struct twobody_f and twobody_f_list **/

struct twobody_f *
twobody_f_init (int nmax, double lambda);
void
twobody_f_free (struct twobody_f *f);

struct twobody_f_list *
twobody_f_list_init (void);
void
twobody_f_list_append (struct twobody_f_list *list,
		       int nmax, double lambda);
void
twobody_f_list_free (struct twobody_f_list *list);


/** far form **/

/* calc XA11 and XA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_far (int n, double l, double s,
		     double *XA11, double *XA12);

/* calc YA11 and YA12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_far (int n, double l, double s,
		     double *YA11, double *YA12);

/* calc YB11 and YB12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_far (int n, double l, double s,
		     double *YB11, double *YB12);

/* calc XC11 and XC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_far (int n, double l, double s,
		     double *XC11, double *XC12);

/* calc YC11 and YC12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_far (int n, double l, double s,
		     double *YC11, double *YC12);

/* calc XG11 and XG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XG11
 *  XG12
 */
void twobody_XG_far (int n, double l, double s,
		     double *XG11, double *XG12);

/* calc YG11 and YG12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YG11
 *  YG12
 */
void twobody_YG_far (int n, double l, double s,
		     double *YG11, double *YG12);

/* calc YH11 and YH12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YH11
 *  YH12
 */
void twobody_YH_far (int n, double l, double s,
		     double *YH11, double *YH12);

/* calc XM11 and XM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XM11
 *  XM12
 */
void twobody_XM_far (int n, double l, double s,
		     double *XM11, double *XM12);

/* calc YM11 and YM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YM11
 *  YM12
 */
void twobody_YM_far (int n, double l, double s,
		     double *YM11, double *YM12);

/* calc ZM11 and ZM12
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  ZM11
 *  ZM12
 */
void twobody_ZM_far (int n, double l, double s,
		     double *ZM11, double *ZM12);

/* calc scalar functions of resistance problem by 1/s expansion
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
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
void twobody_far (int version, int n, double l, double s,
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
void twobody_far_with_f (int version,
			 struct twobody_f *f12,
			 int n, double l, double s,
			 double *far);


/** lubrication form **/

/* calc XA11 and XA12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XA11
 *  XA12
 */
void twobody_XA_lub (int n, double l, double s,
		     double *XA11, double *XA12);

/* calc YA11 and YA12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YA11
 *  YA12
 */
void twobody_YA_lub (int n, double l, double s,
		     double *YA11, double *YA12);

/* calc YB11 and YB12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YB11
 *  YB12
 */
void twobody_YB_lub (int n, double l, double s,
		     double *YB11, double *YB12);

/* calc XC11 and XC12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  XC11
 *  XC12
 */
void twobody_XC_lub (int n, double l, double s,
		     double *XC11, double *XC12);

/* calc YC11 and YC12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YC_lub (int n, double l, double s,
		     double *YC11, double *YC12);

/* calc XG11 and XG12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_XG_lub (int n, double l, double s,
		     double *XG11, double *XG12);

/* calc YG11 and YG12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YG_lub (int n, double l, double s,
		     double *YG11, double *YG12);

/* calc YH11 and YH12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YH_lub (int n, double l, double s,
		     double *YH11, double *YH12);

/* calc XM11 and XM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_XM_lub (int n, double l, double s,
		     double *XM11, double *XM12);

/* calc YM11 and YM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_YM_lub (int n, double l, double s,
		     double *YM11, double *YM12);

/* calc ZM11 and ZM12 for lubrication form
 * INPUT
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
 * OUTPUT
 *  YC11
 *  YC12
 */
void twobody_ZM_lub (int n, double l, double s,
		     double *ZM11, double *ZM12);


/* calc scalar functions of resistance problem by lub form
 * INPUT
 *  version : 0=F, 1=FT, 2=FTS.
 *  n : max order
 *  l := a2 / a1
 *  s := 2 * r / (a1 + a2)
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
void twobody_lub (int version, int n, double l, double s,
		  double *lub);

/* calc scalar functions of resistance problem by lub form
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
void twobody_lub_with_f (int version,
			 struct twobody_f *f2b,
			 int n, double l, double s,
			 double *lub);

/* scale the scalar functions from Jeffrey-Onishi to Stokesian dynamics
 * INPUT
 *  version  : 0=F, 1=FT, 2=FTS.
 *  two [22] : scalar functions in Jeffrey form
 *  l        : lambda = ab / aa,
 *             where aa and ab are radii for particles a(alpha) and b(beta)
 *             Note that the scalar functions are for "a-b" interaction.
 * OUTPUT
 *  two [22] : scalar functions in SD form
 */
void
twobody_scale_SD (int version, double *two, double l);

/* scale the scalar functions from Jeffrey-Onishi to the dimensional form
 * INPUT
 *  version    : 0=F, 1=FT, 2=FTS.
 *  two [22] : scalar functions in Jeffrey form
 *  l        : lambda = ab / aa,
 *             where aa and ab are radii for particles a(alpha) and b(beta)
 *             Note that the scalar functions are for "a-b" interaction.
 * OUTPUT
 *  two [22] : scalar functions in the dimensional form
 */
void
twobody_scale (int version, double *two, double a1, double l);


/* calc scalar functions of two-body exact solution in resistance problem
 * INPUT
 *  version    : 0=F, 1=FT, 2=FTS.
 *  r          : distance between the two := x_b - x_a
 *  aa, ab     : radii for particles a(alpha) and b(beta)
 *  f12        : (struct twobody_f *).
 *               you can give NULL for them.
 *               then, the coefs are calculated on-the-fly (terribly slow).
 *  n          : max order for the coefficients
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 *  flag_scale : 0 no scaling, that is, in Jeffrey form
 *               1 for the dimensional form
 *               2 for the Stokesian dynamics form
 *  res [22]   :
 * OUTPUT
 *  res [22]   : scalar functions. the scaling is given by flag_scale.
 */
void
twobody_scalars_res (int version,
		     double r,
		     double aa, double ab,
		     struct twobody_f *f12,
		     int n, int flag_lub, int flag_scale,
		     double *res);


/***********************/
/* from twobody-slip.h */
/***********************/
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

//#include "stokes.h" // struct stokes
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
calc_lub_f_2b_slip
(struct stokes *sys,
 const double *u1, const double *u2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *f1, double *f2);


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
calc_lub_ft_2b_slip
(struct stokes *sys,
 const double *uo1, const double *uo2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *ft1, double *ft2);


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
calc_lub_fts_2b_slip
(struct stokes *sys,
 const double *uoe1, const double *uoe2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *fts1, double *fts2);




/************************************
 ** miscellaneous routines         **
 ************************************/
/* from bench.h */
/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
long
ptime_ms (void);

/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
double
ptime_ms_d (void);

/* return the current process time in mili-seconds
 * NOTE: not the difference of the times */
long
ptime_micros (void);


/* from matrix.h */

/* D = a * A + b * B . C
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  C [nc1, nc2]
 *  a
 *  b
 *  where (nb2 == nc1) && (na1 = nb1) && (na2 == nc2)
 * OUTPUT
 *  D [na1, na2] = a * A [na1, na2] + b * B [nb1, nb2] . C [nc1, nc2]
 */
void
add_and_mul (const double *A, int na1, int na2,
	     const double *B, int nb1, int nb2,
	     const double *C, int nc1, int nc2,
	     double a, double b,
	     double *D);

/* solve generalized linear set of equations using LU-decomposition
 * INPUT
 *  n1, n2 : dimension
 *  A [n1 * n1] :
 *  B [n1 * n2] :
 *  C [n2 * n1] :
 *  D [n2 * n2] :
 *
 *  E [n1 * n1] :
 *  F [n1 * n2] :
 *  G [n2 * n1] :
 *  H [n2 * n2] :
 *  where the generalized linear set of equations is
 *   [A B](x) = [E F](b)
 *   [C D](y)   [G H](c)
 * OUTPUT
 *  I [n1 * n1] :
 *  J [n1 * n2] :
 *  K [n2 * n1] :
 *  L [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (x) = [I J](b)
 *   (c)   [K L](y)
 *  note that A-D, E-H are destroyed!
 */
void
solve_gen_linear (int n1, int n2,
		  double * A, double * B, double * C, double * D,
		  double * E, double * F, double * G, double * H,
		  double * I, double * J, double * K, double * L);

/* solve linear set of equations using LU-decomposition
 * INPUT
 *  n1, n2 : dimension
 *  A [n1 * n1] :
 *  B [n1 * n2] :
 *  C [n2 * n1] :
 *  D [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (b) = [A B](x)
 *   (c)   [C D](y)
 * OUTPUT
 *  I [n1 * n1] :
 *  J [n1 * n2] :
 *  K [n2 * n1] :
 *  L [n2 * n2] :
 *  where the generalized linear set of equations is
 *   (x) = [I J](b)
 *   (c)   [K L](y)
 *  note that A-D are destroyed!
 */
void
solve_linear (int n1, int n2,
	      double * A, double * B, double * C, double * D,
	      double * I, double * J, double * K, double * L);

/* multiply two matrices (a wrapper to BLAS routine)
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  where (na2 == nb1)
 * OUTPUT
 *  C [na1, nb2] = A [na1, na2] . B [nb1, nb2]
 */
void
mul_matrices (const double * A, int na1, int na2,
	      const double * B, int nb1, int nb2,
	      double * C);
/*
 * INPUT
 *  mat [n1, n2]
 *  x [n2]
 * OUTPUT
 *  y [n1] = mat [n1, n2] . x [n2]
 */
void
dot_prod_matrix (const double * mat, int n1, int n2,
		 const double * x,
		 double * y);
/* 
 * INPUT
 *  alpha
 *  mat [n1, n2]
 *  x [n2]
 *  beta
 * OUTPUT
 *  y [n1] = alpha * mat [n1, n2] . x [n2] + beta * y[]
 */
void
dot_prod_matrix_ (double alpha, const double *mat, int n1, int n2,
		  const double *x,
		  double beta, double *y);

/* utility routine for matrix in the extracted form
 * INPUT
 *  np : # particles (not # elements!)
 *  m [np *11 * np *11] : matrix in the extracted form
 *  x [np *11] : vector in the extracted form
 * INPUT
 *  y [np *11] : output vector in the extracted form (:= m.x)
 */
void
multiply_extmat_with_extvec_3fts (int np, const double * m, const double * x,
				  double * y);


/* from dgetri_c.h */

void lapack_inv (int n, const double *a,
		 double *ai);

/* the version that a[n*n] is input AND output
 */
void lapack_inv_ (int n, double *a);

/* just solve one problem A.x = b
 * INPUT
 *  n : the order of the matrix A
 *  a[n*n] : coefficient matrix
 *  b[n]   : given vector
 * OUTPUT
 *  x[n]   : the solution
 */
void lapack_solve_lin (int n, const double *a, const double *b,
		       double *x);


/* from ewald-3f.h */
/** natural resistance problem **/
/* solve natural resistance problem in F version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f_0
(struct stokes *sys,
 const double *u,
 double *f);
/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f
(struct stokes *sys,
 const double *u,
 double *f);


/** natural mobility problem **/
/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 * OUTPUT
 *  u [np * 3] :
 */
void
solve_mob_3f
(struct stokes * sys,
 const double *f,
 double *u);


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_3f
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *u,
 double *ff);


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in F version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_lub_3f_0
(struct stokes * sys,
 const double *u,
 double *f);
/* solve natural resistance problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_lub_3f
(struct stokes * sys,
 const double *u,
 double *f);


/** mob_lub_3f **/
void
atimes_mob_lub_3f
(int n, const double *x,
 double *y, void *user_data);
/* solve natural mobility problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_lub_3f
(struct stokes * sys,
 const double *f,
 double *u);


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_lub_3f
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *u,
 double *ff);

/* from ewald-3ft.h */
/** natural resistance problem **/
/* solve natural resistance problem in FT version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft_0
(struct stokes * sys,
 const double *u, const double *o,
 double *f, double *t);
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft
(struct stokes * sys,
 const double *u, const double *o,
 double *f, double *t);


/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft
(struct stokes * sys,
 const double *f, const double *t,
 double *u, double *o);


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
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
solve_mix_3ft
(struct stokes * sys,
 const double *f, const double *t,
 const double *uf, const double *of,
 double *u, double *o,
 double *ff, double *tf);


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FT version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft_0
(struct stokes * sys,
 const double *u, const double *o,
 double *f, double *t);
/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft
(struct stokes * sys,
 const double *u, const double *o,
 double *f, double *t);


/** mob_lub_3ft **/
/* solve natural mobility problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft
(struct stokes * sys,
 const double *f, const double *t,
 double *u, double *o);


/** natural mobility problem with lubrication with fixed particles **/
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
solve_mix_lub_3ft
(struct stokes * sys,
 const double *f, const double *t,
 const double *uf, const double *of,
 double *u, double *o,
 double *ff, double *tf);

/* from ewald-3fts.h */
/** natural resistance problem **/
/* solve natural resistance problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts_0
(struct stokes * sys,
 const double *u, const double *o, const double *e,
 double *f, double *t, double *s);
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_3fts
(struct stokes * sys,
 const double *u, const double *o, const double *e,
 double *f, double *t, double *s);


/** natural mobility problem **/
/* solve natural mobility problem in FTS version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys  : system parameters
 *  iter : struct iter (if NULL is given, use sys->it for the solver)
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_3fts_0
(struct stokes *sys, struct iter *iter,
 const double *f, const double *t, const double *e,
 double *u, double *o, double *s);
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : strain tensor     in the labo frame.
 * OUTPUT
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 *  s [np * 5] :
 */
void
solve_mob_3fts
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 double *u, double *o, double *s);


/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  iter : struct iter (if NULL is given, use sys->it for the solver)
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : in the fluid-rest frame
 *  of [nf * 3] : in the fluid-rest frame
 *  ef [nf * 5] : in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts_0
(struct stokes *sys, struct iter *iter,
 const double *f, const double *t, const double *e,
 const double *uf, const double *of, const double *ef,
 double *u, double *o, double *s,
 double *ff, double *tf, double *sf);
/* solve natural mobility problem with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 const double *uf, const double *of, const double *ef,
 double *u, double *o, double *s,
 double *ff, double *tf, double *sf);


/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 *  s [np * 5] :
 */
void
solve_res_lub_3fts_0
(struct stokes * sys,
 const double *u, const double *o, const double *e,
 double *f, double *t, double *s);
/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] : in the labo frame.
 *   o [np * 3] : in the labo frame.
 *   e [np * 5] : in the labo frame.
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_lub_3fts
(struct stokes * sys,
 const double *u, const double *o, const double *e,
 double *f, double *t, double *s);


/** mob_lub_3fts **/
/* solve natural mobility problem with lubrication in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [np * 5] :
 */
void
solve_mob_lub_3fts_0
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 double *u, double *o, double *s);
/* solve natural mobility problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] : in the labo frame.
 * OUTPUT
 *   u [np * 3] : in the labo frame.
 *   o [np * 3] : in the labo frame.
 *   s [np * 5] :
 */
void
solve_mob_lub_3fts
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 double *u, double *o, double *s);


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [np * 5] : = E - E^inf, that is, in the fluid-rest frame
 *  uf [nf * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  of [nf * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  ef [nf * 5] : = E - E^inf, that is, in the fluid-rest frame
 * OUTPUT
 *  u [np * 3] : = U - u^inf, that is, in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, that is, in the fluid-rest frame
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_lub_3fts_0
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 const double *uf, const double *of,
 const double *ef,
 double *u, double *o, double *s,
 double *ff, double *tf, double *sf);
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_lub_3fts
(struct stokes * sys,
 const double *f, const double *t, const double *e,
 const double *uf, const double *of,
 const double *ef,
 double *u, double *o, double *s,
 double *ff, double *tf, double *sf);




/* from ewald-2f.h */
/** natural resistance problem **/
/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_2f (struct stokes * sys,
	      const double *u,
	      double *f);

/** natural mobility problem **/
/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
solve_mob_2f (struct stokes * sys,
	      const double *f,
	      double *u);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_2f (struct stokes * sys,
	      const double *f,
	      const double *uf,
	      double *u,
	      double *ff);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
solve_res_lub_2f (struct stokes * sys,
		  const double *u,
		  double *f);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
solve_mix_lub_2f (struct stokes * sys,
		  const double *f,
		  const double *uf,
		  double *u,
		  double *ff);

/* from ewald-2ft.h */
/** natural resistance problem **/
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_2ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t);

/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 */
void
solve_mob_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       double *u, double *o);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_lub_2ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_lub_2ft (struct stokes * sys,
		   const double *f, const double *t3,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf);

/* from ewald-2fts.h */
/** natural resistance problem **/
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_2fts (struct stokes * sys,
		const double *u, const double *o, const double *e,
		double *f, double *t, double *s);

/** natural mobility problem **/
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
solve_mob_2fts (struct stokes * sys,
		const double *f, const double *t3, const double *e,
		double *u, double *o, double *s);

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   e [nm * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   ef [nf * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_2fts (struct stokes * sys,
		const double *f, const double *t3, const double *e,
		const double *uf, const double *of, const double *ef,
		double *u, double *o, double *s,
		double *ff, double *tf, double *sf);

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
solve_res_lub_2fts (struct stokes * sys,
		    const double *u, const double *o, const double *e,
		    double *f, double *t, double *s);

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   e [nm * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   ef [nf * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
solve_mix_lub_2fts (struct stokes * sys,
		    const double *f, const double *t3,
		    const double *e,
		    const double *uf, const double *of,
		    const double *ef,
		    double *u, double *o, double *s,
		    double *ff, double *tf, double *sf);


/** no-hydrodynamic interaction problems **/
/* from noHI.h */

/* solve natural mobility problem with fixed particles in F version
 * without HI
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  ff [nf * 3] :
 */
void
solve_mix_3f_noHI (struct stokes *sys,
		   const double *f,
		   const double *uf,
		   double *u,
		   double *ff);

/* solve natural mobility problem with fixed particles in FT version
 * without HI
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  o [nm * 3] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 */
void
solve_mix_3ft_noHI (struct stokes *sys,
		    const double *f, const double *t,
		    const double *uf, const double *of,
		    double *u, double *o,
		    double *ff, double *tf);

/* solve natural mobility problem with fixed particles in FTS version
 * without HI
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts_noHI (struct stokes *sys,
		     const double *f, const double *t, const double *e,
		     const double *uf, const double *of, const double *ef,
		     double *u, double *o, double *s,
		     double *ff, double *tf, double *sf);


/** atimes routines **/
/* from ewald.h */

/*
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_ewald_real (int version,
		    double xi, double r,
		    double *xa, double *ya,
		    double *yb,
		    double *xc, double *yc,
		    double *xg, double *yg,
		    double *yh,
		    double *xm, double *ym, double *zm);

/* calculate scalar functions of (12)-interaction for unequal spheres
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_2 - x_1
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_ewald_real_poly (int version,
			 double xi, double r,
			 double aa, double ab,
			 double *xa, double *ya,
			 double *yb,
			 double *xc, double *yc,
			 double *xg, double *yg,
			 double *yh,
			 double *xm, double *ym, double *zm);


/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * this is a wrapper for non-periodic and periodic cases
 * also polydisperse systems for non-periodic
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_3all
(int n, const double *x, double *y, void * user_data);


/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all
(int n, const double *x, double *y, void * user_data);


/* from non-ewald.h */
/* calculate scalar functions under no periodic boundary condition
 * all functions are explicitly shown in Durlofsky-Brady (1987).
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 * OUTPUT
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald (int version,
		  double r,
		  double *scalar);

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald_poly (int version,
		       double r,
		       double aa, double ab,
		       double *scalar);

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [44] : scalar functions in dimensional form!
 *   0, 1, 2, 3 : (xa11, xa12, xa21, xa22)
 *   4, 5, 6, 7 : (ya11, ya12, ya21, ya22)
 *   8, 9,10,11 : (yb11, yb12, yb21, yb22)
 *  12,13,14,15 : (xc11, xc12, xc21, xc22)
 *  16,17,18,19 : (yc11, yc12, yc21, yc22)
 *  20,21,22,23 : (xg11, xg12, xg21, xg22)
 *  24,25,26,27 : (yg11, yg12, yg21, yg22)
 *  28,29,30,31 : (yh11, yh12, yh21, yh22)
 *  32,33,34,35 : (xm11, xm12, xm21, xm22)
 *  36,37,38,39 : (ym11, ym12, ym21, ym22)
 *  40,41,42,43 : (zm11, zm12, zm21, zm22)
 */
void
scalars_nonewald_poly_full (int version,
			    double r,
			    double aa, double ab,
			    double *scalar);


/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all
(int n, const double *x, double *y, void *user_data);


/* from lub.h */
/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 *  lubmax2 : square of the max distance (0 means no limit)
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
int
cond_lub (const double *x1, const double *x2, double lubmax2);

/* condition for lubrication for polydisperse system
 * INPUT
 *  x1 [3], x2 [3] : position
 *  a1, a2         : radii for particles 1 and 2
 *  lubmax2        : square of the max distance (0 means no limit)
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
int
cond_lub_poly (const double *x1, const double *x2,
	       double a1, double a2,
	       double lubmax2);

/* calculate lubrication f by u for all particles
 * for both under the periodic and non-periodic boundary conditions.
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
void
calc_lub_3f
(struct stokes *sys,
 const double *u,
 double *f);

/* calculate lubrication ft by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
void
calc_lub_3ft
(struct stokes * sys,
 const double * uo, double * ft);

/* calculate lubrication fts by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_3fts
(struct stokes * sys,
 const double * uoe, double * fts);




/* from f.h */
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
matrix_f_atimes (const double *x,
		 double *y,
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
	    const double *f2);

/* calc scalar functions of (M^inf)^-1 in F
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *  scalar_f [4] :
 */
void
scalar_minv_f (double s, double * scalar_f);

/* calculate f by u for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin is used.
 *   u1 [3] : velocity
 *   u2 [3] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   f1 [3] : force
 *   f2 [3] :
 */
void
calc_lub_f_2b (struct stokes * sys,
	       const double * u1, const double * u2,
	       const double * x1, const double * x2,
	       double * f1, double * f2);

/* calculate lub-matrix in F version for pair of particles 1 and 2
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
 *            (i,j) are used to assign the results in mat[].
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   n : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b (struct stokes * sys,
		 int i, int j,
		 const double *x1, const double *x2,
		 int n, double * mat);


/* calculate f by u for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin       : min distance for lub calculation.
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
calc_lub_f_2b_poly
(struct stokes *sys,
 const double *u1, const double *u2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *f1, double *f2);




/* pre-process for imposed flow shifting, that is, converting U
 * from the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * to the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  u[np*3] : velocity in the fluid-rest frame
 *            (data is preserved)
 * OUTPUT
 *  u0[np*3] : velocity in the labo frame
 */
void
shift_labo_to_rest_U (struct stokes * sys,
		      int np, const double *u,
		      double *u0);

/* post-process for imposed flow shifting, that is, converting U
 * from the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * to the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  u[np*3] : velocity in the fluid-rest frame
 *            (data is overwritten after the process)
 * OUTPUT
 *  u[np*3] : velocity in the labo frame
 */
void
shift_rest_to_labo_U (struct stokes * sys,
		      int np, double *u);


/* from ft.h */
/* store matrix in FT format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n6 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n6 * n6] :
 */
void
matrix_ft_ij (int i, int j,
	      double ex, double ey, double ez,
	      double xa, double ya,
	      double yb,
	      double xc, double yc,
	      int n6, double *mat);

/* store matrix in B part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   yb : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_B (int n, double *mat,
	     double ex, double ey, double ez,
	     double yb);

/* store matrix in C part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xc, yc : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_C (int n, double *mat,
	     double ex, double ey, double ez,
	     double xc, double yc);

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FT format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [6] : FTS of particle 'i'
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [6] : UOE of particle 'j'
 */
void
matrix_ft_atimes (const double *x,
		  double *y,
		  double ex, double ey, double ez,
		  double xa, double ya,
		  double yb,
		  double xc, double yc);

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FT format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- alpha(i)' interaction is stored.
 * INPUT
 *   x [6] : FTS of particle 'i'
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [6] : UOE of particle 'j'
 */
void
matrix_ft_self_atimes (const double *x,
		       double *y,
		       double ex, double ey, double ez,
		       double xa, double ya,
		       double yb,
		       double xc, double yc);

/* convert ft[] to f[], t[] (this is applicable for UO)
 * INPUT
 *  n : # particles
 *  ft [n * 6] :
 * OUTPUT
 *  f[n * 3] :
 *  t[n * 3] :
 */
void
set_FT_by_ft (int n,
	      double *f, double *t,
	      const double *ft);

/* convert ft[] to f[], t[] (this is applicable for UO)
 * INPUT
 *  n : # particles
 *  f[n * 3] :
 *  t[n * 3] :
 * OUTPUT
 *  ft [n * 6] :
 */
void
set_ft_by_FT (int n,
	      double *ft,
	      const double *f, const double *t);

/* calc scalar functions of (M^inf)^-1 in FT
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *  scalar_ft [10] :
 */
void
scalar_minv_ft (double s, double * scalar_ft);

/* calculate ft by uoe for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin is used.
 *   uo1 [6] : velocity, angular velocity
 *   uo2 [6] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   ft1 [6] : force, torque
 *   ft2 [6] :
 */
void
calc_lub_ft_2b (struct stokes * sys,
		const double *uo1, const double *uo2,
		const double *x1, const double *x2,
		double *ft1, double *ft2);

/* calculate lub-matrix in FT version for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin is used.
 *   i : particle index for '1'
 *   j : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   n : dimension of matrix 'mat' (must be np*6)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_ft_2b (struct stokes * sys,
		  int i, int j,
		  const double *x1, const double *x2,
		  int n, double * mat);

/* calculate ft by uo for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin       : min distance for lub calculation.
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
calc_lub_ft_2b_poly
(struct stokes *sys,
 const double *uo1, const double *uo2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *ft1, double *ft2);




/* pre-process for imposed flow shifting, that is, converting U,O
 * from the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * to the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  o[np*3] : angular velocity in the fluid-rest frame
 *            (data is preserved)
 * OUTPUT
 *  o0[np*3] : angular velocity in the labo frame
 */
void
shift_labo_to_rest_O (struct stokes * sys,
		      int np, const double *o,
		      double *o0);

/* post-process for imposed flow shifting, that is, converting U,O
 * from the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * to the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  o[np*3] : angular velocity in the fluid-rest frame
 *            (data is overwritten after the process)
 * OUTPUT
 *  o[np*3] : angular velocity in the labo frame
 */
void
shift_rest_to_labo_O (struct stokes * sys,
		      int np, double *o);


/* from fts.h */
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
calc_lub_fts_2b_poly
(struct stokes *sys,
 const double *uoe1, const double *uoe2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *fts1, double *fts2);




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




/* from minv-poly.h */
/* calc scalar functions of (M^inf)^-1 in F for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_f [8] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 */
void
scalars_minv_f_poly (double r, double aa, double ab,
		     double *scalar_f);

/* calc scalar functions of (M^inf)^-1 in FT for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_ft [20] : scalar functions in dimensional form!
 *      0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *      4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *      8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *     12,13,14,15 : (XC11, XC12, XC21, XC22)
 *     16,17,18,19 : (YC11, YC12, YC21, YC22)
 */
void
scalars_minv_ft_poly (double r, double aa, double ab,
		      double *scalar_ft);

/* calc scalar functions of (M^inf)^-1 in FTS for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_fts [44] : scalar functions in dimensional form!
 *       0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *       4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *       8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *      12,13,14,15 : (XC11, XC12, XC21, XC22)
 *      16,17,18,19 : (YC11, YC12, YC21, YC22)
 *      20,21,22,23 : (XG11, XG12, XG21, XG22)
 *      24,25,26,27 : (YG11, YG12, YG21, YG22)
 *      28,29,30,31 : (YH11, YH12, YH21, YH22)
 *      32,33,34,35 : (XM11, XM12, XM21, XM22)
 *      36,37,38,39 : (YM11, YM12, YM21, YM22)
 *      40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 */
void
scalars_minv_fts_poly (double r, double aa, double ab,
		       double *scalar_fts);




/** lubrication functions for polydisperse systems **/

//#include "twobody.h" // struct twobody_f

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12    : (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 * OUTPUT
 *  lub [22] : scalar functions in dimensional form!
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
void
scalars_lub_poly (int version,
		  double r, double a1, double a2,
		  struct twobody_f *f12,
		  int n, int flag_lub,
		  double *lub);

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12,f21: (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
 *  n : max order
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
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
scalars_lub_poly_full (int version,
		       double r, double a1, double a2,
		       struct twobody_f *f12, struct twobody_f *f21,
		       int n, int flag_lub,
		       double *lub);

/* from two-body-res.h */
/* calc scalar functions of two-body resistance
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   res [22] : scalar functions
 */
void
scalar_two_body_res (double s, double *res);




/* the following header files are internal use only */

// dgemm_c.h
// dgemv_c.h
// memory-check.h


#endif /* !_LIBSTOKES_CORE_H_ */