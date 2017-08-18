/* fast semi-implicit Brownian dynamics algorithms
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // memcpy()
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_poly.h>
#include <nitsol_c.h> // struct NITSOL
#include <KIrand.h> // struct KIrand

#include <libstokes-core.h>

#include <brownian.h>  // struct BD_params
#include <bd-imp.h>    // struct BD_imp
#include <bonds.h>     // struct BONDS
#include <bonds-groups.h> // struct BONDS_GROUPS

#include "memory-check.h" // CHECK_MALLOC

#include "bd-imp-fast.h"


/* calc imposed flow difference for a spring at the mid-point
 * basically a copy from shift_rest_to_labo_U() in f.c.
 * INPUT
 *  b : struct BD_imp
 *  q[3] : connector vector at which the imposed flow difference is evaluated
 * OUTPUT
 *  u[3] : overwritten
 */
void
fastSI_calc_imposed (struct stokes *sys,
		     const double *q,
		     double *u)
{
  // O\times x
  u[0] = sys->Oi[1] * q[2] - sys->Oi[2] * q[1];
  u[1] = sys->Oi[2] * q[0] - sys->Oi[0] * q[2];
  u[2] = sys->Oi[0] * q[1] - sys->Oi[1] * q[0];

  // E.x
  double Ezz = - sys->Ei[0] - sys->Ei[4];
  u[0] += sys->Ei[0] * q[0] // xx . x
    +     sys->Ei[1] * q[1] // xy . y
    +     sys->Ei[2] * q[2];// xz . z
  u[1] += sys->Ei[1] * q[0] // yx . x
    +     sys->Ei[4] * q[1] // yy . y
    +     sys->Ei[3] * q[2];// yz . z
  u[2] += sys->Ei[2] * q[0] // zx . x
    +     sys->Ei[3] * q[1] // zy . y
    +     Ezz        * q[2];// zz . z
}

/* calc imposed flow difference for a spring at the mid-point
 * basically a copy from shift_rest_to_labo_U() in f.c.
 * INPUT
 *  b : struct BD_imp
 *  q1[3] q2[3] : two connector vectors
 *                the imposed flow difference is evaluated at the mid-point
 *                (1/2)(q1 + q2).
 * OUTPUT
 *  u[3] : overwritten
 */
void
fastSI_calc_imposed_midpoint (struct stokes *sys,
			      const double *q1,
			      const double *q2,
			      double *u)
{
  // calc mid-point
  double q[3];
  q[0] = 0.5 * (q1[0] + q2[0]);
  q[1] = 0.5 * (q1[1] + q2[1]);
  q[2] = 0.5 * (q1[2] + q2[2]);

  fastSI_calc_imposed (sys, q, u);
}


// borrowed from atimes_nonewald_3all() in non-ewald.c
/* calc U_i += M(i,j) . F_j
 * INPUT
 *  sys : struct BD_imp
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  i, j : particle indices
 *  f[3] : force on particle j
 * OUTPUT
 *  u[3] := M(i,j).f[], the velocity of particle i
 */
static void
mobility_F_atimes (struct stokes *sys, int flag_noHI,
		   int i, int j,
		   const double *f,
		   double *u)
{
  // zero clear
  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  double ai;
  double aj;
  if (sys->a == NULL)
    {
      ai = 1.0;
      aj = 1.0;
    }
  else
    {
      ai = sys->a[i];
      aj = sys->a[j];
    }

  if (i == j)
    {
      // self
      double self_a = 1.0;
      matrix_f_atimes (f, u,
		       0.0, 0.0, 0.0,
		       self_a, self_a); // a part
    }
  else if (flag_noHI == 0)
    {
      double mob [22];

      double xx
	= sys->pos [j*3  ]
	- sys->pos [i*3  ];
      double yy
	= sys->pos [j*3+1]
	- sys->pos [i*3+1];
      double zz
	= sys->pos [j*3+2]
	- sys->pos [i*3+2];
      double rr = xx * xx + yy * yy + zz * zz;

      double r = sqrt (rr);
      double rmin = (ai + aj) * sys->rmin;
      if (r < rmin) r = rmin;

      double ex = xx / r;
      double ey = yy / r;
      double ez = zz / r;

      // polydisperse code
      scalars_nonewald_poly (0, // F version
			     r,
			     ai, aj,
			     mob);
      scalars_mob_poly_scale_SD (0, // F version
				 ai,
				 mob);
      // now mob is in the SD form

      // note that interaction (i,j) should be for (U[i], F[j])
      matrix_f_atimes (f, u,
		       ex, ey, ez,
		       mob[0], mob[1]); // xa, ya
    }

  /* convert the resultant velocity into the Global scaling
   */
  u[0] /= ai;
  u[1] /= ai;
  u[2] /= ai;
}

double
calc_B (struct BD_imp *b, int ibond)
{
  double B;
  if (b->BD->sys->a == NULL)
    {
      B = 2.0;
    }
  else
    {
      int ia = b->BD->bonds->ia[ibond];
      int ib = b->BD->bonds->ib[ibond];
      B = 1.0 / b->BD->sys->a[ia]
	+ 1.0 / b->BD->sys->a[ib];
    }
  return (B);
}

/* calculate P = V_a - V_b, where V_a = sum_{j} M_{a,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  bonds : struct BONDS
 *  f[np] : forces on the particles
 *  ib : bond index considering
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar (struct stokes *sys, int flag_noHI,
	    struct BONDS *bonds,
	    const double *f,
	    int ib,
	    double *P)
{
  // zero clear
  P[0] = 0.0;
  P[1] = 0.0;
  P[2] = 0.0;

  int j;
  for (j = 0; j < sys->np; j ++)
    {
      double u[3];
      int j3 = j * 3;

      // from particle "j" to particle "ia[ib]"
      mobility_F_atimes (sys, flag_noHI, 
			 bonds->ia[ib], j,
			 f + j3, u);
      P[0] += u[0];
      P[1] += u[1];
      P[2] += u[2];

      // from particle "j" to particle "ib[ib]"
      mobility_F_atimes (sys, flag_noHI, 
			 bonds->ib[ib], j,
			 f + j3, u);
      P[0] -= u[0];
      P[1] -= u[1];
      P[2] -= u[2];
    }
}

/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[np*3] : initial connector vectors for imposed velocity
 *  q1[np*3] : updated connector vectors for imposed velocity
 *  q2[np*3] : the last-step connector vectors
 *  q [np*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ib       : the bond index for native bond (not COM), that is,
 *             ib = 0, 1, ..., ng-2. (COM is at ib=(ng-1).)
 *  iq       : index for connector divided by 3
 * OUTPUT
 *  r[3]     : the right-hand side vector for the cubic equation
 */
void
fastSI_rhs (struct BD_imp *b,
	    const double *q0,
	    const double *q1,
	    const double *q2,
	    const double *q,
	    int ib,
	    int iq,
	    double *r)
{
  int iq3 = iq * 3;
  fastSI_calc_imposed_midpoint (b->BD->sys,
				q0 + iq3,
				q1 + iq3,
				r);
  // r[] is set by the imposed flow

  double fs[3];
  BONDS_calc_force_spring_i (b->BD->bonds,
			     ib,
			     q2 + iq3,
			     fs);
  double B = calc_B (b, ib);
  r[0] += B * fs[0];
  r[1] += B * fs[1];
  r[2] += B * fs[2];


  // interaction part
  int np = b->BD->sys->np;
  int np3 = np * 3;
  double *qs  = (double *)calloc (np3, sizeof (double));
  double *f   = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (qs,  "fastSI_rhs");
  CHECK_MALLOC (f,   "fastSI_rhs");
  // all arrays are zero-cleared

  // form connector vectors
  int j;
  for (j = 0; j < iq; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q[j3  ];
      qs[j3+1] = q[j3+1];
      qs[j3+2] = q[j3+2];
    }
  for (; j < np; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q2[j3  ];
      qs[j3+1] = q2[j3+1];
      qs[j3+2] = q2[j3+2];
    }
  /*
  int n1 = iq3;
  int n2 = np3 - n1;
  qs = memcpy (qs,      q,       sizeof (double) * n1);
  qs = memcpy (qs + n1, q2 + n1, sizeof (double) * n2);
  */

  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     qs,
		     b->BD->sys->pos);

  // calc force for all particles
  BD_add_FP (b->BD, b->BD->sys->pos, f);

  // calculate P^* = V^*_a - V^*_b
  double P[3] = {0.0, 0.0, 0.0};
  calc_Pstar (b->BD->sys, b->BD->flag_noHI,
	      b->BD->bonds,
	      f,
	      ib,
	      P);

  r[0] += P[0];
  r[1] += P[1];
  r[2] += P[2];

  // house-keeping
  free (qs);
  free (f);

  
  /* **NOTE** b->z[] should be in the mobility mode
   *          and stored in the connector index.
   */
  r[0] = b->dt * r[0] + q0[iq3  ] + b->z[iq3  ];
  r[1] = b->dt * r[1] + q0[iq3+1] + b->z[iq3+1];
  r[2] = b->dt * r[2] + q0[iq3+2] + b->z[iq3+2];
}

/* calculate P^COM = (1/N)sum_i V_i, where V_i = sum_{j} M_{i,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  group : struct BONDS_GROUP for the considering group
 *  f[np] : forces on the particles
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar_COM (struct stokes *sys, int flag_noHI,
		struct BONDS_GROUP *group,
		const double *f,
		double *P)
{
  // zero clear
  P[0] = 0.0;
  P[1] = 0.0;
  P[2] = 0.0;

  int ng = group->np;
  int i;
  for (i = 0; i < ng; i ++)
    {
      int j;
      if (flag_noHI == 1)
	{
	  j = group->ip[i];
	  double u[3];
	  int j3 = j * 3;

	  // from particle "j" to particle "ip[i]"
	  mobility_F_atimes (sys, flag_noHI,
			     group->ip[i], j,
			     f + j3, u);
	  P[0] += u[0];
	  P[1] += u[1];
	  P[2] += u[2];
	}
      else
	{
	  for (j = 0; j < sys->np; j ++)
	    {
	      double u[3];
	      int j3 = j * 3;

	      // from particle "j" to particle "ip[i]"
	      mobility_F_atimes (sys, flag_noHI,
				 group->ip[i], j,
				 f + j3, u);
	      P[0] += u[0];
	      P[1] += u[1];
	      P[2] += u[2];
	    }
	}
    }
  P[0] /= (double)ng;
  P[1] /= (double)ng;
  P[2] /= (double)ng;
}
/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[ng*3] : initial connector vectors for imposed velocity
 *  q1[ng*3] : updated connector vectors for imposed velocity
 *  q2[ng*3] : the last-step connector vectors
 *  q [ng*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ig       : the group index
 *  iq       : index for connector divided by 3
 *             this should be the COM component
 * OUTPUT
 *  qCOM [3] : the updated COM
 */
void
fastSI_update_COM (struct BD_imp *b,
		   const double *q0,
		   const double *q1,
		   const double *q2,
		   const double *q,
		   int ig,
		   int iq,
		   double *qCOM)
{
  int iq3 = iq * 3;
  fastSI_calc_imposed_midpoint (b->BD->sys,
				q0 + iq3,
				q1 + iq3,
				qCOM);

  // interaction part
  int np = b->BD->sys->np;
  int np3 = np * 3;
  double *qs  = (double *)calloc (np3, sizeof (double));
  double *f   = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (qs,  "fastSI_rhs");
  CHECK_MALLOC (f,   "fastSI_rhs");
  // all arrays are zero-cleared

  // form connector vectors
  int j;
  for (j = 0; j < iq; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q[j3  ];
      qs[j3+1] = q[j3+1];
      qs[j3+2] = q[j3+2];
    }
  for (; j < np; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q2[j3  ];
      qs[j3+1] = q2[j3+1];
      qs[j3+2] = q2[j3+2];
    }
  /*
  int n1 = iq3;
  int n2 = np3 - n1;
  qs = memcpy (qs,      q,       sizeof (double) * n1);
  qs = memcpy (qs + n1, q2 + n1, sizeof (double) * n2);
  */

  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     qs,
		     b->BD->sys->pos);

  // calc force for all particles
  BD_add_FP (b->BD, b->BD->sys->pos, f);

  // calculate P^COM = (1/N) sum V
  double P[3] = {0.0, 0.0, 0.0};
  calc_Pstar_COM (b->BD->sys, b->BD->flag_noHI,
		  b->BD->groups->group[ig],
		  f,
		  P);

  qCOM[0] += P[0];
  qCOM[1] += P[1];
  qCOM[2] += P[2];

  // house-keeping
  free (qs);
  free (f);

  
  /* **NOTE** b->z[] should be in the mobility mode
   *          and stored in the connector index.
   */
  qCOM[0] = b->dt * qCOM[0] + q0[iq3  ] + b->z[iq3  ];
  qCOM[1] = b->dt * qCOM[1] + q0[iq3+1] + b->z[iq3+1];
  qCOM[2] = b->dt * qCOM[2] + q0[iq3+2] + b->z[iq3+2];
}

/* calculate V_i = sum_{j} M_{i,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  ip       : particle index of the isolated particle
 *  f[np] : forces on the particles
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar_isolated (struct stokes *sys, int flag_noHI,
		     int ip,
		     const double *f,
		     double *P)
{
  // zero clear
  P[0] = 0.0;
  P[1] = 0.0;
  P[2] = 0.0;

  int j;
  if (flag_noHI == 1)
    {
      j = ip;
      double u[3];
      int j3 = j * 3;

      // from particle "j" to particle "ip"
      mobility_F_atimes (sys, flag_noHI,
			 ip, j,
			 f + j3, u);
      P[0] += u[0];
      P[1] += u[1];
      P[2] += u[2];
    }
  else
    {
      for (j = 0; j < sys->np; j ++)
	{
	  double u[3];
	  int j3 = j * 3;

	  // from particle "j" to particle "ip"
	  mobility_F_atimes (sys, flag_noHI,
			     ip, j,
			     f + j3, u);
	  P[0] += u[0];
	  P[1] += u[1];
	  P[2] += u[2];
	}
    }
}
/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[ng*3] : initial connector vectors for imposed velocity
 *  q1[ng*3] : updated connector vectors for imposed velocity
 *  q2[ng*3] : the last-step connector vectors
 *  q [ng*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ip       : particle index of the isolated particle
 *  iq       : index for connector divided by 3
 *             this should be the component of isolated particle
 * OUTPUT
 *  qp [3]   : the updated connector of the isolated particle
 */
void
fastSI_update_isolated (struct BD_imp *b,
			const double *q0,
			const double *q1,
			const double *q2,
			const double *q,
			int ip,
			int iq,
			double *qp)
{
  int iq3 = iq * 3;
  fastSI_calc_imposed_midpoint (b->BD->sys,
				q0 + iq3,
				q1 + iq3,
				qp);

  // interaction part
  int np = b->BD->sys->np;
  int np3 = np * 3;
  double *qs  = (double *)calloc (np3, sizeof (double));
  double *f   = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (qs,  "fastSI_rhs");
  CHECK_MALLOC (f,   "fastSI_rhs");
  // all arrays are zero-cleared

  // form connector vectors
  int j;
  for (j = 0; j < iq; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q[j3  ];
      qs[j3+1] = q[j3+1];
      qs[j3+2] = q[j3+2];
    }
  for (; j < np; j ++)
    {
      int j3 = j * 3;
      qs[j3  ] = q2[j3  ];
      qs[j3+1] = q2[j3+1];
      qs[j3+2] = q2[j3+2];
    }
  /*
  int n1 = iq3;
  int n2 = np3 - n1;
  qs = memcpy (qs,      q,       sizeof (double) * n1);
  qs = memcpy (qs + n1, q2 + n1, sizeof (double) * n2);
  */

  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     qs,
		     b->BD->sys->pos);

  // calc force for all particles
  BD_add_FP (b->BD, b->BD->sys->pos, f);

  // calculate P_ip = V_ip
  double P[3] = {0.0, 0.0, 0.0};
  // from particle "j" to particle "ip"
  calc_Pstar_isolated (b->BD->sys, b->BD->flag_noHI,
		       ip,
		       f,
		       P);

  qp[0] += P[0];
  qp[1] += P[1];
  qp[2] += P[2];

  // house-keeping
  free (qs);
  free (f);

  
  /* **NOTE** b->z[] should be in the mobility mode
   *          and stored in the connector index.
   */
  qp[0] = b->dt * qp[0] + q0[iq3  ] + b->z[iq3  ];
  qp[1] = b->dt * qp[1] + q0[iq3+1] + b->z[iq3+1];
  qp[2] = b->dt * qp[2] + q0[iq3+2] + b->z[iq3+2];
}

double
fastSI_solve_cubic_in_range (double a, double b, double c, double d,
			     double xmin, double xmax)
{
  double x1, x2, x3;
  int num = gsl_poly_solve_cubic (b/a, c/a, d/a,
				  &x1, &x2, &x3);

  if (num == 3)
    {
      int flag = 0;
      int flag1 = 0;
      if (x1 >= xmin && x1 <= xmax)
	{
	  flag ++;
	  flag1 = 1;
	}
      int flag2 = 0;
      if (x2 >= xmin && x2 <= xmax)
	{
	  flag ++;
	  flag2 = 1;
	}
      int flag3 = 0;
      if (x3 >= xmin && x3 <= xmax)
	{
	  flag ++;
	  flag3 = 1;
	}

      if (flag == 0)
	{
	  fprintf (stdout,
		   "# fastSI_solve_cubic_in_range:"
		   " no solution in the range\n");
	  exit (1);
	}
      /*
      else if (flag > 1)
	{
	  fprintf (stdout,
		   "# fastSI_solve_cubic_in_range:"
		   " too many solutions (%d) in the range\n",
		   flag);
	  if (flag1 == 1 && flag2 == 1)
	    {
	      double d = fabs (x1 - x2);
	      fprintf (stdout,
		       "# fastSI_solve_cubic_in_range:"
		       " diff = %e\n", d);
	    }
	  else if (flag2 == 1 && flag3 == 1)
	    {
	      double d = fabs (x2 - x3);
	      fprintf (stdout,
		       "# fastSI_solve_cubic_in_range:"
		       " diff = %e\n", d);
	    }
	  else
	    {
	      double d = fabs (x3 - x1);
	      fprintf (stdout,
		       "# fastSI_solve_cubic_in_range:"
		       " diff = %e\n", d);
	    }
	  exit (1);
	}
      */
      else
	{
	  if (flag1 == 1)      return (x1);
	  else if (flag2 == 1) return (x2);
	  else                 return (x3);
	}
    }
  else if (num == 1)
    {
      if (x1 < xmin || x1 > xmax)
	{
	  fprintf (stdout, "# fastSI_solve_cubic_in_range:"
		   " no solution in the range\n");
	  exit (1);
	}
      return (x1);
    }
  else
    {
      fprintf (stdout, "# fastSI_solve_cubic_in_range:"
	       " num = %d != 1 or 3\n",
	       num);
      exit (1);
    }
}

/* solve the cubic equation
 * INPUT
 *  b  : struct BD_imp
 *  ib : bond index for (struct BONDS) b->BD->bonds
 *  R  : magnitude of the right-hand side vector
 */
double
fastSI_solve_cubic (struct BD_imp *b,
		    int ib,
		    double R)
{
  double Q0 = b->BD->bonds->p2[ib]; // = Ls / length
  double A  = b->BD->bonds->p1[ib];
  // B  := Asp * dt * (mob_{i} + mob_{i+1})
  double B = A * b->dt * calc_B (b, ib);
  

  double q;
  double aa, bb, cc, dd;
  //int status;

  int type = b->BD->bonds->type[ib];
  double s, s2;
  switch (type)
    {
    case 0: // Hookean (Fraenkel)
    case 6: // dWLC
      //q = (R + B) * Q0 / (Q0 + B);
      q = (R + Q0 * B) / (Q0 + B); // here A is replaced by A/Q0
      break;

    case 1: // WLC
      aa = 4.0 * B + 6.0 * Q0;
      bb = -3.0 * Q0 * (4.0 * Q0 + 2.0 * R + 3.0 * B);
      cc = 6.0 * Q0 * Q0 * (Q0 + 2.0 * R + B);
      dd = - 6.0 * R * Q0 * Q0 * Q0;
      q = fastSI_solve_cubic_in_range (aa, bb, cc, dd,
				       0.0, Q0);
      break;

    case 2: // ILC
      fprintf (stderr,
	       "# fast semi-implicit scheme for ILC is not implemented\n");
      exit (1);
      break;

    case 3: // Cohen
      aa = -(3.0 * Q0 + B);
      bb = 3.0 * R * Q0;
      cc = 3.0 * Q0 * Q0 * (Q0 + B);
      dd = -3.0 * R * Q0 * Q0 * Q0;
      q = fastSI_solve_cubic_in_range (aa, bb, cc, dd,
				       0.0, Q0);
      break;

    case 4: // Werner
      aa = 1.0;
      bb = -R;
      cc = -Q0 * (Q0 + B);
      dd = Q0 * Q0 * R;
      q = fastSI_solve_cubic_in_range (aa, bb, cc, dd,
				       0.0, Q0);
      break;

    case 5: // Hookean
      q = R * Q0 / (B + Q0);
      break;

    case 7: // FENE-Fraenkel
      s = b->BD->bonds->p3[ib]; // = s, tolerance
      s2 = s * s;

      aa = 1.0;
      bb = -2.0 * Q0 - R;
      cc = Q0 * (2.0 * R - B * s2 + Q0 * (1.0 - s2));
      dd = Q0 * Q0 * (R * (s2 - 1.0) + B * s2);
      double qmin = (1.0 - s) * Q0;
      double qmax = (1.0 + s) * Q0;
      q = fastSI_solve_cubic_in_range (aa, bb, cc, dd,
				       qmin, qmax);
      break;

    default:
      fprintf (stderr,
	       "# invalid spring type\n");
      exit (1);
      break;
    }

  return (q);
}

/* 
 * INPUT
 *  b    : struct BD_imp
 *  q0[] : initial connector
 *  qC[] : the predictor for the initial corrector step (flag_first == 1)
 *         the previous corrector for update steps (flag_first == 0)
 *  flag_first : 1 for the initial corrector
 *               0 for the update steps
 * OUTPUT
 *  returned value : residual (only for flag_first == 0)
 *  q[]  : updated corrector
 */
double
fastSI_solve_corrector (struct BD_imp *b,
			const double *q0,
			const double *qC,
			int flag_first,
			double *q)
{
  double eps = 0.0;

  // loop for groups
  int ig0 = 0; // head pointer of the group "ig"
  int ng = b->BD->groups->n;
  int ig;
  for (ig = 0; ig < ng; ig ++)
    {
      struct BONDS_GROUP *g = b->BD->groups->group[ig];

      if (g->np == 1)
	{
	  // isolated particle
	  int iq = ig0;
	  int iq3 = iq * 3;
	  double qp[3];
	  if (flag_first == 1) // first step
	    {
	      // here qC[] is the predictor
	      fastSI_update_isolated
		(b,
		 q0, qC, // for imposed flow
		 q0, q,  // for the interactions
		 g->ip[0], // particle index
		 iq, // index for connector (divided by 3)
		 qp);
	    }
	  else // update steps
	    {
	      fastSI_update_isolated
		(b,
		 q0, qC, // for imposed flow
		 qC, q,  // for the interactions
		 g->ip[0], // particle index
		 iq, // index for connector (divided by 3)
		 qp);
	    }

	  q[iq3  ] = qp[0];
	  q[iq3+1] = qp[1];
	  q[iq3+2] = qp[2];
	  /*
	  if (flag_first == 0) // update steps
	    {
	      double dx = q[iq3  ] - qC[iq3  ];
	      double dy = q[iq3+1] - qC[iq3+1];
	      double dz = q[iq3+2] - qC[iq3+2];
	      eps += dx*dx + dy*dy + dz*dz;
	    }
	  */
	}
      else
	{
	  // loop for bonds in group "ig"
	  int i;
	  for (i = 0; i < g->np - 1; i ++)
	    {
	      int ib = g->bonds[i];
	      int iq = ig0 + i;
	      int iq3 = iq * 3;
	      double r[3];
	      if (flag_first == 1) // first step
		{
		  // here qC[] is the predictor
		  fastSI_rhs
		    (b,
		     q0, qC, // for imposed flow
		     q0, q,  // for the interactions
		     ib, // bond index
		     iq, // index for connector (divided by 3)
		     r);
		}
	      else // update steps
		{
		  fastSI_rhs
		    (b,
		     q0, qC, // for imposed flow
		     qC, q,  // for the interactions
		     ib, // bond index
		     iq, // index for connector (divided by 3)
		     r);
		}
	      double r_norm = sqrt (r[0] * r[0]
				    + r[1] * r[1]
				    + r[2] * r[2]);
	      r[0] /= r_norm;
	      r[1] /= r_norm;
	      r[2] /= r_norm;
	      // now r[] is the unit vector
      
	      double q_norm
		= fastSI_solve_cubic (b,
				      ib, // bond index
				      r_norm);
	      q[iq3  ] = q_norm * r[0];
	      q[iq3+1] = q_norm * r[1];
	      q[iq3+2] = q_norm * r[2];

	      if (flag_first == 0) // update steps
		{
		  double dx = q[iq3  ] - qC[iq3  ];
		  double dy = q[iq3+1] - qC[iq3+1];
		  double dz = q[iq3+2] - qC[iq3+2];
		  eps += dx*dx + dy*dy + dz*dz;
		}
	    }
	  // the COM component
	  int iCOM  = ig0 + i;
	  int iCOM3 = iCOM * 3;
	  double qCOM[3];
	  if (flag_first == 1) // first step
	    {
	      // here qC[] is the predictor
	      fastSI_update_COM
		(b,
		 q0, qC, // for imposed flow
		 q0, q,  // for the interactions
		 ig, // group index
		 iCOM, // index for connector (divided by 3)
		 qCOM);
	    }
	  else // update steps
	    {
	      fastSI_update_COM
		(b,
		 q0, qC, // for imposed flow
		 qC, q,  // for the interactions
		 ig, // group index
		 iCOM, // index for connector (divided by 3)
		 qCOM);
	    }

	  q[iCOM3  ] = qCOM[0];
	  q[iCOM3+1] = qCOM[1];
	  q[iCOM3+2] = qCOM[2];
	  /*
	  if (flag_first == 0) // update steps
	    {
	      double dx = q[iCOM3  ] - qC[iCOM3  ];
	      double dy = q[iCOM3+1] - qC[iCOM3+1];
	      double dz = q[iCOM3+2] - qC[iCOM3+2];
	      eps += dx*dx + dy*dy + dz*dz;
	    }
	  */
	}
      ig0 += g->np;
    }


  return (sqrt (eps));
}

/* update connector q[] by Euler scheme
 * this is used to build predictor in fastSI algorithm
 * INPUT
 *  b : struct BD_imp
 *      **NOTE** b->z[] should be in the mobility mode
 *               and stored in the connector index.
 */
void
fastSI_solve_Euler (struct BD_imp *b,
		    const double *q0,
		    double *q)
{
  int np = b->BD->sys->np;
  int np3 = np * 3;

  double *f = (double *)calloc (np3, sizeof (double));
  double *u = (double *)calloc (np3, sizeof (double));
  double *P = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (f, "fastSI_solve_Euler");
  CHECK_MALLOC (u, "fastSI_solve_Euler");
  CHECK_MALLOC (P, "fastSI_solve_Euler");

  // calc force for all particles for the connector q0[]
  BONDS_conn_to_pos (b->BD->bonds,
		     b->BD->groups,
		     q0,
		     b->BD->sys->pos);
  BD_add_FP (b->BD, b->BD->sys->pos, f);

  if (b->BD->sys->version != 0)
    {
      fprintf (stderr, "# fastSI_solve_Euler: only F version is allowed\n");
      exit (1);
    }
  atimes_3all (np3, f, u, (void *)b->BD->sys);

  BONDS_pos_to_conn (b->BD->bonds,
		     b->BD->groups,
		     u,
		     P);
  free (u);


  // loop for groups
  int ig0 = 0; // head pointer of the group "ig"
  int ng = b->BD->groups->n;
  int ig;
  for (ig = 0; ig < ng; ig ++)
    {
      struct BONDS_GROUP *g = b->BD->groups->group[ig];

      if (g->np == 1)
	{
	  // isolated particle
	  int iq = ig0;
	  int iq3 = iq * 3;
	  fastSI_calc_imposed (b->BD->sys,
			       q0 + iq3,
			       q  + iq3);

	  q[iq3  ] += P[iq3  ];
	  q[iq3+1] += P[iq3+1];
	  q[iq3+2] += P[iq3+2];

	  q[iq3  ] = b->dt * q[iq3  ] + q0[iq3  ] + b->z[iq3  ];
	  q[iq3+1] = b->dt * q[iq3+1] + q0[iq3+1] + b->z[iq3+1];
	  q[iq3+2] = b->dt * q[iq3+2] + q0[iq3+2] + b->z[iq3+2];
	}
      else
	{
	  // loop for bonds in group "ig"
	  int i;
	  for (i = 0; i < g->np - 1; i ++)
	    {
	      int iq = ig0 + i;
	      int iq3 = iq * 3;

	      // calc imposed flow
	      fastSI_calc_imposed (b->BD->sys,
				   q0 + iq3,
				   q  + iq3);

	      q[iq3  ] += P[iq3  ];
	      q[iq3+1] += P[iq3+1];
	      q[iq3+2] += P[iq3+2];

	      // update the connector
	      /* **NOTE** b->z[] should be in the mobility mode
	       *          and stored in the connector index.
	       */
	      q[iq3  ] = b->dt * q[iq3  ] + q0[iq3  ] + b->z[iq3  ];
	      q[iq3+1] = b->dt * q[iq3+1] + q0[iq3+1] + b->z[iq3+1];
	      q[iq3+2] = b->dt * q[iq3+2] + q0[iq3+2] + b->z[iq3+2];
	    }
	  // the COM component
	  int iCOM  = ig0 + i;
	  int iCOM3 = iCOM * 3;

	  fastSI_calc_imposed (b->BD->sys,
			       q0 + iCOM3,
			       q  + iCOM3);

	  q[iCOM3  ] += P[iCOM3  ];
	  q[iCOM3+1] += P[iCOM3+1];
	  q[iCOM3+2] += P[iCOM3+2];

	  q[iCOM3  ] = b->dt * q[iCOM3  ] + q0[iCOM3  ] + b->z[iCOM3  ];
	  q[iCOM3+1] = b->dt * q[iCOM3+1] + q0[iCOM3+1] + b->z[iCOM3+1];
	  q[iCOM3+2] = b->dt * q[iCOM3+2] + q0[iCOM3+2] + b->z[iCOM3+2];
	}
      ig0 += g->np;
    }

  free (f);
  free (P);
}


void
fastSI_solve (struct BD_imp *b,
	      const double *q0,
	      double *q)
{
  int np = b->BD->sys->np;
  int np3 = np * 3;
  double *qC = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (qC, "fastSI_evolve");

  // set the predictor
  fastSI_solve_Euler (b, q0, qC);
  // here qC[] is the predictor

  // solve the initial corrector
  fastSI_solve_corrector (b, q0, qC,
			  1, // first step
			  q);
  // here q[] is the initial corrector

  double eps = 1.0; // just for fail safe
  int iter = 0;
  do
    {
      int i;
      for (i = 0; i < np3; i ++)
	{
	  qC[i] = q[i];
	}
      eps = fastSI_solve_corrector (b, q0, qC,
				    0, // update steps
				    q);
      //fprintf (stdout, "# iter = %d, eps = %e\n", iter, eps);
      iter ++;
    }
  while (eps > b->eps && iter < b->itmax);

  if (eps > b->eps)
    {
      fprintf (stdout, "fastSI_solve : fail to converge. iter = %d, eps = %e\n",
	       iter, eps);
      exit (1);
    }
  /*
  if (b->verbose != 0)
    {
      fprintf (stdout, "fastSI_solve : iter = %d, eps = %e\n",
	       iter, eps);
    }
  */
  fprintf (stdout, "fastSI_solve : iter = %d, eps = %e\n",
	   iter, eps);

  free (qC);
}


/**
 * for nonlinear solvers
 */
/* set the initial connector vector q0[np*3] in b->x0[]
 */
void
fastSI_set_Q0 (struct BD_imp *b,
	       const double *q0)
{
  int i;
  for (i = 0; i < b->BD->sys->np * 3; i ++)
    {
      b->x0[i] = q0[i];
    }
}


/* calc nonlinear equation f[] = 0, where
 * f[i] = q0[i]
 *       + dt * [uimp((q0[i]+q[i])/2)
 *               + (1/a[i+1]) * F(q[i+1])
 *               - (1/a[i+1] + 1/a[i]) * F(q[i])
 *               + (1/a[i]) * F(q[i-1])]
 *       + Y[i]
 *       - q[i]
 * INTPUT
 *  b : struct BD_imp
 *      NOTE that b->x0[np*3] should be set by q0[] before calling
 */
void
fastSI_f (struct BD_imp *b,
	  const double *q,
	  double *f)
{
  // loop for groups
  int ig0 = 0; // head pointer of the group "ig"
  int ng = b->BD->groups->n;
  int ig;
  for (ig = 0; ig < ng; ig ++)
    {
      struct BONDS_GROUP *g = b->BD->groups->group[ig];

      if (g->np == 1)
	{
	  // isolated particle
	  int iq = ig0;
	  int iq3 = iq * 3;
	  double qp[3];
	  fastSI_update_isolated
	    (b,
	     b->x0, q, // for imposed flow
	     q,     q, // for the interactions
	     g->ip[0], // particle index
	     iq, // index for connector (divided by 3)
	     qp);

	  f[iq3  ] = qp[0] - q[iq3  ];
	  f[iq3+1] = qp[1] - q[iq3+1];
	  f[iq3+2] = qp[2] - q[iq3+2];
	}
      else
	{
	  // loop for bonds in group "ig"
	  int i;
	  for (i = 0; i < g->np - 1; i ++)
	    {
	      int ib = g->bonds[i];
	      int iq = ig0 + i;
	      int iq3 = iq * 3;
	      fastSI_rhs (b,
			  b->x0, q, // for imposed flow
			  q,     q, // for the interactions
			  ib, // bond index
			  iq, // index for connector (divided by 3)
			  f + iq3);

	      // subtract the left-hand side term
	      double fs[3];
	      // for (i) spring with q[]
	      BONDS_calc_force_spring_i (b->BD->bonds,
					 ib, // bond index
					 q + iq3,
					 fs);
	      double B = b->dt * calc_B (b, ib);
	      f[iq3  ] -= q[iq3  ] + B * fs[0];
	      f[iq3+1] -= q[iq3+1] + B * fs[1];
	      f[iq3+2] -= q[iq3+2] + B * fs[2];
	    }
	  // the COM component
	  int iCOM  = ig0 + i;
	  int iCOM3 = iCOM * 3;
	  double qCOM[3];
	  fastSI_update_COM
	    (b,
	     b->x0, q, // for imposed flow
	     q,     q, // for the interactions
	     ig, // group index
	     iCOM, // index for connector (divided by 3)
	     qCOM);

	  f[iCOM3  ] = qCOM[0] - q[iCOM3  ];
	  f[iCOM3+1] = qCOM[1] - q[iCOM3+1];
	  f[iCOM3+2] = qCOM[2] - q[iCOM3+2];
	}
      ig0 += g->np;
    }
}


/**
 * NITSOL stuff
 */
/* a wrapper of the nonlinear equations for the semi-implicit algorithms
 * (both JGdP and siPC) for NITSOL
 * INTPUT
 *  xcur[n] : = (x[nm*3])          for F version
 *            = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  rpar    : (struct BD_imp *)BDimp
 *  ipar    : not used.
 * OUTPUT
 *  itrmf   : 0 means success
 *  fcur[n] := x - x0
 *           - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for JGdP, or
 *  fcur[n] := x - x0
 *           - (dt/2) * (U^pr
 *                       + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for siPC, where U^pr is the predictor obtained by Euler scheme as 
 *  U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
void
fastSI_NITSOL_f (int *n, double *xcur, double *fcur,
		 double *rpar, int *ipar, int *itrmf)
{
  // let's say that (double *)rpar corresponds to (struct BD_imp *)BDimp
  struct BD_imp *b = (struct BD_imp *)rpar;

  fastSI_f (b, xcur, fcur);

  *itrmf = 0;
}
void
fastSI_NITSOL_wrap (struct BD_imp *b,
		    const double *q0,
		    double *q)
{
  // set the initial connector in b->x0[]
  fastSI_set_Q0 (b, q0);

  int np = b->BD->sys->np;
  int np3 = np * 3;
  // copy q0[] -> q[]
  int i;
  for (i = 0; i < np3; i ++)
    {
      q[i] = q0[i];
    }
  NITSOL_solve (b->nit, q);
}



/**
 * GSL related stuff
 */
/* wrapper of fastSI_f() for GSL-MULTROOT routine
 */
int
fastSI_GSL_MULTIROOT_func (const gsl_vector *x, void *p,
			   gsl_vector *f)
{
  struct BD_imp *b = (struct BD_imp *)p;

  int np3 = 3 * b->BD->sys->np;
  double *xcur = (double *)malloc (sizeof (double) * np3);
  double *fcur = (double *)malloc (sizeof (double) * np3);
  CHECK_MALLOC (xcur, "fastSI_GSL_MULTIROOT_func");
  CHECK_MALLOC (fcur, "fastSI_GSL_MULTIROOT_func");

  int i;
  for (i = 0; i < np3; i ++)
    {
      xcur[i] = gsl_vector_get (x, i);
    }

  fastSI_f (b, xcur, fcur);

  for (i = 0; i < np3; i ++)
    {
      gsl_vector_set (f, i, fcur[i]);
    }

  free (xcur);
  free (fcur);

  return (GSL_SUCCESS);
}

void
fastSI_GSL_MULTIROOT_wrap (struct BD_imp *b,
			   const double *q0,
			   double *q)
{
  int np3 = 3 * b->BD->sys->np;

  // set the initial connector in b->x0[]
  fastSI_set_Q0 (b, q0);

  /**
   * set the initial guess
   */
  int i;
  for (i = 0; i < np3; i ++)
    {
      gsl_vector_set (b->guess, i, q0[i]);
    }
  gsl_multiroot_fsolver_set (b->S, b->F, b->guess);

  int status;
  int iter = 0;
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (b->S);

      if (status)   /* check if solver is stuck */
	break;

      status = gsl_multiroot_test_residual (b->S->f, b->eps);
    }
  while (status == GSL_CONTINUE && iter < b->itmax);

  if (status != GSL_SUCCESS)
    {
      fprintf (stdout, "status = %s\n", gsl_strerror (status));
    }

  /**
   * retreive the solution
   */
  gsl_vector *root = gsl_multiroot_fsolver_root (b->S);
  for (i = 0; i < np3; i ++)
    {
      q[i] = gsl_vector_get (root, i);
    }
}



/**
 * general routines
 */
void
fastSI_calc_FB (struct BD_imp *b,
		double *X)
{
  if (b->BD->flag_noHI == 1)
    {
      int i;
      for (i = 0; i < b->BD->sys->np; i ++)
	{
	  double factor = sqrt (2.0 / b->BD->peclet
				* b->dt
				/ b->BD->sys->a[i]);
	  fprintf (stderr, "# FB : factor = %e\n", factor);

	  int ix = i * 3;

	  X[ix  ] = KIrand_Gaussian (b->BD->rng) * factor;
	  X[ix+1] = KIrand_Gaussian (b->BD->rng) * factor;
	  X[ix+2] = KIrand_Gaussian (b->BD->rng) * factor;
	}
    }
  else
    {
      fprintf (stdout, "# fastSI_calc_FB() : HI is not implemented.\n");
      exit (1);
    }
}

/* evolve position of particles by semi-implicit predictor-corrector
 * INPUT
 *  t       : current time
 *  b       : struct BD_imp
 *            b->solver == 0 : GSL-multiroot
 *                           == 1 : NITSOL
 *                           == 2 : fastSI
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
fastSI_evolve (double t,
	       struct BD_imp *b,
	       double *x,
	       double dt)
{
  /**
   * update BD_imp data
   */
  //b->dt = dt;
  BD_imp_set_dt (b, dt);

  // calc brownian force for the configuration x[]
  int np = b->BD->sys->np;
  int np3 = np * 3;
  double *X = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (X, "fastSI_evolve");
  fastSI_calc_FB (b, X);

  // set the random vector for the springs
  /* **NOTE** b->z[] is used for it
   *          this should be in the mobility mode
   *          and stored in the connector index.
   */
  BONDS_pos_to_conn (b->BD->bonds,
		     b->BD->groups,
		     X,
		     b->z);
  free (X);

  double *q0 = (double *)calloc (np3, sizeof (double));
  double *q  = (double *)calloc (np3, sizeof (double));
  CHECK_MALLOC (q0, "fastSI_evolve");
  CHECK_MALLOC (q,  "fastSI_evolve");

  // set connector vectors
  BONDS_pos_to_conn (b->BD->bonds,
		     b->BD->groups,
		     x,
		     q0);

  // set the initial connector in b->x0[]
  fastSI_set_Q0 (b, q0);

  /**
   * solve the nonlinear equations
   */
  if (b->solver == 2) // fastSI
    {
      fastSI_solve (b, q0, q);
    }
  else if (b->solver == 1) // NITSOL
    {
      fastSI_NITSOL_wrap (b, q0, q);
    }
  else if (b->solver == 0) // GSL-MULTIROOT
    {
      fastSI_GSL_MULTIROOT_wrap (b, q0, q);
    }
  else
    {
      fprintf (stderr, "# fastSI_evolve() : invalid solver %d\n",
	       b->solver);
    }

  free (q0);
  free (q);

  return (dt);
}
