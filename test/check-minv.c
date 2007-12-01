/* test code for the calculation of minv
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-minv.c,v 1.2 2007/12/01 18:33:09 kichiki Exp $
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
#include "memory-check.h"

#include <non-ewald.h>
#include <dgetri_c.h> // lapack_inv_()
#include <fts.h> // scalar_minv_fts()
#include <minv-poly.h>

#include "check.h" // compare()


/*
 * note that tilde part is not checked.
 * to ensure that, check the symmetry for the whole matrix.
 * INPUT
 *  flag_sym : 0 for the matrix in extracted form
 *             1 for the matrix in symmetric form, that is,
 *             the matrix multiplied by p for M part, where
 *             p = (2,0,0,0,1,
 *                  0,2,0,0,0,
 *                  0,0,2,0,0,
 *                  0,0,0,2,0,
 *                  1,0,0,0,2).
 */
static int
check_matrix_fts_with_ez_alphabeta (int n11, double *mat,
				    int flag_sym,
				    double xa, double ya,
				    double yb,
				    double xc, double yc,
				    double xg, double yg,
				    double yh,
				    double xm, double ym, double zm,
				    int verbose, double tiny, double *max)
{
  int check = 0;

  // A_{alpha,beta}
  check += compare_max (mat[ 0*n11+ 0], ya,  "A (Ux,Fx) ya", verbose, tiny, max);
  check += compare_max (mat[ 0*n11+ 1], 0.0, "A (Ux,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 0*n11+ 2], 0.0, "A (Ux,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 1*n11+ 0], 0.0, "A (Uy,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 1*n11+ 1], ya,  "A (Uy,Fy) ya", verbose, tiny, max);
  check += compare_max (mat[ 1*n11+ 2], 0.0, "A (Uy,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 2*n11+ 0], 0.0, "A (Uz,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 2*n11+ 1], 0.0, "A (Uz,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 2*n11+ 2], xa,  "A (Uz,Fz) xa", verbose, tiny, max);

  // B_{alpha,beta}
  check += compare_max (mat[ 3*n11+ 0], 0.0, "B (Ox,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 3*n11+ 1], yb,  "B (Ox,Fy) yb", verbose, tiny, max);
  check += compare_max (mat[ 3*n11+ 2], 0.0, "B (Ox,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 4*n11+ 0], -yb, "B (Oy,Fx) -yb", verbose, tiny, max);
  check += compare_max (mat[ 4*n11+ 1], 0.0, "B (Oy,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 4*n11+ 2], 0.0, "B (Oy,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 5*n11+ 0], 0.0, "B (Oz,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 5*n11+ 1], 0.0, "B (Oz,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 5*n11+ 2], 0.0, "B (Oz,Fz) 0", verbose, tiny, max);

  // C_{alpha,beta}
  check += compare_max (mat[ 3*n11+ 3], yc,  "C (Ox,Ox) yc", verbose, tiny, max);
  check += compare_max (mat[ 3*n11+ 4], 0.0, "C (Ox,Oy) 0", verbose, tiny, max);
  check += compare_max (mat[ 3*n11+ 5], 0.0, "C (Ox,Oz) 0", verbose, tiny, max);

  check += compare_max (mat[ 4*n11+ 3], 0.0, "C (Oy,Ox) 0", verbose, tiny, max);
  check += compare_max (mat[ 4*n11+ 4], yc,  "C (Oy,Oy) yc", verbose, tiny, max);
  check += compare_max (mat[ 4*n11+ 5], 0.0, "C (Oy,Oz) 0", verbose, tiny, max);

  check += compare_max (mat[ 5*n11+ 3], 0.0, "C (Oz,Ox) 0", verbose, tiny, max);
  check += compare_max (mat[ 5*n11+ 4], 0.0, "C (Oz,Oy) 0", verbose, tiny, max);
  check += compare_max (mat[ 5*n11+ 5], xc,  "C (Oz,Oz) xc", verbose, tiny, max);

  // G_{alpha,beta}
  check += compare_max (mat[ 6*n11+ 0], 0.0, "G (Exx,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 6*n11+ 1], 0.0, "G (Exx,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 6*n11+ 2], -xg/3.0, "G (Exx,Fz) -xg/3", verbose, tiny, max);

  check += compare_max (mat[ 7*n11+ 0], 0.0, "G (Exy,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 7*n11+ 1], 0.0, "G (Exy,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 7*n11+ 2], 0.0, "G (Exy,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 8*n11+ 0], yg,  "G (Exz,Fx) yg", verbose, tiny, max);
  check += compare_max (mat[ 8*n11+ 1], 0.0, "G (Exz,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[ 8*n11+ 2], 0.0, "G (Exz,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[ 9*n11+ 0], 0.0, "G (Eyz,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[ 9*n11+ 1], yg,  "G (Eyz,Fy) yg", verbose, tiny, max);
  check += compare_max (mat[ 9*n11+ 2], 0.0, "G (Eyz,Fz) 0", verbose, tiny, max);

  check += compare_max (mat[10*n11+ 0], 0.0, "G (Eyy,Fx) 0", verbose, tiny, max);
  check += compare_max (mat[10*n11+ 1], 0.0, "G (Eyy,Fy) 0", verbose, tiny, max);
  check += compare_max (mat[10*n11+ 2], -xg/3.0, "G (Eyy,Fz) -xg/3", verbose, tiny, max);
  // H_{alpha,beta}
  check += compare_max (mat[ 6*n11+ 3], 0.0, "H (Exx,Tx) 0", verbose, tiny, max);
  check += compare_max (mat[ 6*n11+ 4], 0.0, "H (Exx,Ty) 0", verbose, tiny, max);
  check += compare_max (mat[ 6*n11+ 5], 0.0, "H (Exx,Tz) 0", verbose, tiny, max);

  check += compare_max (mat[ 7*n11+ 3], 0.0, "H (Exy,Tx) 0", verbose, tiny, max);
  check += compare_max (mat[ 7*n11+ 4], 0.0, "H (Exy,Ty) 0", verbose, tiny, max);
  check += compare_max (mat[ 7*n11+ 5], 0.0, "H (Exy,Tz) 0", verbose, tiny, max);

  check += compare_max (mat[ 8*n11+ 3], 0.0, "H (Exz,Tx) 0", verbose, tiny, max);
  check += compare_max (mat[ 8*n11+ 4], yh,  "H (Exz,Ty) yh", verbose, tiny, max);
  check += compare_max (mat[ 8*n11+ 5], 0.0, "H (Exz,Tz) 0", verbose, tiny, max);

  check += compare_max (mat[ 9*n11+ 3], -yh, "H (Eyz,Tx) -yh", verbose, tiny, max);
  check += compare_max (mat[ 9*n11+ 4], 0.0, "H (Eyz,Ty) 0", verbose, tiny, max);
  check += compare_max (mat[ 9*n11+ 5], 0.0, "H (Eyz,Tz) 0", verbose, tiny, max);

  check += compare_max (mat[10*n11+ 3], 0.0, "H12 (Eyy,Tx) 0", verbose, tiny, max);
  check += compare_max (mat[10*n11+ 4], 0.0, "H12 (Eyy,Ty) 0", verbose, tiny, max);
  check += compare_max (mat[10*n11+ 5], 0.0, "H12 (Eyy,Tz) 0", verbose, tiny, max);

  // M_{alpha,beta}
  if (flag_sym == 0)
    {
      // extracted form
      check += compare_max (mat[ 6*n11+ 6], xm/6.0+zm/2.0, "M (Exx,Sxx) xm/6+zm/2", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 7], 0.0, "M (Exx,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 8], 0.0, "M (Exx,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 9], 0.0, "M (Exx,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+10], xm/6.0-zm/2.0, "M (Exx,Syy) xm/6-zm/2", verbose, tiny, max);

      check += compare_max (mat[ 7*n11+ 6], 0.0, "M (Exy,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 7], zm/2.0, "M (Exy,Sxy) zm/2", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 8], 0.0, "M (Exy,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 9], 0.0, "M (Exy,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+10], 0.0, "M (Exy,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[ 8*n11+ 6], 0.0, "M (Exz,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 7], 0.0, "M (Exz,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 8], ym/2.0, "M (Exz,Sxz) ym/2", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 9], 0.0, "M (Exz,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+10], 0.0, "M (Exz,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[ 9*n11+ 6], 0.0, "M (Eyz,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 7], 0.0, "M (Eyz,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 8], 0.0, "M (Eyz,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 9], ym/2.0, "M (Eyz,Syz) ym/2", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+10], 0.0, "M (Eyz,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[10*n11+ 6], xm/6.0-zm/2.0, "M (Eyy,Sxx) xm/6-zm/2", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 7], 0.0, "M (Eyy,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 8], 0.0, "M (Eyy,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 9], 0.0, "M (Eyy,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+10], xm/6.0+zm/2.0, "M (Eyy,Syy) xm/6+zm/2", verbose, tiny, max);
    }
  else
    {
      // symmetric form
      check += compare_max (mat[ 6*n11+ 6], xm/2.0+zm/2.0, "M (Exx,Sxx) xm/2+zm/2", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 7], 0.0, "M (Exx,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 8], 0.0, "M (Exx,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+ 9], 0.0, "M (Exx,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 6*n11+10], xm/2.0-zm/2.0, "M (Exx,Syy) xm/2-zm/2", verbose, tiny, max);

      check += compare_max (mat[ 7*n11+ 6], 0.0, "M (Exy,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 7], zm,  "M (Exy,Sxy) zm", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 8], 0.0, "M (Exy,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+ 9], 0.0, "M (Exy,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 7*n11+10], 0.0, "M (Exy,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[ 8*n11+ 6], 0.0, "M (Exz,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 7], 0.0, "M (Exz,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 8], ym,  "M (Exz,Sxz) ym", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+ 9], 0.0, "M (Exz,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[ 8*n11+10], 0.0, "M (Exz,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[ 9*n11+ 6], 0.0, "M (Eyz,Sxx) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 7], 0.0, "M (Eyz,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 8], 0.0, "M (Eyz,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+ 9], ym,  "M (Eyz,Syz) ym", verbose, tiny, max);
      check += compare_max (mat[ 9*n11+10], 0.0, "M (Eyz,Syy) 0", verbose, tiny, max);

      check += compare_max (mat[10*n11+ 6], xm/2.0-zm/2.0, "M (Eyy,Sxx) xm/2-zm/2", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 7], 0.0, "M (Eyy,Sxy) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 8], 0.0, "M (Eyy,Sxz) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+ 9], 0.0, "M (Eyy,Syz) 0", verbose, tiny, max);
      check += compare_max (mat[10*n11+10], xm/2.0+zm/2.0, "M (Eyy,Syy) xm/2+zm/2", verbose, tiny, max);
    }

  return (check);
}

/*
 * INPUT
 *  flag_sym : 0 for the matrix in extracted form
 *             1 for the matrix in symmetric form, that is,
 *             the matrix multiplied by p for M part, where
 *             p = (2,0,0,0,1,
 *                  0,2,0,0,0,
 *                  0,0,2,0,0,
 *                  0,0,0,2,0,
 *                  1,0,0,0,2).
 */
static int
get_scalars_from_matrix_fts_with_ez_alphabeta
(int n11, double *mat,
 int flag_sym,
 double *xa, double *ya,
 double *yb,
 double *xc, double *yc,
 double *xg, double *yg,
 double *yh,
 double *xm, double *ym, double *zm,
 int verbose, double tiny, double *max)
{
  int check = 0;


  // A_{alpha,beta}
  *ya = mat[ 0*n11+ 0];
  *xa = mat[ 2*n11+ 2];

  // B_{alpha,beta}
  *yb = mat[ 3*n11+ 1];

  // C_{alpha,beta}
  *yc = mat[ 3*n11+ 3];
  *xc = mat[ 5*n11+ 5];

  // G_{alpha,beta}
  *xg = -3.0*mat[ 6*n11+ 2];
  *yg = mat[ 8*n11+ 0];

  // H_{alpha,beta}
  *yh = mat[ 8*n11+ 4];

  // M_{alpha,beta}
  double mp;
  double mm;
  if (flag_sym == 0)
    {
      // extracted form
      mp = mat[ 6*n11+ 6]; // = xm/6.0+zm/2.0
      mm = mat[ 6*n11+10]; // = xm/6.0-zm/2.0
      *xm = 3.0 * (mp + mm);
      *zm = 2.0 * mat[ 7*n11+ 7];
      *ym = 2.0 * mat[ 8*n11+ 8];
    }
  else
    {
      // symmetric form
      mp = mat[ 6*n11+ 6]; // = xm/2.0+zm/2.0
      mm = mat[ 6*n11+10]; // = xm/2.0-zm/2.0
      *xm = mp + mm;
      *zm = mat[ 7*n11+ 7];
      *ym = mat[ 8*n11+ 8];
    }

  // check the whole matrix
  check = check_matrix_fts_with_ez_alphabeta (n11, mat,
					      flag_sym,
					      *xa, *ya,
					      *yb,
					      *xc, *yc,
					      *xg, *yg,
					      *yh,
					      *xm, *ym, *zm,
					      verbose, tiny, max);
  if (verbose != 0)
    {
      if (check > 0)
	{
	  fprintf (stdout,
		   "get_scalars_from_matrix_fts_with_ez_alphabeta: "
		   "inconsistent!!\n");
	}
    }

  return (check);
}


static int
check_matrix_2B_fts_with_ez (double *mat,
			     double *scalars,
			     int verbose, double tiny, double *max)
{
  int check = 0;

  int n11 = 2 * 11;

  double xa11 = scalars [0];
  double xa12 = scalars [1];
  double ya11 = scalars [2];
  double ya12 = scalars [3];
  double yb11 = scalars [4];
  double yb12 = scalars [5];
  double xc11 = scalars [6];
  double xc12 = scalars [7];
  double yc11 = scalars [8];
  double yc12 = scalars [9];
  double xg11 = scalars[10];
  double xg12 = scalars[11];
  double yg11 = scalars[12];
  double yg12 = scalars[13];
  double yh11 = scalars[14];
  double yh12 = scalars[15];
  double xm11 = scalars[16];
  double xm12 = scalars[17];
  double ym11 = scalars[18];
  double ym12 = scalars[19];
  double zm11 = scalars[20];
  double zm12 = scalars[21];

  double xa22 =  xa11;
  double xa21 =  xa12;
  double ya22 =  ya11;
  double ya21 =  ya12;
  double yb22 = -yb11;
  double yb21 = -yb12;
  double xc22 =  xc11;
  double xc21 =  xc12;
  double yc22 =  yc11;
  double yc21 =  yc12;
  double xg22 = -xg11;
  double xg21 = -xg12;
  double yg22 = -yg11;
  double yg21 = -yg12;
  double yh22 =  yh11;
  double yh21 =  yh12;
  double xm22 =  xm11;
  double xm21 =  xm12;
  double ym22 =  ym11;
  double ym21 =  ym12;
  double zm22 =  zm11;
  double zm21 =  zm12;

  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 0*n11+0, // 11
					       0, // extracted form
					       xa11, ya11,
					       yb11,
					       xc11, yc11,
					       xg11, yg11,
					       yh11,
					       xm11, ym11, zm11,
					       verbose, tiny, max);
  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+11, // 22
					       0, // extracted form
					       xa22, ya22,
					       yb22,
					       xc22, yc22,
					       xg22, yg22,
					       yh22,
					       xm22, ym22, zm22,
					       verbose, tiny, max);

  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 0*n11+11, // 12
					       0, // extracted form
					       xa12, ya12,
					       yb12,
					       xc12, yc12,
					       xg12, yg12,
					       yh12,
					       xm12, ym12, zm12,
					       verbose, tiny, max);
  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+0, // 21
					       0, // extracted form
					       xa21, ya21,
					       yb21,
					       xc21, yc21,
					       xg21, yg21,
					       yh21,
					       xm21, ym21, zm21,
					       verbose, tiny, max);

  // symmetry check for tilde parts
  int i, j;
  char label[80];
  for (i = 0; i < n11; i ++)
    {
      for (j = i+1; j < n11; j ++)
	{
	  sprintf (label, "Symmetry Check (%d,%d)", i, j);
	  check += compare_max (mat[i*n11+j], mat[j*n11+i],
				label, verbose, tiny, max);
	}
    }

  return (check);
}


/* check matrix_mob_nonewald_3all() in FTS version with e=(0,0,1)
 */
int
check_matrix_mob_nonewald_fts (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_matrix_mob_nonewald_fts(r=%f) : start\n", r);
    }

  int check = 0;
  double max = 0.0;


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_matrix_mob_nonewald_fts");
  sys->version = 2; // FTS

  int np = 2;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;
  sys->pos [3] = 0.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = r;
  // so that r12 = r2 - r1 = (0, 0, r)

  int n11 = np * 11;
  double *mat = NULL;
  mat = (double *)malloc (sizeof (double) * n11 * n11);
  CHECK_MALLOC (mat, "check_matrix_mob_nonewald_fts");
  make_matrix_mob_nonewald_3all (sys, mat);


  double mob [11];
  scalars_nonewald (sys->version, r, mob);
  /*  scalar [11]:
   *   0, 1,    : (xa12, ya12) for F version
   *   2,       : (yb12)
   *   3, 4,    : (xc12, yc12) for FT version
   *   5, 6,    : (xg12, yg12)
   *   7,       : (yh12)
   *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
   */
  double xa12 = mob [0];
  double ya12 = mob [1];
  double yb12 = mob [2];
  double xc12 = mob [3];
  double yc12 = mob [4];
  double xg12 = mob [5];
  double yg12 = mob [6];
  double yh12 = mob [7];
  double xm12 = mob [8];
  double ym12 = mob [9];
  double zm12 = mob[10];

  double xa11 = 1.0;
  double ya11 = 1.0;
  double yb11 = 0.0;
  double xc11 = 0.75;
  double yc11 = 0.75;
  double xg11 = 0.0;
  double yg11 = 0.0;
  double yh11 = 0.0;
  double xm11 = 0.9;
  double ym11 = 0.9;
  double zm11 = 0.9;

  double scalars[22];
  scalars [0] = xa11;
  scalars [1] = xa12;
  scalars [2] = ya11;
  scalars [3] = ya12;
  scalars [4] = yb11;
  scalars [5] = yb12;
  scalars [6] = xc11;
  scalars [7] = xc12;
  scalars [8] = yc11;
  scalars [9] = yc12;
  scalars[10] = xg11;
  scalars[11] = xg12;
  scalars[12] = yg11;
  scalars[13] = yg12;
  scalars[14] = yh11;
  scalars[15] = yh12;
  scalars[16] = xm11;
  scalars[17] = xm12;
  scalars[18] = ym11;
  scalars[19] = ym12;
  scalars[20] = zm11;
  scalars[21] = zm12;

  check += check_matrix_2B_fts_with_ez (mat, scalars, verbose, tiny, &max);

  stokes_free (sys);
  free (mat);


  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

void
conv_mat_from_ext_to_e2e_alphabeta (int n, double *mat)
{
  double p[25] = {2.0, 0.0, 0.0, 0.0, 1.0,
		  0.0, 2.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 2.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 2.0, 0.0,
		  1.0, 0.0, 0.0, 0.0, 2.0};
  double tmp[55]; // 11 * 5

  int i, j, k;
  for (i = 0; i < 11; i ++)
    {
      for (j = 0; j < 5; j ++)
	{
	  tmp[i*5+j] = 0.0;
	  for (k = 0; k < 5; k ++)
	    {
	      tmp[i*5+j] += mat[i*n+(k+6)] * p[k*5+j];
	    }
	}
    }

  for (i = 0; i < 11; i ++)
    {
      for (j = 0; j < 5; j ++)
	{
	  mat[i*n+(j+6)] = tmp[i*5+j];
	}
    }
}

/* test the following functions
 *  make_matrix_mob_nonewald_3all() in non-ewald.c,
 *  scalar_minv_fts() in fts.c,
 *  scalars_minv_fts_poly() in minv-poly.c.
 * 1) calculate scalar functions of (M^infty)^-1 in FTS version directly by
 *    make_matrix_mob_nonewald_3all() inversing by lapack_inv_().
 *    check the symmetry on the inverse matrix.
 * 2) calculate scalar functions by the function
 *    scalar_minv_fts() in fts.c,
 *    and compare with the results in 1).
 * 3) calculate scalar functions by the function
 *    scalars_minv_fts_poly() in minv-poly.c
 *    with a1=a2=a (monodisperse case)
 *    and compare with the results in 1).
 */
int
check_minv_fts (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_minv_fts(r=%f) : start\n", r);
    }

  int check = 0;
  double max = 0.0;


  /* calculate the inverse of M^infty directly:
   * M^infty in FTS version is calculated by make_matrix_mob_nonewald_3all()
   * then, inserved by lapack_inv_()
   * then, the scalar functions are extracted based on the fact that e=(0,0,1).
   * the results, xa11, xa12, ..., are the reference to check.
   */
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_minv_fts");
  sys->version = 2; // FTS

  int np = 2;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;
  sys->pos [3] = 0.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = r;
  // so that r12 = r2 - r1 = (0, 0, r)

  int n11 = np * 11;
  double *mat = NULL;
  mat = (double *)malloc (sizeof (double) * n11 * n11);
  CHECK_MALLOC (mat, "check_matrix_mob_nonewald_fts");
  make_matrix_mob_nonewald_3all (sys, mat);
  // mat is M^inf in ext form

  conv_mat_from_ext_to_e2e_alphabeta (n11, mat +  0*n11+ 0); // 11
  conv_mat_from_ext_to_e2e_alphabeta (n11, mat +  0*n11+11); // 12
  conv_mat_from_ext_to_e2e_alphabeta (n11, mat + 11*n11+ 0); // 21
  conv_mat_from_ext_to_e2e_alphabeta (n11, mat + 11*n11+11); // 22

  lapack_inv_ (n11, mat);
  // mat is (M^inf)^-1 in sym form

  double xa11, xa12;
  double ya11, ya12;
  double yb11, yb12;
  double xc11, xc12;
  double yc11, yc12;
  double xg11, xg12;
  double yg11, yg12;
  double yh11, yh12;
  double xm11, xm12;
  double ym11, ym12;
  double zm11, zm12;

  double xa22, xa21;
  double ya22, ya21;
  double yb22, yb21;
  double xc22, xc21;
  double yc22, yc21;
  double xg22, xg21;
  double yg22, yg21;
  double yh22, yh21;
  double xm22, xm21;
  double ym22, ym21;
  double zm22, zm21;

  check +=
    get_scalars_from_matrix_fts_with_ez_alphabeta (n11, mat +  0*n11+ 0, // 11
						   1, // symmetric form
						   &xa11, &ya11,
						   &yb11,
						   &xc11, &yc11,
						   &xg11, &yg11,
						   &yh11,
						   &xm11, &ym11, &zm11,
						   verbose, tiny, &max);
  check +=
    get_scalars_from_matrix_fts_with_ez_alphabeta (n11, mat +  0*n11+11, // 12
						   1, // symmetric form
						   &xa12, &ya12,
						   &yb12,
						   &xc12, &yc12,
						   &xg12, &yg12,
						   &yh12,
						   &xm12, &ym12, &zm12,
						   verbose, tiny, &max);
  check +=
    get_scalars_from_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+ 0, // 21
						   1, // symmetric form
						   &xa21, &ya21,
						   &yb21,
						   &xc21, &yc21,
						   &xg21, &yg21,
						   &yh21,
						   &xm21, &ym21, &zm21,
						   verbose, tiny, &max);
  check +=
    get_scalars_from_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+11, // 22
						   1, // symmetric form
						   &xa22, &ya22,
						   &yb22,
						   &xc22, &yc22,
						   &xg22, &yg22,
						   &yh22,
						   &xm22, &ym22, &zm22,
						   verbose, tiny, &max);
  /* NOTE: the extraction routine get_scalar_from_matrix...() is for
   * e = (0,0,1), not e = (0,0,-1), that is, for (21)-interaction here.
   * that is why B and G parts for (21) are opposite sign to (12).
   */
  yb21 *= -1.0;
  yb22 *= -1.0;
  xg21 *= -1.0;
  xg22 *= -1.0;
  yg21 *= -1.0;
  yg22 *= -1.0;


  // check the symmetry for the monodisperse system
  check += compare_max (xa22, xa11, "XA22", verbose, tiny, &max);
  check += compare_max (ya22, ya11, "YA22", verbose, tiny, &max);
  check += compare_max (yb22, yb11, "YB22", verbose, tiny, &max);
  check += compare_max (xc22, xc11, "XC22", verbose, tiny, &max);
  check += compare_max (yc22, yc11, "YC22", verbose, tiny, &max);
  check += compare_max (xg22, xg11, "XG22", verbose, tiny, &max);
  check += compare_max (yg22, yg11, "YG22", verbose, tiny, &max);
  check += compare_max (yh22, yh11, "YH22", verbose, tiny, &max);
  check += compare_max (xm22, xm11, "XM22", verbose, tiny, &max);
  check += compare_max (ym22, ym11, "YM22", verbose, tiny, &max);
  check += compare_max (zm22, zm11, "ZM22", verbose, tiny, &max);

  check += compare_max (xa21, xa12, "XA21", verbose, tiny, &max);
  check += compare_max (ya21, ya12, "YA21", verbose, tiny, &max);
  check += compare_max (yb21, yb12, "YB21", verbose, tiny, &max);
  check += compare_max (xc21, xc12, "XC21", verbose, tiny, &max);
  check += compare_max (yc21, yc12, "YC21", verbose, tiny, &max);
  check += compare_max (xg21, xg12, "XG21", verbose, tiny, &max);
  check += compare_max (yg21, yg12, "YG21", verbose, tiny, &max);
  check += compare_max (yh21, yh12, "YH21", verbose, tiny, &max);
  check += compare_max (xm21, xm12, "XM21", verbose, tiny, &max);
  check += compare_max (ym21, ym12, "YM21", verbose, tiny, &max);
  check += compare_max (zm21, zm12, "ZM21", verbose, tiny, &max);


  /* Again, here, check_matrix...() is for e = (0,0,1), not e = (0,0,-1),
   * that is, for (21)-interaction here.
   * that is why B and G parts for (21) are opposite sign to (12).
   */
  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 0*n11+0, // 11
					       1, // symmetric form
					       xa11, ya11,
					       yb11,
					       xc11, yc11,
					       xg11, yg11,
					       yh11,
					       xm11, ym11, zm11,
					       verbose, tiny, &max);
  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+11, // 22
					       1, // symmetric form
					       xa22, ya22,
					       -yb22,
					       xc22, yc22,
					       -xg22, -yg22,
					       yh22,
					       xm22, ym22, zm22,
					       verbose, tiny, &max);

  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 0*n11+11, // 12
					       1, // symmetric form
					       xa12, ya12,
					       yb12,
					       xc12, yc12,
					       xg12, yg12,
					       yh12,
					       xm12, ym12, zm12,
					       verbose, tiny, &max);
  check += check_matrix_fts_with_ez_alphabeta (n11, mat + 11*n11+0, // 21
					       1, // symmetric form
					       xa21, ya21,
					       -yb21,
					       xc21, yc21,
					       -xg21, -yg21,
					       yh21,
					       xm21, ym21, zm21,
					       verbose, tiny, &max);

  // compare it with scalar_minv_fts() in fts.c
  double mono[22];
  scalar_minv_fts (r, mono);

  check += compare_max (mono [0], xa11, "minv_fts XA11", verbose, tiny, &max);
  check += compare_max (mono [1], xa12, "minv_fts XA12", verbose, tiny, &max);
  check += compare_max (mono [2], ya11, "minv_fts YA11", verbose, tiny, &max);
  check += compare_max (mono [3], ya12, "minv_fts YA12", verbose, tiny, &max);
  check += compare_max (mono [4], yb11, "minv_fts YB11", verbose, tiny, &max);
  check += compare_max (mono [5], yb12, "minv_fts YB12", verbose, tiny, &max);
  check += compare_max (mono [6], xc11, "minv_fts XC11", verbose, tiny, &max);
  check += compare_max (mono [7], xc12, "minv_fts XC12", verbose, tiny, &max);
  check += compare_max (mono [8], yc11, "minv_fts YC11", verbose, tiny, &max);
  check += compare_max (mono [9], yc12, "minv_fts YC12", verbose, tiny, &max);
  check += compare_max (mono[10], xg11, "minv_fts XG11", verbose, tiny, &max);
  check += compare_max (mono[11], xg12, "minv_fts XG12", verbose, tiny, &max);
  check += compare_max (mono[12], yg11, "minv_fts YG11", verbose, tiny, &max);
  check += compare_max (mono[13], yg12, "minv_fts YG12", verbose, tiny, &max);
  check += compare_max (mono[14], yh11, "minv_fts YH11", verbose, tiny, &max);
  check += compare_max (mono[15], yh12, "minv_fts YH12", verbose, tiny, &max);
  check += compare_max (mono[16], xm11, "minv_fts XM11", verbose, tiny, &max);
  check += compare_max (mono[17], xm12, "minv_fts XM12", verbose, tiny, &max);
  check += compare_max (mono[18], ym11, "minv_fts YM11", verbose, tiny, &max);
  check += compare_max (mono[19], ym12, "minv_fts YM12", verbose, tiny, &max);
  check += compare_max (mono[20], zm11, "minv_fts ZM11", verbose, tiny, &max);
  check += compare_max (mono[21], zm12, "minv_fts ZM12", verbose, tiny, &max);


  double poly[44];
  scalars_minv_fts_poly (r, 1.0, 1.0, poly);
  check += compare_max (poly [0], xa11, "minv_fts_poly XA11", verbose, tiny, &max);
  check += compare_max (poly [1], xa12, "minv_fts_poly XA12", verbose, tiny, &max);
  check += compare_max (poly [2], xa21, "minv_fts_poly XA21", verbose, tiny, &max);
  check += compare_max (poly [3], xa22, "minv_fts_poly XA22", verbose, tiny, &max);
  check += compare_max (poly [4], ya11, "minv_fts_poly YA11", verbose, tiny, &max);
  check += compare_max (poly [5], ya12, "minv_fts_poly YA12", verbose, tiny, &max);
  check += compare_max (poly [6], ya21, "minv_fts_poly YA21", verbose, tiny, &max);
  check += compare_max (poly [7], ya22, "minv_fts_poly YA22", verbose, tiny, &max);
  check += compare_max (poly [8], yb11, "minv_fts_poly YB11", verbose, tiny, &max);
  check += compare_max (poly [9], yb12, "minv_fts_poly YB12", verbose, tiny, &max);
  check += compare_max (poly[10], yb21, "minv_fts_poly YB21", verbose, tiny, &max);
  check += compare_max (poly[11], yb22, "minv_fts_poly YB22", verbose, tiny, &max);
  check += compare_max (poly[12], xc11, "minv_fts_poly XC11", verbose, tiny, &max);
  check += compare_max (poly[13], xc12, "minv_fts_poly XC12", verbose, tiny, &max);
  check += compare_max (poly[14], xc21, "minv_fts_poly XC21", verbose, tiny, &max);
  check += compare_max (poly[15], xc22, "minv_fts_poly XC22", verbose, tiny, &max);
  check += compare_max (poly[16], yc11, "minv_fts_poly YC11", verbose, tiny, &max);
  check += compare_max (poly[17], yc12, "minv_fts_poly YC12", verbose, tiny, &max);
  check += compare_max (poly[18], yc21, "minv_fts_poly YC21", verbose, tiny, &max);
  check += compare_max (poly[19], yc22, "minv_fts_poly YC22", verbose, tiny, &max);
  check += compare_max (poly[20], xg11, "minv_fts_poly XG11", verbose, tiny, &max);
  check += compare_max (poly[21], xg12, "minv_fts_poly XG12", verbose, tiny, &max);
  check += compare_max (poly[22], xg21, "minv_fts_poly XG21", verbose, tiny, &max);
  check += compare_max (poly[23], xg22, "minv_fts_poly XG22", verbose, tiny, &max);
  check += compare_max (poly[24], yg11, "minv_fts_poly YG11", verbose, tiny, &max);
  check += compare_max (poly[25], yg12, "minv_fts_poly YG12", verbose, tiny, &max);
  check += compare_max (poly[26], yg21, "minv_fts_poly YG21", verbose, tiny, &max);
  check += compare_max (poly[27], yg22, "minv_fts_poly YG22", verbose, tiny, &max);
  check += compare_max (poly[28], yh11, "minv_fts_poly YH11", verbose, tiny, &max);
  check += compare_max (poly[29], yh12, "minv_fts_poly YH12", verbose, tiny, &max);
  check += compare_max (poly[30], yh21, "minv_fts_poly YH21", verbose, tiny, &max);
  check += compare_max (poly[31], yh22, "minv_fts_poly YH22", verbose, tiny, &max);
  check += compare_max (poly[32], xm11, "minv_fts_poly XM11", verbose, tiny, &max);
  check += compare_max (poly[33], xm12, "minv_fts_poly XM12", verbose, tiny, &max);
  check += compare_max (poly[34], xm21, "minv_fts_poly XM21", verbose, tiny, &max);
  check += compare_max (poly[35], xm22, "minv_fts_poly XM22", verbose, tiny, &max);
  check += compare_max (poly[36], ym11, "minv_fts_poly YM11", verbose, tiny, &max);
  check += compare_max (poly[37], ym12, "minv_fts_poly YM12", verbose, tiny, &max);
  check += compare_max (poly[38], ym21, "minv_fts_poly YM21", verbose, tiny, &max);
  check += compare_max (poly[39], ym22, "minv_fts_poly YM22", verbose, tiny, &max);
  check += compare_max (poly[40], zm11, "minv_fts_poly ZM11", verbose, tiny, &max);
  check += compare_max (poly[41], zm12, "minv_fts_poly ZM12", verbose, tiny, &max);
  check += compare_max (poly[42], zm21, "minv_fts_poly ZM21", verbose, tiny, &max);
  check += compare_max (poly[43], zm22, "minv_fts_poly ZM22", verbose, tiny, &max);

  check += compare_max (poly [0], mono [0], "poly - mono XA11", verbose, tiny, &max);
  check += compare_max (poly [1], mono [1], "poly - mono XA12", verbose, tiny, &max);
  check += compare_max (poly [4], mono [2], "poly - mono YA11", verbose, tiny, &max);
  check += compare_max (poly [5], mono [3], "poly - mono YA12", verbose, tiny, &max);
  check += compare_max (poly [8], mono [4], "poly - mono YB11", verbose, tiny, &max);
  check += compare_max (poly [9], mono [5], "poly - mono YB12", verbose, tiny, &max);
  check += compare_max (poly[12], mono [6], "poly - mono XC11", verbose, tiny, &max);
  check += compare_max (poly[13], mono [7], "poly - mono XC12", verbose, tiny, &max);
  check += compare_max (poly[16], mono [8], "poly - mono YC11", verbose, tiny, &max);
  check += compare_max (poly[17], mono [9], "poly - mono YC12", verbose, tiny, &max);
  check += compare_max (poly[20], mono[10], "poly - mono XG11", verbose, tiny, &max);
  check += compare_max (poly[21], mono[11], "poly - mono XG12", verbose, tiny, &max);
  check += compare_max (poly[24], mono[12], "poly - mono YG11", verbose, tiny, &max);
  check += compare_max (poly[25], mono[13], "poly - mono YG12", verbose, tiny, &max);
  check += compare_max (poly[28], mono[14], "poly - mono YH11", verbose, tiny, &max);
  check += compare_max (poly[29], mono[15], "poly - mono YH12", verbose, tiny, &max);
  check += compare_max (poly[32], mono[16], "poly - mono XM11", verbose, tiny, &max);
  check += compare_max (poly[33], mono[17], "poly - mono XM12", verbose, tiny, &max);
  check += compare_max (poly[36], mono[18], "poly - mono YM11", verbose, tiny, &max);
  check += compare_max (poly[37], mono[19], "poly - mono YM12", verbose, tiny, &max);
  check += compare_max (poly[40], mono[20], "poly - mono ZM11", verbose, tiny, &max);
  check += compare_max (poly[41], mono[21], "poly - mono ZM12", verbose, tiny, &max);

  stokes_free (sys);
  free (mat);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
