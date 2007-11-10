/* test code for solve_gen_linear() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-solve-gen-linear.c,v 1.4 2007/11/10 23:10:21 kichiki Exp $
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
#include <math.h> // fabs()
#include "check.h" // compare()
#include "memory-check.h" // macro CHECK_MALLOC

#include <dgetri_c.h> // lapack_inv_()
#include <matrix.h> // mul_matrices()


static void
split_matrix (int n, const double *a,
	      int n1, int n2,
	      double *a_ll, double *a_lh, double *a_hl, double *a_hh)
{
  int i, j;
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  a_ll [i*n1+j] = a [i*n+j];
	}
      for (j = n1; j < n; j ++)
	{
	  a_lh [i*n2+(j-n1)] = a [i*n+j];
	}
    }
  for (i = n1; i < n; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  a_hl [(i-n1)*n1+j] = a [i*n+j];
	}
      for (j = n1; j < n; j ++)
	{
	  a_hh [(i-n1)*n2+(j-n1)] = a [i*n+j];
	}
    }
}
static void
merge_matrix (int n1, int n2,
	      double c_ll, const double *a_ll,
	      double c_lh, const double *a_lh,
	      double c_hl, const double *a_hl,
	      double c_hh, const double *a_hh,
	      int n, double *a)
{
  int i, j;
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  a [i*n+j] = c_ll * a_ll [i*n1+j];
	}
      for (j = n1; j < n; j ++)
	{
	  a [i*n+j] = c_lh * a_lh [i*n2+(j-n1)];
	}
    }
  for (i = n1; i < n; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  a [i*n+j] = c_hl * a_hl [(i-n1)*n1+j];
	}
      for (j = n1; j < n; j ++)
	{
	  a [i*n+j] = c_hh * a_hh [(i-n1)*n2+(j-n1)];
	}
    }
}

/* check for some local routines
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_split_merge (int n1, int n2,
		   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_split_merge : start\n");
    }

  int check = 0;


  int n = n1 + n2;

  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_split_merge");
  CHECK_MALLOC (b, "check_split_merge");

  double *a_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *a_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *a_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *a_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (a_ll, "check_split_merge");
  CHECK_MALLOC (a_lh, "check_split_merge");
  CHECK_MALLOC (a_hl, "check_split_merge");
  CHECK_MALLOC (a_hh, "check_split_merge");

  double *b_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *b_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *b_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *b_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (b_ll, "check_split_merge");
  CHECK_MALLOC (b_lh, "check_split_merge");
  CHECK_MALLOC (b_hl, "check_split_merge");
  CHECK_MALLOC (b_hh, "check_split_merge");


  int i;
  srand48(0);
  for (i = 0; i < n * n; i ++)
    {
      a [i] = drand48();
    }

  // a => split
  split_matrix (n, a, n1, n2, a_ll, a_lh, a_hl, a_hh);
  // merge => b
  merge_matrix (n1, n2, 1.0, a_ll, 1.0, a_lh, 1.0, a_hl, 1.0, a_hh, n, b);
  // b => split
  split_matrix (n, b, n1, n2, b_ll, b_lh, b_hl, b_hh);


  // compare a and b
  char label[80];
  int j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "check_split_merge : A-B[%d,%d]", i,j);
	  check += compare (a [i*n+j], b [i*n+j],
			    label, verbose, tiny);
	}
    }

  // compare a_** and b_**
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "check_split_merge : All-Bll[%d,%d]", i,j);
	  check += compare (a_ll [i*n1+j], b_ll [i*n1+j],
			    label, verbose, tiny);
	}
    }
  for (i = 0; i < n1; i ++)
    {
      for (j = 0; j < n2; j ++)
	{
	  sprintf (label, "check_split_merge : Alh-Blh[%d,%d]", i,j);
	  check += compare (a_lh [i*n2+j], b_lh [i*n2+j],
			    label, verbose, tiny);
	}
    }
  for (i = 0; i < n2; i ++)
    {
      for (j = 0; j < n1; j ++)
	{
	  sprintf (label, "check_split_merge : Ahl-Bhl[%d,%d]", i,j);
	  check += compare (a_hl [i*n1+j], b_hl [i*n1+j],
			    label, verbose, tiny);
	}
    }
  for (i = 0; i < n2; i ++)
    {
      for (j = 0; j < n2; j ++)
	{
	  sprintf (label, "check_split_merge : Ahh-Bhh[%d,%d]", i,j);
	  check += compare (a_hh [i*n2+j], b_hh [i*n2+j],
			    label, verbose, tiny);
	}
    }

  free (a);
  free (b);
  free (a_ll);
  free (a_lh);
  free (a_hl);
  free (a_hh);
  free (b_ll);
  free (b_lh);
  free (b_hl);
  free (b_hh);

  if (verbose != 0)
    {
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}




static void
solve_gen_linear_ (int n1, int n2,
		   double * A, double * B, double * C, double * D,
		   double * E, double * F, double * G, double * H,
		   double * I, double * J, double * K, double * L)
{
  int n = n1 + n2;

  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  double *c = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "solve_gen_linear_");
  CHECK_MALLOC (b, "solve_gen_linear_");
  CHECK_MALLOC (c, "solve_gen_linear_");

  merge_matrix (n1, n2, 1.0, A, -1.0, F, 1.0, C, -1.0, H, n, a);

  merge_matrix (n1, n2, 1.0, E, -1.0, B, 1.0, G, -1.0, D, n, b);

  // a = a^-1
  lapack_inv_ (n, a);

  // c = a^-1 . b
  mul_matrices (a, n, n,
		b, n, n,
		c);

  split_matrix (n, c,
		n1, n2, I, J, K, L);

  free (a);
  free (b);
  free (c);
}


static void
inverse_by_sub (int n1, int n2,
		const double *A, const double *B, const double *C, double *D,
		double *W, double *X, double *Y, double *Z)
{
  double *b = (double *)malloc (sizeof (double) * n1 * n2);
  double *c = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (b, "inverse_by_sub");
  CHECK_MALLOC (c, "inverse_by_sub");
  
  /* D := D^-1 */
  lapack_inv_ (n2, D);

  /* c := (D^-1) . C */
  mul_matrices (D, n2, n2, C, n2, n1, c);
  /* W [n1, n1] := A - B.(D^-1).C */
  add_and_mul (A, n1, n1, B, n1, n2, c, n2, n1, 1.0, -1.0, W);
  /* W :=  (A - B.(D^-1).C)^-1 */
  lapack_inv_ (n1, W);

  /* b := B . (D^-1) */
  mul_matrices (B, n1, n2, D, n2, n2, b);
  /* X [n1, n2] := - W . B.(D^-1) */
  mul_matrices (W, n1, n1, b, n1, n2, X);
  int i;
  for (i = 0; i < n1*n2; i ++) X[i] *= -1.0;

  /* Y [n2, n1] := - (D^-1).C . W */
  mul_matrices (c, n2, n1, W, n1, n1, Y);
  for (i = 0; i < n2*n1; i ++) Y[i] *= -1.0;
  /*
  dgemm_wrap (n2, n1, n1,
	      -1.0, c, W,
	      0.0, Y);
  */

  /* Z [n2, n2] := (D^-1) - (D^-1).C . X */
  add_and_mul (D, n2, n2, c, n2, n1, X, n1, n2, 1.0, -1.0, Z);

  free (b);
  free (c);
}


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_inverse_by_sub (int n1, int n2,
		      int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_inverse_by_sub : start\n");
    }

  int check = 0;


  int n = n1 + n2;

  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_inverse_by_sub");
  CHECK_MALLOC (b, "check_inverse_by_sub");

  // make symmetric random matrix
  int i, j;
  for (i = 0; i < n; i ++)
    {
      for (j = i; j < n; j ++)
	{
	  a[i*n+j] = drand48 ();
	}
    }
  for (i = 0; i < n; i ++)
    {
      for (j = i+1; j < n; j ++)
	{
	  a[j*n+i] = a[i*n+j];
	}
    }

  // copy
  for (i = 0; i < n*n; i ++)
    {
      b[i] = a[i];
    }

  // split
  double *a_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *a_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *a_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *a_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (a_ll, "check_inverse_by_sub");
  CHECK_MALLOC (a_lh, "check_inverse_by_sub");
  CHECK_MALLOC (a_hl, "check_inverse_by_sub");
  CHECK_MALLOC (a_hh, "check_inverse_by_sub");

  double *c_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *c_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *c_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *c_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (c_ll, "check_inverse_by_sub");
  CHECK_MALLOC (c_lh, "check_inverse_by_sub");
  CHECK_MALLOC (c_hl, "check_inverse_by_sub");
  CHECK_MALLOC (c_hh, "check_inverse_by_sub");

  split_matrix (n, a, n1, n2, a_ll, a_lh, a_hl, a_hh);
  inverse_by_sub (n1, n2,
		  a_ll, a_lh, a_hl, a_hh,
		  c_ll, c_lh, c_hl, c_hh);


  double *c = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (c, "check_inverse_by_sub");

  merge_matrix (n1, n2, 1.0, c_ll, 1.0, c_lh, 1.0, c_hl, 1.0, c_hh, n, c);

  double *I = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (I, "check_inverse_by_sub");
  mul_matrices (a, n, n, c, n, n, I);

  char label[80];
  double d;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  if (i == j) d = fabs (I[i*n+j] - 1.0);
	  else        d = fabs (I[i*n+j]);

	  sprintf (label, "check_inverse_by_sub : [%d,%d]", i, j);
	  check += compare (d+1.0, 1.0, label, verbose, tiny);
	}
    }

  free (a);
  free (b);
  free (a_ll);
  free (a_lh);
  free (a_hl);
  free (a_hh);
  free (c_ll);
  free (c_lh);
  free (c_hl);
  free (c_hh);
  free (c);
  free (I);

  if (verbose != 0)
    {
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* 
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_solve_gen_linear (int n1, int n2,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_solve_gen_linear : start\n");
    }

  int check = 0;


  int n = n1 + n2;

  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_solve_gen_linear");
  CHECK_MALLOC (b, "check_solve_gen_linear");

  // make symmetric random matrix
  int i, j;
  srand48(0);
  for (i = 0; i < n; i ++)
    {
      for (j = i; j < n; j ++)
	{
	  /*
	  a[i*n+j] = (double)(i*n+j+1);
	  b[i*n+j] = (double)(n*n+i*n+j);
	  */
	  a[i*n+j] = drand48 ();
	  b[i*n+j] = drand48 ();
	}
    }
  for (i = 0; i < n; i ++)
    {
      for (j = i+1; j < n; j ++)
	{
	  a[j*n+i] = a[i*n+j];
	  b[j*n+i] = b[i*n+j];
	}
    }

  // split
  double *a_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *a_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *a_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *a_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (a_ll, "check_solve_gen_linear");
  CHECK_MALLOC (a_lh, "check_solve_gen_linear");
  CHECK_MALLOC (a_hl, "check_solve_gen_linear");
  CHECK_MALLOC (a_hh, "check_solve_gen_linear");

  double *b_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *b_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *b_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *b_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (b_ll, "check_solve_gen_linear");
  CHECK_MALLOC (b_lh, "check_solve_gen_linear");
  CHECK_MALLOC (b_hl, "check_solve_gen_linear");
  CHECK_MALLOC (b_hh, "check_solve_gen_linear");

  double *c_ll = (double *)malloc (sizeof (double) * n1 * n1);
  double *c_lh = (double *)malloc (sizeof (double) * n1 * n2);
  double *c_hl = (double *)malloc (sizeof (double) * n2 * n1);
  double *c_hh = (double *)malloc (sizeof (double) * n2 * n2);
  CHECK_MALLOC (c_ll, "check_solve_gen_linear");
  CHECK_MALLOC (c_lh, "check_solve_gen_linear");
  CHECK_MALLOC (c_hl, "check_solve_gen_linear");
  CHECK_MALLOC (c_hh, "check_solve_gen_linear");

  // case 1
  split_matrix (n, a, n1, n2, a_ll, a_lh, a_hl, a_hh);
  split_matrix (n, b, n1, n2, b_ll, b_lh, b_hl, b_hh);

  solve_gen_linear (n1, n2,
		    a_ll, a_lh, a_hl, a_hh,
		    b_ll, b_lh, b_hl, b_hh,
		    c_ll, c_lh, c_hl, c_hh);
  double *c1 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (c1, "check_solve_gen_linear");
  merge_matrix (n1, n2, 1.0, c_ll, 1.0, c_lh, 1.0, c_hl, 1.0, c_hh, n, c1);


  // case 2
  split_matrix (n, a, n1, n2, a_ll, a_lh, a_hl, a_hh);
  split_matrix (n, b, n1, n2, b_ll, b_lh, b_hl, b_hh);

  solve_gen_linear_ (n1, n2,
		     a_ll, a_lh, a_hl, a_hh,
		     b_ll, b_lh, b_hl, b_hh,
		     c_ll, c_lh, c_hl, c_hh);
  double *c2 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (c2, "check_solve_gen_linear");
  merge_matrix (n1, n2, 1.0, c_ll, 1.0, c_lh, 1.0, c_hl, 1.0, c_hh, n, c2);

  // compare
  char label[80];
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "check_solve_gen_linear : c1-c2[%d,%d]", i, j);
	  check += compare (c1[i*n+j], c2[i*n+j], label, verbose, tiny);
	}
    }


  double *vb = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (vb, "check_solve_gen_linear");
  for (i = 0; i < n; i ++)
    {
      vb[i] = 1.0;
    }

  double *vx = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (vx, "check_solve_gen_linear");

  // use the result by solve_gen_linear_()
  dot_prod_matrix (c2, n, n, vb, vx);


  // check the result
  double *xy = (double *)malloc (sizeof (double) * n);
  double *bc = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (xy, "check_solve_gen_linear");
  CHECK_MALLOC (bc, "check_solve_gen_linear");

  for (i = 0; i < n1; i ++)
    {
      xy [i] = vx [i];
      bc [i] = vb [i];
    }
  for (i = n1; i < n; i ++)
    {
      xy [i] = vb [i];
      bc [i] = vx [i];
    }

  double *left = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (left, "check_solve_gen_linear");
  dot_prod_matrix (a, n, n, xy, left);

  double *right = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (right, "check_solve_gen_linear");
  dot_prod_matrix (b, n, n, bc, right);

  for (i = 0; i < n; i ++)
    {
      sprintf (label, "check_solve_gen_linear : left-right[%d]", i);
      check += compare (left[i], right[i], label, verbose, tiny);
    }

  free (a);
  free (b);
  free (a_ll);
  free (a_lh);
  free (a_hl);
  free (a_hh);
  free (b_ll);
  free (b_lh);
  free (b_hl);
  free (b_hh);
  free (c_ll);
  free (c_lh);
  free (c_hl);
  free (c_hh);
  free (c1);
  free (c2);
  free (vb);
  free (vx);
  free (xy);
  free (bc);
  free (left);
  free (right);

  if (verbose != 0)
    {
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
