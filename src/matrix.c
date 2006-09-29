/* matrix-manipulating routines
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: matrix.c,v 1.2 2006/09/29 03:38:48 ichiki Exp $
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
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
//#include "../../misc/ludcmp.h" /* lu_inv() */
#include "dgetri_c.h" /* lapack_inv_() */

#include "matrix.h"

/* (local) global variable */


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
		  double * I, double * J, double * K, double * L)
{
  /* type 2 */

  /* H := H^-1 */
  //lu_inv (H, n2);
  lapack_inv_ (n2, H);

  /* C := (H^-1) . C */
  mul_left_sq (C, n2, n1, H, n2);
  /* G := (H^-1) . G */
  mul_left_sq (G, n2, n1, H, n2);
  /* D := (H^-1) . D */
  mul_left_sq (D, n2, n2, H, n2);

  /* A [n1, n1] := A - F.(H^-1).C */
  add_and_mul (A, n1, n1, F, n1, n2, C, n2, n1, 1.0, -1.0, A);
  /* E [n1, n1] := E - F.(H^-1).G */
  add_and_mul (E, n1, n1, F, n1, n2, G, n2, n1, 1.0, -1.0, E);
  /* B [n1, n2] := - B + F.(H^-1).D */
  add_and_mul (B, n1, n2, F, n1, n2, D, n2, n2, -1.0, 1.0, B);

  /* A := (A-F.(H^-1).C)^-1 */
  //lu_inv (A, n1);
  lapack_inv_ (n1, A);

  /* I := A.E */
  mul_matrices (A, n1, n1, E, n1, n1, I);
  /* J := A.B */
  mul_matrices (A, n1, n1, B, n1, n2, J);

  /* K := - G + C.I */
  add_and_mul (G, n2, n1, C, n2, n1, I, n1, n1, -1.0, 1.0, K);
  /* L := D + C.J */
  add_and_mul (D, n2, n2, C, n2, n1, J, n1, n2, 1.0, 1.0, L);
}

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
	      double * I, double * J, double * K, double * L)
{
  int i;
  /* type 1 */

  /* A := A^-1 */
  //lu_inv (A, n1);
  lapack_inv_ (n1, A);

  /* B := (A^-1).B */
  mul_left_sq (B, n1, n2, A, n1);

  /* L [n2, n2] := D - C.(A^-1).B */
  add_and_mul (D, n2, n2, C, n2, n1, B, n1, n2, 1.0, -1.0, L);

  /* K := C.A */
  mul_matrices (C, n2, n1, A, n1, n1, K);

  for (i = 0; i < n1 * n1; ++i)
    I [i] = A [i];
  for (i = 0; i < n1 * n2; ++i)
    J [i] = - B [i];
}

/*
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
	      double * C)
{
  int i, j, k;


  if (na2 != nb1)
    {
      fprintf (stderr, "illeagal multiplication in mul_left().\n");
      exit (1);
    }

  for (i = 0; i < na1; ++i)
    {
      for (j = 0; j < nb2; ++j)
	{
	  C [i * nb2 + j] = 0.0;
	  for (k = 0; k < na2; ++k)
	    {
	      C [i * nb2 + j] +=
		A [i * na2 + k] * B [k * nb2 + j];
	    }
	}
    }
}

/*
 * INPUT
 *  A [na1, na2]
 *  B [nb, nb] : square matrix!
 *  where (nb == na1)
 * OUTPUT
 *  A [na1, na2] = B [nb, nb] . A [na1, na2]
 */
void
mul_left_sq (double * A, int na1, int na2,
	     const double * B, int nb)
{
  int i, j, k;
  double * tmp;


  if (nb != na1)
    {
      fprintf (stderr, "illeagal multiplication in mul_left_sq().\n");
      exit (1);
    }

  tmp = (double *) malloc (sizeof (double) * na1 * na2);
  if (tmp == NULL)
    {
      fprintf (stderr, "allocation error in mul_left_sq().\n");
      exit (1);
    }

  for (i = 0; i < na1; ++i)
    {
      for (j = 0; j < na2; ++j)
	{
	  tmp [i * na2 + j] = 0.0;
	  for (k = 0; k < na1; ++k)
	    {
	      tmp [i * na2 + j] +=
		B [i * nb + k]
		* A [k * na2 + j];
	    }
	}
    }

  for (i = 0; i < na1 * na2; ++i)
    A [i] = tmp [i];

  free (tmp);
}

/*
 * INPUT
 *  A [na1, na2]
 *  B [nb1, nb2]
 *  C [nc1, nc2]
 *  a
 *  b
 *  where (nb2 == nc1) && (na1 = nb1) && (na2 == nc2)
 * OUTPUT
 *  D [na1, na2] = a * A [na1, na2] + b * B [nb1, nb2] . C [nc1, nc2]
 *  D could be same to A itself!
 */
void
add_and_mul (const double * A, int na1, int na2,
	     const double * B, int nb1, int nb2,
	     const double * C, int nc1, int nc2,
	     double a, double b,
	     double * D)
{
  int i, j, k;


  if (nb2 != nc1
      || na1 != nb1
      || na2 != nc2)
    {
      fprintf (stderr, "illeagal multiplication in add_and_mul().\n");
      exit (1);
    }

  for (i = 0; i < na1; ++i)
    {
      for (j = 0; j < na2; ++j)
	{
	  D [i * na2 + j] = a * A [i * na2 + j];
	  for (k = 0; k < nc1; ++k)
	    {
	      D [i * na2 + j] +=
		b * B [i * nb2 + k] * C [k * nc2 + j];
	    }
	}
    }
}

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
		 double * y)
{
  int i, j;

  for (i = 0; i < n1; ++i)
    {
      y [i] = 0.0;
      for (j = 0; j < n2; ++j)
	{
	  y [i] += mat [i * n2 + j] * x [j];
	}
    }
}
