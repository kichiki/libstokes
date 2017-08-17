/* header file for matrix.c --
 * matrix-manipulating routines
 * Copyright (C) 2001-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: matrix.h,v 1.6 2007/11/28 03:26:38 kichiki Exp $
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
#ifndef	_MATRIX_H_
#define	_MATRIX_H_

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


#endif /* !_MATRIX_H_ */
