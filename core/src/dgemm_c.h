/* header file for dgemm.c --
 * C wrappers for BLAS' dgemm_()
 * Copyright (C) 2005-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: dgemm_c.h,v 1.1 2006/10/22 22:43:44 kichiki Exp $
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
#ifndef	_DGEMM_C_H_
#define	_DGEMM_C_H_

/* C := alpha*op( A )*op( B ) + beta*C
 * INPUT
 *  m : the number of rows   of the matrix op(A),op(C)=> op(A)[ij],i=0~m.
 *  n : the number of colums of the matrix op(B),op(C)=> op(B)[ij],j=0~n.
 *  k : the number of colums of the matrix op(A) and
 *      the number of rows   of the matrix op(B)      => op(A)[ij],j=0~k.
 *  A[m,k] . B[k,n] = C[m,n]
 *  op(A) is 'm by k' matrix => A[m,k]
 *  op(B)    'k by n' matrix => B[k,n]
 *  op(C)    'm by n' matrix => C[m,n]
 *
 *  A := [a(1,1) a(1,2) ... a(1,n)]
 *       [a(2,1)            a(2,n)]
 *       [ :                      ]
 *       [a(m,1) a(m,2) ... a(m,n)]
 *  A is 'm by n' matrix, where
 *       where n is the number of columns (colum means TATE no NARABI)
 *       and   m is the number of rows.   (row   means 'line')
 * in C, op(A)_{ij} := a[i*k+j], with i=(0~m-1), j=(0~k-1),
 *       op(B)_{ij} := b[i*n+j], with i=(0~k-1), j=(0~n-1),
 *       op(C)_{ij} := c[i*n+j], with i=(0~m-1), j=(0~n-1).
 */
void dgemm_wrap (int m, int n, int k,
		 double alpha, double *a, double *b,
		 double beta, double *c);

#endif /* !_DGEMM_C_H_ */
