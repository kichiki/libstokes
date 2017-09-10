/* header file for f-old.c --
 * backup of bug fixing for polydisperse systems
 * subroutine for the procedure of F version
 * Copyright (C) 2001-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_F_OLD_H_
#define	_F_OLD_H_


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
calc_lub_f_2b_poly_old
(struct stokes *sys,
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
 *            (i,j) are used to assign the results in mat[].
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *            used to access sys->a[] and sys->poly_table[]
 *   n      : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b_poly_old
(struct stokes *sys,
 int i, int j,
 const double *x1, const double *x2,
 int i1, int i2,
 int n, double *mat);


#endif /* !_F_OLD_H_ */
