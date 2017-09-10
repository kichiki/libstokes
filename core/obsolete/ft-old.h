/* header file for ft-old.c --
 * backup of bug fixing for polydisperse systems
 * subroutine for the procedure of FT version
 * Copyright (C) 2000-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_FT_OLD_H_
#define	_FT_OLD_H_


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
calc_lub_ft_2b_poly_old
(struct stokes *sys,
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
matrix_lub_ft_2b_poly_old
(struct stokes *sys,
 int i, int j,
 const double *x1, const double *x2,
 int i1, int i2,
 int n, double *mat);


#endif /* !_FT_OLD_H_ */
