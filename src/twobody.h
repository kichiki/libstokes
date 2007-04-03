/* header file for [XYZ][ABCGHM].c and
 * RYUON-twobody : exact 2-body resistance scalar functions
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: twobody.h,v 1.3 2007/04/03 02:35:01 kichiki Exp $
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
#ifndef	_TWOBODY_H_
#define	_TWOBODY_H_


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
void twobody_far (int n, double l, double s,
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

/* calc scalar functions of resistance problem by lub form
 * INPUT
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
 *     10,11 : (XG11, XG12) not implemented
 *     12,13 : (YG11, YG12) not implemented
 *     14,15 : (YH11, YH12) not implemented
 *     16,17 : (XM11, XM12) not implemented
 *     18,19 : (YM11, YM12) not implemented
 *     20,21 : (ZM11, ZM12) not implemented
 */
void twobody_lub (int n, double l, double s,
		  double *far);


/* scale the scalar functions from Jeffrey-Onishi to Stokesian dynamics
 * INPUT
 *  two [22] : scalar functions in Jeffrey form
 *  l        : lambda = ab / aa,
 *             where aa and ab are radii for particles a(alpha) and b(beta)
 *             Note that the scalar functions are for "a-b" interaction.
 * OUTPUT
 *  two [22] : scalar functions in SD form
 */
void
twobody_scale (double *two, double l);

/* calc scalar functions of two-body exact solution in resistance problem
 * INPUT
 *  r        : distance between the two := x_b - x_a
 *  aa, ab   : radii for particles a(alpha) and b(beta)
 *  n        : max order for the coefficients
 *  flag_lub : 0 to use twobody_far()
 *             1 to use twobody_lub()
 *  res [22] :
 * OUTPUT
 *  res [22] : scalar functions are scaled by Stokesian dynamics form.
 */
void
twobody_scalars_res (double r,
		     double aa, double ab,
		     int n, int flag_lub,
		     double *res);


#endif /* !_TWOBODY_H_ */
