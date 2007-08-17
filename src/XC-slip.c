/* twobody solutions for slip particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: XC-slip.c,v 1.1 2007/08/17 04:31:17 kichiki Exp $
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

double SL_G1 (int m, int n);
double SL_G2 (int m, int n);

/* XC1  */
void twobody_XC_slip (int n, double l, double *f)
{
  int i;
  for (i = 0; i <= n; i ++) f[i] = 0.0;

  f [0] =
    + ( // q = 0
      SL_G1(0,3) 
    ) ;
  if (n == 2) return;
  f [3] =
    + l * l * l * ( // q = 3
      8.0 * SL_G1(0,3) * SL_G2(0,3) 
    ) ;
  if (n == 5) return;
  f [6] =
    + l * l * l * ( // q = 3
      64.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ;
  if (n == 7) return;
  f [8] =
    + l * l * l * l * l * ( // q = 5
      768.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) 
    ) ;
  if (n == 8) return;
  f [9] =
    + l * l * l * l * l * l * ( // q = 6
      512.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    ) ;
  if (n == 9) return;
  f [10] =
    + l * l * l * l * l * l * l * ( // q = 7
      6144.0 * SL_G2(-2,5) * SL_G1(0,3) * SL_G1(0,3) 
    ) ;
  if (n == 10) return;
  f [11] =
    + l * l * l * l * l * l * ( // q = 6
      6144.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      6144.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ) ;
  if (n == 11) return;
  f [12] =
    + l * l * l * l * l * l * ( // q = 6
      4096.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * l * ( // q = 9
      40960.0 * SL_G2(-3,6) * SL_G1(0,3) * SL_G1(0,3) 
    ) ) ;
  if (n == 12) return;
  f [13] =
    + l * l * l * l * l * l * ( // q = 6
      49152.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      98304.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 10
      49152.0 * SL_G2(-2,5) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ) ) ;
  if (n == 13) return;
  f [14] =
    + l * l * l * l * l * l * ( // q = 6
      49152.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      98304.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * l * ( // q = 11
      245760.0 * SL_G2(-4,7) * SL_G1(0,3) * SL_G1(0,3) 
    ) ) ) ;
  if (n == 14) return;
  f [15] =
    + l * l * l * l * l * l * ( // q = 6
      327680.0 * SL_G1(-3,6) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      983040.0 * SL_G1(-2,5) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 9
      32768.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 10
      983040.0 * SL_G2(-2,5) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 12
      327680.0 * SL_G2(-3,6) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ) ) ) ) ;
  if (n == 15) return;
  f [16] =
    + l * l * l * l * l * l * ( // q = 6
      393216.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      1572864.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 10
      196608.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 4.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 3.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
    + l * l * l * ( // q = 13
      1376256.0 * SL_G2(-5,8) * SL_G1(0,3) * SL_G1(0,3) 
    ) ) ) ) ;
  if (n == 16) return;
  f [17] =
    + l * l * l * l * l * l * ( // q = 6
      1966080.0 * SL_G1(-4,7) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      7864320.0 * SL_G1(-3,6) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 9
      786432.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 10
      11796480.0 * SL_G1(-2,5) * SL_G2(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 11
      786432.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 12
      7864320.0 * SL_G2(-3,6) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 14
      1966080.0 * SL_G2(-4,7) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ) ) ) ) ) ) ;
  if (n == 17) return;
  f [18] =
    + l * l * l * l * l * l * ( // q = 6
      2621440.0 * SL_G1(-3,6) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      15728640.0 * SL_G1(-2,5) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 9
      262144.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 10
      3145728.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 4.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
    + l * l * ( // q = 12
      1048576.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-3,6) * SL_G2(0,3) 
      + 9.0 * SL_G2(-2,5) * SL_G2(-1,4) ) 
    + l * l * l * ( // q = 15
      7340032.0 * SL_G2(-6,9) * SL_G1(0,3) * SL_G1(0,3) 
    ) ) ) ) ) ) ;
  if (n == 18) return;
  f [19] =
    + l * l * l * l * l * l * ( // q = 6
      11010048.0 * SL_G1(-5,8) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      55050240.0 * SL_G1(-4,7) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 9
      1572864.0 * SL_G1(0,3) * ( 4.0 * SL_G1(-2,5) * SL_G1(0,3) 
      + 3.0 * SL_G1(-1,4) * SL_G1(-1,4) ) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 10
      110100480.0 * SL_G1(-3,6) * SL_G2(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 11
      23592960.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 12
      110100480.0 * SL_G2(-3,6) * SL_G1(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 13
      1572864.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * ( 4.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 3.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
    + l * ( // q = 14
      55050240.0 * SL_G2(-4,7) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 16
      11010048.0 * SL_G2(-5,8) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    ) ) ) ) ) ) ) ) ) ;
  if (n == 19) return;
  f [20] =
    + l * l * l * l * l * l * ( // q = 6
      15728640.0 * SL_G1(-4,7) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * l * ( // q = 8
      125829120.0 * SL_G1(-3,6) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
    + l * ( // q = 9
      6291456.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 10
      31457280.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G1(0,3) * ( 6.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 5.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
    + l * ( // q = 11
      9437184.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
    + l * ( // q = 12
      125829120.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * ( SL_G2(-3,6) * SL_G2(0,3) 
      + 2.0 * SL_G2(-2,5) * SL_G2(-1,4) ) 
    + l * l * ( // q = 14
      6291456.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-4,7) * SL_G2(0,3) 
      + 10.0 * SL_G2(-3,6) * SL_G2(-1,4) 
      + 6.0 * SL_G2(-2,5) * SL_G2(-2,5) ) 
    + l * l * l * ( // q = 17
      37748736.0 * SL_G2(-7,10) * SL_G1(0,3) * SL_G1(0,3) 
    ) ) ) ) ) ) ) ) ;
}
