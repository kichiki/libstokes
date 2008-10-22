/* twobody solutions for slip particles
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: XC-slip.c,v 1.3 2008/10/22 05:54:25 kichiki Exp $
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
#include <stdio.h>  // fprintf()
#include <stdlib.h> // malloc(), free()
#include "memory-check.h" // CHECK_MALLOC

double SL_G1 (int m, int n);
double SL_G2 (int m, int n);

/* XC1  */
void twobody_XC_slip (int n, double l, double *f)
{
  int i;
  for (i = 0; i <= n; i ++) f[i] = 0.0;
  int q;
  double *fk = (double *)malloc (sizeof (double) * (n+1));
  CHECK_MALLOC (fk, "twobody_XC_slip");

  /* k = 0 */
  for (q = 0; q <= 0; q++) fk[q] = 0.0;
  fk[0] =
      SL_G1(0,3) 
  ;
  f[0] = 0.0;
  for (q = 0; q <= 0; q++) f[0] = f[0] * l + fk[0-q];
  if (n == 0)
    {
      free (fk);
      return;
    }

  /* k = 3 */
  for (q = 0; q <= 3; q++) fk[q] = 0.0;
  fk[3] =
      8.0 * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[3] = 0.0;
  for (q = 0; q <= 3; q++) f[3] = f[3] * l + fk[3-q];
  if (n == 3)
    {
      free (fk);
      return;
    }

  /* k = 6 */
  for (q = 0; q <= 6; q++) fk[q] = 0.0;
  fk[3] =
      64.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[6] = 0.0;
  for (q = 0; q <= 6; q++) f[6] = f[6] * l + fk[6-q];
  if (n == 6)
    {
      free (fk);
      return;
    }

  /* k = 8 */
  for (q = 0; q <= 8; q++) fk[q] = 0.0;
  fk[5] =
      768.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[8] = 0.0;
  for (q = 0; q <= 8; q++) f[8] = f[8] * l + fk[8-q];
  if (n == 8)
    {
      free (fk);
      return;
    }

  /* k = 9 */
  for (q = 0; q <= 9; q++) fk[q] = 0.0;
  fk[6] =
      512.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  f[9] = 0.0;
  for (q = 0; q <= 9; q++) f[9] = f[9] * l + fk[9-q];
  if (n == 9)
    {
      free (fk);
      return;
    }

  /* k = 10 */
  for (q = 0; q <= 10; q++) fk[q] = 0.0;
  fk[7] =
      6144.0 * SL_G2(-2,5) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[10] = 0.0;
  for (q = 0; q <= 10; q++) f[10] = f[10] * l + fk[10-q];
  if (n == 10)
    {
      free (fk);
      return;
    }

  /* k = 11 */
  for (q = 0; q <= 11; q++) fk[q] = 0.0;
  fk[6] =
      6144.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      6144.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[11] = 0.0;
  for (q = 0; q <= 11; q++) f[11] = f[11] * l + fk[11-q];
  if (n == 11)
    {
      free (fk);
      return;
    }

  /* k = 12 */
  for (q = 0; q <= 12; q++) fk[q] = 0.0;
  fk[6] =
      4096.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      40960.0 * SL_G2(-3,6) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[12] = 0.0;
  for (q = 0; q <= 12; q++) f[12] = f[12] * l + fk[12-q];
  if (n == 12)
    {
      free (fk);
      return;
    }

  /* k = 13 */
  for (q = 0; q <= 13; q++) fk[q] = 0.0;
  fk[6] =
      49152.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      98304.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      49152.0 * SL_G2(-2,5) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[13] = 0.0;
  for (q = 0; q <= 13; q++) f[13] = f[13] * l + fk[13-q];
  if (n == 13)
    {
      free (fk);
      return;
    }

  /* k = 14 */
  for (q = 0; q <= 14; q++) fk[q] = 0.0;
  fk[6] =
      49152.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      98304.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[11] =
      245760.0 * SL_G2(-4,7) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[14] = 0.0;
  for (q = 0; q <= 14; q++) f[14] = f[14] * l + fk[14-q];
  if (n == 14)
    {
      free (fk);
      return;
    }

  /* k = 15 */
  for (q = 0; q <= 15; q++) fk[q] = 0.0;
  fk[6] =
      327680.0 * SL_G1(-3,6) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      983040.0 * SL_G1(-2,5) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      32768.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      983040.0 * SL_G2(-2,5) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[12] =
      327680.0 * SL_G2(-3,6) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[15] = 0.0;
  for (q = 0; q <= 15; q++) f[15] = f[15] * l + fk[15-q];
  if (n == 15)
    {
      free (fk);
      return;
    }

  /* k = 16 */
  for (q = 0; q <= 16; q++) fk[q] = 0.0;
  fk[6] =
      393216.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      1572864.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      196608.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 4.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 3.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
  ;
  fk[13] =
      1376256.0 * SL_G2(-5,8) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[16] = 0.0;
  for (q = 0; q <= 16; q++) f[16] = f[16] * l + fk[16-q];
  if (n == 16)
    {
      free (fk);
      return;
    }

  /* k = 17 */
  for (q = 0; q <= 17; q++) fk[q] = 0.0;
  fk[6] =
      1966080.0 * SL_G1(-4,7) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      7864320.0 * SL_G1(-3,6) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      786432.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      11796480.0 * SL_G1(-2,5) * SL_G2(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[11] =
      786432.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[12] =
      7864320.0 * SL_G2(-3,6) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[14] =
      1966080.0 * SL_G2(-4,7) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[17] = 0.0;
  for (q = 0; q <= 17; q++) f[17] = f[17] * l + fk[17-q];
  if (n == 17)
    {
      free (fk);
      return;
    }

  /* k = 18 */
  for (q = 0; q <= 18; q++) fk[q] = 0.0;
  fk[6] =
      2621440.0 * SL_G1(-3,6) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      15728640.0 * SL_G1(-2,5) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      262144.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      3145728.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 4.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
  ;
  fk[12] =
      1048576.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-3,6) * SL_G2(0,3) 
      + 9.0 * SL_G2(-2,5) * SL_G2(-1,4) ) 
  ;
  fk[15] =
      7340032.0 * SL_G2(-6,9) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[18] = 0.0;
  for (q = 0; q <= 18; q++) f[18] = f[18] * l + fk[18-q];
  if (n == 18)
    {
      free (fk);
      return;
    }

  /* k = 19 */
  for (q = 0; q <= 19; q++) fk[q] = 0.0;
  fk[6] =
      11010048.0 * SL_G1(-5,8) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      55050240.0 * SL_G1(-4,7) * SL_G2(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      1572864.0 * SL_G1(0,3) * ( 4.0 * SL_G1(-2,5) * SL_G1(0,3) 
      + 3.0 * SL_G1(-1,4) * SL_G1(-1,4) ) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      110100480.0 * SL_G1(-3,6) * SL_G2(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[11] =
      23592960.0 * SL_G1(-1,4) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[12] =
      110100480.0 * SL_G2(-3,6) * SL_G1(-2,5) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[13] =
      1572864.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * ( 4.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 3.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
  ;
  fk[14] =
      55050240.0 * SL_G2(-4,7) * SL_G1(-1,4) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[16] =
      11010048.0 * SL_G2(-5,8) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  f[19] = 0.0;
  for (q = 0; q <= 19; q++) f[19] = f[19] * l + fk[19-q];
  if (n == 19)
    {
      free (fk);
      return;
    }

  /* k = 20 */
  for (q = 0; q <= 20; q++) fk[q] = 0.0;
  fk[6] =
      15728640.0 * SL_G1(-4,7) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[8] =
      125829120.0 * SL_G1(-3,6) * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) 
  ;
  fk[9] =
      6291456.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[10] =
      31457280.0 * SL_G1(-2,5) * SL_G1(0,3) * SL_G1(0,3) * ( 6.0 * SL_G2(-2,5) * SL_G2(0,3) 
      + 5.0 * SL_G2(-1,4) * SL_G2(-1,4) ) 
  ;
  fk[11] =
      9437184.0 * SL_G2(-1,4) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * SL_G2(0,3) * SL_G2(0,3) 
  ;
  fk[12] =
      125829120.0 * SL_G1(-1,4) * SL_G1(0,3) * SL_G1(0,3) * ( SL_G2(-3,6) * SL_G2(0,3) 
      + 2.0 * SL_G2(-2,5) * SL_G2(-1,4) ) 
  ;
  fk[14] =
      6291456.0 * SL_G1(0,3) * SL_G1(0,3) * SL_G1(0,3) * ( 5.0 * SL_G2(-4,7) * SL_G2(0,3) 
      + 10.0 * SL_G2(-3,6) * SL_G2(-1,4) 
      + 6.0 * SL_G2(-2,5) * SL_G2(-2,5) ) 
  ;
  fk[17] =
      37748736.0 * SL_G2(-7,10) * SL_G1(0,3) * SL_G1(0,3) 
  ;
  f[20] = 0.0;
  for (q = 0; q <= 20; q++) f[20] = f[20] * l + fk[20-q];
  if (n == 20)
    {
      free (fk);
      return;
    }

  fprintf (stderr, "twobody_XC_slip is implemented only up to k=20\n");
  free (fk);
  return;
}
