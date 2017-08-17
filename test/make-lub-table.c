/* make lub table for comparison with SD code
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: $
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

#include "stokes.h" // struct stokes
#include "fts.h" // scalar_minv_fts()
#include "minv-poly.h"
#include "two-body-res.h" // scalar_two_body_res()
#include "twobody.h" // twobody_scalars_res()

#include "memory-check.h"


/* calc scalar functions of (M^inf)^-1 in FTS
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   lub [22] : scalar functions
 */
void
scalar_minv_fts_old (double s,  double * scalar_fts)
{
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

  double s2, s3, s4, s5, s6, s7, s8, s10;
  double gx, gy;


  s2 = s * s;
  s3 = s * s2;
  s4 = s2 * s2;
  s5 = s2 * s3;
  s6 = s2 * s4;
  s7 = s * s6;
  s8 = s4 * s4;
  s10 = s5 * s5;
  
  gx = 
    (2304.0 + s2 *
     (- 21120.0 + s2 *
      (55600.0 + s2 *
       (- 90600.0 + s2 *
	(45945.0 + s2 *
	 (- 800.0 + s2 *
	  (- 1800.0 + s2 *
	   (- 900.0 + s2 *
	    (400.0)))))))));
  gy =
    (256.0 + s2 *
     (640.0 + s2 *
      (- 2000.0 + s2 *
       (- 29624.0 + s2 *
	(16505.0 + s2 *
	 (- 14880.0 + s2 *
	  (112600.0 + s2 *
	   (- 66360 + s2 *
	    (22800.0 + s2 *
	     (3600.0 + s2 *
	      (900.0 + s2 *
	       (- 1600.0))))))))))));
  xa11 =
    20.0 * s6 *
    (- 2880.0 + s2 *
     (2208.0 + s2 *
      (- 260.0 + s2 *
       (- 75.0 + s2  * s2*
	(20.0))))) / gx;
  xa12 =
    20.0 * s3 *
    (- 576.0 + s2 *
     (2880.0 + s2 *
      (- 2000.0 + s2 *
       (375.0 + s2  * s2*
	(20.0 + s2 *
	 (- 30.0)))))) / gx;
  ya11 =
    - 80.0 * s6 *
    (320.0 + s2 *
     (- 544.0 + s2 *
      (140.0 + s2 *
       (- 1130.0 + s2 *
	(736.0 + s2 *
	 (- 280.0 + s2 *
	  (- 15.0 + s2  * s2*
	   (20.0)))))))) / gy;
  ya12 =
    40.0 * s3 *
    (64.0 + s2  * s2*
     (- 400.0 + s2 *
      (119.0 + s2 *
       (- 1440.0 + s2 *
	(1160.0 + s2 *
	 (- 405.0 + s2  * s2*
	  (20.0 + s2 *
	   (30.0)))))))) / gy;
  yb11 =
    - 80.0 * s7 *
    (384.0 + s2 *
     (- 248.0 + s2 *
      (- 190.0 + s2 *
       (150.0 + s2 *
	(- 80.0 + s2 *
	 (- 20.0 + s2 *
	  (- 15.0))))))) / gy;
  yb12 =
    - 20.0 * s4 *
    (- 128.0 + s2 *
     (- 80.0 + s2 *
      (700.0 + s2 *
       (- 2935.0 + s2 *
	(2464.0 + s2 *
	 (- 540.0 + s2 *
	  (- 30.0 + s2  * s2*
	   (80.0)))))))) / gy;
  xc11 =
    4.0 * s6
    / 3.0 / (s6 - 1.0);
  xc12 =
    - 4.0 * s3
    / 3.0 / (s6 - 1.0);
  yc11 =
    - 16.0 * s6 *
    (256.0 + s2 *
     (7840.0 + s2 *
      (1960.0 + s2 *
       (- 31525.0 + s2 *
	(15690.0 + s2 *
	 (- 4100.0 + s2 *
	  (- 600.0 + s2 *
	   (- 225.0 + s2 *
	    (400.0)))))))))
    / 3.0 / gy;
  yc12 =
    - 4.0 * s3 *
    (512.0 + s2 *
     (3680.0 + s2 *
      (200.0 + s2 *
       (- 66950.0 + s2 *
	(80505.0 + s2 *
	 (- 10600.0 + s2  * s2*
	  (450.0 + s2 *
	   (800.0))))))))
    / 3.0 / gy;
  xg11 =
    100.0 * s7 *
    (192.0 + s2 *
     (- 184.0 + s2 *
      (16.0 + s2 *
       (15.0)))) / gx;
  xg12 =
    10.0 * s4 *
    (384.0 + s2 *
     (- 2000.0 + s2 *
      (1700.0 + s2 *
       (- 375.0 + s2 *
	(160.0 + s2 *
	 (- 100.0)))))) / gx;
  yg11 =
    - 100.0 * s7 *
    (128.0 + s2 *
     (320.0 + s2 *
      (- 48.0 + s2 *
       (377.0 + s2 *
	(- 128.0 + s2 *
	 (108.0))))))
    / 3.0 / gy;
  yg12 =
    - 10.0 * s4 *
    (128.0 + s2 *
     (- 80.0 + s2 *
      (- 900.0 + s2 *
       (613.0 + s2 *
	(- 5280.0 + s2 *
	 (2580.0 + s2 *
	  (- 450.0 + s2 *
	   (- 640.0))))))))
    / 3.0 / gy;
  yh11 =
    - 20.0 * s8 *
    (736.0 + s2 *
     (- 1240.0 + s2 *
      (- 1020.0 + s2 *
       (3875.0 + s2 *
	(- 480.0)))))
    / 3.0 / gy;
  yh12 =
    - 10.0 * s5 *
    (96.0 + s2 *
     (- 120.0 + s2 *
      (- 750.0 + s2 *
       (6885.0 + s2 *
	(- 5160.0 + s2 *
	 (- 720.0 + s2 *
	  (- 450.0 + s2 *
	   (800.0))))))))
    / 3.0 / gy;
  xm11 =
    200.0 * s8 *
    (- 192.0 + s2 *
     (220.0 + s2 *
      (- 15.0 + s2 *
       (- 45.0 + s2 *
	(20.0)))))
    / 9.0 / gx;
  xm12 =
    100.0 * s5 *
    (96.0 + s2 *
     (- 584.0 + s2 *
      (810.0 + s2 *
       (- 705.0 + s2 *
	(200.0)))))
    / 9.0 / gx;
  ym11 =
    - 200.0 * s8 *
    (64.0 + s2 *
     (- 208.0 + s2 *
      (75.0 + s2 *
       (- 76.0 + s2 *
	(- 340.0 + s2 *
	 (- 180.0 + s2 *
	  (- 45.0 + s2 *
	   (80.0))))))))
    / 9.0 / gy;
  ym12 =
    - 50.0 * s5 *
    (32.0 + s2 *
     (- 8.0 + s2 *
      (- 210.0 + s2 *
       (- 543.0 + s2 *
	(- 2072.0 + s2 *
	 (360.0 + s2 *
	  (3010.0 + s2 *
	   (- 800.0))))))))
    / 9.0 / gy;
  zm11 =
    10.0 * s10
    / 9.0 / (s10 - 4.0);
  zm12 =
    - 20.0 * s5
    / 9.0 / (s10 - 4.0);

  scalar_fts [ 0] = xa11;
  scalar_fts [ 1] = xa12;
  scalar_fts [ 2] = ya11;
  scalar_fts [ 3] = ya12;
  scalar_fts [ 4] = yb11;
  scalar_fts [ 5] = yb12;
  scalar_fts [ 6] = xc11;
  scalar_fts [ 7] = xc12;
  scalar_fts [ 8] = yc11;
  scalar_fts [ 9] = yc12;
  scalar_fts [10] = xg11;
  scalar_fts [11] = xg12;
  scalar_fts [12] = yg11;
  scalar_fts [13] = yg12;
  scalar_fts [14] = yh11;
  scalar_fts [15] = yh12;
  scalar_fts [16] = xm11;
  scalar_fts [17] = xm12;
  scalar_fts [18] = ym11;
  scalar_fts [19] = ym12;
  scalar_fts [20] = zm11;
  scalar_fts [21] = zm12;
}


static void
print_scalars (FILE *out, char *label, double r, double *scalars)
{
  fprintf (stdout, "%s %6.3f %e %e %e %e %e %e %e %e %e %e"
	   " %e %e %e %e %e %e %e %e %e %e %e %e\n",
	   label, r,
	   scalars [0],
	   scalars [1],
	   scalars [2],
	   scalars [3],
	   scalars [4],
	   scalars [5],
	   scalars [6],
	   scalars [7],
	   scalars [8],
	   scalars [9],
	   scalars[10],
	   scalars[11],
	   scalars[12],
	   scalars[13],
	   scalars[14],
	   scalars[15],
	   scalars[16],
	   scalars[17],
	   scalars[18],
	   scalars[19],
	   scalars[20],
	   scalars[21]
	   );
}

static void
print_scalars_full (FILE *out, char *label, double r, double *scalars)
{
  fprintf (stdout, "%s %6.3f %e %e %e %e %e %e %e %e %e %e"
	   " %e %e %e %e %e %e %e %e %e %e %e %e\n",
	   label, r,
	   scalars [0],
	   scalars [1],
	   scalars [4],
	   scalars [5],
	   scalars [8],
	   scalars [9],
	   scalars[12],
	   scalars[13],
	   scalars[16],
	   scalars[17],
	   scalars[20],
	   scalars[21],
	   scalars[24],
	   scalars[25],
	   scalars[28],
	   scalars[29],
	   scalars[32],
	   scalars[33],
	   scalars[36],
	   scalars[37],
	   scalars[40],
	   scalars[41]
	   );
}


/* main program */
int
main (int argc, char** argv)
{
  double scalars [44];
  double r;
  double a1 = 1.5;
  double a2 = 0.5;

  int i;
  for (i = 0; i < 39; i ++)
    {
      r = 2.1 + (double)i * 0.05;

      scalar_two_body_res (r, scalars);
      print_scalars (stdout, "2B-old", r, scalars);

      twobody_scalars_res (2, // FTS
			   r, 1.0, 1.0, NULL,
			   100, // n
			   1, // lub
			   1, // dimensional
			   scalars);
      print_scalars (stdout, "2B-new", r, scalars);

      scalar_minv_fts (r, scalars);
      print_scalars (stdout, "minv-old", r, scalars);

      scalar_minv_fts_old (r, scalars);
      print_scalars (stdout, "minv-old2", r, scalars);

      scalars_minv_fts_poly (r, 1.0, 1.0, scalars);
      print_scalars_full (stdout, "minv-new", r, scalars);

      scalars_lub_poly (2, // FTS
			r, 1.0, 1.0, NULL,
			100, // n
			1, // lub
			scalars);
      print_scalars (stdout, "lub-fts", r, scalars);


      // polydisperse case
      twobody_scalars_res (2, // FTS
			   r, a1, a2, NULL,
			   100, // n
			   1, // lub
			   1, // dimensional
			   scalars);
      print_scalars (stdout, "2B-poly", r, scalars);

      scalars_minv_fts_poly (r, a1, a2, scalars);
      print_scalars_full (stdout, "minv-poly", r, scalars);

      scalars_lub_poly (2, // FTS
			r, a1, a2, NULL,
			100, // n
			1, // lub
			scalars);
      print_scalars (stdout, "lub-poly", r, scalars);

    }

  return 0;
}
