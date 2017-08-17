/* confinement interaction to the particles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: confinement.c,v 1.1 2008/05/24 05:43:45 kichiki Exp $
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
#include <stdlib.h>
#include <math.h> // M_PI, sqrt()

#include "stokes.h" // struct stokes
#include "memory-check.h" // macro CHECK_MALLOC

#include "confinement.h"


/* initialize struct confinement
 * INPUT
 *  type : confinement type parameter (0 - 4)
 *         0 == sphere
 *               R : radius of cavity centered at (0, 0, 0)
 *         1 == sphere + hole
 *               R : radius of cavity centered at (0, 0, 0)
 *               r : radius of the hole at (0, 0, 1) direction
 *         2 == cylinder
 *               r       : radius of the cylinder
 *               x, y, z : direction vector of the cylinder
 *               the cylinder center goes through (0, 0, 0) and (x, y, z).
 *         3 == dumbbell
 *               R  : left cavity radius centered at (center1, 0, 0)
 *               R2 : right cavity radius centered at (center2, 0, 0)
 *               L  : length of the cylinder
 *               r  : cylinder radius
 *               the origin is at the center of the cylinder
 *         4 == hex2d
 *               R : cavity radius
 *               r : cylinder radius
 *               L : lattice spacing
 *         5 == porous
 *               R : particle radius
 *               L : lattice spacing in x (2R for touching case)
 *  R       : cavity radius
 *  r       : cylinder radius
 *  x, y, z : cylinder direction (for type 2, cylinder)
 *  R2      : right cavity radius (for type 3, dumbbell)
 *  L       : lattice spacing (for type 4, hex2d)
 *  Lennard-Jones parameters
 *  flag_LJ : if 0, e = hat(e) = "e" / kT
 *            if 1, e = hat(e) / (peclet * hat(r0))
 *  e  : the dimensionless number defined by
 *         e = hat(e) for flag_LJ == 0, or
 *         e = hat(e) / (peclet * hat(r0)) for flag_LJ == 1,
 *       where
 *         hat(e)  = "e" / kT
 *         hat(r0) = "r0" / length
 *       note that "e" and "r0" are dimensional numbers.
 *  r0 : the dimensionless number defined by
 *         r0 = hat(r0) = "r0" / length
 *       note that "r0" is a dimensional number.
 */
struct confinement *
CF_init (int type,
	 double R,
	 double r,
	 double x, double y, double z,
	 double R2,
	 double L,
	 int flag_LJ,
	 double e,
	 double r0)
{
  struct confinement *cf
    = (struct confinement *)malloc (sizeof (struct confinement));
  CHECK_MALLOC (cf, "CF_init");

  // some parameter check
  if ((type != 0 && type != 2) && R < r)
    {
      fprintf (stderr, "CF_init : r is larger than R. it is set by 0 now.\n");
      r = 0.0;
    }
  if (type == 3 && R2 < r)
    {
      fprintf (stderr, "CF_init : r is larger than R2. it is set by 0 now.\n");
      r = 0.0;
    }

  cf->type = type;
  cf->R  = R;
  cf->r  = r;
  cf->x  = x;
  cf->y  = y;
  cf->z  = z;
  cf->R2 = R2;
  cf->L  = L;

  cf->flag_LJ = flag_LJ;
  cf->e  = e;
  cf->r0 = r0;

  // set derived parameters
  if (R > 0.0 && R > r)
    {
      cf->theta = asin (r / R);
    }
  else
    {
      cf->theta = 0.0;
    }

  // for dumbbell
  if (type == 3)
    {
      if (R2 > 0.0 && R2 > r)
	{
	  cf->theta2 = asin (r / R2);
	}
      else
	{
	  cf->theta2 = 0.0;
	}

      double L2 = 0.5 * L;
      // x component of the center of left cavity
      cf->center1 = -(L2 + R * cos (cf->theta));

      // x component of the center of right cavity
      cf->center2 = L2 + R2 * cos (cf->theta2);
    }
  else
    {
      cf->theta2 = 0.0;
      cf->center1 = 0.0;
      cf->center2 = 0.0;
    }

  // for hex2d
  if (type == 4)
    {
      cf->Lx = 0.5 * L;
      cf->Ly = sqrt (3.0) * cf->Lx;
      cf->Lz = 0.0;
    }
  // for porous
  else if (type == 5)
    {
      cf->Lx = L;
      cf->Ly = L * sqrt (3.0);
      cf->Lz = L * sqrt (6.0) * 2.0 / 3.0;
    }
  else
    {
      cf->Lx = 0.0;
      cf->Ly = 0.0;
      cf->Lz = 0.0;
    }

  return (cf);
}

void
CF_free (struct confinement *cf)
{
  if (cf != NULL) free (cf);
}

/* set LJ parameters for run
 * INPUT
 *  cf     : struct confinement
 *  peclet : peclet number
 * OUTPUT
 *  cf->e  : defined as 
 *             cf->e = hat(e)/(peclet * hat(r0))
 *           where
 *             hat(e)  = "e" / kT
 *             hat(r0) = "r0" / length
 *             "e" and "r0" are dimensional numbers
 *           therefore, cf->flag-LJ is set by 1.
 */
void
CF_set (struct confinement *cf,
	double peclet)
{
  if (cf->flag_LJ == 0)
    {
      cf->flag_LJ = 1;
      cf->e = cf->e / peclet / cf->r0;
    }
}

/*
 * INPUT
 *  x : depth of overlap
 *      therefore, the LJ distance is given by 
 *        r = r0 - x.
 *      x have to be smaller than r0.
 *      zero is returned if x < 0 (no overlap).
 */
static double
CF_fr (struct confinement *cf,
       double x)
{
  if (x < 0.0) return (0.0);

  // NOTE: cf->r0 is converted to hat(r0) = r0 / length
  double r = cf->r0 - x;
  if (r < 0.0)
    {
      fprintf (stderr, "CF_fr: too much overlap."
	       " the force value is adjusted.\n");
      r = 1.0e-12; // give some small value
    }

  double ri = cf->r0 / r;
  double ri6 = pow (ri, 6.0);
  double ri7 = ri6 * ri;

  // NOTE: cf->e is converted to hat(e)/(peclet * hat(r0))
  double fr = 12.0 * cf->e * ri7 * (ri6 - 1.0);

  return (fr);
}


/*
 * INPUT
 *  cf         : struct confinement
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_sphere_calc_force (struct confinement *cf,
		      struct stokes *sys,
		      double *f,
		      int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      double dx = x;
      double dy = y;
      double dz = z;
      double r2 = dx * dx + dy * dy + dz * dz;

      double Ra = cf->R - a;
      double Ra2 = Ra * Ra;
      if (r2 > Ra2)
	{
	  // particles are touching the cavity wall
	  double r = sqrt (r2);
	  dx = dx / r;
	  dy = dy / r;
	  dz = dz / r;

	  double fr = CF_fr (cf, r - Ra);
	  // where (r - Ra) > 0
	  f[ix] += -fr * dx;
	  f[iy] += -fr * dy;
	  f[iz] += -fr * dz;
	}
    }
}

/*
 * INPUT
 *  CF : struct confinement
 *  sys            : struct stokes (only nm and pos are used)
 *  f [nm * 3]     : force is assigned only for the mobile particles
 *  flag_add       : if 0 is given, zero-clear and set the force
 *                   otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_sphere_hole_calc_force (struct confinement *cf,
			   struct stokes *sys,
			   double *f,
			   int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      double dx = x;
      double dy = y;
      double dz = z;
      double r2 = dx * dx + dy * dy + dz * dz;

      double Ram = cf->R - a;
      double Ram2 = Ram * Ram;
      double Rap = cf->R + a;
      double Rap2 = Rap * Rap;
      if (r2 > Ram2 && r2 < Rap2)
	{
	  // check the contact with the hole edge
	  // first, check the distance from the hole center
	  double hx = x;
	  double hy = y;
	  double hz = z - cf->R * cos (cf->theta);
	  double h2 = hx * hx + hy * hy + hz * hz;
	  double ra = cf->r - a;
	  double ra2 = ra * ra;
	  if (h2 > ra2)
	    {
	      // obtain the distance from the nearest edge point
	      double theta = atan2 (dy, dx);
	      double ex = x - cf->r * cos (theta);
	      double ey = y - cf->r * sin (theta);
	      double ez = hz; // same above
	      double e2 = ex * ex + ey * ey + ez * ez;

	      double a2 = a * a;
	      if (e2 < a2)
		{
		  // touching to the edge
		  double r = sqrt (e2);
		  dx = ex / r;
		  dy = ey / r;
		  dz = ez / r;

		  double fr = CF_fr (cf, a - r);
		  // where (a - r) > 0
		  f[ix] += fr * dx;
		  f[iy] += fr * dy;
		  f[iz] += fr * dz;
		  // vector e is pointing "from" the wall "to" the particle
		}
	      else
		{
		  // particles are touching the cavity wall
		  double r = sqrt (r2);
		  dx = dx / r;
		  dy = dy / r;
		  dz = dz / r;

		  if (r < cf->R)
		    {
		      // inside the cavity
		      double fr = CF_fr (cf, r - Ram);
		      // where (r - Ram) > 0
		      f[ix] += -fr * dx;
		      f[iy] += -fr * dy;
		      f[iz] += -fr * dz;
		    }
		  else
		    {
		      // outside the cavity
		      double fr = CF_fr (cf, Rap - r);
		      // where (Rap - r) > 0
		      f[ix] += fr * dx;
		      f[iy] += fr * dy;
		      f[iz] += fr * dz;
		    }
		}
	    }
	}
    }
}

/*
 * INPUT
 *  (x,y,z)    : position of the particle
 *  a          : radius of the particle
 *  (x1,y1,z1) : one point on the center line of the cylinder
 *  (x2,y2,z2) : another point on the center line of the cylinder
 *  r_cyl      : radius of the cylinder
 * OUTPUT
 *  f[3]       : force is added
 */
static void
CF_cylinder (struct confinement *cf,
	     double x, double y, double z, double a, 
	     double x1, double y1, double z1,
	     double x2, double y2, double z2, double r_cyl,
	     double *f)
{
  // cylinder vector
  double rx = x2 - x1;
  double ry = y2 - y1;
  double rz = z2 - z1;
  double r  = sqrt (rx * rx + ry * ry + rz * rz);
  rx /= r;
  ry /= r;
  rz /= r;
  // now (rx,ry,rz) is the unit vector

  // vector from (x1,y1,z1) to the particle
  double px = x - x1;
  double py = y - y1;
  double pz = z - z1;

  // inner product
  double pr
    = px * rx
    + py * ry
    + pz * rz;

  // tangential vector
  double tx = pr * rx;
  double ty = pr * ry;
  double tz = pr * rz;

  // normal vector
  double nx = px - tx;
  double ny = py - ty;
  double nz = pz - tz;
  double n2
    = nx * nx
    + ny * ny
    + nz * nz;
  // n2 is the squared distance from the center line of the cylinder
  /*
  fprintf (stderr, "# check for (%f,%f,%f), d = %f, x0 = (%f,%f,%f)\n",
	   x, y, z,
	   sqrt(n2),
	   tx, ty, tz);
  */

  double ra = r_cyl - a;
  double ra2 = ra * ra;
  if (n2 <= ra2) return;
  //fprintf (stderr, "# cylinder contact n2 > ra2 : %e, %e\n", n2, ra2);

  double n  = sqrt(n2);
  double dx = nx / n;
  double dy = ny / n;
  double dz = nz / n;
  // now (dx,dy,dz) is the unit vector pointing to the particle

  double fr = CF_fr (cf, n - ra);
  // where (n - ra) > 0
  f[0] += -fr * dx;
  f[1] += -fr * dy;
  f[2] += -fr * dz;
}

/*
 * INPUT
 *  cf_cylinder : struct CF_cylinder
 *  sys         : struct stokes (only nm and pos are used)
 *  f [nm * 3]  : force is assigned only for the mobile particles
 *  flag_add    : if 0 is given, zero-clear and set the force
 *                otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_cylinder_calc_force (struct confinement *cf,
			struct stokes *sys,
			double *f,
			int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      CF_cylinder (cf,
		   x, y, z, a,
		   0.0, 0.0, 0.0,
		   cf->x, cf->y, cf->z,
		   cf->r,
		   f + ix);
    }
}

/*
 * INPUT
 *  cf_hex2d   : struct CF_hex2d
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_dumbbell_calc_force (struct confinement *cf,
			struct stokes *sys,
			double *f,
			int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      if (x < -0.5 * cf->L)
	{
	  // left cavity at (-L/2-R1, 0)
	  double dx = x - cf->center1;
	  double dy = y;
	  double dz = z;
	  double r2 = dx * dx + dy * dy + dz * dz;

	  // cf->R is the radius of the left cavity
	  double Ra = cf->R - a;
	  double Ra2 = Ra * Ra;
	  if (r2 > Ra2)
	    {
	      double theta = atan2 (dy, dx);
	      // cf->theta is for the left cavity
	      if (fabs (theta) < cf->theta)
		{
		  // obtain the distance from the nearest edge point
		  double alpha = atan2 (dz, dy);
		  double ex = x + 0.5 * cf->L; // -L/2 is the position
		  double ey = y - cf->r * cos (alpha);
		  double ez = z - cf->r * sin (alpha);
		  double e2 = ex * ex + ey * ey + ez * ez;

		  double a2 = a * a;
		  if (e2 < a2)
		    {
		      // touching to the edge
		      double r = sqrt (e2);
		      dx = ex / r;
		      dy = ey / r;
		      dz = ez / r;

		      double fr = CF_fr (cf, a - r);
		      // where (a - r) > 0
		      f[ix] += fr * dx;
		      f[iy] += fr * dy;
		      f[iz] += fr * dz;
		      // vector e is pointing "from" the wall "to" the particle
		    }
		}
	      else
		{
		  // touching cavity wall
		  double r = sqrt (r2);
		  dx = dx / r;
		  dy = dy / r;
		  dz = dz / r;

		  double fr = CF_fr (cf, r - Ra);
		  // where (r - Ra) > 0
		  f[ix] += -fr * dx;
		  f[iy] += -fr * dy;
		  f[iz] += -fr * dz;
		}
	    }
	}
      else if (x > 0.5 * cf->L)
	{
	  // right cavity at (+L/2+R2, 0)
	  double dx = x - cf->center2;
	  double dy = y;
	  double dz = z;
	  double r2 = dx * dx + dy * dy + dz * dz;

	  double Ra = cf->R2 - a;
	  double Ra2 = Ra * Ra;
	  if (r2 > Ra2)
	    {
	      double theta = atan2 (dy, -dx);
	      // therefore, the cylinder is around theta = 0
	      if (fabs (theta) < cf->theta2)
		{
		  // obtain the distance from the nearest edge point
		  double alpha = atan2 (dz, dy);
		  double ex = x - 0.5 * cf->L;
		  double ey = y - cf->r * cos (alpha);
		  double ez = z - cf->r * sin (alpha);
		  double e2 = ex * ex + ey * ey + ez * ez;

		  double a2 = a * a;
		  if (e2 < a2)
		    {
		      // touching to the edge
		      double r = sqrt (e2);
		      dx = ex / r;
		      dy = ey / r;
		      dz = ez / r;

		      double fr = CF_fr (cf, a - r);
		      // where (a - r) > 0
		      f[ix] += fr * dx;
		      f[iy] += fr * dy;
		      f[iz] += fr * dz;
		      // vector e is pointing "from" the wall "to" the particle
		    }
		}
	      else
		{
		  // touching cavity wall
		  double r = sqrt (r2);
		  dx = dx / r;
		  dy = dy / r;
		  dz = dz / r;

		  double fr = CF_fr (cf, r - Ra);
		  // where (r - Ra) > 0
		  f[ix] += -fr * dx;
		  f[iy] += -fr * dy;
		  f[iz] += -fr * dz;
		}
	    }
	}
      else
	{
	  // cylinder between (-L/2, 0) to (+L/2, 0)
	  CF_cylinder (cf,
		       x, y, z, a,
		       0.0, 0.0, 0.0,
		       cf->L, 0.0, 0.0,
		       cf->r,
		       f + ix);
	}
    }
}

/* consider a cavity at (0, 0, 0)
 * in the lower left half of the region (0, 0) to (Lx, Ly)
 * INPUT
 *  x, y, z : particle position (should be in the region)
 *  a       : radius of the particle
 * OUTPUT
 *  f[3]    : added (not zero-cleared)
 */
static void
CF_hex2d_1cavity_calc_force (struct confinement *cf,
			     double x, double y, double z,
			     double a,
			     double *f)
{
  double r2 = x * x + y * y + z * z;
  double Ra = cf->R - a;
  double Ra2 = Ra * Ra;
  double theta60 = M_PI / 3.0; // = 60 degree
  
  // consider in cavity at (0,0)
  if (r2 > Ra2)
    {
      //fprintf (stderr, "# r2 > Ra2 : %e, %e\n", r2, Ra2);
      // particles are touching the cavity wall
      double theta = atan2 (y, x);
      // theta is in the range [0, pi/2]

      double E = cf->R * cos (cf->theta);
      if (fabs (theta) < cf->theta)
	{
	  /*
	  fprintf (stderr, "# |theta| < cf->theta : %e, %e\n",
		   theta, cf->theta);
	  */
	  if (x < E)
	    {
	      //fprintf (stderr, "# x < E : %e, %e\n", x, E);
	      // obtain the distance from the nearest edge point
	      double alpha = atan2 (z, y);
	      double ex = x - E;
	      double ey = y - cf->r * cos (alpha);
	      double ez = z - cf->r * sin (alpha);
	      double e2 = ex * ex + ey * ey + ez * ez;

	      double a2 = a * a;
	      if (e2 < a2)
		{
		  fprintf (stderr, "#1 e2 < a2 : %e, %e\n", e2, a2);
		  // touching to the edge
		  double r = sqrt (e2);
		  double dx = ex / r;
		  double dy = ey / r;
		  double dz = ez / r;

		  double fr = CF_fr (cf, a - r);
		  // where (a - r) > 0
		  f[0] += fr * dx;
		  f[1] += fr * dy;
		  f[2] += fr * dz;
		  // vector e is pointing "from" the wall "to" the particle
		}
	    }
	  else
	    {
	      //fprintf (stderr, "# in Cylinder at Lx,0\n");
	      // in the cylinder at (Lx, 0)
	      CF_cylinder (cf,
			   x, y, z, a,
			   0.0, 0.0, 0.0,
			   cf->Lx, 0.0, 0.0,
			   cf->r,
			   f);
	    }
	}
      else if (fabs (theta - theta60) < cf->theta)
	{
	  // transform to (x60,y60) coordinate
	  double x60 = 0.5 * (x + sqrt (3.0) * y);
	  if (x60 < E)
	    {
	      //fprintf (stderr, "# x60 < E : %e, %e\n", x60, E);
	      double y60 = 0.5 * (-sqrt (3.0) * x + y);
	      // obtain the distance from the nearest edge point
	      double alpha = atan2 (z, y60);
	      double ex = x60 - E;
	      double ey = y60 - cf->r * cos (alpha);
	      double ez = z   - cf->r * sin (alpha);
	      double e2 = ex * ex + ey * ey + ez * ez;

	      double a2 = a * a;
	      if (e2 < a2)
		{
		  fprintf (stderr, "#2 e2 < a2 : %e, %e\n", e2, a2);
		  // touching to the edge
		  double r = sqrt (e2);
		  // back to (x,y) coordinate
		  double dx = 0.5 * (ex - sqrt (3.0) * ey) / r;
		  double dy = 0.5 * (sqrt (3.0) * ex + ey) / r;
		  double dz = ez / r;

		  double fr = CF_fr (cf, a - r);
		  // where (a - r) > 0
		  f[0] += fr * dx;
		  f[1] += fr * dy;
		  f[2] += fr * dz;
		  // vector e is pointing "from" the wall "to" the particle
		}
	    }
	  else
	    {
	      //fprintf (stderr, "# in Cylinder at Lx/2,Ly/2\n");
	      // in the cylinder at (Lx/2, Ly/2)
	      CF_cylinder (cf,
			   x, y, z, a,
			   0.0, 0.0, 0.0,
			   cf->Lx, cf->Ly, 0.0,
			   cf->r,
			   f);
	    }
	}
      else
	{
	  fprintf (stderr, "# touching cavity wall r2 > Ra2 : %e, %e\n",
		   r2, Ra2);
	  double r = sqrt (r2);
	  double dx = x / r;
	  double dy = y / r;
	  double dz = z / r;

	  // touching cavity wall
	  double fr = CF_fr (cf, r - Ra);
	  // where (r - Ra) > 0
	  f[0] += -fr * dx;
	  f[1] += -fr * dy;
	  f[2] += -fr * dz;
	}
    }
}

/*
 * INPUT
 *  cf_hex2d   : struct CF_hex2d
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_hex2d_calc_force (struct confinement *cf,
		     struct stokes *sys,
		     double *f,
		     int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double f_local[3] = {0.0, 0.0, 0.0};

      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      int ilx = (int)floor(x / cf->Lx);
      int ily = (int)floor(y / cf->Ly);
      x -= (double)ilx * cf->Lx;
      y -= (double)ily * cf->Ly;
      /* now x is in range [0, Lx)
       * and y is in range [0, Ly)
       */

      if ((ilx % 2) != (ily % 2))
	{
	  // flip x -> -x and shift in [0, Lx)
	  x = -x; // now x is in [-Lx, 0)
	  x += cf->Lx; // now x is in [0, Lx)
	}
      // check
      if (x < 0.0 || x >= cf->Lx ||
	  y < 0.0 || y >= cf->Ly)
	{
	  fprintf (stderr, "CF_hex2d_calc_force"
		   " : out of range (x,y) = (%f,%f) / (Lx,Ly) = (%f,%f)\n",
		   x, y, cf->Lx, cf->Ly);
	}
      //fprintf (stderr, "# shifted => x,y,z = %f, %f, %f\n", x, y, z);

      // consider (0,0)-block

      // determine which cavity the particle is
      // sphere at (0,0)
      double dx_0 = x;
      double dy_0 = y;
      double dz_0 = z;
      double r2_0 = dx_0 * dx_0 + dy_0 * dy_0 + dz_0 * dz_0;
      // sphere at (Lx,Ly)
      double dx_1 = x - cf->Lx;
      double dy_1 = y - cf->Ly;
      double dz_1 = z;
      double r2_1 = dx_1 * dx_1 + dy_1 * dy_1 + dz_1 * dz_1;
      if (r2_0 < r2_1)
	{
	  // consider in cavity at (0,0)
	  CF_hex2d_1cavity_calc_force (cf,
				       x, y, z, a,
				       f_local);
	}
      else
	{
	  // consider in cavity at (Lx,Ly)
	  // convert (x,y) => (x',y') = (Lx-x, Ly-y)
	  // so that the center is at (x'=0, y'=0)
	  double xx = cf->Lx - x;
	  double yy = cf->Ly - y;
	  CF_hex2d_1cavity_calc_force (cf,
				       xx, yy, z, a,
				       f_local);
	  // convert (x',y') => (x,y), that is, f = -ff
	  f_local[0] = -f_local[0];
	  f_local[1] = -f_local[1];
	}

      if ((ilx % 2) != (ily % 2))
	{
	  // flip x component of the force
	  f_local[0] = -f_local[0];
	}

      // for check
      if (fabs (f_local[0]) > 1.0e-12 ||
	  fabs (f_local[1]) > 1.0e-12 ||
	  fabs (f_local[2]) > 1.0e-12)
	{
	  fprintf (stdout, "# contact! f = (%f, %f, %f)\n",
		   f_local[0], f_local[1], f_local[2]);
	}
      f[ix] += f_local[0];
      f[iy] += f_local[1];
      f[iz] += f_local[2];
    }
}


/*
 * INPUT
 *  cf         : struct confinement
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_porous_calc_force (struct confinement *cf,
		      struct stokes *sys,
		      double *f,
		      int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->nm; i ++)
    {
      double a; // bead radius
      if (sys->a == NULL) a = 1.0;
      else                a = sys->a[i];

      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double x = sys->pos[ix];
      double y = sys->pos[iy];
      double z = sys->pos[iz];

      int ilx = (int)floor(x / cf->Lx); // Lx = 2R
      int ily = (int)floor(y / cf->Ly); // Ly = 2R (3^1/2) 
      int ilz = (int)floor(z / cf->Lz); // Lz = 2R (2/3) (6^1/2)
      x -= (double)ilx * cf->Lx;
      y -= (double)ily * cf->Ly;
      z -= (double)ilz * cf->Lz;

      double Ra = cf->R - a;
      double Ra2 = Ra * Ra;

      double fr;

      double dx;
      double dy;
      double dz;
      double r2;
      double r;

      // (0,0,0)
      dx = x - 0.0;
      dy = y - 0.0;
      dz = z - 0.0;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx,0,0)
      dx = x - cf->Lx;
      dy = y - 0.0;
      dz = z - 0.0;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx/2,Ly/2,0)
      dx = x - 0.5 * cf->Lx;
      dy = y - 0.5 * cf->Ly;
      dz = z - 0.0;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (0,Ly,0)
      dx = x - 0.0;
      dy = y - cf->Ly;
      dz = z - 0.0;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx,Ly,0)
      dx = x - cf->Lx;
      dy = y - cf->Ly;
      dz = z - 0.0;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx/2,Ly/6,Lz/2)
      dx = x - 0.5 * cf->Lx;
      dy = y - cf->Ly / 6.0;
      dz = z - 0.5 * cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (0,Ly2/3,Lz/2)
      dx = x - 0.0;
      dy = y - cf->Ly * 2.0 / 3.0;
      dz = z - 0.5 * cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx,Ly2/3,Lz/2)
      dx = x - cf->Lx;
      dy = y - cf->Ly * 2.0 / 3.0;
      dz = z - 0.5 * cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx/2,Ly7/6,Lz/2)
      dx = x - 0.5 * cf->Lx;
      dy = y - cf->Ly * 7.0 / 6.0;
      dz = z - 0.5 * cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (0,0,Lz)
      dx = x - 0.0;
      dy = y - 0.0;
      dz = z - cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx,0,Lz)
      dx = x - cf->Lx;
      dy = y - 0.0;
      dz = z - cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx/2,Ly/2,Lz)
      dx = x - 0.5 * cf->Lx;
      dy = y - 0.5 * cf->Ly;
      dz = z - cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (0,Ly,Lz)
      dx = x - 0.0;
      dy = y - cf->Ly;
      dz = z - cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}

      // (Lx,Ly,Lz)
      dx = x - cf->Lx;
      dy = y - cf->Ly;
      dz = z - cf->Lz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ra2)
	{
	  // contact
	  r = sqrt (r2);
	  dx /= r;
	  dy /= r;
	  dz /= r;
	  fr = CF_fr (cf, Ra - r);
	  // where (Ra - r) > 0
	  f[ix] += fr * dx;
	  f[iy] += fr * dy;
	  f[iz] += fr * dz;
	  // vector d is pointing "from" the center "to" the particle
	}
    }
}


/* wrapper for confinement forces
 * INPUT
 *  c          : struct confinement
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
CF_calc_force (struct confinement *cf,
	       struct stokes *sys,
	       double *f,
	       int flag_add)
{
  switch (cf->type)
    {
    case 0: // sphere
      CF_sphere_calc_force (cf, sys, f, flag_add);
      break;
    case 1: // sphere + hole
      CF_sphere_hole_calc_force (cf, sys, f, flag_add);
      break;
    case 2: // cylinder
      CF_cylinder_calc_force (cf, sys, f, flag_add);
      break;
    case 3: // dumbbell
      CF_dumbbell_calc_force (cf, sys, f, flag_add);
      break;
    case 4: // hex2d
      CF_hex2d_calc_force (cf, sys, f, flag_add);
      break;
    case 5: // porous
      CF_porous_calc_force (cf, sys, f, flag_add);
      break;
    default:
      fprintf (stderr, "CF_calc_force : invalid confinement type %d\n",
	       cf->type);
      break;
    }
}
