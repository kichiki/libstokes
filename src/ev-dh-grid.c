/* excluded-volume interactions by Debye-Huckel with RYUON_grid
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-grid.c,v 1.2 2008/10/31 05:44:15 kichiki Exp $
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
#include <stdlib.h>
#include <math.h> // sqrt()
#include "memory-check.h" // CHECK_MALLOC

#include "stokes.h" // struct stokes
#include "ev-dh.h" // struce EV_DH

#include "grid.h" // struct RYUON_grid


/*
 * for non-periodic system
 * INPUT
 *  ev_dh      : struct EV_DH
 *  g          : struct RYUON_grid
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid (struct EV_DH *ev_dh,
		       struct RYUON_grid *g,
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

  // loop for pairs within the range
  int ic;
  for (ic = 0; ic < g->n; ic ++)
    {
      //fprintf (stdout, "# grid [%d] np = %d\n", ic, g->np[ic]);
      int ipp;
      for (ipp = 0; ipp < g->np[ic]; ipp ++)
	{
	  int ip = g->ip[ic][ipp];
	  if (ip >= sys->nm) continue; // fixed particle
	  if (ev_dh->nu[ip] == 0.0) continue;

	  int ix, iy, iz;
	  GRID_in_to_ixyz (g, ic, &ix, &iy, &iz);
	  int jx;
	  for (jx = ix - 1; jx <= ix + 1; jx++)
	    {
	      if (jx < 0 || jx >= g->nx) continue; // out of range
	      int jy;
	      for (jy = iy - 1; jy <= iy + 1; jy++)
		{
		  if (jy < 0 || jy >= g->ny) continue; // out of range
		  int jz;
		  for (jz = iz - 1; jz <= iz + 1; jz++)
		    {
		      if (jz < 0 || jz >= g->nz) continue; // out of range

		      int jc = GRID_ixyz_to_in (g, jx, jy, jz);
		      int jpp;
		      for (jpp = 0; jpp < g->np[jc]; jpp ++)
			{
			  int jp = g->ip[jc][jpp];
			  // skip the self part and pair of (ip > jp)
			  if (ip >= jp) continue;
			  if (ev_dh->nu[jp] == 0.0) continue;

			  // now ip, jp is the pair within the range
			  int i3 = ip * 3;
			  int j3 = jp * 3;
			  double x = sys->pos [i3+0] - sys->pos [j3+0];
			  double y = sys->pos [i3+1] - sys->pos [j3+1];
			  double z = sys->pos [i3+2] - sys->pos [j3+2];

			  double r2 = x*x + y*y + z*z;
			  if (r2 < ev_dh->r2)
			    {
			      EV_DH_set_force_ij (sys, ev_dh,
						  ip, jp,
						  r2, x, y, z,
						  f);
			    }
			}
		    }
		}
	    }
	}
    }
}

/*
 * for periodic system
 * INPUT
 *  ev_dh      : struct EV_DH
 *  g          : struct RYUON_grid
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid_periodic (struct EV_DH *ev_dh,
				struct RYUON_grid *g,
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

  // loop for pairs within the range
  int ic;
  for (ic = 0; ic < g->n; ic ++)
    {
      int ipp;
      for (ipp = 0; ipp < g->np[ic]; ipp ++)
	{
	  int ip = g->ip[ic][ipp];
	  if (ip >= sys->nm) continue; // fixed particle
	  if (ev_dh->nu[ip] == 0.0) continue;

	  int ix, iy, iz;
	  double llx = 0.0;
	  double lly = 0.0;
	  double llz = 0.0;
	  GRID_in_to_ixyz (g, ic, &ix, &iy, &iz);
	  int jx;
	  for (jx = ix - 1; jx <= ix + 1; jx++)
	    {
	      if (jx < 0)
		{
		  jx += g->nx;
		  llx -= sys->lx;
		}
	      if (jx >= g->nx)
		{
		  jx -= g->nx;
		  llx += sys->lx;
		}
	      int jy;
	      for (jy = iy - 1; jy <= iy + 1; jy++)
		{
		  if (jy < 0)
		    {
		      jy += g->ny;
		      lly -= sys->ly;
		    }
		  if (jy >= g->ny)
		    {
		      jy -= g->ny;
		      lly += sys->ly;
		    }
		  int jz;
		  for (jz = iz - 1; jz <= iz + 1; jz++)
		    {
		      if (jz < 0)
			{
			  jz += g->nz;
			  llz -= sys->lz;
			}
		      if (jz >= g->nz)
			{
			  jz -= g->nz;
			  llz += sys->lz;
			}

		      int jc = GRID_ixyz_to_in (g, jx, jy, jz);
		      int jpp;
		      for (jpp = 0; jpp < g->np[jc]; jpp ++)
			{
			  int jp = g->ip[jc][jpp];
			  // skip the self part and pair of (ip > jp)
			  if (ip >= jp) continue;
			  if (ev_dh->nu[jp] == 0.0) continue;

			  // now ip, jp is the pair within the range
			  int i3 = ip * 3;
			  int j3 = jp * 3;
			  double x = sys->pos [i3+0]
			    - (sys->pos [j3+0] + llx);
			  double y = sys->pos [i3+1]
			    - (sys->pos [j3+1] + lly);
			  double z = sys->pos [i3+2]
			    - (sys->pos [j3+2] + llz);

			  double r2 = x*x + y*y + z*z;
			  if (r2 < ev_dh->r2)
			    {
			      EV_DH_set_force_ij (sys, ev_dh,
						  ip, jp,
						  r2, x, y, z,
						  f);
			    }
			}
		    }
		}
	    }
	  
	}
    }
}
