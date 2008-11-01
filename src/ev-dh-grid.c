/* excluded-volume interactions by Debye-Huckel with RYUON_grid
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-grid.c,v 1.3 2008/11/01 05:45:13 kichiki Exp $
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
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid (struct EV_DH *ev_dh,
		       struct stokes *sys,
		       double *f,
		       int flag_add)
{
  if (ev_dh->flag_grid == 0 ||
      ev_dh->grid == NULL)
    {
      fprintf (stderr, "# EV_DH_calc_force_grid():"
	       " flag_grid == 0 or grid == NULL\n");
      EV_DH_calc_force (ev_dh, sys, f, flag_add);
      return;
    }

  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < sys->np; i ++)
    {
      if (ev_dh->nu[i] == 0.0) continue;
      int i3 = i * 3;

      int k;
      for (k = 0; k < ev_dh->grid->nnn[i]; k ++)
	{
	  int j = ev_dh->grid->nnp[i][k];
	  if (ev_dh->nu[j] == 0.0) continue;

	  int j3 = j * 3;
	  double x = sys->pos [i3  ] - sys->pos [j3  ];
	  double y = sys->pos [i3+1] - sys->pos [j3+1];
	  double z = sys->pos [i3+2] - sys->pos [j3+2];

	  double r2 = x*x + y*y + z*z;
	  EV_DH_set_force_ij (sys, ev_dh, i, j, r2, x, y, z,
			      f);
	}
    }
}

/*
 * for periodic system
 * INPUT
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid_periodic (struct EV_DH *ev_dh,
				struct stokes *sys,
				double *f,
				int flag_add)
{
  if (ev_dh->flag_grid == 0 ||
      ev_dh->grid == NULL)
    {
      fprintf (stderr, "# EV_DH_calc_force_grid():"
	       " flag_grid == 0 or grid == NULL\n");
      EV_DH_calc_force (ev_dh, sys, f, flag_add);
      return;
    }
  fprintf (stderr, "# EV_DH_calc_force_grid_periodic():"
	   " grid version for periodic is incomplete.\n");

  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  /**
   * some remarks (2008/10/31)
   * 1) shear-shift is not taken into account yet.
   * 2) nearest neighbor list version is not implemented yet.
   */

  // loop for pairs within the range
  int ic;
  for (ic = 0; ic < ev_dh->grid->n; ic ++)
    {
      int ipp;
      for (ipp = 0; ipp < ev_dh->grid->np[ic]; ipp ++)
	{
	  int ip = ev_dh->grid->ip[ic][ipp];
	  if (ip >= sys->nm) continue; // fixed particle
	  if (ev_dh->nu[ip] == 0.0) continue;

	  int ix, iy, iz;
	  double llx = 0.0;
	  double lly = 0.0;
	  double llz = 0.0;
	  GRID_in_to_ixyz (ev_dh->grid, ic, &ix, &iy, &iz);
	  int jx;
	  for (jx = ix - 1; jx <= ix + 1; jx++)
	    {
	      if (jx < 0)
		{
		  jx += ev_dh->grid->nx;
		  llx -= sys->lx;
		}
	      if (jx >= ev_dh->grid->nx)
		{
		  jx -= ev_dh->grid->nx;
		  llx += sys->lx;
		}
	      int jy;
	      for (jy = iy - 1; jy <= iy + 1; jy++)
		{
		  if (jy < 0)
		    {
		      jy += ev_dh->grid->ny;
		      lly -= sys->ly;
		    }
		  if (jy >= ev_dh->grid->ny)
		    {
		      jy -= ev_dh->grid->ny;
		      lly += sys->ly;
		    }
		  int jz;
		  for (jz = iz - 1; jz <= iz + 1; jz++)
		    {
		      if (jz < 0)
			{
			  jz += ev_dh->grid->nz;
			  llz -= sys->lz;
			}
		      if (jz >= ev_dh->grid->nz)
			{
			  jz -= ev_dh->grid->nz;
			  llz += sys->lz;
			}

		      int jc = GRID_ixyz_to_in (ev_dh->grid, jx, jy, jz);
		      int jpp;
		      for (jpp = 0; jpp < ev_dh->grid->np[jc]; jpp ++)
			{
			  int jp = ev_dh->grid->ip[jc][jpp];
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
