/* bug fixing for polydisperse systems
 * lubrication routines -- atimes procedure
 * Copyright (C) 1993-2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include "stokes.h"

#include "f.h" // calc_lub_f_2b()
#include "ft.h" // calc_lub_ft_2b()
#include "fts.h" // calc_lub_fts_2b()

#include "lub.h"

#include "f-new.h"
#include "ft-new.h"
#include "fts-new.h"
#include "twobody-slip-new.h" // calc_lub_[f,ft,fts]_2b_slip_new()

#include "lub-new.h"


// copy from lub.c
/*
 * INPUT
 *  sys : struct stokes *sys.
 *        lx, ly, lz, and lubmax are used.
 * OUTPUT
 *  *imax[xyz] : range for the image-cell loop
 */
static void
set_imax_lub_periodic (struct stokes *sys,
		       int *imaxx, int *imaxy, int *imaxz)
{
  double lubmax;
  if (sys->lubmax <= 0.0)
    {
      // if lubmax2 is not define, at least, take symmetric in x, y, z.
      // set lubmax by max (lx, ly, lz)
      lubmax = sys->lx;
      if (sys->ly > lubmax) lubmax = sys->ly;
      if (sys->lz > lubmax) lubmax = sys->lz;
    }
  else
    {
      lubmax = sys->lubmax;
    }

  *imaxx = (int)(lubmax / sys->lx) + 1;
  *imaxy = (int)(lubmax / sys->ly) + 1;
  *imaxz = (int)(lubmax / sys->lz) + 1;
}


/* calculate lubrication f by u for all particles
 * for both under the periodic and non-periodic boundary conditions.
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
void
calc_lub_3f_new
(struct stokes *sys,
 const double *u,
 double *f)
{
  int i, j;
  int i3;
  int j3;

  double tmp_pos[3];

  int np = sys->np;
  /* clear f [np * 3] */
  for (i = 0; i < np * 3; ++i)
    {
      f [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  /*
	  if (list_ex_check (sys->ex_lub, i, j) == 1)
	    {
	      // j is in the exclusion list for i
	      continue;
	    }
	  */

	  j3 = j * 3;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  // mono noslip is OK
			  calc_lub_f_2b
			    (sys,
			     u + i3, u + j3,
			     sys->pos + i3, sys->pos + j3,
			     f + i3, f + j3);
			}
		      else
			{
			  // slip
			  calc_lub_f_2b_slip_new
			    (sys,
			     u + i3, u + j3,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     f + i3, f + j3);
			}
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  calc_lub_f_2b_poly_new
			    (sys,
			     u + i3, u + j3,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     f + i3, f + j3);
			}
		      else
			{
			  // slip
			  calc_lub_f_2b_slip_new
			    (sys,
			     u + i3, u + j3,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     f + i3, f + j3);
			}
		    }
		}
	    }
	  else
	    {
	      // periodic system
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  /* shift for shear */
			  if (sys->shear_mode == 1)
			    {
			      tmp_pos[0] += (double)iy * sys->shear_shift;
			    }
			  else if (sys->shear_mode == 2)
			    {
			      tmp_pos[0] += (double)iz * sys->shear_shift;
			    }

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      // mono noslip is OK
				      calc_lub_f_2b
					(sys,
					 u + i3, u + j3,
					 sys->pos + i3, tmp_pos,
					 f + i3, f + j3);
				    }
				  else
				    {
				      // slip
				      calc_lub_f_2b_slip_new
					(sys,
					 u + i3, u + j3,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 f + i3, f + j3);
				    }
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, tmp_pos,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      calc_lub_f_2b_poly_new
					(sys,
					 u + i3, u + j3,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 f + i3, f + j3);
				    }
				  else
				    {
				      // slip
				      calc_lub_f_2b_slip_new
					(sys,
					 u + i3, u + j3,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 f + i3, f + j3);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

/* calculate lubrication ft by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
void
calc_lub_3ft_new
(struct stokes * sys,
 const double * uo, double * ft)
{
  int i, j;
  int i3, i6;
  int j3, j6;

  double tmp_pos[3];

  int np = sys->np;
  /* clear ft [np * 6] */
  for (i = 0; i < np * 6; ++i)
    {
      ft [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = i; j < np; ++j)
	{
	  /*
	  if (list_ex_check (sys->ex_lub, i, j) == 1)
	    {
	      // j is in the exclusion list for i
	      continue;
	    }
	  */

	  j3 = j * 3;
	  j6 = j * 6;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  // mono noslip is OK
			  calc_lub_ft_2b
			    (sys,
			     uo + i6, uo + j6,
			     sys->pos + i3, sys->pos + j3,
			     ft + i6, ft + j6);
			}
		      else
			{
			  // slip
			  calc_lub_ft_2b_slip_new
			    (sys,
			     uo + i6, uo + j6,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     ft + i6, ft + j6);
			}
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  calc_lub_ft_2b_poly_new
			    (sys,
			     uo + i6, uo + j6,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     ft + i6, ft + j6);
			}
		      else
			{
			  // slip
			  calc_lub_ft_2b_slip_new
			    (sys,
			     uo + i6, uo + j6,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     ft + i6, ft + j6);
			}
		    }
		}
	    }
	  else
	    {
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  /* shift for shear */
			  if (sys->shear_mode == 1)
			    {
			      tmp_pos[0] += (double)iy * sys->shear_shift;
			    }
			  else if (sys->shear_mode == 2)
			    {
			      tmp_pos[0] += (double)iz * sys->shear_shift;
			    }

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      // mono noslip is OK
				      calc_lub_ft_2b
					(sys,
					 uo + i6, uo + j6,
					 sys->pos + i3, tmp_pos,
					 ft + i6, ft + j6);
				    }
				  else
				    {
				      // slip
				      calc_lub_ft_2b_slip_new
					(sys,
					 uo + i6, uo + j6,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 ft + i6, ft + j6);
				    }
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, tmp_pos,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      calc_lub_ft_2b_poly_new
					(sys,
					 uo + i6, uo + j6,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 ft + i6, ft + j6);
				    }
				  else
				    {
				      // slip
				      calc_lub_ft_2b_slip_new
					(sys,
					 uo + i6, uo + j6,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 ft + i6, ft + j6);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

/* calculate lubrication fts by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse and slip systems can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_3fts_new
(struct stokes * sys,
 const double * uoe, double * fts)
{
  int i, j;
  int i3, i11;
  int j3, j11;

  double tmp_pos[3];

  int np = sys->np;
  /* clear fts [np * 11] */
  for (i = 0; i < np * 11; ++i)
    {
      fts [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i; j < np; ++j)
	{
	  /*
	  if (list_ex_check (sys->ex_lub, i, j) == 1)
	    {
	      // j is in the exclusion list for i
	      continue;
	    }
	  */

	  j3 = j * 3;
	  j11 = j * 11;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  // mono noslip is OK
			  calc_lub_fts_2b
			    (sys,
			     uoe + i11, uoe + j11,
			     sys->pos + i3, sys->pos + j3,
			     fts + i11, fts + j11);
			}
		      else
			{
			  // slip
			  calc_lub_fts_2b_slip_new
			    (sys,
			     uoe + i11, uoe + j11,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     fts + i11, fts + j11);
			}
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      if (sys->slip == NULL)
			{
			  // noslip
			  calc_lub_fts_2b_poly_new
			    (sys,
			     uoe + i11, uoe + j11,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     fts + i11, fts + j11);
			}
		      else
			{
			  // slip
			  calc_lub_fts_2b_slip_new
			    (sys,
			     uoe + i11, uoe + j11,
			     sys->pos + i3, sys->pos + j3,
			     i, j,
			     fts + i11, fts + j11);
			}
		    }
		}
	    }
	  else
	    {
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  /* shift for shear */
			  if (sys->shear_mode == 1)
			    {
			      tmp_pos[0] += (double)iy * sys->shear_shift;
			    }
			  else if (sys->shear_mode == 2)
			    {
			      tmp_pos[0] += (double)iz * sys->shear_shift;
			    }

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      // mono noslip is OK
				      calc_lub_fts_2b
					(sys,
					 uoe + i11, uoe + j11,
					 sys->pos + i3, tmp_pos,
					 fts + i11, fts + j11);
				    }
				  else
				    {
				      // slip
				      calc_lub_fts_2b_slip_new
					(sys,
					 uoe + i11, uoe + j11,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 fts + i11, fts + j11);
				    }
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, tmp_pos,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  if (sys->slip == NULL)
				    {
				      // noslip
				      calc_lub_fts_2b_poly_new
					(sys,
					 uoe + i11, uoe + j11,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 fts + i11, fts + j11);
				    }
				  else
				    {
				      // slip
				      calc_lub_fts_2b_slip_new
					(sys,
					 uoe + i11, uoe + j11,
					 sys->pos + i3, tmp_pos,
					 i, j,
					 fts + i11, fts + j11);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
