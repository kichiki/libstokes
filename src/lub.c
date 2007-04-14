/* lubrication routines -- atimes procedure
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: lub.c,v 5.5 2007/04/14 00:33:28 kichiki Exp $
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

#include "stokes.h"
#include "f.h"
#include "ft.h"
#include "fts.h"

#include "lub.h"


/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 *  lubmax2 : square of the max distance (0 means no limit)
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub (const double * x1, const double * x2, double lubmax2)
{
  double x, y, z;
  double r2;


  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  r2 = x * x
    + y * y
    + z * z;

  if (r2 == 0.0)
    {
      // self part
      return 1; // false
    }
  else if (lubmax2 <= 0.0)
    {
      // no limit
      return 0; // true
    }
  else if (r2 < lubmax2)
    {
      // within the limit
      return 0; // true
    }
  else
    {
      // out of the limit
      return 1; // false
    }
}

/* calculate lubrication f by u for all particles
 * for both under the periodic and non-periodic boundary conditions.
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
void
calc_lub_3f (struct stokes *sys,
	     const double *u,
	     double *f)
{
  int i, j, k;
  int i3;
  int j3;

  double tmp_pos[3];

  int np = sys->np;
  /* clear f [np * 3] */
  for (i = 0; i < np * 3; ++i)
    {
      f [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (cond_lub (sys->pos + i3, sys->pos + j3,
			    sys->lubmax2) == 0)
		{
		  if (sys->a == NULL)
		    {
		      // monodisperse
		      calc_lub_f_2b (sys,
				     u + i3, u + j3,
				     sys->pos + i3, sys->pos + j3,
				     f + i3, f + j3);
		    }
		  else
		    {
		      // polydisperse
		      calc_lub_f_2b_poly (sys,
					  u + i3, u + j3,
					  sys->pos + i3, sys->pos + j3,
					  sys->a[i], sys->a[j],
					  f + i3, f + j3);
		    }
		}
	    }
	  else
	    {
	      /* all image cells */
	      for (k = 0; k < 27; ++k)
		{
		  tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
		  tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
		  tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
		  if (cond_lub (sys->pos + i3, tmp_pos,
				sys->lubmax2) == 0)
		    {
		      if (sys->a == NULL)
			{
			  // monodisperse
			  calc_lub_f_2b (sys,
					 u + i3, u + j3,
					 sys->pos + i3, tmp_pos,
					 f + i3, f + j3);
			}
		      else
			{
			  // polydisperse
			  calc_lub_f_2b_poly (sys,
					      u + i3, u + j3,
					      sys->pos + i3, tmp_pos,
					      sys->a[i], sys->a[j],
					      f + i3, f + j3);
			}
		    }
		}
	    }
	}
    }
}

/* calculate lubrication ft by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
void
calc_lub_3ft (struct stokes * sys,
	      const double * uo, double * ft)
{
  int i, j, k;
  int i3, i6;
  int j3, j6;

  double tmp_pos[3];

  int np = sys->np;
  /* clear ft [np * 6] */
  for (i = 0; i < np * 6; ++i)
    {
      ft [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j6 = j * 6;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (cond_lub (sys->pos + i3, sys->pos + j3,
			    sys->lubmax2) == 0)
		{
		  if (sys->a == NULL)
		    {
		      // monodisperse
		      calc_lub_ft_2b (sys,
				      uo + i6, uo + j6,
				      sys->pos + i3, sys->pos + j3,
				      ft + i6, ft + j6);
		    }
		  else
		    {
		      // polydisperse
		      calc_lub_ft_2b_poly (sys,
					   uo + i6, uo + j6,
					   sys->pos + i3, sys->pos + j3,
					   sys->a[i], sys->a[j],
					   ft + i6, ft + j6);
		    }
		}
	    }
	  else
	    {
	      /* all image cells */
	      for (k = 0; k < 27; ++k)
		{
		  tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
		  tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
		  tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
		  if (cond_lub (sys->pos + i3, tmp_pos,
				sys->lubmax2) == 0)
		    {
		      if (sys->a == NULL)
			{
			  // monodisperse
			  calc_lub_ft_2b (sys,
					  uo + i6, uo + j6,
					  sys->pos + i3, tmp_pos,
					  ft + i6, ft + j6);
			}
		      else
			{
			  // polydisperse
			  calc_lub_ft_2b_poly (sys,
					       uo + i6, uo + j6,
					       sys->pos + i3, tmp_pos,
					       sys->a[i], sys->a[j],
					       ft + i6, ft + j6);
			}
		    }
		}
	    }
	}
    }
}

/* calculate lubrication fts by uoe for all particles
 * for both under the periodic and non-periodic boundary conditions
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_3fts (struct stokes * sys,
	       const double * uoe, double * fts)
{
  int i, j, k;
  int i3, i11;
  int j3, j11;

  double tmp_pos[3];

  int np = sys->np;
  /* clear fts [np * 11] */
  for (i = 0; i < np * 11; ++i)
    {
      fts [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j11 = j * 11;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (cond_lub (sys->pos + i3, sys->pos + j3,
			    sys->lubmax2) == 0)
		{
		  if (sys->a == NULL)
		    {
		      // monodisperse
		      calc_lub_fts_2b (sys,
				       uoe + i11, uoe + j11,
				       sys->pos + i3, sys->pos + j3,
				       fts + i11, fts + j11);
		    }
		  else
		    {
		      // polydisperse
		      calc_lub_fts_2b_poly (sys,
					    uoe + i11, uoe + j11,
					    sys->pos + i3, sys->pos + j3,
					    sys->a[i], sys->a[j],
					    fts + i11, fts + j11);
		    }
		}
	    }
	  else
	    {
	      /* all image cells */
	      for (k = 0; k < 27; ++k)
		{
		  tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
		  tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
		  tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
		  if (cond_lub (sys->pos + i3, tmp_pos,
				sys->lubmax2) == 0)
		    {
		      if (sys->a == NULL)
			{
			  // monodisperse
			  calc_lub_fts_2b (sys,
					   uoe + i11, uoe + j11,
					   sys->pos + i3, tmp_pos,
					   fts + i11, fts + j11);
			}
		      else
			{
			  // polydisperse
			  calc_lub_fts_2b_poly (sys,
						uoe + i11, uoe + j11,
						sys->pos + i3, tmp_pos,
						sys->a[i], sys->a[j],
						fts + i11, fts + j11);
			}
		    }
		}
	    }
	}
    }
}
