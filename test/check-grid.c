/* test for grid.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-grid.c,v 1.3 2008/11/02 06:16:14 kichiki Exp $
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
#include <math.h>
#include <stdio.h> /* printf() fprintf() */
#include <stdlib.h> /* exit() */
#include <string.h> /* strcmp() */

#include "memory-check.h" // CHECK_MALLOC
#include "check.h" // compare_max

#include "grid.h"

/*
 * INPUT
 */
int
check_GRID_ixyz_to_in_to_ixyz (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_GRID_ixyz_to_in_to_ixyz: start\n");
    }

  int check = 0;
  double max = 0.0;


  double x0 = 0.0;
  double x1 = 10.0;
  double n  = 10;
  struct RYUON_grid *g = GRID_init (x0, x1,
				    x0, x1,
				    x0, x1,
				    n, n, n);
  CHECK_MALLOC (g, "GRID_ixyz_to_in_to_ixyz");

  char label [80];
  int ixx, iyy, izz;
  int ix, iy, iz;
  for (ix = 0; ix < n; ix ++)
    {
      for (iy = 0; iy < n; iy ++)
	{
	  for (iz = 0; iz < n; iz ++)
	    {
	      int in = GRID_ixyz_to_in (g, ix, iy, iz);
	      GRID_in_to_ixyz (g, in,
			       &ixx, &iyy, &izz);
	      sprintf (label, " ix for in=%d", in);
	      check += compare_max ((double)ix, (double)ixx,
				    label, verbose, tiny, &max);
	      sprintf (label, " iy for in=%d", in);
	      check += compare_max ((double)iy, (double)iyy,
				    label, verbose, tiny, &max);
	      sprintf (label, " iz for in=%d", in);
	      check += compare_max ((double)iz, (double)izz,
				    label, verbose, tiny, &max);
	    }
	}
    }

  GRID_free (g);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/*
 * INPUT
 *  n       : number of elements
 *  list[n] : list of numbers
 *  flag    : 0 == ascending
 *            1 == descending
 */
static void
sort_list (int n, int *list, int flag)
{
  if (flag == 0)
    {
      // ascending order
      int i;
      for (i = 1; i < n; i ++)
	{
	  if (list[i] < list[i-1])
	    {
	      // swap the elements at i-1 and i
	      int l = list[i];
	      list[i] = list[i-1];
	      list[i-1] = l;
	      int j;
	      for (j = i-1; j >= 1; j --)
		{
		  if (list[j] < list[j-1])
		    {
		      // swap the elements at i-1 and i
		      l = list[j];
		      list[j] = list[j-1];
		      list[j-1] = l;
		    }
		  else
		    {
		      break;
		    }
		}
	    }
	}
    }
  else
    {
      // descending order
      int i;
      for (i = 1; i < n; i ++)
	{
	  if (list[i] > list[i-1])
	    {
	      // swap the elements at i-1 and i
	      int l = list[i];
	      list[i] = list[i-1];
	      list[i-1] = l;
	      int j;
	      for (j = i-1; j >= 1; j --)
		{
		  if (list[j] > list[j-1])
		    {
		      // swap the elements at i-1 and i
		      l = list[j];
		      list[j] = list[j-1];
		      list[j-1] = l;
		    }
		  else
		    {
		      break;
		    }
		}
	    }
	}
    }
}


/*
 * INPUT
 *  np : number of particles
 *  flag_conf : 0 == random configuration
 *              1 == regular array configuration
 */
int
check_GRID_init_all_by_cutoff (int np,
			       int flag_conf,
			       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_GRID_init_all_by_cutoff(np=%d,conf=%d): start\n",
	       np, flag_conf);
    }

  int check = 0;
  double max = 0.0;


  double r_cutoff = 1.0; // cut-off distance
  double r_cutoff2 = r_cutoff * r_cutoff;

  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_GRID_init_all_by_cutoff");
  sys->version = 0; // F version

  stokes_set_np (sys, np, np);


  double r0 = 0.5 * r_cutoff; // (mean) lattice distance
  int nx = (int)(pow ((double)np, 1.0 / 3.0)) + 1;

  int i;
  if (flag_conf == 0)
    {
      // random configuration
      double L = r0 * (double)nx;
      srand48 (0);
      for (i = 0; i < np; i ++)
	{
	  int ix = i*3;
	  sys->pos[ix  ] = L * drand48 ();
	  sys->pos[ix+1] = L * drand48 ();
	  sys->pos[ix+2] = L * drand48 ();
	}
    }
  else
    {
      // array configuration
      i = 0;
      int ix, iy, iz;
      for (ix = 0; ix < nx; ix ++)
	{
	  double x = (double)ix * r0;
	  for (iy = 0; iy < nx; iy ++)
	    {
	      double y = (double)iy * r0;
	      for (iz = 0; iz < nx; iz ++)
		{
		  double z = (double)iz * r0;

		  if (i >= np) break;
		  // now (i < np)

		  int ip = i*3;
		  sys->pos [ip  ] = x;
		  sys->pos [ip+1] = y;
		  sys->pos [ip+2] = z;

		  i++;
		}
	      if (i >= np) break;
	    }
	  if (i >= np) break;
	}
    }


  struct RYUON_grid *g = GRID_init_all_by_cutoff (sys, r_cutoff);
  CHECK_MALLOC (g, "check_GRID_init_all_by_cutoff");

  for (i = 0; i < np; i ++)
    {
      int ix = i * 3;

      // conventional approach
      int in1 = 0;
      int *ip1 = NULL;

      int j;
      for (j = i+1; j < np; j ++)
	{
	  int jx = j * 3;
	  double x = sys->pos [ix  ] - sys->pos [jx  ];
	  double y = sys->pos [ix+1] - sys->pos [jx+1];
	  double z = sys->pos [ix+2] - sys->pos [jx+2];

	  double r2 = x*x + y*y + z*z;
	  if (r2 < r_cutoff2)
	    {
	      in1 ++;
	      ip1 = (int *)realloc (ip1, sizeof (int) * in1);
	      CHECK_MALLOC (ip1, "check_GRID_init_all_by_cutoff");
	      ip1[in1-1] = j;
	    }
	}

      // grid approach
      int in2 = 0;
      int *ip2 = NULL;

      int k;
      for (k = 0; k < g->nnn[i]; k ++)
	{
	  int j = g->nnp[i][k];
	  // skip the self part and pair of (i > j)
	  if (i >= j) continue;

	  // j is in the nearest neighbor cell
	  in2 ++;
	  ip2 = (int *)realloc (ip2, sizeof (int) * in2);
	  CHECK_MALLOC (ip2, "check_GRID_init_all_by_cutoff");
	  ip2[in2-1] = j;
	}


      // compare between the conventional and grid approaches
      char label [80];
      sprintf (label, " n [%d]", i);
      check += compare_max ((double)in1, (double)in2,
			    label, verbose, tiny, &max);
      if (in1 == in2 && in1 > 0)
	{
	  // sort the lists in ascending order
	  sort_list (in1, ip1, 0);
	  sort_list (in2, ip2, 0);

	  for (k = 0; k < in1; k ++)
	    {
	      sprintf (label, " ip[%d] for %d", k, i);
	      check += compare_max ((double)ip1[k], (double)ip2[k],
				    label, verbose, tiny, &max);
	    }
	}

      // house-keeping
      if (ip1 != NULL) free (ip1);
      if (ip2 != NULL) free (ip2);
    }
  // end of for (i) -- particle loop


  GRID_free (g);
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
