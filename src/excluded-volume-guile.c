/* guile interface for struct EV
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume-guile.c,v 1.5 2008/07/17 05:35:19 kichiki Exp $
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
#include <libguile.h>
#include <math.h> // M_PI, sqrt()
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes-guile.h" // guile_load()
#include "excluded-volume.h" // struct EV

#include "excluded-volume-guile.h"


/* get ev from SCM and set struct EV
 * in SCM, ev are given by something like
 *  (define ev '(
 *   5.0     ; max distance [nm] (or in the same dimension of "length")
 *   ( ; for the EV 1
 *    0.0012 ; v [nm^3] (or in the same dimension of "length")
 *    0      ; fene
 *    1.0    ; p1 = A^{sp}, scaled spring const
 *    2.1    ; p2 = L_{s} / length, scaled max extension
 *    (0 1 2); list of particles belongs to the EV parameters
 *   )
 *   ( ; for the EV 2
 *    0.002  ; v [nm^3] (or in the same dimension of "length")
 *    1      ; fene
 *    19.8   ; p1 = N_{K,s}, the Kuhn steps for a spring
 *    106.0  ; p2 = b_{K} [nm], the Kuhn length
 *    (3 4)  ; list of particles belongs to the EV parameters
 *   )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-v".
 *  np     : number of particles (beads)
 *  length : unit length given by "length" in SCM (dimensional value)
 *  peclet : peclet number
 * OUTPUT
 *  returned value : struct EV
 *                   if NULL is returned, it failed (not defined)
 */
struct EV *
EV_guile_get (const char *var,
	      int np,
	      double length, double peclet)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "EV_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_EVs
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_EVs)))
    {
      fprintf (stderr, "EV_guile_get: %s is not a list\n", var);
      return (NULL);
    }

  unsigned long len
    = scm_num2ulong (scm_length (scm_EVs),
		     0, "EV_guile_get");
  if (len < 2)
    {
      // null list is given ==> no excluded-volume interaction
      return (NULL);
    }

  double *v = (double *)calloc (np, sizeof (double));
  double *N = (double *)calloc (np, sizeof (double));
  double *b = (double *)calloc (np, sizeof (double));
  CHECK_MALLOC (v, "EV_guile_get");
  CHECK_MALLOC (N, "EV_guile_get");
  CHECK_MALLOC (b, "EV_guile_get");

  // 1st element
  double rmax = scm_num2dbl (scm_list_ref (scm_EVs, scm_int2num (0)),
			     "EV_guile_get");
  rmax /= length; // dimensionless
  double r2 = rmax * rmax;

  int i;
  for (i = 1; i < len; i ++)
    {
      SCM scm_EV = scm_list_ref (scm_EVs, scm_int2num (i));
      if (!SCM_NFALSEP (scm_list_p (scm_EV)))
	{
	  // scm_EV is not a list
	  fprintf (stderr, "EV_guile_get:"
		   " %d-th EV of %s is not a list\n",
		   i-1, var);
	  return (NULL); // failed
	}
      unsigned long EV_len
	= scm_num2ulong (scm_length (scm_EV),
			 0, "EV_guile_get");
      if (EV_len != 5)
	{
	  fprintf (stderr, "EV_guile_get:"
		   " length of %d-th EV of %s is not 5\n",
		   i-1, var);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_EV
      double v0 = scm_num2dbl (scm_list_ref (scm_EV, scm_int2num (0)),
			       "EV_guile_get");
      // 2nd element (1) of the list scm_EV
      int fene = scm_num2int (scm_list_ref (scm_EV, scm_int2num (1)),
			      0,
			      "EV_guile_get");
      // 3rd element (2) of the list scm_EV
      double p1 = scm_num2dbl (scm_list_ref (scm_EV, scm_int2num (2)),
			       "EV_guile_get");
      // 4th element (3) of the list scm_EV
      double p2 = scm_num2dbl (scm_list_ref (scm_EV, scm_int2num (3)),
			       "EV_guile_get");

      double N_Ks, b_K;
      if (fene == 1)
	{
	  // given parameters
	  N_Ks = p1;
	  b_K  = p2; // dimensional
	}
      else
	{
	  // given parameters
	  double Asp = p1;
	  double Ls  = p2;

	  b_K = 3.0 / (peclet * Asp); // dimensionless
	  N_Ks = Ls / b_K;

	  b_K *= length; // dimensional
	}

      // 5th element (4) of the list scm_EV
      SCM scm_particles = scm_list_ref (scm_EV, scm_int2num (4));
      if (!SCM_NFALSEP (scm_list_p (scm_particles)))
	{
	  // scm_params is not a list
	  fprintf (stderr, "EV_guile_get:"
		   " params of %d-th EV in %s is not a list\n",
		   i-1, var);
	  return (NULL); // failed
	}
      unsigned long particles_len
	= scm_num2ulong (scm_length (scm_particles), 0, "EV_guile_get");
      int j;
      for (j = 0; j < particles_len; j ++)
	{
	  int ip = scm_num2int (scm_list_ref (scm_particles, scm_int2num (j)),
				0,
				"EV_guile_get");
	  if (ip < 0 || ip >= np)
	    {
	      fprintf (stderr, "EV_guile_get:"
		       " invalid particle index %d (np=%d)\n",
		       ip, np);
	      return (NULL); // failed
	    }

	  v[ip] = v0;
	  N[ip] = N_Ks;
	  b[ip] = b_K;
	}
    }

  struct EV *ev = EV_init (np, length, peclet, r2, v, N, b);
  CHECK_MALLOC (ev, "EV_guile_get");

  free (v);
  free (N);
  free (b);

  return (ev); // success
}
