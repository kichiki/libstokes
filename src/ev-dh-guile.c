/* guile interface for struct EV_DH.
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-guile.c,v 1.5 2008/06/13 03:10:32 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes-guile.h" // guile_load()
#include "ev-dh.h"

#include "ev-dh-guile.h"


/* get ev_dh from SCM
 * in SCM, "ev-dh" are given by something like
 *  (define ev-dh '(
 *    ; system parameters
 *    1.0e-6   ; 1) epsilon for the cut-off distance of EV_DH interaction
 *    298.0    ; 2) temperature [K]
 *    80.0     ; 3) dielectric constant of the solution
 *    3.07     ; 4) Debye length [nm]
 *    (        ; 5) list of chain types
 *     (; chain type 1
 *      2.43    ; 1) nu [e/nm]
 *      5.00    ; 2) l0 [nm]
 *      (0 1 2) ; 3) list of particles
 *     )
 *     (; chain type 2
 *      2.00    ; 1) nu [e/nm]
 *      4.00    ; 2) l0 [nm]
 *      (3 4)   ; 3) list of particles
 *     )
 *    )
 *  ))
 * INPUT
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  length : the characteristic length (dimensional, usually in nm)
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_DH
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_DH *
EV_DH_guile_get (const char *var,
		 double length, double peclet, int np)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "EV_DH_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_ev_dh
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_ev_dh)))
    {
      fprintf (stderr, "EV_DH_guile_get: %s is not a list\n", var);
      return (NULL);
    }

  unsigned long len
    = scm_num2ulong (scm_length (scm_ev_dh),
		     0, "EV_DH_guile_get");
  if (len == 0)
    {
      // null is given
      return (NULL);
    }
  else if (len != 5)
    {
      fprintf (stderr, "EV_DH_guile_get:"
	       " length of %s is not 5\n", var);
      return (NULL); // failed
    }

  // 1st element (0) of the list scm_angle
  double eps = scm_num2dbl (scm_list_ref (scm_ev_dh, scm_int2num (0)),
			    "EV_DH_guile_get");
  // 2nd element (1) of the list scm_angle
  double T = scm_num2dbl (scm_list_ref (scm_ev_dh, scm_int2num (1)),
			  "EV_DH_guile_get");
  // 3rd element (2) of the list scm_ev_dh
  double e = scm_num2dbl (scm_list_ref (scm_ev_dh, scm_int2num (2)),
			  "EV_DH_guile_get");
  // 4th element (3) of the list scm_ev_dh
  double rd = scm_num2dbl (scm_list_ref (scm_ev_dh, scm_int2num (3)),
			   "EV_DH_guile_get");
  // note that rd should be dimensional for EV_DH_init()

  struct EV_DH *ev_dh = EV_DH_init (length, peclet, rd, T, e, eps, np);
  CHECK_MALLOC (ev_dh, "EV_DH_guile_get");

  // 5th element (4) of the list scm_ev_dh
  SCM scm_chains
    = scm_list_ref (scm_ev_dh, scm_int2num (4));
  if (!SCM_NFALSEP (scm_list_p (scm_chains)))
    {
      // scm_chains is not a list
      fprintf (stderr, "EV_DH_guile_get:"
	       " 5th element of %s is not a list\n",
	       var);
      EV_DH_free (ev_dh);
      return (NULL); // failed
    }
  unsigned long chains_len
    = scm_num2ulong (scm_length (scm_chains),
		     0, "EV_DH_guile_get");

  int i;
  for (i = 0; i < chains_len; i ++)
    {
      SCM scm_chain
	= scm_list_ref (scm_chains,
			scm_int2num (i));
      if (!SCM_NFALSEP (scm_list_p (scm_chain)))
	{
	  // scm_chain is not a list
	  fprintf (stderr, "EV_DH_guile_get:"
		   " %d-th chain of %s is not a list\n",
		   i, var);
	  EV_DH_free (ev_dh);
	  return (NULL); // failed
	}
      unsigned long chain_len
	= scm_num2ulong (scm_length (scm_chain),
			 0, "EV_DH_guile_get");
      if (chain_len != 3)
	{
	  fprintf (stderr, "EV_DH_guile_get:"
		   " length of %d-th chain of %s is not 3\n",
		   i, var);
	  EV_DH_free (ev_dh);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_chain
      double nu = scm_num2dbl (scm_list_ref (scm_chain, scm_int2num (0)),
			       "EV_DH_guile_get");

      // 2nd element (1) of the list scm_chain
      double l0 = scm_num2dbl (scm_list_ref (scm_chain, scm_int2num (1)),
			       "EV_DH_guile_get");

      // 3rd element (2) of the list scm_chain
      SCM scm_particles
	= scm_list_ref (scm_chain, scm_int2num (2));
      if (!SCM_NFALSEP (scm_list_p (scm_particles)))
	{
	  // scm_particles is not a list
	  fprintf (stderr, "EV_DH_guile_get:"
		   " particle-list of %d-th chain of %s is not a list\n",
		   i, var);
	  EV_DH_free (ev_dh);
	  return (NULL); // failed
	}
      unsigned long particles_len
	= scm_num2ulong (scm_length (scm_particles),
			 0, "EV_DH_guile_get");
      int j;
      for (j = 0; j < particles_len; j++)
	{
	  int p = scm_num2int (scm_list_ref (scm_particles, scm_int2num (j)),
			       0,
			       "EV_DH_guile_get");
	  if (p < 0 || p >= np)
	    {
	      fprintf (stderr, "EV_DH_guile_get:"
		       " %d-th particle index %d of %d-th chain of %s"
		       " is out of range [0,%d)\n",
		       j, p, i, var, np);
	      EV_DH_free (ev_dh);
	      return (NULL); // failed
	    }
	  ev_dh->nu[p] = nu * l0;
	}
    }

  return (ev_dh); // success
}
