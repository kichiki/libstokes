/* guile interface for struct EV_LJ.
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-LJ-guile.c,v 1.1 2008/05/24 05:47:33 kichiki Exp $
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
#include "ev-LJ.h"

#include "ev-LJ-guile.h"


/* get ev-LJ from SCM
 * in SCM, angles are given by something like
 *  (define ev-LJ '(
 *   (; LJ type 1
 *    10.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     0 1 2
 *    )
 *   )
 *   (; LJ type 2
 *    8.0 ; 1) LJ parameter epsilon in kT (so this is dimensionless value)
 *    2.0  ; 2) LJ parameter r0 in "length" (so this is dimensionless value)
 *    (    ; 3) list of particles
 *     3 4
 *    )
 *   )
 *  ))
 * INPUT
 *  var    : name of the variable.
 *           in the above example, set "ev-dh".
 *  peclet : peclet number
 *  np  : number of particles used for ev_dh_init()
 * OUTPUT
 *  returned value : struct EV_LJ
 *                   if NULL is returned, it failed (not defined)
 */
struct EV_LJ *
EV_LJ_guile_get (const char *var, int np)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "EV_LJ_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_LJs
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_LJs)))
    {
      fprintf (stderr, "EV_LJ_guile_get: %s is not a list\n", var);
      return (NULL);
    }

  unsigned long len
    = scm_num2ulong (scm_length (scm_LJs),
		     0, "EV_LJ_guile_get");
  if (len == 0)
    {
      // null is given
      return (NULL);
    }

  struct EV_LJ *ev_LJ = EV_LJ_init (np);
  CHECK_MALLOC (ev_LJ, "EV_LJ_guile_get");

  int i;
  for (i = 0; i < len; i ++)
    {
      SCM scm_LJ;
      scm_LJ = scm_list_ref (scm_LJs,
			     scm_int2num (i));
      if (!SCM_NFALSEP (scm_list_p (scm_LJ)))
	{
	  // scm_LJ is not a list
	  fprintf (stderr, "EV_LJ_guile_get:"
		   " %d-th LJ of %s is not a list\n",
		   i, var);
	  EV_LJ_free (ev_LJ);
	  return (NULL); // failed
	}
      unsigned long LJ_len;
      LJ_len = scm_num2ulong (scm_length (scm_LJ),
				0, "EV_LJ_guile_get");
      if (LJ_len != 3)
	{
	  fprintf (stderr, "EV_LJ_guile_get:"
		   " length of %d-th LJ of %s is not 3\n",
		   i, var);
	  EV_LJ_free (ev_LJ);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_LJ is (double)e
      double e  = scm_num2dbl (scm_list_ref (scm_LJ, scm_int2num (0)),
			       "EV_LJ_guile_get");
      // 2nd element (1) of the list scm_LJ is (double)r0
      double r0 = scm_num2dbl (scm_list_ref (scm_LJ, scm_int2num (1)),
			       "EV_LJ_guile_get");
      // 3rd element (2) of the list scm_LJ is list
      SCM scm_plist = scm_list_ref (scm_LJ, scm_int2num (2));
      if (!SCM_NFALSEP (scm_list_p (scm_plist)))
	{
	  // scm_params is not a list
	  fprintf (stderr, "EV_LJ_guile_get:"
		   " 3rd element of %d-th LJ in %s is not a list\n",
		   i, var);
	  EV_LJ_free (ev_LJ);
	  return (NULL); // failed
	}
      unsigned long plist_len
	= scm_num2ulong (scm_length (scm_plist), 0, "EV_LJ_guile_get");
      int j;
      for (j = 0; j < plist_len; j ++)
	{
	  int p = scm_num2int (scm_list_ref (scm_plist, scm_int2num (j)),
			       0,
			       "EV_LJ_guile_get");
	  if (p < 0 || p >= np)
	    {
	      fprintf (stderr, "EV_LJ_guile_get:"
		       " %d-th particle index %d of %d-th LJ of %s"
		       " is out of range [0,%d)\n",
		       j, p, i, var, np);
	      EV_LJ_free (ev_LJ);
	      return (NULL); // failed
	    }
	  ev_LJ->e[p]  = e;
	  ev_LJ->r0[p] = r0;
	}
    }

  return (ev_LJ); // success
}
