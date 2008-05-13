/* guile interface for struct EV
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume-guile.c,v 1.1 2008/05/13 01:11:54 kichiki Exp $
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
#include "bonds.h"           // struct bonds
#include "excluded-volume.h" // struct EV

#include "excluded-volume-guile.h"


/* get ev-v from SCM and set struct EV
 * in SCM, ev-v is a list of parameter v [nm^3] or [micro m^3]
 * (depending on the dimension of the parameter "length")
 * for each spring:
 *  (define ev-v '(
 *   0.0012 ; for the spring 1
 *   0.002  ; for the spring 2
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "ev-v".
 *  bonds : struct bonds
 *  length : unit length given by "length" in SCM (dimensional value)
 *  peclet : peclet number
 *  ev_r2  : square of max distance for EV interaction
 *  np     : number of particles (beads)
 * OUTPUT
 *  returned value : struct EV
 *                   if NULL is returned, it failed (not defined)
 */
struct EV *
guile_get_ev_v (const char *var,
		const struct bonds *bonds,
		double length, double peclet,
		double ev_r2, int np)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "guile_get_ev_v: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_ev
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_ev)))
    {
      fprintf (stderr, "guile_get_ev_v: %s is not a list\n", var);
      return (NULL);
    }

  unsigned long len
    = scm_num2ulong (scm_length (scm_ev),
		     0, "guile_get_ev_v");
  if (len == 0)
    {
      // null list is given ==> no excluded-volume interaction
      return (NULL);
    }
  else if (len != bonds->n)
    {
      fprintf (stderr, "guile_get_ev_v:"
	       " length of %s %ld != %d number of springs\n",
	       var, len, bonds->n);
      return (NULL);
    }

  double *ev_v = (double *)malloc (sizeof (double) * len);
  CHECK_MALLOC (ev_v, "guile_get_ev_v");

  if (guile_get_doubles (var, len, ev_v) != 1) // FALSE
    {
      fprintf (stderr, "guile_get_ev_v: fail to parse %s\n", var);
      return (NULL);
    }

  struct EV *ev = EV_init (bonds, length, peclet, ev_r2, ev_v, np);
  CHECK_MALLOC (ev, "guile_get_ev_v");

  return (ev); // success
}
