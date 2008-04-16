/* guile interface for struct angles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: angles-guile.c,v 1.2 2008/04/16 00:38:43 kichiki Exp $
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
#include <math.h> // M_PI, sqrt()
#include "memory-check.h" // macro CHECK_MALLOC

#include <libguile.h>
#include <guile/gh.h>
#include "stokes-guile.h" // guile_load()
#include "angles.h"

#include "angles-guile.h"


/* get angles from SCM
 * in SCM, angles are given by something like
 *  (define angles '(
 *    (; angle type 1
 *     10.0    ; 1) constant (k^{angle})
 *     0.0     ; 2) angle in degree (theta_0)
 *     ((0 1 2); 3) list of triplets
 *      (1 2 3)
 *      (2 3 4)
 *     )
 *    )
 *    (; angle type 2
 *     20.0    ; 1) constant (k^{angle})
 *     90.0    ; 2) angle in degree (theta_0)
 *     ((3 4 5); 3) list of triplets
 *      (4 5 6)
 *     )
 *    )
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "angles".
 * OUTPUT
 *  returned value : struct angles
 *                   if NULL is returned, it failed (not defined)
 */
struct angles *
guile_get_angles (const char *var)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "guile_get_angles: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_angles
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_angles)))
    {
      fprintf (stderr, "guile_get_angles: %s is not a list\n", var);
      return (NULL);
    }


  struct angles *angles = angles_init ();
  CHECK_MALLOC (angles, "guile_get_angles");

  unsigned long len
    = scm_num2ulong (scm_length (scm_angles),
		     0, "guile_get_angles");
  int i;
  for (i = 0; i < len; i ++)
    {
      SCM scm_angle
	= scm_list_ref (scm_angles,
			scm_int2num (i));
      if (!SCM_NFALSEP (scm_list_p (scm_angle)))
	{
	  // scm_angle is not a list
	  fprintf (stderr, "guile_get_angles:"
		   " %d-th bond of %s is not a list\n",
		   i, var);
	  angles_free (angles);
	  return (NULL); // failed
	}
      unsigned long angle_len
	= scm_num2ulong (scm_length (scm_angle),
			 0, "guile_get_angles");
      if (angle_len != 3)
	{
	  fprintf (stderr, "guile_get_angles:"
		   " length of %d-th bond of %s is not 3\n",
		   i, var);
	  angles_free (angles);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_angle
      double k = scm_num2dbl (scm_list_ref (scm_angle, scm_int2num (0)),
			      "guile_get_angles");

      // 2nd element (1) of the list scm_angle
      double t0 = scm_num2dbl (scm_list_ref (scm_angle, scm_int2num (1)),
			       "guile_get_angles");
      // convert degree to radian
      t0 = t0 * M_PI / 180.0;

      // 3rd element (2) of the list scm_angle
      SCM scm_triplets = scm_list_ref (scm_angle, scm_int2num (2));
      if (!SCM_NFALSEP (scm_list_p (scm_triplets)))
	{
	  // scm_triplets is not a list
	  fprintf (stderr, "guile_get_angles:"
		   " pairs of %d-th bond in %s is not a list\n",
		   i, var);
	  angles_free (angles);
	  return (NULL); // failed
	}

      unsigned long triplets_len;
      triplets_len = scm_num2ulong (scm_length (scm_triplets),
				    0, "guile_get_angles");
      int j;
      for (j = 0; j < triplets_len; j ++)
	{
	  SCM scm_triplet = scm_list_ref (scm_triplets,
					  scm_int2num (j));
	  if (!SCM_NFALSEP (scm_list_p (scm_triplet)))
	    {
	      // scm_triplet is not a list
	      fprintf (stderr, "guile_get_angles:"
		       " %d-th pair of %d-th bond of %s is not a list\n",
		       j, i, var);
	      angles_free (angles);
	      return (NULL); // failed
	    }

	  unsigned long triplet_len;
	  triplet_len = scm_num2ulong (scm_length (scm_triplet),
				       0, "guile_get_angles");
	  if (triplet_len != 3)
	    {
	      fprintf (stderr, "guile_get_angles:"
		       " length of %d-th triplet of %d-th bond of %s"
		       " is not 3\n",
		       j, i, var);
	      angles_free (angles);
	      return (NULL); // failed
	    }
	  
	  int ia = scm_num2int (scm_list_ref (scm_triplet,
					      scm_int2num (0)),
				0,
				"guile_get_angles");
	  int ib = scm_num2int (scm_list_ref (scm_triplet,
					      scm_int2num (1)),
				0,
				"guile_get_angles");
	  int ic = scm_num2int (scm_list_ref (scm_triplet,
					      scm_int2num (2)),
				0,
				"guile_get_angles");

	  angles_add (angles, ia, ib, ic, k, t0);
	}
    }

  return (angles); // success
}
