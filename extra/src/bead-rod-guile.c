/* guile interface for struct BeadRod
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <string.h> // strcmp()
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes-guile.h" // guile_load()
#include "bead-rod.h"     // struct BeadRod

#include "bead-rod-guile.h"


/* get constraints from SCM and set struct BeadRod
 * "constraints" is given in SCM as 
 *  (define constraints '(
 *   ; system parameters
 *   1.0e-6    ; 1) tolerance
 *   "NITSOL"  ; 2) scheme for solving nonlinear equations
 *             ;    "linear" for iterative scheme in linear approximation
 *             ;    "NITSOL" for Newton-GMRES scheme by NITSOL library
 *   ; the following is for each constraint
 *   (         ; 3) constraint type 1
 *    5.0      ; 3-1) distance [nm]
 *    (        ; 3-2) list of particle-pairs
 *     (0 1)
 *     (1 2)
 *     (2 3)
 *   ))
 *   (         ; 4) constraint type 2
 *    10.0     ; 4-1) distance [nm]
 *    (        ; 4-2) list of particle-pairs
 *     (3 4)
 *     (4 5)
 *   ))
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "constraints".
 *  sys : struct stokes
 *  length : unit length in the simulation
 * OUTPUT
 *  returned value : struct BeadRod
 *                   if NULL is returned, it failed (not defined)
 */
struct BeadRod *
BeadRod_guile_get (const char *var,
		   struct stokes *sys,
		   double length)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "BeadRod_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_brs
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_brs)))
    {
      fprintf (stderr, "BeadRod_guile_get: %s is not a list\n", var);
      return (NULL);
    }

  unsigned long len
    //= scm_num2ulong (scm_length (scm_brs),
    //		     0, "BeadRod_guile_get");
    = scm_to_uint64 (scm_length (scm_brs));
  if (len <= 2)
    {
      // null list is given ==> no constraint
      return (NULL);
    }

  // 1st element (0) of the list scm_brs
  //double eps = scm_num2dbl (scm_list_ref (scm_brs, scm_int2num (0)),
  //			    "BeadRod_guile_get");
  double eps = scm_to_double (scm_list_ref (scm_brs, scm_from_int32 (0)));
  // 2nd element (1) of the list scm_brs
  // get the string
  char *str_scheme = NULL;
  //SCM scm_scheme = scm_list_ref (scm_brs, scm_int2num (1));
  SCM scm_scheme = scm_list_ref (scm_brs, scm_from_int32 (1));
  if (scm_is_string (scm_scheme))
    {
      str_scheme = scm_to_locale_string (scm_scheme);
    }

  int scheme;
  if (strcmp (str_scheme, "linear") == 0)
    {
      scheme = 0;
    }
  else if (strcmp (str_scheme, "nitsol") == 0 ||
	   strcmp (str_scheme, "NITSOL") == 0)
    {
      scheme = 1;
    }
  else
    {
      fprintf (stderr, "BeadRod_guile_get: invalid scheme %s\n", str_scheme);
      return (NULL);
    }


  double *a = NULL;
  int *ia = NULL;
  int *ib = NULL;
  int nc = 0;

  // (len - 2) is the number of types of constraints
  int i;
  for (i = 2; i < len; i ++)
    {
      // (i-1)-th constraint
      SCM scm_br
	//= scm_list_ref (scm_brs,
	//		scm_int2num (i));
	= scm_list_ref (scm_brs,
			scm_from_int32 (i));
      if (!SCM_NFALSEP (scm_list_p (scm_br)))
	{
	  // scm_br is not a list
	  fprintf (stderr, "BeadRod_guile_get:"
		   " %d-th constraint of %s is not a list\n",
		   i-1, var);
	  if (a  != NULL) free (a);
	  if (ia != NULL) free (ia);
	  if (ib != NULL) free (ib);
	  return (NULL); // failed
	}
      unsigned long br_len
	//= scm_num2ulong (scm_length (scm_br),
	//		 0, "BeadRod_guile_get");
	= scm_to_uint64 (scm_length (scm_br));
      if (br_len != 2)
	{
	  fprintf (stderr, "BeadRod_guile_get:"
		   " length of %d-th constraint of %s is not 2\n",
		   i-1, var);
	  if (a  != NULL) free (a);
	  if (ia != NULL) free (ia);
	  if (ib != NULL) free (ib);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_br
      //double a0 = scm_num2dbl (scm_list_ref (scm_br, scm_int2num (0)),
      //		       "BeadRod_guile_get");
      double a0 = scm_to_double (scm_list_ref (scm_br, scm_from_int32 (0)));
      a0 /= length; // scale by the unit length

      // 2nd element (1) of the list scm_br
      SCM scm_pairs
	//= scm_list_ref (scm_br,	scm_int2num (1));
	= scm_list_ref (scm_br,	scm_from_int32 (1));

      if (!SCM_NFALSEP (scm_list_p (scm_pairs)))
	{
	  // scm_pairs is not a list
	  fprintf (stderr, "BeadRod_guile_get:"
		   " 2nd element of %d-th constraint of %s is not a list\n",
		   i-1, var);
	  if (a  != NULL) free (a);
	  if (ia != NULL) free (ia);
	  if (ib != NULL) free (ib);
	  return (NULL); // failed
	}
      unsigned long pairs_len
	//= scm_num2ulong (scm_length (scm_pairs),
	//		 0, "BeadRod_guile_get");
	= scm_to_uint64 (scm_length (scm_pairs));
      int j;
      for (j = 0; j < pairs_len; j ++)
	{
	  SCM scm_pair
	    //= scm_list_ref (scm_pairs,
	    //		    scm_int2num (j));
	    = scm_list_ref (scm_pairs,
	    		    scm_from_int32 (j));
	  unsigned long pair_len
	    //= scm_num2ulong (scm_length (scm_pair),
	    //		     0, "BeadRod_guile_get");
	    = scm_to_uint64 (scm_length (scm_pair));
	  if (pair_len != 2)
	    {
	      fprintf (stderr, "BeadRod_guile_get:"
		       " length of %d-th pair for %d-th constraint of %s"
		       " is not 2\n",
		       j, i-1, var);
	      if (a  != NULL) free (a);
	      if (ia != NULL) free (ia);
	      if (ib != NULL) free (ib);
	      return (NULL); // failed
	    }
	  // 1st element (0) of the list scm_pair
	  //int ia0 = scm_num2int (scm_list_ref (scm_pair, scm_int2num (0)),
	  //			 0,
	  //			 "BeadRod_guile_get");
	  int ia0 = scm_to_int32 (scm_list_ref (scm_pair, scm_from_int32 (0)));
	  // 2nd element (1) of the list scm_pair
	  //int ib0 = scm_num2int (scm_list_ref (scm_pair, scm_int2num (1)),
	  //			 0,
	  //			 "BeadRod_guile_get");
	  int ib0 = scm_to_int32 (scm_list_ref (scm_pair, scm_from_int32 (1)));
	  nc ++;
	  a  = (double *)realloc (a, sizeof (double) * nc);
	  ia = (int *)realloc (ia, sizeof (int) * nc);
	  ib = (int *)realloc (ib, sizeof (int) * nc);
	  a [nc-1] = a0;
	  ia[nc-1] = ia0;
	  ib[nc-1] = ib0;
	}

    }

  // now we have nc (number of constraints) and tables a[nc], ia[nc], ib[nc]
  struct BeadRod *br = BeadRod_init (sys, nc, a, ia, ib);
  CHECK_MALLOC (br, "BeadRod_guile_get");
  free (a);
  free (ia);
  free (ib);

  BeadRod_set_scheme (br, scheme, eps);

  return (br); // success
}
