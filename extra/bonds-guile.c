/* guile interface for struct BONDS
 * Copyright (C) 2007-2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include "stokes.h" // struct stokes
#include "stokes-guile.h" // ...
#include "bonds.h"

#include "bonds-guile.h"


/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     0       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      0      ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
 *      1.0    ;    p1   = A^{sp}, scaled spring constant
 *      2.1)   ;    p2   = L_{s} / length, scaled max extension
 *     ((0 1)  ; 3) list of pairs
 *      (1 2)
 *      (2 3))
 *      -1)    ; 4) number of exclusion for lubrication
 *             ;    negative means all particles in the chain is excluded.
 *    (; bond 2
 *     2       ; 1) spring type
 *     (       ; 2) spring parameters (list with 3 elements)
 *      1      ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
 *      19.8   ;    p1 = N_{K,s}, the Kuhn steps for a spring
 *      106.0) ;    p2 = b_{K} [nm], the Kuhn length
 *             ;    note that, for dWLC (type == 6),
 *             ;    (p1, p2) = (k, r0 [nm]), where the potential is
 *             ;    (k/2) * (kT / r0^2) * (r-r0)^2
 *     ((4 5)  ; 3) list of pairs
 *      (5 6)
 *      (6 7))
 *      1)     ; 4) number of exclusion for lubrication
 *    (; bond 3
 *     7       ; 1) spring type (FENE-Fraenkel)
 *     (       ; 2) spring parameters (list with 4 elements)
 *      0      ;    fene = 0 means (p1, p2, p3) = (H, r0 [nm], tol)
 *      1.0e6  ;    p1 = H, the spring constant
 *      0.5    ;    p2 = r0 [nm], the natural length of the spring
 *      0.01)  ;    p3 = tol, the tolerance parameter "s"
 *             ;    note that, for FENE-Fraenkel (type == 7),
 *             ;    the scalar part of the force is
 *             ;    fr = H * (r/hat(r0) - 1.0) / (1 - ((1-r/hat(r0))/tol)^2)
 *             ;    where hat(r0) = r0 / L0 (L0 is given by "length" [nm])
 *     ((8 9)  ; 3) list of pairs
 *      (9 10))
 *      0)     ; 4) number of exclusion for lubrication
 *   ))
 * where spring types are
 *   0 : Hookean spring (Asp * (r - Ls))
 *   1 : wormlike chain (WLC)
 *   2 : inverse Langevin chain (ILC)
 *   3 : Cohen's Pade approximation
 *   4 : Warner spring
 *   5 : Hookean spring (Asp * r / Ls)
 *   6 : Hookean spring for dWLC
 *   7 : FENE-Fraenkel
 * OUTPUT
 *  returned value : struct BONDS
 *                   if NULL is returned, it failed (not defined)
 */
struct BONDS *
BONDS_guile_get (const char * var)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "BONDS_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol = scm_c_lookup (var);

  SCM scm_bonds = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_bonds)))
    {
      fprintf (stderr, "BONDS_guile_get: %s is not a list\n", var);
      return (NULL);
    }

  //unsigned long len = scm_num2ulong (scm_length (scm_bonds),
  //				     0, "BONDS_guile_get");
  unsigned long len = scm_to_uint64 (scm_length (scm_bonds));
  if (len == 0)
    {
      // null is given
      return (NULL);
    }

  struct BONDS *b = BONDS_init ();
  CHECK_MALLOC (b, "BONDS_guile_get");

  int i;
  for (i = 0; i < len; i ++)
    {
      //SCM scm_bond = scm_list_ref (scm_bonds,
      //			   scm_int2num (i));
      SCM scm_bond = scm_list_ref (scm_bonds,
				   scm_from_int32 (i));
      if (!SCM_NFALSEP (scm_list_p (scm_bond)))
	{
	  // scm_bond is not a list
	  fprintf (stderr, "BONDS_guile_get:"
		   " %d-th bond of %s is not a list\n",
		   i, var);
	  BONDS_free (b);
	  return (NULL); // failed
	}
      unsigned long bond_len
	//= scm_num2ulong (scm_length (scm_bond),
	//		 0, "BONDS_guile_get");
	= scm_to_uint64 (scm_length (scm_bond));
      if (bond_len != 4)
	{
	  fprintf (stderr, "BONDS_guile_get:"
		   " length of %d-th bond of %s is not 4\n",
		   i, var);
	  BONDS_free (b);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_bond
      //int type = scm_num2int (scm_list_ref (scm_bond, scm_int2num (0)),
      //		      0,
      //		      "BONDS_guile_get");
      int type = scm_to_int32 (scm_list_ref (scm_bond, scm_from_int32 (0)));

      // 2nd element (1) of the list scm_bond
      //SCM scm_params = scm_list_ref (scm_bond, scm_int2num (1));
      SCM scm_params = scm_list_ref (scm_bond, scm_from_int32 (1));
      if (!SCM_NFALSEP (scm_list_p (scm_params)))
	{
	  // scm_params is not a list
	  fprintf (stderr, "BONDS_guile_get:"
		   " params of %d-th bond in %s is not a list\n",
		   i, var);
	  BONDS_free (b);
	  return (NULL); // failed
	}
      unsigned long params_len
	//= scm_num2ulong (scm_length (scm_params), 0, "BONDS_guile_get");
	= scm_to_uint64 (scm_length (scm_params));
      if (params_len != 3 &&
	  params_len != 4)
	{
	  fprintf (stderr, "BONDS_guile_get:"
		   " params has wrong length %ld != 3 nor 4\n",
		   params_len);
	  BONDS_free (b);
	  return (NULL); // failed
	}
      //int fene = scm_num2int (scm_list_ref (scm_params, scm_int2num (0)),
      //		      0,
      //		      "BONDS_guile_get");
      int fene = scm_to_int32 (scm_list_ref (scm_params, scm_from_int32 (0)));
      //double p1 = scm_num2dbl (scm_list_ref (scm_params, scm_int2num (1)),
      //			       "BONDS_guile_get");
      double p1 = scm_to_double (scm_list_ref (scm_params, scm_from_int32 (1)));
      //double p2 = scm_num2dbl (scm_list_ref (scm_params, scm_int2num (2)),
      //		       "BONDS_guile_get");
      double p2 = scm_to_double (scm_list_ref (scm_params, scm_from_int32 (2)));
      double p3 = 0.0;
      if (params_len == 4)
	{
	  // params_len == 4
	  //p3 = scm_num2dbl (scm_list_ref (scm_params, scm_int2num (3)),
	  //		    "BONDS_guile_get");
	  p3 = scm_to_double (scm_list_ref (scm_params, scm_from_int32 (3)));
	}
      // end of parsing scm_params


      // 4th element (3) of the list scm_bond
      // (3rd element will be taken later)
      //int nex = scm_num2int (scm_list_ref (scm_bond, scm_int2num (3)),
      //		     0,
      //		     "BONDS_guile_get");
      int nex = scm_to_int32 (scm_list_ref (scm_bond, scm_from_int32 (3)));


      // 3rd element (2) of the list scm_bond
      //SCM scm_pairs = scm_list_ref (scm_bond, scm_int2num (2));
      SCM scm_pairs = scm_list_ref (scm_bond, scm_from_int32 (2));
      if (!SCM_NFALSEP (scm_list_p (scm_pairs)))
	{
	  // scm_pairs is not a list
	  fprintf (stderr, "BONDS_guile_get:"
		   " pairs of %d-th bond in %s is not a list\n",
		   i, var);
	  BONDS_free (b);
	  return (NULL); // failed
	}

      unsigned long pairs_len
	//= scm_num2ulong (scm_length (scm_pairs),
	//			 0, "BONDS_guile_get");
	= scm_to_uint64 (scm_length (scm_pairs));
      int j;
      for (j = 0; j < pairs_len; j ++)
	{
	  //SCM scm_pair = scm_list_ref (scm_pairs,
	  //			       scm_int2num (j));
	  SCM scm_pair = scm_list_ref (scm_pairs,
	  			       scm_from_int32 (j));
	  if (!SCM_NFALSEP (scm_list_p (scm_pair)))
	    {
	      // scm_pair is not a list
	      fprintf (stderr, "BONDS_guile_get:"
		       " %d-th pair of %d-th bond of %s is not a list\n",
		       j, i, var);
	      BONDS_free (b);
	      return (NULL); // failed
	    }

	  unsigned long pair_len
	    //= scm_num2ulong (scm_length (scm_pair),
	    //			    0, "BONDS_guile_get");
	    = scm_to_uint64 (scm_length (scm_pair));
	  if (pair_len != 2)
	    {
	      fprintf (stderr, "BONDS_guile_get:"
		       " length of %d-th pair of %d-th bond of %s is not 2\n",
		       j, i, var);
	      BONDS_free (b);
	      return (NULL); // failed
	    }
	  
	  int ia, ib;
	  //ia = scm_num2int (scm_list_ref (scm_pair,
	  //				  scm_int2num (0)),
	  //		    0,
	  //		    "BONDS_guile_get");
	  //ib = scm_num2int (scm_list_ref (scm_pair,
	  //				  scm_int2num (1)),
	  //		    0,
	  //		    "BONDS_guile_get");
	  ia = scm_to_int32 (scm_list_ref (scm_pair,
					   scm_from_int32 (0)));
	  ib = scm_to_int32 (scm_list_ref (scm_pair,
					   scm_from_int32 (1)));

	  BONDS_append (b, type, fene, p1, p2, p3,
			ia, ib);
	}
    }

  return (b); // success
}
