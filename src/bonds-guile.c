/* guile interface for struct bonds
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-guile.c,v 1.6 2008/04/16 03:07:15 kichiki Exp $
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
#include "stokes.h" // struct stokes
#include "stokes-guile.h" // ...
#include "bonds.h"

#include "bonds-guile.h"


/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     0         ; 1) spring type
 *     (         ; 2) spring parameters (list with 3 elements)
 *      0        ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
 *      1.0      ;    p1   = A^{sp}, scaled spring constant  (for fene == 0)
 *      2.1)     ;    p2   = L_{s} / a, scaled max extension (for fene == 0)
 *     ((0 1)    ; 3) list of pairs
 *      (1 2)
 *      (2 3))
 *      -1)      ; 4) number of exclusion for lubrication
 *               ;    negative means all particles in the chain is excluded.
 *    (; bond 2
 *     2         ; 1) spring type
 *     (         ; 2) spring parameters (list with 3 elements)
 *      1        ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
 *      19.8     ;    p1 = N_{K,s}, the Kuhn steps for a spring (for fene = 1)
 *      106.0)   ;    p2 = b_{K} [nm], the Kuhn length          (for fene = 1)
 *     ((4 5)    ; 3) list of pairs
 *      (5 6)
 *      (6 7))
 *       1)      ; 4) number of exclusion for lubrication
 *   ))
 * where spring types are
 *   0 : Hookean spring (Asp * (r - Ls)
 *   1 : wormlike chain (WLC)
 *   2 : inverse Langevin chain (ILC)
 *   3 : Cohen's Pade approximation
 *   4 : Warner spring
 *   5 : Hookean spring (Asp * r / Ls)
 * OUTPUT
 *  returned value : struct bonds
 *                   if NULL is returned, it failed (not defined)
 */
struct bonds *
guile_get_bonds (const char * var)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "guile_get_bonds: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol;
  scm_symbol = scm_c_lookup (var);

  SCM scm_bonds;
  scm_bonds = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_bonds)))
    {
      fprintf (stderr, "guile_get_bonds: %s is not a list\n", var);
      return (NULL);
    }


  struct bonds *bonds = NULL;
  bonds = bonds_init ();

  unsigned long len;
  len = scm_num2ulong (scm_length (scm_bonds),
		       0, "guile_get_bonds");
  int i;
  for (i = 0; i < len; i ++)
    {
      SCM scm_bond;
      scm_bond = scm_list_ref (scm_bonds,
			       scm_int2num (i));
      if (!SCM_NFALSEP (scm_list_p (scm_bond)))
	{
	  // scm_bond is not a list
	  fprintf (stderr, "guile_get_bonds:"
		   " %d-th bond of %s is not a list\n",
		   i, var);
	  bonds_free (bonds);
	  return (NULL); // failed
	}
      unsigned long bond_len;
      bond_len = scm_num2ulong (scm_length (scm_bond),
				0, "guile_get_bonds");
      if (bond_len != 4)
	{
	  fprintf (stderr, "guile_get_bonds:"
		   " length of %d-th bond of %s is not 4\n",
		   i, var);
	  bonds_free (bonds);
	  return (NULL); // failed
	}

      // 1st element (0) of the list scm_bond
      int type = scm_num2int (scm_list_ref (scm_bond, scm_int2num (0)),
			      0,
			      "guile_get_bonds");

      // 2nd element (1) of the list scm_bond
      SCM scm_params = scm_list_ref (scm_bond, scm_int2num (1));
      if (!SCM_NFALSEP (scm_list_p (scm_params)))
	{
	  // scm_params is not a list
	  fprintf (stderr, "guile_get_bonds:"
		   " params of %d-th bond in %s is not a list\n",
		   i, var);
	  bonds_free (bonds);
	  return (NULL); // failed
	}
      unsigned long params_len
	= scm_num2ulong (scm_length (scm_params), 0, "guile_get_bonds");
      if (params_len != 3)
	{
	  fprintf (stderr, "guile_get_bonds:"
		   " params has wrong length %ld != 3\n",
		   params_len);
	  bonds_free (bonds);
	  return (NULL); // failed
	}
      int fene = scm_num2int (scm_list_ref (scm_params, scm_int2num (0)),
			      0,
			      "guile_get_bonds");
      double p1 = scm_num2dbl (scm_list_ref (scm_params, scm_int2num (1)),
			       "guile_get_bonds");
      double p2 = scm_num2dbl (scm_list_ref (scm_params, scm_int2num (2)),
			       "guile_get_bonds");

      // 4th element (3) of the list scm_bond
      // (3rd element will be taken later)
      int nex = scm_num2int (scm_list_ref (scm_bond, scm_int2num (3)),
			     0,
			     "guile_get_bonds");

      bonds_add_type (bonds, type, fene, p1, p2, nex);


      // 3rd element (2) of the list scm_bond
      SCM scm_pairs = scm_list_ref (scm_bond, scm_int2num (2));
      if (!SCM_NFALSEP (scm_list_p (scm_pairs)))
	{
	  // scm_pairs is not a list
	  fprintf (stderr, "guile_get_bonds:"
		   " pairs of %d-th bond in %s is not a list\n",
		   i, var);
	  bonds_free (bonds);
	  return (NULL); // failed
	}

      unsigned long pairs_len;
      pairs_len = scm_num2ulong (scm_length (scm_pairs),
				 0, "guile_get_bonds");
      int j;
      for (j = 0; j < pairs_len; j ++)
	{
	  SCM scm_pair;
	  scm_pair = scm_list_ref (scm_pairs,
				   scm_int2num (j));
	  if (!SCM_NFALSEP (scm_list_p (scm_pair)))
	    {
	      // scm_pair is not a list
	      fprintf (stderr, "guile_get_bonds:"
		       " %d-th pair of %d-th bond of %s is not a list\n",
		       j, i, var);
	      bonds_free (bonds);
	      return (NULL); // failed
	    }

	  unsigned long pair_len;
	  pair_len = scm_num2ulong (scm_length (scm_pair),
				    0, "guile_get_bonds");
	  if (pair_len != 2)
	    {
	      fprintf (stderr, "guile_get_bonds:"
		       " length of %d-th pair of %d-th bond of %s is not 2\n",
		       j, i, var);
	      bonds_free (bonds);
	      return (NULL); // failed
	    }
	  
	  int ia, ib;
	  ia = scm_num2int (scm_list_ref (scm_pair,
					  scm_int2num (0)),
			    0,
			    "guile_get_bonds");
	  ib = scm_num2int (scm_list_ref (scm_pair,
					  scm_int2num (1)),
			    0,
			    "guile_get_bonds");

	  bond_pairs_add (bonds->pairs [i], ia, ib);
	}
    }

  return (bonds); // success
}
