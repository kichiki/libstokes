/* guile interface for struct bonds
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-guile.c,v 1.1 2007/03/27 00:53:57 kichiki Exp $
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
#include <guile/gh.h>
#include "stokes.h" // struct stokes
#include "stokes-guile.h" // ...
#include "bonds.h"

#include "bonds-guile.h"


/* get bonds from SCM
 * in SCM, bonds are something like
 *  (define bonds '(
 *    (; bond 1
 *     1.0 ; spring const
 *     2.1 ; natural distance
 *     ((0 1) ; list of pairs
 *      (1 2)
 *      (2 3)))
 *    (; bond 2
 *     1.0 ; spring const
 *     2.5 ; natural distance
 *     ((4 5) ; list of pairs
 *      (5 6)
 *      (6 7)))
 *   ))
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
      if (bond_len != 3)
	{
	  fprintf (stderr, "guile_get_bonds:"
		   " length of %d-th bond of %s is not 3\n",
		   i, var);
	  bonds_free (bonds);
	  return (NULL); // failed
	}

      double k, r0;
      k  = scm_num2dbl (scm_list_ref (scm_bond,
				      scm_int2num (0)),
			"guile_get_bonds");
      r0 = scm_num2dbl (scm_list_ref (scm_bond,
				      scm_int2num (1)),
			"guile_get_bonds");
      bonds_add_type (bonds, k, r0);


      SCM scm_pairs;
      scm_pairs = scm_list_ref (scm_bond,
				scm_int2num (2));
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
