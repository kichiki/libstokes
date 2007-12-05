/* test code for struct list_ex in bonds.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-list-ex.c,v 1.1 2007/12/05 03:53:06 kichiki Exp $
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
#include <stdlib.h>
#include <math.h>
#include "memory-check.h"
#include "check.h" // compare()

#include <bonds.h>


static void
list_ex_print (struct list_ex *ex)
{
  fprintf (stdout, " ex->np = %d\n", ex->np);
  int i;
  for (i = 0; i < ex->np; i ++)
    {
      fprintf (stdout, " ex->n[%d] = %d :", i, ex->n[i]);
      int j;
      for (j = 0; j < ex->n[i]; j ++)
	{
	  fprintf (stdout, " %d", ex->i[i][j]);
	}
      fprintf (stdout, "\n");
    }
}


int
check_list_ex (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_list_ex : start\n");
    }

  int check = 0;
  //double max = 0.0;


  struct bonds *b = bonds_init ();
  CHECK_MALLOC (b, "check_list_ex");

  /* set struct bonds *b by the following parameters:
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
   */
  // linear chain (0)-(1)-(2)-(3) with nex = -1
  bonds_add_type (b,
		  0,   // type
		  0,   // fene
		  1.0, // p1 = A^{sp}
		  2.1, // p2 = L_{s}
		  -1); // nex
  bond_pairs_add (b->pairs[0], 0, 1);
  bond_pairs_add (b->pairs[0], 1, 2);
  bond_pairs_add (b->pairs[0], 2, 3);

  // linear chain (4)-(5)-(6)-(7) with nex = 1
  bonds_add_type (b,
		  2,    // type
		  1,    // fene
		  19.8, // p1 = N_{K,s}
		  106.0,// p2 = b_{K}
		  1);   // nex
  bond_pairs_add (b->pairs[1], 4, 5);
  bond_pairs_add (b->pairs[1], 5, 6);
  bond_pairs_add (b->pairs[1], 6, 7);

  // linear chain (8)-(9)-(10)-(11) with nex = 2
  bonds_add_type (b,
		  2,    // type
		  1,    // fene
		  19.8, // p1 = N_{K,s}
		  106.0,// p2 = b_{K}
		  2);   // nex
  bond_pairs_add (b->pairs[2], 8,  9);
  bond_pairs_add (b->pairs[2], 9,  10);
  bond_pairs_add (b->pairs[2], 10, 11);

  /* complex bonds with nex = 2
   *      (12)
   *       |
   * (13)-(14)-(15)
   *       |
   *      (16)-(17)
   *       |    |
   *      (18)-(19)
   */
  bonds_add_type (b,
		  2,    // type
		  1,    // fene
		  19.8, // p1 = N_{K,s}
		  106.0,// p2 = b_{K}
		  2);   // nex
  bond_pairs_add (b->pairs[3], 12, 14);
  bond_pairs_add (b->pairs[3], 13, 14);
  bond_pairs_add (b->pairs[3], 14, 15);
  bond_pairs_add (b->pairs[3], 14, 16);
  bond_pairs_add (b->pairs[3], 16, 17);
  bond_pairs_add (b->pairs[3], 17, 19);
  bond_pairs_add (b->pairs[3], 19, 18);
  bond_pairs_add (b->pairs[3], 18, 16);

  int np = 20;
  struct list_ex *ex = list_ex_init (np);
  CHECK_MALLOC (ex, "check_list_ex");

  list_ex_set_by_bonds (ex, b);


  /**
   * expected results
   */
  /* linear chain (0)-(1)-(2)-(3) with nex = -1
    ex->n[0] = 3 : 1 2 3
    ex->n[1] = 3 : 0 2 3
    ex->n[2] = 3 : 1 3 0
    ex->n[3] = 3 : 2 1 0
  */
  if (ex->n[0] != 3) check ++;
  if (ex->n[1] != 3) check ++;
  if (ex->n[2] != 3) check ++;
  if (ex->n[3] != 3) check ++;

  /* linear chain (4)-(5)-(6)-(7) with nex = 1
    ex->n[4] = 1 : 5
    ex->n[5] = 2 : 4 6
    ex->n[6] = 2 : 5 7
    ex->n[7] = 1 : 6
  */
  if (ex->n[4] != 1) check ++;
  if (ex->n[5] != 2) check ++;
  if (ex->n[6] != 2) check ++;
  if (ex->n[7] != 1) check ++;

  /* linear chain (8)-(9)-(10)-(11) with nex = 2
    ex->n[8] = 2 : 9 10
    ex->n[9] = 3 : 8 10 11
    ex->n[10] = 3 : 9 11 8
    ex->n[11] = 2 : 10 9
  */
  if (ex->n[8]  != 2) check ++;
  if (ex->n[9]  != 3) check ++;
  if (ex->n[10] != 3) check ++;
  if (ex->n[11] != 2) check ++;

  /* complex bonds with nex = 2
   *      (12)
   *       |
   * (13)-(14)-(15)
   *       |
   *      (16)-(17)
   *       |    |
   *      (18)-(19)
   *
    ex->n[12] = 4 : 14 13 15 16
    ex->n[13] = 4 : 14 12 15 16
    ex->n[14] = 6 : 12 13 15 16 17 18
    ex->n[15] = 4 : 14 12 13 16
    ex->n[16] = 7 : 14 17 18 12 13 15 19
    ex->n[17] = 4 : 16 19 14 18
    ex->n[18] = 4 : 19 16 17 14
    ex->n[19] = 3 : 17 18 16
  */
  if (ex->n[12] != 4) check ++;
  if (ex->n[13] != 4) check ++;
  if (ex->n[14] != 6) check ++;
  if (ex->n[15] != 4) check ++;
  if (ex->n[16] != 7) check ++;
  if (ex->n[17] != 4) check ++;
  if (ex->n[18] != 4) check ++;
  if (ex->n[19] != 3) check ++;

  if (verbose != 0)
    {
      list_ex_print (ex);
    }


  bonds_free (b);
  list_ex_free (ex);

  if (verbose != 0)
    {
      //fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
