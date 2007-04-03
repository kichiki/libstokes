/* test code for libstokes
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-all.c,v 1.2 2007/04/03 02:38:04 kichiki Exp $
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
#include "check-lapack-inv.h" // check_lapack_inv_()
#include "check-mul-matrices.h" // check_mul_matrices()
#include "check-solve-gen-linear.h" // check_solve_gen_linear()
#include "check-poly.h" // check_scalars_nonewald_poly()

#include "memory-check.h"


/* main program */
int
main (int argc, char** argv)
{

  check_lapack_inv_ (100, 1);
  check_split_merge (200, 200, 1);
  check_mul_matrices (200, 1);
  check_inverse_by_sub (100, 100, 1);

  check_solve_gen_linear (200, 200, 1);

  check_scalars_nonewald_poly (1, 1.0e-10);

  return 0;
}
