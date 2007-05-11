/* test code for libstokes
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-all.c,v 1.8 2007/05/11 02:08:10 kichiki Exp $
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

#include "check-lapack-inv.h" // check_lapack_inv_()
#include "check-mul-matrices.h" // check_mul_matrices()
#include "check-solve-gen-linear.h" // check_solve_gen_linear()
#include "check-twobody.h" // check_twobody_scalars_res_with_equal()
#include "check-poly.h" // check_scalars_nonewald_poly()
#include "check-minv.h" // check_matrix_mob_nonewald_fts(), check_minv_fts()
#include "check-minv-poly.h" // check_scalars_minv_[f,ft,fts]_poly_...()
#include "check-lub-poly.h" // check_lub_fts_2b_poly()
#include "check-ewald.h" // check_atimes_ewald_3all_SC ()
#include "check-ewald-poly.h" // check_atimes_ewald_3all_poly_SC_...()

#include "memory-check.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;
  check += check_lapack_inv_ (100, 1);
  check += check_split_merge (200, 200, 1);
  check += check_mul_matrices (200, 1);
  check += check_inverse_by_sub (100, 100, 1);

  check += check_solve_gen_linear (200, 200, 1);


  double r = 2.5;

  double a1 = 1.9;
  double a2 = 0.1;

  // check-twobody.c
  check += check_twobody_lub(2, /* FTS */ r, a1, a2, 100, 1, 2.0e-16);
  check += check_twobody_scalars_res_with_equal(r, 50, 1, 8.0e-7);

  // check-poly.c
  check += check_scalars_nonewald_poly (r, 1, 1.0e-17);
  check += check_scalars_nonewald_poly_symmetry (2.5, a1, a2, 1, 1.0e-17);
  check += check_scalars_ewald_real_poly (2.5, 1.0, 1, 1.0e-15);

  check += check_matrix_mob_nonewald_fts (r, 1, 1.0e-14);
  check += check_minv_fts (2.5, 1, 1.0e-14);

  check += check_scalars_minv_f_poly_with_equal   (r, 1, 1.0e-15);
  check += check_scalars_minv_ft_poly_with_equal  (r, 1, 1.0e-15);
  check += check_scalars_minv_fts_poly_with_equal (r, 1, 1.0e-14);

  // check-minv-poly.c
  check += check_scalars_minv_f_poly_ana   (r, a1, a2, 1, 1.0e-15);
  check += check_scalars_minv_ft_poly_ana  (r, a1, a2, 1, 1.0e-12);
  check += check_scalars_minv_fts_poly_ana (r, a1, a2, 1, 1.0e-14);

  // check-lub-poly.c
  // the reason of such big tiny is because of the accuracy of two-body-res.c
  check += check_lub_scalars_poly (r, 50, 1, 3.0e-6);
  check += check_lub_fts_2b_poly (r, 1, 1.0e-4);
  check += check_matrix_lub_fts_2b_poly (r, 1, 1.0e-5);
  check += check_atimes_matrix_lub_fts_2b_poly (r, a1, a2, 1, 1.0e-15);

  // check-ewald.c
  int version = 2; // FTS
  double phi = 0.2;
  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;
  check += check_ewald_3all_atimes_matrix_SC
    (version, phi, ewald_tr, ewald_eps, 1, 5.0e-14);

  // check-ewald-poly.c
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      4.0, // ewald_tr
				     ewald_eps, 1, 7.0e-15);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      1.0, // ewald_tr
					      ewald_eps, 1, 4.0e-15);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      10.0, // ewald_tr
					      ewald_eps, 1, 9.0e-15);

  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 4.0, // ewald_tr
     ewald_eps, 1, 2.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 1.0, // ewald_tr
     ewald_eps, 1, 8.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 10.0, // ewald_tr
     ewald_eps, 1, 6.0e-15);

  check += check_atimes_ewald_3all_poly_SC_2
    (version, 0, // x dir
     phi, ewald_tr, ewald_eps, 1, 4.0e-15);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 1, // y dir
     phi, ewald_tr, ewald_eps, 1, 5.0e-15);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 2, // z dir
     phi, ewald_tr, ewald_eps, 1, 4.0e-15);

  check += check_make_matrix_mob_ewald_3all_poly_SC_2
    (version, 0, // x dir
     phi, ewald_tr, ewald_eps,
     1, 2.0e-14);
  check += check_make_matrix_mob_ewald_3all_poly_SC_2
    (version, 1, // y dir
     phi, ewald_tr, ewald_eps,
     1, 2.0e-14);
  check += check_make_matrix_mob_ewald_3all_poly_SC_2
    (version, 2, // z dir
     phi, ewald_tr, ewald_eps,
     1, 1.0e-14);

  fprintf (stdout,
	   "==================================================\n"
	   "TOTAL ERRORS : %d\n", check);
  if (check == 0)
    {
      fprintf (stdout, "Conglaturation!! ALL TESTS PASSED\n");
    }
  else
    {
      fprintf (stdout, "sorry. some test(s) FAILED\n");
    }

  return 0;
}
