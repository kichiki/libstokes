/* test code for libstokes-core
 * Copyright (C) 2007-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include "check-lapack-solve-lin.h" // check_lapack_solve_lin()
#include "check-mul-matrices.h" // check_mul_matrices()
#include "check-solve-gen-linear.h" // check_solve_gen_linear()

#include "check-twobody.h" // check_twobody_scalars_res_with_equal()
#include "check-poly.h" // check_scalars_nonewald_poly()

#include "check-minv-poly.h" // check_scalars_minv_[f,ft,fts]_poly_...()
#include "check-lub-poly.h" // check_lub_fts_2b_poly()

#include "check-ewald-poly.h" // check_atimes_ewald_3all_poly_SC_...()
#include "check-ewald-shear.h" // check_atimes_3all_ewald_shear()
#include "check-twobody-slip.h" // check_twobody_slip_with_noslip()

#include "check-matrix.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;

  // check-mul-matrices.c
  check += check_mul_matrices (200, 1, 0.0);

  // check-matrix.c
  check += benchmark_dgemm (200, 1, 2.1e-15);
  check += benchmark_mul_matrices (200, 1, 2.0e-15);

  // check-solve-gen-linear.c
  check += check_split_merge (200, 200, 1, 0.0);
  check += check_inverse_by_sub (100, 100, 1, 4.0e-13);
  check += check_solve_gen_linear (200, 200, 1, 1.9e-7);

  // check-lapack-inv.c
  check += check_lapack_inv_ (100, 1, 1.1e-13);
  // check-lapack-solve-lin.c
  check += check_lapack_solve_lin (100, 1, 1.3e-12);
  check += check_lapack_solve_lin (1000, 1, 3.1e-9);


  double r = 2.5;

  double a1 = 1.9;
  double a2 = 0.1;

  // check-twobody.c
  check += check_twobody_lub(2, /* FTS */ r, a1, a2, 100, 1, 4.0e-13);
  check += check_twobody_scalars_res_with_equal(r, 50, 1, 8.0e-7);

  // check-poly.c
  check += check_scalars_nonewald_poly (r, 1, 2.0e-16);
  check += check_scalars_nonewald_poly_symmetry (2.5, a1, a2, 1, 1.0e-17);
  check += check_scalars_ewald_real_poly (2.5, 1.0, 1, 1.3e-15);

  check += check_scalars_minv_f_poly_with_equal   (r, 1, 1.0e-15);
  check += check_scalars_minv_ft_poly_with_equal  (r, 1, 1.0e-15);
  check += check_scalars_minv_fts_poly_with_equal (r, 1, 1.0e-14);

  // check-minv-poly.c
  check += check_scalars_minv_f_poly_ana   (r, a1, a2, 1, 1.0e-15);
  check += check_scalars_minv_ft_poly_ana  (r, a1, a2, 1, 1.0e-12);
  check += check_scalars_minv_fts_poly_ana (r, a1, a2, 1, 1.0e-14);

  // check-lub-poly.c
  check += check_lub_fts_2b_poly (r, 1, 1.0e-4);


  int version = 2; // FTS
  double phi = 0.2;
  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;

  // check-ewald-poly.c
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      4.0, // ewald_tr
					      ewald_eps, 1, 1.1e-14);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      1.0, // ewald_tr
					      ewald_eps, 1, 7.0e-15);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      10.0, // ewald_tr
					      ewald_eps, 1, 2.3e-14);


  check += check_atimes_ewald_3all_poly_SC_2
    (version, 0, // x dir
     phi, ewald_tr, ewald_eps, 1, 8.0e-15);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 1, // y dir
     phi, ewald_tr, ewald_eps, 1, 1.0e-14);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 2, // z dir
     phi, ewald_tr, ewald_eps, 1, 1.0e-14);


  // check-ewald-shear.c
  check += check_atimes_3all_ewald_shear
    (0, // F version
     1, // (x=flow, y=grad)
     0.2, // phi
     1.0, // ewald_tr
     1.0e-12, // ewald_eps,
     1, 6.0e-14);
  check += check_atimes_3all_ewald_shear
    (1, // FT version
     1, // (x=flow, y=grad)
     0.4, // phi
     1.0, // ewald_tr
     1.0e-12, // ewald_eps,
     1, 2.0e-13);
  check += check_atimes_3all_ewald_shear
    (2, // FT version
     2, // (x=flow, z=grad)
     0.5, // phi
     1.0, // ewald_tr
     1.0e-12, // ewald_eps,
     1, 5.0e-12);

  // check-twobody-slip.c
  check += check_twobody_slip_with_noslip (2.5, 1.0,  1, 0.0);
  check += check_twobody_slip_with_noslip (2.5, 10.0, 1, 2.0e-16);


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
