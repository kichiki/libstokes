/* test code for libstokes
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-all.c,v 1.27 2008/06/13 03:13:17 kichiki Exp $
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
#include "check-minv.h" // check_matrix_mob_nonewald_fts(), check_minv_fts()
#include "check-minv-poly.h" // check_scalars_minv_[f,ft,fts]_poly_...()
#include "check-lub-poly.h" // check_lub_fts_2b_poly()
#include "check-mob-fts.h" // check_mob_fts ()
#include "check-ewald.h" // check_atimes_ewald_3all_SC ()
#include "check-ewald-poly.h" // check_atimes_ewald_3all_poly_SC_...()
#include "check-ewald-shear.h" // check_atimes_3all_ewald_shear()
#include "check-ode-quaternion.h" // check_quaternion_Winv()
#include "check-twobody-slip.h" // check_twobody_slip_with_noslip()

// Brownian dynamics stuff
#include "check-KIrand.h"
#include "check-chebyshev.h"
#include "check-dsaupd_c.h"
#include "check-dnaupd_c.h"
#include "check-bd.h"
#include "check-dgeev_c.h"

#include "check-dpotrf_c.h"
#include "check-brownian.h"
#include "check-bd-imp.h"
#include "check-bd-imp-nitsol.h"
#include "check-sqrt-dgeev.h"

#include "check-matrix.h"

#include "check-list-ex.h"
#include "check-angles.h"
#include "check-ev-dh.h"
#include "check-ev-dh-guile.h"
#include "check-ev-LJ.h"
#include "check-ev-LJ-guile.h"

#include "check-confinement.h"
#include "check-confinement-guile.h"

#include "check-bead-rod.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;


  // Brownian dynamics stuff
  check += check_chebyshev (1, 1e-10);
  check += check_dsaupd_c (100, 1, 3.1e-15);
  check += check_dnaupd_c (10, 1, 1.5e-15);
  check += check_dnaupd_c (100, 1, 1.2e-15);
  check += check_bd (50, 100, 2.005, 1, 8.0e-10);
  check += check_bd (100, 100, 2.005, 1, 2.0e-15);
  check += check_dgeev_c (1, 6.0e-15);

  // check-mul-matrices.c
  check += check_mul_matrices (200, 1, 0.0);

  // check-matrix.c
  check += benchmark_dgemm (200, 1, 2.1e-15);
  check += benchmark_mul_matrices (200, 1, 0.0);

  // check-solve-gen-linear.c
  check += check_split_merge (200, 200, 1, 0.0);
  check += check_inverse_by_sub (100, 100, 1, 2.9e-13);
  check += check_solve_gen_linear (200, 200, 1, 1.8e-7);

  // check-lapack-inv.c
  check += check_lapack_inv_ (100, 1, 7.1e-14);
  // check-lapack-solve-lin.c
  check += check_lapack_solve_lin (100, 1, 1.3e-12);
  check += check_lapack_solve_lin (1000, 1, 5.0e-10);


  double r = 2.5;

  double a1 = 1.9;
  double a2 = 0.1;

  // check-twobody.c
  check += check_twobody_lub(2, /* FTS */ r, a1, a2, 100, 1, 4.0e-13);
  check += check_twobody_scalars_res_with_equal(r, 50, 1, 8.0e-7);

  // check-poly.c
  check += check_scalars_nonewald_poly (r, 1, 2.0e-16);
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

  // check-mob-fts.c
  check += check_mob_fts (1, 7.0e-12);
  check += check_mob_lub_fts (1, 3.3e-11);

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
					      ewald_eps, 1, 8.5e-15);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      1.0, // ewald_tr
					      ewald_eps, 1, 4.0e-15);
  check += check_atimes_ewald_3all_poly_SC_1 (version, phi,
					      10.0, // ewald_tr
					      ewald_eps, 1, 2.3e-14);

  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 4.0, // ewald_tr
     ewald_eps, 1, 7.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 1.0, // ewald_tr
     ewald_eps, 1, 8.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 10.0, // ewald_tr
     ewald_eps, 1, 2.0e-14);

  check += check_atimes_ewald_3all_poly_SC_2
    (version, 0, // x dir
     phi, ewald_tr, ewald_eps, 1, 8.0e-15);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 1, // y dir
     phi, ewald_tr, ewald_eps, 1, 1.0e-14);
  check += check_atimes_ewald_3all_poly_SC_2
    (version, 2, // z dir
     phi, ewald_tr, ewald_eps, 1, 1.0e-14);

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
     1, 2.0e-14);

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

  // check-ode-quaternion.c
  check += check_quaternion_Winv (1, 3.0e-16);

  // check-twobody-slip.c
  check += check_twobody_slip_with_noslip (2.5, 1.0,  1, 0.0);
  check += check_twobody_slip_with_noslip (2.5, 10.0, 1, 2.0e-16);

  // check-dpotrf_c.c
  check += check_dpotrf_c (1000, 1, 5.2e-15);
  check += check_dpotf2_c (1000, 1, 5.0e-15);

  // check-brownian.c
  check += check_cheb_minv (10, 1, 2.0e-14);
  check += check_cheb_lub  (10, 1, 6.4e-15);

  // check-brownian.c -- serious tests for calc_brownian_force()
  check += check_minv_FU (1, 1.6e-10);
  check += check_lub_FU (1, 2.6e-13);
  check += benchmark_BD_minv_FU_in_FTS (50, 1, 6.0e-9); // a bit large...
  check += check_inv_by_submatrices (50, 50, 1, 1.4e-10);

  // check-bd-imp.c
  // for BD_evolve_JGdP00()
  // dt dependence
  check += check_BD_evolve_JGdP00 (0,   // F version
				   0,   // flag_lub
				   0,   // flag_mat
				   0,   // flag_Q
				   0.1, // dt
				   1, 8.0e-6);
  check += check_BD_evolve_JGdP00 (0,    // F version
				   0,   // flag_lub
				   0,   // flag_mat
				   0,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);
  check += check_BD_evolve_JGdP00 (0,    // F version
				   0,   // flag_lub
				   0,   // flag_mat
				   0,    // flag_Q
				   0.001,// dt
				   1, 8.0e-10);
  // higher versions
  check += check_BD_evolve_JGdP00 (1,    // FT version
				   0,   // flag_lub
				   0,   // flag_mat
				   0,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);
  check += check_BD_evolve_JGdP00 (2,    // FTS version
				   0,   // flag_lub
				   0,   // flag_mat
				   0,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);
  // lub
  check += check_BD_evolve_JGdP00 (2,    // FTS version
				   1,    // flag_lub
				   1,    // flag_mat
				   0,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);
  // quaternion
  check += check_BD_evolve_JGdP00 (1,    // FTS version
				   0,    // flag_lub
				   1,    // flag_mat
				   1,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);
  check += check_BD_evolve_JGdP00 (2,    // FT version
				   1,    // flag_lub
				   1,    // flag_mat
				   1,    // flag_Q
				   0.01, // dt
				   1, 8.0e-8);

  // for BD_imp_ode_evolve()
  // t_out dependence
  check += check_BD_imp_ode_evolve (0,   // F version
				    0,   // flag_lub
				    0,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);
  check += check_BD_imp_ode_evolve (0,   // F version
				    0,   // flag_lub
				    0,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    10.0,// t_out
				    1, 9.0e-5);
  check += check_BD_imp_ode_evolve (0,   // F version
				    0,   // flag_lub
				    0,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    100.0,// t_out
				    1, 3.0e-3);
  // higher versions
  check += check_BD_imp_ode_evolve (1,   // FT version
				    0,   // flag_lub
				    0,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);
  check += check_BD_imp_ode_evolve (2,   // FTS version
				    0,   // flag_lub
				    0,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);
  // lub
  check += check_BD_imp_ode_evolve (2,   // FTS version
				    1,   // flag_lub
				    1,   // flag_mat
				    0,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);
  // quaternion
  check += check_BD_imp_ode_evolve (1,   // FT version
				    0,   // flag_lub
				    1,   // flag_mat
				    1,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);
  check += check_BD_imp_ode_evolve (2,   // FTS version
				    1,   // flag_lub
				    1,   // flag_mat
				    1,   // flag_Q
				    0.01,// h
				    1.0, // t_out
				    1, 8.0e-6);

  // check-bd-imp-nitsol.c
  check += check_BD_imp_NITSOL (0,   // F version
				20,  // N
				0,   // flag_lub
				0,   // flag_mat
				0,   // flag_Q
				0.1, // dt
				1, 3.1e-4);

  check += check_BD_imp_NITSOL (1,   // FT version
				20,  // N
				0,   // flag_lub
				0,   // flag_mat
				1,   // flag_Q
				0.1, // dt
				1, 3.1e-4);

  check += check_BD_imp_NITSOL (2,   // FTS version
				20,  // N
				0,   // flag_lub
				0,   // flag_mat
				1,   // flag_Q
				0.1, // dt
				1, 3.1e-4);


  // check-sqrt-dgeev.c
  check += check_BD_sqrt_by_dgeev (100, 1, 4.1e-13);

  // check-KIrand.c
  check += check_KIrand_Gaussian (1, 5.0e-4);
  check += check_KIrand_Gaussian_cont (1, 0.0);


  // check-list-ex.c
  check += check_list_ex (1, 0.0);

  // check-angles.c
  check += check_angles_guile_get (1, 0.0);
  check += check_angles_calc_force (1.0, 180.0, 1, 2.0e-13);
  check += check_angles_calc_force (2.0, 170.0, 1, 2.0e-13);
  check += check_angles_calc_force (3.0, 90.0,  1, 2.0e-13);

  // check-ev-dh.c
  check += check_EV_DH_calc_force (3.0, 1, 0.0);
  // check-ev-dh-guile.c
  check += check_EV_DH_guile_get (1, 0.0);

  // check-ev-LJ.c
  check += check_EV_LJ_calc_force (3.0, 1, 0.0);
  // check-ev-LJ-guile.c
  check += check_EV_LJ_guile_get (1, 0.0);


  // check-confinement.c
  check += check_CF_init (1, 1.0e-15);
  check += check_CF_sphere_calc_force (10.0, 1, 4.0e-14);
  check += check_CF_sphere_calc_force (100.0, 1, 2.0e-13);
  // check-confinement-guile.c
  check += check_CF_guile_get (1, 0.0);


  // check-bead-rod.c
  check += check_BeadRod_constraint_displacement (10, 1, 2.0e-15);
  check += check_BeadRod_solve_iter_gamma (10, 1.0e-8, 1, 4.0e-10);
  check += check_BeadRod_solve_gamma_by_NITSOL (10, 1.0e-8, 1, 4.0e-10);


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
