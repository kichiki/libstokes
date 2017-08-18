/* test code for libstokes-extra
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

// Brownian dynamics stuff
#include "check-KIrand.h"
#include "check-chebyshev.h"
#include "check-dsaupd_c.h"
#include "check-dnaupd_c.h"
#include "check-bd.h"
#include "check-dgeev_c.h"

#include "check-ode-quaternion.h" // check_quaternion_Winv()

#include "check-dpotrf_c.h"
#include "check-brownian.h"
#include "check-bd-imp.h"
#include "check-bd-imp-nitsol.h"
#include "check-bd-imp-fast.h"
#include "check-sqrt-dgeev.h"

//#include "check-list-ex.h"
#include "check-bonds-guile.h"
#include "check-angles.h"
#include "check-excluded-volume-guile.h"
#include "check-ev-dh.h"
#include "check-ev-dh-guile.h"
#include "check-ev-LJ.h"
#include "check-ev-LJ-guile.h"

#include "check-confinement.h"
#include "check-confinement-guile.h"

#include "check-bead-rod.h"
#include "check-bead-rod-guile.h"

#include "check-solve-cubic.h"

#include "check-grid.h"
#include "check-ev-dh-grid.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;

  //check_BeadRod_calc_dr (10, 1.0e-8, 1, 0.0);
  //exit (1);


  // Brownian dynamics stuff
  check += check_chebyshev (1, 1e-10);
  check += check_dsaupd_c (100, 1, 3.1e-15);
  check += check_dnaupd_c (10, 1, 2.0e-15);
  check += check_dnaupd_c (100, 1, 1.2e-15);
  check += check_bd (50, 100, 2.005, 1, 8.0e-10);
  check += check_bd (100, 100, 2.005, 1, 2.0e-15);
  check += check_dgeev_c (1, 6.0e-15);

  // check-ode-quaternion.c
  check += check_quaternion_Winv (1, 3.0e-16);

  // check-dpotrf_c.c
  check += check_dpotrf_c (1000, 1, 5.2e-15);
  check += check_dpotf2_c (1000, 1, 5.0e-15);

  // check-brownian.c
  check += check_cheb_minv (10, 1, 2.0e-14);
  check += check_cheb_lub  (10, 1, 6.4e-15);

  // check-brownian.c -- serious tests for calc_brownian_force()
  check += check_minv_FU (1, 1.6e-10);
  check += check_lub_FU (1, 2.6e-13);
  check += benchmark_BD_minv_FU_in_FTS (50, 1, 8.0e-9); // a bit large...
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

  // check-bd-imp-fast.c
  // noHI
  check += check_fastSI_rhs (0, 3, 1, 1, 3.0e-16); // Fraenkel
  check += check_fastSI_rhs (1, 3, 1, 1, 0.0);     // WLC
  check += check_fastSI_rhs (3, 3, 1, 1, 0.0);     // Cohen
  check += check_fastSI_rhs (4, 3, 1, 1, 0.0);     // Werner
  check += check_fastSI_rhs (5, 3, 1, 1, 0.0); // Hook
  check += check_fastSI_rhs (6, 3, 1, 1, 3.0e-16); // Fraenkel for dWLC
  check += check_fastSI_rhs (7, 3, 1, 1, 3.0e-16); // FENE-Fraenkel

  // HI
  check += check_fastSI_rhs (0, 3, 0, 1, 0.0);     // Fraenkel
  check += check_fastSI_rhs (1, 3, 0, 1, 2.0e-16); // WLC
  check += check_fastSI_rhs (3, 3, 0, 1, 2.0e-16); // Cohen
  check += check_fastSI_rhs (4, 3, 0, 1, 2.0e-16); // Werner
  check += check_fastSI_rhs (5, 3, 0, 1, 2.0e-16); // Hook
  check += check_fastSI_rhs (6, 3, 0, 1, 0.0);     // Fraenkel for dWLC
  check += check_fastSI_rhs (7, 3, 0, 1, 0.0);     // FENE-Fraenkel


  // noHI
  check += check_fastSI_solve_cubic (0, 3, 1, 1, 3.6e-16); // Fraenkel
  check += check_fastSI_solve_cubic (1, 3, 1, 1, 2.0e-15); // WLC
  check += check_fastSI_solve_cubic (3, 3, 1, 1, 3.0e-15); // Cohen
  check += check_fastSI_solve_cubic (4, 3, 1, 1, 7.5e-16); // Werner
  check += check_fastSI_solve_cubic (5, 3, 1, 1, 2.1e-16); // Hook
  check += check_fastSI_solve_cubic (6, 3, 1, 1, 3.6e-16); // Fraenkel for dWLC
  check += check_fastSI_solve_cubic (7, 3, 1, 1, 5.0e-12); // FENE-Fraenkel

  // HI
  check += check_fastSI_solve_cubic (0, 3, 0, 1, 2.2e-16); // Fraenkel
  check += check_fastSI_solve_cubic (1, 3, 0, 1, 1.6e-15); // WLC
  check += check_fastSI_solve_cubic (3, 3, 0, 1, 4.0e-16); // Cohen
  check += check_fastSI_solve_cubic (4, 3, 0, 1, 2.0e-15); // Werner
  check += check_fastSI_solve_cubic (5, 3, 0, 1, 3.0e-16); // Hook
  check += check_fastSI_solve_cubic (6, 3, 0, 1, 2.2e-16); // Fraenkel for dWLC
  check += check_fastSI_solve_cubic (7, 3, 0, 1, 4.4e-12); // FENE-Fraenkel

  // noHI
  check += check_fastSI_solve (0, 3, 1, 1.0e-2, 1, 3.0e-9);  // Fraenkel
  check += check_fastSI_solve (1, 3, 1, 1.0e-2, 1, 1.0e-9);  // WLC
  check += check_fastSI_solve (3, 3, 1, 1.0e-2, 1, 2.0e-10); // Cohen
  check += check_fastSI_solve (4, 3, 1, 1.0e-2, 1, 4.0e-9);  // Werner
  check += check_fastSI_solve (5, 3, 1, 1.0e-2, 1, 5.0e-10); // Hook
  check += check_fastSI_solve (6, 3, 1, 1.0e-2, 1, 3.0e-9);  // Fraenkel
  check += check_fastSI_solve (7, 3, 1, 1.0e-2, 1, 5.0e-9);  // FENE-Fraenkel

  // noHI
  check += check_fastSI_solve (0, 10, 1, 1.0e-2, 1, 2.0e-9);  // Fraenkel
  check += check_fastSI_solve (1, 10, 1, 1.0e-2, 1, 2.0e-9);  // WLC
  check += check_fastSI_solve (3, 10, 1, 1.0e-2, 1, 8.0e-10); // Cohen
  check += check_fastSI_solve (4, 10, 1, 1.0e-2, 1, 3.0e-9);  // Werner
  check += check_fastSI_solve (5, 10, 1, 1.0e-2, 1, 3.0e-9);  // Hook
  check += check_fastSI_solve (6, 10, 1, 1.0e-2, 1, 2.0e-9);  // Fraenkel
  check += check_fastSI_solve (7, 10, 1, 1.0e-2, 1, 3.0e-9);  // FENE-Fraenkel

  // noHI
  check += check_fastSI_solve (0, 100, 1, 1.0e-2, 1, 5.0e-9);  // Fraenkel
  check += check_fastSI_solve (1, 100, 1, 1.0e-2, 1, 3.0e-10); // WLC
  check += check_fastSI_solve (3, 100, 1, 1.0e-2, 1, 2.0e-9);  // Cohen
  check += check_fastSI_solve (4, 100, 1, 1.0e-2, 1, 8.0e-10); // Werner
  check += check_fastSI_solve (5, 100, 1, 1.0e-2, 1, 3.0e-10); // Hook
  check += check_fastSI_solve (6, 100, 1, 1.0e-2, 1, 5.0e-9);  // Fraenkel
  check += check_fastSI_solve (7, 100, 1, 1.0e-2, 1, 2.0e-9);  // FENE-Fraenkel

  // HI -- dt is reduced
  check += check_fastSI_solve (0, 3, 0, 1.0e-3, 1, 7.0e-9);  // Fraenkel
  check += check_fastSI_solve (1, 3, 0, 1.0e-3, 1, 3.0e-10); // WLC
  check += check_fastSI_solve (3, 3, 0, 1.0e-3, 1, 5.0e-10); // Cohen
  check += check_fastSI_solve (4, 3, 0, 1.0e-3, 1, 3.0e-10); // Werner
  check += check_fastSI_solve (5, 3, 0, 1.0e-3, 1, 2.0e-9);  // Hook
  check += check_fastSI_solve (6, 3, 0, 1.0e-3, 1, 7.0e-9);  // Fraenkel
  check += check_fastSI_solve (7, 3, 0, 1.0e-3, 1, 4.0e-9);  // FENE-Fraenkel

  // some tweaks for dt are needed...
  check += check_fastSI_solve (0, 10, 0, 1.0e-3, 1, 7.0e-10); // Fraenkel
  check += check_fastSI_solve (1, 10, 0, 1.0e-4, 1, 3.0e-10); // WLC
  check += check_fastSI_solve (3, 10, 0, 1.0e-4, 1, 5.0e-11); // Cohen
  check += check_fastSI_solve (4, 10, 0, 1.0e-4, 1, 2.0e-10); // Werner
  check += check_fastSI_solve (5, 10, 0, 1.0e-3, 1, 7.0e-9);  // Hook
  check += check_fastSI_solve (6, 10, 0, 1.0e-3, 1, 7.0e-10); // Fraenkel
  check += check_fastSI_solve (7, 10, 0, 1.0e-3, 1, 3.0e-9);  // FENE-Fraenkel

  // further tweaks for dt are needed...
  check += check_fastSI_solve (0, 100, 0, 1.0e-3, 1, 9.0e-11); // Fraenkel
  check += check_fastSI_solve (1, 100, 0, 1.0e-6, 1, 4.0e-9);  // WLC
  check += check_fastSI_solve (3, 100, 0, 1.0e-6, 1, 4.0e-9);  // Cohen
  check += check_fastSI_solve (4, 100, 0, 1.0e-6, 1, 2.0e-9);  // Werner
  check += check_fastSI_solve (5, 100, 0, 1.0e-6, 1, 9.0e-11); // Hook
  check += check_fastSI_solve (6, 100, 0, 1.0e-3, 1, 9.0e-11); // Fraenkel
  check += check_fastSI_solve (7, 100, 0, 1.0e-3, 1, 3.0e-9);  // FENE-Fraenkel


  // check-sqrt-dgeev.c
  check += check_BD_sqrt_by_dgeev (100, 1, 4.1e-13);

  // check-KIrand.c
  check += check_KIrand_Gaussian (1, 5.0e-4);
  check += check_KIrand_Gaussian_cont (1, 0.0);


  // check-list-ex.c
  //check += check_list_ex (1, 0.0);

  // check-bonds-guile.c
  check += check_bonds_guile_get (1, 0.0);

  // check-excluded-volume-guile.c
  check += check_EV_guile_get (1, 1.0e-15);

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
  check += check_CF_sphere_calc_force (10.0, 1, 7.8e-14);
  check += check_CF_sphere_calc_force (100.0, 1, 2.0e-13);
  // check-confinement-guile.c
  check += check_CF_guile_get (1, 0.0);


  // check-bead-rod.c
  check += check_BeadRod_constraint_displacement (10, 1, 2.0e-15);
  check += check_BeadRod_solve_iter_gamma (10, 1.0e-8, 1, 3.0e-9);
  check += check_BeadRod_solve_gamma_by_NITSOL (10, 1.0e-8, 1, 5.0e-9);
  check += check_BeadRod_solve_gamma (10, 1.0e-8, 1, 2.0e-9);
  // check-bead-rod-guile.c
  check += check_BeadRod_guile_get (1, 0.0);


  // check-solve-cubic.c
  check += check_solve_cubic (1, 3.0e-16);


  // check-grid.c
  check += check_GRID_ixyz_to_in_to_ixyz (1, 0.0);
  // randon config
  check += check_GRID_init_all_by_cutoff (100, 0, 1, 0.0);
  check += check_GRID_init_all_by_cutoff (1000, 0, 1, 0.0);
  check += check_GRID_init_all_by_cutoff (10000, 0, 1, 0.0);
  // regular array config
  check += check_GRID_init_all_by_cutoff (100, 1, 1, 0.0);
  check += check_GRID_init_all_by_cutoff (1000, 1, 1, 0.0);
  check += check_GRID_init_all_by_cutoff (10000, 1, 1, 0.0);

  // check-ev-dh-grid.c
  check += check_EV_DH_calc_force_grid (100, 1, 7.0e-16);
  check += check_EV_DH_calc_force_grid (1000, 1, 2.0e-15);
  check += check_EV_DH_calc_force_grid (10000, 1, 2.0e-15);


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
