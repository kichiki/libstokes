/* fragment of test/test-all.c for obsolete codes
 * test code for libstokes-core
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

#include "check-minv.h" // check_matrix_mob_nonewald_fts(), check_minv_fts()
#include "check-lub-poly-old.h" // check_lub_fts_2b_poly()
#include "check-mob-fts.h" // check_mob_fts ()
#include "check-ewald.h" // check_atimes_ewald_3all_SC ()
#include "check-ewald-poly-old.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;

  double r = 2.5;

  double a1 = 1.9;
  double a2 = 0.1;

  // check-minv.c
  check += check_matrix_mob_nonewald_fts (r, 1, 1.0e-14);
  check += check_minv_fts (2.5, 1, 1.0e-14);

  // check-lub-poly-old.c
  // the reason of such big tiny is because of the accuracy of two-body-res.c
  check += check_lub_scalars_poly (r, 50, 1, 3.0e-6);
  check += check_matrix_lub_fts_2b_poly (r, 1, 1.0e-5);
  check += check_atimes_matrix_lub_fts_2b_poly (r, a1, a2, 1, 1.0e-15);

  // check-mob-fts.c
  check += check_mob_fts (1, 7.6e-12);
  check += check_mob_lub_fts (1, 3.4e-11);


  int version = 2; // FTS
  double phi = 0.2;
  double ewald_tr = 1.0;
  double ewald_eps = 1.0e-12;
  // check-ewald.c
  check += check_ewald_3all_atimes_matrix_SC
    (version, phi, ewald_tr, ewald_eps, 1, 5.0e-14);

  // check-ewald-poly.c

  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 4.0, // ewald_tr
     ewald_eps, 1, 7.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 1.0, // ewald_tr
     ewald_eps, 1, 8.0e-15);
  check += check_make_matrix_mob_ewald_3all_poly_SC_1
    (version, phi, 10.0, // ewald_tr
     ewald_eps, 1, 2.0e-14);

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
