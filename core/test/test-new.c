/* test code for polydisperse bugs, regarding to
 * non-ewald-new.c and ewald-new.c
 * Copyright (C) 2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include "check-non-ewald-new.h" // check_atimes_nonewald_3all_new ()
#include "check-ewald-new.h" // check_atimes_ewald_3all_new ()

#include "check-ewald-3fts-new_res0.h" // check_solve_res_3fts_0_a ()

#include "check-ewald-3fts-new_res.h" // check_solve_res_3fts_a ()
#include "check-ewald-3fts-new_mob.h" // check_solve_mob_3fts_a ()
#include "check-ewald-3fts-new_mix.h" // check_solve_mix_3fts_a ()

#include "check-ewald-3fts-new_res_lub.h" // check_solve_res_lub_3fts_a ()
#include "check-ewald-3fts-new_mob_lub.h" // check_solve_mob_lub_3fts_a ()
#include "check-ewald-3fts-new_mix_lub.h" // check_solve_mix_lub_3fts_a ()

/* main program */
int
main (int argc, char** argv)
{
  int check = 0;

  int version = 2; // FTS


  // check-non-ewald-new.c
  // compare with monodisperse systems
  check += check_atimes_nonewald_3all_new_0 (version, 1, 0.0);
  // compare with polydisperse systems with a = 1.0 (monodisperse)
  check += check_atimes_nonewald_3all_new_1 (version, 1, 0.0);
  // compare with polydisperse systems scaled by the factor "scale"
  check += check_atimes_nonewald_3all_new_2 (version, 10.0, 1, 1.0e-14);


  // check-ewald-new.c
  // compare with monodisperse systems
  check += check_atimes_ewald_3all_new_0 (version, 1, 0.0);
  // compare with polydisperse systems with a = 1.0 (monodisperse)
  check += check_atimes_ewald_3all_new_1 (version, 1, 0.0);
  // compare with polydisperse systems scaled by the factor "scale"
  check += check_atimes_ewald_3all_new_2 (version, 10.0, 0.25, 1, 3.0e-12);

  // compare mono (without a) with polydisperse systems with a = 1.0
  check += check_atimes_ewald_3all_new_1b (version, 1, 1.3e-14);


  // check-ewald-3fts-new_res0.c
  // compare new (mono) with old (mono)
  check += check_solve_res_3fts_0_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_3fts_0_b (1, 0.0);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_res_3fts_0_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_3fts_0_d (1, 3.9e-11);


  // check-ewald-3fts-new_res.c
  // compare new (mono) with old (mono)
  check += check_solve_res_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_3fts_b (1, 0.0);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_res_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_3fts_d (1, 8.0e-12);

  // check-ewald-3fts-new_mob.c
  // compare new (mono) with old (mono)
  check += check_solve_mob_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mob_3fts_b (1, 0.0);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_mob_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mob_3fts_d (1, 9.2e-14);

  // check-ewald-3fts-new_mix.c
  // compare new (mono) with old (mono)
  check += check_solve_mix_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mix_3fts_b (1, 0.0);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_mix_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mix_3fts_d (1, 1.3e-08);


  // check-ewald-3fts-new_res_lub.c
  // compare new (mono) with old (mono)
  check += check_solve_res_lub_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_lub_3fts_b (1, 1.0e-14);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_res_lub_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_res_lub_3fts_d (1, 2.8e-06);

  // check-ewald-3fts-new_mob_lub.c
  // compare new (mono) with old (mono)
  check += check_solve_mob_lub_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mob_lub_3fts_b (1, 1.2e-14);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_mob_lub_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mob_lub_3fts_d (1, 7.2e-06);

  // check-ewald-3fts-new_mix_lub.c
  // compare new (mono) with old (mono)
  check += check_solve_mix_lub_3fts_a (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mix_lub_3fts_b (1, 6.3e-14);
  // peroidic systems
  // compare new (mono) with old (mono)
  check += check_solve_mix_lub_3fts_c (1, 0.0);
  // compare with new (mono) and new (poly a=1)
  check += check_solve_mix_lub_3fts_d (1, 7.6e-06);


  // check-ewald-3fts-new_res.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_res_3fts_b2 (10.0, 1, 8.7e-13);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_res_3fts_d2 (10.0, 1, 3.6e-08);

  // check-ewald-3fts-new_mob.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_mob_3fts_b2 (10.0, 1, 8.6e-13);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_mob_3fts_d2 (10.0, 1, 3.3e-11);

  // check-ewald-3fts-new_mix.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_mix_3fts_b2 (10.0, 1, 2.5e-08);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_mix_3fts_d2 (10.0, 1, 2.8e-07);

  // check-ewald-3fts-new_res_lub.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_res_lub_3fts_b2 (10.0, 1, 7.5e-13);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_res_lub_3fts_d2 (10.0, 1, 4.2e-08);

  // check-ewald-3fts-new_mob_lub.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_mob_lub_3fts_b2 (10.0, 1, 1.4e-12);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_mob_lub_3fts_d2 (10.0, 1, 1.7e-06);

  // check-ewald-3fts-new_mix_lub.c PART II
  // compare with mono and scaled poly (a=10)
  check += check_solve_mix_lub_3fts_b2 (10.0, 1, 2.5e-06);
  // ** peroidic systems **
  // compare with mono and scaled poly (a=10)
  check += check_solve_mix_lub_3fts_d2 (10.0, 1, 1.9e-05);


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
