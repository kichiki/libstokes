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
