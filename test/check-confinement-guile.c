/* test code for confinement-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-confinement-guile.c,v 1.1 2008/05/24 06:02:54 kichiki Exp $
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

#include "stokes-guile.h" // guile_load()
#include "confinement-guile.h"
#include "confinement.h"


/* check reading SCM script
 */
int
check_CF_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_CF_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;

  struct confinement *cf = NULL;

  char *file0 = "check-confinement-guile-sphere.scm";
  char *file1 = "check-confinement-guile-sphere-hole.scm";
  char *file2 = "check-confinement-guile-cylinder.scm";
  char *file3 = "check-confinement-guile-dumbbell.scm";
  char *file4 = "check-confinement-guile-hex2d.scm";
  char *file5 = "check-confinement-guile-porous.scm";

  // sphere
  guile_load (file0);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " sphere : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " sphere : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " sphere : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 0.0, " sphere : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 10.0, " sphere : R", verbose, tiny, &max);
  check += compare_max (cf->r, 0.0, " sphere : r", verbose, tiny, &max);
  check += compare_max (cf->x, 0.0, " sphere : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " sphere : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " sphere : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 0.0, " sphere : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 0.0, " sphere : L", verbose, tiny, &max);

  CF_free (cf);

  // sphere+hole
  guile_load (file1);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " sphere+hole : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " sphere+hole : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " sphere+hole : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 1.0, " sphere+hole : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 10.0, " sphere+hole : R", verbose, tiny, &max);
  check += compare_max (cf->r, 1.0, " sphere+hole : r", verbose, tiny, &max);
  check += compare_max (cf->x, 0.0, " sphere+hole : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " sphere+hole : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " sphere+hole : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 0.0, " sphere+hole : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 0.0, " sphere+hole : L", verbose, tiny, &max);

  CF_free (cf);

  // cylinder
  guile_load (file2);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " cylinder : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " cylinder : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " cylinder : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 2.0, " cylinder : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 0.0, " cylinder : R", verbose, tiny, &max);
  check += compare_max (cf->r, 10.0, " cylinder : r", verbose, tiny, &max);
  check += compare_max (cf->x, 1.0, " cylinder : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " cylinder : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " cylinder : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 0.0, " cylinder : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 0.0, " cylinder : L", verbose, tiny, &max);

  CF_free (cf);

  // dumbbell
  guile_load (file3);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " dumbbell : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " dumbbell : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " dumbbell : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 3.0, " dumbbell : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 10.0, " dumbbell : R", verbose, tiny, &max);
  check += compare_max (cf->r, 1.0, " dumbbell : r", verbose, tiny, &max);
  check += compare_max (cf->x, 0.0, " dumbbell : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " dumbbell : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " dumbbell : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 10.0, " dumbbell : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 2.0, " dumbbell : L", verbose, tiny, &max);

  CF_free (cf);

  // hex2d
  guile_load (file4);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " hex2d : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " hex2d : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " hex2d : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 4.0, " hex2d : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 10.0, " hex2d : R", verbose, tiny, &max);
  check += compare_max (cf->r, 1.0, " hex2d : r", verbose, tiny, &max);
  check += compare_max (cf->x, 0.0, " hex2d : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " hex2d : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " hex2d : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 0.0, " hex2d : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 12.0, " hex2d : L", verbose, tiny, &max);

  CF_free (cf);

  // porous
  guile_load (file5);
  cf = CF_guile_get ("confinement");
  CHECK_MALLOC (cf, "check_CF_guile_get");

  // check the result
  check += compare_max ((double)cf->flag_LJ, 0.0, " porous : LJ-flag",
			verbose, tiny, &max);
  check += compare_max (cf->e, 10.0, " porous : LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, 1.0, " porous : LJ-r0", verbose, tiny, &max);
  check += compare_max ((double)cf->type, 5.0, " porous : type",
			verbose, tiny, &max);
  check += compare_max (cf->R, 10.0, " porous : R", verbose, tiny, &max);
  check += compare_max (cf->r, 0.0, " porous : r", verbose, tiny, &max);
  check += compare_max (cf->x, 0.0, " porous : x", verbose, tiny, &max);
  check += compare_max (cf->y, 0.0, " porous : y", verbose, tiny, &max);
  check += compare_max (cf->z, 0.0, " porous : z", verbose, tiny, &max);
  check += compare_max (cf->R2, 0.0, " porous : R2", verbose, tiny, &max);
  check += compare_max (cf->L, 20.0, " porous : L", verbose, tiny, &max);

  CF_free (cf);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
