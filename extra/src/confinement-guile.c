/* guile interface for struct confinement
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <libguile.h>
#include <string.h> // strcmp()
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes-guile.h" // guile_load()
#include "confinement.h"

#include "confinement-guile.h"


/* get confinement from SCM
 * in SCM, confinement is given by one of these:
 * for spherical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *  ))
 * for spherical confinement with a hole,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "sphere+hole"
 *    10.0 ;; radius of the cavity at (0, 0, 0)
 *    1.0  ;; radius of the hole at (0, 0, 1) direction
 *  ))
 * for cylindrical confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "cylinder"    ;; the cylinder center goes through (0, 0, 0) and (x, y, z).
 *    10.0          ;; radius of the cylinder
 *    1.0  0.0  0.0 ;; direction vector (x, y, z) of the cylinder
 *  ))
 * for dumbbell confinement,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "dumbbell" ;; the origin is at the center of the cylinder
 *    10.0       ;; left cavity radius centered at (center1, 0, 0)
 *    10.0       ;; right cavity radius centered at (center2, 0, 0)
 *    2.0        ;; length of the cylinder
 *    1.0        ;; cylinder radius
 *  ))
 * for 2D hexagonal confinement with cylinder pipe,
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "hex2d"
 *    10.0    ;; cavity radius
 *    1.0     ;; cylinder radius
 *    12.0    ;; lattice spacing
 *  ))
 * for porous media (outside of the 3D hexagonal particle array)
 *  (define confinement '(
 *    10.0 ;; LJ parameter epsilon in kT (so this is dimensionless value)
 *    1.0  ;; LJ parameter r0 in "length" (so this is dimensionless value)
 *    "porous"
 *    10.0    ;; particle radius
 *    20.0    ;; lattice spacing in x (2R for touching case)
 *  ))
 * INPUT
 *  var : name of the variable.
 *        in the above example, set "confinement".
 * OUTPUT
 *  returned value : struct confinement
 *                   if NULL is returned, it failed (not defined)
 */
struct confinement *
CF_guile_get (const char *var)
{
  if (guile_check_symbol (var) == 0)
    {
      fprintf (stderr, "CF_guile_get: %s is not defined\n", var);
      return (NULL);
    }

  SCM scm_symbol
    = scm_c_lookup (var);

  SCM scm_confinement
    = scm_variable_ref (scm_symbol);

  if (!SCM_NFALSEP (scm_list_p (scm_confinement)))
    {
      fprintf (stderr, "CF_guile_get: %s is not a list\n", var);
      return (NULL);
    }


  struct confinement *cf = NULL;

  unsigned long len
    //= scm_num2ulong (scm_length (scm_confinement),
    //		     0, "CF_guile_get");
    = scm_to_uint64 (scm_length (scm_confinement));
  if (len == 0)
    {
      // no confinement
      return (cf);
    }
  else if (len < 4)
    {
      fprintf (stderr, "CF_guile_get: %s is too short\n", var);
      return (NULL);
    }

  double epsilon
    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (0)),
    //		   "CF_guile_get");
    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (0)));
  double r0
    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (1)),
    //		   "CF_guile_get");
    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (1)));

  // get the string
  char *str_cf = NULL;
  //SCM scm_conf = scm_list_ref (scm_confinement, scm_int2num (2));
  SCM scm_conf = scm_list_ref (scm_confinement, scm_from_int32 (2));
  if (scm_is_string (scm_conf))
    {
      str_cf = scm_to_locale_string (scm_conf);
    }

  if (strcmp (str_cf, "sphere") == 0)
    {
      if (len != 4)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for sphere, number of parameter must be 1\n");
	}
      else
	{
	  double R
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  cf = CF_init (0, // sphere
			R,
			0.0, // r
			0.0, 0.0, 0.0, // x, y, z
			0.0, // R2
			0.0, // L
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else if (strcmp (str_cf, "sphere+hole") == 0)
    {
      if (len != 5)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for sphere+hole, number of parameter must be 2\n");
	}
      else
	{
	  double R
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  double r
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (4)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (4)));
	  cf = CF_init (1, // sphere+hole
			R,
			r,
			0.0, 0.0, 0.0, // x, y, z
			0.0, // R2
			0.0, // L
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else if (strcmp (str_cf, "cylinder") == 0)
    {
      if (len != 7)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for cylinder, number of parameter must be 4\n");
	}
      else
	{
	  double r
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  double x
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (4)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (4)));
	  double y
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (5)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (5)));
	  double z
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (6)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (6)));
	  cf = CF_init (2, // cylinder
			0.0, // R,
			r,
			x, y, z,
			0.0, // R2
			0.0, // L
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else if (strcmp (str_cf, "dumbbell") == 0)
    {
      if (len != 7)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for dumbbell, number of parameter must be 4\n");
	}
      else
	{
	  double R
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  double R2
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (4)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (4)));
	  double L
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (5)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (5)));
	  double r
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (6)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (6)));
	  cf = CF_init (3, // dumbbell
			R,
			r,
			0.0, 0.0, 0.0, // x, y, z
			R2,
			L,
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else if (strcmp (str_cf, "hex2d") == 0)
    {
      if (len != 6)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for hex2d, number of parameter must be 3\n");
	}
      else
	{
	  double R
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  double r
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (4)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (4)));
	  double L
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (5)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (5)));
	  cf = CF_init (4, // hex2d
			R,
			r,
			0.0, 0.0, 0.0, // x, y, z
			0.0, // R2
			L,
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else if (strcmp (str_cf, "porous") == 0)
    {
      if (len != 5)
	{
	  fprintf (stderr, "CF_guile_get:"
		   " for hex2d, number of parameter must be 2\n");
	}
      else
	{
	  double R
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (3)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (3)));
	  double L
	    //= scm_num2dbl (scm_list_ref (scm_confinement, scm_int2num (4)),
	    //		   "CF_guile_get");
	    = scm_to_double (scm_list_ref (scm_confinement, scm_from_int32 (4)));
	  cf = CF_init (5, // porous
			R,
			0.0,
			0.0, 0.0, 0.0, // x, y, z
			0.0, // R2
			L,
			0, // flag_LJ
			epsilon,
			r0);
	  CHECK_MALLOC (cf, "CF_guile_get");
	}
    }
  else
    {
      fprintf (stderr, "CF_guile_get: invalid confinement %s\n",
	       str_cf);
    }
  free (str_cf);

  return (cf); // success
}
