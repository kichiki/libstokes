/* guile interface for libstokes
 * Copyright (C) 2006-2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifdef GUILE16
#include <guile/gh.h>
#endif // GUILE16

/* this utility function (original name is does_scm_symbol_exist)
 * written by Michael Gran is taken from
 * http://www.lonelycactus.com/guilebook/c319.html
 */
int
guile_check_symbol (const char *name)
{ 
  SCM sym; 
  SCM var;

  //sym = scm_str2symbol (name);
  sym = scm_string_to_symbol (scm_from_locale_string (name));

  /* Check to see if the symbol exists */
  //var = scm_sym2var (sym, 
  //		     scm_current_module_lookup_closure (), 
  //		     SCM_BOOL_F);
  var = scm_symbol_interned_p (sym);

  if (var != SCM_BOOL_F) 
    {
      return 1; // true
    }
  
  return 0; // false
}

/* check boolean
 * INPUT
 *  var : string of the SCM variable
 * OUTPUT
 *  true  : if var is not nil
 *  false : if var is nil or even not defined
 */
int
guile_get_bool (const char * var)
{
  SCM scm_symbol;
  SCM scm_param;
  int i;

  i = 0; // FALSE
  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (SCM_NFALSEP (scm_param)) // TRUE (!FALSE)
	{
	  i = 1; // TRUE
	}
    }

  return (i);
}

/*
 * INPUT
 *  var : string of the SCM variable
 *  i0  : default value (for the undefined case)
 * OUTPUT
 */
int
guile_get_int (const char * var, int i0)
{
  SCM scm_symbol;
  SCM scm_param;
  int i;

  i = i0;
  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (scm_number_p (scm_param))
	{
	  //i = scm_num2int (scm_param, 0, "guile_get_int");
	  i = scm_to_int32 (scm_param);
	}
    }

  return (i);
}

/*
 * INPUT
 *  var : string of the SCM variable
 *  d0  : default value (for the undefined case)
 * OUTPUT
 */
double
guile_get_double (const char * var, double d0)
{
  SCM scm_symbol;
  SCM scm_param;
  double d;

  d = d0;
  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (scm_number_p (scm_param))
	{
	  //d = scm_num2dbl (scm_param, "guile_get_double");
	  d = scm_to_double (scm_param);
	}
    }

  return (d);
}

/* get doubles from SCM list or vector with length check
 * OUTPUT
 *  returned value : 0 = failed (not defined)
 *                   1 = success
 */
int
guile_get_doubles (const char * var, int n, double * x)
{
  SCM scm_symbol;
  SCM scm_param;
  unsigned long len;
  int i;


  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (SCM_NFALSEP (scm_list_p (scm_param)))
	{
	  //len = scm_num2ulong (scm_length (scm_param),
	  //		       0, "guile_get_doubles");
	  len = scm_to_uint64 (scm_length (scm_param));
	  if (len != n)
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	  for (i = 0; i < n; i ++)
	    {
	      //x[i] = scm_num2dbl (scm_list_ref (scm_param,
	      //					scm_int2num (i)),
	      //			  "guile_get_doubles");
	      x[i] = scm_to_double (scm_list_ref (scm_param,
						  scm_from_int32 (i)));
	    }
	  return (1); // success
	}
#ifdef GUILE16
      else if (SCM_VECTORP(scm_param))
	{
	  len = SCM_VECTOR_LENGTH (scm_param);
#else // !GUILE16
      else if (scm_is_vector(scm_param))
	{
	  len = scm_c_vector_length (scm_param);
#endif // GUILE16
	  if (len != n)
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	  for (i = 0; i < n; i ++)
	    {
	      //x[i] = scm_num2dbl (scm_vector_ref (scm_param,
	      //				  scm_int2num (i)),
	      //		  "guile_get_doubles");
	      x[i] = scm_to_double (scm_vector_ref (scm_param,
						    scm_from_int32 (i)));
	    }
	  return (1); // success
	}
    }

  return (0); // failed
}

/* get doubles from SCM list or vector with unknown length
 * OUTPUT
 *  returned value : NULL = failed (not defined)
 */
double *
guile_get_doubles_ (const char * var)
{
  SCM scm_symbol;
  SCM scm_param;
  unsigned long len;
  int i;
  double *x = NULL;


  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (SCM_NFALSEP (scm_list_p (scm_param)))
	{
	  //len = scm_num2ulong (scm_length (scm_param),
	  //		       0, "guile_get_doubles_");
	  len = scm_to_uint64 (scm_length (scm_param));
	  x = (double *) malloc (sizeof (double) * len);
	  for (i = 0; i < len; i ++)
	    {
	      //x[i] = scm_num2dbl (scm_list_ref (scm_param,
	      //				scm_int2num (i)),
	      //		  "guile_get_doubles_");
	      x[i] = scm_to_double (scm_list_ref (scm_param,
						  scm_from_int32 (i)));
	    }
	}
#ifdef GUILE16
      else if (SCM_VECTORP(scm_param))
	{
	  len = SCM_VECTOR_LENGTH (scm_param);
#else // !GUILE16
      else if (scm_is_vector(scm_param))
	{
	  len = scm_c_vector_length (scm_param);
#endif // GUILE16
	  x = (double *) malloc (sizeof (double) * len);
	  for (i = 0; i < len; i ++)
	    {
	      //x[i] = scm_num2dbl (scm_vector_ref (scm_param,
	      //				  scm_int2num (i)),
	      //		  "guile_get_doubles_");
	      x[i] = scm_to_double (scm_vector_ref (scm_param,
						    scm_from_int32 (i)));
	    }
	}
    }

  return (x);
}

/* get length of SCM list or vector
 * OUTPUT
 *  returned value : length (not defined, 0 is returned)
 */
int
guile_get_length (const char * var)
{
  SCM scm_symbol;
  SCM scm_param;
  unsigned long len;


  len = 0;
  if (guile_check_symbol (var) != 0) // true == found
    {
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (SCM_NFALSEP (scm_list_p (scm_param)))
	{
	  //len = scm_num2ulong (scm_length (scm_param),
	  //		       0, "guile_get_length");
	  len = scm_to_uint64 (scm_length (scm_param));
	}
#ifdef GUILE16
      else if (SCM_VECTORP(scm_param))
	{
	  len = SCM_VECTOR_LENGTH (scm_param);
#else // !GUILE16
      else if (scm_is_vector(scm_param))
	{
	  len = scm_c_vector_length (scm_param);
#endif // GUILE16
	}
    }

  return (len);
}

/*
 */
char *
guile_get_string (const char * var)
{
  SCM scm_param;
  char * str = NULL;
#ifdef GUILE16
  size_t str_len;

  scm_param = gh_eval_str (var);
  if (gh_string_p (scm_param))
    {
      str = gh_scm2newstr (scm_param, &str_len);
    }
#else // !GUILE16
  SCM scm_symbol;

  scm_symbol = scm_c_lookup (var);
  scm_param = scm_variable_ref (scm_symbol);
  if (scm_is_string (scm_param))
    {
      str = scm_to_locale_string (scm_param);
    }
#endif // GUILE16

  return str;
}

/*
 */
FILE *
guile_open_file (const char * var, const char * mode)
{
  char * filename = NULL;
  FILE * f = NULL;


  filename = guile_get_string (var);
  if (filename != NULL)
    {
      f = fopen (filename, mode);
      if (f == NULL)
	{
	  fprintf (stderr, "Cannot open file %s!\n", filename);
	  exit (1);
	}
    }
  free (filename);

  return f;
}


/* load scm script file into the libstokes C system
 */
void
guile_load (const char *file)
{
  scm_init_guile ();
  scm_c_primitive_load (file);
}
