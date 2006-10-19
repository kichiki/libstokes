/* guile interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-guile.c,v 5.1 2006/10/19 03:11:30 ichiki Exp $
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
#include <guile/gh.h>


/* this utility function (original name is does_scm_symbol_exist)
 * written by Michael Gran is taken from
 * http://www.lonelycactus.com/guilebook/c319.html
 */
int
guile_check_symbol (const char *name)
{ 
  SCM sym; 
  SCM var;

  sym = scm_str2symbol (name);

  /* Check to see if the symbol exists */
  var = scm_sym2var (sym, 
		     scm_current_module_lookup_closure (), 
		     SCM_BOOL_F);
  
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
      /*
      scm_param = gh_eval_str (var);
      if (gh_number_p (scm_param))
	{
	  i = gh_scm2int (scm_param);
	}
      */
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (scm_number_p (scm_param))
	{
	  i = scm_num2int (scm_param, 0, "guile_get_int");
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
      /*
      scm_param = gh_eval_str (var);
      if (gh_number_p (scm_param))
	{
	  d = gh_scm2double (scm_param);
	}
      */
      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      if (scm_number_p (scm_param))
	{
	  d = scm_num2dbl (scm_param, "guile_get_double");
	}
    }

  return (d);
}

/*
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
      /*
      scm_param = gh_eval_str (var);
      if (gh_vector_p (scm_param))
	{
	  if (gh_vector_length (scm_param) == n)
	    {
	      gh_scm2doubles (scm_param, x);
	      return (1); // success
	    }
	  else
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	}
      else if (gh_list_p (scm_param))
	{
	  if (gh_length (scm_param) == n)
	    {
	      gh_scm2doubles (scm_param, x);
	      return (1); // success
	    }
	  else
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	}
      */

      scm_symbol = scm_c_lookup (var);
      scm_param = scm_variable_ref (scm_symbol);
      //if (scm_list_p (scm_param))
      if (SCM_NFALSEP (scm_list_p (scm_param)))
	{
	  len = scm_num2ulong (scm_length (scm_param),
			       0, "guile_get_doubles");
	  if (len != n)
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	  for (i = 0; i < n; i ++)
	    {
	      x[i] = scm_num2dbl (scm_list_ref (scm_param,
						scm_int2num (i)),
				  "guile_get_doubles");
	    }
	  return (1); // success
	}
      else if (SCM_VECTORP(scm_param)/*scm_vector_p (scm_param)*/)
	{
	  len = SCM_VECTOR_LENGTH (scm_param);
	  if (len != n)
	    {
	      fprintf (stderr, "wrong size of %s in guile_get_doubles().\n",
		       var);
	      return (0); // failed
	    }
	  for (i = 0; i < n; i ++)
	    {
	      x[i] = scm_num2dbl (scm_vector_ref (scm_param,
						  scm_int2num (i)),
				  "guile_get_doubles");
	    }
	  return (1); // success
	}
    }

  return (0); // failed
}

/*
 */
char *
guile_get_string (const char * var)
{
  SCM scm_param;
  char * str = NULL;
  int str_len;


  scm_param = gh_eval_str (var);
  if (gh_string_p (scm_param))
    {
      str = gh_scm2newstr (scm_param, &str_len);
    }

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