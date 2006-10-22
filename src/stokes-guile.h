/* header file for stokes-guile.c --
 * guile interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-guile.h,v 5.2 2006/10/22 04:17:51 kichiki Exp $
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
#ifndef	_STOKES_GUILE_H_
#define	_STOKES_GUILE_H_


/* this utility function (original name is does_scm_symbol_exist)
 * written by Michael Gran is taken from
 * http://www.lonelycactus.com/guilebook/c319.html
 */
int
guile_check_symbol (const char *name);

/* check boolean
 * INPUT
 *  var : string of the SCM variable
 * OUTPUT
 *  true  : if var is not nil
 *  false : if var is nil or even not defined
 */
int
guile_get_bool (const char * var);

/*
 * INPUT
 *  var : string of the SCM variable
 *  i0  : default value (for the undefined case)
 * OUTPUT
 */
int
guile_get_int (const char * var, int i0);

/*
 * INPUT
 *  var : string of the SCM variable
 *  d0  : default value (for the undefined case)
 * OUTPUT
 */
double
guile_get_double (const char * var, double d0);

/* get doubles from SCM list or vector with length check
 * OUTPUT
 *  returned value : 0 = failed (not defined)
 *                   1 = success
 */
int
guile_get_doubles (const char * var, int n, double * x);
/* get doubles from SCM list or vector with unknown length
 * OUTPUT
 *  returned value : NULL = failed (not defined)
 */
double *
guile_get_doubles_ (const char * var);
/* get length of SCM list or vector
 * OUTPUT
 *  returned value : length (not defined, 0 is returned)
 */
int
guile_get_length (const char * var);

/*
 */
char *
guile_get_string (const char * var);

/*
 */
FILE *
guile_open_file (const char * var, const char * mode);


#endif /* !_STOKES_GUILE_H_ */
