#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.5 2006/10/26 01:52:35 kichiki Exp $

AC_PREREQ(2.59)
AC_INIT(libstokes, 0.4.1, kichiki@users.sourceforge.net)
AC_CONFIG_SRCDIR([src/bench.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_LIBTOOL
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h])
AC_CHECK_HEADERS([cblas.h])
AC_CHECK_HEADERS([blas.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
