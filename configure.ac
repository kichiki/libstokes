#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.2 2006/10/12 17:26:54 ichiki Exp $

AC_PREREQ(2.59)
AC_INIT(libstokes, 0.3, kichiki@users.sourceforge.net)
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

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
