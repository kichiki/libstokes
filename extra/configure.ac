#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Copyright (C) 1993-2006,2017 Kengo Ichiki <kengoichiki@gmail.com>

AC_PREREQ(2.59)
AC_INIT(libstokes-extra, 0.1.0, kengoichiki@gmail.com)
AC_CONFIG_SRCDIR([src/KIrand.c])
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
