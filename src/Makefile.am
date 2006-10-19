# for libstokes
# $Id: Makefile.am,v 1.4 2006/10/19 04:18:59 ichiki Exp $

lib_LTLIBRARIES = libstokes.la

libstokes_la_SOURCES = \
	stokes.c		\
	ewald-3f.c		\
	ewald-3ft.c		\
	ewald-3fts.c		\
	ewald-2f.c		\
	ewald-2ft.c		\
	ewald-2fts.c		\
	f.c			\
	ft.c			\
	fts.c			\
	two-body-res.c		\
	lub.c			\
	ewald.c			\
	ewald-3f-matrix.c	\
	ewald-3ft-matrix.c	\
	ewald-3fts-matrix.c	\
	lub-matrix.c		\
	matrix.c		\
	dgetri_c.c		\
	bench.c			\
	stokes-nc.c		\
	stokes-nc-read.c	\
	stokes-guile.c

libstokes_la_LDFLAGS = -version-info 3:0:0

include_HEADERS = libstokes.h
