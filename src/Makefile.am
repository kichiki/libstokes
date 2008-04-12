# for libstokes
# $Id: Makefile.am,v 1.21 2008/04/12 18:19:05 kichiki Exp $

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
        XA.c			\
        YA.c			\
        YB.c			\
        XC.c			\
        YC.c			\
        XG.c			\
        YG.c			\
        YH.c			\
        XM.c			\
        YM.c			\
        ZM.c			\
	XP.c			\
	XQ.c			\
	twobody.c		\
        XA-slip.c		\
        YA-slip.c		\
        YB-slip.c		\
        XC-slip.c		\
        YC-slip.c		\
        XG-slip.c		\
        YG-slip.c		\
        YH-slip.c		\
        XM-slip.c		\
        YM-slip.c		\
        ZM-slip.c		\
	twobody-slip.c		\
	lub.c			\
	minv-poly.c		\
	ewald.c			\
	non-ewald.c		\
	ewald-3f-matrix.c	\
	ewald-3ft-matrix.c	\
	ewald-3fts-matrix.c	\
	lub-matrix.c		\
	matrix.c		\
	dgetri_c.c		\
	dgemv_c.c		\
	dgemm_c.c		\
	bench.c			\
	stokes-nc.c		\
	stokes-nc-read.c	\
	stokes-guile.c		\
	coll.c			\
	periodicity.c		\
	bonds.c			\
	bonds-guile.c		\
	angles.c		\
	angles-guile.c		\
	ode.c			\
	ode-quaternion.c	\
	chebyshev.c		\
	dsaupd_c.c		\
	dnaupd_c.c		\
	dgeev_c.c		\
	dpotrf_c.c		\
	KIrand.c		\
	excluded-volume.c	\
	brownian.c		\
	bd-imp.c


libstokes_la_LDFLAGS = -version-info 6:0:0

include_HEADERS = libstokes.h
