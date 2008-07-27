# $Id: Makefile,v 1.26 2008/07/27 01:11:28 kichiki Exp $

CC           = cc
#LD           = libtool --mode=link $(CC)
LD           = libtool --mode=link --tag=CC $(CC)

OPTS = \
	-march=prescott -mfpmath=sse \
	-pipe -O3 -fomit-frame-pointer -ffast-math
# for SUSE-Linux
#OPTS = \
#	-mtune=opteron -m64 \
#	-pipe -O3 -fomit-frame-pointer -ffast-math
# for i1200 (openSUSE)
#OPTS = \
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math


# for header file besides libstokes.h
RYUON_SRC = ../src

RYUON_INC = -I/usr/local/include
RYUON_LIB = -L/usr/local/lib
# for SUSE-Linux
#RYUON_INC = -I/home/ichiki/RYUON/include
#RYUON_LIB = -L/home/ichiki/RYUON/lib

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.5 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
# for openSUSE
#ATLAS_LIBS	= \
#	-lgfortran -lblas -llapack\
#	-larpack
#ATLAS_LIBS	= \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack

CFLAGS	     = \
	-Wall		\
	$(OPTS)		\
	-I$(RYUON_SRC)	\
	$(RYUON_INC)	\
	$(GUILE_INC)	\
	$(GSL_INC)	\
	$(NITSOL_INC)

LDFLAGS      =

LIBS         = \
	../src/libstokes.la	\
	-lm			\
	$(RYUON_LIB) -liter	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)		\
	$(NITSOL_LIB)

all:	test-all	\
	test-twobody	\
	make-lub-table	\
	test-twobody-slip


# test-all
test-all-OBJ =\
	test-all.o \
	check.o	\
	check-mul-matrices.o \
	check-lapack-inv.o \
	check-lapack-solve-lin.o \
	check-solve-gen-linear.o \
	check-twobody.o	\
	check-poly.o	\
	check-minv.o	\
	check-minv-poly.o	\
	check-lub-poly.o	\
	check-mob-fts.o		\
	check-ewald.o		\
	check-ewald-poly.o	\
	check-ewald-shear.o	\
	check-ode-quaternion.o	\
	check-twobody-slip.o	\
	check-dsaupd_c.o	\
	check-dnaupd_c.o	\
	check-chebyshev.o	\
	check-bd.o		\
	check-dgeev_c.o		\
	check-dpotrf_c.o	\
	check-brownian.o	\
	check-bd-imp.o		\
	check-bd-imp-nitsol.o	\
	check-bd-imp-fast.o	\
	check-sqrt-dgeev.o	\
	check-matrix.o		\
	check-bonds-guile.o	\
	check-angles.o		\
	check-excluded-volume-guile.o	\
	check-ev-dh.o		\
	check-ev-dh-guile.o	\
	check-ev-LJ.o		\
	check-ev-LJ-guile.o	\
	check-confinement.o	\
	check-confinement-guile.o	\
	check-KIrand.o		\
	check-bead-rod.o	\
	check-bead-rod-guile.o	\
	check-solve-cubic.o
# below the check-matrix.o
#	check-list-ex.o		\

test-all:	$(test-all-OBJ)
	$(LD) $(LDFLAGS) -o test-all	\
	$(test-all-OBJ)		\
	$(LIBS)

# test-twobody
test-twobody-OBJ =\
	test-twobody.o

test-twobody:	$(test-twobody-OBJ)
	$(LD) $(LDFLAGS) -o test-twobody	\
	$(test-twobody-OBJ)		\
	$(LIBS)

# make-lub-table
make-lub-table-OBJ =\
	make-lub-table.o

make-lub-table:	$(make-lub-table-OBJ)
	$(LD) $(LDFLAGS) -o make-lub-table	\
	$(make-lub-table-OBJ)	\
	$(LIBS)


# test-twobody-slip
test-twobody-slip-OBJ =\
	test-twobody-slip.o

test-twobody-slip:	$(test-twobody-slip-OBJ)
	$(LD) $(LDFLAGS) -o test-twobody-slip	\
	$(test-twobody-slip-OBJ)		\
	$(LIBS)

clean:	
	rm -f *~ *.o	\
	test-all	\
	test-twobody	\
	make-lub-table	\
	test-twobody-slip
