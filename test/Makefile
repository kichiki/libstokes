# $Id: Makefile,v 1.28 2008/10/31 06:07:45 kichiki Exp $

# for FreeBSD
CC = cc
RYUON_PREFIX = /usr/local
F2C_LIBS = -L/usr/local/lib/gcc-4.2.5 -lgfortran
# plain g2c library
#F2C_LIBS = -lg2c

#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
# for FreeBSD
OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math\
	-march=prescott -mfpmath=sse

# for SUSE-Linux
#RYUON_PREFIX = /home/ichiki/RYUON
#F2C_LIBS = -lgfortran
#OPTS = \
#	-mtune=opteron -m64 \
#	-pipe -O3 -fomit-frame-pointer -ffast-math

# for i1200 (openSUSE)
#OPTS = \
#	-march=pentium3m \
#	-pipe -O3 -fomit-frame-pointer -ffast-math

# for glacier (with intel compiler)
#CC = /global/software/intel-9/cc/bin/icc
#RYUON_PREFIX = /global/home/ichiki/RYUON
#F2C_LIBS = -L/global/software/intel-9/fc/lib -lifcore
#OPTS = -O3

# for robson gcc version 3 (32bit)
#CC = /usr/bin/gcc
#RYUON_PREFIX = /home/ichiki/RYUON3-32

# for robson gcc version 3 (64bit)
#CC = /usr/bin/gcc -m64
#RYUON_PREFIX = /home/ichiki/RYUON3

# for robson gcc version 4
#CC = /usr/local/bin/gcc4
#RYUON_PREFIX = /home/ichiki/RYUON4-32

# for robson gcc version 4 (64bit)
#CC = /usr/local/bin/gcc4 -m64
#RYUON_PREFIX = /home/ichiki/RYUON4

# for robson gcc3/g77
#F2C_LIBS = -lg2c
# for robson gcc4/gfortran
#F2C_LIBS = -L/usr/local/gcc-4.2.2/lib/gcc/powerpc64-suse-linux/4.2.2 -lgfortran

# for robson
#OPTS = -O3 -mcpu=970 -mtune=970 -maltivec


RYUON_BIN = $(RYUON_PREFIX)/bin
RYUON_INC = -I$(RYUON_PREFIX)/include
RYUON_LIB = -L$(RYUON_PREFIX)/lib

LD = libtool --mode=link --tag=CC $(CC)
#LD = $(CC)

# for header file besides libstokes.h
RYUON_SRC = ../src

LIBSTOKES = ../src/libstokes.la
#LIBSTOKES = $(RYUON_LIB) -lstokes


#MATH_LIBS	= \
#	$(F2C_LIBS)\
#	-lblas -llapack\
#	-lcblas -lf77blas -latlas\
#	-larpack
# for FreeBSD
MATH_LIBS = \
	$(F2C_LIBS)\
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack
# for SUSE-Linux
# for openSUSE
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-lblas -llapack\
#	-larpack
#MATH_LIBS = \
#	-llapack -lptcblas -lptf77blas -latlas\
#	-larpack
# for glacier
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-L/usr/lib -llapack -lblas \
#	$(RYUON_LIB) -larpack
# for robson (gcc/g77, 64bit)
#MATH_LIBS = $(F2C_LIBS) -L/usr/local/lib64 -lcblas -latlas -L/usr/local/lib64 -larpack
# for robson (gcc/g77, 64bit or 32bit)
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	$(RYUON_LIB) -larpack -llapack -lblas
# for robson (gcc4/gfortran, 64bit)
#MATH_LIBS = -L/usr/local/lib64 -larpack -lessl -lxlf90_r
# for robson (gcc4/gfortran, 64bit or 32bit)
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	$(RYUON_LIB) -larpack -llapack -lblas

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`
# for glacier
#GSL_INC = `gsl-config --cflags`
#GSL_LIB = `gsl-config --libs`
# for robson
#GSL_INC = `$(RYUON_BIN)/gsl-config --cflags`
#GSL_LIB = `$(RYUON_BIN)/gsl-config --libs`

NITSOL_INC = -I/usr/local/include
NITSOL_LIB = -L/usr/local/lib -lnitsol
# for SUSE-Linux
#NITSOL_INC = -I/home/ichiki/RYUON/include
#NITSOL_LIB = -L/home/ichiki/RYUON/lib -lnitsol
# for glacier
# for robson
#NITSOL_INC = $(RYUON_INC)
#NITSOL_LIB = $(RYUON_LIB) -lnitsol

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`
# for robson
#GUILE_INC = `$(RYUON_BIN)/guile-config compile`
#GUILE_LIB = `$(RYUON_BIN)/guile-config link`

NETCDF_INC = 
NETCDF_LIB = -lnetcdf
# for glacier
# for robson
#NETCDF_INC = $(RYUON_INC)
#NETCDF_LIB = $(RYUON_LIB) -lnetcdf


CFLAGS = \
	-Wall		\
	$(OPTS)		\
	-I$(RYUON_SRC)	\
	$(RYUON_INC)	\
	$(NETCDF_INC)	\
	$(GUILE_INC)	\
	$(GSL_INC)	\
	$(NITSOL_INC)

LDFLAGS =

LIBS = \
	$(LIBSTOKES)		\
	-lm			\
	$(RYUON_LIB) -liter	\
	$(NETCDF_LIB)		\
	$(GUILE_LIB)		\
	$(GSL_LIB)		\
	$(NITSOL_LIB)		\
	$(MATH_LIBS)

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
	check-solve-cubic.o	\
	check-grid.o		\
	check-ev-dh-grid.o
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
