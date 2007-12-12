# $Id: Makefile,v 1.15 2007/12/12 06:32:31 kichiki Exp $

CC           = cc
#LD           = libtool --mode=link $(CC)
LD           = libtool --mode=link --tag=CC $(CC)

SRC_DIR = ../src

GSL_INC = `pkg-config --cflags gsl`
GSL_LIB = `pkg-config --libs gsl`

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

CFLAGS	     = \
	-Wall		\
	-I$(SRC_DIR)	\
	-pipe -O3 -march=prescott -fomit-frame-pointer\
	-mfpmath=sse -ffast-math\
	$(GUILE_INC)	\
	$(GSL_INC)

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.3 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
	-larpack

LDFLAGS      =
LIBS         = \
	../src/libstokes.la	\
	-lm			\
	-L/usr/local/lib -liter	\
	$(ATLAS_LIBS)		\
	-lnetcdf		\
	$(GUILE_LIB)		\
	$(GSL_LIB)

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
	check-ewald.o		\
	check-ewald-poly.o	\
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
	check-sqrt-dgeev.o	\
	check-matrix.o		\
	check-list-ex.o


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
