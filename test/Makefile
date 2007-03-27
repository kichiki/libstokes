# $Id: Makefile,v 1.1 2007/03/19 04:10:37 kichiki Exp $

CC           = cc
LD           = libtool --mode=link $(CC)

SRC_DIR = ../src

GUILE_INC = `guile-config compile`
GUILE_LIB = `guile-config link`

#CFLAGS	     = -Wall -I/usr/local/include -O3 -mtune=pentium
CFLAGS	     = -Wall\
		-I$(SRC_DIR)\
		-pipe -O3 -march=prescott -fomit-frame-pointer\
		-mfpmath=sse -ffast-math\
		$(GUILE_INC)

#ATLAS_LIBS	= \
#	-lg2c -lblas -llapack\
#	-lcblas -lf77blas -latlas
ATLAS_LIBS	= \
	-L/usr/local/lib/gcc-4.2.0 -lgfortran \
	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas

LDFLAGS      =
LIBS         = ../src/libstokes.la\
		-lm\
		-L/usr/local/lib -liter\
		$(ATLAS_LIBS)\
		-lnetcdf\
		$(GUILE_LIB)

all:	test-all


# test-all
test-all-OBJ =\
	test-all.o \
	check-mul-matrices.o \
	check-lapack-inv.o \
	check-solve-gen-linear.o

test-all:	$(test-all-OBJ)
	$(LD) $(LDFLAGS) -o test-all	\
	$(test-all-OBJ)		\
	$(LIBS)


clean:	
	rm -f *~ *.o	\
	test-all