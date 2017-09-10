# Makefile for libstokes-core/test
# Copyright (C) 1993-2008,2017 Kengo Ichiki <kengoichiki@gmail.com>

# for Mac
CC = gcc
RYUON_PREFIX = /Users/ichiki/WORK/201708/RYUON
OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
F2C_LIBS = -L/opt/local/lib/gcc49 -lgfortran

# for FreeBSD
#CC = cc
#RYUON_PREFIX = /usr/local
#F2C_LIBS = -L/usr/local/lib/gcc-4.2.5 -lgfortran
# plain g2c library
#F2C_LIBS = -lg2c

#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math
# for FreeBSD
#OPTS = -pipe -O3 -fomit-frame-pointer -ffast-math\
#	-march=prescott -mfpmath=sse

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

#LD = libtool --mode=link --tag=CC $(CC)
LD = $(CC)

# for header file besides libstokes.h
RYUON_SRC = ../src

#LIBSTOKES = ../src/libstokes.la
LIBSTOKES_CORE = $(RYUON_LIB) -lstokes-core


# for Mac
MATH_LIBS = \
	$(F2C_LIBS)\
	-L/opt/local/lib -llapack -lblas\
	-larpack

#MATH_LIBS	= \
#	$(F2C_LIBS)\
#	-lblas -llapack\
#	-lcblas -lf77blas -latlas\
#	-larpack
# for FreeBSD
#MATH_LIBS = \
#	$(F2C_LIBS)\
#	-L/usr/local/lib -lalapack -lf77blas -lcblas -latlas\
#	-larpack
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


CFLAGS = \
	-Wall		\
	$(OPTS)		\
	-I$(RYUON_SRC)	\
	$(RYUON_INC)

LDFLAGS =

LIBS = \
	$(LIBSTOKES_CORE)	\
	-lm			\
	$(RYUON_LIB) -liter	\
	$(MATH_LIBS)

all:	test-all	\
	test-twobody	\
	make-lub-table	\
	test-twobody-slip	\
	test-new


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
	check-minv-poly.o	\
	check-lub-poly.o	\
	check-ewald-poly.o	\
	check-ewald-shear.o	\
	check-twobody-slip.o	\
	check-matrix.o

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


# test-new
test-new-OBJ =\
	test-new.o	\
	check.o	\
	check-non-ewald-new.o	\
	check-ewald-new.o	\
	check-ewald-3fts-new_res0.o	\
	check-ewald-3fts-new_res.o	\
	check-ewald-3fts-new_mob.o	\
	check-ewald-3fts-new_mix.o	\
	check-ewald-3fts-new_res_lub.o	\
	check-ewald-3fts-new_mob_lub.o	\
	check-ewald-3fts-new_mix_lub.o

test-new:	$(test-new-OBJ)
	$(LD) $(LDFLAGS) -o test-new	\
	$(test-new-OBJ)		\
	$(LIBS)


clean:	
	rm -f *~ *.o	\
	test-all	\
	test-twobody	\
	make-lub-table	\
	test-twobody-slip	\
	test-new
