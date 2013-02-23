#!/bin/sh
rm -rf aclocal.m4 autom4te.cache\
 config.* configure configure.scan\
 depcomp install-sh libtool ltmain.sh\
 Makefile Makefile.in missing mkinstalldirs stamp-h1\
 src/Makefile.in

# for FreeBSD
#aclocal -I /usr/local/share/aclocal
aclocal
#libtoolize
glibtoolize
autoheader
autoconf
automake -a
