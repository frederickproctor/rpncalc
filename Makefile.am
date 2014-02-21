bin_PROGRAMS = rpn variate

rpn_SOURCES = src/rpnmain.c
rpn_LDADD = -L. -lrpncalc
rpn_DEPENDENCIES = librpncalc.a

variate_SOURCES = src/variates.c src/variates.h
variate_CFLAGS = -DMAIN

lib_LIBRARIES = librpncalc.a
librpncalc_a_SOURCES = src/rpncalc.c src/rpncalc.h src/variates.c src/variates.h src/ptime.c src/ptime.h

include_HEADERS = src/rpncalc.h src/variates.h src/ptime.h
