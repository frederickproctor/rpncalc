AC_DEFUN([ACX_RPNCALC],
	[AC_HAVE_LIBRARY(m)]
	[AC_CHECK_FUNCS([pow sqrt])]
	[AC_SEARCH_LIBS([clock_gettime], [rt])]
	[AC_MSG_CHECKING([for RPN calculator functions])]
	[AC_ARG_WITH(rpncalc,
		[ --with-rpncalc=<path to rpncalc>  Specify path to rpncalc directory],	dirs=$withval,dirs="/usr/local")]
	for dir in $dirs ; do
		if test -f $dir/include/rpncalc.h ; then RPNCALC_DIR=$dir ; break ; fi
	done
	if test x$RPNCALC_DIR = x ; then
	[AC_MSG_ERROR([not found, specify using --with-rpncalc=<path to rpncalc>])]
	else
	[AC_MSG_RESULT([$RPNCALC_DIR])]
	RPNCALC_CFLAGS="-I$RPNCALC_DIR/include"
	RPNCALC_LIBS="-L$RPNCALC_DIR/lib -lrpncalc"
dnl put HAVE_RPNCALC in config.h
	[AC_DEFINE(HAVE_RPNCALC, 1, [Define non-zero if you have rpncalc.])]
dnl put RPNCALC_DIR in Makefile
	[AC_SUBST(RPNCALC_DIR)]
	[AC_SUBST(RPNCALC_CFLAGS)]
	[AC_SUBST(RPNCALC_LIBS)]
	fi
dnl enable HAVE_RPNCALC test in Makefile
	[AM_CONDITIONAL(HAVE_RPNCALC, test x$RPNCALC_DIR != x)]
)
