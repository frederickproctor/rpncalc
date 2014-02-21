AC_INIT(rpn, 1.0)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE(rpn, 1.0)
AC_CONFIG_SRCDIR([src/rpncalc.c])
AC_CONFIG_HEADERS([config/config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([readline/readline.h readline/history.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_HAVE_LIBRARY(m)
AC_CHECK_FUNCS([pow sqrt])
AC_SEARCH_LIBS([clock_gettime], [rt])
AC_HAVE_LIBRARY(history)
AC_HAVE_LIBRARY(curses)
AC_HAVE_LIBRARY(readline, , , -lcurses)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
