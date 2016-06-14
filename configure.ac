#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(CCA, 1.0, davidgil@usc.edu)
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE()

# Checks for programs.
AC_PROG_CC
AC_PROG_FC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_ARG_WITH([etree-include-path],
  [AS_HELP_STRING([--with-etree-include-path],
    [location of the Etree headers])],
  [ETREE_INCL="-I$withval"],
  [ETREE_INCL=""])
AC_ARG_WITH([etree-lib-path],
  [AS_HELP_STRING([--with-etree-lib-path], [location of the Etree libraries])],
  [ETREE_LIB="-L$withval -letree"],
  [ETREE_LIB="-letree"])

# Get required Proj.4 library and include locations
AC_ARG_WITH([proj4-include-path],
  [AS_HELP_STRING([--with-proj4-include-path],
    [location of the Proj.4 headers])],
  [PROJ4_INCL="-I$withval"],
  [PROJ4_INCL=""])
AC_ARG_WITH([proj4-lib-path],
  [AS_HELP_STRING([--with-proj4-lib-path], [location of the Proj.4 libraries])],
  [PROJ4_LIB="-L$withval -lproj -pthread"],
  [PROJ4_LIB="-lproj -pthread"])

# Check required ETREE installation
CFLAGS="$ETREE_INCL $CHECK_CFLAGS"
LDFLAGS="$CHECK_LDFLAGS $ETREE_LIB"
AC_CHECK_LIB(etree, etree_open, [AC_CHECK_HEADER([etree.h], [], [AC_MSG_ERROR(["Etree headers not found in $ETREE_INCL; use --with-etree-include-path"])], [AC_INCLUDES_DEFAULT])],[AC_MSG_ERROR(["Etree library not found; use --with-etree-lib-path"])])

# Check required PROJ4 installation
CFLAGS="$PROJ4_INCL $CHECK_CFLAGS"
LDFLAGS="$CHECK_LDFLAGS $PROJ4_LIB"
AC_CHECK_LIB(proj, pj_init_plus, [AC_CHECK_HEADER([proj_api.h], [], [AC_MSG_ERROR([Proj4 header not found in $PROJ4_INCL; use --with-proj4-include-path"])
  ], [AC_INCLUDES_DEFAULT])],[AC_MSG_ERROR(["Proj4 library not found; use --with-proj4-lib-path"])], [-pthread -lm])

# Set final CFLAGS and LDFLAGS
CFLAGS="$CHECK_CFLAGS $ETREE_INCL $PROJ4_INCL"
LDFLAGS="$CHECK_LDFLAGS $ETREE_LIB $PROJ4_LIB"
LDFLAGS="$LDFLAGS -lm"

AC_CONFIG_FILES([Makefile
				 data/Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
