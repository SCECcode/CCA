# -*- Autoconf -*-

# ----------------------------------------------------------------------
# SCEC_SQLITE_HEADER
# ----------------------------------------------------------------------
AC_DEFUN([SCEC_SQLITE3_HEADER], [
  scec_save_cppflags=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $SQLITE3_INCLUDES"
  AC_LANG(C)
  AC_REQUIRE_CPP
  AC_CHECK_HEADER([sqlite3.h], [], [
    AC_MSG_ERROR([sqlite3 header not found; try --with-sqlite-incdir=<sqlite3 include dir>"])
  ])dnl
  CPPFLAGS=$scec_save_cppflags
])dnl SCEC_SQLITE3_HEADER


# ----------------------------------------------------------------------
# SCEC_SQLITE_LIB
# ----------------------------------------------------------------------
AC_DEFUN([SCEC_SQLITE3_LIB], [
  scec_save_CPPFLAGS=$CPPFLAGS
  scec_save_LDFLAGS=$LDFLAGS
  scec_save_libs=$LIBS
  CPPFLAGS="$CPPFLAGS $SQLITE3_INCLUDES"
  LDFLAGS="$LDFLAGS $SQLITE3_LDFLAGS"
  AC_LANG(C)
  AC_REQUIRE_CPP
  AC_CHECK_LIB(sqlite3, sqlite3_open, [],[
    AC_MSG_ERROR([sqlite3 library not found; try --with-sqlite-libdir=<sqlite3 lib dir>])
  ])dnl
  CPPFLAGS=$scec_save_cppflags
  LDFLAGS=$scec_save_ldflags
  LIBS=$scec_save_libs
])dnl SCEC_SQLITE3_LIB


dnl end of file
