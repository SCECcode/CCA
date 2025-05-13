# -*- Autoconf -*-

# ----------------------------------------------------------------------
# SCEC_TIFF_HEADER
# ----------------------------------------------------------------------
AC_DEFUN([SCEC_TIFF_HEADER], [
  scec_save_cppflags=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $TIFF_INCLUDES"
  AC_LANG(C)
  AC_REQUIRE_CPP
  AC_CHECK_HEADER([tiff.h], [], [
    AC_MSG_ERROR([tiff header not found; try --with-tiff-incdir=<tiff include dir>"])
  ])dnl
  CPPFLAGS=$scec_save_cppflags
])dnl SCEC_TIFF_HEADER


# ----------------------------------------------------------------------
# SCEC_TIFF_LIB
# ----------------------------------------------------------------------
AC_DEFUN([SCEC_TIFF_LIB], [
  scec_save_CPPFLAGS=$CPPFLAGS
  scec_save_LDFLAGS=$LDFLAGS
  scec_save_libs=$LIBS
  CPPFLAGS="$CPPFLAGS $TIFF_INCLUDES"
  LDFLAGS="$LDFLAGS $TIFF_LDFLAGS"
  AC_LANG(C)
  AC_REQUIRE_CPP
  AC_CHECK_LIB(tiff, TIFFOpen, [],[
    AC_MSG_ERROR([tiff library not found; try --with-tiff-libdir=<tiff lib dir>])
  ])dnl
  CPPFLAGS=$scec_save_cppflags
  LDFLAGS=$scec_save_ldflags
  LIBS=$scec_save_libs
])dnl SCEC_TIFF_LIB


dnl end of file
