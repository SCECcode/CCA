dnl -*- Autoconf -*-


AC_DEFUN([SCEC_ETREE_HEADER], [
  scec_save_cppflags=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $ETREE_INCLUDES"
  AC_LANG(C)
  AC_CHECK_HEADER([euclid/etree.h], [], [
    AC_MSG_ERROR([Etree (Euclid) header not found; try --with-etree-incdir="<Etree include dir>"])
  ])
  CPPFLAGS=$scec_save_cppflags
])


AC_DEFUN([SCEC_ETREE_LIB], [
  scec_save_cppflags=$CPPFLAGS
  scec_save_ldflags=$LDFLAGS
  CPPFLAGS="$CPPFLAGS $ETREE_INCLUDES"
  LDFLAGS="$LDFLAGS $ETREE_LDFLAGS"
  AC_CHECK_LIB(etree, etree_open, [],
      [AC_MSG_ERROR("Etree library not found; use --with-etree-libdir")])
  CPPFLAGS=$scec_save_cppflags
  LDFLAGS=$scec_save_ldflags
])


dnl end of file
