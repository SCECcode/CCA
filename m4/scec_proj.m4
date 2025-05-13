dnl -*- Autoconf -*-


AC_DEFUN([SCEC_PROJ_HEADER], [
  scec_save_cppflags=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $PROJ_INCLUDES"
  AC_LANG(C)
  AC_CHECK_HEADER([proj.h], [], [
    AC_MSG_ERROR([Proj header not found; try --with-proj-incdir="<Proj include dir>"])
  ])
  CPPFLAGS=$scec_save_cppflags
])


AC_DEFUN([SCEC_PROJ_LIB], [
  scec_save_cppflags=$CPPFLAGS
  scec_save_ldflags=$LDFLAGS
  CPPFLAGS="$CPPFLAGS $PROJ_INCLUDES"
  LDFLAGS="$LDFLAGS $PROJ_LDFLAGS"
  AC_CHECK_LIB(proj, proj_create, [],
      [AC_MSG_ERROR("Proj library not found; use --with-proj-libdir")])
  CPPFLAGS=$scec_save_cppflags
  LDFLAGS=$scec_save_ldflags
])


dnl end of file
