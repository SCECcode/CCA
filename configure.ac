dnl -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(CCA, 1.1.0, software@scec.org)
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE(foreign tar-pax)

dnl
dnl Dependencies
dnl

dnl Etree (Euclid)
AC_ARG_ENABLE([etree],
    [AC_HELP_STRING([--enable-etree],
        [install Etree library @<:@default=yes@:>@])],
	[if test "$enableval" = yes ; then install_etree=yes; else install_etree=no; fi],
	[install_etree=yes])
AM_CONDITIONAL([INSTALL_ETREE], [test "$install_etree" = yes])
AC_ARG_WITH([etree-incdir],
  [AS_HELP_STRING([--with-etree-incdir], [directory containing Etree header files])],
  [with_etree_incdir="$withval"],
  [with_etree_incdir=no])
AC_SUBST(with_etree_incdir)
AC_ARG_WITH([etree-libdir],
  [AS_HELP_STRING([--with-etree-libdir], [directory containing Etree library])],
  [with_etree_libdir="$withval"],
  [with_etree_libdir=no])
AC_SUBST(with_etree_libdir)

dnl SQLITE
AC_ARG_ENABLE([sqlite],
    [AC_HELP_STRING([--enable-sqlite],
        [install Sqlite @<:@default=no@:>@])],
	[if test "$enableval" = "yes" ; then install_sqlite=yes; else install_sqlite=no; fi],
	[install_sqlite=no])
AM_CONDITIONAL([INSTALL_SQLITE], [test "$install_sqlite" = yes])
AC_ARG_WITH([sqlite-incdir],
    [AC_HELP_STRING([--with-sqlite-incdir],
        [location of sqlite header files @<:@default=no@:>@])],
	[with_sqlite_incdir=$withval],
	[with_sqlite_incdir=no])
AC_SUBST(with_sqlite_incdir)
AC_ARG_WITH([sqlite-libdir],
    [AC_HELP_STRING([--with-sqlite-libdir],
        [location of sqlite library @<:@default=no@:>@])],
	[with_sqlite_libdir=$withval],
	[with_sqlite_libdir=no])
AC_SUBST(with_sqlite_libdir)

dnl TIFF
AC_ARG_ENABLE([tiff],
    [AC_HELP_STRING([--enable-tiff],
        [install Tiff @<:@default=no@:>@])],
	[if test "$enableval" = "yes" ; then install_tiff=yes; else install_tiff=no; fi],
	[install_tiff=no])
AM_CONDITIONAL([INSTALL_TIFF], [test "$install_tiff" = yes])
AC_ARG_WITH([tiff-incdir],
    [AC_HELP_STRING([--with-tiff-incdir],
        [location of tiff header files @<:@default=no@:>@])],
	[with_tiff_incdir=$withval],
	[with_tiff_incdir=no])
AC_SUBST(with_tiff_incdir)
AC_ARG_WITH([tiff-libdir],
    [AC_HELP_STRING([--with-tiff-libdir],
        [location of tiff library @<:@default=no@:>@])],
	[with_tiff_libdir=$withval],
	[with_tiff_libdir=no])
AC_SUBST(with_tiff_libdir)

dnl PROJ
AC_ARG_ENABLE([proj],
    [AC_HELP_STRING([--enable-proj],
        [install Proj @<:@default=yes@:>@])],
	[if test "$enableval" = "yes" ; then install_proj=yes; else install_proj=no; fi],
	[install_proj=yes])
AM_CONDITIONAL([INSTALL_PROJ], [test "$install_proj" = yes])
AC_ARG_WITH([proj-incdir],
    [AC_HELP_STRING([--with-proj-incdir], [directory containing Proj header files @<:@default=no@:>@])],
	[with_proj_incdir="$withval"],
	[with_proj_incdir=no])
AC_SUBST(with_proj_incdir)
AC_ARG_WITH([proj-libdir],
    [AC_HELP_STRING([--with-proj-libdir], [directory containing Proj library @<:@default=no@:>@])],
	[with_proj_libdir="$withval"],
	[with_proj_libdir=no])
AC_SUBST(with_proj_libdir)


##check optional large data path 
if test x"$CVM_LARGEDATA_DIR" != x; then
   # test directory existence
   CCA_LARGEDATA_DIR=$CVM_LARGEDATA_DIR/model/cca_i06
   AC_CHECK_FILE($CVM_LARGEDATA_DIR/model/cca_i06/density.dat, [], [AC_MSG_ERROR(["CCA data not found in cvm_largedata_dir"])])
   AM_CONDITIONAL([WITH_CCA_LARGEDATA_DIR], true)
   AC_SUBST(CCA_LARGEDATA_DIR)
else
   AM_CONDITIONAL(WITH_CCA_LARGEDATA_DIR, false)
fi

dnl
dnl Command line tools
dnl

AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_RANLIB
AC_PROG_LN_S

dnl CURL (used to download model files)
AC_PATH_PROG(CURL, curl)
if test -z "$CURL" ; then
  AC_MSG_FAILURE([cannot find 'curl' program.])
  CURL=`echo "Error: curl is not installed." ; false`
fi

dnl TAR
AC_PATH_PROG(TAR, tar)
if test -z "$TAR" ; then
  AC_MSG_FAILURE([cannot find 'tar' program.])
  TAR=`echo "Error: tar is not installed." ; false`
fi

dnl
dnl Setup environment so dependencies are used in build
dnl

CFLAGS="$CFLAGS"
CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib -L$prefix/lib64 $LDFLAGS"

dnl
dnl Verify configuration
dnl

dnl ETREE (required)
if test "$install_etree" = yes ; then
  ETREE_INCLUDES=""
  ETREE_LDFLAGS="-letree"
else
  SCEC_ETREE_HEADER
  SCEC_ETREE_LIB
fi
if test "$with_etree_incdir" != no ; then
  ETREE_INCLUDES="-I$with_etree_incdir"
fi
if test "$with_etree_libdir" != no ; then
  ETREE_LDFLAGS="-L$with_etree_libdir -letree"
fi

AC_SUBST(ETREE_INCLUDES)
AC_SUBST(ETREE_LDFLAGS)

dnl proj (required)
if test "$install_proj" = yes ; then
  dnl sqlite (required by proj)
  if test "$install_sqlite" = yes ; then
    SQLITE3_INCLUDES="-I$prefix/include"
    SQLITE3_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_sqlite_incdir" != no; then
      SQLITE3_INCLUDES="-I$with_sqlite_incdir"
    fi
    if test "$with_sqlite_libdir" != no; then
      SQLITE3_LDFLAGS="-L$with_sqlite_libdir"
    fi
    SCEC_SQLITE3_HEADER
    SCEC_SQLITE3_LIB
  fi
  AC_SUBST(SQLITE3_INCLUDES)
  AC_SUBST(SQLITE3_LDFLAGS)

  dnl tiff (required by proj)
  if test "$install_tiff" = yes ; then
    TIFF_INCLUDES="-I$prefix/include"
    TIFF_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_tiff_incdir" != no; then
      TIFF_INCLUDES="-I$with_tiff_incdir"
    fi
    if test "$with_tiff_libdir" != no; then
      TIFF_LDFLAGS="-L$with_tiff_libdir"
    fi
    SCEC_TIFF_HEADER
    SCEC_TIFF_LIB
  fi
  AC_SUBST(TIFF_INCLUDES)
  AC_SUBST(TIFF_LDFLAGS)

  PROJ_INCLUDES="-I$prefix/include"
  PROJ_LDFLAGS="-L$prefix/lib"
else
  SCEC_PROJ_HEADER
  SCEC_PROJ_LIB
fi
if test "$with_proj_incdir" != no ; then
  PROJ_INCLUDES="-I$with_proj_incdir"
fi
if test "$with_proj_libdir" != no ; then
  PROJ_LDFLAGS="-L$with_proj_libdir -lproj"
fi
AC_SUBST(PROJ_INCLUDES)
AC_SUBST(PROJ_LDFLAGS)




AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  dependencies/Makefile
  data/Makefile
  ])

AC_OUTPUT

AM_CONDITIONAL([BUILD_SHARED_LIBRARY], [test "$enable_shared" = yes])
if test "$enable_shared" = yes; then
  AC_MSG_NOTICE([Building shared libraries. UCVM must be built with shared libraries for plugins to work.])
else
  AC_MSG_WARN([Building static libraries. Shared libraries are not enabled.])
fi


dnl End of file
