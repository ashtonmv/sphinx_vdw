# ---------------------------------------------------------------------------
#
#      The ab-initio based multiscale library
#
#                  S / P H I / n X
#
#      Copyright:  Max-Planck-Institute for Iron Research
#                  40237 Duesseldorf, Germany
#
#      Contact:    https://sxlib.mpie.de
#      Authors:    see sphinx/AUTHORS
#      License:    see sphinx/LICENSE
#
# ---------------------------------------------------------------------------

AC_INIT([SPHInX],[2.6.1], [https://sxrepo.mpie.de])
SHARED_VERSION_INFO="1:7:0"

AC_CANONICAL_HOST

# --- we are (ab)using the autoconf install mechanism to build sxdist
ac_default_prefix=`pwd`

   
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([geom/Makefile])
AC_CONFIG_FILES([geom/examples/Makefile])
AC_CONFIG_FILES([dirac/Makefile])
AC_CONFIG_FILES([dirac/examples/Makefile])
AC_CONFIG_FILES([dft/Makefile])
AC_CONFIG_FILES([dft/examples/Makefile])
AC_CONFIG_FILES([struct/Makefile])
AC_CONFIG_FILES([struct/examples/Makefile])
AC_CONFIG_FILES([exx/Makefile])
AC_CONFIG_FILES([classic/Makefile])
AC_CONFIG_FILES([add-ons/Makefile])
AC_CONFIG_FILES([tools/Makefile])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([share/sphinx/Makefile])
AC_CONFIG_FILES([share/sphinx/std/Makefile])
AC_CONFIG_FILES([share/sphinx/structures/Makefile])
AC_CONFIG_FILES([add-ons/sxget:add-ons/sxget.pl.in], [chmod +x add-ons/sxget])
AM_CONDITIONAL(COND_PLAYGROUND,[ls $srcdir/playground | grep -q cpp])
AM_COND_IF([COND_PLAYGROUND],AC_CONFIG_FILES([playground/Makefile]))

# --- Package options
SX_TITLE([S/PHI/nX Physics])
SX_ARG_WITH(  [addons], [$srcdir/add-ons],  [WITH_ADDONS], [yes], 
              [compile S/PHI/nX add-ons executables. libsxext.la will always be compiled.])
SX_ARG_WITH(  [tools], [$srcdir/tools],  [WITH_TOOLS], [yes], 
              [compile S/PHI/nX tools])
SX_ARG_WITH(  [vimsyntax], [$srcdir/system],  [WITH_VIM_SYNTAX], [no], 
              [generate the vim syntax files for S/PHI/nX source code])
SX_ARG_WITH(  [examples],  [$srcdir/dft/examples], [WITH_EXAMPLES], [no], 
              [compile S/PHI/nX with source code examples])
SX_ARG_WITH(  [tests],  [$srcdir/tests], [WITH_TESTS], [no], 
              [compile S/PHI/nX with source code tests])
SX_ARG_ENABLE([sxgemmm],      [USE_SXGEMMM], [no], 
              [compile package with mmm algorithm for projectors])

# dummy-config will be written by autoheader and protects SxDFTConfig.h
AC_CONFIG_HEADERS([dummy-config dft/SxDFTConfig.h])
# --- check environment and operating system
SX_ENVIRONMENT
SX_CHECK_SXACCELERATE
SXDIR=`cd ..; pwd`
AC_SUBST([SXDIR])

AM_CONDITIONAL([COND_BUILD_ADDONS],   [test x"$ac_cv_with_addons"   = x"yes"])
AM_CONDITIONAL([COND_BUILD_TOOLS],    [test x"$ac_cv_with_tools"    = x"yes"])
AM_CONDITIONAL([COND_BUILD_EXAMPLES], [test x"$ac_cv_with_examples" = x"yes"])
AM_CONDITIONAL([COND_BUILD_DOC], [test x"$ac_cv_with_doc" = x"yes"])

# --- Checks for programs and tools
AC_PROG_LN_S
AC_PROG_AWK
AM_PROG_LEX
AC_PROG_YACC
AC_CHECK_TOOL(DOXYGEN, [doxygen])

# --- setup compiler and determine suitable compiler flags
SX_COMPFLAGS

# --- static or dynamic linking?
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
AM_CONDITIONAL([BUILD_SHAREDLIBS], [test x"$enable_shared" = x"yes"])

# --- numerical libraries
SX_CHECK_NUMLIBS

# --- part from system/sxsetup.sh.in
AC_PATH_TOOL(PERL,    [perl])
AC_PATH_TOOL(GREP,    [grep])
AC_PATH_TOOL(FIND,    [find])
AC_PATH_TOOL(RM,      [rm])
AC_PATH_PROGS(ZCAT,   [gzcat zcat])
AC_PATH_PROGS(TAIL,   [gtail tail])
AC_PATH_PROGS(AWK,    [gawk awk])
AC_PATH_PROGS(SED,    [gsed sed])

# --- Check if required packages were found
SX_REQUIRED(CC,           GNU C compiler)   
SX_REQUIRED(CXX,          GNU C++ compiler)
SX_REQUIRED(AWK,          awk)
SX_REQUIRED(YACC,         bison)
SX_REQUIRED(LEX,          flex)

# --- Check for available functions
AC_CHECK_LIB(m, sincos, [AC_DEFINE([HAVE_SINCOS],[1],  
             [Define to 1 if your system has 'sincos'.])])

# --- avoid lazy linking of dynamic objects
test x"$enable_shared" != x"yes" && LIBS="$LIBS -all-static"

# --- Substitute variables in makefile [system/make.conf.in]
AC_SUBST(SHARED_VERSION_INFO)
AC_SUBST(DOXYGEN)
AC_SUBST(SUBDIRS)

AC_CONFIG_COMMANDS([helper],    [test -f GNUmakefile || rm -f GNUmakefile])
if test "$WITH_VIM_SYNTAX" = "yes" ; then
  AC_CONFIG_COMMANDS([vimsyntax], [system/sxvimsyntax.sh .])
fi
if test "$ABS_TOP_BUILDDIR" != `cd $srcdir;pwd` ; then
  AC_CONFIG_COMMANDS([sharelink],[(cd $srcdir/ ; pwd) > share/.srcdir])
fi

AC_OUTPUT
