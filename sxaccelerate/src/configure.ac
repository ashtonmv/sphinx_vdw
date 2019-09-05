# ----------------------------------------------------------------------------
#
#           The general purpose cross platform C/C++ framework
#
#                       S x A c c e l e r a t e
#
#           Home:       https://www.sxlib.de
#           License:    Apache 2
#           Authors:    see src/AUTHORS
#
# ----------------------------------------------------------------------------

AC_INIT([S/PHI/nX Accelerate Libraries], [2.10.0], [https://www.sxlib.de], [sxaccelerate])
SHARED_VERSION_INFO="1:0:0"

AC_CANONICAL_HOST

# --- we are (ab)using the autoconf install mechanism to build sxdist
ac_default_prefix=`pwd`

   
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([nonstd/Makefile])
AC_CONFIG_FILES([util/Makefile])
AC_CONFIG_FILES([util/examples/Makefile])
AC_CONFIG_FILES([util/tests/Makefile])
AC_CONFIG_FILES([graph/Makefile])
#AC_CONFIG_FILES([graph/examples/Makefile])
AC_CONFIG_FILES([math/Makefile])
AC_CONFIG_FILES([math/examples/Makefile])
AC_CONFIG_FILES([math/tests/Makefile])
AC_CONFIG_FILES([math/benchmarks/Makefile])
AC_CONFIG_FILES([ipc/Makefile])
AC_CONFIG_FILES([ipc/tests/Makefile])
AC_CONFIG_FILES([ipc/examples/Makefile])
AC_CONFIG_FILES([ipc/examples/plugins/Makefile])
AC_CONFIG_FILES([fs/Makefile])
AC_CONFIG_FILES([fs/examples/Makefile])
AC_CONFIG_FILES([fs/tests/Makefile])
AC_CONFIG_FILES([parserkit/Makefile])
AC_CONFIG_FILES([parserkit/examples/Makefile])
AC_CONFIG_FILES([parserkit/examples/01_hello/Makefile])
AC_CONFIG_FILES([parserkit/examples/02_stateful_scanner/Makefile])
AC_CONFIG_FILES([parserkit/examples/03_unicode/Makefile])
AC_CONFIG_FILES([parserkit/examples/04_include_files/Makefile])
AC_CONFIG_FILES([parserkit/examples/05_parser_stack/Makefile])
AC_CONFIG_FILES([parserkit/examples/06_schema/Makefile])
AC_CONFIG_FILES([json/Makefile])
AC_CONFIG_FILES([json/examples/Makefile])
AC_CONFIG_FILES([io/Makefile])
AC_CONFIG_FILES([io/examples/Makefile])
AC_CONFIG_FILES([io/examples/plugins/Makefile])
AC_CONFIG_FILES([mpi/Makefile])
AC_CONFIG_FILES([mpi/examples/Makefile])
AC_CONFIG_FILES([mpi/tests/Makefile])
AC_CONFIG_FILES([apps/Makefile])
AC_CONFIG_FILES([apps/libexec/Makefile])
AC_CONFIG_FILES([apps/tests/Makefile])
AC_CONFIG_FILES([apps/tests/bin/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([share/manual/Makefile])
AC_CONFIG_FILES([share/manual/manual.sx])
AC_CONFIG_FILES([share/manual/preamble.a4paper.tex])
AC_CONFIG_FILES([share/manual/preamble.a5paper.tex])
AC_CONFIG_FILES([share/manual/title.a4paper.tex])
AC_CONFIG_FILES([share/manual/title.a5paper.tex])
AC_CONFIG_FILES([system/Makefile])
AC_CONFIG_FILES([system/sxwarnstat.pl])
AC_CONFIG_FILES([system/Doxyfile])
#AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([gpu/Makefile])


# --- SxConfig.h ------------------------------------------------------------
AC_CONFIG_HEADERS([util/SxUnixConfig.h])
AH_TOP([
#ifndef _SX_UNIX_CONFIG_H_
#define _SX_UNIX_CONFIG_H_
])
AH_BOTTOM([
#define SX_REGISTER_CLASS(name) void name##_register() { name var; } \
                                extern char SX_NO_OPERATION
#endif /* _SX_UNIX_CONFIG_H_ */
])
# ---------------------------------------------------------------------------


# --- Package options
SX_TITLE([S/PHI/nX Accelerate Libraries])
SX_ARG_WITH(  [examples],  [$srcdir/util/examples], [WITH_EXAMPLES], [no], 
              [compile with source code examples])
SX_ARG_WITH(  [apps],      [$srcdir/apps], [WITH_APPS], [yes], 
              [compile with tools])
SX_ARG_WITH(  [sxmath], [.], [USE_SXMATH], [no], [whether to compile sxmath])
SX_ARG_WITH(  [tests],  [$srcdir/fs/tests], [WITH_TESTS], [no], 
              [compile S/PHI/nX with source code tests])
AM_CONDITIONAL([COND_BUILD_EXAMPLES], [test x"$ac_cv_with_examples" = x"yes"])
AM_CONDITIONAL([COND_BUILD_SXMATH],   [test x"$ac_cv_with_sxmath"   = x"yes"])
AM_CONDITIONAL([COND_BUILD_TESTS],    [test x"$ac_cv_with_tests"    = x"yes"])
AM_CONDITIONAL([COND_BUILD_APPS],     [test x"$ac_cv_with_apps"     = x"yes"])

# --- if math is needed, get settings for optimzed number-crunching
if test x"$ac_cv_with_sxmath"   = x"yes" ; then
   SX_PARAM_CPUSPECIFIC=yes
else
   SX_PARAM_CPUSPECIFIC=no
fi

# --- check environment and operating system
SX_ENVIRONMENT

# --- older versions of autoconf don't provide abs_top_srcdir/abs_top_builddir
SXBASE="\$(abs_top_srcdir)"
if test x"$ac_abs_top_srcdir" = x; then
   abs_top_srcdir=`cd $top_srcdir && pwd`
   AC_SUBST(abs_top_srcdir)
fi
if test x"$ac_abs_top_builddir" = x; then
   abs_top_builddir=`cd $top_builddir && pwd`
   AC_SUBST(abs_top_builddir)
fi
AC_SUBST([SXBASE])

# --- Checks for programs and tools
AC_PROG_LN_S
AC_PROG_AWK
AM_PROG_LEX
AC_PROG_YACC
AC_CHECK_TOOL(DOXYGEN, [doxygen])

# --- sxpublisher
SXQ_CHECK_DOC

# --- setup compiler and determine suitable compiler flags
SX_COMPFLAGS

## --- setup optional packages
#SX_OPT_DIRS=""
#if test x"${ac_cv_with_sxlm}" = x"yes"; then
#   SX_OPT_DIRS="$SX_OPT_DIRS sxkeygen"
#   AC_CONFIG_FILES([sxkeygen/Makefile])
#fi
#AC_SUBST([SX_OPT_DIRS], [$SX_OPT_DIRS])


# --- part from system/sxsetup.sh.in
AC_SUBST(SXGET_RC,    [`pwd`/add-ons/sxget.rc])
AC_PATH_TOOL(PERL,    [perl])
AC_PATH_TOOL(GREP,    [grep])
AC_PATH_TOOL(FIND,    [find])
AC_PATH_TOOL(RM,      [rm])
AC_PATH_PROGS(ZCAT,   [gzcat zcat])
AC_PATH_PROGS(TAIL,   [gtail tail])
AC_PATH_PROGS(AWK,    [gawk awk])
AC_PATH_PROGS(SED,    [gsed sed])

# --- static or dynamic linking?
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
AM_CONDITIONAL([BUILD_SHAREDLIBS], [test x"$enable_shared" = x"yes"])

## --- md5 or md5sum?
#AC_CHECK_TOOL(SX_MD5,  md5sum)
#if test -z ${ac_cv_prog_ac_ct_SX_MD5}; then
#   AC_CHECK_TOOL(SX_MD5, md5)
#fi

# --- Check if required packages were found
SX_REQUIRED(CC,           GNU C compiler)   
SX_REQUIRED(CXX,          GNU C++ compiler)
SX_REQUIRED(AWK,          awk)
SX_REQUIRED(YACC,         bison)
SX_REQUIRED(LEX,          flex)


# --- Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([arpa/inet.h fcntl.h netdb.h netinet/in.h stdlib.h signal.h string.h sys/socket.h sys/types.h sys/time.h sys/wait.h sys/resource.h unistd.h winsock.h winsock2.h])

# --- activate numlibs support
SX_CHECK_NUMLIBS
AM_CONDITIONAL([COND_BUILD_SXMATH], [test x"$ac_cv_with_sxmath" = x"yes"])

# --- Check for available functions
AC_CHECK_FUNCS(mkstemp posix_memalign)
AC_FUNC_STRERROR_R

# --- avoid lazy linking of dynamic objects
test x"$enable_shared" != x"yes" && LIBS="$LIBS -all-static"

# --- Substitute variables in makefile [system/make.conf.in]
AC_SUBST(SHARED_VERSION_INFO)
AC_SUBST(VERSION)
AC_SUBST(MODE)
AC_SUBST(DOXYGEN)
AC_SUBST(DOXYWIZARD)
AC_SUBST(SUBDIRS)

AC_CONFIG_COMMANDS([chmodwarn], [chmod +x system/sxwarnstat.pl])
AC_CONFIG_COMMANDS([helper],    [test -f GNUmakefile || rm -f GNUmakefile])
AC_CONFIG_FILES([system/sxmake])
AC_CONFIG_COMMANDS([sxmake], [chmod +x system/sxmake])

sxaccel_src=$(cd $srcdir;pwd)
sxaccel_build=$(pwd)
AC_SUBST([SXACCELERATE_BUILD],[$sxaccel_build])
AC_SUBST([SXACCELERATE_SRC],[$sxaccel_src])
AC_SUBST([SXACCELERATE],[$prefix])

mkdir -p system/$SXBINTARGET
cat > system/$SXBINTARGET/sxconfig.dat <<CONFIG_END
ac_cv_with_sxaccelerate_src=${sxaccel_src}
ac_cv_with_sxaccelerate_build=${sxaccel_build}
ac_cv_with_mklpath=${ac_cv_with_mklpath}
ac_cv_with_numlibs=${ac_cv_with_numlibs:-no}
ac_cv_with_sxmath=${ac_cv_with_sxmath:-no}
ac_cv_enable_acml=${ac_cv_enable_acml:-no}
ac_cv_enable_acmlfft=${ac_cv_enable_acmlfft:-no}
ac_cv_enable_atlas=${ac_cv_enable_atlas:-no}
ac_cv_enable_cpuspecific=${ac_cv_enable_cpuspecific:-no}
ac_cv_enable_debug=${ac_cv_enable_debug}
ac_cv_enable_fftw=${ac_cv_enable_fftw:-no}
ac_cv_enable_goto=${ac_cv_enable_goto:-no}
ac_cv_enable_mkl=${ac_cv_enable_mkl:-no}
ac_cv_enable_mklfft=${ac_cv_enable_mklfft:-no}
ac_cv_enable_mpi=${ac_cv_enable_mpi:-no}
ac_cv_enable_hdf5=${ac_cv_enable_hdf5:-no}
ac_cv_enable_openmp=${ac_cv_enable_openmp:-no}
ac_cv_enable_profile=${ac_cv_enable_profile:-no}
CONFIG_END

AC_OUTPUT
