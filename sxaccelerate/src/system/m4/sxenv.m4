AC_DEFUN([SX_ENVIRONMENT],
[


SX_ARG_WITH(  [bintarget],   [$srcdir], [WITH_DBGDIR], [], 
              [Append directory to lib, libexec, bin, sbin.])

SXBINTARGET=$ac_cv_with_bintarget
if test x"$SXBINTARGET" = x""; then :; else
   libdir="$libdir/$SXBINTARGET"
   libexecdir="$libexecdir/$SXBINTARGET"
   bindir="$bindir/$SXBINTARGET"
   sbindir="$sbindir/$SXBINTARGET"
fi

AC_ARG_VAR([DIST_OS],
           AS_HELP_STRING([Overwrite OS / distribution when cross-compiling]))


UNAME=`uname`
WHOAMI=`whoami`
TODAY=`date +%D`
PROCESSOR=`uname -m`

ABS_TOP_SRCDIR=`cd $srcdir; pwd`
ABS_TOP_BUILDDIR=`pwd`
SXLIBTOOL="$ABS_TOP_BUILDDIR/libtool"

isMobile="no"
isAndroid="no"
iIOS="no"

# --- determine distribution
if test -n "$DIST_OS"; then
   if test -z "$TARGET_ARCH"; then
      AC_MSG_ERROR([Environment variable TARGET_ARCH not set.])
   fi

   case "$DIST_OS" in
      "Android"*)

         DIST_OS=Android
         DIST_NAME=Android
         DIST_VERSION=0
         DIST_MAJOR=0
         DIST_MINOR=0
         DIST_PATCH=0
         DIST_REPO_VER=0
         PROCESSOR=$TARGET_ARCH
         isMobile="yes"
         isAndroid="yes"
         AC_DEFINE(SX_ANDROID, "1", [Setup Android Compilation])
         AC_DEFINE(SX_MOBILE,  "1", [Setup mobile compilation])
         ;;
      "iOS"*)
         DIST_OS=iOS
         DIST_NAME=iOS
         DIST_VERSION=0
         DIST_MAJOR=0
         DIST_MINOR=0
         DIST_PATCH=0
         DIST_REPO_VER=0
         PROCESSOR=$TARGET_ARCH
         isMobile="yes"
         isIOS="yes"
         AC_DEFINE(SX_IOS,    "1", [Setup iOS Compilation])
         AC_DEFINE(SX_MOBILE, "1", [Setup mobile compilation])
         ;;
      *)
         AC_MSG_ERROR([Environment variable DIST_OS=$DIST_OS is invalid (Android, iOS).])
   esac
elif test -f /etc/redhat-release; then
   PKG_TYPE=rpm
   case `cat /etc/redhat-release` in
   "CentOS release"*)
      DIST_OS=Linux
      DIST_NAME=CentOS
      DIST_VERSION=`cut -f 3 -d ' ' /etc/redhat-release`
      DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
      DIST_MINOR=`echo $DIST_VERSION | cut -f 2 -d '.'`
      DIST_PATCH=0
#     DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      DIST_REPO_VER="${DIST_MAJOR}"
      ;;
   "CentOS Linux release"*)
      DIST_OS=Linux
      DIST_NAME=CentOS
      DIST_VERSION=`cut -f 4 -d ' ' /etc/redhat-release`
      DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
      DIST_MINOR=`echo $DIST_VERSION | cut -f 2 -d '.'`
      DIST_PATCH=0
#     DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      DIST_REPO_VER="${DIST_MAJOR}"
      ;;
   "Scientific Linux release"*)
      DIST_OS=Linux
      DIST_NAME=SL
      DIST_VERSION=`cut -f 4 -d ' ' /etc/redhat-release`
      DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
      DIST_MINOR=`echo $DIST_VERSION | cut -f 2 -d '.'`
      DIST_PATCH=0
#     DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      ;;
   "Red Hat Enterprise Linux Server release 6.5"*)
      DIST_OS=Linux
      DIST_NAME=RHEL
      DIST_VERSION=`cut -f 7 -d ' ' /etc/redhat-release`
      DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
      DIST_MINOR=`echo $DIST_VERSION | cut -f 2 -d '.'`
      DIST_PATCH=0
#     DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
      ;;
   "Fedora release"*)
      DIST_OS=Linux
      DIST_NAME=Fedora
      DIST_VERSION=`cut -f 3 -d ' ' /etc/redhat-release`
      DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
      DIST_MINOR=0
      DIST_PATCH=0
      DIST_REPO_VER="${DIST_MAJOR}"
      ;;
   esac
elif test -f /etc/debian_version; then
   DIST_OS=Linux
   DIST_NAME=Debian
   DIST_VERSION=`cat /etc/debian_version`
   DIST_MAJOR=`cut -f 1 -d '.' /etc/debian_version`
   DIST_MINOR=`cut -f 2 -d '.' /etc/debian_version`
   DIST_PATCH=`cut -f 3 -d '.' /etc/debian_version`
   DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
   PKG_TYPE=deb
elif test -f /etc/SuSE-release; then
   DIST_OS=Linux
   DIST_NAME=SuSE
   DIST_VERSION=`cat /etc/SuSE-release|grep VERSION|awk '{print $3}'`
   PKG_TYPE=rpm
elif test -f /etc/os-release; then
   DIST_OS=Linux
   DIST_NAME=`sed -ne 's/^ID="\?\([^"]*\)"\?/\1/p' /etc/os-release`
   DIST_VERSION=`sed -ne 's/^VERSION_ID="\?\([^"]*\)"\?/\1/p' /etc/os-release`
   # Can we figure out the package type? maybe via DIST_NAME or ID_LIKE...
   #PKG_TYPE=
elif test -d /Library/Preferences; then
   PKG_TYPE=pkg
   DIST_OS=Darwin
   DIST_NAME="MacOSX"
   DIST_VERSION=`/usr/bin/sw_vers|grep ProductVersion|cut -f 2`
   DIST_MAJOR=`echo $DIST_VERSION | cut -f 1 -d '.'`
   DIST_MINOR=`echo $DIST_VERSION | cut -f 2 -d '.'`
   DIST_PATCH=`echo $DIST_VERSION | cut -f 3 -d '.'`
   DIST_REPO_VER="${DIST_MAJOR}.${DIST_MINOR}"
   if test x"$DIST_PATCH" = x""; then
      DIST_PATCH="0"
   fi
   sysctl=`/usr/sbin/sysctl hw | grep hw.cpu64bit_capable | cut -f 2 -d ':' | sed s/\ //g`
   if test x"$sysctl" = x"1"; then
      PROCESSOR="x86_64"
   fi
else
   echo "WARNING: sxenv.m4 could not determine distribution"
fi
DIST_VERSION_L=`echo ${DIST_MAJOR}${DIST_MINOR}${DIST_PATCH}L`

REPOPROC=$PROCESSOR
if test x"$REPOPROC" = x"x86_64"; then
   if test x"$PKG_TYPE" = x"deb"; then
      REPOPROC="amd64"
   fi
fi

AC_CHECK_TOOLS([GDB], [ddd gdb])
AC_CHECK_PROG([MEMTRACER], [valgrind], [`which valgrind` --leak-check=full --show-reachable=yes --run-libc-freeres=no --log-fd=1])
GDB=`which $GDB`
case "${host}" in
   *-mingw32)
      isWindowsOS="yes"
      LIBEXT="dll"
      DLLLIBEXT="-0.dll"
      IMPLIBEXT=".dll.a"
      CP_A=`which cp`" -a"
      AC_CHECK_PROG([LDD], [ldd], [ldd])
      LDD_AWK_COL="3"
      ;;
   *-darwin*)
      isWindowsOS="no"
      buildMacBundles="yes"
      LIBEXT="dylib"
      DLLLIBEXT=".$LIBEXT"
      IMPLIBEXT=".$LIBEXT"
      CP_A=`which cp`" -pR"
      AC_CHECK_PROG([LDD], [otool], [otool -L])
      LDD_AWK_COL="1"
      ;;
   *)
      isWindowsOS="no"
      LIBEXT="so"
      DLLLIBEXT=".$LIBEXT"
      IMPLIBEXT=".$LIBEXT"
      CP_A=`which cp`" -a"
      AC_CHECK_PROG([LDD], [ldd], [ldd])
      LDD_AWK_COL="3"
      ;;
esac

AM_CONDITIONAL([BUILD_LINUX],        [test -n "`uname | grep Linux`"])
AM_CONDITIONAL([BUILD_DARWIN],       [test -n "`uname | grep Darwin`"])
AM_CONDITIONAL([COND_BUILD_WIN32],   [test x"$isWindowsOS" = x"yes"])
AM_CONDITIONAL([COND_BUILD_MOBILE],  [test x"$isMobile"    = x"yes"])
AM_CONDITIONAL([COND_BUILD_ANDROID], [test x"$isAndroid"   = x"yes"])
AM_CONDITIONAL([COND_BUILD_IOS],     [test x"$isIOS"       = x"yes"])



# --- 32 or 64 bit?
AC_CHECK_SIZEOF(void*)
if test x"$ac_cv_sizeof_voidp" = x"8"; then
   AC_DEFINE_UNQUOTED(HAS_64BIT, "1",   [64bit environment])
#   LIB3264=lib64
else
   AC_DEFINE_UNQUOTED(HAS_32BIT, "1",   [32bit environment])
#   LIB3264=lib
fi

# --- Check endian type
AC_C_BIGENDIAN(
[ AC_DEFINE_UNQUOTED(HAS_BIG_ENDIAN, 1, [System byte order]) ], [],
[ AC_MSG_ERROR([Cannot determine byte order of this system.])])

# --- Define variables in header file
AC_DEFINE_UNQUOTED(UNAME,         "$UNAME",         [Output of uname(1)])
AC_DEFINE_UNQUOTED(PROCESSOR,     "$PROCESSOR",     [CPU type])
AC_DEFINE_UNQUOTED(WHOAMI,        "$WHOAMI",        [User who compiled package])
AC_DEFINE_UNQUOTED(GDB,           "$GDB",           [Debugger])
AC_DEFINE_UNQUOTED(MEMTRACER,     "$MEMTRACER",     [Memory tracer tool])
AC_DEFINE_UNQUOTED(DIST_OS,       "$DIST_OS",       [Operating system ID])
AC_DEFINE_UNQUOTED(DIST_NAME,     "$DIST_NAME",     [Distribution ID])
AC_DEFINE_UNQUOTED(DIST_VERSION,  "$DIST_VERSION",  [Distribution Version])
AC_DEFINE_UNQUOTED(DIST_VERSION_L, $DIST_VERSION_L ,[Distribution Version Integer Number])
AC_DEFINE_UNQUOTED(DIST_MAJOR,    "$DIST_MAJOR",    [Distribution Major Version])
AC_DEFINE_UNQUOTED(DIST_MINOR,    "$DIST_MINOR",    [Distribution Minor Version])
AC_DEFINE_UNQUOTED(DIST_PATCH,    "$DIST_PATCH",    [Distribution Patch Version])
AC_DEFINE_UNQUOTED(DIST_REPO_VER, "$DIST_REPO_VER", [Distribution Repository Version])
AC_DEFINE_UNQUOTED(PKG_TYPE,      "$PKG_TYPE",      [Package type])
AC_DEFINE_UNQUOTED(REPOPROC,      "$REPOPROC",      [Repository processor id])
AC_DEFINE_UNQUOTED(SXBINTARGET,   "$SXBINTARGET",   [Binary target directory])
AC_DEFINE_UNQUOTED(SXACCELERATE_SRC,"$ABS_TOP_SRCDIR", [Path to source tree])
AC_DEFINE_UNQUOTED(SXACCELERATE_BUILD,"$ABS_TOP_BUILDDIR", [Path to build tree])

# --- Substitute variables in makefile
AC_SUBST([UNAME])
AC_SUBST([PROCESSOR])
AC_SUBST([WHOAMI])
AC_SUBST([TODAY])
AC_SUBST([LIBEXT])
AC_SUBST([DLLLIBEXT])
AC_SUBST([IMPLIBEXT])
AC_SUBST([CP_A])
AC_SUBST([GDB])
AC_SUBST([LDD])
AC_SUBST([LDD_AWK_COL])
AC_SUBST([SXLIBTOOL])
AC_SUBST([DIST_OS])
AC_SUBST([DIST_NAME])
AC_SUBST([DIST_VERSION])
AC_SUBST([DIST_MAJOR])
AC_SUBST([DIST_MINOR])
AC_SUBST([DIST_PATCH])
AC_SUBST([DIST_REPO_VER])
AC_SUBST([PKG_TYPE])
AC_SUBST([REPOPROC])
AC_SUBST([SXBINTARGET])
])
