# --- Define macro for --enable-FEATURE
#     1: FEATURE name,          addons for --enable-atlas
#     2: variable name,         USE_ATLAS
#     3: default value,         yes, no
#     4: help message,          compile with ATLAS support
AC_DEFUN([SX_ARG_ENABLE],
[AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1], [$4 (default=$3)]),
ac_cv_enable_$1=$enableval, ) 
AC_CACHE_CHECK(whether to enable $1, ac_cv_enable_$1, ac_cv_enable_$1=$3)
AC_SUBST([$2],$ac_cv_enable_$1,[myDescription])
if test "${ac_cv_enable_$1}" = "yes"; then
  AC_DEFINE([$2],[],[$4])
fi
])

# --- Define macro for --with-FEATURE
#     1: FEATURE name,          addons for --with-addons
#     2: available directory,   add-ons for ../src/add-ons
#     3: variable name,         WITH_ADDONS
#     4: default value,         yes, no
#     5: help message,          compile with add-ons
AC_DEFUN([SX_ARG_WITH],
[AC_ARG_WITH([$1],
AC_HELP_STRING([--with-$1], [$5 (default=$4)]),
ac_cv_with_$1=$withval, )
if test "$4" = "opt"; then
   if test -d "[$2]"; then
     with_$1="yes"
   else
     with_$1="no"
   fi
else
   with_$1="$4"
fi
AC_CACHE_CHECK(whether to use $1, ac_cv_with_$1, ac_cv_with_$1=$with_$1)
if test "${ac_cv_with_$1}" = "yes"; then
  AC_DEFINE([$3],[],[$5])
  if test -d "[$2]"; then
     AC_DEFINE([$3],[],[$5])
  else
     echo
     echo "Module [$1] is not available. Please run"
     echo "  configure --without-[$1]"
     exit 1
  fi
fi
AC_SUBST([$3],$ac_cv_with_$1)
])


# --- Define macro for required packages
#     1: variable               CXX
#     2: package name           GNU C++ compiler
AC_DEFUN([SX_REQUIRED],
[ if test -z "$[$1]"; then
     cat <<MSGEND
$SEPARATOR
ERROR: No '[$2]' found.
$SEPARATOR
       The $PACKAGE_NAME package cannot be compiled without [$2].
       Please make sure that the corresponding package is installed
       and located in the search path.
$SEPARATOR
MSGEND
     exit 1
   fi
])

# --- Define macro for crosschecking modules and headers
#     1:    feature name        enable_tb
#        or module name         with_tb
#        or program name        proc_ac_ct_LATEX
#     2:    headers             mysql/mysql.h => header_mysql_mysql_h
#        or functions           sqrt          => func_sqrt
AC_DEFUN([SX_VALIDATE_MODULE],
[
   if test "${ac_cv_$1}" = "yes"; then
      if test ! "${ac_cv_$2}" = "yes"; then
         id=`echo $2 | sed -e 's/_h$/.h/g;s/^header_//g;s/^func_//g;s/prog_ac_ct_//g;s/_/\//g'`
         mod=`echo $1 | sed -e 's/enable_//g;s/with_//g'`
         echo
         echo "$id not found. Cannot support module $mod"
         echo "Please deactive it by running"
         echo "  configure --without-[$1]"
         exit 1
      fi
   fi
])

# --- Define macro for call configure for SUBDIR
#     with extra configure settings as requested
#     (similar to AC_CONFIG_SUBDIRS, but with specific configure flags)
#     1: SUBDIR name,         e.g. sxaccelerate
#     2: extra flags,         e.g. --with-sxmath
#     if prefix is not set, it is set to current path
AC_DEFUN([SX_CONFIG_SUBDIR],[
if false ; then
  # never go here, but trigger autotools recursions
  AC_CONFIG_SUBDIRS([$1])
fi
if echo "$ac_configure_args" | grep -q -- "--prefix" ; then : ; else
   if test x"$prefix" = x"NONE" ; then
      prefix=$(pwd)
   fi
   ac_configure_args="--prefix=$prefix $ac_configure_args"
fi
sx_topsrc=$(cd $srcdir && pwd)
mkdir -p $1
AC_MSG_NOTICE([configure $1 with $sx_topsrc/$1/configure --disable-option-checking $2 $ac_configure_args])
(cd $1 && eval "$sx_topsrc/$1/configure --disable-option-checking $2 $ac_configure_args" || AC_MSG_ERROR([Failed to configure $1]))
])

