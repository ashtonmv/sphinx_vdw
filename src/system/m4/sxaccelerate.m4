AC_DEFUN([SX_CHECK_SXACCELERATE], [

   SX_ARG_WITH(  [sxaccelerate], [.], [SXACCELERATE], [],
                 [Absolute path to the SxAccelerate/ top level folder])

   sxaccel_build=true
   if test x"$ac_cv_with_sxaccelerate" = x; then
      if test -f ../sxaccelerate/src/system/$ac_cv_with_bintarget/sxconfig.dat; then
         # --- found easy build: sphinx + sxaccelerate in combined source tree
         . ../sxaccelerate/src/system/$ac_cv_with_bintarget/sxconfig.dat
      else
         AC_MSG_ERROR([Could not find SxAccelerate Framework])
      fi
   else
      if test -f ${ac_cv_with_sxaccelerate}/system/$ac_cv_with_bintarget/sxconfig.dat; then
         # --- found separate sxaccelerate build tree
         . ${ac_cv_with_sxaccelerate}/system/$ac_cv_with_bintarget/sxconfig.dat
      else
         if test -f ${ac_cv_with_sxaccelerate}/include/SxString.h; then :; else
            AC_MSG_ERROR([Could not find SxString.h in ${ac_cv_with_sxaccelerate}/include])
         fi
         if test -f ${ac_cv_with_sxaccelerate}/lib/libsxutil.la; then :; else
            AC_MSG_ERROR([Could not find libsxutil.la in ${ac_cv_with_sxaccelerate}/lib])
         fi
         sxaccel_build=false
      fi
      if ${sxaccel_build} ; then : ; else
        abs_sxaccelerate=$(cd $ac_cv_with_sxaccelerate;pwd)
        ac_cv_with_sxaccelerate_src=$abs_sxaccelerate
        ac_cv_with_sxaccelerate_build=$abs_sxaccelerate
        CPPFLAGS="$CPPFLAGS -I ${abs_sxaccelerate}/include"
      fi
   fi

   AC_SUBST([SXACCELERATE_BUILD],[$ac_cv_with_sxaccelerate_build])
   AC_SUBST([SXACCELERATE_SRC],[$ac_cv_with_sxaccelerate_src])
   AC_SUBST([SXACCELERATE],[$ac_cv_with_sxaccelerate])
   AM_CONDITIONAL([COND_BUILD_ACCELERATE],[$sxaccel_build])

   # --- older versions of autoconf don't provide abs_top_srcdir/abs_top_builddir
   if test x"$ac_abs_top_srcdir" = x; then
      abs_top_srcdir=`cd $top_srcdir && pwd`
      AC_SUBST(abs_top_srcdir)
   fi
   if test x"$ac_abs_top_builddir" = x; then
      abs_top_builddir=`cd $top_builddir && pwd`
      AC_SUBST(abs_top_builddir)
   fi

])
