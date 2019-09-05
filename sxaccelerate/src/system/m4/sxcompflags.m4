AC_DEFUN([SX_COMPFLAGS], [


SX_ARG_ENABLE([debug],    [DBG_MODE], [yes], 
              [compile package in DEBUG mode - RELEASE mode otherwise])
SX_ARG_ENABLE([logging],  [USE_SX_LOG], [yes],
              [compile package with logging support])
SX_ARG_ENABLE([profile],  [PRO_MODE], [no], 
              [compile package in the profile mode])
SX_ARG_ENABLE([openmp],   [USE_OPENMP], [no], 
              [compile package with symmetric multiprocessing support])
SX_ARG_ENABLE([mpi],      [USE_MPI], [no], 
              [compile package with Message-Passing Interface support])
# --- control CPU-specific compilation in some libraries, notably sxmath
SX_ARG_ENABLE([cpuspecific], [MACHINE_SPECIFIC], [$SX_PARAM_CPUSPECIFIC],
              [compile code to exploit all features of the current CPU])
SX_ARG_ENABLE([errorlimit], [GCCERRORLIMIT], [1],
              [limit number of compiler errors before aborting])
AC_ARG_VAR([CXX_M_ARCH],
           AS_HELP_STRING([compiler flags to select CPU features, e.g. -march=...]))

AM_CONDITIONAL([COND_MPI], [test x"$ac_cv_enable_mpi" = x"yes"])

# --- GPU support
SX_ARG_ENABLE([gpu],      [USE_GPU], [no], 
              [compile package with NVIDIA CUDA support])
AM_CONDITIONAL([COND_GPU], [test x"$ac_cv_enable_gpu" = x"yes"])
SX_PARAM_CUDA="/usr/local/cuda"
SX_ARG_WITH(  [cuda], [.], [CUDA], [$SX_PARAM_CUDA],
              [absolute path to the NVIDIA CUDA toolkit])
# ---
AC_MSG_CHECKING(operating system)


AC_MSG_NOTICE([Determine compiler flags])

# --- Search for C/C++ compilers
if test x"$DIST_OS" = x"Darwin"; then
   # MacOS: prefer gcc over clang
   AC_PROG_CC(clang gcc)              # GNU C compiler for dependencies
   if test "$ac_cv_enable_debug" = "yes"; then
   AC_PROG_CXX(clang++ g++)              # GNU C++ compiler for DEBUG mode
   else
   AC_PROG_CXX(clang++ g++ icc ipcp xlC aCC)  # C++ compilers for RELEASE mode
   fi
   AC_PROG_OBJCXX(clang++ g++)           # Objective C++ compiler
else
   # Linux: prefer gcc over clang
   AC_PROG_CC(gcc clang)              # GNU C compiler for dependencies
   if test "$ac_cv_enable_debug" = "yes"; then
   AC_PROG_CXX(g++ clang++)              # GNU C++ compiler for DEBUG mode
   else
   AC_PROG_CXX(g++ clang++ icc icpc xlC aCC)  # C++ compilers for RELEASE mode
   fi
   AC_PROG_OBJCXX(g++ clang++)           # Objective C++ compiler
fi

cxxname=`echo "$CXX" | sed -e's#.*/##'`
if test x"$cxxname" != x"g++" -a $ac_cv_cxx_compiler_gnu = yes ; then
   # CXX not "g++", but seems to be a GNU C++ compiler. Try to verify...
   if $CXX --version | head -n 1 | grep -q -e 'GCC' -e "g++" ; then
      cxxname="g++"
      AC_MSG_NOTICE([$CXX has been identified as GNU C++ (g++)])
   elif $CXX --version | head -n 1 | grep -q -e '^icpc (ICC)' ; then
      AC_MSG_NOTICE([$CXX has been identified as Intel C++ (g++)])
      cxxname="icpc"
   else
      AC_MSG_NOTICE([$CXX pretends to be GCC, but fails to prove it in --version])
   fi
fi

if test x"$ac_cv_enable_cpuspecific" = x"" ; then
   if test x"${CXX_M_ARCH:-empty}" = x"empty" ; then
     ac_cv_enable_cpuspecific="no"
   else
     ac_cv_enable_cpuspecific="yes"
   fi
fi
if test x"$ac_cv_enable_cpuspecific" = x"no" ; then
   if test x"${CXX_M_ARCH:-empty}" = x"empty" ; then : ; else
      AC_MSG_ERROR([Conflict: --disable-cpuspecific, but CXX_M_ARCH='$CXX_M_ARCH'])
   fi
   CXX_M_ARCH=""
fi

# --- Set OS dependent flags
if test x"$DIST_OS" = x"Linux"; then
   AC_MSG_RESULT(Linux)
   AC_DEFINE(LINUX, "1", [Setup Linux Compilation])
   MAKE_SILENT_ARGS="-s --no-print-directory"
   MODEL=`cat /proc/cpuinfo | grep 'model name' | uniq`
   FLAGS=`cat /proc/cpuinfo | grep 'flags' | uniq`
elif test x"$DIST_OS" = x"Darwin"; then
   AC_MSG_RESULT(MacOSX)
   AC_DEFINE(MACOSX, "1", [Setup MACOSX Compilation])
   MAKE_SILENT_ARGS="-s"
   LIBS="$LIBS -framework CoreFoundation"
elif test x"$DIST_OS" = x"HP-UX"; then
   AC_MSG_RESULT(HP-UX)
   AC_DEFINE(HPUX, "1", [Setup HPUX Compilation])
   MAKE_SILENT_ARGS="-s --no-print-directory"
elif test x"$DIST_OS" = x"SunOS"; then
   AC_MSG_RESULT(Athlon XP)
   AC_DEFINE(SUNOS, "1", [Setup SunOS Compilation])
   MAKE_SILENT_ARGS="-s --no-print-directory"
else   
   AC_MSG_ERROR(unknown)
fi   



# --- Set compiler flags for RELEASE mode
found_compiler=yes
case "$cxxname" in
   g++)  # GNU C/C++ Compiler
      GCC_M_ARCH="${CXX_M_ARCH--march=native}"
      AC_MSG_CHECKING([GCC_M_ARCH])
      AC_MSG_RESULT([$GCC_M_ARCH])
      CXX_VER=`$CXX --version | head -1`
      # -----------------------------------------------------------------------
      CXX_REL="-Wall -O3 $GCC_M_ARCH"
      CXX_REL="$CXX_REL -fpermissive"
      #
      #bl: needed for vector conversion in Accelerate Framwork in MacOsX 10.9
      if test ! -z "`uname | grep Darwin`"; then
         CXX_REL="$CXX_REL -flax-vector-conversions"
      fi
      CXX_REL="$CXX_REL -DNDEBUG -pipe"
      CXX_REL="$CXX_REL -Wno-variadic-macros"
      # -----------------------------------------------------------------------
      CXX_DBG="-g -O2 -W -Wall -Wcast-align -D_REENTRANT"
      if test x"$ac_cv_enable_errorlimit" = x"no" ; then : ;
      elif test x"$ac_cv_enable_errorlimit" = x"yes" ; then
        CXX_DBG="$CXX_DBG -fmax-errors=1"
      else
        CXX_DBG="$CXX_DBG -fmax-errors=$ac_cv_enable_errorlimit"
      fi
      CXX_DBG="$CXX_DBG -Wcast-qual -Wchar-subscripts"
      CXX_DBG="$CXX_DBG -Wpointer-arith"
      CXX_DBG="$CXX_DBG -Wshadow"
      CXX_DBG="$CXX_DBG -Wwrite-strings"
      CXX_DBG="$CXX_DBG -pipe"
      CXX_DBG="$CXX_DBG -Wno-variadic-macros"
      #bl: needed for vector conversion in Accelerate Framwork in MacOsX 10.9
      if test ! -z "`uname | grep Darwin`"; then
         CXX_DBG="$CXX_DBG -flax-vector-conversions"  
      fi
      CXX_PEDANTIC="-std=c++14 -pedantic -Wconversion -Wredundant-decls -Wno-long-long"
      CXX_WARNING_IS_ERROR="-Werror -Wfatal-errors"
      LEX_CXXFLAGS="-Wno-conversion -Wno-redundant-decls -Wno-unused"
      LEX_CXXFLAGS="$LEX_CXXFLAGS -Wno-unused-parameter -Wno-sign-compare"
      # -----------------------------------------------------------------------
      CXX_M_ARCH=$GCC_M_ARCH
      # -----------------------------------------------------------------------
      CXX_OPENMP="-fopenmp"
      LDC_OPENMP="-fopenmp"
      # -----------------------------------------------------------------------
      CXX_PRO="-g -pg -DNDEBUG"
      LDF_PRO=""
      # -----------------------------------------------------------------------
      AR_FLAGS="ruv"
      # -----------------------------------------------------------------------
      ;;
   clang++)  # GNU C/C++ Compiler
      GCC_M_ARCH="${CXX_M_ARCH--march=native}"
      AC_MSG_CHECKING([GCC_M_ARCH])
      AC_MSG_RESULT([$GCC_M_ARCH])
      CXX_VER=`$CXX --version | head -1`
      # -----------------------------------------------------------------------
      CXX_REL="-Wall -O3 $GCC_M_ARCH"
      CXX_REL="$CXX_REL -fpermissive"
      #
      #bl: needed for vector conversion in Accelerate Framwork in MacOsX 10.9
      if test ! -z "`uname | grep Darwin`"; then
         CXX_REL="$CXX_REL -flax-vector-conversions"
      fi
      CXX_REL="$CXX_REL -DNDEBUG -pipe"
      CXX_REL="$CXX_REL -Wno-variadic-macros"
      # -----------------------------------------------------------------------
      CXX_DBG="-g -O2 -W -Wall -Wcast-align -D_REENTRANT"
#      if test x"$ac_cv_enable_errorlimit" = x"no" ; then : ;
#      elif test x"$ac_cv_enable_errorlimit" = x"yes" ; then
#        CXX_DBG="$CXX_DBG -fmax-limit=1"
#      else
#        CXX_DBG="$CXX_DBG -fmax-limit=$ac_cv_enable_errorlimit"
#      fi
      CXX_DBG="$CXX_DBG -Wcast-qual -Wchar-subscripts"
      CXX_DBG="$CXX_DBG -Wpointer-arith"
      CXX_DBG="$CXX_DBG -Wshadow"
      CXX_DBG="$CXX_DBG -Wwrite-strings"
      CXX_DBG="$CXX_DBG -pipe"
      CXX_DBG="$CXX_DBG -Wno-variadic-macros"
      #bl: needed for vector conversion in Accelerate Framwork in MacOsX 10.9
      if test ! -z "`uname | grep Darwin`"; then
         CXX_DBG="$CXX_DBG -flax-vector-conversions"  
      fi
      CXX_PEDANTIC="-std=c++14 -pedantic -Wconversion -Wredundant-decls -Wno-long-long"
      CXX_WARNING_IS_ERROR="-Werror -Wfatal-errors"
      LEX_CXXFLAGS="-Wno-conversion -Wno-redundant-decls -Wno-unused"
      LEX_CXXFLAGS="$LEX_CXXFLAGS -Wno-unused-parameter -Wno-sign-compare"
      # -----------------------------------------------------------------------
      CXX_M_ARCH=$GCC_M_ARCH
      # -----------------------------------------------------------------------
      CXX_OPENMP="-fopenmp"
      LDC_OPENMP="-fopenmp"
      # -----------------------------------------------------------------------
      CXX_PRO="-g -pg -DNDEBUG"
      LDF_PRO=""
      # -----------------------------------------------------------------------
      AR_FLAGS="ruv"
      # -----------------------------------------------------------------------
      ;;
   icpc) # Intel C/C++ Compiler
      ICC_M_ARCH="${CXX_M_ARCH--xHost}"
      CXX_VER=`$CXX --version | head -1`
      # -----------------------------------------------------------------------
      CXX_REL="-O3 $ICC_M_ARCH"
      CXX_REL="$CXX_REL -DNDEBUG"
      # -----------------------------------------------------------------------
      CXX_DBG="-O0 -no-ip -g -Wall -Wcheck -w1 $ICC_M_ARCH"
      # -----------------------------------------------------------------------
      CXX_OPENMP="-qopenmp"
      LDC_OPENMP="-qopenmp"
      # -----------------------------------------------------------------------
      CXX_PRO="-g -O -DNDEBUG"
      LDF_PRO=""
      # -----------------------------------------------------------------------
      CXX_M_ARCH=$ICC_M_ARCH
      # -----------------------------------------------------------------------
      AR_FLAGS="ruv"
      # -----------------------------------------------------------------------
      ;;
   aCC)  # HP-UX aCC Compiler
      CXX_VER=`$CXX -V`
      # -----------------------------------------------------------------------
      CXX_REL="-ext -mt -w +p -AA -D_HPUX"
      CXX_REL="$CXX_REL -D_XOPEN_SOURCE_EXTENDED" # force aCC to deal with UNIX sockets
      CXX_REL="$CXX_REL +W749"                    # suppress 'reinterpret_cast' warning
      CXX_REL="$CXX_REL -D__HP_NO_MATH_OVERLOADS" # ignore HP-UX cmath overloading
      CXX_REL="$CXX_REL +O2 +DA2.0W -DNDEBUG"
      CXX_REL="$CXX_REL -DUSE_VECLIB -DUSE_VECLIB_FFT -D_REENTRANT"
      CXX_REL="$CXX_REL -D_POSIX_C_SOURCE=199506L"
      CXX_REL="$CXX_REL -I/usr/local/include/"
      CXX_REL="$CXX_REL -I/usr/local/include/pa20_64"
      # -----------------------------------------------------------------------
      CXX_OPENMP=""
      LDC_OPENMP=""
      # -----------------------------------------------------------------------
      CXX_PRO="-g -O -DNDEBUG"
      LDF_PRO=""
      # -----------------------------------------------------------------------
      AR_FLAGS="ruvl"
      # -----------------------------------------------------------------------
      ;;
   *)    # --- No Compiler found
      found_compiler=no
esac
if test $found_compiler = no -a "$ac_test_CXXFLAGS" != "set"; then
   cat << ERRMSG
$SEPARATOR
ERROR: Compiler not supported.
$SEPARATOR
       This package does not support the compiler "$CXX". 
       You can either provide a C++ compiler supported by SPHInX or
       set the compiler flags (CXXFLAGS) by hand.
$SEPARATOR
ERRMSG
   exit 1
fi

# --- GPU: NVIDIA CUDA compiler wrapper and flags
if test x"${ac_cv_enable_gpu}"  = x"yes"; then
   AC_MSG_CHECKING(CUDA toolkit)
   # ---
   NVCC="$ac_cv_with_cuda/bin/nvcc"
   CUDA_CFLAGS="-O3 -I$ac_cv_with_cuda/include"
   CUDA_LIBS="-lcudart"
   CUDA_LDFLAGS="-L$ac_cv_with_cuda/lib64"
   # --- check if nvcc is actually available
   SX_REQUIRED(NVCC, nvcc)
   # ---
   NVCC_VER=`$NVCC --version | tail -1`
   AC_MSG_RESULT($NVCC_VER)
fi

CXXLD="$CXX"



# --- Set up compiler flags according to OS, compiler, and mode   
if test  "$ac_test_CXXFLAGS" != "set"; then
   # remove standard CXXFLAGS provided by autoconf
   CXXFLAGS=""
   if test "$ac_cv_enable_debug" = "yes"; then
      AM_CXXFLAGS="$CXX_DBG"
   else
      AM_CXXFLAGS="$CXX_REL"
   fi
   if test "$ac_cv_enable_profile" = "yes"; then
      AM_CXXFLAGS="$AM_CXXFLAGS $CXX_PRO"
   fi
   if test x"${ac_cv_enable_openmp}" = x"yes"; then
      AM_CXXFLAGS="$AM_CXXFLAGS $CXX_OPENMP"
   fi
   if test x"${ac_cv_enable_mpi}" = x"yes"; then
      AM_CXXFLAGS="$AM_CXXFLAGS -DMPICH_IGNORE_CXX_SEEK"
   fi

   AC_SUBST(AM_CXXFLAGS,  "$AM_CXXFLAGS")
   AC_SUBST(CXX_PEDANTIC, "$CXX_PEDANTIC")
   AC_SUBST(LEX_CXXFLAGS, "$LEX_CXXFLAGS")
   AC_SUBST(CXX_WARNING_IS_ERROR,"$CXX_WARNING_IS_ERROR")
   CXXFLAGS="$AM_CXXFLAGS"
   # bug in libtool: postdeps_CXX is set from CFLAGS rather than CXXFLAGS
   CFLAGS="$AM_CXXFLAGS"
else
   AC_SUBST(AM_CXXFLAGS, "$CXXFLAGS")
fi


# --- check extra math functions (needs -lm)
ORIG_LIBS="$LIBS"
LIBS="$LIBS -lm"
AC_CHECK_FUNCS(round lround sincos)
LIBS="$ORIG_LIBS"


# --- check OS dependent libraries
AC_CHECK_LIB([dl],[dlopen], [LIBS="$LIBS -ldl"])
AC_CHECK_LIB([pthreadGCE2], [pthread_create], [LIBS="$LIBS -lpthreadGCE2"],
   [AC_CHECK_LIB([pthread], [pthread_create], [LIBS="$LIBS -lpthread"])]
)

if test x"${ac_cv_enable_openmp}" = x"yes"; then
   LIBS="$LIBS $LDC_OPENMP"
fi

AC_MSG_CHECKING([for suitable compiler flags])
AC_MSG_RESULT([$AM_CXXFLAGS])


# --- SxConfig.h
AC_DEFINE_UNQUOTED(CXX,         "$CXX",         [C++ compiler])
AC_DEFINE_UNQUOTED(CXXFLAGS,    "$AM_CXXFLAGS", [C++ compiler flags])
AC_DEFINE_UNQUOTED(CXXVERSION,  "$CXX_VER",     [C++ compiler version])
AC_DEFINE_UNQUOTED(F77,         "$F77",         [FORTRAN77 compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      "$FFLAGS",      [FORTRAN77 compiler flags])


# --- Makefile & Co
AC_SUBST(CC)
AC_SUBST(CXXLD)
AC_SUBST(MAKE_SILENT_ARGS)


# --- GPU/CUDA compiler and flags
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_LDFLAGS)
AC_SUBST(NVCC)


])
