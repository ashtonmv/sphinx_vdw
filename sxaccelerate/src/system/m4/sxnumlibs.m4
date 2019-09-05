AC_DEFUN([SX_IF_CHECKNUMLIBS], [
if test x"$ac_cv_enable_numlibschecks" = x"yes"; then
   $1
fi
])
AC_DEFUN([SX_CHECK_NUMLIBS], [


AC_MSG_NOTICE([Setting up numlibs])


ORIG_LIBS=$LIBS
LIBS="-lm $LIBS"

# --- some numeric libs require FORTRAN
AC_CHECK_LIB(gfortran, rand)
AC_CHECK_LIB(gfortran, _gfortran_copy_string, 
             [AC_DEFINE([HAVE_GFORTRAN_COPY_STRING],[1],
             [Define to 1 if gfortran does not come with gcc bug 33647])])
SX_BLAS_LIBS="$SX_BLAS_LIBS $LIBS"

SX_PARAM_ATLAS="no"
SX_PARAM_ACML="no"
SX_PARAM_MKL="no"
SX_PARAM_GOTO="no"
SX_PARAM_MKLPATH="/opt/intel/mkl"
SX_PARAM_FFTW="no"
SX_PARAM_ACMLFFT="no"
SX_PARAM_MKLFFT="no"

if test x"$ac_cv_with_sxmath" = x"yes"; then :; else
   SX_PARAM_ATLAS="no"
   SX_PARAM_FFTW="no"
fi

configFile=$1
if test -n "$configFile"; then
   if test -f "$configFile"; then
      . $configFile
   fi
fi


SX_ARG_WITH(  [numlibs], [.], [NUMLIBS], [$SX_PARAM_NUMLIBS],
              [absolute paths to top level folder for external libraries])
SX_ARG_ENABLE([numlibschecks], [NUMLIBS_CHECKS], [no], 
              [check that numlibs can be used])
if echo x"$ac_cv_enable_numlibschecks" | grep -q '^x *-' ; then
  AC_MSG_NOTICE([Appending $ac_cv_enable_numlibschecks to linker flags for checking numlibs])
   LIBS="$ac_cv_enable_numlibschecks $LIBS"
   ac_cv_enable_numlibschecks="yes"
fi
SX_ARG_ENABLE([atlas],    [USE_ATLAS], [$SX_PARAM_ATLAS], 
              [compile package with ATLAS support])
SX_ARG_ENABLE([acml],     [USE_ACML], [$SX_PARAM_ACML], 
              [compile package with AMD Core Math Library support])
SX_ARG_WITH(  [mklpath], [.], [MKLPATH], [$SX_PARAM_MKLPATH],
              [absolute path to the Intel MKL top level folder])
SX_ARG_ENABLE([mkl],      [USE_INTEL_MKL], [$SX_PARAM_MKL], 
              [compile package with Intel Math Kernel Library support])
SX_ARG_ENABLE([goto],     [USE_GOTO], [$SX_PARAM_GOTO], 
              [compile package with GotoBLAS support])
SX_ARG_ENABLE([fftw],     [USE_FFTW],  [$SX_PARAM_FFTW], 
              [compile package with FFTW support])
SX_ARG_ENABLE([mklfft],   [USE_MKL_FFT],  [$SX_PARAM_MKLFFT], 
              [compile package with MKL's FFT support])
SX_ARG_ENABLE([acmlfft],  [USE_ACML_FFT],  [$SX_PARAM_ACMLFFT], 
              [compile package with ACML's FFT support])
SX_ARG_ENABLE([mpi],      [USE_MPI], [no], 
              [compile package with Message-Passing Interface support])
SX_ARG_ENABLE([netcdf4], [USE_NETCDF4], [no], 
              [compile package with support for NetCDF4 IO])
SX_ARG_ENABLE([parnetcdf4], [USE_PARALLEL_NETCDF4], [no], 
              [compile package with support for parallel NetCDF4 IO])
SX_ARG_ENABLE([openmp],   [USE_OPENMP], [no], 
              [compile package with OpenMP support])
SX_ARG_ENABLE([hdf5], [USE_HDF5], [no], [compile package with HDF5 support])
SX_ARG_ENABLE([pcre2], [USE_PCRE2], [auto], [compile package with PCRE2 support])


# --- check mutual exclusions
if test x"${ac_cv_enable_acml}"  = x"yes" \
     -a x"${ac_cv_enable_atlas}" = x"yes"; then
   AC_MSG_ERROR([Cannot use ACML and ATLAS simultaneously.])
fi
if test x"${ac_cv_enable_acml}"  = x"yes" \
     -a x"${ac_cv_enable_mkl}" = x"yes"; then
   AC_MSG_ERROR([Cannot use ACML and MKL simultaneously.])
fi
if test x"${ac_cv_enable_acml}"  = x"yes" \
     -a x"${ac_cv_enable_goto}" = x"yes"; then
   AC_MSG_ERROR([Cannot use ACML and GotoBLAS simultaneously.])
fi
if test x"${ac_cv_enable_mkl}"  = x"yes" \
     -a x"${ac_cv_enable_goto}" = x"yes"; then
   AC_MSG_ERROR([Cannot use MKL and GotoBLAS simultaneously.])
fi
if test x"${ac_cv_enable_atlas}"  = x"yes" \
     -a x"${ac_cv_enable_goto}" = x"yes"; then
   AC_MSG_ERROR([Cannot use ATLAS and GotoBLAS simultaneously.])
fi
if test x"${ac_cv_enable_mkl}"  = x"yes" \
     -a x"${ac_cv_enable_atlas}" = x"yes"; then
   AC_MSG_ERROR([Cannot use MKL and ATLAS simultaneously.])
fi
if test x"${ac_cv_enable_acmlfft}" = x"yes" \
     -a x"${ac_cv_enable_fftw}"    = x"yes"; then
   AC_MSG_ERROR([Cannot use ACML's FFT and FFTW interface simultaneously.])
fi
if test x"${ac_cv_enable_parnetcdf4}" = x"yes" \
     -a x"${ac_cv_enable_mpi}"              = x"no"; then
   AC_MSG_ERROR([Parallel NetCDF4 IO requires MPI.])
fi
if test x"${ac_cv_enable_parnetcdf4}" = x"yes" \
     -a x"${ac_cv_enable_netcdf4}"          = x"no"; then
   AC_MSG_NOTICE([Enabling NetCDF4 also for serial IO.])
   ac_cv_enable_netcdf4="yes"
fi

# Copyrights and licenses
SXCOPY="| S/PHI/nX utilizes the following 3rd-party libraries:\n"
if test x"${ac_cv_enable_acml}"  = x"yes"; then
   SXCOPY="${SXCOPY}|    - ACML            http://developer.amd.com/acml\n"
fi
if test x"${ac_cv_enable_atlas}"  = x"yes"; then
   SXCOPY="${SXCOPY}|    - Atlas           http://math-atlas.sourceforge.net\n"
fi
if test x"${ac_cv_enable_mkl}"  = x"yes"; then
   SXCOPY="${SXCOPY}|    - Intel MKL       http://software.intel.com/en-us/intel-mkl\n"
fi
if test x"${ac_cv_enable_fftw}"  = x"yes"; then
   SXCOPY="${SXCOPY}|    - FFTW            http://www.fftw.org\n"
fi
SXCOPY="${SXCOPY}|    - Flex            http://flex.sourceforge.net\n"
SXCOPY="${SXCOPY}|    - NetCDF          http://www.unidata.ucar.edu/software/netcdf\n"
AC_DEFINE_UNQUOTED(SXCOPYRIGHT, "$SXCOPY",[Copyright information])

if test x"$ac_cv_with_numlibs" = x"no"; then
   ac_cv_with_numlibs=""
fi

# add numlibs top levels to library and include search path
for dir in `echo "$ac_cv_with_numlibs" | sed -e's/:/ /g'` ; do
   if test -d $dir/include ; then
      AC_MSG_NOTICE([Adding include path $dir/include])
      CPPFLAGS="$CPPFLAGS `cd $dir/include && pwd | sed -e's/^/-I/'`"
   else
      AC_MSG_NOTICE([No $dir/include])
   fi
   if test -d $dir/lib; then
      AC_MSG_NOTICE([Adding library path $dir/lib])
      LDFLAGS="$LDFLAGS `cd $dir/lib && pwd | sed -e's/^/-L/'`"
   else
      AC_MSG_NOTICE([No $dir/lib])
   fi
done
   
# MKL include and library paths
if test x"${ac_cv_enable_mkl}"  = x"yes"; then
   if test x"$ac_cv_with_mklpath" = x"null"; then
      # we require that the path to MKL is explicitly given
      AC_MSG_ERROR([Cannot find Intel MKL.  Use the flag --with-mklpath to point to an MKL installation.])
   fi
   AC_CHECK_FILE([$ac_cv_with_mklpath/include/mkl.h], [],
                 [AC_MSG_ERROR([Cannot find Intel MKL.  Use the flag --with-mklpath to point to an MKL installation.])]
   )
   if test -d "$ac_cv_with_mklpath" ; then
     ac_cv_with_mklpath=`cd $ac_cv_with_mklpath; pwd`
   fi
   CPPFLAGS="$CPPFLAGS -I$ac_cv_with_mklpath/include"

   # find the actual lib path
   if test -d ${ac_cv_with_mklpath}/lib/em64t ; then
      mkllibpath=`cd ${ac_cv_with_mklpath}/lib/em64t && pwd`
   elif test -d ${ac_cv_with_mklpath}/lib/intel64 ; then
      mkllibpath=`cd ${ac_cv_with_mklpath}/lib/intel64 && pwd`
   elif test -f ${ac_cv_with_mklpath}/lib/libmkl_sequential.so ; then
      mkllibpath=`cd ${ac_cv_with_mklpath}/lib && pwd`
   else
      AC_MSG_ERROR([Cannot find Intel MKL library path. Expected ${ac_cv_with_mklpath}/lib/intel64 or  ${ac_cv_with_mklpath}/lib/em64t ])
   fi
     
   if test x"${enable_shared}" = x"yes"; then
      LDFLAGS="-L${mkllibpath} -Wl,-rpath,${mkllibpath} ${LDFLAGS}"
   fi
   
   # cxxname is set in sxcompflags.m4
   case "$cxxname" in
      g++)
         mkl_thread=mkl_gnu_thread
         ;;
      icpc)
         mkl_thread=mkl_intel_thread
         ;;
      *)
         AC_MSG_ERROR([Intel MKL is not supported for your compiler $CXX.])
   esac
fi

# --- Determine suitable set of numlibs
if test x"$ac_cv_with_sxmath" = x"yes"; then
   if test x"$ac_cv_enable_fftw" = x"yes"; then
      SX_FFT_LIBS="-lfftw3"
      SX_IF_CHECKNUMLIBS([
         AC_CHECK_HEADER([fftw3.h],,[AC_MSG_ERROR([Cannot find fftw3.h])] )
         AC_CHECK_LIB([fftw3], [fftw_plan_dft], 
                      [],
                      [AC_MSG_ERROR([Cannot link against FFTW library])])
       ])

      if test x"$ac_cv_enable_openmp" = x"yes" ; then
         OLD_LIBS="$LIBS"
         # try to find the right FFTW threaded library
         if test x"$SX_FFTW_OMP" = x ; then 
            AC_SEARCH_LIBS([fftw_init_threads], [fftw3_omp fftw3_threads], 
                           [SX_FFTW_OMP=`echo $ac_cv_search_fftw_init_threads | sed -e's/\-lfftw3_//'`],
                           [SX_IF_CHECKNUMLIBS([AC_MSG_ERROR([
============================================================================
Cannot link against threaded FFTW library.
with FFTW 3.2 it should be fftw3_threads
with FFTW 3.3 it should be fftw3_omp
============================================================================])])
                           SX_FFTW_OMP=omp
                           ], [-lfftw3])
         else
            SX_IF_CHECKNUMLIBS([
               AC_SEARCH_LIBS([fftw_init_threads], [fftw3_${SX_FFTW_OMP}], [],
                              [AC_MSG_ERROR([
============================================================================
Cannot link against openMP FFTW library fftw3_${SX_FFTW_OMP}.
To change the suffix, set SX_FFTW_OMP (e.g. SX_FFTW_OMP=omp for fftw3_omp)
with FFTW 3.2 it should be SX_FFTW_OMP=threads
with FFTW 3.3 it should be SX_FFTW_OMP=omp
============================================================================])
                              ])
            ])
         fi
         SX_FFT_LIBS="-lfftw3_${SX_FFTW_OMP} -lfftw3"
         LIBS="$OLD_LIBS"
      fi
   fi
fi

SX_NETCDF_LIBS="-lnetcdf"
SX_IF_CHECKNUMLIBS([
   AC_CHECK_LIB([netcdf], [nc_get_vars], [],
                [AC_MSG_ERROR([Cannot link against netCDF library])])
])

if test x"$ac_cv_enable_acml"     = x"yes" \
     -o x"$ac_cv_enable_acmlfft"  = x"yes"; then

   # does ACML provide fast math support?
   AC_CHECK_LIB([acml_mv], [fastpow], [SX_BLAS_LIBS="-lacml_mv $SX_BLAS_LIBS"])

   case "${host}" in
      *-mingw32*)
         acmllib="acml"
         test x"${enable_shared}" = x"yes" && acmllib="acml_dll"
         ;;
      *)
         echo "ACML-normal: ${enable_shared}"
         acmllib="acml"
         test x"${ac_cv_enable_openmp}" = x"yes" && acmllib="acml_mp"
         ;;
   esac
   SX_BLAS_LIBS="-l${acmllib} $SX_BLAS_LIBS"
   SX_IF_CHECKNUMLIBS([
      AC_CHECK_LIB([$acmllib],   [zfft1mx],       
                   [SX_BLAS_LIBS="-l${acmllib} $SX_BLAS_LIBS"],
                   [AC_MSG_ERROR([Cannot link against ACML library $acmllib, $SX_BLAS_LIBS])])
   ])
fi

# --- check if ATLAS provides BLAS/LAPACK routines
if test x"$ac_cv_enable_atlas" = x"yes"; then
   case "${host}" in
      *-darwin*)
         AC_MSG_CHECKING([for ATLAS support])
         SX_BLAS_LIBS="-framework Accelerate $SX_BLAS_LIBS"
         AC_DEFINE_UNQUOTED(USE_ACCELERATE_FRAMEWORK, "1", [Darwin ATLAS])
         AC_MSG_RESULT([Apple Accelerate Framework])
         ;;
      *)
         SX_IF_CHECKNUMLIBS([
         AC_CHECK_HEADER([cblas.h],,[AC_MSG_ERROR([Cannot find cblas.h])] )
         AC_CHECK_HEADER([f2c.h],,[AC_MSG_ERROR([Cannot find f2c.h])] )
         AC_CHECK_HEADER([clapack.h],,[AC_MSG_ERROR([Cannot find clapack.h])],
                         [#include <f2c.h>])
 
         AC_CHECK_LIB([f2c], [dtime_], [],
                      [AC_MSG_ERROR([Cannot link against F2C library])]
         )

         AC_CHECK_LIB([atlas], [ATL_dnrm2], [],
                      [AC_MSG_ERROR([Cannot link against ATLAS library])],
                      [-lf2c]
         )
         ])
         SX_BLAS_LIBS="-latlas -lf2c $SX_BLAS_LIBS"

         # by setting pt, use threaded ATLAS (which is not open MP!)
         # test x"${ac_cv_enable_openmp}" = x"yes" && pt="pt"
         SX_BLAS_LIBS="-l${pt}lapack -l${pt}cblas -l${pt}f77blas $SX_BLAS_LIBS"

         SX_IF_CHECKNUMLIBS([
         AC_CHECK_LIB([${pt}f77blas], [dnrm2_], [],
                      [AC_MSG_ERROR([Cannot link against ${pt}F77BLAS wrapper])],
                      [-latlas -lf2c]
         )
         AC_CHECK_LIB([${pt}cblas], [cblas_dnrm2], [],
                      [AC_MSG_ERROR([Cannot link against ${pt}CBLAS library])],
                      [-latlas -lf2c]
         )
         AC_SEARCH_LIBS([zhpev_], [${pt}lapack lapack], [],
                        [AC_MSG_ERROR([Cannot link against LAPACK library])],
                        [-lf77blas -latlas -lf2c]
         )
         LIBS="$SX_BLAS_LIBS $LIBS"
         AC_CHECK_FUNC([dsptri_],        
                      [],
                      [AC_MSG_ERROR([ATLAS does not fully support LAPACK])]
         )
         ])
         ;;
   esac
fi



# --- check if INTEL MKL provides BLAS/LAPACK routines
if test x"$ac_cv_enable_mkl" = x"yes"; then
   if test x"$ac_cv_enable_openmp" = x"yes" ; then
      AC_MSG_NOTICE([Using ${mkl_thread} for MKL threading])
   else
      mkl_thread=mkl_sequential
   fi
   
   if test x"${enable_shared}" = x"yes"; then
      SX_MKL_LIBS="-lmkl_intel_lp64 -lmkl_core -l${mkl_thread} -ldl"
   else
      mkdir -p mkllibs
      rm -f mkllibs/*.a 2>/dev/null
      MKLCORE1=`pwd`/mkllibs/libmkl_core_1.a
      MKLCORE2=`pwd`/mkllibs/libmkl_core_2.a
      MKLTHRD1=`pwd`/mkllibs/lib${mkl_thread}_1.a
      ln -s ${mkllibpath}/libmkl_core.a ${MKLCORE1}
      ln -s ${mkllibpath}/libmkl_core.a ${MKLCORE2}
      ln -s ${mkllibpath}/lib${mkl_thread}.a ${MKLTHRD1}
      #SX_MKL_LIBS="-L`pwd`/mkllibs -lmkl_intel_lp64 -lmkl_core -l${mkl_thread} -lmkl_core_1 -l${mkl_thread}_1 -lmkl_core_2"
      #
      # khr: appears to be broken and is not trivially fixed using Intel suggestions
      #
# GCC/ICC sequential: -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group
# GCC OMP: -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group
# ICC OMP: -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group
#      SX_MKL_LIBS="-Wl,--start-group ${mkllibpath}/libmkl_intel_lp64.a ${mkllibpath}/libmkl_core.a ${mkllibpath}/lib${mkl_thread}.a -Wl,--end-group"    
      SX_MKL_LIBS="-L${mkllibpath} -L`pwd`/mkllibs -lmkl_intel_lp64 -lmkl_core -l${mkl_thread} -lmkl_core_1 -l${mkl_thread}_1 -lmkl_core_2"
   fi
   
   SX_MKL_LIBS="$SX_MKL_LIBS -ldl -lpthread -lm"

   SX_BLAS_LIBS="$SX_MKL_LIBS $SX_BLAS_LIBS"
   LIBS="$SX_BLAS_LIBS $LIBS"

   AC_CHECK_FUNC(dnrm2,,AC_MSG_ERROR([Cannot link the MKL]))
fi


# --- check for GotoBLAS/OpenBLAS
if test x"$ac_cv_enable_goto" = x"yes"; then

   AC_CHECK_HEADER([cblas.h],,[AC_MSG_ERROR([Cannot find cblas.h])] )
   AC_CHECK_HEADER([f2c.h],,[AC_MSG_ERROR([Cannot find f2c.h])] )
   AC_CHECK_HEADER([clapack.h],,[AC_MSG_ERROR([Cannot find clapack.h])],
                   [#include <f2c.h>])

   # If you're using OpenBLAS ("libopenblas.a/so") create a symbolic link named
   # "libgoto2.a/so" or change the name of the library below:
   #
   SX_GOTO_LIBS="-lgoto2 -lgfortran -lpthread -m64"
   SX_BLAS_LIBS="$SX_GOTO_LIBS $SX_BLAS_LIBS"
   
   if test  "$ac_test_LDFLAGS" != "set"; then
      # cxxname is set in sxcompflags.m4
      case "$cxxname" in
         g++)
            LDFLAGS="-Wl,-rpath=${ac_cv_with_numlibs}/lib ${LDFLAGS}"
            ;;
         icpc)
            LDFLAGS="-Xlinker -rpath=${ac_cv_with_numlibs}/lib ${LDFLAGS}"
            ;;
         *)
            LDFLAGS=""
            AC_MSG_ERROR([Intel MKL is not supported for your compiler.])
      esac
   fi
      
fi


# --- set up the MPI compiler ---
# (the initialization of CXXLD below is necessary)
CXXLD="${CXX}"
if test x"$ac_cv_enable_mpi" = x"yes"; then
   if test x"$cxxname" = x"icpc"; then
      AC_CHECK_TOOLS(MPICXX, [mpiicpc mpic++ mpicxx])
      AC_CHECK_TOOLS(MPICC, [mpiicc mpicc])
   else
      AC_CHECK_TOOLS(MPICXX, [mpigxx mpic++ mpicxx])
      AC_CHECK_TOOLS(MPICC, [mpigcc mpicc])
   fi
   MPICXX=${ac_cv_prog_ac_ct_MPICXX}
   MPICC=${ac_cv_prog_ac_ct_MPICC}
   if test x"${MPICXX}" = x; then
      AC_MSG_ERROR([MPI C++ compiler wrapper not found.])
   fi
   AC_MSG_CHECKING([whether ${MPICXX} wraps ${CXX}])
   ${CXX}    --version > config.out.cxx     2>&1
   ${MPICXX} --version > config.out.mpicxx  2>&1
   if test -f config.out.cxx; then
      if test -f config.out.mpicxx; then
         cmp config.out.cxx config.out.mpicxx > /dev/null
         if test $? -ne 0; then
            AC_MSG_ERROR([${MPICXX} is not a wrapper for ${CXX}.])
         fi
      else
         AC_MSG_ERROR([Could not execute ${MPICXX} -v])
      fi
   else
      AC_MSG_ERROR([Could not execute ${CXX} -v])
   fi
   rm -f config.out.cxx config.out.mpicxx
   AC_MSG_RESULT([yes])
   CXXLD='${MPICXX}'
   if test \( x"$MPICXX" = x"mpigxx" \) -o \( x"$MPICXX" = x"mpiicpc" \); then
      if test x"${ac_cv_enable_openmp}" = x"yes"; then
         # Intel MPI needs to see -fopenmp/-openmp to select proper MPI library
         for sxcxxldflag in "$CXX_OPENMP" ; do
            CXXLD="$CXXLD -Xcompiler $sxcxxldflag"
         done
      fi
   fi
   AC_SUBST(CXX, "$MPICXX")
   AC_SUBST(CXXLD)
fi


LIBS=$ORIG_LIBS

# HDF5 (if requested and working)
if test -d "$ac_cv_enable_hdf5"/lib ; then
   LDFLAGS="-L$ac_cv_enable_hdf5/lib $LDFLAGS"
fi
if test -d "$ac_cv_enable_hdf5"/include ; then
   CPPFLAGS="-I$ac_cv_enable_hdf5/include $CPPFLAGS"
fi
if test ! x"$ac_cv_enable_hdf5" = x"no"; then
   AC_DEFINE([USE_HDF5])
   AC_CHECK_LIB([hdf5], [H5Fopen], [],
                AC_MSG_ERROR([Cannot link to HDF5]), [-lz -lm])
   LIBS="-lhdf5_cpp -lhdf5 -lz -lm $LIBS"
fi

          

SX_PARAM_ATLAS="$ac_cv_enable_atlas"
SX_PARAM_MKL="$ac_cv_enable_mkl"
SX_PARAM_ACML="$ac_cv_enable_acml"
SX_PARAM_FFTW="$ac_cv_enable_fftw"
SX_PARAM_ACMLFFT="$ac_cv_enable_acmlfft"
SX_PARAM_MKLFFT="$ac_cv_enable_mklfft"

SX_SHORTCUT_LIBS=""


# --- replace libraries with libtool archives
case "${host}" in
   *-mingw32*)
      AC_CHECK_LIB([shortcut], [sxnumlibs_createShortcut],
                   [SX_SHORTCUT_LIBS="-lshortcut"],
                   [AC_MSG_ERROR([Cannot link against WinShortCut helper])])

      SX_BLAS_LIBS=`echo $SX_BLAS_LIBS | sed s#-lacml_dll#${ac_cv_with_numlibs}/lib/libacml_dll.la#`
      SX_BLAS_LIBS=`echo $SX_BLAS_LIBS | sed s#-lacml#${ac_cv_with_numlibs}/lib/libacml.la#`
      SX_NETCDF_LIBS=`echo $SX_NETCDF_LIBS | sed s#-lnetcdf#${ac_cv_with_numlibs}/lib/libnetcdf.la#`
      ;;
esac

# --- PCRE2
if test x"$ac_cv_enable_pcre2" = x"yes"; then
   AC_CHECK_HEADER([pcre2.h],,[AC_MSG_ERROR([Cannot find PCRE2])],[#define PCRE2_CODE_UNIT_WIDTH 8])
   AC_DEFINE([USE_PCRE2],[],[compile package with PCRE2 support])
   LIBS="-lpcre2-8 ${LIBS}"
elif test x"$ac_cv_enable_pcre2" = x"auto" ; then
   AC_CHECK_HEADER([pcre2.h],[
      AC_DEFINE([USE_PCRE2],[],[compile package with PCRE2 support])
      LIBS="-lpcre2-8 ${LIBS}"
      AC_MSG_NOTICE([auto-enabled pcre2])
   ],AC_MSG_NOTICE([pcre2 will not be used]),[#define PCRE2_CODE_UNIT_WIDTH 8])
fi


AC_SUBST([SX_FFT_LIBS])
AC_SUBST([SX_NETCDF_LIBS])
AC_SUBST([SX_SHORTCUT_LIBS])
AC_SUBST([SX_BLAS_LIBS])
AC_SUBST([SX_PARAM_NUMLIBS])
AC_SUBST([SX_PARAM_SXMATH])
AC_SUBST([SX_PARAM_ATLAS])
AC_SUBST([SX_PARAM_MKL])
AC_SUBST([SX_PARAM_FFTW])
AC_SUBST([SX_PARAM_ACML])
AC_SUBST([SX_PARAM_ACMLFFT])
AC_SUBST([SX_PARAM_MKLFFT])
])
