// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------
// ref1 - Comp. Phys. Comm (128), 1-45 (2000)

//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <iostream>
#include <math.h>
#include <SxBlasLib.h>
#include <SxError.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

//#include <SxRandom.h>



// --- BLAS, LAPACK
#if defined(USE_VECLIB)
//#   include <veclib.h>   // --- don't use HP header files!
//#   include <lapack.h>   //     more details in SxMLIB.h
#     include <SxMLIB.h>
#elif defined (ESSE_ESSL)
    // --- before including ESSL we have to do a messy definition
    //     otherwise it doesn't compile. i don't know a better way.
    //     see also: SxFFT.h
#   define _ESV_COMPLEX_
#   include <complex>
#   include <essl.h>
#elif defined (USE_INTEL_MKL)
// mkl uses 'long long', which is not ISO C++. Disable errors here.
//#pragma GCC diagnostic push
//#pragma GCC diagnostic warning "-Wlong-long"
    extern "C"  {
#      include <mkl.h>
    }
//#pragma GCC diagnostic pop
#elif defined (USE_ACML)
    extern "C"  {
#      include <acml.h>
    }
#elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
#  define blasint int
       extern "C"  {
#         include <f2c.h>
#         include <cblas.h>
#         include <clapack.h>
       }
#elif defined (USE_ATLAS)
#  if defined (USE_ACCELERATE_FRAMEWORK)
#     include <Accelerate/Accelerate.h>
      typedef __CLPK_integer integer;
      typedef __CLPK_logical logical;
      typedef __CLPK_real real;
      typedef __CLPK_doublereal doublereal;
      typedef __CLPK_ftnlen ftnlen;
      typedef __CLPK_complex complex;
      typedef __CLPK_doublecomplex doublecomplex;
#  else
       extern "C"  { 
#         include <f2c.h>
#         include <cblas.h>
#         include <clapack.h>
       }
#  endif /* USE_ACCELERATE_FRAMEWORK */
#else
#   error "No numeric library specified"
// make sure that we do not drown in error messages
#define SX_IGNORE_THE_REST_OF_THE_FILE
#endif

#ifndef SX_IGNORE_THE_REST_OF_THE_FILE
//------------------------------------------------------------------------------
// BLAS/LAPACK error handling   
//------------------------------------------------------------------------------
#if   defined (USE_VECLIB)
   // error handling not yet supported
#elif defined (USE_ESSL)
   // error handling not yet supported
#elif defined (USE_INTEL_MKL)
   // error handling not yet supported
#elif defined (USE_ACML)
   // error handling not yet supported
#elif defined (USE_GOTO)
   // error handling not yet supported
#else
#  include <stdarg.h>

#ifdef MACOSX
#  if ( DIST_VERSION_L >= 1070L )
      extern "C" void cblas_xerbla (int p, char *rout, char *form, ...)
#  else
      extern "C" void cblas_xerbla (int p, const char *rout, const char *form, ...)
#  endif /* DIST_VERSION_L */
#else
   extern "C" void cblas_xerbla (int p, const char *rout, const char *form, ...)
#endif /* MACOSX */
   { 
      //sxprintf ("\nA BLAS/LAPACK error has occured!\n");
      std::cout << "\nA BLAS/LAPACK error has occured!\n";
      // --- original code from ATLAS/interfaces/blas/C/src/cblas_xerbla.c
      va_list argptr;
      va_start(argptr, form);
#     ifdef GCCWIN
         //if (p)  sxprintf("Parameter %d to routine %s was incorrect\n", p, rout);
         if (p) std::cout << "Parameter " << p << " to routine "
                          << rout << " was incorrect\n";
         vprintf(form, argptr);
#     else
         if (p)
            fprintf(stderr, "Parameter %d to routine %s was incorrect\n", 
                    p, rout);
         vfprintf(stderr, form, argptr);
#     endif /* CYGWIN */
      va_end(argptr);
      // --- end of original ATLAS code
      SX_EXIT;
   }

#endif /* USE_VECLIB */


//------------------------------------------------------------------------------
// norm of vectors
//------------------------------------------------------------------------------
float norm2 (const float *vec, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      return snrm2 (&n, (float *)vec, &incx);
#  elif defined (USE_ESSL)
      return snrm2 (n, vec, incx);
#  elif defined (USE_INTEL_MKL)
      return snrm2 (&n, const_cast<float *>(vec), &incx);
#  elif defined (USE_ACML)
      return snrm2 (n, const_cast<float *>(vec), incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      return cblas_snrm2 (n, const_cast<float *>(vec), incx);
#  else 
      return cblas_snrm2 (n, vec, incx);
#  endif   
}

double norm2 (const double *vec, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      return dnrm2 (&n, (double *)vec, &incx);
#  elif defined (USE_ESSL)
      return dnrm2 (n, vec, incx);
#  elif defined (USE_INTEL_MKL)
      return dnrm2 (&n, const_cast<double *>(vec), &incx);
#  elif defined (USE_ACML)
      return dnrm2 (n, const_cast<double *>(vec), incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      return cblas_dnrm2 (n, const_cast<double*>(vec), incx);
#  else 
      return cblas_dnrm2 (n, vec, incx);
#  endif      
}


float norm2 (const SxComplex8 *vec, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      return scnrm2 (&n, (complex8_t *)vec, &incx);
#  elif defined (USE_ESSL)
      return scnrm2 (n, (const complex<float> *)vec,  incx);
#  elif defined (USE_INTEL_MKL)
      return scnrm2 (&n, (MKL_Complex8 *)const_cast<SxComplex8 *>(vec), &incx);
#  elif defined (USE_ACML)
      return scnrm2 (n, (complex *)const_cast<SxComplex8 *>(vec), incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      return cblas_scnrm2 (n, (float*)vec, incx);
#  else 
      return cblas_scnrm2 (n, vec, incx);
#  endif      
}


double norm2 (const SxComplex16 *vec, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      return dznrm2 (&n, (complex16_t *)vec, &incx);
#  elif defined (USE_ESSL)
      return dznrm2 (n, (const complex<double> *)vec, incx);
#  elif defined (USE_INTEL_MKL)
      return dznrm2 (&n, (MKL_Complex16 *)const_cast<SxComplex16*>(vec), &incx);
#  elif defined (USE_ACML)
      return dznrm2 (n, (doublecomplex *)const_cast<SxComplex16 *>(vec), 
                     incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      return cblas_dznrm2 (n, (double*)vec, incx);
#  else
      return cblas_dznrm2 (n, vec, incx);
#  endif      
}



//------------------------------------------------------------------------------
// scale vectors
//------------------------------------------------------------------------------
void scale (float *vec,  const float alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      sscal (&n, (float *)&alpha, (float *)vec, &incx);
#  elif defined (USE_ESSL)
      sscal (n, alpha, vec, incx);
#  elif defined (USE_INTEL_MKL)
      sscal (&n, const_cast<float *>(&alpha), (float *)vec, &incx);
#  elif defined (USE_ACML)
      sscal (n, alpha, vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_sscal (n, alpha, vec, incx);
#  else
      cblas_sscal (n, alpha, vec, incx);
#  endif      
}


void scale (double *vec, const double alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      dscal (&n, (double *)&alpha, (double *)vec, &incx);
#  elif defined (USE_ESSL)
      dscal (n, alpha, vec, incx);
#  elif defined (USE_INTEL_MKL)
      dscal (&n, const_cast<double *>(&alpha), (double *)vec, &incx);
#  elif defined (USE_ACML)
      dscal (n, alpha, vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_dscal (n, alpha, vec, incx);
#  else
      cblas_dscal (n, alpha, vec, incx);
#  endif      
}

void scale (SxComplex8 *vec,  const SxComplex8 &alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      cscal (&n, (complex8_t *)&alpha, (complex8_t *)vec, &incx);
#  elif defined (USE_ESSL)
      const complex<float> alphaTmp (alpha.re, alpha.im);
      cscal (n, alphaTmp, (complex<float> *)vec, incx);
#  elif defined (USE_INTEL_MKL)
      cscal (&n, (MKL_Complex8 *)const_cast<SxComplex8*>(&alpha),
             (MKL_Complex8 *)vec, &incx);
#  elif defined (USE_ACML)
      cscal (n, (complex *)const_cast<SxComplex8 *>(&alpha), 
             (complex *)vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_cscal (n, (float*)&alpha, (float*)vec, incx);
#  else
      cblas_cscal (n, &alpha, vec, incx);
#  endif      
}


void scale (SxComplex16 *vec, const SxComplex16 &alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      zscal (&n, (complex16_t *)&alpha, (complex16_t *)vec, &incx);
#  elif defined (USE_ESSL)
      const complex<double> alphaTmp (alpha.re, alpha.im);
      zscal (n, alphaTmp, (complex<double> *)vec, incx);
#  elif defined (USE_INTEL_MKL)
      zscal (&n, (MKL_Complex16 *)const_cast<SxComplex16*>(&alpha),
             (MKL_Complex16 *)vec, &incx);
#  elif defined (USE_ACML)
      zscal (n, (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
                (doublecomplex *)vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_zscal (n, (double*)&alpha, (double*)vec, incx);
#  else
      cblas_zscal (n, &alpha, vec, incx);
#  endif      
}


//------------------------------------------------------------------------------
// Y += a*X 
//------------------------------------------------------------------------------
void axpy (float *yOut, const float &alpha, const float *xIn, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      saxpy (&n, (float *)&alpha, 
             (float *)xIn, &incx, (float *)yOut, &incx);
#  elif defined (USE_ESSL)
      saxpy (n, a, 
            (float *)xIn, incx, (float *)yOut, incx);
#  elif defined (USE_INTEL_MKL)
      saxpy (&n, const_cast<float *>(&alpha),
            const_cast<float *>(xIn), &incx, yOut, &incx);
#  elif defined (USE_ACML)
      saxpy (n, alpha, (float *)const_cast<float *>(xIn), 
             incx, (float *)yOut, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_saxpy (n, alpha, const_cast<float*>(xIn), incx, yOut, incx);
#  else
      cblas_saxpy (n, alpha, xIn, incx, yOut, incx);
#  endif      
}

void axpy (double *yOut, const double &alpha, const double *xIn, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      daxpy (&n, (double *)&alpha, 
             (double *)xIn, &incx, (double *)yOut, &incx);
#  elif defined (USE_ESSL)
      daxpy (n, a, 
            (double *)xIn, incx, (double *)yOut, incx);
#  elif defined (USE_INTEL_MKL)
      daxpy (&n, const_cast<double *>(&alpha),
            const_cast<double *>(xIn), &incx, (double *)yOut, &incx);
#  elif defined (USE_ACML)
      daxpy (n, alpha, (double *)const_cast<double *>(xIn), 
             incx, (double *)yOut, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_daxpy (n, alpha, const_cast<double*>(xIn), incx, yOut, incx);
#  else
      cblas_daxpy (n, alpha, xIn, incx, yOut, incx);
#  endif      
}

void axpy (SxComplex8 *yOut, 
           const SxComplex8 &alpha, const SxComplex8 *xIn, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      caxpy (&n, (complex8_t *)&alpha, 
             (complex8_t *)xIn, &incx, (complex8_t *)yOut, &incx);
#  elif defined (USE_ESSL)
      const complex<double> alphaTmp (alpha.re, alpha.im);
      caxpy (n, alphaTmp, 
            (complex<float> *)xIn, incx, (complex<float> *)yOut, incx);
#  elif defined (USE_INTEL_MKL)
      caxpy (&n, (MKL_Complex8 *)const_cast<SxComplex8 *>(&alpha),
            (MKL_Complex8 *)const_cast<SxComplex8 *>(xIn), &incx,
            (MKL_Complex8 *)yOut, &incx);
#  elif defined (USE_ACML)
      caxpy (n, (complex *)const_cast<SxComplex8 *>(&alpha), 
                (complex *)const_cast<SxComplex8 *>(xIn), incx, 
                (complex *)yOut, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_caxpy (n, (float*)&alpha, (float*)xIn, incx, (float*)yOut, incx);
#  else
      cblas_caxpy (n, &alpha, xIn, incx, yOut, incx);
#  endif      
}

void axpy (SxComplex16 *yOut, 
           const SxComplex16 &alpha, const SxComplex16 *xIn, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      zaxpy (&n, (complex16_t *)&alpha, 
             (complex16_t *)xIn, &incx, (complex16_t *)yOut, &incx);
#  elif defined (USE_ESSL)
      const complex<double> alphaTmp (alpha.re, alpha.im);
      zaxpy (n, alphaTmp, 
            (complex<double> *)xIn, incx, (complex<double> *)yOut, incx);
#  elif defined (USE_INTEL_MKL)
      zaxpy (&n, (MKL_Complex16 *)const_cast<SxComplex16 *>(&alpha),
            (MKL_Complex16 *)const_cast<SxComplex16 *>(xIn), &incx,
            (MKL_Complex16 *)yOut, &incx);
#  elif defined (USE_ACML)
      zaxpy (n, (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
            (doublecomplex *)const_cast<SxComplex16 *>(xIn), 
            incx, (doublecomplex *)yOut, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_zaxpy (n, (double*)&alpha, (double*)xIn, incx, (double*)yOut, incx);
#  else
      cblas_zaxpy (n, &alpha, xIn, incx, yOut, incx);
#  endif      
}


//------------------------------------------------------------------------------
// scalar product
//------------------------------------------------------------------------------

float scalarProduct (const float *aVec, const float *bVec, int n)
{
   const float *aPtr = aVec;
   const float *bPtr = bVec;
   float res = 0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += *aPtr * *bPtr;
   return res;
}



double scalarProduct (const double *aVec, const double *bVec, int n)
{
   const double *aPtr = aVec;
   const double *bPtr = bVec;
   double res = 0;
   for (int i=0; i < n; i++, aPtr++, bPtr++)  
      res += *aPtr * *bPtr;
   return res;
}


SxComplex8 scalarProduct (const SxComplex8 *aVec, const SxComplex8 *bVec, int n)
{
   const SxComplex8 *aPtr = aVec;
   const SxComplex8 *bPtr = bVec;
   SxComplex8 res = (SxComplex8)0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += aPtr->conj() * *bPtr;
   return res;
}


SxComplex16 scalarProduct (const SxComplex16 *aVec, const SxComplex16 *bVec, int n)
{
   const SxComplex16 *aPtr = aVec;
   const SxComplex16 *bPtr = bVec;
   SxComplex16 res = (SxComplex8)0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += aPtr->conj() * *bPtr;
   return res;
}


//------------------------------------------------------------------------------
// general matrix-matrix multiplication
//------------------------------------------------------------------------------
void matmult (float *resMat, const float *aMat, const float *bMat, 
              int aMatRows, int aMatCols, int bMatCols)
{
   float alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      sgemm (&noTrans, &noTrans, 
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, (float *)aMat, &aMatRows,
             (float *)bMat, &aMatCols, &beta, (float *)resMat, &aMatRows, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      sgemm (&noTrans, &noTrans, 
             aMatRows, bMatCols, aMatCols,
             alpha, aMat, aMatRows,
             bMat, aMatCols, beta, resMat, aMatRows);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      sgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             &alpha, const_cast<float *>(aMat), &aMatRows,
             const_cast<float *>(bMat), &aMatCols, &beta, resMat, &aMatRows);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      sgemm (noTrans, noTrans, aMatRows, bMatCols, aMatCols,
             alpha, const_cast<float *>(aMat), aMatRows,
             const_cast<float *>(bMat), aMatCols, beta, resMat, aMatRows);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                         aMatRows, bMatCols, aMatCols,
                         alpha, const_cast<float*>(aMat), aMatRows,
                         const_cast<float*>(bMat), aMatCols, beta, resMat,
                         aMatRows);
#  else  
      cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols, 
                   alpha, aMat, aMatRows, 
                   bMat, aMatCols, beta, resMat, aMatRows);
#  endif
}


void matmult (double *resMat, const double *aMat, const double *bMat, 
              int aMatRows, int aMatCols, int bMatCols)
{
   double alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      dgemm (&noTrans, &noTrans, 
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, (double *)aMat, &aMatRows,
             (double *)bMat, &aMatCols, &beta, (double *)resMat, &aMatRows, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      dgemm (&noTrans, &noTrans, 
             aMatRows, bMatCols, aMatCols,
             alpha, aMat, aMatRows,
             bMat, aMatCols, beta, resMat, aMatRows);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      dgemm (&noTrans, &noTrans, 
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, const_cast<double *>(aMat), &aMatRows,
             const_cast<double *>(bMat), &aMatCols, &beta, resMat, &aMatRows);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      dgemm (noTrans, noTrans, 
             aMatRows, bMatCols, aMatCols,
             alpha, const_cast<double *>(aMat), aMatRows,
             const_cast<double *>(bMat), aMatCols, beta, resMat, aMatRows);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                         aMatRows, bMatCols, aMatCols,
                         alpha, const_cast<double*>(aMat), aMatRows,
                         const_cast<double*>(bMat), aMatCols, beta, resMat, aMatRows);
#  else  
      cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols, 
                   alpha, aMat, aMatRows, 
                   bMat, aMatCols, beta, resMat, aMatRows);
#  endif      
}


void matmult (SxComplex8 *resMat, 
              const SxComplex8 *aMat, const SxComplex8 *bMat, 
              int aMatRows, int aMatCols, int bMatCols)
{
   SxComplex8 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      cgemm (&noTrans, &noTrans, 
             &aMatRows, &bMatCols, &aMatCols,
             (complex8_t *)&alpha, (complex8_t *)aMat, &aMatRows,
             (complex8_t *)bMat,   &aMatCols, (complex8_t *)&beta, 
             (complex8_t *)resMat, &aMatRows, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      complex<float> alphaT (alpha.re, alpha.im);
      complex<float> betaT  (beta.re,  beta.im);
      cgemm (&noTrans, &noTrans, 
             aMatRows, bMatCols, aMatCols,
             alphaT, (complex<float> *)aMat, aMatRows,
             (complex<float> *)bMat,   aMatCols, betaT, 
             (complex<float> *)resMat, aMatRows);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      cgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             (MKL_Complex8 *)&alpha,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(aMat), &aMatRows,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(bMat), &aMatCols,
             (MKL_Complex8 *)&beta,
             (MKL_Complex8 *)resMat, &aMatRows);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      cgemm (noTrans, noTrans, 
             aMatRows, bMatCols, aMatCols,
             (complex *)const_cast<SxComplex8 *>(&alpha), 
             (complex *)const_cast<SxComplex8 *>(aMat), aMatRows,
             (complex *)const_cast<SxComplex8 *>(bMat), aMatCols, 
             (complex *)const_cast<SxComplex8 *>(&beta), 
             (complex *)resMat, aMatRows);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   (float*)&alpha, (float*)aMat, aMatRows,
                   (float*)bMat, aMatCols, (float*)&beta,
                   (float*)resMat, aMatRows);
#  else  
      cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols, 
                   &alpha, aMat, aMatRows, 
                   bMat, aMatCols, &beta, resMat, aMatRows);
#  endif      
}


void matmult (SxComplex16 *resMat, 
              const SxComplex16 *aMat, const SxComplex16 *bMat, 
              int aMatRows, int aMatCols, int bMatCols)
{
   SxComplex16 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      zgemm (&noTrans, &noTrans, 
             &aMatRows, &bMatCols, &aMatCols,
             (complex16_t *)&alpha, (complex16_t *)aMat, &aMatRows,
             (complex16_t *)bMat,   &aMatCols, (complex16_t *)&beta, 
             (complex16_t *)resMat, &aMatRows, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      complex<double> alphaT (alpha.re, alpha.im);
      complex<double> betaT  (beta.re,  beta.im);
      zgemm (&noTrans, &noTrans, 
             aMatRows, bMatCols, aMatCols,
             alphaT, (complex<double> *)aMat, aMatRows,
             (complex<double> *)bMat,   aMatCols, betaT, 
             (complex<double> *)resMat, aMatRows);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      zgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             (MKL_Complex16 *)&alpha,
             (MKL_Complex16 *)const_cast<SxComplex16 *>(aMat), &aMatRows,
             (MKL_Complex16 *)const_cast<SxComplex16 *>(bMat), &aMatCols,
             (MKL_Complex16 *)&beta, (MKL_Complex16 *)resMat, &aMatRows);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      zgemm (noTrans, noTrans, 
             aMatRows, bMatCols, aMatCols,
             (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
             (doublecomplex *)const_cast<SxComplex16 *>(aMat), aMatRows,
             (doublecomplex *)const_cast<SxComplex16 *>(bMat), aMatCols, 
             (doublecomplex *)const_cast<SxComplex16 *>(&beta), 
             (doublecomplex *)resMat, aMatRows);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasNoTrans, aMatRows, aMatCols,
                      (double*)&alpha, (double*)aMat, aMatRows,
                      (double*)bMat, 1, (double*)&beta, (double*)resMat, 1);
      } else if (aMatRows == 1)  {
         cblas_zgemv (CblasColMajor, CblasTrans, aMatCols, bMatCols,
                      (double*)&alpha, (double*)bMat, aMatCols,
                      (double*)aMat, 1, (double*)&beta, (double*)resMat, 1);
      } else {
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      aMatRows, bMatCols, aMatCols,
                      (double*)&alpha, (double*)aMat, aMatRows,
                      (double*)bMat, aMatCols, (double*)&beta, (double*)resMat,
                      aMatRows);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasNoTrans, aMatRows, aMatCols,
                      &alpha, aMat, aMatRows,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatRows == 1)  {
         cblas_zgemv (CblasColMajor, CblasTrans, aMatCols, bMatCols,
                      &alpha, bMat, aMatCols,
                      aMat, 1, &beta, resMat, 1);
      } else {
#if defined(USE_OPENMP) && defined(USE_ATLAS)
         if (aMatRows > 1024)  {
#           pragma omp parallel
            {
               int nThreads = omp_get_num_threads ();
               int nb = aMatRows / nThreads;
               // distribute remaining elements over all threads
               if (nb * nThreads < aMatRows) nb++;
               int offset = omp_get_thread_num () * nb;
               // reduce nb for last thread
               if (offset + nb > aMatRows) nb = aMatRows - offset;
               cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nb, bMatCols, aMatCols, 
                            &alpha, aMat + offset, aMatRows, 
                            bMat, aMatCols, &beta, resMat + offset, aMatRows);
            }
         } else // no openMP parallelism ...
#endif
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      aMatRows, bMatCols, aMatCols, 
                      &alpha, aMat, aMatRows, 
                      bMat, aMatCols, &beta, resMat, aMatRows);
      }
#  endif      
}

//------------------------------------------------------------------------------
// overlap matrices
//------------------------------------------------------------------------------
void matovlp (float *resMat, 
              const float *aMat, const float *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   float alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, (float *)aMat, &aMatRows,
             (float *)bMat,   &bMatRows, &beta, 
             (float *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, aMat, aMatRows,
             bMat, bMatRows, beta, resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
//      SX_EXIT; // not tested
      // change by khr
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, const_cast<float *>(aMat), &aMatRows,
             const_cast<float *>(bMat), &bMatRows, &beta, resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, const_cast<float *>(aMat), aMatRows,
             const_cast<float *>(bMat), bMatRows, 
             beta, resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
            cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                         alpha, const_cast<float*>(aMat), aMatRows,
                         const_cast<float*>(bMat), 1, beta, resMat, 1);
         } else if (aMatCols == 1)  {
            cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                         alpha, const_cast<float*>(bMat), bMatRows,
                         const_cast<float*>(aMat), 1, beta, resMat, 1);
         } else {
            cblas_sgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize,
                         alpha, const_cast<float *>(aMat), aMatRows,
                         const_cast<float*>(bMat), bMatRows, beta, resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, aMat, aMatRows,
                      bMat, 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, bMat, bMatRows,
                      aMat, 1, beta, resMat, 1);
      } else {
         cblas_sgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      alpha, aMat, aMatRows, 
                      bMat, bMatRows, beta, resMat, aMatCols);
      }
#  endif      
}

void matovlp (double *resMat, 
              const double *aMat, const double *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   double alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, (double *)aMat, &aMatRows,
             (double *)bMat,   &bMatRows, &beta, 
             (double *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, aMat, aMatRows,
             bMat, bMatRows, beta, resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
//      SX_EXIT; // not tested
      // change by khr
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, const_cast<double *>(aMat), &aMatRows,
             const_cast<double *>(bMat), &bMatRows, &beta, resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, const_cast<double *>(aMat), aMatRows,
             const_cast<double *>(bMat), bMatRows, 
             beta, resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, const_cast<double*>(aMat), aMatRows,
                      const_cast<double*>(bMat), 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, const_cast<double*>(bMat), bMatRows,
                      const_cast<double*>(aMat), 1, beta, resMat, 1);
      } else {
         cblas_dgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize,
                      alpha, const_cast<double*>(aMat), aMatRows,
                      const_cast<double*>(bMat), bMatRows, beta, resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, aMat, aMatRows,
                      bMat, 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, bMat, bMatRows,
                      aMat, 1, beta, resMat, 1);
      } else {
         cblas_dgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      alpha, aMat, aMatRows, 
                      bMat, bMatRows, beta, resMat, aMatCols);
      }
#  endif      
}

void matovlp (SxComplex8 *resMat, 
              const SxComplex8 *aMat, const SxComplex8 *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   SxComplex8 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      cgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             (complex8_t *)&alpha, (complex8_t *)aMat, &aMatRows,
             (complex8_t *)bMat,   &bMatRows, (complex8_t *)&beta, 
             (complex8_t *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      complex<float> alphaT (alpha.re, alpha.im);
      complex<float> betaT  (beta.re,  beta.im);
      cgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alphaT, (complex<float> *)aMat, aMatRows,
             (complex<float> *)bMat,   bMatRows, betaT, 
             (complex<float> *)resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N', conjTrans = 'C';
      cgemm (&conjTrans, &noTrans, &aMatCols, &bMatCols, &sumSize,
             (MKL_Complex8 *)&alpha,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(aMat), &aMatRows,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(bMat), &bMatRows,
             (MKL_Complex8 *)&beta, (MKL_Complex8 *)resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      cgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             (complex *)const_cast<SxComplex8 *>(&alpha), 
             (complex *)const_cast<SxComplex8 *>(aMat), aMatRows,
             (complex *)const_cast<SxComplex8 *>(bMat), bMatRows, 
             (complex *)const_cast<SxComplex8 *>(&beta), 
             (complex *)resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      (float*)&alpha, (float*)aMat, aMatRows,
                      (float*)bMat, 1, (float*)&beta, (float*)resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      (float*)&alpha, (float*)bMat, bMatRows,
                      (float*)aMat, 1, (float*)&beta, (float*)resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else {
         cblas_cgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize,
                      (float*)&alpha, (float*)aMat, aMatRows,
                      (float*)bMat, bMatRows, (float*)&beta, (float*)resMat,
                      aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      &alpha, aMat, aMatRows,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      &alpha, bMat, bMatRows,
                      aMat, 1, &beta, resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else {
         cblas_cgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      &alpha, aMat, aMatRows, 
                      bMat, bMatRows, &beta, resMat, aMatCols);
      }
#  endif      
}

/*
#include <SxTimer.h>
enum BlasLibTimer { Matovlp };
REGISTER_TIMERS (BlasLibTimer)
{
   regTimer (Matovlp, "matovlp");
}
*/

void matovlp (SxComplex16 *resMat, 
              const SxComplex16 *aMat, const SxComplex16 *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   SxComplex16 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      zgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             (complex16_t *)&alpha, (complex16_t *)aMat, &aMatRows,
             (complex16_t *)bMat,   &bMatRows, (complex16_t *)&beta, 
             (complex16_t *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      complex<double> alphaT (alpha.re, alpha.im);
      complex<double> betaT  (beta.re,  beta.im);
      zgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alphaT, (complex<double> *)aMat, aMatRows,
             (complex<double> *)bMat,   bMatRows, betaT, 
             (complex<double> *)resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N', conjTrans = 'C';
      if (aMat == bMat && aMatRows == bMatRows && aMatCols == bMatCols)  {
         // special case: A.overlap (A) => use Hermitean routine ...
         char uplo = 'U';
         zherk (&uplo, &conjTrans, &aMatCols, &sumSize,
                &alpha.re,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(aMat), &aMatRows,
                &beta.re, (MKL_Complex16 *)resMat, &aMatCols);
         // ... and copy upper half to lower half
         for (int i = 0; i < aMatCols; i++)
            for (int j = 0; j < i; j++)
               resMat[i + aMatCols * j] = resMat[j + aMatCols * i].conj ();
      } else {
         zgemm (&conjTrans, &noTrans, &aMatCols, &bMatCols, &sumSize,
                (MKL_Complex16 *)&alpha,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(aMat), &aMatRows,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(bMat), &bMatRows,
                (MKL_Complex16 *)&beta, (MKL_Complex16 *)resMat, &aMatCols);
      }
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      zgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
             (doublecomplex *)const_cast<SxComplex16 *>(aMat), aMatRows,
             (doublecomplex *)const_cast<SxComplex16 *>(bMat), bMatRows, 
             (doublecomplex *)const_cast<SxComplex16 *>(&beta), 
             (doublecomplex *)resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
            cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                         (double*)&alpha, (double*)aMat, aMatRows,
                         (double*)bMat, 1, (double*)&beta, (double*)resMat, 1);
         } else if (aMatCols == 1)  {
            cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                         (double*)&alpha, (double*)bMat, bMatRows,
                         (double*)aMat, 1, (double*)&beta, (double*)resMat, 1);
            // now resMat contains (B.adjoint () ^ A)
            // perform conjugate (note: res is a vector)
            for (ssize_t i = 0; i < bMatCols; ++i)
               resMat[i].im = -resMat[i].im;
         } else {
            cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize,
                         (double*)&alpha, (double*)aMat, aMatRows,
                         (double*)bMat, bMatRows, (double*)&beta,
                         (double*)resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      &alpha, aMat, aMatRows,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      &alpha, bMat, bMatRows,
                      aMat, 1, &beta, resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else {
#ifndef USE_OPENMP
         cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      &alpha, aMat, aMatRows, 
                      bMat, bMatRows, &beta, resMat, aMatCols);
#else
         //CLOCK (Matovlp);
         int nb = 128;
         int np = sumSize / nb;

         // --- compute contribution from rest elements
         if (sumSize > nb * np)  {
            cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize - nb * np, 
                         &alpha, aMat + nb * np, aMatRows, 
                         bMat + nb * np, bMatRows, &beta, resMat, aMatCols);
         } else {
            for (int i = 0; i < aMatCols * bMatCols; ++i)
               resMat[i].re = resMat[i].im = 0.;

         }
#        pragma omp parallel
         {
            SxComplex16 *part = NULL;
            // --- compute partial results (in parallel)
#           pragma omp for
            for (int ip = 0; ip < np; ip++) {
               if (!part)  {
                  part = new SxComplex16[aMatCols * bMatCols];
                  cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                               aMatCols, bMatCols, nb, 
                               &alpha, aMat + ip * nb, aMatRows, 
                               bMat + ip * nb, bMatRows, &beta, part, aMatCols);
               } else {
                  cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                               aMatCols, bMatCols, nb, 
                               &alpha, aMat + ip * nb, aMatRows, 
                               bMat + ip * nb, bMatRows, &alpha, part, aMatCols);
               }
            }
            // --- sum partial results
            if (part)  {
#              pragma omp critical
               cblas_zaxpy (aMatCols * bMatCols, &alpha, part, 1, resMat, 1);
               delete [] part;
            }
         }
#     endif
      }
#  endif      
}

//------------------------------------------------------------------------------
// Matrix decompositions
//------------------------------------------------------------------------------
void cholesky (float *resMat, enum UPLO uplo, float * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      spotrf (&uploChar, &n, (float *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      spotrf (&uploChar, n, resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      spotrf (&uploChar, &n, resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      spotrf (uploChar, n, resMat, n, &err);
#  else 
      integer rank  = (integer)n, err = 0;
      spotrf_ (&uploChar, &rank, (real *)resMat, &rank, &err);
      if ( err )  { // TODO: throw exception
         std::cout << "cholesky err=" << err << std::endl;
         SX_EXIT;
      }
#   endif      
}


void cholesky (double *resMat, enum UPLO uplo, double * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      dpotrf (&uploChar, &n, (double *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      dpotrf (&uploChar, n, resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      dpotrf (&uploChar, &n, resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      dpotrf (uploChar, n, resMat, n, &err);
#  else 
      integer rank  = (integer)n, err = 0;
      dpotrf_ (&uploChar, &rank, (doublereal *)resMat, &rank, &err);
      if ( err )  { // TODO: throw exception
         std::cout << "cholesky err=" << err << std::endl;
         SX_EXIT;
      }
#  endif      
}


void cholesky (SxComplex8 *resMat, enum UPLO uplo, 
               SxComplex8 * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      cpotrf (&uploChar, &n, (complex8_t *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      cpotrf (&uploChar, n, (complex<float> *)resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      cpotrf (&uploChar, &n, (MKL_Complex8 *)resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      cpotrf (uploChar, n, (complex *)resMat, n, &err);
#  else 
      integer rank  = (integer)n, err = 0;
      cpotrf_ (&uploChar, &rank, (complex *)resMat, &rank, &err);
      if ( err )  { // TODO: throw exception
         std::cout << "cholesky err=" << err << std::endl;
         SX_EXIT;
      }
#  endif      
}


void cholesky (SxComplex16 *resMat, enum UPLO uplo, 
               SxComplex16 * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      zpotrf (&uploChar, &n, (complex16_t *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      zpotrf (&uploChar, n, (complex<double> *)resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      zpotrf(&uploChar, &n, (MKL_Complex16 *)resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      zpotrf(uploChar, n, (doublecomplex *)resMat, n, &err);
#  else 
      integer rank  = (integer)n, err = 0;
      zpotrf_ (&uploChar, &rank, (doublecomplex *)resMat, &rank, &err);
      if ( err )  { // TODO: throw exception
         std::cout << "cholesky err=" << err << std::endl;
         SX_EXIT;
      }
#  endif   
}

void singularValueDecomp (float *mat, int nRows, int nCols,
                                 float *vals,
                                 float *left,
                                 float *right, // V^H
                                 bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
   int ldvt = zeroSpace ? nCols : minMN;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   // workspace query:
   int *iwork = new int[8 * minMN];
   int err = 0;
   int lwork = -1;
   float optwork;
   sgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (int)optwork;
   float* work = new float[lwork];
   // --- actual compute
   sgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  else
   integer *iwork = new integer[8 * minMN];
   integer lwork = -1;
   integer err = 0;
   // workspace query:
   float optwork;
   integer nr = nRows, nc = nCols, ldvt_ = ldvt;
   sgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork;
   float* work = new float[lwork];
   // --- actual compute
   sgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   
   delete[] iwork;
}

void singularValueDecomp (double *mat, int nRows, int nCols,
                                 double *vals,
                                 double *left,
                                 double *right, // V^H
                                 bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
   int ldvt = zeroSpace ? nCols : minMN;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   // workspace query:
   int err = 0;
   int *iwork = new int[8 * minMN];
   int lwork = -1;
   double optwork;
   dgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork));
   double* work = new double[lwork];
   // --- actual compute
   dgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  else
   integer err = 0;
   integer lwork = -1;
   integer *iwork = new integer[8 * minMN];
   // workspace query:
   double optwork;
   integer nr = nRows, nc = nCols, ldvt_ = ldvt;
   dgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork;
   double* work = new double[lwork];
   // --- actual compute
   dgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   
   delete[] iwork;
}


void singularValueDecomp (SxComplex8 *mat, int nRows, int nCols,
                                 float *vals,
                                 SxComplex8 *left,
                                 SxComplex8 *right, // V^H
                                 bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   float *rwork = new float[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   MKL_Complex8 optwork;
   cgesdd(&jobz, &nRows, &nCols, (MKL_Complex8*)mat, &nRows,
          vals, (MKL_Complex8*)left, &nRows, (MKL_Complex8*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.real));
   MKL_Complex8* work = new MKL_Complex8[lwork];
   // --- actual compute
   cgesdd(&jobz, &nRows, &nCols, (MKL_Complex8*)mat, &nRows,
          vals, (MKL_Complex8*)left, &nRows, (MKL_Complex8*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  else
   integer ldvt = zeroSpace ? nCols : minMN;
   integer lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   integer err = 0;
   integer lwork = -1;
   float *rwork = new float[lrwork];
   integer *iwork = new integer[8 * minMN];
   // workspace query:
   integer nr = nRows, nc = nCols;
   complex optwork;
   SX_EXIT;
   /* missing prototype
   cgesdd_(&jobz, &nr, &nc, (complex *)mat, &nr,
          vals, (complex*)left, &nr, (complex*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
          */
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork.r;
   complex* work = new complex[lwork];
   // --- actual compute
   SX_EXIT;
   /* missing prototype
   cgesdd_(&jobz, &nr, &nc, (complex *)mat, &nr,
          vals, (complex*)left, &nr, (complex*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
          */
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   

}

void singularValueDecomp (SxComplex16 *mat, int nRows, int nCols,
                                 double *vals,
                                 SxComplex16 *left,
                                 SxComplex16 *right, // V^H
                                 bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   double *rwork = new double[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   MKL_Complex16 optwork;
   zgesdd(&jobz, &nRows, &nCols, (MKL_Complex16*)mat, &nRows,
          vals, (MKL_Complex16*)left, &nRows, (MKL_Complex16*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.real));
   MKL_Complex16* work = new MKL_Complex16[lwork];
   // --- actual compute
   zgesdd(&jobz, &nRows, &nCols, (MKL_Complex16*)mat, &nRows,
          vals, (MKL_Complex16*)left, &nRows, (MKL_Complex16*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  else
   integer ldvt = zeroSpace ? nCols : minMN;
   integer lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   integer err = 0;
   integer lwork = -1;
   double *rwork = new double[lrwork];
   integer *iwork = new integer[8 * minMN];
   // workspace query:
   doublecomplex optwork;
   integer nr = nRows, nc = nCols;
   SX_EXIT;
   /* missing prototype:
   zgesdd_(&jobz, &nr, &nc, (doublecomplex *)mat, &nr,
          vals, (doublecomplex*)left, &nr, (doublecomplex*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
          */
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork.r;
   doublecomplex* work = new doublecomplex[lwork];
   // --- actual compute
   /* missing prototype:
   zgesdd_(&jobz, &nr, &nc, (doublecomplex *)mat, &nr,
          vals, (doublecomplex*)left, &nr, (doublecomplex*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
          */
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   

}


//------------------------------------------------------------------------------
// Matrix inversion
//------------------------------------------------------------------------------
void matInverse (float *mat, int nRows, int nCols)
{
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (&nRows, &nCols, (float *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (nRows, nCols, mat, r, pivots, err);
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (&nRows, &nCols, mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (nRows, nCols, mat, r, pivots, &err);
#  else 
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      sgetrf_ (&r, &c, (real *)mat, &r, pivots, &err);
#  endif      
   if (err)  {   // TODO: throw execption
      std::cout << "SxMatrix<T>::inverse: Error in SGETRF: " << err <<std::endl;
      delete [] pivots;
      SX_EXIT;
   }

// --- (2) diagonalization
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (&r, (float *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (r, mat, r, pivots, work, lWork, err);
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (&r, mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      float *work = NULL; // ACML doesn't use work
      sgetri (r, mat, r, pivots, &err);
#  else
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      real *work = new real [lWork];
      sgetri_ (&r, (real *)mat, &r, pivots, work, &lWork, &err);      
#  endif      

   if ( err )  {  // TODO: throw execption
      std::cout << "SxMatrix<T>::inverse: Error in SGETRI: " << err <<std::endl;
      delete [] work; delete [] pivots;
      SX_EXIT;
   }

   delete [] work; delete [] pivots;
}


void matInverse (double *mat, int nRows, int nCols)
{
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (&nRows, &nCols, (double *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (nRows, nCols, mat, r, pivots, err);
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (&nRows, &nCols, mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (nRows, nCols, mat, r, pivots, &err);
#  else 
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      dgetrf_ (&r, &c, (doublereal *)mat, &r, pivots, &err);
#  endif      

   if (err)  {   // TODO: throw execption
      std::cout << "SxMatrix<T>::inverse: Error in DGETRF: " << err <<std::endl;
      delete [] pivots;
      SX_EXIT;
   }

   // --- (2) diagonalization
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (&nRows, (double *)mat, &nRows, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (nRows, mat, nRows, pivots, work, lWork, err);
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (&nRows, mat, &nRows, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      double *work = NULL;  // ACML doesn't use work
      dgetri (nRows, mat, nRows, pivots, &err);
#  else 
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      doublereal *work = new doublereal [lWork];
      dgetri_ (&r, (doublereal *)mat, &r, pivots, work, &lWork, &err);
#  endif      

   if ( err )  {  // TODO: throw execption
      std::cout << "SxMatrix<T>::inverse: Error in DGETRI: " << err <<std::endl;
      delete [] work; delete [] pivots;
      SX_EXIT;
   }

   delete [] work; delete [] pivots;
}


void matInverse (SxComplex8 *mat, int nRows, int nCols)
{
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (&r, &c, (complex8_t *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r, c, err=0;
      // map complex inversion to real inversion
      //             +----+----+     +-----+-----+
      // +----+      | Re |-Im |     | re1 | im1 |
      // |cmlx|  =>  +----+----+  =  +-----+-----+
      // +----+      | Im | Re |     | im2 | re2 |
      //             +----+----+     +-----+-----+
      SxComplex8 *matPtr = mat;
      float *realMat = new float [ (2*nRows)*(2*nCols) ];
      float *re1Ptr  =  realMat;
      float *im1Ptr  = &realMat[nCols];
      float *im2Ptr  = &realMat[2*nRows*nCols];
      float *re2Ptr  = &realMat[2*nRows*nCols + nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, re2Ptr+=nCols,
                              im1Ptr+=nCols, im2Ptr+=nCols)  
      {
         for (c=0; c<nCols; c++, matPtr++)  {
            *re1Ptr++ = *re2Ptr++ = matPtr->re;
            *im1Ptr++ =            -matPtr->im;
            *im2Ptr++ =             matPtr->im;
         }
      }
      // --- compute inverse of real helper matrix
      matInverse (realMat, 2*nRows, 2*nCols);
      // construct complex result
      // +----+----+     +-----+-----+
      // | Re | *  |     | re1 |  *  |    +----+
      // +----+----+  =  +-----+-----+ => |cmlx|
      // | Im | *  |     | im2 |  *  |    +----+
      // +----+----+     +-----+-----+
      matPtr = mat;
      re1Ptr =  realMat;
      im2Ptr = &realMat[2*nRows*nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, im2Ptr+=nCols)  {
         for (c=0; c<nCols; c++, matPtr++)  {
            matPtr->re = *re1Ptr++;
            matPtr->im = *im2Ptr++;
         }
      }
      delete [] realMat;
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (&r, &c, (MKL_Complex8 *)mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (r, c, (complex *)mat, r, pivots, &err);
#  else  
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      cgetrf_ (&r, &c, (complex *)mat, &r, pivots, &err);
#  endif

#  ifndef USE_ESSL      
      if (err)  {   // TODO: throw execption
         std::cout << "SxMatrix<T>::inverse: Error in CGETRF: "<<err<<std::endl;
         delete [] pivots;
         SX_EXIT;
      }
#  endif 

   // --- (2) diagonalization
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex8_t *work = new complex8_t [lWork];
      cgetri (&r, (complex8_t *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      // empty
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      MKL_Complex8 *work = new MKL_Complex8 [lWork];
      cgetri (&r, (MKL_Complex8 *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      complex *work = NULL;  // ACML doesn't use work
      cgetri (r, (complex *)mat, r, pivots, &err);
#  else  
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex *work = new complex [lWork];
      cgetri_ (&r, (complex *)mat, &r, pivots, work, &lWork, &err);
#  endif      

#  ifndef USE_ESSL
      if ( err )  {  // TODO: throw execption
         std::cout << "SxMatrix<T>::inverse: Error in CGETRI: "<<err<<std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}


void matInverse (SxComplex16 *mat, int nRows, int nCols)
{
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (&r, &c, (complex16_t *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r, c, err=0;
      // map complex inversion to real inversion
      //             +----+----+     +-----+-----+
      // +----+      | Re |-Im |     | re1 | im1 |
      // |cmlx|  =>  +----+----+  =  +-----+-----+
      // +----+      | Im | Re |     | im2 | re2 |
      //             +----+----+     +-----+-----+
      SxComplex16 *matPtr = mat;
      double *realMat = new double [ (2*nRows)*(2*nCols) ];
      double *re1Ptr  =  realMat;
      double *im1Ptr  = &realMat[nCols];
      double *im2Ptr  = &realMat[2*nRows*nCols];
      double *re2Ptr  = &realMat[2*nRows*nCols + nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, re2Ptr+=nCols,
                              im1Ptr+=nCols, im2Ptr+=nCols)  
      {
         for (c=0; c<nCols; c++, matPtr++)  {
            *re1Ptr++ = *re2Ptr++ = matPtr->re;
            *im1Ptr++ =            -matPtr->im;
            *im2Ptr++ =             matPtr->im;
         }
      }
      // --- compute inverse of real helper matrix
      matInverse (realMat, 2*nRows, 2*nCols);
      // construct complex result
      // +----+----+     +-----+-----+
      // | Re | *  |     | re1 |  *  |    +----+
      // +----+----+  =  +-----+-----+ => |cmlx|
      // | Im | *  |     | im2 |  *  |    +----+
      // +----+----+     +-----+-----+
      matPtr = mat;
      re1Ptr =  realMat;
      im2Ptr = &realMat[2*nRows*nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, im2Ptr+=nCols)  {
         for (c=0; c<nCols; c++, matPtr++)  {
            matPtr->re = *re1Ptr++;
            matPtr->im = *im2Ptr++;
         }
      }
      delete [] realMat;
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (&r, &c, (MKL_Complex16 *)mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (r, c, (doublecomplex *)mat, r, pivots, &err);
#  else 
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      zgetrf_ (&r, &c, (doublecomplex *)mat, &r, pivots, &err);
#  endif   

#  ifndef USE_ESSL      
      if (err)  {   // TODO: throw execption
         std::cout << "SxMatrix<T>::inverse: Error in ZGETRF: "<<err<<std::endl;
         delete [] pivots;
         SX_EXIT;
      }
#  endif

   // --- (2) diagonalization
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex16_t *work = new complex16_t [lWork];
      zgetri (&r, (complex16_t *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      // empty
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      MKL_Complex16 *work = new MKL_Complex16 [lWork];
      zgetri (&r, (MKL_Complex16 *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      doublecomplex *work = NULL; // ACML doesn't use work
      zgetri (r, (doublecomplex *)mat, r, pivots, &err);
#  else 
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      doublecomplex *work = new doublecomplex [lWork];
      zgetri_ (&r, (doublecomplex *)mat, &r, pivots, work, &lWork, &err);
#  endif   

#  ifndef USE_ESSL      
      if ( err )  {  // TODO: throw execption
         std::cout << "SxMatrix<T>::inverse: Error in ZGETRI: "<<err<<std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}



void matInverseTri (float * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

void matInverseTri (double *mat, int nRows, enum UPLO uplo)
{
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
   // --- (1) Factorization: A = U*D*U^t 
#  if   defined (USE_VECLIB)
      SX_EXIT; // not yet implemented
#  elif defined ( USE_ESSL)
      SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
//     int err = 0;
//     SX_EXIT; // not yet implemented
      // change by khr
      int r = nRows, err = 0;
      int *pivots = new int[r];
      dsptrf (&uploChar, &r, mat, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, err = 0;
      int *pivots = new int[r];
      dsptrf (uploChar, r, mat, pivots, &err);
#  else
      integer r = (integer)nRows, err = 0;
      integer *pivots = new integer[r];
      dsptrf_ (&uploChar, &r, (doublereal *)mat, pivots, &err);
#  endif    
      if ( err )  {  // TODO: throw execption
         std::cout << "SxMatrix<T>::inverseTri: Error in DSPTRF: "
                   << err << std::endl;
         SX_EXIT;
      }

   // --- (2) Inverse of A
#  if   defined (USE_VECLIB)
       SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
       SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
//       double *work = NULL, *pivots=NULL;
//       SX_EXIT; // not yet implemented
      // change by khr
      double *work = new double [r];
      dsptri (&uploChar, &r, mat, pivots, work, &err);
#  elif defined (USE_ACML)
      double *work = NULL; // ACML doesn't use work
      dsptri (uploChar, r, mat, pivots, &err);
#  else
      doublereal *work = new doublereal [r];
      dsptri_ (&uploChar, &r, (doublereal *)mat, pivots, work, &err);
#  endif      

#  ifndef USE_ESSL      
      if ( err )  {  // TODO: throw execption
         std::cout << "SxMatrix<T>::inverseTri: Error in DSPTRI: "
                   << err << std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}

void matInverseTri (SxComplex8 * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

void matInverseTri (SxComplex16 * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

//------------------------------------------------------------------------------
// Linear equation solver (least square based)
//------------------------------------------------------------------------------
void solveLinEq (float *mat, int nRows, int nCols, float *b, int bCols)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      int iopt = 2; //compute singulare Values, V and U^TB
      int m = nRows;
      int n = nCols;
      int nrhs = bCols;
      int lda = m;
      int ldb = m;
      char transa = 'N';
      int naux = 0; //ESSL choose size of work array dynamically
      float *aux = NULL; //workarray is ignored
      int info;
      float *s = new float [n];
      float *x = new float [n*nrhs];
      float tau = 1e-10; // error tolerance for zero
      sgesvf(iopt,mat,lda,b,ldb,nrhs,s,m,n,aux,naux);
      sgesvs (mat,n,b,n,nrhs,s,x,n,m,n,tau);
      // results are now in x; copy to b
      for(int i = 0; i < n*nrhs; i++)   {
         b[i] = x[i];
      }

      delete [] s; delete [] x;

      //TODO TEST
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      // for ilaenv
      MKL_INT iSpec      = 9;
      char *name   = const_cast<char*>("SGELSD");
      char *opts   = const_cast<char*>(" ");
      MKL_INT n1         = 0;
      MKL_INT n2         = 0;
      MKL_INT n3         = 0;
      MKL_INT n4         = 0;
      // for dgelsd
      MKL_INT m          = nRows;
      MKL_INT n          = nCols;
      MKL_INT minmn      = n < m ? n : m;
      MKL_INT nrhs       = bCols;
      MKL_INT lda        = m;
      MKL_INT ldb        = m;
      float *s           = new float [minmn];
      float rcond        = -1.;
      MKL_INT rank;
      MKL_INT SMLSIZ     = ilaenv(&iSpec,name,opts,&n1,&n2,&n3,&n4);
      MKL_INT logVal     = MKL_INT(log( 1.0*n/(SMLSIZ+1))/log(2.0)) + 1;
      MKL_INT NLVL       = 0 < logVal ? logVal : 0;
      MKL_INT lwork      = -1;
      float autolwork;
      MKL_INT liwork     = 3 * minmn * NLVL + 11 * minmn;
      MKL_INT *iwork     = new MKL_INT [liwork];
      MKL_INT info;
      // determine optimal lwork
//      dgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      sgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = MKL_INT(autolwork+0.5);
      float *work        = new float [lwork];

      sgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
      //TODO TEST
      SX_EXIT;
#  elif defined (USE_ACML)
      int m        = nRows;
      int n        = nCols;
      int nrhs     = bCols;
      int lda      = m;
      int ldb      = m;
      int minmn    = n < m ? n : m;
      float *s    = new float [minmn];
      float rcond = -1.;
      int rank;
      int info;
      sgelsd (m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&info);

      delete [] s; 
      //TODO TEST
      SX_EXIT;
#  else 
      // for ilaenv
      integer iSpec = 9;
      char *name    = const_cast<char*>("SGELSD");
      char *opts    = const_cast<char*>(" ");
      integer n1    = 0;
      integer n2    = 0;
      integer n3    = 0;
      integer n4    = 0;
      ftnlen lname  = 6;
      ftnlen lopts  = 1;
      // for sgelsd
      integer m          = nRows;
      integer n          = nCols;
      integer minmn      = n < m ? n : m;
      integer nrhs       = bCols;
      integer lda        = m;
      integer ldb        = m;
      float *s          = new float [minmn];
      float rcond       = -1.;
      integer rank;
#ifdef MACOSX
#  if ( DIST_VERSION_L >= 1070L )
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4);
#  else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#  endif /* DIST_VERSION_L */
#else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#endif /* MACOSX */
      integer logVal     = integer(log( double(n)/double(SMLSIZ+1))/log(2.0)) + 1;
      integer NLVL       = 0 < logVal ? logVal : 0;
      integer lwork      = -1;
      float autolwork;
      integer liwork     = 3 * minmn * NLVL + 11 * minmn;
      integer *iwork     = new integer [liwork];
      integer info;
      // determine optimal lwork
      sgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = integer(autolwork+0.5);
      float *work       = new float [lwork];

      sgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  endif      
}

void solveLinEq (double *mat, int nRows, int nCols, double *b, int bCols)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      int iopt = 2; //compute singulare Values, V and U^TB
      int m = nRows;
      int n = nCols;
      int nrhs = bCols;
      int lda = m;
      int ldb = m;
      char transa = 'N';
      int naux = 0; //ESSL choose size of work array dynamically
      double *aux = NULL; //workarray is ignored
      int info;
      double *s = new double [n];
      double *x = new double [n*nrhs];
      double tau = 1e-10; // error tolerance for zero
      dgesvf(iopt,mat,lda,b,ldb,nrhs,s,m,n,aux,naux);
      dgesvs (mat,n,b,n,nrhs,s,x,n,m,n,tau);
      // results are now in x; copy to b
      for(int i = 0; i < n*nrhs; i++)   {
         b[i] = x[i];
      }

      delete [] s; delete [] x;

      //TODO TEST
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      // for ilaenv
      MKL_INT iSpec = 9;
      char *name    = const_cast<char*>("DGELSD");
      char *opts    = const_cast<char*>(" ");
      MKL_INT n1    = 0;
      MKL_INT n2    = 0;
      MKL_INT n3    = 0;
      MKL_INT n4    = 0;
      // for dgelsd
      MKL_INT m          = nRows;
      MKL_INT n          = nCols;
      MKL_INT minmn      = n < m ? n : m;
      MKL_INT nrhs       = bCols;
      MKL_INT lda        = m;
      MKL_INT ldb        = m;
      double *s          = new double [minmn];
      double rcond       = -1.;
      MKL_INT rank;
      MKL_INT SMLSIZ     = ilaenv(&iSpec,name,opts,&n1,&n2,&n3,&n4);
      MKL_INT logVal     = MKL_INT(log( 1.0*n/(SMLSIZ+1))/log(2.0)) + 1;
      MKL_INT NLVL       = 0 < logVal ? logVal : 0;
      MKL_INT lwork      = -1;
      double  autolwork;
      MKL_INT liwork     = 3 * minmn * NLVL + 11 * minmn;
      MKL_INT *iwork     = new MKL_INT [liwork];
      MKL_INT info;
      // determine optimal lwork
      dgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = MKL_INT(autolwork+0.5);
      double *work       = new double [lwork];

      dgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  elif defined (USE_ACML)
      int m        = nRows;
      int n        = nCols;
      int nrhs     = bCols;
      int lda      = m;
      int ldb      = m;
      int minmn    = n < m ? n : m;
      double *s    = new double [minmn];
      double rcond = -1.;
      int rank;
      int info;
      dgelsd (m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&info);

      delete [] s; 
      //TODO TEST
      SX_EXIT;
#  else 
      // for ilaenv
      integer iSpec = 9;
      char *name    = const_cast<char*>("DGELSD");
      char *opts    = const_cast<char*>(" ");
      integer n1    = 0;
      integer n2    = 0;
      integer n3    = 0;
      integer n4    = 0;
      ftnlen lname  = 6;
      ftnlen lopts  = 1;
      // for dgelsd
      integer m          = nRows;
      integer n          = nCols;
      integer minmn      = n < m ? n : m;
      integer nrhs       = bCols;
      integer lda        = m;
      integer ldb        = m;
      double *s          = new double [minmn];
      double rcond       = -1.;
      integer rank;
#ifdef MACOSX      
#  if ( DIST_VERSION_L >= 1070L )
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4);
#  else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#  endif /* DIST_VERSION_L */
#else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#endif /* MACOSX */
      integer logVal     = integer(log( double(n)/double(SMLSIZ+1))/log(2.0)) + 1;
      integer NLVL       = 0 < logVal ? logVal : 0;
      integer lwork      = -1;
      double  autolwork;
      integer liwork     = 3 * minmn * NLVL + 11 * minmn;
      integer *iwork     = new integer [liwork];
      integer info;
      // determine optimal lwork
      dgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = integer(autolwork+0.5);
      double *work       = new double [lwork];

      dgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  endif      
}

void solveLinEq (SxComplex8 * /*mat*/, int /*nRows*/, int /*nCols*/,
                 SxComplex8 * /*b*/, int /*bCols*/)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_ACML)
      //TODO
      SX_EXIT;
#  endif      
}

void solveLinEq (SxComplex16 * /*mat*/, int /*nRows*/, int /*nCols*/,
                 SxComplex16 * /*b*/, int /*bCols*/)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_ACML)
      //TODO 
      SX_EXIT;
#  else 
      //TODO
      SX_EXIT;
#  endif      
}

//------------------------------------------------------------------------------
// Eigensolver
//------------------------------------------------------------------------------
int matEigensolver (SxComplex8 *eigVals, float *eigVecs,
                    float *inMat, int n, EIGCMD cmd, int size)
{
#  if   defined (USE_VECLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_INTEL_MKL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else 
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 4*n : size;
      integer workDim = lWork;
#  endif   
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      float *work = new float [workDim];
      float *epsRe = new float [n], *epsIm = new float [n];
      float *eigVecLeft = new float [ldVecLeft];
      sgeev (&jobvl, &jobvr,
             &rank, (float *)inMat, &rank, epsRe, epsIm, eigVecLeft, 
             &ldVecLeft, (float *)eigVecs, &rank, work, &lWork, &info, 0, 0);
#  elif defined (USE_ESSL)
      complex<float> *vecs = new complex<float> [rank*rank];
      int iOpt = 1;   // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         sgeev (iOpt, inMat, rank, 
                (complex<float> *)eigVals, vecs, rank, 
                NULL, rank, NULL, 0);
         complex<float> *srcPtr = vecs;
         float *dstPtr = eigVecs;
         int i, len = n*n;
         for (i=0; i < len; i++, srcPtr++)  {
            *dstPtr++ = real(*srcPtr);
         }
         // --- normalize eigenvectors
         float c;
         for (i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
      delete [] vecs;
      return 0;
#  elif defined (USE_INTEL_MKL)
      float *work = new float [workDim];
      float *epsRe = new float [n], *epsIm = new float [n];
      float *eigVecLeft = new float [ldVecLeft];
      sgeev (&jobvl, &jobvr, 
             &rank, inMat, &rank, epsRe, epsIm, eigVecLeft, 
             &ldVecLeft, eigVecs, &rank, work, &lWork, &info);
#  elif defined (USE_ACML)
      float *work = NULL;  // ACML doesn't use work
      float *epsRe = new float [n], *epsIm = new float [n];
      float *eigVecLeft = new float [ldVecLeft];
      sgeev (jobvl, jobvr, 
             rank, inMat, rank, epsRe, epsIm, eigVecLeft, 
             ldVecLeft, eigVecs, rank, &info);
#  else
      real *work = new real [workDim];
      real *epsRe = new real[n], *epsIm = new real[n];
      real *eigVecLeft = new real [ldVecLeft];
      sgeev_ (&jobvl, &jobvr,
              &rank, (real *)inMat, &rank, epsRe, epsIm, eigVecLeft, 
              &ldVecLeft, (real *)eigVecs, &rank, work, &lWork, &info);
#  endif      

#  ifndef USE_ESSL      
      float *ptr=(float *)eigVals;
      float *rePtr=(float *)epsRe, *imPtr=(float *)epsIm;
      for (int i=0; i < n; i++)  {
         *ptr++ = *rePtr++;  // eigVals[i].re = epsRe[i];
         *ptr++ = *imPtr++;  // eigVals[i].im = epsIm[i];
      }
      delete [] eigVecLeft; delete [] epsIm; delete [] epsRe; delete [] work;

      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in SGEEV: " << info << std::endl;
         SX_EXIT;
      }
      return (cmd == OptSize)  ?  (int)work[0] : 0;
#  endif  

}


int matEigensolver (SxComplex16 *eigVals, double *eigVecs,
                    double *inMat, int n, EIGCMD cmd, int size)
{
#  if   defined (USE_VECLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_INTEL_MKL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else 
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 4*n : size;
      integer workDim = lWork;
#  endif      
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      double *work = new double [workDim];
      double *epsRe = new double[n], *epsIm = new double[n];
      double *eigVecLeft = new double [ldVecLeft];
      dgeev (&jobvl, &jobvr,
              &rank, (double *)inMat, &rank, epsRe, epsIm, eigVecLeft, 
              &ldVecLeft, (double *)eigVecs, &rank, work, &lWork, &info, 0, 0);
#  elif defined (USE_ESSL)
      complex<double> *vecs = new complex<double> [rank*rank];
      int iOpt = 1;   // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         dgeev (iOpt, inMat, rank, 
                (complex<double> *)eigVals, vecs, rank, 
                NULL, rank, NULL, 0);
         complex<double> *srcPtr = vecs;
         double *dstPtr = eigVecs;
         int i, len = n*n;
         for (i=0; i < len; i++, srcPtr++)  {
            *dstPtr++ = real(*srcPtr);
         }
         // --- normalize eigenvectors
         double c;
         for (i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
      delete [] vecs;
      return 0;
#  elif defined (USE_INTEL_MKL)
      double *work       = new double [workDim];
      double *epsRe      = new double [n], *epsIm = new double [n];
      double *eigVecLeft = new double [ldVecLeft];
      dgeev (&jobvl, &jobvr, 
             &rank, inMat, &rank, epsRe, epsIm, eigVecLeft, 
             &ldVecLeft, eigVecs, &rank, work, &lWork, &info);      
#  elif defined (USE_ACML)
      double *work       = NULL; // ACML doesn't use work
      double *epsRe      = new double [n], *epsIm = new double [n];
      double *eigVecLeft = new double [ldVecLeft];
      dgeev (jobvl, jobvr, 
             rank, inMat, rank, epsRe, epsIm, eigVecLeft, 
             ldVecLeft, eigVecs, rank, &info);      
#  else  
      doublereal *work = new doublereal [workDim];
      doublereal *epsRe = new doublereal[n], *epsIm = new doublereal[n];
      doublereal *eigVecLeft = new doublereal [ldVecLeft];
      dgeev_ (&jobvl, &jobvr,
              &rank, (doublereal *)inMat, &rank, epsRe, epsIm, eigVecLeft, 
              &ldVecLeft, (doublereal *)eigVecs, &rank, work, &lWork, &info);
#  endif      

#  ifndef USE_ESSL      
      double *ptr=(double *)eigVals;
      double *rePtr=(double *)epsRe, *imPtr=(double *)epsIm;
      for (int i=0; i < n; i++)  {
         *ptr++ = *rePtr++;  // eigVals[i].re = epsRe[i];
         *ptr++ = *imPtr++;  // eigVals[i].im = epsIm[i];
      }
      delete [] eigVecLeft; delete [] epsIm; delete [] epsRe; delete [] work;

      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in DGEEV: " << info << std::endl;
         SX_EXIT;
      }
      return (cmd == OptSize)  ?  (int)work[0] : 0;
#  endif

}


int matEigensolver (SxComplex8 *eigVals, SxComplex8 *eigVecs, 
                    SxComplex8 *inMat, int n, EIGCMD cmd, int size)
{
#  if   defined (USE_VECLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_INTEL_MKL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else  
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 2*n : size;
      integer workDim = lWork;
#  endif   
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      complex8_t *eigVecLeft = new complex8_t [ldVecLeft];
      complex8_t *work       = new complex8_t [workDim];
      float      *rWork      = new float      [workDim];
      cgeev (&jobvl, &jobvr,
             &rank, (complex8_t *)inMat, &rank, 
             (complex8_t *)eigVals, eigVecLeft, 
             &ldVecLeft, (complex8_t *)eigVecs, &rank, work, &lWork, 
             rWork, &info, 0, 0);
      if (cmd == OptSize) lWork = (int)work[0].re;
      delete [] rWork;  delete [] work; delete [] eigVecLeft; 
#  elif defined (USE_ESSL)
      int iOpt = 1;  // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         cgeev (iOpt, (complex<float> *)inMat, rank, 
                (complex<float> *)eigVals, (complex<float> *)eigVecs, rank, 
                NULL, rank, NULL, 0);
         // --- normalize eigenvectors
         double c;
         for (int i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
#  elif defined (USE_INTEL_MKL)
      MKL_Complex8 *eigVecLeft = new MKL_Complex8 [ldVecLeft];
      MKL_Complex8 *work       = new MKL_Complex8 [workDim];
      float        *rWork      = new float        [workDim];
      cgeev (&jobvl, &jobvr, 
             &rank, (MKL_Complex8 *)inMat, &rank, 
             (MKL_Complex8 *)eigVals, eigVecLeft, 
             &ldVecLeft, (MKL_Complex8 *)eigVecs, &rank, work, &lWork, 
             rWork, &info);      
      if (cmd == OptSize) lWork = (int)work[0].real;
      delete [] rWork;  delete [] work; delete [] eigVecLeft; 
#  elif defined (USE_ACML)
      complex *eigVecLeft = new complex [ldVecLeft];
      cgeev (jobvl, jobvr, 
             rank, (complex *)inMat, rank, 
             (complex *)eigVals, eigVecLeft, 
             ldVecLeft, (complex *)eigVecs, rank, &info);      
      delete eigVecLeft;
#  else  
      complex *eigVecLeft = new complex [ldVecLeft];
      complex *work       = new complex [workDim];
      real    *rWork      = new real    [workDim];
      cgeev_ (&jobvl, &jobvr,
              &rank, (complex *)inMat, &rank, 
              (complex *)eigVals, eigVecLeft, 
              &ldVecLeft, (complex *)eigVecs, &rank, work, &lWork, 
              rWork, &info);
      if (cmd == OptSize) lWork = (int)work[0].r;
      delete [] rWork;  delete [] work; delete [] eigVecLeft; 
#  endif      

#  ifndef USE_ESSL
      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in CGEEV: " << info << std::endl;
         SX_EXIT;
      }
#  endif

#  if   defined (USE_VECLIB)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ESSL)
      return 0;
#  elif defined (USE_INTEL_MKL)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ACML)
      return 0;
#  else  
      return (cmd == OptSize)  ?  (int)lWork : 0;
#  endif      
}


int matEigensolver (SxComplex16 *eigVals, SxComplex16 *eigVecs, 
                    SxComplex16 *inMat, int n, EIGCMD cmd, int size)
{
#  if   defined (USE_VECLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_INTEL_MKL)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else 
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 2*n : size;
      integer workDim = !size ? 2*n : size;
#  endif
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      complex16_t *eigVecLeft = new complex16_t [ldVecLeft];
      complex16_t *work       = new complex16_t [workDim];
      double      *rWork      = new double      [workDim];
      zgeev (&jobvl, &jobvr,
             &rank, (complex16_t *)inMat, &rank, 
             (complex16_t *)eigVals, eigVecLeft, 
             &ldVecLeft, (complex16_t *)eigVecs, &rank, work, &lWork, 
             rWork, &info, 0, 0);
      if (cmd == OptSize) lWork = (int)work[0].re;
      delete [] rWork;  delete [] work; delete [] eigVecLeft; 
#  elif defined (USE_ESSL)
      int iOpt = 1;   // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         zgeev (iOpt, (complex<double> *)inMat, rank, 
                (complex<double> *)eigVals, (complex<double> *)eigVecs, rank, 
                NULL, rank, NULL, 0);
         // --- normalize eigenvectors
         double c;
         for (int i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
#  elif defined (USE_INTEL_MKL)
      MKL_Complex16 *eigVecLeft = new MKL_Complex16 [ldVecLeft];
      MKL_Complex16 *work       = new MKL_Complex16 [workDim];
      double        *rWork      = new double        [workDim];
      zgeev (&jobvl, &jobvr, 
             &rank, (MKL_Complex16 *)inMat, &rank, 
             (MKL_Complex16 *)eigVals, eigVecLeft, 
             &ldVecLeft, (MKL_Complex16 *)eigVecs, &rank, work, &lWork, 
             rWork, &info);
      if (cmd == OptSize) lWork = (int)work[0].real;
      delete [] rWork;
      delete [] work;
      delete [] eigVecLeft;
#  elif defined (USE_ACML)
      doublecomplex *eigVecLeft = new doublecomplex [ldVecLeft];
      zgeev (jobvl, jobvr, 
             rank, (doublecomplex *)inMat, rank, 
             (doublecomplex *)eigVals, eigVecLeft, 
             ldVecLeft, (doublecomplex *)eigVecs, rank, &info);           
      delete [] eigVecLeft;
#  else
      doublecomplex *eigVecLeft = new doublecomplex [ldVecLeft];
      doublecomplex *work       = new doublecomplex [workDim];
      doublereal    *rWork      = new doublereal    [workDim];
      zgeev_ (&jobvl, &jobvr,
              &rank, (doublecomplex *)inMat, &rank, 
              (doublecomplex *)eigVals, eigVecLeft, 
              &ldVecLeft, (doublecomplex *)eigVecs, &rank, work, &lWork, 
              rWork, &info);
      if (cmd == OptSize) lWork = (int)work[0].r;
      delete [] rWork;  delete [] work; delete [] eigVecLeft; 
#  endif      

#  ifndef USE_ESSL
     if ( info )  {  // TODO: throw exception
        std::cout << "matEigensolver: Error in DGEEV: " << info << std::endl;
        SX_EXIT;
     }
#  endif

#  if   defined (USE_VECLIB)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ESSL)
      return 0;
#  elif defined (USE_INTEL_MKL)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ACML)
      return 0;
#  else  
      return (cmd == OptSize)  ?  (int)lWork : 0;
#  endif      
}


//------------------------------------------------------------------------------
// Eigensolver - tridiagonal matrices
//------------------------------------------------------------------------------
void matEigensolverTri (float *eigVals, float *eigVecs,
                        float *inMat, int n, enum UPLO uplo, 
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 3*n;
      float *work  = new float [workDim];
      sspev (&jobz, &uploChar, 
             &rank, (float *)inMat, (float *)eigVals, 
             (float *)eigVecs, &ldEigVecs,
              work, &info, 0, 0);
      delete [] work;  
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      sspev (iOpt, inMat, eigVals, 
             eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 3*n;
      float *work  = new float [workDim];
      sspev (&jobz, &uploChar, 
             &rank, (float *)inMat, (float *)eigVals, 
             (float *)eigVecs, &ldEigVecs,
              work, &info);
      delete [] work;  
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      sspev (jobz, uploChar, 
             rank, inMat, eigVals, 
             eigVecs, ldEigVecs, &info);
#  else  
      integer rank = (integer)n;
      integer ldEigVecs = rank; 
      integer info = 0;
      integer workDim  = 3*n;
      real *work  = new real [workDim];
      sspev_ (&jobz, &uploChar, 
              &rank, (real *)inMat, (real *)eigVals, 
              (real *)eigVecs, &ldEigVecs,
               work, &info);
      delete [] work;  
#  endif      
   
   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in SSPEV: " << info << std::endl;
      SX_EXIT;
   }
}


void matEigensolverTri (double *eigVals, double *eigVecs,
                        double *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 3*n;
      double *work  = new double [workDim];
      dspev (&jobz, &uploChar, 
             &rank, (double *)inMat, (double *)eigVals, 
             (double *)eigVecs, &ldEigVecs,
              work, &info, 0, 0);
      delete [] work;  
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      dspev (iOpt, inMat, eigVals, 
             eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 3*n;
      double *work  = new double [workDim];
      dspev (&jobz, &uploChar, 
             &rank, (double *)inMat, (double *)eigVals, 
             (double *)eigVecs, &ldEigVecs,
              work, &info);
      delete [] work;  
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      dspev (jobz, uploChar, 
             rank, inMat, eigVals, 
             eigVecs, ldEigVecs, &info);
#  else
      integer rank = (integer)n;
      integer ldEigVecs = rank; 
      integer info = 0;
      integer workDim  = 3*n;
      doublereal *work  = new doublereal [workDim];
      dspev_ (&jobz, &uploChar, 
              &rank, (doublereal *)inMat, (doublereal *)eigVals, 
              (doublereal *)eigVecs, &ldEigVecs,
               work, &info);
      delete [] work;  
#  endif      
   
   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in DSPEV: " << info << std::endl;
      SX_EXIT;
   }
}


void matEigensolverTri (float *eigVals, SxComplex8 *eigVecs,
                        SxComplex8 *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      complex8_t *work = new complex8_t [workDim];
      float    *rWork  = new float    [rWorkDim];
      chpev (&jobz, &uploChar, 
             &rank, (complex8_t *)inMat, (float *)eigVals, 
             (complex8_t *)eigVecs, &ldEigVecs,
              work, rWork, &info, 0, 0);
      delete [] rWork; delete [] work;  
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      chpev (iOpt, (complex<float> *)inMat, eigVals, 
             (complex<float> *)eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      MKL_Complex8 *work = new MKL_Complex8 [workDim];
      float       *rWork = new float       [rWorkDim];
      chpev (&jobz, &uploChar, 
             &rank, (MKL_Complex8 *)inMat, (float *)eigVals, 
             (MKL_Complex8 *)eigVecs, &ldEigVecs,
              work, rWork, &info);
      delete [] rWork; delete [] work;  
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      chpev (jobz, uploChar, 
             rank, (complex *)inMat, (float *)eigVals, 
             (complex *)eigVecs, ldEigVecs, &info);
#  else  
      integer rank = (integer)n;
      integer ldEigVecs = rank; 
      integer info = 0;
      integer workDim  = 2*n-1;
      integer rWorkDim = 3*n-2;
      complex *work       = new complex [workDim];
      real    *rWork      = new real    [rWorkDim];
      chpev_ (&jobz, &uploChar, 
              &rank, (complex *)inMat, (real *)eigVals, 
              (complex *)eigVecs, &ldEigVecs,
               work, rWork, &info);
      delete [] rWork; delete [] work;  
#  endif      
   
   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in CHPEV: " << info << std::endl;
      SX_EXIT;
   }
}


void matEigensolverTri (double *eigVals, SxComplex16 *eigVecs,
                        SxComplex16 *inMat, int n, enum UPLO uplo, 
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      complex16_t *work = new complex16_t [workDim];
      double    *rWork  = new double[rWorkDim];
      zhpev (&jobz, &uploChar, 
             &rank, (complex16_t *)inMat, (double *)eigVals, 
             (complex16_t *)eigVecs, &ldEigVecs,
              work, rWork, &info, 0, 0);
      delete [] rWork; delete [] work;  
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      zhpev (iOpt, (complex<double> *)inMat, eigVals, 
             (complex<double> *)eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      MKL_Complex16 *work = new MKL_Complex16 [workDim];
      double       *rWork = new double[rWorkDim];
      zhpev (&jobz, &uploChar, 
             &rank, (MKL_Complex16 *)inMat, (double *)eigVals, 
             (MKL_Complex16 *)eigVecs, &ldEigVecs,
              work, rWork, &info);
      delete [] rWork; delete [] work;  
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank; 
      int info = 0;
      zhpev (jobz, uploChar, 
             rank, (doublecomplex *)inMat, (double *)eigVals, 
             (doublecomplex *)eigVecs, ldEigVecs, &info);
#  else  
      integer rank = (integer)n;
      integer ldEigVecs = rank; 
      integer info = 0;
      integer workDim  = 2*n-1;
      integer rWorkDim = 3*n-2;
      doublecomplex *work       = new doublecomplex [workDim];
      doublereal    *rWork      = new doublereal    [rWorkDim];
      zhpev_ (&jobz, &uploChar, 
              &rank, (doublecomplex *)inMat, (doublereal *)eigVals, 
              (doublecomplex *)eigVecs, &ldEigVecs,
               work, rWork, &info);
      delete [] rWork; delete [] work;  
#  endif   
   
   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in ZHPEV: " << info << std::endl;
      SX_EXIT;
   }
}

// -------------------------- SX_IGNORE_THE_REST_OF_THE_FILE 
# endif
// -------------------------- 
