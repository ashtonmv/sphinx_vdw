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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <SxMathLib.h>
#include <SxError.h>
#include <SxRandom.h>
#include <SxVector.h>
#ifdef USE_FFTW
# include <fftw3.h>
#endif


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
    extern "C"  {
#      include <mkl.h>
    }
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
#endif

//template<class T> void  *SxMemory<T>::buffer[2];
//template<class T> size_t SxMemory<T>::nBytes[2];
//template<class T> long   SxMemory<T>::accessCounter[2];

//template class SX_EXPORT_MATH SxMemory<SxComplex8>;
//template class SX_EXPORT_MATH SxMemory<SxComplex16>;

void initSPHInXMath ()
{
   // --- setup memory manager
//   SxMemory<int>::init ();
//   SxMemory<float>::init ();
//   SxMemory<double>::init ();
//   SxMemory<SxComplex8>::init ();
//   SxMemory<SxComplex16>::init ();

   initRandom (23, 43, 34, 54);
// gethostname (name, sizeof (name));
// struct hostent *info = gethostbyname (name);
// initVectorClass (info->h_name, gethostid());
#  if defined USE_FFTW && defined USE_OPENMP
      if (fftw_init_threads () == 0)  {
         std::cout << "Error initializing FFTW threads" << std::endl;
         SX_EXIT;
      }
      fftw_plan_with_nthreads (SxUtil::getGlobalObj ().nProcs);
#  endif
}   

void destroySFHIngXMath ()
{
//   SxMemory<SxComplex16>::destroy ();
//   SxMemory<SxComplex8>::destroy ();
//   SxMemory<double>::destroy ();
//   SxMemory<float>::destroy ();
//   SxMemory<int>::destroy ();
}

//------------------------------------------------------------------------------
// in-place matrix multiplication
//------------------------------------------------------------------------------
void inPlaceRot (float *mat,  const float *rotMat, 
                 int nRows, int nCols)
{
   if (nRows <= 0 || nCols <=0) return;
   // make sure that mat and rotMat do not overlap
   SX_CHECK (   mat    >= rotMat + nCols * nCols
             || rotMat >= mat    + nRows * nCols,
             (mat-rotMat), nRows * nCols);
   int blockSize = 1024;
#if defined USE_ATLAS || defined USE_ACML || defined USE_INTEL_MKL
   const int minSize = 2 * blockSize;
   if (nRows < minSize || sizeof(float) * nRows * nCols < 102400)
#endif
   {
      float *workspace = new float[nRows * nCols];
      memcpy (workspace, mat, sizeof(float) * nRows * nCols);
      matmult (mat, workspace, rotMat, nRows, nCols, nCols);
      delete [] workspace;
      return;
   }
   
   SX_EXIT; // needs to be tested against implementation above
   
   const float one = 1.0, zero = 0.;
#if ! defined USE_ATLAS
   char noTrans = 'N';
#endif
// sketch for omp: #  pragma omp parallel
   {
      float *workspace;
      // catching the bad_alloc is essential for omp parallelism, because exceptions from omp
      // threads cannot be caught from the outside
      try {
         workspace = new float[blockSize * nCols];
      } catch (std::bad_alloc &err) {
         cout << "Memory allocation failed in inPlaceRot." << endl;
         cout << err.what () << endl;
         SX_EXIT;
      }
      int r, c, myBlockSize = blockSize;
// sketch for omp: #     pragma omp for
      for (r = 0; r < nRows; r += blockSize)  {
         if (r + blockSize > nRows)  {
            // reduce block size for last block
            myBlockSize = nRows - r;
         }
         // --- apply rotation matrix to current block, output to workspace
// --- critical point for not having omp:
//     * how do we prevent parallelism inside zgemm???
//     * does it hurt if zgemm makes more threads??
//     * maybe if omp is used inside zgemm, omp can deal with it...
#ifdef USE_ATLAS
         cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols, one, mat + r, nRows,
                      rotMat, nCols, zero, workspace, myBlockSize);
#elif defined USE_GOTO // khr: GotoBLAS, experimental!
         cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols, one, mat + r, nRows,
                      const_cast<float*>(rotMat), nCols, zero, workspace,
                      myBlockSize);
#elif defined USE_ACML
         sgemm (noTrans,noTrans, myBlockSize, nCols, nCols, 
                one, (mat + r), nRows, const_cast<float*>(rotMat), nCols,
                zero, workspace, myBlockSize);
#elif defined USE_INTEL_MKL
         sgemm (&noTrans, &noTrans, 
                &myBlockSize, &nCols, &nCols,
                &one, (mat + r), &nRows,
                rotMat, &nCols, &zero, workspace, &myBlockSize);
#else
         SX_EXIT;
#endif
         // copy result into mat
         for (c = 0; c < nCols; ++c)
            memcpy (mat + r + nRows * c, workspace + myBlockSize * c,
                    sizeof(float) * myBlockSize);
      }
      delete [] workspace;
   }
}

void inPlaceRot (double *mat,  const double *rotMat, int nRows, int nCols)
{
   if (nRows <= 0 || nCols <=0) return;
   // make sure that mat and rotMat do not overlap
   SX_CHECK (   mat    >= rotMat + nCols * nCols
             || rotMat >= mat    + nRows * nCols,
             (mat-rotMat), nRows * nCols);
   int blockSize = 1024;
#if defined USE_ATLAS || defined USE_ACML || defined USE_INTEL_MKL
   const int minSize = 2 * blockSize;
   if (nRows < minSize || sizeof(double) * nRows * nCols < 102400)
#endif
   {
      double *workspace = new double[nRows * nCols];
      memcpy (workspace, mat, sizeof(double) * nRows * nCols);
      matmult (mat, workspace, rotMat, nRows, nCols, nCols);
      delete [] workspace;
      return;
   }
   
   SX_EXIT; // needs to be tested against implementation above

   const double one = 1., zero = 0.;
#if ! defined USE_ATLAS
   char noTrans = 'N';
#endif
// sketch for omp: #  pragma omp parallel
   {
      double *workspace;
      // catching the bad_alloc is essential for omp parallelism, because exceptions from omp
      // threads cannot be caught from the outside
      try {
         workspace = new double[blockSize * nCols];
      } catch (std::bad_alloc &err) {
         cout << "Memory allocation failed in inPlaceRot." << endl;
         cout << err.what () << endl;
         SX_EXIT;
      }
      int r, c, myBlockSize = blockSize;
// sketch for omp: #     pragma omp for
      for (r = 0; r < nRows; r += blockSize)  {
         if (r + blockSize > nRows)  {
            // reduce block size for last block
            myBlockSize = nRows - r;
         }
         // --- apply rotation matrix to current block, output to workspace
// --- critical point for not having omp:
//     * how do we prevent parallelism inside zgemm???
//     * does it hurt if zgemm makes more threads??
//     * maybe if omp is used inside zgemm, omp can deal with it...
#ifdef USE_ATLAS
         cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols, 
                      one, mat + r, nRows, 
                      rotMat, nCols, zero, workspace, myBlockSize);
#elif defined USE_GOTO // khr: GotoBLAS, experimental!
         cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols,
                      one, mat + r, nRows,
                      const_cast<double*>(rotMat), nCols, zero, workspace,
                      myBlockSize);
#elif defined USE_ACML
         dgemm (noTrans,noTrans, myBlockSize, nCols, nCols, 
                one, (mat + r), nRows, 
                const_cast<double*>(rotMat), nCols, 
                zero, workspace, myBlockSize);
#elif defined USE_INTEL_MKL
         dgemm (&noTrans, &noTrans, &myBlockSize, &nCols, &nCols,
                &one, (mat + r), &nRows, rotMat, &nCols, &zero, 
                workspace, &myBlockSize);
#else
         SX_EXIT;
#endif
         // copy result into mat
         for (c = 0; c < nCols; ++c)
            memcpy (mat + r + nRows * c, workspace + myBlockSize * c,
                    sizeof(double) * myBlockSize);
      }
      delete [] workspace;
   }
}

void inPlaceRot (SxComplex8 *mat,  const SxComplex8 *rotMat, 
                 int nRows, int nCols)
{
   if (nRows <= 0 || nCols <=0) return;
   // make sure that mat and rotMat do not overlap
   SX_CHECK (   mat    >= rotMat + nCols * nCols
             || rotMat >= mat    + nRows * nCols,
             (mat-rotMat), nRows * nCols);
   int blockSize = 1024;
#if defined USE_ATLAS || defined USE_ACML || defined USE_INTEL_MKL
   const int minSize = 2 * blockSize;
   if (nRows < minSize || sizeof(SxComplex8) * nRows * nCols < 102400)
#endif
   {
      SxComplex8 *workspace = new SxComplex8[nRows * nCols];
      memcpy (workspace, mat, sizeof(SxComplex8) * nRows * nCols);
      matmult (mat, workspace, rotMat, nRows, nCols, nCols);
      delete [] workspace;
      return;
   }
   
   SX_EXIT; // needs to be tested against implementation above

   const SxComplex8 one (1.0, 0.0), zero (0.0, 0.0);
#if ! defined USE_ATLAS
   char noTrans = 'N';
#endif
// sketch for omp: #  pragma omp parallel
   {
      SxComplex8 *workspace;
      // catching the bad_alloc is essential for omp parallelism, because exceptions from omp
      // threads cannot be caught from the outside
      try {
         workspace = new SxComplex8[blockSize * nCols];
      } catch (std::bad_alloc &err) {
         cout << "Memory allocation failed in inPlaceRot." << endl;
         cout << err.what () << endl;
         SX_EXIT;
      }
      int r, c, myBlockSize = blockSize;
// sketch for omp: #     pragma omp for
      for (r = 0; r < nRows; r += blockSize)  {
         if (r + blockSize > nRows)  {
            // reduce block size for last block
            myBlockSize = nRows - r;
         }
         // --- apply rotation matrix to current block, output to workspace
// --- critical point for not having omp:
//     * how do we prevent parallelism inside zgemm???
//     * does it hurt if zgemm makes more threads??
//     * maybe if omp is used inside zgemm, omp can deal with it...
#ifdef USE_ATLAS
         cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols, 
                      &one, mat + r, nRows, 
                      rotMat, nCols, &zero, workspace, myBlockSize);
#elif defined USE_GOTO // khr: GotoBLAS, experimental!
         cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols,
                      (float*)&one, (float*)(mat + r), nRows,
                      (float*)rotMat, nCols, (float*)&zero,
                      (float*)workspace, myBlockSize);
#elif defined USE_ACML
         cgemm (noTrans,noTrans, myBlockSize, nCols, nCols, 
                (complex*)const_cast<SxComplex8*>(&one), 
                (complex*)(mat + r), nRows, 
                (complex*)const_cast<SxComplex8*>(rotMat), nCols, 
                (complex*)const_cast<SxComplex8*>(&zero),
                (complex*)workspace, myBlockSize);
#elif defined USE_INTEL_MKL
         cgemm (&noTrans, &noTrans, &myBlockSize, &nCols, &nCols,
                (MKL_Complex8 *)const_cast<SxComplex8*>(&one),
                (MKL_Complex8 *)const_cast<SxComplex8*>(mat + r), &nRows,
                (MKL_Complex8 *)const_cast<SxComplex8*>(rotMat), &nCols,
                (MKL_Complex8 *)const_cast<SxComplex8*>(&zero),
                (MKL_Complex8 *)workspace, &myBlockSize);
#else
         SX_EXIT;
#endif
         // copy result into mat
         for (c = 0; c < nCols; ++c)
            memcpy (mat + r + nRows * c, workspace + myBlockSize * c,
                    sizeof(SxComplex8) * myBlockSize);
      }
      delete [] workspace;
   }
}

void inPlaceRot (SxComplex16 *mat,  const SxComplex16 *rotMat, 
                 int nRows, int nCols)
{
   if (nRows <= 0 || nCols <=0) return;
   // make sure that mat and rotMat do not overlap
   SX_CHECK (   mat    >= rotMat + nCols * nCols
             || rotMat >= mat    + nRows * nCols,
             (mat-rotMat), nRows * nCols);
#ifdef USE_INTEL_MKL
   int blockSize = 128;
#else
   int blockSize = 1024;
#endif
   const int minSize = 2 * blockSize;
#if defined USE_ATLAS || defined USE_ACML || defined USE_INTEL_MKL
   if (nRows < minSize || sizeof(SxComplex16) * nRows * nCols < 1024000)
#endif
   {
      SxComplex16 *workspace = new SxComplex16[nRows * nCols];
      memcpy (workspace, mat, sizeof(SxComplex16) * nRows * nCols);
      matmult (mat, workspace, rotMat, nRows, nCols, nCols);
      delete [] workspace;
      return;
   }
   const SxComplex16 one (1.0, 0.0), zero (0.0, 0.0);
#if ! defined USE_ATLAS
   char noTrans = 'N';
#endif
#ifdef USE_INTEL_MKL
   if (rotMat[nCols-1].absSqr () < 1e-24 * rotMat[0].absSqr ())  {
      // check for upper triangular matrix
      bool upT = true;
      for (int ic = 0; ic < nCols -1; ic++)  {
         double upNorm = norm2(rotMat + ic * nCols, ic + 1);
         double lowNorm = norm2(rotMat + ic * nCols + ic + 1, nCols - ic -1);
         if (lowNorm >= 1e-24 * upNorm) { upT = false; break; }
      }
      if (upT)  {
         char side = 'R', uplo = 'U', diag = 'N';
         ztrmm (&side, &uplo, &noTrans, &diag, &nRows, &nCols, 
                (MKL_Complex16 *)const_cast<SxComplex16*>(&one),
                (MKL_Complex16 *)const_cast<SxComplex16*>(rotMat), &nCols,
                (MKL_Complex16 *)const_cast<SxComplex16*>(mat), &nRows);
         return;
      }
   }
#endif
#ifdef USE_OPENMP
#  pragma omp parallel
#endif
   {
      SxComplex16 *workspace;
      // catching the bad_alloc is essential for omp parallelism, because exceptions from omp
      // threads cannot be caught from the outside
      try {
         workspace = new SxComplex16[blockSize * nCols];
      } catch (std::bad_alloc &err) {
         cout << "Memory allocation failed in inPlaceRot." << endl;
         cout << err.what () << endl;
         SX_EXIT;
      }
      int r, c, myBlockSize = blockSize;
#ifdef USE_OPENMP
#     pragma omp for
#endif
      for (r = 0; r < nRows; r += blockSize)  {
         if (r + blockSize > nRows)  {
            // reduce block size for last block
            myBlockSize = nRows - r;
         }
         // --- apply rotation matrix to current block, output to workspace
// --- critical point for not having omp:
//     * how do we prevent parallelism inside zgemm???
//     * does it hurt if zgemm makes more threads??
//     * maybe if omp is used inside zgemm, omp can deal with it...
#ifdef USE_ATLAS
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols, 
                      &one, mat + r, nRows, 
                      rotMat, nCols, &zero, workspace, myBlockSize);
#elif defined USE_GOTO // khr: GotoBLAS, experimental!
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      myBlockSize, nCols, nCols,
                      (double*)&one, (double*)(mat + r), nRows,
                      (double*)rotMat, nCols, (double*)&zero,
                      (double*)workspace, myBlockSize);
#elif defined USE_ACML
         zgemm (noTrans,noTrans, myBlockSize, nCols, nCols, 
                (doublecomplex*)const_cast<SxComplex16*>(&one), 
                (doublecomplex*)(mat + r), nRows, 
                (doublecomplex*)const_cast<SxComplex16*>(rotMat), nCols, 
                (doublecomplex*)const_cast<SxComplex16*>(&zero),
                (doublecomplex*)workspace, myBlockSize);
#elif defined USE_INTEL_MKL
         zgemm (&noTrans, &noTrans, &myBlockSize, &nCols, &nCols,
                (MKL_Complex16 *)const_cast<SxComplex16*>(&one),
                (MKL_Complex16 *)const_cast<SxComplex16*>(mat + r), &nRows,
                (MKL_Complex16 *)const_cast<SxComplex16*>(rotMat), &nCols,
                (MKL_Complex16 *)const_cast<SxComplex16*>(&zero),
                (MKL_Complex16 *)workspace, &myBlockSize);
#else
         SX_EXIT;
#endif
         // copy result into mat
         for (c = 0; c < nCols; ++c)
            memcpy (mat + r + nRows * c, workspace + myBlockSize * c,
                    sizeof(SxComplex16) * myBlockSize);
      }
      delete [] workspace;
   }
}




//------------------------------------------------------------------------------
// Erf and Erfc
//------------------------------------------------------------------------------
double derf (double x)
{
   return x < 0.0 ? -gammp (0.5, x*x) : gammp (0.5, x*x);
}


double derfc (double x)
{
	return x < 0.0 ? 1.0+gammp (0.5,x*x) : gammq (0.5,x*x);
}

double gammp (double a, double x)
{
   SX_CHECK (x >= 0.0 && a > 0.0, x, a);
   double gamser, gammcf, gln;

   if (x < (a+1.0))  {
      gser (&gamser, a, x, &gln);
      return gamser;
   }  else  {
      gcf (&gammcf, a, x, &gln);
      return 1.0 - gammcf;
   }
}



double gammq (double a, double x)
{
   SX_CHECK (x >= 0.0 && a > 0.0, x, a);
   double gamser, gammcf, gln;

   if (x < (a+1.0))  {
      gser (&gamser, a, x, &gln);
      return 1.0 - gamser;
   }  else  {
      gcf (&gammcf, a, x, &gln);
      return gammcf;
   }
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double gammln (double xx)
{
   double x, y, tmp, ser;
   static double cof[6] = { 76.18009172947146,     -86.50532032941677,
                            24.01409824083091,      -1.231739572450155,
                             0.1208650973866179e-2, -0.5395239384953e-5 };
   int j;

   if (xx == 0.5) return log(SQRT_PI); // dominant use
   y    = x = xx;
   tmp  = x + 5.5;
   tmp -= (x+0.5) * log(tmp);
   ser  = 1.000000000190015;
   for (j=0; j <= 5; j++)  ser += cof[j] / ++y;
   return -tmp + log (2.5066282746310005 * ser/x);
}



void gser (double *gamser, double a, double x, double *gln)
{
   SX_CHECK (x >= 0.0, x);

   static int ITMAX = 100;
   static double EPS = 3.0e-16;

   int n;
   double sum, del, ap;

   *gln = gammln (a);
   if ( x == 0.0 )  {
      *gamser = 0.0;
      return;
   }  else  {
      ap = a;
      del = sum = 1.0 / a;
      for (n=1; n <= ITMAX; n++)  {
         ap++;
         del *= x / ap;
         sum += del;
         if (fabs (del) < fabs(sum)*EPS)  {
            *gamser = sum * exp (-x+a*log(x)-(*gln));
            return;
         }
      }
      SX_EXIT;
      return;
   }
}

double gcf_sum (double a, double x)
{
   static int    ITMAX = 100;
   static double EPS   = 1.0e-15;
   static double FPMIN = 1.0e-100;
   
   int i;

   double an, b, c, d, del, h;
   b    = x + 1.0 - a;
   c    = 1.0 / FPMIN;
   d    = 1.0 / b;
   h    = d;
   for (i=1; i <= ITMAX; i++)  {
      an  = -i * ( i - a );
      b  += 2.0;
      d   = an * d + b;
      if (fabs (d) < FPMIN) d = FPMIN;
      c   = b + an / c;
      if (fabs (c) < FPMIN) c = FPMIN;
      d   = 1.0 / d;
      del = d * c;
      h  *= del;
      if (fabs (del-1.0) < EPS) break;
   }

   SX_CHECK (i <= ITMAX, i, ITMAX); // a too large, ITMAX to small
   return h;
}

void gcf (double *gammcf, double a, double x, double *gln)
{
   *gln = gammln(a);
   *gammcf = exp (-x + a * log(x) - (*gln) ) * gcf_sum(a, x);
}

double erfcexp (double x)
{
   double x2 = x * x;
   if (x2 < 1.5) return exp(x2) * gammq (0.5, x2);
   return x / SQRT_PI * gcf_sum(0.5, x2);
}


double lineMinimization (double y, double xT, double yT, double dYdX,
                         double *curvature)
{
   SX_CHECK ( fabs(xT) > 1e-50, xT);
   // see ref1, fig4
   double xMin, curv;
   curv = (yT - (y + xT*dYdX)) / (xT*xT);
   xMin = (fabs(curv) > 1e-50) 
        ? -dYdX / (2. * curv)
        : 0.;

   if (curvature)  *curvature = curv;
   return xMin;
}


#ifdef WIN32
double cbrt (double x)
{
   return pow (x, 0.333333333333333333333333333333);
}
#endif


#ifndef HAVE_SINCOS
void sincos (double x, double *s, double *c)
{
   *s = sin(x); 
   *c = cos(x);
}

void sincosf (float x, float *s, float *c)
{
   *s = sinf(x);
   *c = cosf(x);
}
#endif

// work-around for gcc bug 33647
#ifndef HAVE_GFORTRAN_COPY_STRING
/* see also 
   http://www.nabble.com/Patch:-optimized-copy_string-t717416.html */
extern "C" {
   void _gfortran_copy_string (int destlen, char * dest,
                               int srclen, const char * src)
   {
      if (srclen >= destlen)
      {
         /* This will truncate if too long.  */
         memmove (dest, src, destlen);
      }
      else
      {
         memmove (dest, src, srclen);
         /* Pad with spaces.  */
         memset (&dest[srclen], ' ', destlen - srclen);
      }
      const int copy = srclen <= destlen ? srclen : destlen;
      const int fill = srclen <= destlen ? destlen - srclen : 0;
      memmove (dest, src, copy);
      memset (dest + srclen, ' ', fill);
   }
}
#endif /* GFORTRAN_HAVE_COPY_STRING */
