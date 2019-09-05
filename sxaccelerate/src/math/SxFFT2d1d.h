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


#ifndef _SX_FFT_2D1d_H_
#define _SX_FFT_1D1D_H_

#include <SxPrecision.h>
#include <SxMath.h>
#include <SxError.h>
#include <SxMesh3D.h>
#include <SxComplex.h>
#include <SxList.h>
#include <SxVector.h>
#include <SxString.h>
#include <SxConfig.h>
#include <SxFFT.h>

#if defined(USE_FFTW) && defined(HAVE_POSIX_MEMALIGN) && defined(__AVX__)
// note: MKL gets very slow in innerloop 1D FFT when run in parallel with openMP
#  if ! defined(USE_MKL_FFT) || ! defined(USE_OPENMP)
#     define USE_FFT2d1d
#  endif
#else
#  if defined(USE_FFT2d1d)
#    error "SxFFT2d1d needs FFTW API (USE_FFTW) and posix_memalign and AVX"
#  endif
#endif

#ifndef USE_FFTW
typedef void* fftw_plan;
#endif

/** \brief 2+1-dimensional FFT

    \b SxFFT3d = SFHIngX Fast Fourier Transformation in 2+1 dimensions

    Some algorithm that use FFTs perform rather trivial transformations on the
    output of the FFT. For instance, the convolution in the application of the
    DFT real-space potential is FFT + multiply with potential + reverse FFT.

    In this case, it is advantageous to intertwine the trivial
    part with the FFTs to make use of data locality. This is achieved by doing
    transformations on two dimensions first for all data, and then combine the
    remaining 1-dimensional transforms with the operation in question.

    Moreover, by splitting the Fourier transforms in two steps, we can make use
    of the fact that many Fourier coefficients are zero/ignored. We can therefore
    avoid the 2D transforms on planes that contain no non-zero coefficients.

    In practice,


    \ingroup group_num
    \author  Christoph Freysoldt, freysoldt@mpie.de
  */
class SX_EXPORT_MATH SxFFT2d1d
{
   public:

      /** \brief Constructor
          \param nx,ny,nz  dimensions of the FFT
          \param nxNonZero nonZero elements go from -nxNonZero ... nxNonZero - 1
       */
      SxFFT2d1d (int nx, int ny, int nz, int nxNonZero);
      /// Destructor
      ~SxFFT2d1d ();

      /// Copy constructor (not implemented)
      SxFFT2d1d (const SxFFT2d1d &); // not implemented
      /// copy assignment (not implemented)
      void operator= (const SxFFT2d1d &); // not implemented

      /** Physical dimension of the FFT mesh */
      SxMesh3D realMesh;

   protected:
      /// Number of elements in the Nx x Ny plane, with padding
      ssize_t N12p;

      /// Nonzero elements on Nz dimension (from -nonZero ... nonZero-1)
      int nonZero;

      /// Data for cleaning mesh where it is not used by G+k-basis.
      SxArray<ssize_t> cleanCode;

      /// Openmp entry points
      SxArray<ssize_t> cleanEntry;

#ifdef USE_FFTW
      fftw_plan plan1F, plan1B, plan2F, plan2B;
#endif

      /// The internal mesh (N12p x 2*nonZero)
      SxComplex16 *meshData;

      /// Utility: plan a 1D or 2D transform
      fftw_plan plan (int nDim, int *dim, int dir);

      /// Destructor utility: unregister from FFT plans
      void destroyPlan (fftw_plan plan);

      /// Allocate mesh
      void allocate ();

   public:
      /// Get the mesh data pointer
      inline SxComplex16* getMesh () {
         if (!meshData) allocate ();
         return meshData;
      }

      /// Free the mesh
      inline void freeMesh () {
         if (meshData) free(meshData);
         meshData = NULL;
      }

      /// Clean the mesh data
      void clean ();

      /// Convolute the mesh data with the real-space potential V
      void convolute (const double *V);

      /** \brief Transform to real space, take absSqr, and add to rho
        @param weight prefactor for this absSqr
        @param rho the density
        */
      void addToRho (double weight, double *rho);

      /** \brief Transform to real space */
      void fftForward (SxComplex16 *out);

      /** \brief Transform from real space */
      void fftBackward (const SxComplex16 *in);

      /// Translate the standard FFT index to the compact mesh index
      SxVector<TPrecFFTIdx> getN231 (const SxVector<TPrecFFTIdx> &n123);
};

#endif /* _SX_FFT_2D1D_H_ */
