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

#ifndef _SX_FFT_1D_H_
#define _SX_FFT_1D_H_

#include <SxPrecision.h>
#include <SxMath.h>
#include <SxError.h>
#include <SxComplex.h>
#include <SxList.h>
#include <SxString.h>
#include <SxConfig.h>
#include <SxFFT.h>


/** \brief 1-dimensional FFT

    \b SxFFT1d = SFHIngX Fast Fourier Transformation in 1 dimension

    This is the implementation of a cross-platform 1D Fast Fourier
    Transformation.

    \ingroup group_num
    \author  Sixten Boeck
  */
class SX_EXPORT_MATH SxFFT1d : public SxFFT
{
   public:

      SxFFT1d ();
      SxFFT1d (Directions dir, int nx, double omega=1.0, int strideIn=1,
               bool symmetric=true);
      SxFFT1d (const SxFFT1d &);
      ~SxFFT1d ();

      SxFFT1d &operator= (const SxFFT1d &);
      
      void setMesh (int nx, double omega, int stride);

      void init (int nx, double omega=1.0);

      void fftReverse (int n, const void *in, void *out);
      void fftReverse ()  {
         fftReverse (meshSize, inArray, outArray);
      }

      void fftForward (const int n, const void *in, void *out);
      inline void fftForward () {
         fftForward (meshSize, inArray, outArray);
      }

      inline int getFFTIdx (int x)  const {
         SX_CHECK (x>=0 && x < mesh, x, mesh);

         return x;
      }

      void print ();


      /** Physical dimension of the FFT mesh */
      int mesh;
#     if defined (USE_FFTSG)
         /** The optimal mesh size. Usually it is equal the physical mesh size 
             ::mesh. On special machines (RISC or vector machines) the size 
             can differ to increase performance. It will be setup in the 
             constructor depending on the used FFT library. */
         int perfMesh;
         int fftCacheSize;  // TODO: should go to include/SxConfig.h
#     endif /* USE_FFTSG */      

      static bool supportsStridedFFT () {
#        if defined (USE_ESSL)
            SX_EXIT; // not yet implemented
            return false;
#        elif defined (USE_ACML_FFT)
            SX_EXIT; // not yet implemented
            return false;
#        elif defined (USE_VECLIB_FFT)
            SX_EXIT; // not yet implemented
            return false;
#        elif defined (USE_FFTSG)
            SX_EXIT; // not yet implemented
            return false;
#        else /* USE_FFTW */
            return true;
#        endif     
      }       

      void clean () { SxFFT::clean (inArray); }
     
      void print (FFT_COMPLEX *);

// protected:
      int stride;
};



#endif /* _SX_FFT_1D_H_ */
