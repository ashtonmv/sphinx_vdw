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

#ifndef _SX_FFT_2D_H_
#define _SX_FFT_2D_H_

#include <SxPrecision.h>
#include <SxMath.h>
#include <SxError.h>
#include <SxComplex.h>
#include <SxList.h>
#include <SxString.h>
#include <SxConfig.h>
#include <SxFFT.h>
#include <SxVector3.h>


/** \brief 1-dimensional FFT

    \b SxFFT2d = SFHIngX Fast Fourier Transformation in 2 dimensions

    This is the implementation of a cross-platform 2D Fast Fourier
    Transformation.

    \ingroup group_num
    \author  Sixten Boeck
  */
class SX_EXPORT_MATH SxFFT2d : public SxFFT
{
   public:

      /// Constructor
      SxFFT2d ();
      /// Constructor
      SxFFT2d (Directions dir, int nx, int ny, double omega=1.0,
               bool symmetric=true);
      /// Copy constructor
      SxFFT2d (const SxFFT2d &);
      /// Destructor
      ~SxFFT2d ();

      /// Copy assignment
      SxFFT2d &operator= (const SxFFT2d &);
      
      void setMesh (int nx, int ny, double omega);

      void fftReverse (int n, const void *in, void *out);
      void fftReverse ()  {
         fftReverse (meshSize, inArray, outArray);
      }

      void fftForward (const int n, const void *in, void *out);
      inline void fftForward () {
         fftForward (meshSize, inArray, outArray);
      }

      /// Turn (x,y) into linear mesh index 
      inline int getFFTIdx (int x, int y)  const {
         if (x < 0) x += mesh[0];
         if (y < 0) y += mesh[1];
         SX_CHECK (x>=0 && x < mesh[0], x, mesh[0]);
         SX_CHECK (y>=0 && y < mesh[1], y, mesh[1]);

         return y + mesh[1] * x;
      }

      /** \brief Turn linear mesh index into vector

        In lack of a SxVector2<Int>, we use SxVector3<Int> and
        set the 3rd value to 0.

        This functions wraps nx/2 < ix < nx to -nx/2<ix<0.
        This functions wraps ny/2 < iy < ny to -ny/2<iy<0.
        Unwrapped indices can be obtained via ix = idx / ny, iy = idx % ny.
      */
      inline SxVector3<Int> getMeshVecOrigin(ssize_t idx) const
      {
         SX_CHECK (meshSize > 0);
         SX_CHECK(idx >= 0 && idx < meshSize, idx, meshSize);
         SxVector3<Int> res;
         res(2) = 0;
         res(1) = int(idx/mesh[1]);
         res(0) = int(idx - res(1) * mesh[1]);
         if (2 * res(0) > mesh[0]) res(0) -= mesh[0];
         if (2 * res(1) > mesh[1]) res(1) -= mesh[1];
         return res;
      }

      void print ();


      /** Physical dimension of the FFT mesh */
      int mesh[2];
#     if defined (USE_FFTSG)
         /** The optimal mesh size. Usually it is equal the physical mesh size 
             ::mesh. On special machines (RISC or vector machines) the size 
             can differ to increase performance. It will be setup in the 
             constructor depending on the used FFT library. */
         int perfMesh[2];
         int fftCacheSize;  // TODO: should go to include/SxConfig.h
#     endif /* USE_FFTSG */      

      void clean () { SxFFT::clean (inArray); }
     
      void print (FFT_COMPLEX *);

};



#endif /* _SX_FFT_2D_H_ */
