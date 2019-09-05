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


#ifndef _SX_FFT_3D_H_
#define _SX_FFT_3D_H_

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


/** \brief 3-dimensional FFT

    \b SxFFT3d = SFHIngX Fast Fourier Transformation in 3 dimensions

    This is the implementation of a cross-platform 3D Fast Fourier 
    Transformation.

    \ingroup group_num
    \author  Sixten Boeck
  */
class SX_EXPORT_MATH SxFFT3d : public SxFFT
{
   public:

      /// Constructor
      SxFFT3d ();
      /** \brief Constructor
          \param dir       direction of the FFT \sa Directions
          \param nx,ny,nz  dimensions of the FFT
          \param omega     real space unit cell volume
          \param symmetric FFTW only: if set to false, ignore volume scaling
                           factor.
       */
      SxFFT3d (Directions dir, int nx, int ny, int nz, 
               double omega=1.0, bool symmetric=true);
      /** \brief Constructor
          \param dir       direction of the FFT \sa Directions
          \param meshSize  dimensions of the FFT
          \param omega     real space unit cell volume
          \param symmetric FFTW only: if set to false, ignore volume scaling
                           factor.
        */
      SxFFT3d (Directions dir, const SxVector3<Int> meshSize, 
               double omega=1.0, bool symmetric=true);
      SxFFT3d (const SxFFT3d &);
      ~SxFFT3d ();
      
      /** \brief Set the mesh and initialize the FFT according to #dir
        \param nx
        \param ny
        \param nz the mesh dimensions
        \param omega the cell volume
        */
      void setMesh (int nx, int ny, int nz, double omega);

      /** \brief Set the mesh and initialize the FFT according to #dir
        \param mesh the mesh dimensions
        \param omega the cell volume
        */
      void setMesh (const SxVector3<Int> &mesh_, double omega)
      {
         setMesh (mesh_(0), mesh_(1), mesh_(2), omega);
      }
      
      /** \brief reverse Fourier transformation (K -> R)

          The reverse Fourier transformation is usually the transformation
          of a mesh given in frequency space to the realspace.
       
          \f[
              Y(\mathbf{r})
                 = 
              \sum_{j=0}^{n-1} 
                 X_j(\mathbf{k}) e^{+i \mathbf{k}\cdot\mathbf{r}}
          \f]  */
      void fftReverse (int n, const void *in, void *out);
      void fftReverse ()  {
         fftReverse (meshSize, inArray, outArray);
      }

      /** \brief forward Fourier transformation (R -> K)

          The forward Fourier transformation is usually the transformation
          of a mesh given in realspace to the frequency space.
       
           \f[
              Y(\mathbf{k}) 
                 = 
              \sum_{j=0}^{n-1} 
                 X_j(\mathbf{r}) e^{-i \mathbf{k}\cdot\mathbf{r}}
          \f]  */
      void fftForward (const int n, const void *in, void *out);
      inline void fftForward () {
         fftForward (meshSize, inArray, outArray);
      }

      void print ();


      /** Physical dimension of the FFT mesh */
      SxMesh3D mesh;
#     if defined (USE_FFTSG)
         /** The optimal mesh size. Usually it is equal the physical mesh size 
             ::mesh. On special machines (RISC or vector machines) the size 
             can differ to increase performance. It will be setup in the 
             constructor depending on the used FFT library. */
         SxVector3<Int> perfMesh;
         int fftCacheSize;  // TODO: should go to include/SxConfig.h
#     endif /* USE_FFTSG */      

      void clean () { SxFFT::clean (inArray); }

      SxFFT3d &operator= (const SxFFT3d&);
     
      void print (FFT_COMPLEX *);

      /** \brief Switch to trigger verbosity of debug messages in the 3D FFTW
       *         routines */
      static const bool verbdebug=false;


#if defined (USE_FFTW)
#  if defined (USE_FFTW_PARALLEL)

         /**  \brief Temporary array #1 needed for the MPI parallel transform */
         FFT_COMPLEX *fftwTmp1;
         /**  \brief Temporary array #2 needed for the MPI parallel transform */
         FFT_COMPLEX *fftwTmp2;


         /** \brief decomposition in x:
          *         array which holds the number of points
          *         to work on for each MPI rank */
         int *nxMPI;
         /** \brief decomposition in x:
          *         array which holds the start index
          *         to work from for each MPI rank */
         int *x0MPI;
         /** \brief decomposition in x:
          *         array which holds the stop index
          *         to work up to for each MPI rank */
         int *x1MPI;

         /** \brief decomposition in KY x KZ:
          *         array which holds the number of points
          *         to work on for each MPI rank */
         int *nKyKzMPI;
         /** \brief decomposition in KY x KZ:
          *         array which holds the start index
          *         to work from for each MPI rank */
         int *KyKz0MPI;
         /** \brief decomposition in KY x KZ:
          *         array which holds the stop index
          *         to work up to for each MPI rank */
         int *KyKz1MPI;


         /** \brief allocate memory for temporary arrays */
         void createMPIArrays(int nx, int ny, int nz);
         /** \brief free memory of temporary arrays */
         void deleteMPIArrays();
         /** \brief fill integer arrays with information in the domain decomposition */
         void initMPIArrays(int nx, int ny, int nz);

         /** \brief actually execute the parallel FFT */
         void doParallelFFTW(FFT_COMPLEX *input,
                             FFT_COMPLEX *output,
                             int dir);

#  endif // USE_FFTW_PARALLEL
#endif // USE_FFTW

};

#endif /* _SX_FFT_3D_H_ */
