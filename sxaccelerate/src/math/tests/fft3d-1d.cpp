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

#include <SxUtil.h>
#include <SxFFT1d.h>
#include <SxFFT3d.h>
//#include <SxDirac.h>
#include <SxVector.h>
#include <SxMatrix3.h>
#include <SxConstants.h>  /* for complex I */
#include <SxComplex.h>
#include <SxTimer.h>
//#include <SxBinIO.h>
#include <stdio.h>
#include <iostream>

SxVector<SX_T_FFT_COMPLEX> multiple1d (      SxFFT1d &, 
                                               SxFFT1d &,
                                               SxFFT1d &,
                                         const SxFFT3d &,
                                         const SxVector<SX_T_FFT_COMPLEX> &);
/** 
  \example fft3d-1d.cpp

  This test validates the compatibility of the SxFFT1d and the SxFFT3d by 
  checking that the results of a single 3d-FFT call yields the same as
  multiple 1d-FFT calls.
  \author Sixten Boeck
  */
int main ()
{
   SxTimer timer (2);

   // --- initialize FFT library
   FFT_REAL omega;
   //SxVector3<Int> mesh (14, 16, 20);  omega = 100;
   SxVector3<Int> mesh (200, 50, 50); omega = 123.456;
   //SxVector3<Int> mesh (64, 64, 64); omega = 123.456;

   SxFFT3d fft3d (SxFFT3d::Forward, mesh(0), mesh(1), mesh(2), omega);

   double scaleFor = fft3d.scaleFor;
   double scaleRev = fft3d.scaleRev;
   double nrm      = 3;
   double dOmega   = omega / fft3d.meshSize;

   printf ("FFT parameters\n");
   printf ("     norm:     %g\n", nrm);
   printf ("     omega:    %g\n", omega);
   printf ("   1/omega:    %g\n", 1./omega);
   printf ("    dOmega:    %g\n", dOmega);
   printf ("     mesh:     %d x %d x %d (%d)\n", mesh(0), mesh(1), mesh(2),
                                                 mesh.product());
   printf ("     scaleFor: %g\n", scaleFor);
   printf ("     scaleRev: %g\n", scaleRev);
   printf ("   1/scaleFor: %g\n", 1./scaleFor);
   printf ("   1/scaleRev: %g\n", 1./scaleRev);
   printf ("\n");

   // --- only type SX_T_FFT_COMPLEX is allowed (defined at link time)!
   //     via 3d-FFT
   SxVector<SX_T_FFT_COMPLEX> meshW(fft3d.meshSize); // src in frequency space
   SxVector<SX_T_FFT_COMPLEX> meshR(fft3d.meshSize); // dest. in real space
   //     via multiple 1d-FFT
   SxVector<SX_T_FFT_COMPLEX> meshR1(fft3d.meshSize); // dest. in real space

   // --- performing multiple 1d FFTs (w -> R)
   //SxFFT1d fftX (SxFFT1d::Forward, mesh(0), 1.0);
   //SxFFT1d fftY (SxFFT1d::Forward, mesh(1), 1.0);
   //SxFFT1d fftZ (SxFFT1d::Forward, mesh(2), 1.0);
   SxFFT1d fftX (SxFFT1d::Forward, mesh(0), 1.0, mesh(1)*mesh(2));
   SxFFT1d fftY (SxFFT1d::Forward, mesh(1), 1.0, mesh(2) );
   SxFFT1d fftZ (SxFFT1d::Forward, mesh(2), 1.0);

   fftX.autoscale = fftY.autoscale = fftZ.autoscale = false;
   
   // --- Checking forward FFT  (w -> R)
   printf ("Checking forward FFT\n");
   meshW.randomize ();
   printf ("3d-fft\n"); fflush (stdout);
   for (int it=0; it < 10; it++)  {
      printf ("it=%d\n", it); fflush (stdout);
      timer.start(0);
      fft3d.fftForward (fft3d.meshSize, meshW.elements, meshR.elements);
      timer.stop(0);


 
      printf ("1d-ffts\n"); fflush (stdout);
      timer.start(1);
      meshR1 = multiple1d (fftX, fftY, fftZ, fft3d, meshW);
      timer.stop(1);
   }
   printf ("1d-ffts done\n"); fflush (stdout);
   SxMatrix3<Double> aMat (10.,0.,0., 0.,10.,0., 0.,0.,10.);
   /*
   try  {
      SxBinIO io ("meshR.sxb", SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (meshR, aMat, fft3d.mesh);
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (meshR, aMat, fft3d.mesh);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   try  {
      SxBinIO io ("meshR1.sxb", SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (meshR1, aMat, fft3d.mesh);
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (meshR1, aMat, fft3d.mesh);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   */

   printf (">>> %20.18f\n", (meshR1 - meshR).absSqr().sum() );
   timer.print ();

   return 0;
}


SxVector<SX_T_FFT_COMPLEX> multiple1d (      SxFFT1d &fftX,
                                               SxFFT1d &fftY,
                                               SxFFT1d &fftZ,
                                         const SxFFT3d &fft3d,
                                         const SxVector<SX_T_FFT_COMPLEX> &meshIn)
{
   SxTimer timer (10);
   printf ("1d-ffts preparation\n"); fflush (stdout);
   int x, y, z;
   SxVector3<Int> mesh (fft3d.mesh(0), fft3d.mesh(1), fft3d.mesh(2));

   SxComplex<FFT_REAL> sum;
   SxVector<SX_T_FFT_COMPLEX::TReal> arg;
   SxVector<SX_T_FFT_COMPLEX> meshOut (mesh.product()), meshOut1(mesh.product()), meshOut2(mesh.product());
   SxVector<SX_T_FFT_COMPLEX> lineInX, lineInY, lineInZ, lineOut;

   SxIdx range;
   int idx, idx0, stride;

   printf ("1d-ffts (I)\n"); fflush (stdout);
   // --- 1d FFTs along z axis
   timer.start(0);
   lineInZ = SxVector<SX_T_FFT_COMPLEX> (mesh(2), 0.);
   lineOut = SxVector<SX_T_FFT_COMPLEX> (mesh(2));
   for (x=0; x < mesh(0); x++)  {
      for (y=0; y < mesh(1); y++)  {
         
         idx     = fft3d.mesh.getMeshIdx(x,y,0, SxMesh3D::Unknown);
         //range   = SxIdx (idx, idx + mesh(2)-1);
         //lineInZ = meshIn(range);
         //fftZ.fftForward (mesh(2), lineInZ.elements, meshOut.elements + idx);
         fftZ.fftForward (mesh(2), meshIn.elements+idx, meshOut.elements + idx);
      }
   }
   timer.stop(0);

   printf ("1d-ffts (II)\n"); fflush (stdout);
   // --- 1d FFTs along y axis
   timer.start(1);
   lineInY = SxVector<SX_T_FFT_COMPLEX> (mesh(1), 0.);
   lineOut = SxVector<SX_T_FFT_COMPLEX> (mesh(1));
   stride  = mesh(2);
   for (x=0; x < mesh(0); x++)  {
      for (z=0; z < mesh(2); z++)  {
         
         idx = idx0 = fft3d.mesh.getMeshIdx (x,0,z, SxMesh3D::Unknown);
         //for (y=0; y < mesh(1); y++, idx+=stride)  lineInY(y) = meshOut(idx);
         //fftY.fftForward (mesh(1), lineInY.elements, lineOut.elements);
         fftY.fftForward (mesh(1), meshOut.elements+idx, meshOut1.elements+idx);

         //idx = idx0;
         //for (y=0; y < mesh(1); y++, idx+=stride)  meshOut1(idx) = lineOut(y);
      }
   }
   timer.stop (1);

   printf ("1d-ffts (III)\n"); fflush (stdout);
   timer.start (2);
   // --- 1d FFTs along x axis
   lineInX = SxVector<SX_T_FFT_COMPLEX> (mesh(0), 0.);
   lineOut = SxVector<SX_T_FFT_COMPLEX> (mesh(0));
   stride  = mesh(1) * mesh(2);
   idx = 0;
   for (y=0; y < mesh(1); y++)  {
      for (z=0; z < mesh(2); z++)  {

         idx = idx0 = fft3d.mesh.getMeshIdx(0,y,z, SxMesh3D::Unknown);
         //for (x=0; x < mesh(0); x++, idx+=stride)  lineInX(x) = meshOut1(idx);
         //fftX.fftForward (mesh(0), lineInX.elements, lineOut.elements);
         fftX.fftForward (mesh(0), meshOut1.elements+idx, meshOut2.elements+idx);

         //idx = idx0;
         //for (x=0; x < mesh(0); x++, idx+=stride)  meshOut2(idx) = lineOut(x);
      }
   }
   timer.stop (2);

   meshOut2 *= fft3d.scaleFor;

   timer.print ();
   return meshOut2;
}

