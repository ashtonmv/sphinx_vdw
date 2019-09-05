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
#include <SxFFT3d.h>
#include <SxVector.h>
#include <SxConstants.h>  /* for complex I */
#include <SxComplex.h>
#include <SxCLI.h>
#include <stdio.h>
#include <iostream>
#include <SxTimer.h>

enum LocalTimer { TimeFFT};

SX_REGISTER_TIMERS(LocalTimer)
{
   regTimer (TimeFFT, "FFT");
}


SxVector<SX_T_FFT_COMPLEX> discreteFT (const SxFFT3d &,
                                       const SxVector<SX_T_FFT_COMPLEX> &,
                                       FFT_REAL);
/** 
  \example fft3d.cpp

  This small routine is not supposed to be an example how to call FFT.
  It is rather a tool to port the SxFFT3d class to new platforms (tuning
  of prefactors etc.).
  Note, you should never ever call FFT routines directly. 
  Memory handling and data type setting is tricky!
  Use FFT interfaces of SxGBasis or SxGkBasis instead!!!
  \author Sixten Boeck
  */
int main (int argc, char **argv)
{
   {
      SxCLI cli(argc, argv);
      SxVector3<Int> mesh (20, 12, 18);
      if (cli.option ("--mesh", "mesh", "FFT mesh size").exists ())  {
         mesh = SxVector3<Int> (cli.last ().toIntList3 ());
      }
      cli.last ().optional = true;
      cli.last ().defaultValue = "default: 20 x 12 x 18";
      SxString fftWisdom;
      if (SxFFT::hasWisdom ())
         fftWisdom = cli.option("--wisdom","file","FFT wisdom file")
                     .toString ("");
      cli.finalize ();

      // load FFT wisdom to speed up FFT setup
      if (fftWisdom.getSize () > 0) SxFFT::loadWisdom (fftWisdom);
   
      // --- initialize vector class, setting up random number generator etc.
      initSPHInXMath ();

      // --- initialize FFT library
      FFT_REAL omega;
      //SxVector3<Int> mesh (10, 10, 10);  omega = 100;
      omega = 123.456;
      //SxVector3<Int> mesh (64, 64, 64); omega = 123.456;

      SxFFT3d fft (SxFFT3d::Both, mesh(0), mesh(1), mesh(2), omega);

      double scaleFor = fft.scaleFor;
      double scaleRev = fft.scaleRev;
      double nrm      = 3;
      double dOmega   = omega / fft.meshSize;

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
      SxVector<SX_T_FFT_COMPLEX> meshW(fft.meshSize); // src in frequency space
      SxVector<SX_T_FFT_COMPLEX> meshR(fft.meshSize); // dest. in real space

      
      // --- Checking forward FFT  (w -> R)
      printf ("Checking forward FFT\n");
      printf ("   Checking scaleForFFT factor\n");
      meshW.set (0.); meshW(0) = 1.0;
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      cout << "   sum meshR  = " << meshR.sum()
           << " (should be " << fft.meshSize / sqrt(omega) << ") " << endl;
      cout << endl;

      // --- Checking scaleFor factor
      printf ("   Checking scaleFor factor\n");
      meshW.set (0.); meshW(0) = nrm/omega;
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      meshR /= scaleFor;
      cout << "   sum meshR  = " << meshR.sum() * dOmega
           << " (should be " << nrm << ") " << endl;
      cout << endl;

      // --- Checking reverse FFT (R -> w)
      printf ("Checking reverse FFT\n");
      printf ("   Checking scaleRev factor\n");
      meshR.set (nrm / omega / fft.meshSize ); 
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      meshW /= scaleRev;
      cout << "   sum meshW(0)  = " << meshW(0)
           << " (should be "        << nrm/omega << ")\n";
      cout << endl;

      // --- Checking scaleRev factor
      printf ("   Checking scaleRevFFT factor\n");
      meshR.set (0.); meshR(0) = 1.;
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      cout << "   sum meshW  = " << meshW.sum()
           << " (should be " << sqrt(omega) << ")\n";
      cout << endl;


      // --- Checking reverse[forward] == forward[reverse]
      printf ("Checking identity\n");
      SxVector<SX_T_FFT_COMPLEX> meshT(fft.meshSize);
      meshT.randomize();
      meshW.copy (meshT);
      SX_START_TIMER (TimeFFT);
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      SX_STOP_TIMER (TimeFFT);
      cout << "   error = " << (meshW - meshT).absSqr().sum() 
           << " (should be 0)\n";

      meshR.copy (meshT);
      SX_START_TIMER (TimeFFT);
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      SX_STOP_TIMER (TimeFFT);
      cout << "   error = " << (meshR - meshT).absSqr().sum()
           << " (should be 0)\n";
   }

   printTiming ();

   return 0;
}


/** Computes a discrete Fourier transformation.
    \param meshIn  input mesh
    \param fftSign Either +1 or -1. This is the sign in the expontent
                   \f$ e^{+2\pi/n} \f$ or \f$ e^{-2\pi/n}\f$.
    \return transformed mesh  */
SxVector<SX_T_FFT_COMPLEX> discreteFT (const SxFFT3d &fft,
                                       const SxVector<SX_T_FFT_COMPLEX> &meshIn,
                                       FFT_REAL fftSign)
{
   int x, y, z, i, j, k;
   SxVector3<Int> mesh (fft.mesh(0), fft.mesh(1), fft.mesh(2));

   SxComplex<FFT_REAL> sum;
   SxVector3<SX_T_FFT_COMPLEX::TReal> arg;
   SxVector<SX_T_FFT_COMPLEX> meshOut (mesh.product());
   SxComplex<FFT_REAL> I_fftSign (fftSign, 1.);

   for (x=0; x<mesh(0); x++)  {
      printf ("               \r");
      printf ("%d of %d", x, fft.mesh(0)); fflush (stdout);
      for (y=0; y<mesh(1); y++)  {
         for (z=0; z<mesh(2); z++)  {

            sum = 0.;
            for (i=0; i<mesh(0); i++)  {
               for (j=0; j<mesh(1); j++)  {
                  for (k=0; k<mesh(2); k++)  {
                     arg(0)  = TWO_PI/fft.mesh(0) * i*x;
                     arg(1)  = TWO_PI/fft.mesh(1) * j*y;
                     arg(2)  = TWO_PI/fft.mesh(2) * k*z;
                     sum += meshIn(fft.mesh.getMeshIdx(i,j,k, SxMesh3D::Positive)) 
                          * (cos(arg(0)) + I_fftSign*sin(arg(0)))
                          * (cos(arg(1)) + I_fftSign*sin(arg(1)))
                          * (cos(arg(2)) + I_fftSign*sin(arg(2)));
                  }
               }
            }
            
            meshOut(fft.mesh.getMeshIdx(x,y,z, SxMesh3D::Positive)) = sum;

         }
      }
   }
   printf ("               \r");
   return meshOut;
}

