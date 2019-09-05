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
#include <SxMathLib.h>
#include <SxVector.h>
#include <SxConstants.h>  /* for complex I */
#include <SxComplex.h>
#include <stdio.h>
#include <iostream>

SxVector<SX_T_FFT_COMPLEX> discreteFT (const SxFFT1d &,
                                       const SxVector<SX_T_FFT_COMPLEX> &,
                                       FFT_REAL);

/** 
  \example fft1d.cpp

  This small routine is not supposed to be an example how to call FFT.
  It is rather a tool to port the SxFFT1d class to new platforms (tuning
  of prefactors etc.).
  Note, you should never ever call FFT routines directly. 
  Memory handling and data type setting is tricky!
  Use FFT interfaces of SxGBasis or SxGkBasis instead!!!
  \author Sixten Boeck
  */
int main ()
{
   {
      // --- initialize vector class, setting up random number generator etc.
      initSPHInXMath ();

      // --- initialize FFT library
      FFT_REAL omega;
      //int mesh = 10;  omega = 100;
      int mesh = 20; omega = 123.456;

      SxFFT1d fft ;
      fft =  SxFFT1d (SxFFT1d::Both, mesh, omega);

      double scaleFor = fft.scaleFor;
      double scaleRev = fft.scaleRev;
      double nrm      = 3;
      double dOmega   = omega / fft.meshSize;

      printf ("FFT parameters\n");
      printf ("            norm: %g\n", nrm);
      printf ("           omega: %g\n", omega);
      printf ("         1/omega: %g\n", 1./omega);
      printf ("          dOmega: %g\n", dOmega);
      printf ("            mesh: %d\n", mesh);
      printf ("mesh/sqrt(omega): %g\n", mesh/sqrt(omega));
      printf ("        scaleFor: %g\n", scaleFor);
      printf ("        scaleRev: %g\n", scaleRev);
      printf ("      1/scaleFor: %g\n", 1./scaleFor);
      printf ("      1/scaleRev: %g\n", 1./scaleRev);
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
      meshR.set(-1.);
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      cout << "   error = " << (meshW - meshT).absSqr().sum() 
           << " (should be 0)\n";

      meshR.copy (meshT);
      fft.fftReverse (fft.meshSize, meshR.elements, meshW.elements);
      fft.fftForward (fft.meshSize, meshW.elements, meshR.elements);
      cout << "   error = " << (meshR - meshT).absSqr().sum()
           << " (should be 0)\n";
   }
   return 0;
}


/** Computes a discrete Fourier transformation.
    \param meshIn  input mesh
    \param fftSign Either +1 or -1. This is the sign in the expontent
                   \f$ e^{+2\pi/n} \f$ or \f$ e^{-2\pi/n}\f$.
    \return transformed mesh  */
SxVector<SX_T_FFT_COMPLEX> discreteFT (const SxFFT1d &fft,
                                       const SxVector<SX_T_FFT_COMPLEX> &meshIn,
                                       FFT_REAL fftSign)
{
   int x, i;
   int mesh = fft.mesh;

   SxComplex<FFT_REAL> sum;
   SX_T_FFT_COMPLEX::Real arg;
   SxVector<SX_T_FFT_COMPLEX> meshOut (mesh);
   SxComplex<FFT_REAL> I_fftSign (fftSign, 1.);

   for (x=0; x < mesh; x++)  {
      printf ("               \r");
      printf ("%d of %d", x, fft.mesh); fflush (stdout);
      sum = 0.;
      for (i=0; i<mesh; i++)  {
         arg  = TWO_PI/fft.mesh * i*x;
         sum += meshIn(fft.getFFTIdx(i)) * (cos(arg) + I_fftSign*sin(arg));
      }
      meshOut(fft.getFFTIdx(x)) = sum;
   }

   printf ("               \r");
   return meshOut;
}

