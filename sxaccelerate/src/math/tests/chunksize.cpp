#include <stdio.h>
#include <SxTime.h>
#include <SxTimer.h>
#include <SxUtil.h>
#include <SxMathLib.h>
#include <SxVector.h>
#include <SxList.h>
#include <SxFFT3d.h>
#include <SxMatrix.h>
#ifdef USE_OPENMP
#  include <omp.h>
#endif /* USE_OPENMP */

void DISABLE_SMP ()
{
#  ifdef USE_OPENMP   
      omp_set_num_threads(1);
      sxChunkSize = 10000000000;
#ifdef USE_FFTW
      fftw_plan_with_nthreads (1);
#endif
#  endif /* USE_OPENMP */
}


void ENABLE_SMP ()
{
#  ifdef USE_OPENMP   
     omp_set_num_threads(omp_get_num_procs());
     //omp_set_num_threads(2);
     sxChunkSize = 1;  // enable always SMP
#ifdef USE_FFTW
      fftw_plan_with_nthreads (omp_get_num_procs());
#endif
#  endif /* USE_OPENMP */
}

/*
double benchmark (const SxVector<Double> &a, const SxVector<Double> &b)
{
   double r = 0.;
   int nIt = 1000;
   int n = a.getSize();
   SxVector<Double> c(n);
   for (int it=0; it < nIt; ++it)  {


      double *src1Ptr  = a.elements;
      double *src2Ptr  = b.elements;
      double *destPtr  = c.elements;

      if (n > chunkSize)  {
//#        pragma omp parallel shared(destPtr, src1Ptr, src2Ptr) private (smpIdx)
//         {
//#           pragma omp for schedule(dynamic,chunk) nowait
#           pragma omp parallel for shared(destPtr,src1Ptr,src2Ptr) \
            private(smpIdx) schedule(static)
            for (smpIdx=0; smpIdx < n; smpIdx++)
               destPtr[smpIdx] = src1Ptr[smpIdx] + src2Ptr[smpIdx]; 
//       }
      } else {
//       for (int i=0; i < n; i++, ++src1Ptr, ++src2Ptr)  
//          *destPtr++ = *src1Ptr + *src2Ptr; 
         for (int i=0; i < n; i++) destPtr[i] = src1Ptr[i] + src2Ptr[i];
      }

      c = a + b;

      r += c.absSqr().sum();
   }
   return r;
}
*/

double benchmark (const SxVector<Double> &a, const SxVector<Double> &b)
{
   double r = 0.;
   double t1;
   double t2;
   ssize_t i;
   ssize_t nLoops = 32;
   
   //SxVector<Double> c;
   
   // --- configure to 1 sec
   do {
      nLoops = nLoops << 1;
      t1 = SxTime::getRealTime ();
      for (i=0; i < nLoops; ++i)  {
         //c = a + b;
         //c = a * b;
         //r += a.sum ();
         r += (a * b + 2.*a * sin(a)).absSqr().sum ();
      }
      t2 = SxTime::getRealTime ();
   } while (t2 - t1 < 0.5);
   sxprintf ("%12.2f  ", (float)r);
   t1 = (t2 - t1) / (double)nLoops;   
   return t1;
}

int main()
{

   size_t memSize = 2048 *size_t(1024*1024);
   size_t n;
   size_t nMaxVec = 11000000;
   //size_t nMaxMat = 10000000;
   //size_t nMaxFFT = 200L*200L*200L;
   double maxTime = 120.;
   //double maxTimeFFT = 2.;
#ifdef USE_FFTW
   //fftw_set_timelimit (50. * maxTimeFFT);
#endif
   initSPHInXMath ();
   SxTimer timer (6, /* realTime: */true);
   timer.setName (0, "seq vector");
   timer.setName (1, "smp vector");
   timer.setName (2, "seq matrix");
   timer.setName (3, "smp matrix");
   timer.setName (4, "seq FFT");
   timer.setName (5, "smp FFT");

   SxList<size_t> nElemVec, nElemMat, nElemFFT;
   SxList<double> tVec0, tVec1, tMat0, tMat1, tFFT0, tFFT1;

   //double r = 0.;

   DISABLE_SMP ();
   printf ("Testing sequential vector operations...\n");
   for (n=10; n < nMaxVec; n*=2)  {
      if (size_t(n) * sizeof(double) * size_t(4) > memSize) break;
      SxVector<Double> a(n), b(n);
      a.randomize (); b.randomize ();
      tVec0 << benchmark (a, b);
      nElemVec << n;
      printf ("%8ld  %12.6f\n", n, tVec0.last ());
      if (timer.getTime(0) > maxTime / 2.) break;
   }
/*
   printf ("Testing sequential matrix operations...\n");
   for (n=10; n < ::sqrt(nMaxMat); n*=2)  {
      SxMatrix<Complex16> m1(n,n), m2(n,n);
      m1.randomize(); m2.randomize ();
      timer.reset (2); timer.start (2); 
      m1 ^ m2;
      timer.stop (2);
      tMat0 << timer.getTime(2);
      nElemMat << n*n;
      printf ("%8ld  %12.6f\n", n*n, timer.getTime(2));
      if (timer.getTime(2) > maxTime / 8.) break;
   }

   printf ("Testing sequential FFT...\n");
   for (n=10; n < cbrt(nMaxFFT); n*=2)  {
      SxFFT3d fft (SxFFT::Both, n, n, n);
      SxVector<Complex16> meshR (fft.meshSize), meshW(fft.meshSize);
      meshR.randomize ();
      timer.reset (4); timer.start (4); 
      fft.fftForward (fft.meshSize, meshR.elements, meshW.elements);
      fft.fftReverse (fft.meshSize, meshW.elements, meshR.elements);
      timer.stop (4);
      tFFT0 << timer.getTime(4);
      nElemFFT << n*n*n;
      printf ("%8ld  %12.6f\n", n*n*n, timer.getTime(4));
      if (timer.getTime(4) > maxTime / 8.) break;
   }
*/
#  ifndef USE_OPENMP
      printf ("ERROR: Cannot execute SMP benchmark.\n");
      printf ("       Please recompile S/PHI/nX with --enable-openmp.\n");
      exit(1);
#  else

   ENABLE_SMP ();
   printf ("Testing SMP vector operations...\n");
   for (n=10; n < nMaxVec; n*=2)  {
      SxVector<Double> a(n), b(n);
      a.randomize (); b.randomize ();
      tVec1 << benchmark (a, b);
      printf ("%8ld  %12.6f\n", n, tVec1.last ());
      if (tVec1.getSize () >= tVec0.getSize ()) break;
   }
/*
   printf ("Testing SMP matrix operations...\n");
   for (n=10; n < sqrt(nMaxMat); n*=2)  {
      SxMatrix<Complex16> m1(n,n), m2(n,n);
      m1.randomize(); m2.randomize ();
      timer.reset (2); timer.start (2); 
      m1 ^ m2;
      timer.stop (2);
      tMat1 << timer.getTime(2);
      printf ("%8ld  %12.6f\n", n*n, timer.getTime(2));
      if (tMat1.getSize () >= tMat0.getSize ()) break;
   }

   printf ("Testing SMP FFT...\n");
   for (n=10; n < cbrt(nMaxFFT); n*=2)  {
      SxFFT3d fft (SxFFT::Both, n, n, n);
      SxVector<Complex16> meshR (fft.meshSize), meshW(fft.meshSize);
      meshR.randomize ();
      timer.reset (4); timer.start (4); 
      fft.fftForward (fft.meshSize, meshR.elements, meshW.elements);
      fft.fftReverse (fft.meshSize, meshW.elements, meshR.elements);
      timer.stop (4);
      tFFT1 << timer.getTime(4);
      printf ("%8ld  %12.6f\n", n*n*n, timer.getTime(4));
      if (tFFT1.getSize () >= tFFT0.getSize ()) break;
   }
*/

   // --- print statistics
   int nProcs = omp_get_num_procs();
   printf ("\nResults: vector operations\n");
   printf ("    size      serial[ms]     smp(%d)[ms]    serial/smp\n", nProcs);
   for (int i=0; i < tVec0.getSize(); ++i)  {
      printf ("%8ld  %12.6f   %12.6f   %12.6f\n", nElemVec(i), 
              tVec0(i) * 1e3,  tVec1(i) * 1e3, 
              tVec1(i) > 1e-10 ? tVec0(i) / tVec1(i) : 0.);
   }
/*
   printf ("\nResults: matrix operations\n");
   printf ("    size      serial         smp(%d)        serial/smp\n", nProcs);
   for (int i=0; i < tMat0.getSize(); ++i)  {
      printf ("%8ld  %12.6f   %12.6f   %12.6f\n", nElemMat(i), 
              tMat0(i),  tMat1(i), 
              tMat1(i) > 1e-10 ? tMat0(i) / tMat1(i) : 0.);
   }

   printf ("\nResults: FFT\n");
   printf ("    size      serial         smp(%d)        serial/smp\n", nProcs);
   for (int i=0; i < tFFT0.getSize(); ++i)  {
      printf ("%8ld  %12.6f   %12.6f   %12.6f\n", nElemFFT(i), 
              tFFT0(i),  tFFT1(i), 
              tFFT1(i) > 1e-10 ? tFFT0(i) / tFFT1(i): 0.);
   }
   */
#  endif /* USE_OPENMP */
      
   return 0;
}
