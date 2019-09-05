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

#include <SxFFT.h>
#include <SxRandom.h>
#include <SxLoopMPI.h>
#include <SxCLI.h>


SxList<SxString> fftPlanSpecs;
SxList<void *>   fftPlans;
SxList<int>      fftPlanCounter;
#ifdef USE_FFTW
#  if defined (USE_FFTW_PARALLEL)
      SxList<SxString>    fftPlanSpecsParA;
      SxList<SxString>    fftPlanSpecsParB;
      SxList<void *>   fftPlansParA;
      SxList<void *>   fftPlansParB;
      SxList<int>         fftPlanCounterParA;
      SxList<int>         fftPlanCounterParB;
#  endif // USE_FFTW_PARALLEL
   unsigned SxFFT::fftPlanMode(FFTW_MEASURE); // default planning mode
#endif


SxFFT::SxFFT (Directions dirIn, bool autoscaleIn, bool symmetricIn) 
   : meshSize (0),
     inArray(NULL),
     outArray(NULL),
     scaleFor (1.),
     scaleRev (1.),
#    if defined (USE_FFTW)
     iPlanFor (-1),
     iPlanRev (-1),
#    if defined (USE_FFTW_PARALLEL)
        iPlanForParA (-1),
        iPlanRevParA (-1),
        iPlanForParB (-1),
        iPlanRevParB (-1),
#    endif // USE_FFTW_PARALLEL
#    elif defined (USE_ACML_FFT)
     iPlanFor  (-1),
     iPlanRev (-1),
#    else
     FFT_FORWARD (0),
     FFT_REVERSE (0),     
#    endif     
     dir (dirIn),
     autoscale (autoscaleIn),
     symmetric (symmetricIn),
     scaleForFFT (1.),
     scaleRevFFT (1.)
{
   // empty
}


SxFFT::SxFFT (const SxFFT &in)
   : meshSize (in.meshSize),
     inArray(NULL),
     outArray(NULL),
     scaleFor (in.scaleFor),
     scaleRev (in.scaleRev),
#    if defined (USE_FFTW)
     plan_fwd (in.plan_fwd),
     plan_rev (in.plan_rev),
     iPlanFor (in.iPlanFor),
     iPlanRev (in.iPlanRev),
#    if defined (USE_FFTW_PARALLEL)
        plan_fwd_par_a (in.plan_fwd_par_a),
        plan_rev_par_a (in.plan_rev_par_a),
        iPlanForParA   (in.iPlanForParA),
        iPlanRevParA   (in.iPlanRevParA),
        plan_fwd_par_b (in.plan_fwd_par_b),
        plan_rev_par_b (in.plan_rev_par_b),
        iPlanForParB   (in.iPlanForParB),
        iPlanRevParB   (in.iPlanRevParB),
#    endif // USE_FFTW_PARALLEL
#    elif defined (USE_ACML_FFT)
     plan_fwd (in.plan_fwd),
     plan_rev (in.plan_rev),
     iPlanFor (in.iPlanFor),
     iPlanRev (in.iPlanRev),
#    else
     FFT_FORWARD (in.FFT_FORWARD),
     FFT_REVERSE (in.FFT_REVERSE),     
#    endif     
     dir (in.dir),
     autoscale (in.autoscale),
     symmetric (in.symmetric),
     scaleForFFT (in.scaleForFFT),
     scaleRevFFT (in.scaleRevFFT)
{
   if (dir != None)  {
      if (in.inArray) createArrays (InArray);
      if (in.outArray) createArrays (OutArray);
   }
   
#  if defined (USE_FFTW)
      if (iPlanFor != -1)  fftPlanCounter(iPlanFor)++;
      if (iPlanRev != -1)  fftPlanCounter(iPlanRev)++;
#     if defined (USE_FFTW_PARALLEL)
         if (iPlanForParA != -1)  fftPlanCounterParA(iPlanForParA)++;
         if (iPlanRevParA != -1)  fftPlanCounterParA(iPlanRevParA)++;
         if (iPlanForParB != -1)  fftPlanCounterParB(iPlanForParB)++;
         if (iPlanRevParB != -1)  fftPlanCounterParB(iPlanRevParB)++;
#     endif // USE_FFTW_PARALLEL
#  elif defined (USE_ACML_FFT)
      if (iPlanFor != -1)  fftPlanCounter(iPlanFor)++;
      if (iPlanRev != -1)  fftPlanCounter(iPlanRev)++;
#  endif

}


SxFFT::~SxFFT ()
{
   destroy ();
}

void SxFFT::destroy ()
{

#  if defined (USE_FFTW)

      if (iPlanFor != -1 && (dir == Forward || dir == Both))  {
         int i = iPlanFor;
         fftPlanCounter(i)--;
         if (fftPlanCounter(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecs(i) = "";
            fftw_destroy_plan ((fftw_plan)plan_fwd);
            fftPlans(i) = NULL;
         }
         iPlanFor = -1;
         plan_fwd = NULL;
      }

      if (iPlanRev != -1 && (dir == Reverse || dir == Both))  {
         int i = iPlanRev;
         fftPlanCounter(i)--;
         if (fftPlanCounter(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecs(i) = "";
            fftw_destroy_plan ((fftw_plan)plan_rev);
            fftPlans(i) = NULL;
         }
         iPlanRev = -1;
         plan_rev = NULL;
      }

#     if defined (USE_FFTW_PARALLEL)

      if (iPlanForParA != -1 && (dir == Forward || dir == Both))  {
         int i = iPlanForParA;
         fftPlanCounterParA(i)--;
         if (fftPlanCounterParA(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecsParA(i) = "";
            fftw_destroy_plan (plan_fwd_par_a);
            fftPlansParA(i) = NULL;
         }
         iPlanForParA  = -1;
         plan_fwd_par_a = NULL;
      }

      if (iPlanRevParA != -1 && (dir == Reverse || dir == Both))  {
         int i = iPlanRevParA;
         fftPlanCounterParA(i)--;
         if (fftPlanCounterParA(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecsParA(i) = "";
            fftw_destroy_plan (plan_rev_par_a);
            fftPlansParA(i) = NULL;
         }
         iPlanRevParA = -1;
         plan_rev_par_a = NULL;
      }

      if (iPlanForParB != -1 && (dir == Forward || dir == Both))  {
         int i = iPlanForParB;
         fftPlanCounterParB(i)--;
         if (fftPlanCounterParB(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecsParB(i) = "";
            fftw_destroy_plan (plan_fwd_par_b);
            fftPlansParB(i) = NULL;
         }
         iPlanForParB  = -1;
         plan_fwd_par_b = NULL;
      }

      if (iPlanRevParB != -1 && (dir == Reverse || dir == Both))  {
         int i = iPlanRevParB;
         fftPlanCounterParB(i)--;
         if (fftPlanCounterParB(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecsParB(i) = "";
            fftw_destroy_plan (plan_rev_par_b);
            fftPlansParB(i) = NULL;
         }
         iPlanRevParB = -1;
         plan_rev_par_b = NULL;
      }
#     endif // USE_FFTW_PARALLEL

#  elif defined (USE_ACML_FFT)

      // --- destroy ACML plans
      //     note that ACML uses one plan for both directions
      if (iPlanFor != -1)  {
         int i = iPlanFor;
         fftPlanCounter(i)--;
         if (fftPlanCounter(i) == 0)  {
            //cout << "destroying " << fftPlanSpecs(i) << "..." << endl;
            fftPlanSpecs(i) = "";
            fftPlans(i) = NULL;
            delete [] plan_fwd;
            plan_rev = NULL;  // plan_rev just refers to plan_fwd!
         }
         iPlanFor = -1;
      }
#  endif /* USE_FFTW */

   if (dir != None) {
      destroyArrays ();
   }
   inArray = outArray = NULL;
   dir = None;
}


SxFFT &SxFFT::operator= (const SxFFT &in)
{
   if (&in == this)  return *this;

   destroy ();

   meshSize  = in.meshSize;
   if (in.inArray) createArrays (InArray);
   if (in.outArray) createArrays (OutArray);
   scaleFor  = in.scaleFor;
   scaleRev  = in.scaleRev;
#  if defined (USE_FFTW)
      plan_fwd      = in.plan_fwd;
      plan_rev      = in.plan_rev;
      iPlanFor      = in.iPlanFor; 
      iPlanRev      = in.iPlanRev;
      if (iPlanFor != -1)  fftPlanCounter(iPlanFor)++;
      if (iPlanRev != -1)  fftPlanCounter(iPlanRev)++;
#  if defined (USE_FFTW_PARALLEL)
      plan_fwd_par_a    = in.plan_fwd_par_a;
      plan_rev_par_a    = in.plan_rev_par_a;
      iPlanForParA      = in.iPlanForParA;
      iPlanRevParA      = in.iPlanRevParA;
      if (iPlanForParA != -1)  fftPlanCounterParA(iPlanForParA)++;
      if (iPlanRevParA != -1)  fftPlanCounterParA(iPlanRevParA)++;
      plan_fwd_par_b    = in.plan_fwd_par_b;
      plan_rev_par_b    = in.plan_rev_par_b;
      iPlanForParB      = in.iPlanForParB;
      iPlanRevParB      = in.iPlanRevParB;
      if (iPlanForParB != -1)  fftPlanCounterParB(iPlanForParB)++;
      if (iPlanRevParB != -1)  fftPlanCounterParB(iPlanRevParB)++;
#  endif // USE_FFTW_PARALLEL
#  elif defined (USE_ACML_FFT)
      plan_fwd      = in.plan_fwd;
      plan_rev      = in.plan_rev;
      iPlanFor      = in.iPlanFor; 
      iPlanRev      = in.iPlanRev;
      if (iPlanFor != -1)  fftPlanCounter(iPlanFor)++;
      if (iPlanRev != -1)  fftPlanCounter(iPlanRev)++;
#  else
      FFT_FORWARD = in.FFT_FORWARD;
      FFT_REVERSE = in.FFT_REVERSE;
#  endif     
   dir = in.dir;
   autoscale = in.autoscale;
   symmetric = in.symmetric;
   scaleForFFT = in.scaleForFFT;
   scaleRevFFT = in.scaleRevFFT;

   return *this;
}


void SxFFT::clean (FFT_COMPLEX *arrayPtr)
{
   if (!arrayPtr) arrayPtr = inArray;
   SX_CHECK (arrayPtr == inArray || arrayPtr == outArray);
   SX_CHECK (arrayPtr);

#  if   defined (USE_VECLIB_FFT)
      FFT_COMPLEX zero;
      zero.re = zero.im = 0.;
      for (int i=0; i<meshSize; i++)  *arrayPtr++ = zero;
#  elif defined (USE_ACML_FFT)
      FFT_COMPLEX zero; zero.real = zero.imag = 0.;
      for (int i=0; i<meshSize; i++)  *arrayPtr++ = zero;
#  elif defined (USE_ESSL)
      FFT_COMPLEX zero = FFT_COMPLEX (0., 0.);
      for (int i=0; i<meshSize; i++)  *arrayPtr++ = zero;
#  elif defined (USE_FFTSG)
      FFT_COMPLEX zero = FFT_COMPLEX (0., 0.);
      for (int i=0; i<meshSize; i++)  *arrayPtr++ = zero;
#  elif defined (USE_MKL_FFT)
      FFT_COMPLEX zero; zero.re = zero.im = 0.;
      for (int i=0; i<meshSize; i++)  *arrayPtr++ = zero;
#  else  /* FFTW */      
#ifdef USE_OPENMP
      if (meshSize > sxChunkSize)  {
#pragma omp parallel for
         for (int i=0; i<meshSize; i++)
            arrayPtr[i][0] = arrayPtr[i][1] = 0.;
      } else
#endif
      {
         for (int i=0; i<meshSize; i++, ++arrayPtr)  
            (*arrayPtr)[0] = (*arrayPtr)[1] = 0.;
      }
#  endif      
}



void SxFFT::randomize (FFT_COMPLEX *arrayPtr)
{
   if (!arrayPtr)  arrayPtr = inArray;
   SX_CHECK (arrayPtr == inArray || arrayPtr == outArray);
   SX_CHECK (arrayPtr);

#  if   defined (USE_VECLIB_FFT)
      FFT_COMPLEX rand;
      for (int i=0; i<meshSize; i++)  {
         rand.re = SxRandom::get();
         rand.im = SxRandom::get();
         *arrayPtr++ = rand;
      }
#  elif defined (USE_ACML_FFT)
      FFT_COMPLEX rand;
      for (int i=0; i<meshSize; i++)  {
         rand.real = SxRandom::get();
         rand.imag = SxRandom::get();
         *arrayPtr++ = rand;
      }
#  elif defined (USE_ESSL)
      FFT_COMPLEX rand;
      for (int i=0; i<meshSize; i++)  {
         rand = FFT_COMPLEX (SxRandom::get(), SxRandom::get());
         *arrayPtr++ = rand;
      }
#  elif defined (USE_FFTSG)
      FFT_COMPLEX rand;
      for (int i=0; i<meshSize; i++)  {
         rand = FFT_COMPLEX (SxRandom::get(), SxRandom::get());
         *arrayPtr++ = rand;
      }
#  elif defined (USE_MKL_FFT)
      FFT_COMPLEX rand;
      for (int i=0; i<meshSize; i++)  {
         rand.re = SxRandom::get();
         rand.im = SxRandom::get();
         *arrayPtr++ = rand;
      }
#  else  /* FFTW */      
      for (int i=0; i<meshSize; i++)  {
         (*arrayPtr)[0]   = SxRandom::get();
         (*arrayPtr++)[1] = SxRandom::get();
      }
#  endif      
}



void SxFFT::checkMeshSize (int checkSize)
{
   double size;

#  if   defined (USE_VECLIB_FFT)
      int a, b, c;
      double pow2a, pow3b;
      for (a=0; a<=25; a++)
         for (b=0, pow2a=::pow(2.,a); b<=25; b++)
            for (c=0, pow3b=::pow(3.,b); c<=25; c++)  {
               size = pow2a * pow3b  * ::pow(5.,c);
               if (fabs(checkSize - size) < 1e-3) {
                  // --- size is okay
                  return;
               }
            }
#  elif defined (USE_ACML_FFT)
      int a, b, c, d, e, f;
      double pow2a, pow3b, pow5c, pow7d, pow11e;
      for (a=0; a<15; a++)
         for (b=0, pow2a=::pow(2.,a); b<15; b++)
            for (c=0, pow3b=::pow(3.,b); c<15; c++)
               for (d=0, pow5c=::pow(5.,c); d<15; d++)
                  for (e=0, pow7d=::pow(7.,d); e<2; e++)
                     for (f=0, pow11e=::pow(11.,e); f<2; f++)
                        if (e+f < 2)  {
                           size = pow2a * pow3b  * pow5c
                                * pow7d * pow11e * ::pow(13.,f);
                           if (fabs(checkSize - size) < 1e-3) {
                              // --- size is okay
                              return;
                           }
                        }

      // --- no commensurable mesh found
#  elif defined (USE_ESSL)
      int a, b, c, d, e, f;
      double pow2a, pow3b, pow5c, pow7d, pow11e;
      for (a=1; a<=25; a++)
         for (b=0, pow2a=::pow(2.,a); b<=2; b++)
            for (c=0, pow3b=::pow(3.,b); c<=1 ; c++)
               for (d=0, pow5c=::pow(5.,c); d<=1 ; d++)
                  for (e=0, pow7d=::pow(7.,d); e<=1; e++)
                     for (f=0, pow11e=::pow(11.,e); f<=1; f++)  {
                        size = pow2a * pow3b  * pow5c
                             * pow7d * pow11e;
                        if (fabs(checkSize - size) < 1e-3) {
                           // --- size is okay
                           return;
                        }
                     }
#  elif defined (USE_FFTSG)
      int a, b, c;
      double pow2a, pow3b;
      for (a=0; a<=25; a++)
         for (b=0, pow2a=::pow(2.,a); b<=25; b++)
            for (c=0, pow3b=::pow(3.,b); c<=25; c++)  {
               size = pow2a * pow3b  * ::pow(5.,c);
               if (fabs(checkSize - size) < 1e-3) {
                  // --- size is okay
                  return;
               }
            }
 
#  else /* FFTW */
      int a, b, c, d, e, f;
      double pow2a, pow3b, pow5c, pow7d, pow11e;
      for (a=0; a<15; a++)
         for (b=0, pow2a=::pow(2.,a); b<15; b++)
            for (c=0, pow3b=::pow(3.,b); c<15; c++)
               for (d=0, pow5c=::pow(5.,c); d<15; d++)
                  for (e=0, pow7d=::pow(7.,d); e<2; e++)
                     for (f=0, pow11e=::pow(11.,e); f<2; f++)
                        if (e+f < 2)  {
                           size = pow2a * pow3b  * pow5c
                                * pow7d * pow11e * ::pow(13.,f);
                           if (fabs(checkSize - size) < 1e-3) {
                              // --- size is okay
                              return;
                           }
                        }

      sxprintf ("WARNING: One of the provided mesh dimensions (%d) may "
                "not be well suited for FFT (large prime factors)\n.",
                checkSize);
      return;
#  endif 
      sxprintf ("One of the provided mesh dimensions (%d) is uncommensurable "
                "to FFT.\n",
               checkSize);
      SX_EXIT;
}



// --------------------------------------------------------------------------
// How to generate these lists, e.g. for FFTW:
// see SxFFT::checkMeshSize
// --------- abc.cpp -------
// int main ()
// {
//    int a, b, c, d, e, f;
//    double pow2a, pow3b, pow5c, pow7d, pow11e;
//    double size;
//    for (a=0; a<15; a++)
//       for (b=0, pow2a=pow(2.,a); b<15; b++)
//          for (c=0, pow3b=pow(3.,b); c<15; c++)
//             for (d=0, pow5c=pow(5.,c); d<15; d++)
//                for (e=0, pow7d=pow(7.,d); e<2; e++)
//                   for (f=0, pow11e=pow(11.,e); f<2; f++)
//                      if (e+f < 2)  {
//                         size = pow2a * pow3b  * pow5c
//                            * pow7d * pow11e * pow(13.,f);
//                         if (size >= 6 && size <= 1000)
//                            sxprintf ("%g\n", size);
//                      }
//    return 0;
// }
// ----------- 8< ---------
//
// g++ list.cpp; a.out | sort -n  | uniq
//
SxList<int> SxFFT::getStdMeshSizes ()
{
   SxList<int> list;
#  if defined (USE_VECLIB)
   SX_EXIT;
      list = SxList<int> ()
           << 6 << 8 << 9 << 10 << 12 << 15 << 16 << 18 << 20 << 24 << 25
           << 27 << 30 << 32 << 36 << 40 << 45 << 48 << 50 << 54 << 60 << 64
           << 72 << 75 << 80 << 81 << 90 << 96 << 100 << 108 << 120 << 125 
           << 128 << 135 << 144 << 150 << 160 << 162 << 180 << 192 << 200 
           << 216 << 225 << 240 << 243 << 250 << 256 << 270 << 288 << 300
           << 320 << 324 << 360 << 375 << 384 << 400 << 405 << 432 << 450 
           << 480 << 486 << 500 << 512 << 540 << 576 << 600 << 625 << 640
           << 648 << 675 << 720 << 729 << 750 << 768 << 800 << 810 << 864 
           << 900 << 960 << 972 << 1000;
#  elif defined (USE_ESSL)
      list = SxList<int> ()
           << 6 << 8 << 10 << 12 << 14 << 16 << 18 << 20 << 22 << 24 << 28 
           << 30 << 32 << 36 << 40 << 42 << 44 << 48 << 56 << 60 << 64 << 66 
           << 70 << 72 << 80 << 84 << 88 << 90 << 96 << 110 << 112 << 120 
           << 126 << 128 << 132 << 140 << 144 << 154 << 160 << 168 << 176 
           << 180 << 192 << 198 << 210 << 220 << 224 << 240 << 252 << 256 
           << 264 << 280 << 288 << 308 << 320 << 330 << 336 << 352 << 360 
           << 384 << 396 << 420 << 440 << 448 << 462 << 480 << 504 << 512 
           << 528 << 560 << 576 << 616 << 630 << 640 << 660 << 672 << 704 
           << 720 << 768 << 770 << 792 << 840 << 880 << 896 << 924 << 960 
           << 990;
#  elif defined (USE_FFTSG)
      // taken from the fft3d-sg.f90, subroutine: CTRIG
      list = SxList<int> ()
       << 3 << 4 << 5 << 6 << 8 << 9 << 10 << 12 << 15 << 16 << 18 << 20 << 24
       << 25 << 27 << 30 << 32 << 36 << 40 << 45 << 48 << 50 << 54 << 60 << 64
       << 72 << 75 << 80 << 81 << 90 << 96 << 100 << 108 << 120 << 125 << 128
       << 135 << 144 << 150 << 160 << 162 << 180 << 192 << 200 << 216 << 225 
       << 240 << 243 << 250 << 256 << 270 << 288 << 300 << 320 << 324 << 360 
       << 375 << 384 << 400 << 405 << 432 << 450 << 480 << 486 << 500 << 512
       << 540 << 576 << 600 << 625 << 640 << 648 << 675 << 720 << 729 << 750 
       << 768 << 800 << 810 << 864 << 900 << 960 << 972 << 1000 << 1024
       << 1080 << 1152 << 1200 << 1250 << 1280 << 1296 << 1350 << 1440 
       << 1458 << 1500 << 1536 << 1600 << 1620 << 1728 << 1800 << 1920 
       << 1944 << 2000 << 2048;
#  else /* FFTW */
      list = SxList<int> ()
           << 6 << 8 << 10 << 12 << 14 << 16 << 18 << 20 << 22 << 24 << 28 
           << 30 << 32 << 36 << 40 << 42 << 44 << 48 << 56 << 60 << 64 << 66 
           << 70 << 72 << 80 << 84 << 88 << 90 << 96 << 110 << 112 << 120 
           << 126 << 128 << 132 << 140 << 144 << 154 << 160 << 168 << 176 
           << 180 << 192 << 198 << 210 << 220 << 224 << 240 << 252 << 256 
           << 264 << 280 << 288 << 308 << 320 << 330 << 336 << 352 << 360 
           << 384 << 396 << 420 << 440 << 448 << 462 << 480 << 504 << 512 
           << 528 << 560 << 576 << 616 << 630 << 640 << 660 << 672 << 704 
           << 720 << 768 << 770 << 792 << 840 << 880 << 896 << 924 << 960 
           << 990;
#  endif 
   return list;     
}


int SxFFT::getNextMeshSize (int minSize)
{
   int maxFFTSize = 0;
   SxList<int> fftList = getStdMeshSizes ();
   SxList<int>::Iterator it;
   for (it = fftList.begin(); it != fftList.end(); it++)  
      if ( *it > maxFFTSize )  maxFFTSize = *it;
   
   if (minSize > maxFFTSize)  {
      sxprintf ("ERROR: Needed mesh size is larger than largest ");
      sxprintf ("predefined mesh data.\n");
      sxprintf ("Should be > %d\n", minSize);
      sxprintf ("Automatic search... ");
      for (int res = minSize + 1; res < minSize * 2; ++res)  {
         int p = res;
         while (p % 2 == 0) { p /= 2; }
         while (p % 3 == 0) { p /= 3; }
         while (p % 5 == 0) { p /= 5; }
#ifdef USE_FFTW
         while (p % 7 == 0) { p /= 7; }
         while (p % 11 == 0) { p /= 11; }
#endif
         if (p == 1)  {
            // this number has only the above prime factors
            sxprintf("%d\n", res);
            return res;
         }
      }
      sxprintf ("failed.\n");
      SX_EXIT;
   }

   int res = -1;
   for (it = fftList.begin(); it != fftList.end(); ++it)  {
      if ( minSize <= *it )  {
         res = *it;
         break;
      }
   }
   return res;
}


void SxFFT::saveWisdom (const SxString &fileName)
{
#ifdef USE_FFTW
   if (SxLoopMPI::me () == 0) {
      FILE *outFile = fopen (fileName.ascii (), "w");
      if (outFile)  {
         fftw_export_wisdom_to_file(outFile);
         fclose (outFile);
      } else {
         cout << "Cannot open '" << fileName << "' for saving FFTW wisdom." 
              << endl;
      }
   }
#else
   SX_UNUSED (fileName);
#endif
}

void SxFFT::loadWisdom (const SxString &fileName, bool ignoreErrors)
{
#ifdef USE_FFTW
   FILE *fp = fopen(fileName.ascii (), "r");
   (cout << "Read FFTW wisdom from '" << fileName << "'... ").flush ();
   if (fp)  {
      bool ok = (fftw_import_wisdom_from_file (fp) != 0);
      fclose (fp);
      if (ok)  {
         (cout << "done." << endl).flush ();
      } else {
         (cout << "failed in FFTW." << endl).flush ();
         if (!ignoreErrors) { SX_QUIT; }
      }
   } else {
      cout << "failed." << endl;
      cout << "Cannot open '" << fileName << "' for reading FFTW wisdom." 
           << endl;
      if (!ignoreErrors) { SX_QUIT; }
   }
#else
   SX_UNUSED (fileName);
   SX_UNUSED (ignoreErrors);
#endif
}

void SxFFT::createArrays (ArrayType type)
{
   // cout << "SxFFT::createArrays called" << endl << flush;
   if (type & (InArray | InArrayZero))  {
      SX_CHECK (inArray == NULL);
#ifdef USE_FFTW
      inArray = (FFT_COMPLEX*)fftw_malloc(meshSize * sizeof(FFT_COMPLEX));
      if (!inArray) sxOutOfMemoryHandler ();
#else
      inArray = new FFT_COMPLEX [meshSize];
#endif
      if (type & InArrayZero) clean (inArray);
   }
   if (type & (OutArray | OutArrayZero))  {
      SX_CHECK (outArray == NULL);
#ifdef USE_FFTW
      outArray = (FFT_COMPLEX*)fftw_malloc(meshSize * sizeof(FFT_COMPLEX));
      if (!outArray) sxOutOfMemoryHandler ();
#else
      outArray = new FFT_COMPLEX [meshSize];
#endif
      if (type & OutArrayZero) clean (outArray);
   }
}

void SxFFT::destroyArrays ()
{
   // cout << "SxFFT::destroyArrays called" << endl << flush;
#ifdef USE_FFTW
   if (inArray) fftw_free (inArray);
   if (outArray) fftw_free (outArray);
#else
   if (inArray) delete [] inArray;
   if (outArray) delete [] outArray;
#endif
   inArray = outArray = NULL;
}

#ifdef USE_FFTW
void SxFFT::broadcastPlans()
{
#ifdef USE_LOOPMPI
   char * planString = NULL;
   int planStringSize = 0;
   if (SxLoopMPI::me() == 0) {
      planString = fftw_export_wisdom_to_string();
      planStringSize = (int)SxString(planString).getSize();
   }
   planStringSize = SxLoopMPI::bcast(planStringSize, 0);
   if (! planStringSize) {
      return;
   }
   planStringSize++; // '\0' at the end
   if (SxLoopMPI::me() != 0) {
      planString = (char*) malloc(planStringSize * sizeof(char));
   }
   SxLoopMPI::bcast(planString, planStringSize, 0);
   if (SxLoopMPI::me() != 0) {
      fftw_forget_wisdom();
      if (!fftw_import_wisdom_from_string(planString))  {
         cout << "Failed to import wisdom from MPI master!" << endl;
         SX_EXIT;
      }
   }
   free(planString);
#endif
}
#endif

void SxFFT::plannerCLI (SxCLI &cli)
{
#if defined USE_FFTW && ! defined USE_MKL_FFT
   if (cli.option ("--fastfft", "mode",
                   "request fast FFT (FFTW_ESTIMATE) planning.\n Optional mode (=estimate,=measure,=patient) sets FFTW planning mode").toBool ())  {
      if (cli.last ().hasValue ())  {
         SxString mode = cli.last ().getValue ().toLower ();
         if (mode == "estimate")
            SxFFT::quickFFTPlanner (Estimate);
         else if (mode == "measure")
            SxFFT::quickFFTPlanner (Measure);
         else if (mode == "patient")
            SxFFT::quickFFTPlanner (Patient);
         else  {
            cout << "Unknown FFT planning mode: " << mode << endl;
            cli.setError ();
         }
      } else {
         SxFFT::quickFFTPlanner (Estimate);
      }
   }
   SxString fftWisdom = cli.option("--wisdom","file","FFT wisdom file")
                        .toString ("");
   if (fftWisdom.getSize () > 0 && ! cli.error)
      loadWisdom (fftWisdom);
#else
   if (cli.option ("--fastfft", "mode",
                   "(option not available for FFT library)").toBool ())
   {
      cout << "WARNING: --fastfft option is not available for FFT library"
           << endl;
   }
#endif
}
