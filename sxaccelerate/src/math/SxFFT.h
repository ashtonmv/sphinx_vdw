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


// use the parallel implementation of the 3D FFT in case FFTW is used
// #define USE_FFTW   // satisfy Eclipse ...


// #define USE_FFTW_PARALLEL
#undef USE_FFTW_PARALLEL


#ifndef _SX_FFT_H_
#define _SX_FFT_H_

#include <SxString.h>
#include <SxList.h>
#include <SxPrecision.h>
#include <SxMath.h>
#include <SxConfig.h>
#include <SxTimer.h>
class SxCLI;


// --- FFTW
#if defined (USE_MKL_FFT)
#   include <fftw/fftw3.h>
#   define USE_FFTW

#   define FFT_REAL     double
   typedef struct { double re, im; } FFT_COMPLEX;
#   define SX_REAL(x)   x.re
#   define SX_IMAG(x)   x.im

#   define FFT_ROW_MAJOR  true


#elif defined (USE_FFTW)
#   include <fftw3.h>

#   define FFT_REAL     double
#   define FFT_COMPLEX  fftw_complex
#   define SX_REAL(x)   x.re
#   define SX_IMAG(x)   x.im

#   define FFT_ROW_MAJOR  true

// --- Stefan Goedecker FFT
#elif defined (USE_FFTSG)
//#   include <SPHInX.h>
#   define FFT_REAL     double
#   define FFT_COMPLEX  SxComplex16
#   define SX_REAL(x)   x.re
#   define SX_IMAG(x)   x.im

#   define FFT_ROW_MAJOR  false

// --- HP VECLIB
#elif defined (USE_VECLIB_FFT)

//#   include <veclib.h>  // --- don't use HP header files!
                          //     more details in SxMLIB.h
#     include <SxMLIB.h>
#   ifdef FFT_SINGLE_PREC
#      define FFT_REAL     float
#      define FFT_COMPLEX  complex8_t
#   else
#      define FFT_REAL     double
#      define FFT_COMPLEX  complex16_t
#   endif  /* _FFT_SINGLE_PROC */
#   define SX_REAL(x)   x.re
#   define SX_IMAG(x)   x.im

#   define FFT_ROW_MAJOR  false

// --- AMD Core Math Library
#elif defined (USE_ACML_FFT)

#   include <acml.h>
#   ifdef FFT_SINGLE_PREC
#      define FFT_REAL     float
#      define FFT_COMPLEX  complex
#   else
#      define FFT_REAL     double
#      define FFT_COMPLEX  doublecomplex
#   endif  /* _FFT_SINGLE_PROC */
#   define SX_REAL(x)   x.real
#   define SX_IMAG(x)   x.imag

#   define FFT_ROW_MAJOR  false


// --- IBM ESSL
#elif defined (USE_ESSL)
    // --- before including ESSL we have to do a messy definition
    //     otherwise it doesn't compile. i don't know a better way.
    //     see also: SxMath.cpp
#   define _ESV_COMPLEX_
#   include <complex>
#   include <essl.h>
#   ifdef FFT_SINGLE_PREC
#      define FFT_REAL     float
#      define FFT_COMPLEX  complex<float>
#      define SX_REAL(x)   sreal(x)
#      define SX_IMAG(x)   simag(x)
#   else 
#      define FFT_REAL     double
#      define FFT_COMPLEX  complex<double>
#      define SX_REAL(x)   real(x)
#      define SX_IMAG(x)   imag(x)
#   endif  /* FFT_SINGLE_PROC */

#   define FFT_ROW_MAJOR  false
#else
#   error "SxFFT.h: No FFT library specified."
#endif  /* USE_FFTW */



/** \brief FFT Basis class

    \b SxFFT = SPHInX Fast Fourier Transformation basis

    This class contains the functions and macros needed for all 
    Fast Fourier Transformations regardless the dimensionality of the
    transformation. The actual transformations are done in the derived
    classes.

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_MATH SxFFT
{
   public:
      enum Directions { None, Forward, Reverse, Both };
      int meshSize;

      /**  \brief Pointer to local input data array */
      FFT_COMPLEX *inArray;
      /**  \brief Pointer to local output data array */
      FFT_COMPLEX *outArray;

      /** \brief Forward FFT scaling factor */
      SX_FFT_REAL  scaleFor;
      /** \brief Reverse FFT scaling factor */
      SX_FFT_REAL  scaleRev;


#     if defined (USE_FFTW)

            /** \brief FFTW plan for forward FFT

                   \sa http://www.fftw.org/fftw3_doc/Using-Plans.html#Using%20Plans*/
            void *plan_fwd;
            /** \brief FFTW plan for forward FFT

                   \sa http://www.fftw.org/fftw3_doc/Using-Plans.html#Using%20Plans*/
            void *plan_rev;
            /** \brief Index of the FFTW plan number in the list fftPlans
             */
            int iPlanFor;
            /** \brief Index of the FFTW plan number in the list fftPlans
             */
            int iPlanRev;

            /** \brief Push FFTW plans from rank 0 to all higher ranks */
            void broadcastPlans();

#        if defined (USE_FFTW_PARALLEL)
            /** \brief First FFTW plan for the forward FFT */
            void *plan_fwd_par_a;
            /** \brief First FFTW plan for the reverse FFT */
            void *plan_rev_par_a;
            /** \brief Index of the first FFTW plan number in the list */
            int iPlanForParA;
            /** \brief Index of the first FFTW plan number in the list */
            int iPlanRevParA;
            /** \brief Second FFTW plan for the forward FFT */
            void *plan_fwd_par_b;
            /** \brief Second FFTW plan for the reverse FFT */
            void *plan_rev_par_b;
            /** \brief Index of the second FFTW plan number in the list fftPlans */
            int iPlanForParB;
            /** \brief Index of the second FFTW plan number in the list fftPlans */
            int iPlanRevParB;
#        endif


#     elif defined (USE_ACML_FFT)
         /** \brief Forward communication array.

              The ACML communication arrays are comparible to FFTW's plans. This is
              the forward plan.*/
         FFT_COMPLEX *plan_fwd;

         /** \brief Reverse communication array.

              The ACML communication arrays are comparible to FFTW's plans. This is
              the reverse plan.*/
       	 FFT_COMPLEX *plan_rev;

         /** \brief Index of the plan number in the list fftPlans 
          */         
         int iPlanFor;
         /** \brief Index of the plan number in the list fftPlans 
           */         
         int iPlanRev;
#     else      
         int FFT_FORWARD;
         int FFT_REVERSE;
#     endif /* USE_FFTW */      

      
      SxFFT (Directions dirIn=None, 
             bool autoscaleIn=true, 
             bool symmetricIn=true);
      SxFFT (const SxFFT &);
      ~SxFFT ();

      /** Destroy all allocated arrays  and plans */
      void destroy ();

      SxFFT &operator= (const SxFFT &);

      /** \brief Set an element of the input data array

          This function modifies the \em i-th element of the input
          data array (::inArray) with a provided value. */
      inline void setElement (int i, const SxComplex16 &v)  {
         SX_CHECK (i>=0 && i < meshSize, i, meshSize);
#        if   defined (USE_ESSL)
            inArray[i] = FFT_COMPLEX ((FFT_REAL)v.re, (FFT_REAL)v.im);
#        elif defined (USE_VECLIB_FFT)
            inArray[i].re = (FFT_REAL)v.re;
            inArray[i].im = (FFT_REAL)v.im;
#        elif defined (USE_ACML_FFT)
            inArray[i].real = (FFT_REAL)v.re;
            inArray[i].imag = (FFT_REAL)v.im;
#        elif defined (USE_FFTSG)
            inArray[i] = v;
            inArray[i].re = (FFT_REAL)v.re;
            inArray[i].im = (FFT_REAL)v.im;
#        elif defined (USE_MKL_FFT)
            inArray[i].re = (FFT_REAL)v.re;
            inArray[i].im = (FFT_REAL)v.im;
#        else /* USE_FFTW */            
            inArray[i][0] = (FFT_REAL)v.re;
            inArray[i][1] = (FFT_REAL)v.im;
#        endif            
      }


      /** \brief Initialize a provided FFT mesh with zeros. */
      void clean (FFT_COMPLEX *);

      /** \brief for benchmarks only */
      void randomize (FFT_COMPLEX *arrayPtr=NULL);

      /** \brief Creates list of allowed FFT mesh sizes

          Different FFT libraries can deal with different FFT mesh sizes.
          This function returns the proper list of allowed mesh sizes on the
          currently linked FFT library. 
          \sa checkMeshSize */
      static SxList<int> getStdMeshSizes ();
      /** \brief Get the next larger supported FFT mesh size 
       
          FFT libraries support only certain grid sizes. If the demanded
          size (minSize) is supported by the currently linked FFT library
          minSize is returned. Otherwise the next larger allowed FFT size
          is returned.  */
      static int getNextMeshSize (int minSize);

      /// \name Saving wisdom/plans
      /// @{
      /** \brief FFT library supports savable wisdom
          \todo This is presently only used for FFTW. Other libraries need
                to be checked.
        */
      static bool hasWisdom ()
      {
#        ifdef USE_FFTW
            return true;
#        else
            return false;
#        endif
      }
         
      /// Save accumulated FFT wisdom
      static void saveWisdom (const SxString &fileName);
      /// Load FFT wisdom
      static void loadWisdom (const SxString &fileName, 
                              bool ignoreErrors = false);

      enum PlanMode { Estimate, Measure, Patient };
#ifdef USE_FFTW
      static unsigned fftPlanMode;
      static void
      quickFFTPlanner (PlanMode mode = Estimate, unsigned *oldMode = NULL)
      {
         if (oldMode) *oldMode = fftPlanMode;
         if (mode == Estimate) fftPlanMode = FFTW_ESTIMATE;
         else if (mode == Measure) fftPlanMode = FFTW_MEASURE;
         else if (mode == Patient) fftPlanMode = FFTW_PATIENT;
         else { cout << "Unknown FFT planning mode" << endl; SX_EXIT; }
      }
      static void restorePlannerMode (unsigned fftw_mode) {
         fftPlanMode = fftw_mode;
      }
#else
      static void quickFFTPlanner (PlanMode, unsigned * = NULL)  { /* empty */ }
      static void restorePlannerMode (unsigned) { /* empty */ }
#endif
      // --- set planning mode/wisdom from CLI
      static void plannerCLI (SxCLI &);
      ///@}

// protected:

      Directions dir;
      bool autoscale;
      bool symmetric;


      /** \brief Validate the FFT mesh size.

          This function validates whether the currently linked FFT library
          can treat the provided FFT mesh size. 
          \sa getStdMeshSizes */
      static void checkMeshSize (int);

      /** \brief Internal scaling factors to achieve consitency with all
          FFT libraries. */
      SX_FFT_REAL  scaleForFFT;
      /** \brief Internal scaling factors to achieve consitency with all
          FFT libraries. */
      SX_FFT_REAL  scaleRevFFT;

      enum ArrayType { InArray = 0x01, InArrayZero = 0x02, OutArray=0x10, OutArrayZero = 0x20,
                       InOutArray=0x11, InOutZero = 0x22 };
      void createArrays (ArrayType type);
      void destroyArrays ();

};


extern SxList<SxString>   fftPlanSpecs;
extern SxList<void *>     fftPlans;
extern SxList<int>        fftPlanCounter;
#  if defined (USE_FFTW_PARALLEL)
   extern SxList<SxString>   fftPlanSpecsParA;
   extern SxList<SxString>   fftPlanSpecsParB;
   extern SxList<void *>     fftPlansParA;
   extern SxList<void *>     fftPlansParB;
   extern SxList<int>        fftPlanCounterParA;
   extern SxList<int>        fftPlanCounterParB;
#  endif

namespace Timer {
   enum FFTTimer { FFTPlanning } ;
}

SX_REGISTER_TIMERS(Timer::FFTTimer)
{
   using namespace Timer;
#ifdef USE_FFTW
   if (SxFFT::fftPlanMode == FFTW_ESTIMATE)
      regTimer (FFTPlanning, "FFT planning (estim)");
   else if (SxFFT::fftPlanMode == FFTW_MEASURE)
      regTimer (FFTPlanning, "FFT planning (meas)");
   else if (SxFFT::fftPlanMode == FFTW_PATIENT)
      regTimer (FFTPlanning, "FFT planning (pat)");
   else
#endif
   {
      regTimer (FFTPlanning, "FFT planning");
   }
}

#endif /* _SX_FFT_H_ */
