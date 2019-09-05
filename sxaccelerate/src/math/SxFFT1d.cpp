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

#include <SxFFT1d.h>
#include <SxVector.h>
#include <SxLoopMPI.h>

#if defined (USE_MKL_FFT)
#  include <fftw/fftw3.h>
#elif defined (USE_FFTW)
#  include <fftw3.h>
#endif /* USE_MKL_FFTW */


SxFFT1d::SxFFT1d () : SxFFT (None, true, true)
{
   // empty
}


SxFFT1d::SxFFT1d (Directions dir_, int nx, double omega, int strideIn, 
                  bool symmetric_)
   : SxFFT (dir_, true, symmetric_),
     stride (strideIn)
{
   setMesh (nx, omega, stride);
}


SxFFT1d::SxFFT1d (const SxFFT1d &in)
   : SxFFT (in),
     mesh (in.mesh),
#    if defined (USE_FFTSG)
        perfMesh (in.perfMesh),
        fftCacheSize (in.fftCacheSize),
#    endif   
     stride (in.stride)
{
   // empty
} 


SxFFT1d::~SxFFT1d ()
{
   // FFTW plans as well as inArray and outArray are destroyed
   // by SxFFT::SxFFT
}


SxFFT1d &SxFFT1d::operator= (const SxFFT1d &in)
{
   if (&in == this)  return *this;

   SxFFT::operator= (in);

   mesh = in.mesh;
#  if defined (USE_FFTSG)
      perfMesh     = in.perfMesh;
      fftCacheSize = in.fftCacheSize;
#  endif
   stride = in.stride;

   return *this;
}


void SxFFT1d::setMesh (int nx, double omega, int stride_)
{
   SX_CHECK (stride_ == 1 || (stride_ > 1 && supportsStridedFFT()));
   stride = stride_;
   
   //sxprintf ("FFT mesh: %d\n", nx);
   meshSize = nx;
#  if defined (USE_VECLIB_FFT)
      SX_CHECK (symmetric); 
      scaleRevFFT        = sqrt(omega) / (double)meshSize;
      scaleForFFT        = (double)meshSize / sqrt(omega);
      scaleFor           = 1. / sqrt(omega);
      scaleRev           = scaleRevFFT;
      FFT_FORWARD        = -1;
      FFT_REVERSE        = +1;
#  elif defined (USE_ACML_FFT)
      if (symmetric)  {
         scaleForFFT     = 1./sqrt(omega);
         scaleRevFFT     = sqrt(omega)/double(meshSize);
         scaleFor        = scaleForFFT;
         scaleRev        = scaleRevFFT;
      }  else  {
         scaleRevFFT     =  1. / ((double)meshSize);
         scaleForFFT     =  1.;
         scaleFor        = scaleForFFT;
         scaleRev        = scaleRevFFT;
      }
#  elif defined (USE_ESSL)
      SX_CHECK (symmetric); 
      scaleRevFFT        =  sqrt(omega) / ((double)meshSize);
      scaleForFFT        =  1. / sqrt(omega);
      scaleFor           = scaleForFFT;
      scaleRev           = scaleRevFFT;
      FFT_FORWARD        = -1;
      FFT_REVERSE        = +1;
#  elif defined (USE_FFTSG)
      SX_CHECK (symmetric); 
      fftCacheSize       = 4*1024;
      sxprintf ("TODO: FFTSG - NCACHE not yet optimized\n");
      scaleRevFFT        =  sqrt(omega) / ((double)meshSize);
      scaleForFFT        =  1. / sqrt(omega);
      scaleFor           = scaleForFFT;
      scaleRev           = scaleRevFFT;
      FFT_FORWARD        = -1;
      FFT_REVERSE        = +1;
      // FIXME: perfMesh = SxVector1<Int> (nx, ny, nz);
      perfMesh = nx;
      sxprintf (">>> TODO: perfMesh not yet optimized!\n");
#  else  /* FFTW */
      if (symmetric)  {
         scaleRevFFT     =  sqrt(omega) / ((double)meshSize);
         scaleForFFT     =  1. / sqrt(omega);
         scaleFor        = scaleForFFT;
         scaleRev        = scaleRevFFT;
      }  else  {
         scaleRevFFT     =  1. / ((double)meshSize);
         scaleForFFT     =  1.;
         scaleFor        = scaleForFFT;
         scaleRev        = scaleRevFFT;
      }
#  endif
   mesh = nx;
   destroyArrays ();
   meshSize = stride * nx; // for allocating arrays
   createArrays (InOutArray);
   meshSize = nx; 
   
   checkMeshSize (nx);

#  if defined (USE_FFTW)
      int i;
      SxString spec = SxString(nx) + "-" + stride;

      if (dir == Forward || dir == Both)  {
         i = iPlanFor = int(fftPlanSpecs.findPos (spec + ".fwd"));
         if (i >= 0)  {
            fftPlanCounter(i)++;
            plan_fwd = fftPlans(i);
         }  else  {
            SX_CLOCK (Timer::FFTPlanning);
//          plan_fwd = fftw_plan_dft_1d (nx,
//                                       inArray, outArray,
//                                       FFTW_BACKWARD,// see FFTW definition
//                                       fftPlanMode | FFTW_PRESERVE_INPUT);
            plan_fwd = fftw_plan_many_dft (1, &nx, 1, 
                                           (fftw_complex*)inArray,  NULL, stride, 0,
                                           (fftw_complex*)outArray, NULL, stride, 0,
                                           FFTW_BACKWARD,// see FFTW definition
                                           fftPlanMode | FFTW_PRESERVE_INPUT);
            if (!plan_fwd)  { sxprintf ("Couldn't create FFT plan.\n"); SX_EXIT; }
            fftPlanSpecs   << (spec + ".fwd");
            fftPlans       << plan_fwd;
            fftPlanCounter << 1;
            iPlanFor = int(fftPlanCounter.getSize() - 1);
         }

      } else {
         plan_fwd = NULL;
      }

      if (dir == Reverse || dir == Both)  {
         i = iPlanRev = int (fftPlanSpecs.findPos (spec + ".rev"));
         if (i >= 0)  {
            fftPlanCounter(i)++;
            plan_rev = fftPlans(i);
         }  else  {
            SX_CLOCK (Timer::FFTPlanning);
            plan_rev = fftw_plan_many_dft (1, &nx, 1, 
                                           (fftw_complex*)inArray,  NULL, stride, 0,
                                           (fftw_complex*)outArray, NULL, stride, 0,
                                           FFTW_FORWARD,// see FFTW definition
                                           fftPlanMode | FFTW_PRESERVE_INPUT);
            if (!plan_rev)  { sxprintf ("Couldn't create FFT plan.\n"); SX_EXIT; }
            fftPlanSpecs   << (spec + ".rev");
            fftPlans       << plan_rev;
            fftPlanCounter << 1;
            iPlanRev = int (fftPlanCounter.getSize() - 1);
         }
      } else {
         plan_rev = NULL;
      }
#  elif defined (USE_ACML_FFT)
      int i;
      SxString spec = SxString(nx) + "-" + stride;
      
      // --- create plan for both forward and reverse FFT
      //     note that ACML uses one plan for both directions
      i = iPlanFor = iPlanRev = int(fftPlanSpecs.findPos (spec + ".fwd"));
      if (i >= 0)  {
         fftPlanCounter(i)++;
         plan_fwd = (FFT_COMPLEX*)fftPlans(i);
      }  else  {
         int error = 0; // createPlan = 100; // see ACML docu
         plan_fwd = new FFT_COMPLEX [3*nx+100];      // see ACML manual
#        if defined (FFT_SINGLE_PREC)
            cfft1dx (100, 1.0, false, nx, inArray, 1, outArray, 1, 
                     plan_fwd, &error);
#        else
            zfft1dx (100, 1.0, false, nx, inArray, 1, outArray, 1, 
                     plan_fwd, &error);
#        endif /* FFT_SINGLE_PREC */
         if (error)  { sxprintf ("Couldn't create FFT plan. \n"); SX_EXIT; }
         fftPlanSpecs   << spec;
         fftPlans       << (void*)plan_fwd;
         fftPlanCounter << 1;
         iPlanRev = iPlanFor = int (fftPlanCounter.getSize() - 1);
      }
      plan_rev = plan_fwd;

#  endif  /* USE_FFTW */
      if (stride > 1) destroyArrays ();

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT1d::fftReverse (int n, const void *in, void *out)
{
   SX_CHECK (in != out);  // arrays *MUST* be different
#  if defined (USE_VECLIB_FFT)
      SX_EXIT;  // to yet implemented
//    int info;
//    // --- copy input array -> working array
//    FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//    FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
//    for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
//       c1dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
//               &mesh.v[2], &mesh.v[1], &FFT_REVERSE, &info);
#     else
//       z1dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
//               &mesh.v[2], &mesh.v[1], &FFT_REVERSE, &info);
#     endif  /* FFT_SINGLE_PREC */         
//    if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
#  elif defined (USE_ACML_FFT)
      SX_CHECK (plan_rev);
      int error = 0;
      FFT_REAL fac = (autoscale) ? scaleRevFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft1dx (-1, fac, false, n, 
                  (FFT_COMPLEX *)const_cast<void *>(in), 1, 
                  (FFT_COMPLEX *)out, 1, plan_rev, &error);
#     else
         zfft1dx (-1, fac, false, n, 
                  (FFT_COMPLEX *)const_cast<void *>(in), 1, 
                  (FFT_COMPLEX *)out, 1, plan_rev, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
#  elif defined (USE_ESSL)
      SX_EXIT; // not yet implemented
#     if defined (FFT_SINGLE_PREC)
//       scft1 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
//              (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
//              mesh.v[2], mesh.v[1], mesh.v[0],
//              FFT_REVERSE, scaleRevFFT, NULL, 0);
#     else         
//       dcft1 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
//              (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
//              mesh.v[2], mesh.v[1], mesh.v[0],
//              FFT_REVERSE, scaleRevFFT, NULL, 0);
#     endif  /* FFT_SINGLE_PREC */

#  elif defined (USE_FFTSG)
      SX_EXIT; // not yet implemented

//    sxprintf ("TODO: FFTSG still copies vectors!!!\n");
//    FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//    FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out;
//    SxVector<Complex16> inOut (2*meshSize);
//    int i;
//    for (i=0; i < n; i++)  inOut(i) = srcPtr[i];
//
//    int inzee =  1;
//    int isign = -1;
//    FCALL (fft) (FINT(mesh.v[2]),  FINT(mesh.v[1]),  FINT(mesh.v[0]),
//              FINT(perfMesh.v[2]), FINT(perfMesh.v[1]), FINT(perfMesh.v[0]),
//              FCOMPLEX16VEC(inOut),
//              FINT(isign), FINT(inzee), FINT(fftCacheSize));
//
//    for (i=0; i < n; i++)  destPtr[i]  = inOut(i+n);
//    if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
         
#  else  /* FFTW */
      SX_CHECK (plan_rev);
      fftw_execute_dft ((fftw_plan)plan_rev, 
                        static_cast<fftw_complex *>(const_cast<void *>(in)),
                        (fftw_complex *)out);
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
#  endif

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT1d::fftForward (const int n, const void *in, void *out)
{
   SX_CHECK (in != out);  // arrays *MUST* be different

#  if defined (USE_VECLIB_FFT)
      SX_EXIT; // not yet implemented
//    int info;
//    // --- copy input array -> working array
//    FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//    FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
//    for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
//       c1dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
//               &mesh.v[2], &mesh.v[1], &FFT_FORWARD, &info);
#     else
//       z1dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
//               &mesh.v[2], &mesh.v[1], &FFT_FORWARD, &info);
#     endif  /* FFT_SINGLE_PREC */         
//    if (autoscale && symmetric)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);
#  elif defined (USE_ACML_FFT)
      SX_CHECK (plan_fwd);
      int error = 0;
      FFT_REAL fac = (autoscale) ? scaleForFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft1dx (1, fac, false, n, 
                  (FFT_COMPLEX *)const_cast<void *>(in), 1, 
                  (FFT_COMPLEX *)out, 1, plan_fwd, &error);
#     else
         zfft1dx (1, fac, false, n, 
                  (FFT_COMPLEX *)const_cast<void *>(in), 1, 
                  (FFT_COMPLEX *)out, 1, plan_fwd, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
#  elif USE_ESSL
      SX_EXIT; // not yet implemented
#     if defined (FFT_SINGLE_PREC)
//       scft1 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
//              (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.[v1],
//              mesh.v[2], mesh.v[1], mesh.v[0],
//              FFT_FORWARD, scaleForFFT, NULL, 0);
#     else 
//       dcft1 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
//              (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
//              mesh.v[2], mesh.v[1], mesh.v[0],
//              FFT_FORWARD, scaleForFFT, NULL, 0);
#     endif  /* FFT_SINGLE_PREC */         

#  elif defined (USE_FFTSG)
      SX_EXIT; // not yet implemented
//    SX_CHECK (out-in == meshSize, out-in, meshSize);
//    FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//    FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out;
//    SxVector<Complex16> inOut (2*meshSize);
//    int i;
//    for (i=0; i < n; i++)  inOut(i) = srcPtr[i];
//
//    int inzee = 1;
//    int isign = 1;
//    FCALL (fft) (FINT(mesh.v[2]),  FINT(mesh.v[1]),  FINT(mesh.v[0]),
//              FINT(perfMesh.v[2]), FINT(perfMesh.v[1]), FINT(perfMesh.v[0]),
//              FCOMPLEX16VEC(inOut),
//              FINT(isign), FINT(inzee), FINT(fftCacheSize));
//
//    for (i=0; i < n; i++)  destPtr[i]  = inOut(i+n);
//    if (autoscale && symmetric)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);
//
#  else /* FFTW */
      SX_CHECK (plan_fwd);
      fftw_execute_dft ((fftw_plan)plan_fwd, 
                       static_cast<fftw_complex *>(const_cast<void *>(in)),
                        (fftw_complex *)out);
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);
#  endif  /* USE_FFTW */
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT1d::print ()
{
   print (inArray);
}

void SxFFT1d::print (FFT_COMPLEX *array)
{
   int i;
#  ifndef USE_FFTW   
      FFT_COMPLEX c;
#  endif
#  ifdef USE_MKL_FFT
      FFT_COMPLEX c;
#  endif
   for (i=0; i<mesh; i++)  {
#     if defined (USE_MKL_FFT)
         c = array[getFFTIdx(i)];
         sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#     elif defined (USE_FFTW)
       sxprintf ("(%g,%g)  ", array[getFFTIdx(i)][0],
                              array[getFFTIdx(i)][1]);
#     else               
         c = array[getFFTIdx(i)];
         sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#     endif               
   }
   sxprintf ("\n");
}

