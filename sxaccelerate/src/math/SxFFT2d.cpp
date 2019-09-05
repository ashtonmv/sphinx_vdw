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

#include <SxFFT2d.h>
// #include <SxVector.h>
#include <SxLoopMPI.h>
#include <SxBlasLib.h>

SxFFT2d::SxFFT2d () : SxFFT (None, true, true)
{
   // empty
}


SxFFT2d::SxFFT2d (Directions dir_, int nx, int ny, double omega,
                  bool symmetric_)
   : SxFFT (dir_, true, symmetric_)
{
   setMesh (nx, ny, omega);
}

SxFFT2d::SxFFT2d (const SxFFT2d &in)
   : SxFFT (in)
#    if defined (USE_FFTSG)
        ,fftCacheSize (in.fftCacheSize)
#    endif
{
     mesh[0] = in.mesh[0];
     mesh[1] = in.mesh[1];
#    if defined (USE_FFTSG)
       perfMesh[0] = in.perfMesh[0];
       perfMesh[1] = in.perfMesh[1];
#    endif
}

SxFFT2d::~SxFFT2d ()
{
   // FFTW plans as well as inArray and outArray are destroyed
   // by SxFFT::SxFFT
}



// among others, creates plans for FFTW
void SxFFT2d::setMesh (int nx, int ny, double omega)
{

   bool createdNewPlan = false;
   meshSize = nx * ny;

#  if   defined (USE_VECLIB_FFT)
      SX_CHECK (symmetric); 
      scaleRevFFT        = sqrt(omega) / meshSize;
      scaleForFFT        = meshSize / sqrt(omega);
      scaleFor           = 1. / sqrt(omega);
      scaleRev           = scaleRevFFT;
      FFT_FORWARD        = -1;
      FFT_REVERSE        = +1;
#  elif defined (USE_ACML_FFT)
      if (symmetric)  {
         scaleRevFFT     = sqrt(omega)/(double)meshSize;
         scaleForFFT     = 1./sqrt(omega);
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
      perfMesh[0] = nx;
      perfMesh[1] = ny;
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
      
   destroyArrays ();
   createArrays (InOutArray);

#ifndef NDEBUG
   {
      // --- 'invalidate' arrays
      double realNan = sqrt(-1.);
      SxComplex16 nan(realNan, realNan);
      // inArray
      for (int i = 0; i < meshSize; i++)
         setElement (i, nan);
      // outArray via inArray
      FFT_COMPLEX *saveInArray = inArray;
      inArray = outArray;
      for (int i = 0; i < meshSize; i++)
         setElement (i, nan);
      inArray = saveInArray;
   }
#endif

   mesh[0] = nx;
   mesh[1] = ny;
   
   checkMeshSize (nx);
   checkMeshSize (ny);

#  if defined (USE_FFTW)

      int i;
      SxString spec = SxString(nx) + "x" + ny;

      if (dir == Forward || dir == Both)  {

         SxString specDir = spec + "FWD";

         // initialization of the 2D FFT, forward direction
         i = iPlanFor = int(fftPlanSpecs.findPos (specDir));
         if (i >= 0)  {
            plan_fwd = fftPlans(i);
            fftPlanCounter(i)++;
         }  else  {
            SX_CLOCK (Timer::FFTPlanning);
            plan_fwd = fftw_plan_dft_2d (nx, ny,
                                         (fftw_complex*)inArray,
                                         (fftw_complex*)outArray,
                                         FFTW_BACKWARD,// see FFTW definition
                                         fftPlanMode | FFTW_PRESERVE_INPUT);
            if (!plan_fwd)  SX_EXIT;
            fftPlanSpecs   << (specDir);
            fftPlans       << plan_fwd;
            fftPlanCounter << 1;
            iPlanFor = int(fftPlanCounter.getSize() - 1);
            createdNewPlan = true;
         }
      }


      if (dir == Reverse || dir == Both)  {

         SxString specDir = spec + "REV";

         i = iPlanRev = int(fftPlanSpecs.findPos (specDir));
         if (i >= 0)  {
            plan_rev = fftPlans(i);
            fftPlanCounter(i)++;
         }  else  {
            SX_CLOCK (Timer::FFTPlanning);
            plan_rev = fftw_plan_dft_2d (nx, ny,
                                         (fftw_complex*)inArray,
                                         (fftw_complex*)outArray,
                                         FFTW_FORWARD, // see FFT definition
                                         fftPlanMode | FFTW_PRESERVE_INPUT);
            if (!plan_rev)  SX_EXIT;
            fftPlanSpecs   << (specDir);
            fftPlans       << plan_rev;
            fftPlanCounter << 1;
            iPlanRev = int(fftPlanCounter.getSize() - 1);
            createdNewPlan = true;
         }
      }

#  elif defined (USE_ACML_FFT)
      int i;
      SxString spec = SxString(nx) + "x" + ny;

      // --- create plan for both forward and reverse FFT
      //     note that ACML uses one plan for both directions
      i = iPlanFor = iPlanRev = int(fftPlanSpecs.findPos (spec));
      if (i >= 0)  {
         plan_fwd = plan_rev = (FFT_COMPLEX*)fftPlans(i);
         fftPlanCounter(i)++;
      }  else  {
         int error = 0, createPlan = 100;  // see ACML docu
         plan_fwd = new FFT_COMPLEX [nx*ny + 3*(nx+ny)+200];
#           if defined (FFT_SINGLE_PREC)
            cfft2dx (createPlan, 1.0, true, false, ny, nx, 
                     inArray, outArray,  plan_fwd, &error);
#           else
            zfft2dx (createPlan, 1.0, true, false, ny, nx, 
                     inArray, outArray,  plan_fwd, &error);
#           endif /* FFT_SINGLE_PREC */
         if (error)  SX_EXIT;
         fftPlanSpecs   << spec;
         fftPlans       << (void *)plan_fwd;
         fftPlanCounter << 1;
         iPlanRev = iPlanFor =  int(fftPlanCounter.getSize()) - 1;
      }
#  endif  /* USE_FFTW */

   if (createdNewPlan && (SxLoopMPI::me()==0))
      SxFFT::saveWisdom ("fftwisdom.dat");

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT2d::fftReverse (int n, const void *in, void *out)
{
   SX_CHECK (in);
   SX_CHECK (out);
   SX_CHECK (in != out);  // arrays *MUST* be different
#  if defined (USE_VECLIB_FFT)
      SX_EXIT;  // yet to be implemented
      int info;
      // --- copy input array -> working array
      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
      for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
         c2dfft ((FFT_COMPLEX *)out, &mesh[1], &mesh[0], 
                 &mesh[1], &FFT_REVERSE, &info);
#     else
         z2dfft ((FFT_COMPLEX *)out, &mesh[1], &mesh[0], 
                 &mesh[1], &FFT_REVERSE, &info);
#     endif  /* FFT_SINGLE_PREC */         
#  elif defined (USE_ACML_FFT)
      SX_UNUSED (n);
      int error = 0;
      SX_FFT_REAL fac = (autoscale) ? scaleRevFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft2dx (-1, fac, true, false, mesh[1], mesh[0], 
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_rev, &error);
#     else
         zfft2dx (-1, fac, true, false, mesh[1], mesh[0],
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_rev, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
      
#  elif USE_ESSL
#     if defined (FFT_SINGLE_PREC)
         scft3 ((FFT_COMPLEX *)in, mesh[1],  
                (FFT_COMPLEX *)out,mesh[1],
                mesh[1], mesh[0],
                FFT_REVERSE, scaleRevFFT, NULL, 0);
#     else         
         dcft3 ((FFT_COMPLEX *)in, mesh[1],  
                (FFT_COMPLEX *)out,mesh[1],
                mesh[1], mesh[0],
                FFT_REVERSE, scaleRevFFT, NULL, 0);
#     endif  /* FFT_SINGLE_PREC */

#  elif defined (USE_FFTSG)
      SX_EXIT;
//      sxprintf ("TODO: FFTSG still copies vectors!!!\n");
//      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out;
//      SxVector<Complex16> inOut (2*meshSize);
//      int i;
//      for (i=0; i < n; i++)  inOut(i) = srcPtr[i];
//
//      int inzee =  1;
//      int isign = -1;
      
//      FCALL (fft) (FINT(mesh.v[2]),  FINT(mesh.v[1]),  FINT(mesh.v[0]),
//                FINT(perfMesh.v[2]), FINT(perfMesh.v[1]), FINT(perfMesh.v[0]),
//                FCOMPLEX16VEC(inOut),
//                FINT(isign), FINT(inzee), FINT(fftCacheSize));

//      for (i=0; i < n; i++)  destPtr[i]  = inOut(i+n);
//      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
         
#  else  /* FFTW */

      fftw_execute_dft ((fftw_plan)plan_rev,
                        static_cast<fftw_complex *>(const_cast<void *>(in)),
                        (fftw_complex *)out);
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
#  endif // USE_FFTW

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT2d::fftForward (const int n, const void *in, void *out)
{
   SX_CHECK (in);
   SX_CHECK (out);
   SX_CHECK (in != out);  // arrays *MUST* be different

#  if defined (USE_VECLIB_FFT)
   SX_EXIT; // not tested
      int info;
      // --- copy input array -> working array
      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
      for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
         c2dfft ((FFT_COMPLEX *)out, &mesh[1], &mesh[0], 
                 &mesh[1], &FFT_FORWARD, &info);
#     else
         z2dfft ((FFT_COMPLEX *)out, &mesh[1], &mesh[0], 
                 &mesh[1], &FFT_FORWARD, &info);
#     endif  /* FFT_SINGLE_PREC */         
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);

#  elif defined (USE_ACML_FFT)
      SX_UNUSED (n);
      int error = 0;
      SX_FFT_REAL fac = (autoscale) ? scaleForFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft2dx (1, fac, true, false, mesh[1], mesh[0], 
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_fwd, &error);
#     else
         zfft2dx (1, fac, true, false, mesh[1], mesh[0],
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_fwd, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
      


#  elif defined (USE_ESSL)
#     if defined (FFT_SINGLE_PREC)
         scft3 ((FFT_COMPLEX *)in, mesh[1],  
                (FFT_COMPLEX *)out,mesh[1],
                mesh[1], mesh[0],
                FFT_FORWARD, scaleForFFT, NULL, 0);
#     else 
         dcft3 ((FFT_COMPLEX *)in, mesh[1],  
                (FFT_COMPLEX *)out,mesh[1],
                mesh[1], mesh[0],
                FFT_FORWARD, scaleForFFT, NULL, 0);
#     endif  /* FFT_SINGLE_PREC */         

#  elif defined (USE_FFTSG)
      SX_EXIT;
//      SX_CHECK (out-in == meshSize, out-in, meshSize);
//      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
//      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out;
//      SxVector<Complex16> inOut (2*meshSize);
//      int i;
//      for (i=0; i < n; i++)  inOut(i) = srcPtr[i];
//
//      int inzee = 1;
//      int isign = 1;
      
//      FCALL (fft) (FINT(mesh.v[2]),  FINT(mesh.v[1]),  FINT(mesh.v[0]),
//                FINT(perfMesh.v[2]), FINT(perfMesh.v[1]), FINT(perfMesh.v[0]),
//                FCOMPLEX16VEC(inOut),
//                FINT(isign), FINT(inzee), FINT(fftCacheSize));

//      for (i=0; i < n; i++)  destPtr[i]  = inOut(i+n);
//      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);

#  else /* FFTW */

      fftw_execute_dft( (fftw_plan)plan_fwd,
                        (fftw_complex *)(const_cast<void *>(in)),
                        (fftw_complex *)out);
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);
#  endif  /* USE_FFTW */
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT2d::print ()
{
   SX_CHECK (inArray);
   print (inArray);
}

void SxFFT2d::print (FFT_COMPLEX *array)
{
#  ifndef USE_FFTW   
      FFT_COMPLEX c;
#  endif      
#  ifdef USE_MKL_FFT
      FFT_COMPLEX c;
#  endif
   for (int i=0; i<mesh[0]; i++)  {
      for (int j=0; j<mesh[1]; j++)  {
         ssize_t idx = getFFTIdx(i, j);
#     if defined (USE_MKL_FFT)
         c = array[idx];
         sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#     elif defined (USE_FFTW)
         sxprintf ("(%g,%g)  ", array[idx][0],
               array[idx][1]);
#     else               
         c = array[idx];
         sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#     endif               
         }
      sxprintf ("\n");
   }
}

SxFFT2d &SxFFT2d::operator= (const SxFFT2d &in)
{
   if (&in == this)  return *this;

   SxFFT::operator= (in);

   mesh[0] = in.mesh[0];
   mesh[1] = in.mesh[1];
#  if defined (USE_FFTSG)
      perfMesh[0]     = in.perfMesh[0];
      perfMesh[1]     = in.perfMesh[1];
      fftCacheSize = in.fftCacheSize;
#  endif      
   return *this; 
}

