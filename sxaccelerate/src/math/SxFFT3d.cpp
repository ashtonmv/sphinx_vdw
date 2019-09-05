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

#include <SxFFT3d.h>
// #include <SxVector.h>
#include <SxLoopMPI.h>

SxFFT3d::SxFFT3d () : SxFFT (None, true, true)
{
   // empty
}


SxFFT3d::SxFFT3d (Directions dir_, int nx, int ny, int nz, double omega,
                  bool symmetric_)
   : SxFFT (dir_, true, symmetric_)
{
   setMesh (nx, ny, nz, omega);
}

SxFFT3d::SxFFT3d (Directions dir_, const SxVector3<Int> meshSize_, double omega,
                  bool symmetric_)
   : SxFFT (dir_, true, symmetric_)
{
   setMesh (meshSize_(0), meshSize_(1), meshSize_(2), omega);
}

SxFFT3d::SxFFT3d (const SxFFT3d &in)
   : SxFFT (in),
     mesh (in.mesh)
#    if defined (USE_FFTSG)
       ,perfMesh (in.perfMesh),
        fftCacheSize (in.fftCacheSize)
#    endif
{
   // empty
}

SxFFT3d::~SxFFT3d ()
{
   // FFTW plans as well as inArray and outArray are destroyed
   // by SxFFT::SxFFT
}



// among others, creates plans for FFTW
void SxFFT3d::setMesh (int nx, int ny, int nz, double omega)
{

   if (verbdebug)
      cout << "SxFFT3d::setMesh called with nx, ny, nz: " << nx
           << ", " << ny << ", " << nz << endl << flush;

   #if defined (USE_FFTW)
   #  if defined (USE_FFTW_PARALLEL)
//      if(SxLoopMPI::parallelFft ())
//      {
         createMPIArrays(nx, ny, nz);
         initMPIArrays(nx, ny, nz);
//      }
   #  endif // USE_FFTW_PARALLEL
   #endif // USE_FFTW

   bool createdNewPlan = false;
   meshSize = nx * ny * nz;

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
      perfMesh = SxVector3<Int> (nx, ny, nz);
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
      SxComplex16 nan = sqrt(-1.);
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

   mesh(0) = nx;
   mesh(1) = ny;
   mesh(2) = nz;
   
   checkMeshSize (nx);
   checkMeshSize (ny);
   checkMeshSize (nz);

#  if defined (USE_FFTW)

      int i;
      SxString spec = SxString(nx) + "x" + ny + "x" + nz;

      if (dir == Forward || dir == Both)  {

         SxString specDir = spec + "FWD";
         if (fftPlanMode == FFTW_ESTIMATE) specDir += "-E";

#        if defined (USE_FFTW_PARALLEL)
         // initialization of the 2D+1D=3D FFT, forward direction
//         if(SxLoopMPI::parallelFft ())
         {

            // initialization of the 2D FFT over y and z
            i = iPlanForParA = fftPlanSpecsParA.findPos (specDir);
            if (verbdebug) cout << specDir << " 2D plan ";
            if (i >= 0)  {
               plan_fwd_par_a = fftPlansParA(i);
               fftPlanCounterParA(i)++;
               if (verbdebug) cout << "loaded, ";
            }
            else
            {
               SX_CLOCK(Timer::FFTPlanning);
               plan_fwd_par_a = fftw_plan_dft_2d(ny, nz,
                                    inArray, fftwTmp1,
                                    FFTW_BACKWARD,// see FFTW definition
                                    fftPlanMode | FFTW_PRESERVE_INPUT);
               if (!plan_fwd_par_a)  SX_EXIT;
               //
               fftPlanSpecsParA   << specDir;
               fftPlansParA       << plan_fwd_par_a;
               fftPlanCounterParA << 1;
               iPlanForParA = fftPlanCounterParA.getSize() - 1;
               createdNewPlan = true;
               if (verbdebug) cout << "created, ";
            }
            if (verbdebug) cout << "plan=" << plan_fwd_par_a << endl;


            // initialization of the 1D FFT over x
            i = iPlanForParB = fftPlanSpecsParB.findPos (specDir);
            if (verbdebug) cout << specDir << " 1D plan ";
            if (i >= 0)  {
               plan_fwd_par_b = fftPlansParB(i);
               fftPlanCounterParB(i)++;
               if (verbdebug) cout << "loaded, ";
            }
            else
            {
               SX_CLOCK(Timer::FFTPlanning);
               plan_fwd_par_b = fftw_plan_dft(1, &nx,
                                    inArray, fftwTmp1,
                                    FFTW_BACKWARD,// see FFTW definition
                                    fftPlanMode | FFTW_PRESERVE_INPUT);
               if (!plan_fwd_par_b)  SX_EXIT;
               //
               fftPlanSpecsParB   << specDir;
               fftPlansParB       << plan_fwd_par_b;
               fftPlanCounterParB << 1;
               iPlanForParB = fftPlanCounterParB.getSize() - 1;
               createdNewPlan = true;
               if (verbdebug) cout << "created, ";
            }
            if (verbdebug) cout << "plan=" << plan_fwd_par_b << endl;

         }
#        endif // USE_FFTW_PARALLEL

         // initialization of the conventional 3D FFT, forward direction
         i = iPlanFor = int(fftPlanSpecs.findPos (specDir));
         if (verbdebug) cout << specDir << " 3D plan ";
         if (i >= 0)  {
            plan_fwd = fftPlans(i);
            fftPlanCounter(i)++;
            if (verbdebug) cout << "loaded, ";
         }  else  {
            SX_CLOCK(Timer::FFTPlanning);
            plan_fwd = fftw_plan_dft_3d (nx, ny, nz,
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
            if (verbdebug) cout << "created, ";
         }
         if (verbdebug) cout << "plan=" << plan_fwd << endl;

      }


      if (dir == Reverse || dir == Both)  {

         SxString specDir = spec + "REV";
         if (fftPlanMode == FFTW_ESTIMATE) specDir += "-E";

#        if defined (USE_FFTW_PARALLEL)
         // initialization of the 2D+1D=3D FFT, reverse direction
//         if(SxLoopMPI::parallelFft ())
         {
            i = iPlanRevParA = int(fftPlanSpecsParA.findPos (specDir));
            if (verbdebug) cout << specDir << " 2D plan ";
            if (i >= 0)  {
               plan_rev_par_a = fftPlansParA(i);
               fftPlanCounterParA(i)++;
               if (verbdebug) cout << "loaded, ";
            }
            else
            {
               SX_CLOCK(Timer::FFTPlanning);
               plan_rev_par_a = fftw_plan_dft_2d(ny, nz,
                                              inArray, fftwTmp1,
                                              FFTW_FORWARD,// see FFTW definition
                                              fftPlanMode | FFTW_PRESERVE_INPUT);
               if (!plan_rev_par_a)  SX_EXIT;
               fftPlanSpecsParA   << specDir;
               fftPlansParA       << plan_rev_par_a;
               fftPlanCounterParA << 1;
               iPlanRevParA = fftPlanCounterParA.getSize() - 1;
               createdNewPlan = true;
               if (verbdebug) cout << "created, ";
            }
            if (verbdebug) cout << "plan=" << plan_rev_par_a << endl;


            i = iPlanRevParB = int(fftPlanSpecsParB.findPos (specDir));
            if (verbdebug) cout << specDir << " 1D plan ";
            if (i >= 0)  {
               plan_rev_par_b = fftPlansParB(i);
               fftPlanCounterParB(i)++;
               if (verbdebug) cout << "loaded, ";
            }
            else
            {
               SX_CLOCK(Timer::FFTPlanning);
               plan_rev_par_b = fftw_plan_dft(1, &nx,
                                             inArray, fftwTmp1,
                                             FFTW_FORWARD,// see FFTW definition
                                             fftPlanMode | FFTW_PRESERVE_INPUT);
               if (!plan_rev_par_b)  SX_EXIT;
               fftPlanSpecsParB   << specDir;
               fftPlansParB       << plan_rev_par_b;
               fftPlanCounterParB << 1;
               iPlanRevParB = fftPlanCounterParB.getSize() - 1;
               createdNewPlan = true;
               if (verbdebug) cout << "created, ";
            }
            if (verbdebug) cout << "plan=" << plan_rev_par_b << endl;

         }
#        endif // USE_FFTW_PARALLEL

         i = iPlanRev = int(fftPlanSpecs.findPos (specDir));
         if (verbdebug) cout << specDir << " 3D plan ";
         if (i >= 0)  {
            plan_rev = fftPlans(i);
            fftPlanCounter(i)++;
            if (verbdebug) cout << "loaded, ";
         }  else  {
            SX_CLOCK(Timer::FFTPlanning);
            plan_rev = fftw_plan_dft_3d (nx, ny, nz,
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
            if (verbdebug) cout << "created, ";
         }
         if (verbdebug) cout << "plan=" << plan_rev << endl;

      }

#  elif defined (USE_ACML_FFT)
      int i;
      SxString spec = SxString(nx) + "x" + ny + "x" + nz;

      // --- create plan for both forward and reverse FFT
      //     note that ACML uses one plan for both directions
      i = iPlanFor = iPlanRev = static_cast<int>(fftPlanSpecs.findPos (spec));
      if (i >= 0)  {
         plan_fwd = plan_rev = (FFT_COMPLEX*)fftPlans(i);
         fftPlanCounter(i)++;
      }  else  {
         int error = 0, createPlan = 100;  // see ACML docu
         plan_fwd = new FFT_COMPLEX [nx*ny*nz + 3*(nx+ny+nz)+300];
#           if defined (FFT_SINGLE_PREC)
            cfft3dx (createPlan, 1.0, true, false, nz, ny, nx, 
                     inArray, outArray,  plan_fwd, &error);
#           else
            zfft3dx (createPlan, 1.0, true, false, nz, ny, nx, 
                     inArray, outArray,  plan_fwd, &error);
#           endif /* FFT_SINGLE_PREC */
         if (error)  SX_EXIT;
         fftPlanSpecs   << spec;
         fftPlans       << (void *)plan_fwd;
         fftPlanCounter << 1;
         iPlanRev = iPlanFor =  static_cast<int>(fftPlanCounter.getSize()) - 1;
      }
#  endif  /* USE_FFTW */

   #if defined (USE_FFTW)
   #  if defined (USE_FFTW_PARALLEL)
//      if(SxLoopMPI::parallelFft ())
//      {
         deleteMPIArrays();
//      }
   #  endif // USE_FFTW_PARALLEL
   #endif // USE_FFTW

   if (createdNewPlan && (SxLoopMPI::me()==0))
      SxFFT::saveWisdom ("fftwisdom.dat");

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT3d::fftReverse (int n, const void *in, void *out)
{
   SX_CHECK (in);
   SX_CHECK (out);
   SX_CHECK (in != out);  // arrays *MUST* be different
#  if defined (USE_VECLIB_FFT)
      int info;
      // --- copy input array -> working array
      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
      for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
         c3dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
                 &mesh.v[2], &mesh.v[1], &FFT_REVERSE, &info);
#     else
         z3dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
                 &mesh.v[2], &mesh.v[1], &FFT_REVERSE, &info);
#     endif  /* FFT_SINGLE_PREC */         
#  elif defined (USE_ACML_FFT)
      SX_UNUSED (n);
      int error = 0;
      SX_FFT_REAL fac = (autoscale) ? scaleRevFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft3dx (-1, fac, true, false, mesh.v[2], mesh.v[1], mesh.v[0], 
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_rev, &error);
#     else
         zfft3dx (-1, fac, true, false, mesh.v[2], mesh.v[1], mesh.v[0],
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_rev, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
      
#  elif USE_ESSL
#     if defined (FFT_SINGLE_PREC)
         scft3 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
                (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
                mesh.v[2], mesh.v[1], mesh.v[0],
                FFT_REVERSE, scaleRevFFT, NULL, 0);
#     else         
         dcft3 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
                (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
                mesh.v[2], mesh.v[1], mesh.v[0],
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

#     if defined (USE_FFTW_PARALLEL)
         if(SxLoopMPI::parallelFft ())
         {
            createMPIArrays (mesh(0), mesh(1), mesh(2));
            initMPIArrays (mesh(0), mesh(1), mesh(2));
            // the flag "-1" means by definition "reverse FFT"
            doParallelFFTW ((FFT_COMPLEX *)in, (FFT_COMPLEX *)out, -1);
            deleteMPIArrays ();
         }
         else
         {
#     endif // USE_FFTW_PARALLEL
         fftw_execute_dft ((fftw_plan)plan_rev,
                           static_cast<fftw_complex *>(const_cast<void *>(in)),
                           (fftw_complex *)out);
#     if defined (USE_FFTW_PARALLEL)
         }
#     endif // USE_FFTW_PARALLEL

      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleRevFFT, n);
#  endif // USE_FFTW

}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT3d::fftForward (const int n, const void *in, void *out)
{
   SX_CHECK (in);
   SX_CHECK (out);
   SX_CHECK (in != out);  // arrays *MUST* be different

#  if defined (USE_VECLIB_FFT)
      int info;
      // --- copy input array -> working array
      FFT_COMPLEX *srcPtr  = (FFT_COMPLEX *)in;
      FFT_COMPLEX *destPtr = (FFT_COMPLEX *)out; 
      for (int i=0; i < n; i++)  *destPtr++ = *srcPtr++;
#     if defined (FFT_SINGLE_PREC)
         c3dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
                 &mesh.v[2], &mesh.v[1], &FFT_FORWARD, &info);
#     else
         z3dfft ((FFT_COMPLEX *)out, &mesh.v[2], &mesh.v[1], &mesh.v[0], 
                 &mesh.v[2], &mesh.v[1], &FFT_FORWARD, &info);
#     endif  /* FFT_SINGLE_PREC */         
      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);

#  elif defined (USE_ACML_FFT)
      SX_UNUSED (n);
      int error = 0;
      SX_FFT_REAL fac = (autoscale) ? scaleForFFT : 1.0;
#     if defined (FFT_SINGLE_PREC)
         cfft3dx (1, fac, true, false, mesh.v[2], mesh.v[1], mesh.v[0], 
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_fwd, &error);
#     else
         zfft3dx (1, fac, true, false, mesh.v[2], mesh.v[1], mesh.v[0],
                  (FFT_COMPLEX *)const_cast<void *>(in),  
                  (FFT_COMPLEX *)out, plan_fwd, &error);
#     endif  /* FFT_SINGLE_PREC */         
      SX_CHECK (!error, error);
      


#  elif defined (USE_ESSL)
#     if defined (FFT_SINGLE_PREC)
         scft3 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
                (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.[v1],
                mesh.v[2], mesh.v[1], mesh.v[0],
                FFT_FORWARD, scaleForFFT, NULL, 0);
#     else 
         dcft3 ((FFT_COMPLEX *)in, mesh.v[2], mesh.v[2]*mesh.v[1],  
                (FFT_COMPLEX *)out,mesh.v[2], mesh.v[2]*mesh.v[1],
                mesh.v[2], mesh.v[1], mesh.v[0],
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

#     if defined (USE_FFTW_PARALLEL)
         if(SxLoopMPI::parallelFft ())
         {
            createMPIArrays (mesh(0), mesh(1), mesh(2));
            initMPIArrays (mesh(0), mesh(1), mesh(2));
            // the flag "-1" means by definition "reverse FFT"
            doParallelFFTW ((FFT_COMPLEX *)in, (FFT_COMPLEX *)out, +1);
            deleteMPIArrays ();
         }
         else
         {
#     endif // USE_FFTW_PARALLEL
            fftw_execute_dft( (fftw_plan)plan_fwd,
                              (fftw_complex *)(const_cast<void *>(in)),
                              (fftw_complex *)out);
#     if defined (USE_FFTW_PARALLEL)
         }
#     endif // USE_FFTW_PARALLEL

      if (autoscale)  scale ((SX_FFT_COMPLEX *)out, scaleForFFT, n);
#  endif  /* USE_FFTW */
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxFFT3d::print ()
{
   SX_CHECK (inArray);
   print (inArray);
}

void SxFFT3d::print (FFT_COMPLEX *array)
{
   int i, j, k;
#  ifndef USE_FFTW   
      FFT_COMPLEX c;
#  endif      
#  ifdef USE_MKL_FFT
      FFT_COMPLEX c;
#  endif
   for (i=0; i<mesh.v[0]; i++)  {
      for (j=0; j<mesh.v[1]; j++)  {
         for (k=0; k<mesh.v[2]; k++)  {
            ssize_t idx = mesh.getMeshIdx(i, j, k, SxMesh3D::Positive);
#           if defined (USE_MKL_FFT)
               c = array[idx];
               sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#           elif defined (USE_FFTW)
               sxprintf ("(%g,%g)  ", array[idx][0],
                                    array[idx][1]);
#           else               
               c = array[idx];
               sxprintf ("(%g,%g)  ", SX_REAL(c), SX_IMAG(c));
#           endif               
         }
         sxprintf ("\n");
      }
      sxprintf ("\n");
   }
}

SxFFT3d &SxFFT3d::operator= (const SxFFT3d &in)
{
   if (&in == this)  return *this;

   SxFFT::operator= (in);

   mesh = in.mesh;
#  if defined (USE_FFTSG)
      perfMesh     = in.perfMesh;
      fftCacheSize = in.fftCacheSize;
#  endif      
   return *this; 
}


#if defined (USE_FFTW)
#if defined (USE_FFTW_PARALLEL)

void SxFFT3d::createMPIArrays(int nx, int ny, int nz)
{
   if (verbdebug)
      cout << "SxFFT3d::createMPIArrays called with nx, ny, nz: " << nx
            << ", " << ny << ", " << nz << endl << flush;

   // allocate space for temporary arrays
   fftwTmp1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
   fftwTmp2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);

   int np = SxLoopMPI::nr();

//   // Layout from test code:
//   nzMPI = new int[np];
//   z0MPI = new int[np];
//   z1MPI = new int[np];
//
//   nKxKyMPI = new int[np];
//   KxKy0MPI = new int[np];
//   KxKy1MPI = new int[np];

   // Layout adapted to the SPHINX code:
   nxMPI = new int[np];
   x0MPI = new int[np];
   x1MPI = new int[np];

   nKyKzMPI = new int[np];
   KyKz0MPI = new int[np];
   KyKz1MPI = new int[np];

   return;
}

void SxFFT3d::deleteMPIArrays()
{
   if (verbdebug)
      cout << "SxFFT3d::deleteMPIArrays called with nx, ny, nz: " << mesh(0)
            << ", " << mesh(1) << ", " << mesh(2) << endl << flush;

   if (fftwTmp1)   fftw_free(fftwTmp1);
   if (fftwTmp2)   fftw_free(fftwTmp2);

//   if (nzMPI)      delete [] nzMPI;
//   if (z0MPI)      delete [] z0MPI;
//   if (z1MPI)      delete [] z1MPI;
//   if (nKxKyMPI)   delete [] nKxKyMPI;
//   if (KxKy0MPI)   delete [] KxKy0MPI;
//   if (KxKy1MPI)   delete [] KxKy1MPI;

   if (nxMPI)      delete [] nxMPI;
   if (x0MPI)      delete [] x0MPI;
   if (x1MPI)      delete [] x1MPI;

   if (nKyKzMPI)   delete [] nKyKzMPI;
   if (KyKz0MPI)   delete [] KyKz0MPI;
   if (KyKz1MPI)   delete [] KyKz1MPI;

   return;
}


void SxFFT3d::initMPIArrays(int nx, int ny, int nz)
{
   if (verbdebug)
      cout << "SxFFT3d::initMPIArrays called with nx, ny, nz: " << nx
            << ", " << ny << ", " << nz << endl << flush;

//   // allocate space for temporary arrays
//   fftwTmp1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
//   fftwTmp2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);

   // make sure that this routine is and can be called only once
   int np = SxLoopMPI::nr();

   // --- initialize the domain decomposition in z ---

//   nzMPI = new int[np];
//   z0MPI = new int[np];
//   z1MPI = new int[np];

   // calculate the number of z points that remain when distributing
   // these points equally to the ranks
   int rxMPI=nx%np;
   // calculate the actual number of z points each rank has to work on
   for (int ip=0; ip<np; ip++)
   {
      nxMPI[ip] = nx/np;
      if (rxMPI>0)
      {
         nxMPI[ip]++;
         rxMPI--;
      }
   }
   // calculate the start and stop indices in z for each rank
   for (int ip=0; ip<np; ip++)
   {
      if(ip>0)
      {
         x0MPI[ip] = x1MPI[ip-1]+1;
         x1MPI[ip] = x1MPI[ip-1]+nxMPI[ip];
      }
      else
      {
         x0MPI[ip] = 0;
         x1MPI[ip] = nxMPI[ip]-1;
      }
   }
//   if(SxLoopMPI::me()==0)
//   {
//      sxprintf("rank  nzMPI  z0MPI  z1MPI\n");
//      for (int ip=0; ip<np; ip++)
//      {
//         sxprintf("   %d     %d      %d     %d\n",
//                  ip, nzMPI[ip], z0MPI[ip], z1MPI[ip]);
//      }
//   }


   // --- initialize the domain decomposition in (KX x KY) ---

//   nKxKyMPI = new int[np];
//   KxKy0MPI = new int[np];
//   KxKy1MPI = new int[np];

   // total number of 1D transforms along z
   int nKyKz = ny*nz;
   // calculate the number of (KX x KY) points that remain
   // when distributing them equally to the ranks
   int rKyKzMPI = nKyKz%np;
   // calculate the actual number of (KX x KY) points each rank has to work on
   for (int ip=0; ip<np; ip++)
   {
      nKyKzMPI[ip] = nKyKz/np;
      if (rKyKzMPI>0)
      {
         nKyKzMPI[ip]++;
         rKyKzMPI--;
      }
   }
   // calculate the start and stop indices in (KX x KY) for each rank
   for (int ip=0; ip<np; ip++)
   {
      if(ip>0)
      {
         KyKz0MPI[ip] = KyKz1MPI[ip-1]+1;
         KyKz1MPI[ip] = KyKz1MPI[ip-1]+nKyKzMPI[ip];
      }
      else
      {
         KyKz0MPI[ip] = 0;
         KyKz1MPI[ip] = nKyKzMPI[ip]-1;
      }
   }
//   if(me==0)
//   {
//      sxprintf("rank  nKxKyMPI  KxKy0MPI  KxKy1MPI\n");
//      for (int ip=0; ip<np; ip++)
//      {
//         sxprintf("   %d     %d      %d     %d\n",
//               ip, nKxKyMPI[ip], KxKy0MPI[ip], KxKy1MPI[ip]);
//      }
//   }

   return;
}


void SxFFT3d::doParallelFFTW(FFT_COMPLEX *input,
                    FFT_COMPLEX *output,
                    int dir)
{
   int np = SxLoopMPI::nr ();
   int me = SxLoopMPI::me ();

   int nx = mesh(0);
   int ny = mesh(1);
   int nz = mesh(2);

   if (verbdebug)
      cout << "SxFFT3d::doParallelFFTW called for nx, ny, nz: " << nx
            << ", " << ny << ", " << nz << endl << flush;

   fftw_plan plan1;
   fftw_plan plan2;

   if(dir==-1)
   {
      plan1=plan_rev_par_a;
      plan2=plan_rev_par_b;
   }
   else if(dir==1)
   {
      plan1=plan_fwd_par_a;
      plan2=plan_fwd_par_b;
   }
   else
   {
      SX_EXIT;
   }


   // implement a parallel 3D FFT using consecutive 2D and 1D FFTs

   // (1) do 2D FFTs over x and y
   for (int ix=x0MPI[me]; ix<=x1MPI[me]; ix++)
   {
      int off=ix*(ny*nz);
      fftw_execute_dft((fftw_plan)plan1, &input[off], &fftwTmp1[off]);
      // cout << "iz=" << iz << endl << flush;
   }

   // (1a) synchronize the arrays between the ranks
   for (int ip=0; ip<np; ip++)
   {
      int nelem =ny*nz*nxMPI[ip];
      int offset=ny*nz*x0MPI[ip];
      // MPI_Bcast(&fftwTmp1[offset], nelem, MPI_DOUBLE_COMPLEX, ip, MPI_COMM_WORLD);
      SxLoopMPI::bcast( (double*)&fftwTmp1[offset], 2*nelem, ip);
   }

   // (2) rearrange the array
   {
//      int offout=0;
//      for (int ix=0; ix<nx; ix++)
//      {
//         for (int iy=0; iy<ny; iy++)
//         {
//            for (int iz=0; iz<nz; iz++)
//            {
//               int offin  = ix+iy*nx+iz*(nx*ny);
//               fftwTmp2[offout][0] = fftwTmp1[offin][0];
//               fftwTmp2[offout][1] = fftwTmp1[offin][1];
//               offout++;
//            }
//         }
//      }
      int offout=0;
      for (int iz=0; iz<nz; iz++)
      {
         for (int iy=0; iy<ny; iy++)
         {
            for (int ix=0; ix<nx; ix++)
            {
               int offin = iz+iy*nz+ix*(ny*nz);
               fftwTmp2[offout][0] = fftwTmp1[offin][0];
               fftwTmp2[offout][1] = fftwTmp1[offin][1];
               offout++;
            }
         }
      }
   }

   // (3) do a 1D FFT along x
   // for (int iKxKy=0; iKxKy<nKxKy; iKxKy++)
   for (int iKyKz=KyKz0MPI[me]; iKyKz<=KyKz1MPI[me]; iKyKz++)
   {
      int off=iKyKz*nx;
      fftw_execute_dft((fftw_plan)plan2, &fftwTmp2[off], &fftwTmp1[off]);
   }

   // (3a) synchronize the arrays between the ranks
   for (int ip=0; ip<np; ip++)
   {
      int nelem =nx*nKyKzMPI[ip];
      int offset=nx*KyKz0MPI[ip];
      // MPI_Bcast(&fftwtmp2[offset], nelem, MPI_DOUBLE_COMPLEX, ip, MPI_COMM_WORLD);
      // TODO:  implement the method with the native datatype
      SxLoopMPI::bcast( (double*)&fftwTmp1[offset], 2*nelem, ip);
   }

   // (4) rearrange again to restore the desired order
   {
      int offin=0;
//      for (int ix=0; ix<nx; ix++)
//      {
//         for (int iy=0; iy<ny; iy++)
//         {
//            for (int iz=0; iz<nz; iz++)
//            {
//               int offout = ix+iy*nx+iz*(nx*ny);
//               output[offout][0] = fftwTmp1[offin][0];
//               output[offout][1] = fftwTmp1[offin][1];
//               offin++;
//            }
//         }
//      }
      for (int iz=0; iz<nz; iz++)
      {
         for (int iy=0; iy<ny; iy++)
         {
            for (int ix=0; ix<nx; ix++)
            {
               int offout = iz+iy*nz+ix*(ny*nz);
               output[offout][0] = fftwTmp1[offin][0];
               output[offout][1] = fftwTmp1[offin][1];
               offin++;
            }
         }
      }
   }

   return;
}

#endif // USE_FFTW_PARALLEL
#endif // USE_FFTW
