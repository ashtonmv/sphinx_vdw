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

#include <SxFFT2d1d.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

#if defined(__GNUG__) && defined(__AVX__)
#define RESTRICT __restrict__
#include <immintrin.h>
typedef __m128d v2df;
typedef __m256d v4df;
namespace {

   // --- load a duplicated double
   inline v2df load2x(const double *in)
   {
      return _mm_loaddup_pd(in);
   }

   // --- load complex number
   inline v2df  load(const SxComplex16 *in)
   {
      return _mm_loadu_pd((const double*)in);
   }

   // --- store complex number
   inline void store(SxComplex16 *in, v2df x)
   {
      _mm_storeu_pd((double*)in, x);
   }

   // --- load complex number
   inline v2df load2(const double *in)
   {
      return _mm_loadu_pd(in);
   }
   // --- store 2 doubles
   inline void store2(double *in, v2df x)
   {
      _mm_storeu_pd(in, x);
   }

   // --- copy 4 subsequent rows of a matrix into 4 arrays
   //     and prefetch the next 4 rows into level 2 cache
   //     the input should be aligned to 64-byte cache line boundary
   void fourStreamCopy(ssize_t N,
                       const SxComplex16 * RESTRICT src,
                       ssize_t strideIn,
                       SxComplex16 * RESTRICT out1,
                       SxComplex16 * RESTRICT out2,
                       SxComplex16 * RESTRICT out3,
                       SxComplex16 * RESTRICT out4)
   {
      for (ssize_t i = 0; i < N; ++i, src += strideIn)  {
         //__builtin_prefetch(src + 4, 0, 1);
         __builtin_prefetch(src + 8, 0, 2);
         __builtin_prefetch(src + 4 * strideIn, 0, 3);
         store(out1++, load(src + 0));
         store(out2++, load(src + 1));
         store(out3++, load(src + 2));
         store(out4++, load(src + 3));
      }
   }

   // --- copy 4 arrays into 4 subsequent rows of a matrix
   //     the output should be aligned to 64-byte cache line boundary
   void fourStreamCopyOut(ssize_t N,
                          SxComplex16 * RESTRICT dest,
                          ssize_t strideOut,
                          const SxComplex16 * RESTRICT in1,
                          const SxComplex16 * RESTRICT in2,
                          const SxComplex16 * RESTRICT in3,
                          const SxComplex16 * RESTRICT in4)
   {
      for (ssize_t i = 0; i < N; ++i, dest += strideOut)  {
         __builtin_prefetch(dest + 4 * strideOut, 1, 3);
         store(dest + 0, load(in1++));
         store(dest + 1, load(in2++));
         store(dest + 2, load(in3++));
         store(dest + 3, load(in4++));
      }
   }
}
#else
namespace {
   // --- copy 4 subsequent rows of a matrix into 4 arrays
   //     and prefetch the next 4 rows into level 2 cache
   //     the input should be aligned to 64-byte cache line boundary
   void fourStreamCopy(ssize_t ,
                       const SxComplex16 *,
                       ssize_t ,
                       SxComplex16 * ,
                       SxComplex16 * ,
                       SxComplex16 * ,
                       SxComplex16 * )
   {
      SX_EXIT;
   }
   // --- copy 4 arrays into 4 subsequent rows of a matrix
   //     the output should be aligned to 64-byte cache line boundary
   void fourStreamCopyOut(ssize_t ,
                          SxComplex16 * ,
                          ssize_t ,
                          const SxComplex16 * ,
                          const SxComplex16 * ,
                          const SxComplex16 * ,
                          const SxComplex16 * )
   {
      SX_EXIT;
   }
}
#endif

SxFFT2d1d::SxFFT2d1d (int nx, int ny, int nz, int nxNonZero)
   : realMesh (nx, ny, nz), nonZero(nxNonZero), meshData(NULL)
{
   SxFFT::checkMeshSize (nx);
   SxFFT::checkMeshSize (ny);
   SxFFT::checkMeshSize (nz);

   N12p = nx * ny;
   // padd N12 to cache line boundary, but make sure that N12p * (0...nx)
   // spans full cache size (avoid early cache eviction from limited cache
   // associativity)
   while (N12p % 8 != 4) N12p++;

   // allocate array
   getMesh ();
   // clean, in order to avoid nan's
   clean ();
   int dim[3];
   dim[0] = nx;
   dim[1] = ny;
   dim[2] = nz;

   // --- plan the 1D transforms (not parallel)
#ifdef USE_OPENMP
   fftw_plan_with_nthreads (1);
#endif
   plan1F = plan(1, dim, +1);
   plan1B = plan(1, dim, -1);

   // --- plan the 2D transforms (possibly parallel)
#ifdef USE_OPENMP
   // cf. SxMathLib.cpp
   fftw_plan_with_nthreads (SxUtil::getGlobalObj ().nProcs);
#endif
   plan2F = plan(2, dim, +1);
   plan2B = plan(2, dim, -1);
#ifdef USE_FFTW
   SxFFT::saveWisdom ("fftwisdom.dat");
#endif
   freeMesh ();
}

fftw_plan SxFFT2d1d::plan (int nDim, int *dim, int dir)
{
   SX_CHECK (dir == 1 || dir == -1, dir);
   SX_CHECK (nDim == 1 || nDim == 2, nDim);
   int howMany, dist;
   int offset = 0;
   SxString spec;
   // settings depending on dimensionality
   if (nDim == 1)  {
      howMany = 4;
      dist    = dim[2];
      spec    = SxString(dim[2]) + " (4x)";
      dim+=2; // use dim[2]
      offset = howMany * dist; // out-of-place
   } else {
      howMany = 2 * nonZero;
      dist    = (int)N12p;
      spec    = SxString(dim[0]) + "x" + SxString(dim[1])
              + " (" + howMany + "x)";
   }

   // spec: Forward or Backward?
   if (dir == 1)
      spec += " FWD";
   else
      spec += " REV";

   // --- check for existing plan
   ssize_t iPlan = fftPlanSpecs.findPos (spec);
   if (iPlan >= 0)  {
      fftPlanCounter(iPlan)++;
      return (fftw_plan)fftPlans(iPlan);
   }

   // --- create plan
#ifdef USE_FFTW
   SX_CLOCK (Timer::FFTPlanning);
#ifdef NDEBUG
   int mode =  FFTW_PATIENT | FFTW_DESTROY_INPUT;// always plan patient in release!
#else
   int mode = SxFFT::fftPlanMode | FFTW_DESTROY_INPUT;
#endif
   fftw_plan newPlan
      = fftw_plan_many_dft(nDim, dim, howMany,
                           (fftw_complex*)meshData, NULL, 1, dist,
                           (fftw_complex*)meshData + offset, NULL, 1, dist,
                           dir, mode);
#else
   SX_EXIT;
#endif

   fftPlanSpecs << spec;
   fftPlans     << newPlan;
   fftPlanCounter << 1;

   return newPlan;
}

void SxFFT2d1d::allocate ()
{
   SX_CHECK (nonZero > 0);
   SX_CHECK (N12p > 0);
#ifdef HAVE_POSIX_MEMALIGN
   void *memory = NULL;
   ssize_t size = N12p * 2 * nonZero;
   // align meshData at 64-byte cache-line boundary
   if (posix_memalign(&memory, 64, size * sizeof(SxComplex16)) == 0) {
      meshData = (SxComplex16*)memory;
   } else {
      cout << SX_SEPARATOR;
      cout << "|Insufficient memory.\n";
      cout << "| Can't allocate further " << size << " elements. "
         << "(" << size * sizeof(SxComplex16) << " bytes).";
      SX_EXIT;
   }
#else
   SX_EXIT;
#endif
}

void SxFFT2d1d::clean ()
{
   if (!meshData) return;
   SX_CHECK(N12p);
   SX_CHECK(nonZero);

#if defined (__GNUG__) && defined(__AVX__)
   if (cleanCode.getSize () == 0)  {
      // --- clean up complete mesh
      v4df zero = { 0., 0., 0., 0. };
      // note: N12p is multiple of 4, so no cleanup necessary
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < 2 * nonZero * N12p; i+=8)  {
         *(v4df*)(meshData + i + 0) = zero;
         *(v4df*)(meshData + i + 2) = zero;
         *(v4df*)(meshData + i + 4) = zero;
         *(v4df*)(meshData + i + 6) = zero;
      }
      return;
   }

   // --- clean up unused parts of the mesh
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
#ifdef USE_OPENMP
      SX_CHECK(omp_get_num_threads () == cleanEntry.getSize ()/2,
               omp_get_num_threads (), cleanEntry.getSize ()/2);
	   int iThread = omp_get_thread_num ();
      ssize_t start = cleanEntry(2 * iThread);
      ssize_t end = (iThread == cleanEntry.getSize ()/2 - 1)
                  ? cleanCode.getSize ()
                  : cleanEntry(2 * iThread + 2);
      ssize_t idx = cleanEntry(2 * iThread + 1);
#else
      ssize_t start = 0;
      ssize_t end = cleanCode.getSize ();
      ssize_t idx = 0;
#endif
      v2df zero  = { 0., 0. };
      v4df zero2 = { 0., 0., 0., 0. };
      for (ssize_t iCode = start; iCode < end; iCode += 2)  {
         // [iCode] is number of elements to clean out
         // [iCode + 1] is number of elements to jump over
         ssize_t nClean = cleanCode(iCode);
         // --- uneven element zero-out
         if (idx & 1)  {
            *(v2df *)(meshData + idx) = zero;
           idx++;
           nClean--;
         }
         // --- now zero out element pairs
         v4df *dest = (v4df*)(meshData + idx);
         // advance idx
         idx += nClean + cleanCode(iCode + 1);
         // get rest for unrolled loop
         ssize_t nClean2 = nClean & 7;
         nClean -= nClean2;
         SX_CHECK(nClean % 8 == 0, nClean);
         SX_CHECK((ssize_t)(dest) % sizeof(v4df) == 0);
         // unrolled loop
         for ( ; nClean; nClean -=8, dest += 4)  {
            dest[0] = zero2;
            dest[1] = zero2;
            dest[2] = zero2;
            dest[3] = zero2;
         }
         // unrolled loop trailing pairs
         while (nClean2 & 6) { *dest++ = zero2; nClean2 -= 2; }
         // unrolled loop trailing single element
         if (nClean2 & 1) *(v2df*)(dest) = zero;
      }
   }
#else
   SX_EXIT;
#endif
}

void SxFFT2d1d::destroyPlan (fftw_plan thePlan)
{
   ssize_t i = fftPlans.findPos (thePlan);
   SX_CHECK(i >= 0);
   if (--fftPlanCounter(i) == 0)  {
      fftPlanSpecs(i) = "";
      fftw_destroy_plan ((fftw_plan)thePlan);
      fftPlans(i) = NULL;
   }
}

SxFFT2d1d::~SxFFT2d1d ()
{
   freeMesh ();
   destroyPlan (plan1F);
   destroyPlan (plan1B);
   destroyPlan (plan2F);
   destroyPlan (plan2B);
}

void SxFFT2d1d::convolute (const double *V)
{
#ifdef USE_FFTW
   fftw_execute_dft(plan2F, (fftw_complex*)meshData, (fftw_complex*)meshData);
#else
   SX_EXIT;
#endif

   ssize_t N3  = realMesh(2);
   ssize_t N12 = realMesh(0) * realMesh(1);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxComplex16 *meshes;
#ifdef USE_FFTW
      meshes  = (SxComplex16*)fftw_malloc(N3 * 4 * 2 * sizeof(SxComplex16));
#else
   SX_EXIT;
#endif

#ifdef USE_OPENMP
#pragma omp for
#endif
      for (int i12 = 0; i12 < N12; i12 += 4)  {
         // copy in lower non-zero part
         fourStreamCopy(nonZero, meshData + i12, N12p,
                        meshes         , meshes +     N3,
                        meshes + 2 * N3, meshes + 3 * N3);
#if defined(__GNUG__) && defined(__AVX__)
         // --- set middle part to zero
         v2df zero = {0.,0.};
         for (int j12 = 0; j12 < 4; ++j12)  {
            SxComplex16 *dest = meshes + j12 * N3 + nonZero;
            for (int i3 = nonZero; i3 < N3 - nonZero; i3++)
               store(dest++, zero);
         }
#else
         SX_EXIT;
#endif

         // copy in upper non-zero part
         fourStreamCopy(nonZero, meshData + i12 + N12p * nonZero, N12p,
                        meshes +     N3 - nonZero, meshes + 2 * N3 - nonZero,
                        meshes + 3 * N3 - nonZero, meshes + 4 * N3 - nonZero);

         // do the 1D FFT on the current 4 meshes
#ifdef USE_FFTW
         fftw_execute_dft (plan1F, (fftw_complex*)meshes,
                                   (fftw_complex*)meshes + 4 * N3);
#else
         SX_EXIT;
#endif

         // multiply with potential V
#if defined(__GNUG__) && defined(__AVX__)
         {
            ssize_t offset = N3 * i12;
            ssize_t iv;
            SxComplex16 *workPlace = meshes + 4 * N3;
            const double *VinWork = V + offset;
            // unrolled loop with potential prefetching
            for (iv = 0;
                 iv < N3 * 4 - 8;
                 iv+=8, workPlace += 8, VinWork += 8)
            {
               // future bit of potential into level 2 cache...
               __builtin_prefetch (VinWork + N3 * 4, 0, 2);
               // next bit of potential into level 1 cache
               __builtin_prefetch (VinWork + 16, 0, 3);
               store(workPlace+0, load(workPlace + 0) * load2x(VinWork + 0));
               store(workPlace+1, load(workPlace + 1) * load2x(VinWork + 1));
               store(workPlace+2, load(workPlace + 2) * load2x(VinWork + 2));
               store(workPlace+3, load(workPlace + 3) * load2x(VinWork + 3));
               store(workPlace+4, load(workPlace + 4) * load2x(VinWork + 4));
               store(workPlace+5, load(workPlace + 5) * load2x(VinWork + 5));
               store(workPlace+6, load(workPlace + 6) * load2x(VinWork + 6));
               store(workPlace+7, load(workPlace + 7) * load2x(VinWork + 7));
            }
            // --- cleanup
            for (; iv < N3 * 4; iv++, workPlace++, VinWork++)
               store(workPlace, load(workPlace) * load2x(VinWork));
         }
#else
         (void)V;
         SX_EXIT;
#endif

         // --- backward transform
#ifdef USE_FFTW
         fftw_execute_dft(plan1B, (fftw_complex*)meshes + 4 * N3, (fftw_complex*)meshes);
#else
         SX_EXIT;
#endif
         // copy out lower non-zero part
         fourStreamCopyOut(nonZero, meshData + i12, N12p,
                           meshes         , meshes +     N3,
                           meshes + 2 * N3, meshes + 3 * N3);
         // copy out upper non-zero part
         fourStreamCopyOut(nonZero, meshData + i12 + N12p * nonZero, N12p,
                           meshes +     N3 - nonZero,
                           meshes + 2 * N3 - nonZero,
                           meshes + 3 * N3 - nonZero,
                           meshes + 4 * N3 - nonZero);

      }
#ifdef USE_FFTW
      fftw_free(meshes);
#else
      SX_EXIT;
#endif
   }

   // 2D back transform
#ifdef USE_FFTW
   fftw_execute_dft(plan2B, (fftw_complex*)meshData, (fftw_complex*)meshData);
#else
   SX_EXIT;
#endif
}

void SxFFT2d1d::addToRho (double weight, double *rho)
{
#ifdef USE_FFTW
   fftw_execute_dft(plan2F, (fftw_complex*)meshData, (fftw_complex*)meshData);
#else
   SX_EXIT;
#endif

   ssize_t N3  = realMesh(2);
   ssize_t N12 = realMesh(0) * realMesh(1);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxComplex16 *meshes;
#ifdef USE_FFTW
      meshes  = (SxComplex16*)fftw_malloc(N3 * 4 * 2 * sizeof(SxComplex16));
#else
   SX_EXIT;
#endif

#ifdef USE_OPENMP
#pragma omp for
#endif
      for (int i12 = 0; i12 < N12; i12 += 4)  {
         // copy in lower non-zero part
         fourStreamCopy(nonZero, meshData + i12, N12p,
                        meshes         , meshes +     N3,
                        meshes + 2 * N3, meshes + 3 * N3);
#if defined(__GNUG__) && defined(__AVX__)
         // --- set middle part to zero
         v2df zero = {0.,0.};
         for (int j12 = 0; j12 < 4; ++j12)  {
            SxComplex16 *dest = meshes + j12 * N3 + nonZero;
            for (int i3 = nonZero; i3 < N3 - nonZero; i3++)
               store(dest++, zero);
         }
#else
         SX_EXIT;
#endif

         // copy in upper non-zero part
         fourStreamCopy(nonZero, meshData + i12 + N12p * nonZero, N12p,
                        meshes +     N3 - nonZero, meshes + 2 * N3 - nonZero,
                        meshes + 3 * N3 - nonZero, meshes + 4 * N3 - nonZero);

         // do the 1D FFT on the current 4 meshes
#ifdef USE_FFTW
         fftw_execute_dft (plan1F, (fftw_complex*)meshes,
                                   (fftw_complex*)meshes + 4 * N3);
#else
         SX_EXIT;
#endif

         // --- take absSqr and add to rho with weight
#if defined(__GNUG__) && defined(__AVX__)
         {
            ssize_t offset = N3 * i12;
            ssize_t iv;
            SxComplex16 *workPlace = meshes + 4 * N3;
            double *rhoInWork = rho + offset;
            ssize_t ivMax = min(4L, N12 - i12) * N3;
            v2df w2 = load2x(&weight);
            // unrolled loop with rho prefetching
            for (iv = 0;
                 iv < ivMax - 8;
                 iv+=8, workPlace += 8, rhoInWork += 8)
            {
               // future bit of rho into level 2 cache...
               __builtin_prefetch (rhoInWork + N3 * 4, 0, 2);
               // next bit of rho into level 1 cache
               __builtin_prefetch (rhoInWork + 16, 0, 3);
               v2df psi0 = load(workPlace + 0);
               v2df psi1 = load(workPlace + 1);
               v2df psi2 = load(workPlace + 2);
               v2df psi3 = load(workPlace + 3);
               v2df psi4 = load(workPlace + 4);
               v2df psi5 = load(workPlace + 5);
               v2df psi6 = load(workPlace + 6);
               v2df psi7 = load(workPlace + 7);
               store2(rhoInWork + 0, load2(rhoInWork + 0)
                                 + w2 * _mm_hadd_pd(psi0 * psi0, psi1 * psi1));
               store2(rhoInWork + 2, load2(rhoInWork + 2)
                                 + w2 * _mm_hadd_pd(psi2 * psi2, psi3 * psi3));
               store2(rhoInWork + 4, load2(rhoInWork + 4)
                                 + w2 * _mm_hadd_pd(psi4 * psi4, psi5 * psi5));
               store2(rhoInWork + 6, load2(rhoInWork + 6)
                                 + w2 * _mm_hadd_pd(psi6 * psi6, psi7 * psi7));

            }
            // --- cleanup
            for (; iv < ivMax; iv++, workPlace++, rhoInWork++)
               rhoInWork[0] += weight * workPlace[0].absSqr ();
         }
#else
         (void)rho;
         SX_EXIT;
#endif

      }
#ifdef USE_FFTW
      fftw_free(meshes);
#else
      SX_EXIT;
#endif
   }

}

SxVector<TPrecFFTIdx>
SxFFT2d1d::getN231 (const SxVector<TPrecFFTIdx> &n123)
{
   SxVector<TPrecFFTIdx> res(n123.getSize ());
   ssize_t meshSize = 2 * N12p * nonZero;
   // --- setup data for cleaning unused mesh entries
   SxArray<char> used(meshSize);
   used.set (0);
   for (ssize_t i123 = 0; i123 < n123.getSize (); ++i123)  {
      SxVector3<Int> vec = realMesh.getMeshVec (n123(i123),
                                                SxMesh3D::Positive);
      SX_CHECK (   vec(2) < nonZero
                || realMesh(2) - vec(2) <= nonZero,
                vec(2), realMesh(2) - vec(2), nonZero);
      if (vec(2) * 2 > realMesh(2))  {
         vec(2) += 2 * nonZero - realMesh(2);
      }
      res(i123) = vec(1) + realMesh(1) * vec(0)
                + (TPrecFFTIdx::Type)N12p * vec(2);
      used(res(i123)) = 1;
   }
   // --- setup the cleanCode: stripes of data to be zeroed
   //     and stripes of data to keep
   SxStack<ssize_t> sizes;
   ssize_t n = 0;
   char useVal = 0;
   ssize_t nZero = 0;
   for (int i = 0; i < meshSize; i++)  {
      if (used(i) == useVal) {
         n++;
      } else {
         sizes << n;
         nZero += (1-useVal) * n;
         n = 1;
         useVal = used(i);
      }
   }
   sizes << n;
   // add final stripe of data to keep
   if (useVal == 0) sizes << 0;
   cleanCode = sizes;

#ifdef USE_OPENMP
   // --- setup entry points for openMP
   int nOmp = omp_get_max_threads ();
   ssize_t entryPoint = 0;
   ssize_t targetLoad = nZero / nOmp;
   cleanEntry.resize (2 * nOmp);
   cleanEntry(0) = 0;
   cleanEntry(1) = 0;
   ssize_t idx0 = 0;
   //cout << "nClean = " << nZero << " / " << meshSize << endl;
   //cout << "nCleanStripes = " << cleanCode.getSize () / 2 << endl;
   for (int iOmp = 1; iOmp < nOmp; ++iOmp)  {
      ssize_t threadNow = 0;
      for ( ; entryPoint < cleanCode.getSize (); entryPoint += 2)  {
         if (threadNow >= targetLoad) break;
         threadNow += cleanCode(entryPoint);
         idx0 += cleanCode(entryPoint) + cleanCode(entryPoint + 1);
      }
      cleanEntry(2 * iOmp) = entryPoint;
      cleanEntry(2 * iOmp + 1) = idx0;
   }
#endif
   return res;
}

void SxFFT2d1d::fftForward (SxComplex16 *out)
{
#ifdef USE_FFTW
   fftw_execute_dft(plan2F, (fftw_complex*)meshData, (fftw_complex*)meshData);
#else
   SX_EXIT;
#endif

   ssize_t N3  = realMesh(2);
   ssize_t N12 = realMesh(0) * realMesh(1);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxComplex16 *meshes;
#ifdef USE_FFTW
      meshes  = (SxComplex16*)fftw_malloc(N3 * 4 * 2 * sizeof(SxComplex16));
#else
   SX_EXIT;
#endif

#ifdef USE_OPENMP
#pragma omp for
#endif
      for (int i12 = 0; i12 < N12; i12 += 4)  {
         // copy in lower non-zero part
         fourStreamCopy(nonZero, meshData + i12, N12p,
                        meshes         , meshes +     N3,
                        meshes + 2 * N3, meshes + 3 * N3);
#if defined(__GNUG__) && defined(__AVX__)
         // --- set middle part to zero
         v2df zero = {0.,0.};
         for (int j12 = 0; j12 < 4; ++j12)  {
            SxComplex16 *dest = meshes + j12 * N3 + nonZero;
            for (int i3 = nonZero; i3 < N3 - nonZero; i3++)
               store(dest++, zero);
         }
#else
         SX_EXIT;
#endif

         // copy in upper non-zero part
         fourStreamCopy(nonZero, meshData + i12 + N12p * nonZero, N12p,
                        meshes +     N3 - nonZero, meshes + 2 * N3 - nonZero,
                        meshes + 3 * N3 - nonZero, meshes + 4 * N3 - nonZero);

         // do the 1D FFT on the current 4 meshes
#ifdef USE_FFTW
         fftw_execute_dft (plan1F, (fftw_complex*)meshes,
                                   (fftw_complex*)meshes + 4 * N3);
#else
         SX_EXIT;
#endif

         // copy result to "out"
#if defined(__GNUG__) && defined(__AVX__)
         {
            ssize_t offset = N3 * i12;
            ssize_t iv;
            SxComplex16 *workPlace = meshes + 4 * N3;
            SxComplex16 *outWork = out + offset;
            // unrolled loop with potential prefetching
            for (iv = 0;
                 iv < N3 * 4 - 8;
                 iv+=8, workPlace += 8, outWork += 8)
            {
               // future bit of potential into level 2 cache...
               __builtin_prefetch (outWork + N3 * 4, 0, 2);
               // next bit of potential into level 1 cache
               __builtin_prefetch (outWork + 16, 0, 3);
               store(outWork+0, load(workPlace + 0));
               store(outWork+1, load(workPlace + 1));
               store(outWork+2, load(workPlace + 2));
               store(outWork+3, load(workPlace + 3));
               store(outWork+4, load(workPlace + 4));
               store(outWork+5, load(workPlace + 5));
               store(outWork+6, load(workPlace + 6));
               store(outWork+7, load(workPlace + 7));
            }
            // --- cleanup
            for (; iv < N3 * 4; iv++, workPlace++, outWork++)
               store(outWork, load(workPlace));
         }
#else
         (void)out;
         SX_EXIT;
#endif

      }
#ifdef USE_FFTW
      fftw_free(meshes);
#else
      SX_EXIT;
#endif
   }
}

void SxFFT2d1d::fftBackward (const SxComplex16 *in)
{
   ssize_t N3  = realMesh(2);
   ssize_t N12 = realMesh(0) * realMesh(1);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxComplex16 *meshes;
#ifdef USE_FFTW
      meshes  = (SxComplex16*)fftw_malloc(N3 * 4 * 2 * sizeof(SxComplex16));
#else
   SX_EXIT;
#endif

#ifdef USE_OPENMP
#pragma omp for
#endif
      for (int i12 = 0; i12 < N12; i12 += 4)  {
         // copy data from "in"
#if defined(__GNUG__) && defined(__AVX__)
         {
            ssize_t offset = N3 * i12;
            ssize_t iv;
            SxComplex16 *workPlace = meshes + 4 * N3;
            const SxComplex16 *inWork = in + offset;
            // unrolled loop with potential prefetching
            for (iv = 0;
                 iv < N3 * 4 - 8;
                 iv+=8, workPlace += 8, inWork += 8)
            {
               // future bit of potential into level 2 cache...
               __builtin_prefetch (inWork + N3 * 4, 0, 2);
               // next bit of potential into level 1 cache
               __builtin_prefetch (inWork + 16, 0, 3);
               store(workPlace+0, load(inWork + 0));
               store(workPlace+1, load(inWork + 1));
               store(workPlace+2, load(inWork + 2));
               store(workPlace+3, load(inWork + 3));
               store(workPlace+4, load(inWork + 4));
               store(workPlace+5, load(inWork + 5));
               store(workPlace+6, load(inWork + 6));
               store(workPlace+7, load(inWork + 7));
            }
            // --- cleanup
            for (; iv < N3 * 4; iv++, workPlace++, inWork++)
               store(workPlace, load(inWork));
         }
#else
         (void)out;
         SX_EXIT;
#endif

         // --- backward transform
#ifdef USE_FFTW
         fftw_execute_dft(plan1B, (fftw_complex*)meshes + 4 * N3, (fftw_complex*)meshes);
#else
         SX_EXIT;
#endif
         // copy out lower non-zero part
         fourStreamCopyOut(nonZero, meshData + i12, N12p,
                           meshes         , meshes +     N3,
                           meshes + 2 * N3, meshes + 3 * N3);
         // copy out upper non-zero part
         fourStreamCopyOut(nonZero, meshData + i12 + N12p * nonZero, N12p,
                           meshes +     N3 - nonZero,
                           meshes + 2 * N3 - nonZero,
                           meshes + 3 * N3 - nonZero,
                           meshes + 4 * N3 - nonZero);

      }
#ifdef USE_FFTW
      fftw_free(meshes);
#else
      SX_EXIT;
#endif
   }

   // 2D back transform
#ifdef USE_FFTW
   fftw_execute_dft(plan2B, (fftw_complex*)meshData, (fftw_complex*)meshData);
#else
   SX_EXIT;
#endif
}


