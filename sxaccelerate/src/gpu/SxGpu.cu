#include <cstdio>

#ifdef STANDALONE
#include <cstdlib>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <assert.h>
#else
#include <SxComplex.h>
#include "SxGpu.h"
#endif

#include <cuda.h>
#include <cublas_v2.h>

#undef CUDA_DEBUG
//#define CUDA_DEBUG
#ifdef CUDA_DEBUG
#define CU_CHECK( CU_CALL )                                \
{                                                           \
   CU_CALL;                                                 \
   cudaDeviceSynchronize();                                 \
   cudaError_t status = cudaGetLastError();                 \
   if (status != cudaSuccess)                               \
   {                                                        \
      printf("CUDA error at %s:%d\n", __FILE__, __LINE__);  \
      /*exit(1);*/                                          \
   }                                                        \
}
#else
#define CU_CHECK( CU_CALL )                                \
{                                                           \
   CU_CALL;                                                 \
}
#endif


namespace gpu
{
   int dev = 0;
   namespace mmm3vars
   {
      bool initialized = false;
      // --- CUDA device pointers
      cuDoubleComplex *A_d, *B_d, *R_d, *BC_d;
      double *C_d;
      // --- CUDA grid parameters
      dim3 block;
      dim3 grid;
      // --- CUDA hardware parameters
      int warpsize;
      int maxThreadsPerBlock;
      // --- handle for CUBLAS ZGEMM
      cublasHandle_t cublasHandle;
      // --- factors for CUBLAS ZGEMM
      const cuDoubleComplex alpha = {1.0, 0.0};
      const cuDoubleComplex beta  = {0.0, 0.0};

      __global__
      void knl_mmm3_calcBC(int nb, int nc, int nd,
            cuDoubleComplex *B, double *C, cuDoubleComplex *BC)
      {
         int id = blockIdx.x*blockDim.x + threadIdx.x;
         int ib = blockIdx.y*blockDim.y + threadIdx.y;
         int ic = blockIdx.z*blockDim.z + threadIdx.z;

         if ( (id < nd) && (ic < nc) && (ib < nb) )
         {
            int idxb = id + ib*nd;
            int idxc = id + ic*nd;
            int idxr = id + ic*nd + ib*nd*nc;
            BC[idxr] = make_cuDoubleComplex( cuCreal(B[idxb]) * C[idxc],
                                             cuCimag(B[idxb]) * C[idxc] );
         }
      }
   } // namespace mmm3vars

   void mmm3(int Na, int Nb, int Nc, int Nd,
              const cuDoubleComplex *A, const cuDoubleComplex *B, const double *C, cuDoubleComplex *R,
              int rldc, int rldbc)
   {
      using namespace mmm3vars;
      if (!initialized)
      {
         cudaSetDevice(dev);
         cudaDeviceProp prop;
         CU_CHECK( cudaGetDeviceProperties(&prop, dev) );
         warpsize = prop.warpSize;
         maxThreadsPerBlock = prop.maxThreadsPerBlock;

         // --- set up grid configuration for BC calculation
         int blockSize = warpsize;
         block.x = blockSize;
         block.y = blockSize;
         block.z = 1;
         grid.x  = (unsigned)ceil((Nd)/(double)(block.x));
         grid.y  = (unsigned)ceil((Nb)/(double)(block.y));
         grid.z  = (unsigned)ceil((Nc)/(double)(block.z));

         printf("MMM3 CUDA launch configuration : block(%d,%d,%d), grid(%d,%d,%d)\n",
                block.x,block.y,block.z, grid.x,grid.y,grid.z);

         CU_CHECK( cudaMalloc( (void**)&A_d,     Na*Nd*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&B_d,     Nb*Nd*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&C_d,     Nc*Nd*sizeof(double)) );
         CU_CHECK( cudaMalloc( (void**)&BC_d, Nd*Nc*Nb*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&R_d,  Na*rldbc*sizeof(cuDoubleComplex)) );

         cublasCreate(&cublasHandle);

         initialized = true;
      }

      {
         CU_CHECK( cudaMemcpy( A_d, A, Na*Nd*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) );
         CU_CHECK( cudaMemcpy( B_d, B, Nb*Nd*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) );
         CU_CHECK( cudaMemcpy( C_d, C, Nc*Nd*sizeof(double), cudaMemcpyHostToDevice) );
      }

      // --- STEP 1 : calculate BC[d,c,b] = B[d,b] * C[d,c]
      {
         knl_mmm3_calcBC<<<grid,block>>>(Nb, Nc, Nd, B_d, C_d, BC_d);
         CU_CHECK();
      }

      // --- STEP 2 : calculate R[c,b,a] = sum_d BC[d,c,b] A[d,a] using zgemm
      {
         cublasZgemm       // --- documentation for input parameters, cf CUBLAS guide
         (
            cublasHandle,  // input handle to the cuBLAS library context
            CUBLAS_OP_T,   // input operation op(A) that is non- or (conj.) transpose
            CUBLAS_OP_N,   // input operation op(B) that is non- or (conj.) transpose
            rldbc,         // m -- input number of rows of matrix op(A) and C
            Na,            // n -- input number of columns of matrix op(B) and C
            Nd,            // k -- input number of columns of op(A) and rows of op(B)
            &alpha,        // scalar used for multiplication
            BC_d,          // A -- array of dimensions lda x k with lda>=max(1,m) if transa == CUBLAS_OP_N and lda x m with lda>=max(1,k) otherwise
            Nd,            // lda -- input leading dimension of two-dimensional array used to store the matrix A
            A_d,           // B -- array of dimension ldb x n with ldb>=max(1,k) if transa == CUBLAS_OP_N and ldb x k with ldb>=max(1,n) otherwise
            Nd,            // ldb -- input leading dimension of two-dimensional array used to store the matrix B
            &beta,         // scalar used for multiplication
            R_d,           // C -- array of dimensions ldc x n with ldc>=max(1,m)
            rldbc          // ldc -- input leading dimension of a two-dimensional array used to store the matrix C
         );
         CU_CHECK();
      }

      {
         CU_CHECK( cudaMemcpy( R, R_d, Na*rldbc*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost) );
      }
   } // mmm3


   namespace mm3mvars
   {
      bool initialized = false;
      // --- CUDA device pointers
      cuDoubleComplex *B_d, *BC_d, *X_d, *R_d;
      double *C_d;
      // --- CUDA grid parameters
      dim3 block;
      dim3 grid;
      // --- CUDA hardware parameters
      int warpsize;
      int maxThreadsPerBlock;
      // --- handle for CUBLAS ZGEMM
      cublasHandle_t cublasHandle;
      // --- factors for CUBLAS ZGEMM
      const cuDoubleComplex alpha = {1.0, 0.0};
      const cuDoubleComplex beta  = {0.0, 0.0};

      __global__
      void knl_mm3m_calcBC(int nb, int nc, int nd,
            cuDoubleComplex *B, double *C, cuDoubleComplex *BC)
      {
         int id = blockIdx.x*blockDim.x + threadIdx.x;
         int ib = blockIdx.y*blockDim.y + threadIdx.y;
         int ic = blockIdx.z*blockDim.z + threadIdx.z;

         if ( (id < nd) && (ic < nc) && (ib < nb) )
         {
            int idxb = id + ib*nd;
            int idxc = id + ic*nd;
            int idxr = id + ic*nd + ib*nd*nc;
            BC[idxr] = make_cuDoubleComplex( cuCreal(B[idxb]) * C[idxc],
                                             cuCimag(B[idxb]) * C[idxc] );
         }
      }
   } // namespace mm3mvars


   void mm3m(int Na, int Nb, int Nc, int Nd,
              const cuDoubleComplex *B, const double *C, const cuDoubleComplex *X, cuDoubleComplex *R,
              int xldc, int xldbc)
   {
      using namespace mm3mvars;
      if (!initialized)
      {
         cudaSetDevice(dev);
         cudaDeviceProp prop;
         CU_CHECK( cudaGetDeviceProperties(&prop, dev) );
         warpsize = prop.warpSize;
         maxThreadsPerBlock = prop.maxThreadsPerBlock;

         // --- set up grid configuration for BC calculation
         int blockSize = warpsize;
         block.x = blockSize;
         block.y = blockSize;
         block.z = 1;
         grid.x  = (unsigned)ceil((Nd)/(double)(block.x));
         grid.y  = (unsigned)ceil((Nb)/(double)(block.y));
         grid.z  = (unsigned)ceil((Nc)/(double)(block.z));

         printf("MM3M CUDA launch configuration : block(%d,%d,%d), grid(%d,%d,%d)\n",
                block.x,block.y,block.z, grid.x,grid.y,grid.z);

         CU_CHECK( cudaMalloc( (void**)&B_d,     Nb*Nd*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&C_d,     Nc*Nd*sizeof(double)) );
         CU_CHECK( cudaMalloc( (void**)&BC_d, Nd*Nc*Nb*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&X_d,  Na*xldbc*sizeof(cuDoubleComplex)) );
         CU_CHECK( cudaMalloc( (void**)&R_d,     Na*Nd*sizeof(cuDoubleComplex)) );

         cublasCreate(&cublasHandle);

         initialized = true;
      }

      {
         CU_CHECK( cudaMemcpy( B_d, B,    Nb*Nd*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) );
         CU_CHECK( cudaMemcpy( C_d, C,    Nc*Nd*sizeof(double),         cudaMemcpyHostToDevice) );
         CU_CHECK( cudaMemcpy( X_d, X, Na*xldbc*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) );
      }

      // --- STEP 1 : calculate BC[d,c,b] = B[d,b] * C[d,c]
      {
         knl_mm3m_calcBC<<<grid,block>>>(Nb, Nc, Nd, B_d, C_d, BC_d);
         CU_CHECK();
      }

      // --- STEP 2 : calculate R[d,a] = sum_{b,c} BC[d,c,b] X[c,b,a] using zgemm
      {
         cublasZgemm       // --- documentation for input parameters, cf CUBLAS guide
         (
            cublasHandle,  // input handle to the cuBLAS library context
            CUBLAS_OP_N,   // input operation op(A) that is non- or (conj.) transpose
            CUBLAS_OP_N,   // input operation op(B) that is non- or (conj.) transpose
            Nd,            // m -- input number of rows of matrix op(A) and C
            Na,            // n -- input number of columns of matrix op(B) and C
            xldbc,         // k -- input number of columns of op(A) and rows of op(B)
            &alpha,        // scalar used for multiplication
            BC_d,          // A -- array of dimensions lda x k with lda>=max(1,m) if transa == CUBLAS_OP_N and lda x m with lda>=max(1,k) otherwise
            Nd,            // lda -- input leading dimension of two-dimensional array used to store the matrix A
            X_d,           // B -- array of dimension ldb x n with ldb>=max(1,k) if transa == CUBLAS_OP_N and ldb x k with ldb>=max(1,n) otherwise
            xldbc,         // ldb -- input leading dimension of two-dimensional array used to store the matrix B
            &beta,         // scalar used for multiplication
            R_d,           // C -- array of dimensions ldc x n with ldc>=max(1,m)
            Nd             // ldc -- input leading dimension of a two-dimensional array used to store the matrix C
         );
         CU_CHECK();
      }

      {
         CU_CHECK( cudaMemcpy( R, R_d, Na*Nd*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost) );
      }
   } // mm3m

   void free()
   {
      {
         using namespace mmm3vars;
         if (initialized)
         {
            CU_CHECK( cudaFree(A_d)  );
            CU_CHECK( cudaFree(B_d)  );
            CU_CHECK( cudaFree(C_d)  );
            CU_CHECK( cudaFree(BC_d) );
            CU_CHECK( cudaFree(R_d)  );
            cublasDestroy(cublasHandle);
            initialized = false;
         }
      }
      {
         using namespace mm3mvars;
         if (initialized)
         {
            CU_CHECK( cudaFree(B_d)  );
            CU_CHECK( cudaFree(C_d)  );
            CU_CHECK( cudaFree(BC_d) );
            CU_CHECK( cudaFree(X_d)  );
            CU_CHECK( cudaFree(R_d)  );
            cublasDestroy(cublasHandle);
            initialized = false;
         }
      }
   }

} // namespace gpu





#ifdef STANDALONE

using namespace std;

typedef struct {
   double re;
   double im;
} CMPLX;

namespace cpu {
   #if defined __GNUG__
      #define RESTRICT __restrict__
   #else
      #define RESTRICT
   #endif

   inline void mmm3_d_loop(const int Nd, CMPLX * RESTRICT x, CMPLX * RESTRICT y, double * RESTRICT z, CMPLX * RESTRICT r)
   {
      CMPLX sum = {0.0, 0.0};
      for (int d=0; d<Nd; ++d)
      {
         CMPLX xy;
         // --- X*Y
         xy.re   = x[d].re * y[d].re - x[d].im * y[d].im;
         xy.im   = x[d].re * y[d].im + x[d].im * y[d].re;
         // --- XY*Z
         sum.re += xy.re * z[d];
         sum.im += xy.im * z[d];
      }
      *r = sum;
   }

   void mmm3(int Na, int Nb, int Nc, int Nd,
              CMPLX * RESTRICT x, CMPLX * RESTRICT y, double * RESTRICT z, CMPLX * RESTRICT r)
   {
      int a, idxx, b, idxy, c, idxz, idxr;

#pragma omp parallel for \
   default(none) \
   private(a,idxx,b,idxy,c,idxz,idxr) \
   shared(Na,Nb,Nc,Nd,r,x,y,z)
      for (a=0; a<Na; ++a)
      {
         idxx = a*Nd;
         for (b=0; b<Nb; ++b)
         {
            idxy = b*Nd;
            for (c=0; c<Nc; ++c)
            {
               idxz = c*Nd;
               idxr = a*Nb*Nc + b*Nc + c;
               mmm3_d_loop(Nd, &x[idxx], &y[idxy], &z[idxz], &r[idxr]);
            }
         }
      }
   }

   inline CMPLX mm3m_3mul_add(const CMPLX &B, const double &C, const CMPLX &X)
   {
      CMPLX BC;
      // --- B*C
      BC.re = B.re * C;
      BC.im = B.im * C;
      // --- BC*X
      CMPLX R;
      R.re = BC.re * X.re - BC.im * X.im;
      R.im = BC.re * X.im + BC.im * X.re;
      // ---
      return R;
   }

   // --- naive mm3m CPU implementation
   void mm3m(int Na, int Nb, int Nc, int Nd,
            CMPLX * RESTRICT B, double * RESTRICT C, CMPLX * RESTRICT X, CMPLX * RESTRICT R)
   {
      int a, b, c, d;
      int idxb, idxc, idxx, idxr;
      CMPLX BCX;

#pragma omp parallel for \
   default(none) \
   private(a,b,c,d,idxb,idxc,idxx,idxr,BCX) \
   shared(Na,Nb,Nc,Nd,B,C,X,R)
      for (a=0; a<Na; ++a)
      {
         for (d=0; d<Nd; ++d)
         {
            idxr = a*Nd + d;
            R[idxr].re = 0.0;
            R[idxr].im = 0.0;
            //idxr *= 2;
            for (c=0; c<Nc; ++c)
            {
               idxc = c*Nd + d;
               //idxc *= 2;
               for (b=0; b<Nb; ++b)
               {
                  idxb = b*Nd + d;
                  //idxb *= 2;
                  idxx = a*Nb*Nc + b*Nc + c;
                  //idxx *= 2;
                  // --- definitely not optimized for speed
                  BCX = mm3m_3mul_add( B[idxb], C[idxc], X[idxx] );
                  R[idxr].re += BCX.re;
                  R[idxr].im += BCX.im;
               }
            }
         }
      }

   }


} // namespace cpu

   class stopwatch
   {
      struct timeval start;
      struct timeval stop;
      char label[128];
      double duration;
   public:
      stopwatch() {
         strcpy(label, "(no label)");
         gettimeofday(&start, NULL);
      }
      stopwatch(const char * label_) {
         strcpy(label, label_);
         gettimeofday(&start, NULL);
      }
      ~stopwatch() {
         gettimeofday(&stop, NULL);
         double dstart, dstop;
         dstart = (double) start.tv_sec + ((double) start.tv_usec)*1.e-6;
         dstop  = (double)  stop.tv_sec + ((double)  stop.tv_usec)*1.e-6;
         duration = dstop - dstart;
         printf("%s : %g s\n", label, duration);
      }

   };

   int main (int argc, char ** argv)
   {
      double *A, *B, *C, *X, *RC, *RG;

      int Na = 32;
      int Nb = 32;
      int Nc = 8;
      int Nd = 5000;

      // --- allowed relative difference between GPU and CPU results during the tests
      const double eps = 1.e-9;

      if (argc == 1)
      {
         /* use default parameters */
      }
      else if (argc == 5)
      {
         sscanf(argv[1], "%d", &Na);
         sscanf(argv[2], "%d", &Nb);
         sscanf(argv[3], "%d", &Nc);
         sscanf(argv[4], "%d", &Nd);
      }
      else
      {
         printf("Usage: %s Na Nb Nc Nd\n", argv[0]);
         return 1;
      }

      printf("parameters : Na=%d  Nb=%d  Nc=%d  Nd=%d\n",
                           Na,    Nb,    Nc,    Nd);

      // --- align memory to 64 bit boundaries (formerly to make MIC happy)
      assert(posix_memalign((void**)&A,  64,    Na*Nd*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&B,  64,    Nb*Nd*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&C,  64,    Nc*Nd*  sizeof(double)) == 0);
      assert(posix_memalign((void**)&RC, 64, Na*Nb*Nc*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&RG, 64, Na*Nb*Nc*2*sizeof(double)) == 0);

      srand(time(NULL));
      for (int i=0; i<Na*Nd*2; ++i) A[i] = double(rand()) / double(RAND_MAX);
      for (int i=0; i<Nb*Nd*2; ++i) B[i] = double(rand()) / double(RAND_MAX);
      for (int i=0; i<Nc*Nd  ; ++i) C[i] = double(rand()) / double(RAND_MAX);

      printf("\n");


      // (A) --- run MMM3 tests

      // (A1) --- run on the CPU
      {
         for (int j=0; j<Na*Nb*Nc*2; ++j) RC[j] = 0.;
         stopwatch sw("mmm3-cpu");
         cpu::mmm3(Na, Nb, Nc, Nd, (CMPLX*)A, (CMPLX*)B, C, (CMPLX*)RC);
      }
      printf("\n");

      // (A2) --- run on the GPU
      for (int i=0; i<2; ++i)
      {
         for (int j=0; j<Na*Nb*Nc*2; ++j) RG[j] = 0.;
         stopwatch sw("mmm3-gpu");
         gpu::mmm3 (Na, Nb, Nc, Nd, (cuDoubleComplex*)A, (cuDoubleComplex*)B, C, (cuDoubleComplex*)RG, Nc, Nb*Nc);
      }

      // (A3) --- compare CPU and accelerator results
      {
         for (int a=0; a<Na; ++a)
         {
            for (int b=0; b<Nb; ++b)
            {
               for (int c=0; c<Nc; ++c)
               {
                  int idxr = a*Nb*Nc + b*Nc + c;
                  idxr *= 2;
                  for (int i=0; i<2; ++i)
                  {
                     bool ok = fabs((RG[idxr] - RC[idxr])/RC[idxr]) < eps;
                     if (!ok)
                        printf("%d %d %d %d : %f %f\n", a, b, c, idxr, RG[idxr], RC[idxr]);
                     assert(ok);
                     ++idxr;
                  }
               }
            }
         }
      }
      printf("\nMMM3 : CPU and accelerator results match!\n");

      free(A);
      free(B);
      free(C);
      free(RC);
      free(RG);
      gpu::free();

      printf("\n");
      printf("\n");



      // (B) --- run MM3M tests

      assert(posix_memalign((void**)&B,  64,    Nb*Nd*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&C,  64,    Nc*Nd*  sizeof(double)) == 0);
      assert(posix_memalign((void**)&X,  64, Na*Nb*Nc*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&RC, 64,    Na*Nd*2*sizeof(double)) == 0);
      assert(posix_memalign((void**)&RG, 64,    Na*Nd*2*sizeof(double)) == 0);

      for (int i=0;    i<Nb*Nd*2; ++i) B[i] = double(rand()) / double(RAND_MAX);
      for (int i=0;    i<Nc*Nd  ; ++i) C[i] = double(rand()) / double(RAND_MAX);
      for (int i=0; i<Na*Nb*Nc*2; ++i) X[i] = double(rand()) / double(RAND_MAX);


      // (B1) --- run contraction on the CPU
      {
         for (int j=0; j<Na*Nd*2; ++j) RC[j] = 0.;
         stopwatch sw("mm3m-cpu");
         cpu::mm3m(Na, Nb, Nc, Nd, (CMPLX*)B, C, (CMPLX*)X, (CMPLX*)RC);
      }
      printf("\n");

      // (B2) --- run on the GPU
      for (int i=0; i<2; ++i)
      {
         for (int j=0; j<Na*Nd*2; ++j) RG[j] = 0.;
         stopwatch sw("mm3m-gpu");
         gpu::mm3m (Na, Nb, Nc, Nd, (cuDoubleComplex*)B, C, (cuDoubleComplex*)X, (cuDoubleComplex*)RG, Nc, Nb*Nc);
      }

      // (B3) --- compare CPU and accelerator results
      {
         for (int a=0; a<Na; ++a)
         {
            for (int d=0; d<Nd; ++d)
            {
               int idxr = a*Nd + d;
               idxr *= 2;
               for (int i=0; i<2; ++i)
               {
                  bool ok = fabs((RG[idxr] - RC[idxr])/RC[idxr]) < eps;
                  if (!ok)
                     printf("%d %d %d : %f %f\n", a, d, idxr, RG[idxr], RC[idxr]);
                  assert(ok);
                  ++idxr;
               }
            }
         }
      }
      printf("\nMM3M : CPU and accelerator results match!\n");

      free(B);
      free(C);
      free(X);
      free(RC);
      free(RG);
      gpu::free();

      return 0;
   }
#else // STANDALONE

void sx_gpu_gemmm3 (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *A,
              const SxComplex16 *B,
              const double *C,
              SxComplex16 *res, ssize_t rldc, ssize_t rldbc)
{
   gpu::mmm3 (int(Na), int(Nb), int(Nc), int(Nd),
         (cuDoubleComplex*)A, (cuDoubleComplex*)B, C,
         (cuDoubleComplex*)res,
         int(Nc), int(Nb*Nc));
}

void sx_gpu_gemm3m (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *B,
              const double *C,
              const SxComplex16 *X,
              SxComplex16 *res, ssize_t xldc, ssize_t xldbc)
{
   gpu::mm3m (int(Na), int(Nb), int(Nc), int(Nd),
         (cuDoubleComplex*)B, C, (cuDoubleComplex*)X,
         (cuDoubleComplex*)res,
         int(Nc), int(Nb*Nc));
}

void sx_gpu_gemmm_free ()
{
   gpu::free ();
}

void sx_gpu_set_device (int id)
{
   gpu::dev = id;
}

#endif // STANDALONE
