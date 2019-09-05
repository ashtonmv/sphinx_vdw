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
#include <sxavx.h>

#if (USE_AVX == 1)
#ifndef BA
#  define BA 3
#endif
#ifndef BB
#  define BB 1
#endif
#ifndef BD
#   define BD 8
#endif

#define PASTE4(A,B,C,D) A ## _ ## B ## _ ## C ## _ ## D
#define GEMM3MKERNEL_AVX(BA,BB,BD) PASTE4(gemm3m_kernel_avx,BA,BB,BD)
#ifndef GEMM3MKERNEL
#  define GEMM3MKERNEL GEMM3MKERNEL_AVX(BA,BB,BD)
#endif
#ifndef GEMM3MKERNEL1
#  define GEMM3MKERNEL1 GEMM3MKERNEL_AVX(1,1,BD)
#endif

#endif

#if ! (defined (BA) && defined (BB) && defined (BD))
#error You must define the block sizes BA BB and BD
#endif

#ifndef CRTYPE
#define CRTYPE SxComplex16
#define CRTYPE_C double
#endif
#ifndef CRTYPE_A
#define CRTYPE_A CRTYPE
#endif
#ifndef CRTYPE_B
#define CRTYPE_B CRTYPE
#endif
#ifndef CRTYPE_C
#define CRTYPE_C CRTYPE
#endif

#define SXGEMM3M_TUNE
#define CRTYPE_D TM<CRTYPE_A,TM<CRTYPE_B,CRTYPE_C>::Res>::Res
#include <../SxGemmm.cpp>

#include <SxCLI.h>


inline double sqr(double x) { return x*x; }
inline double sqr(SxComplex16 x) { return x.re*x.re + x.im * x.im; }

void initIn (int Na, int Nd, double *A)
{
   for (int ia = 0; ia < Na; ++ia)
      for (int id = 0; id < Nd; ++id)
         A[id + Nd * ia] = 1./double(ia + id + 1);
}

void initIn (int Na, int Nd, SxComplex16 *A)
{
   for (int ia = 0; ia < Na; ++ia)
      for (int id = 0; id < Nd; ++id)
         A[id + Nd * ia] = SxComplex16(1./double(ia + id + 1), 
                                       1./double(abs(ia-id) + 1));
}

void init1 (int Na, int Nd, double *A)
{
   for (int ia = 0; ia < Na; ++ia)
      for (int id = 0; id < Nd; ++id)
         A[id + Nd * ia] = 1.;
}

void init1 (int Na, int Nd, SxComplex16 *A)
{
   for (int ia = 0; ia < Na; ++ia)
      for (int id = 0; id < Nd; ++id)
         A[id + Nd * ia] = 1.;
}

template<class T>
double diff (int Na, int Nb, int Nc, T *res, T *res2)
{
   double s = 0.;
   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         for (int ic = 0; ic < Nc; ++ic)  {
            s += sqr( res[ic + Nc * (ib + Nb * ia)] 
                     -res2[ic + Nc * (ib + Nb * ia)]);
	    //cout << ia << ' ' << ib << ' ' << ic << endl;
	    //cout << res[ic + Nc * (ib + Nb * ia)]  << "=";
	    //cout << res2[ic + Nc * (ib + Nb * ia)]  << endl;
         }
      }
   }
   return s;
}

#define STRING2(x) #x
#define STRING(x) STRING2(x)

#include <SxTimer.h>
enum GemmmTimers { Trivial, Default, Tuned } ;
SX_REGISTER_TIMERS(GemmmTimers)
{
   regTimer (Trivial, "trivial");
   regTimer (Default, "default");
   regTimer (Tuned, "tuned");
}

void setZero (CRTYPE_A *res, int Na, int Nd)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
   for (ssize_t i = 0; i < ssize_t(Na) * ssize_t(Nd); ++i)
      res[i] = 0.;
}

int main (int argc, char **argv)
{
   SxCLI cli(argc, argv);
   int Na = cli.option ("--Na|-A","number","dimension Na (usually states)")
            .toInt(96,1);
   int Nb = cli.option ("--Nb|-B","number","dimension Nb (usually atoms)")
            .toInt(96,1);
   int Nc = cli.option ("--Nc|-C","number","dimension Nc (usually projectors)")
            .toInt(10,1);
   int Nd = cli.option ("--Nd|-D","number","dimension Nd (usually G vectors)")
            .toInt(10000,1);

   bool check = cli.option ("--check", "check blocked vs. trivial algorithm")
                .toBool ();
   bool timeDefault = cli.option ("--default", "time default blocked algorithm")
                .toBool ();

   long minOps = long(1e9 * 
            cli.option ("--ops","number","min. number of operations (in 10^9))")
                       .toDouble (3));
   cli.finalize ();

   long ops = long(Na) * long(Nb) * long(Nc) * long(Nd);
   int nCycle = int(minOps / ops) + 1;
   cout << "A: " << STRING(CRTYPE_A) << endl;
   cout << "B: " << STRING(CRTYPE_B) << endl;
   cout << "C: " << STRING(CRTYPE_C) << endl;

#if (USE_AVX == 1)
   cout << "using AVX kernel: ";
#endif
   cout << "BA=" << STRING(BA)
        << " BB=" << STRING(BB)
        << " BD=" << STRING(BD) << endl;
#if (USE_AVX == 1)
   cout << STRING((GEMM3MKERNEL for kernel)) << endl;
   cout << STRING((GEMM3MKERNEL1 for cleanup)) << endl;
#endif
   cout << "Na=" << Na << " Nb=" << Nb << " Nc=" << Nc << endl;
   cout << "ops = " << ops * nCycle 
        << " = " << double(ops * nCycle) * 1e-9 << " Gops" << endl;
#ifdef USE_OPENMP
   int nThreads = omp_get_max_threads();
#else
   int nThreads = 1;
#endif
   if (nThreads > 1)
      cout << "using " << nThreads << " openMP threads" << endl;

   CRTYPE_B *B = new CRTYPE_B[Nb * Nd];
   CRTYPE_C *C = new CRTYPE_C[Nc * Nd];
   CRTYPE_D *X = new CRTYPE_D[Na * Nb * Nc];
   initIn (Na, Nb * Nc, X);
   //init1 (Na, Nb * Nc, X);

   initIn (Nb, Nd, B);
   //init1 (Nb, Nd, B);
   initIn (Nc, Nd, C);
   //init1 (Nc, Nd, C);

   CRTYPE_A *res = new CRTYPE_A[Na * Nd];
   CRTYPE_A *res2 = NULL;
   if (check)  {
      res2 = new CRTYPE_A[Na * Nd];
      SX_CLOCK (Trivial);
      gemm3m(Na, Nb, Nc, Nd, B, Nd, C, Nd, X, Nc, Nc * Nb, res2, Nd);
   }

   if (timeDefault)  {
      for (int iCycle = 0; iCycle < nCycle; ++iCycle) {
         {
            SX_CLOCK(Default);
            setZero (res, Na, Nd);
            sxpgemm3m (Na, Nb, Nc, Nd, B, C, X, Nc, Nc * Nb, res);
         }
      }
   }
   if (timeDefault && check)
      cout << "Error default kernel: " << diff(Na, Nd, 1, res, res2) << endl;

   for (int iCycle = 0; iCycle < nCycle; ++iCycle) {
      SX_CLOCK(Tuned);
      setZero (res, Na, Nd);
#ifdef USE_OPENMP
      pgemm3m_b_omp<BA,BB,BD>
#else
      pgemm3m_b<BA,BB,BD>
#endif
         (Na, Nb, Nc, Nd, B, Nd, C, Nd, X, Nc, Nc * Nb, res, Nd);
   }

   if (check)
      cout << "Error tuned kernel: " << diff(Na, Nd, 1, res, res2) << endl;

   double tDefault = GETTIME(Default)/double(ops * nCycle) * 1e9,
          tTuned   = GETTIME(Tuned) / double(ops * nCycle) * 1e9;
   cout.precision (3);
   double flop = (2. + 2. * (Nc - 1.)/Nc + 8. / Nc);
   if (timeDefault)  {
      cout << "t per op (default): " << tDefault 
           << " ns (" << flop/tDefault << " Gflop/s";
      if (nThreads > 1)
         cout << ", per thread: " << flop/tDefault/nThreads;
      cout << ")" << endl;
   }
   cout << "t per op (tuned):   " << tTuned
        << " ns (" << flop/tTuned << " Gflop/s";
   if (nThreads > 1)
      cout << ", per thread: " << flop/tTuned/nThreads;
   cout << ")" << endl;
   long data = (sizeof(*res) * Na + sizeof(*B)*Nb + sizeof(*C)*Nc) * Nd
             + sizeof(*X) * (Na * Nb * Nc);
   data *= nCycle;
   if (timeDefault)
      cout << "global data rate (default): " 
           << double(data) / GETTIME(Default) * 1e-9 << " GB/s" << endl;
   cout << "global data rate (tuned):   " << double(data)/GETTIME(Tuned)*1e-9
        << " GB/s" << endl;

   delete [] res;
   if (res2) delete [] res2;
   delete [] X;
   delete [] B;
   delete [] C;
   printTiming ();
   return 0;
}
