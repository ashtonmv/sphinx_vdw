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
#pragma GCC optimize "3"

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

#define SXGEMMM_TUNE
#define CRTYPE_D TM<CRTYPE_A,TM<CRTYPE_B,CRTYPE_C>::Res>::Res
#include <../SxGemmm.cpp>
#if ! (defined (BA) && defined (BB) && defined (BD))
#error You must define the block sizes BA BB and BD
#endif

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
                                       0.5/double(abs(ia-id) + 1));
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
//            cout << res[ic + Nc * (ib + Nb * ia)] << ' ' << res2[ic + Nc * (ib + Nb * ia)] << endl;
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

   long ops = size_t(Na) * size_t(Nb) * size_t(Nc) * size_t(Nd);
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

   CRTYPE_A *A = new CRTYPE_A[Na * Nd];
   CRTYPE_B *B = new CRTYPE_B[Nb * Nd];
   CRTYPE_C *C = new CRTYPE_C[Nc * Nd];

   initIn (Na, Nd, A);
   initIn (Nb, Nd, B);
   initIn (Nc, Nd, C);

   CRTYPE_D *res = new CRTYPE_D[Na * Nb * Nc];
   CRTYPE_D *res2 = NULL;
   if (check)  {
      res2 = new CRTYPE_D[Na * Nb * Nc];
      SX_CLOCK (Trivial);
      gemmm(Na, Nb, Nc, Nd, A, Nd, B, Nd, C, Nd, res2, Nc, Nc * Nb);
   }

   if (timeDefault)  {
      for (ssize_t i = 0; i < ssize_t(Na) * ssize_t(Nb) * ssize_t(Nc); i++)
         res[i] = -1.;
      for (int iCycle = 0; iCycle < nCycle; ++iCycle) {
         {
            SX_CLOCK(Default);
            sxgemmm(Na, Nb, Nc, Nd, A, B, C, res, Nc, Nc * Nb);
         }
      }
   }
   if (timeDefault && check)
      cout << "Error default kernel: " << diff(Na, Nb, Nc, res, res2) << endl;


   for (ssize_t i = 0; i < ssize_t(Na) * ssize_t(Nb) * ssize_t(Nc); i++)
      res[i] = -1.;
   for (int iCycle = 0; iCycle < nCycle; ++iCycle) {
      SX_CLOCK(Tuned);
#ifdef USE_OPENMP
      gemmm_b2_omp<BA,BB,BD,CRTYPE_A,CRTYPE_B,CRTYPE_C>
#else
      gemmm_b2<BA,BB,BD,CRTYPE_A,CRTYPE_B,CRTYPE_C>
#endif
         (Na, Nb, Nc, Nd, A, Nd, B, Nd, C, Nd, res, Nc, Nc * Nb);
   }

   if (check)
      cout << "Error tuned kernel: " << diff(Na, Nb, Nc, res, res2) << endl;

   double tDefault = GETTIME(Default)/double(ops * nCycle) * 1e9,
          tTuned   = GETTIME(Tuned) / double(ops * nCycle) * 1e9;
   cout.precision (3);
   double flop = double(sizeof(*res) / sizeof(double))
            * (2. + 3. / Nc);
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
   long data = (sizeof(*A) * Na + sizeof(*B)*Nb + sizeof(*C)*Nc) * Nd
             + sizeof(*res) * (Na * Nb * Nc);
   data *= nCycle;
   if (timeDefault)
      cout << "global data rate (default): " 
           << double(data) / GETTIME(Default) * 1e-9 << " GB/s" << endl;
   cout << "global data rate (tuned):   " << double(data)/GETTIME(Tuned)*1e-9
        << " GB/s" << endl;

   long opsbf = sizeof(*A) * Na * Nd  // load A (negligible)
                + sizeof(*B) * Nb * Nd * max(Na/BA,1) // load B (negligible)
                + sizeof(*C) * max(Na/BA,1) * Nb * Nc * Nd // load C
                + 2*sizeof(*res) * max(Nd/BD,1) * Na * Nb * Nc;// load/store res
   opsbf *= nCycle;
   if (timeDefault)
      cout << "kernel data rate (default): " 
           << double(opsbf) / GETTIME(Default) * 1e-9
           << " GB/s (for all cores)" << endl;
   cout << "kernel data rate (tuned):   " 
        << double(opsbf) / GETTIME(Tuned) * 1e-9
        << " GB/s (for all cores)" << endl;

   delete [] res;
   if (res2) delete [] res2;
   delete [] A;
   delete [] B;
   delete [] C;
   printTiming ();
   return 0;
}
