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
#include <SxGemmm.h>

#ifndef WIN32


#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <sxavx.h>
#if !defined __GNUG__
//#error The sxgemmm codes needs special compiler flags, which are only provided for gcc at present.
// for other compilers, please
// - provide the optimization flags
// - modify __restrict__ attribute, if necessary
// - adjust AVX intrinsics
void sxgemmm (ssize_t, ssize_t, ssize_t, ssize_t,
              const SxComplex16 *, ssize_t,
              const SxComplex16 *,
              const double      *,
              SxComplex16 *, ssize_t, ssize_t)
{
   cout << "The sxgemmm codes needs special compiler flags, which are only provided for gcc at present." << endl;
   SX_EXIT;
}
void sxpgemm3m(ssize_t, ssize_t, ssize_t, ssize_t,
              const SxComplex16*,
              const double     *,
              const SxComplex16*, ssize_t, ssize_t,
                    SxComplex16*)
{
   cout << "The sxgemmm codes needs special compiler flags, which are only provided for gcc at present." << endl;
   SX_EXIT;
}
#else
#pragma GCC optimize "3"
#include <immintrin.h>
// note: if USE_AVX=1, we require TA=TB=SxComplex16, TC=double

namespace {
// --- result type mapper
template <class T1, class T2> class TM;
template<> class TM<double, double> { public: typedef double Res; };
template<> class TM<double, SxComplex16> { public: typedef SxComplex16 Res;};
template<> class TM<SxComplex16, double> { public: typedef SxComplex16 Res;};
template<> class TM<SxComplex16, SxComplex16> { public: typedef SxComplex16 Res;};
}

template<>
inline void SxComplex<double>::operator+= (const SxComplex<double> &in)
{
   re += in.re;
   im += in.im;
}

template<>
inline SxComplex<double> SxComplex<double>::operator* (double x) const
{
   return SxComplex<double> (re * x, im * x);
}

// we hide all internals in an anonymous namespace
// all symbols should be compiled away anyway (inline + templates)
namespace {
// --- adds a contribution, used for cleanup
template<class TA, class TB, class TC>
inline void pgemmm (int Na, int Nb, int Nc, int Nd,
           const TA *__restrict__ A, int lda,
           const TB *__restrict__ B, int ldb,
           const TC *__restrict__ C, int ldc,
           typename TM<TA,typename TM<TB, TC>::Res>::Res *__restrict__ res,
           int rldc, int rldbc)
{
   if (Na == 0 || Nb == 0 || Nc == 0 || Nd == 0) return;
   const int bd = 64;
   typename TM<TA,TB>::Res ab[bd];
   if (Nd <= bd)  {
      for (int ia = 0; ia < Na; ++ia)  {
         for (int ib = 0; ib < Nb; ++ib)  {
            // compute A * B
            for (int id = 0; id < Nd; id++)  {
               ab[id] = A[id + lda * ia]
                      * B[id + ldb * ib];
            }
            // inner loop: sum over d
            for (int ic = 0; ic < Nc; ++ic)  {
               typename TM<TA,typename TM<TB, TC>::Res>::Res S = 0.;
               for (int id = 0; id < Nd; ++id)  {
                  S += ab[id] * C[id + ldc * ic];
               }
               res[ic + rldc * ib + rldbc * ia] += S;
            }
         }
      }
      return;
   }
   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         // block over id
         for (int idb = 0; idb < Nd; idb += bd)  {
            int Ndd = (idb + bd > Nd) ? (Nd - idb) : bd;
            // compute A * B
            for (int id = 0; id < Ndd; id++)  {
               ab[id] = A[idb + id + lda * ia]
                      * B[idb + id + ldb * ib];
            }
            // inner loop: sum over d
            for (int ic = 0; ic < Nc; ++ic)  {
               typename TM<TA,typename TM<TB, TC>::Res>::Res S = 0.;
               for (int id = 0; id < Ndd; ++id)  {
                  S += ab[id] * C[idb + id + ldc * ic];
               }
               res[ic + rldc * ib + rldbc * ia] += S;
            }
         }
      }
   }
}
template<>
inline void pgemmm (int Na, int Nb, int Nc, int Nd,
           const SxComplex16 *__restrict__ A, int lda,
           const SxComplex16 *__restrict__ B, int ldb,
           const double *__restrict__ C, int ldc,
           SxComplex16 *__restrict__ res,
           int rldc, int rldbc)
{
   if (Na == 0 || Nb == 0 || Nc == 0 || Nd == 0) return;
   const int bd = 64;
   SxComplex16 ab[bd];
   int Nc2 = (Nc / 3)*3;
   if (Nd <= bd)  {
      for (int ia = 0; ia < Na; ++ia)  {
         for (int ib = 0; ib < Nb; ++ib)  {
            // compute A * B
            for (int id = 0; id < Nd; id++)  {
               ab[id].re = A[id + lda * ia].re * B[id + ldb * ib].re
                         - A[id + lda * ia].im * B[id + ldb * ib].im;
               ab[id].im = A[id + lda * ia].re * B[id + ldb * ib].im
                         + A[id + lda * ia].im * B[id + ldb * ib].re;
            }
            // inner loop: sum over d
            for (int ic = 0; ic < Nc2; ic+=3)  {
               double Sre = 0., Sim = 0.;
               double Tre = 0., Tim = 0.;
               double Ure = 0., Uim = 0.;
               for (int id = 0; id < Nd; id++)  {
                  Sre += ab[id].re * C[id + ldc * ic];
                  Sim += ab[id].im * C[id + ldc * ic];
                  Tre += ab[id].re * C[id + ldc * (ic+1)];
                  Tim += ab[id].im * C[id + ldc * (ic+1)];
                  Ure += ab[id].re * C[id + ldc * (ic+2)];
                  Uim += ab[id].im * C[id + ldc * (ic+2)];
               }
               res[ic + rldc * ib + rldbc * ia].re += Sre;
               res[ic + rldc * ib + rldbc * ia].im += Sim;
               res[ic+1 + rldc * ib + rldbc * ia].re += Tre;
               res[ic+1 + rldc * ib + rldbc * ia].im += Tim;
               res[ic+2 + rldc * ib + rldbc * ia].re += Ure;
               res[ic+2 + rldc * ib + rldbc * ia].im += Uim;
            }
            for (int ic = Nc2; ic < Nc; ic++)  {
               double Sre = 0., Sim = 0.;
               for (int id = 0; id < Nd; id++)  {
                  Sre += ab[id].re * C[id + ldc * ic];
                  Sim += ab[id].im * C[id + ldc * ic];
               }
               res[ic + rldc * ib + rldbc * ia].re += Sre;
               res[ic + rldc * ib + rldbc * ia].im += Sim;
            }
         }
      }
      return;
   }
   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         // block over id
         for (int idb = 0; idb < Nd; idb += bd)  {
            int Ndd = (idb + bd > Nd) ? (Nd - idb) : bd;
            // compute A * B
            for (int id = 0; id < Ndd; id++)  {
               ab[id] = A[idb + id + lda * ia]
                      * B[idb + id + ldb * ib];
            }
            // inner loop: sum over d
            for (int ic = 0; ic < Nc; ++ic)  {
               SxComplex16 S = 0.;
               for (int id = 0; id < Ndd; ++id)  {
                  S += ab[id] * C[idb + id + ldc * ic];
               }
               res[ic + rldc * ib + rldbc * ia] += S;
            }
         }
      }
   }
}

// --- copy data into contiguous array
template<int ba, int bd, class TA, class TB>
inline void copy(TA *__restrict__ A, int lda, TB *__restrict__ res)
{
   for (int ia = 0; ia < ba; ++ia)
      for (int id = 0; id < bd; ++id)  {
         res[id + bd * ia] = A[id + lda * ia];
         _mm_prefetch(A + id + lda * ia + bd, _MM_HINT_T2); // L3 cache
      }
}

template<int ba, int bb, int bd, class TA, class TB, class TC>
inline
void gemmm_b_kernel2 (int Nc,
                     const TA *__restrict__ Ab,
                     const TB *__restrict__ Bb,
                     const TC *__restrict__ C, int ldc,
      typename TM<TA, typename TM<TB,TC>::Res>::Res *__restrict__ res,
                     int rldc, int rldbc)
{
   typename TM<TA,TB>::Res AB[ba * bb * bd];
   // --- calculate AB[ab] = A[a] * B[b] for d-block
   for (int ia = 0; ia < ba; ++ia)
      for (int ib = 0; ib < bb; ++ib)
         for (int id = 0; id < bd; ++id)
            AB[id + bd * (ib + bb * ia)] = Ab[id + bd * ia]
                                         * Bb[id + bd * ib];
   // --- multiply AB[ab] with C[c] for d-block, and sum into output
   for (int ic = 0; ic < Nc; ++ic)  {
      for (int ia = 0, iab = 0; ia < ba; ++ia) {
         for (int ib = 0; ib < bb; ++ib, ++iab)  {
            typename TM<TA, typename TM<TB,TC>::Res>::Res S = 0.;
            for (int id = 0; id < bd; ++id)
               S += AB[id + bd * iab] * C[id + ldc * ic];
            res[ic + rldc * ib + rldbc * ia] += S;
         }
      }
   }
}

// adds a contribution, used for cleanup
template<class TB, class TC, class TD>
inline void pgemm3m(int Na, int Nb, int Nc, int Nd,
            const TB *__restrict__ B, int ldb,
            const TC *__restrict__ C, int ldc,
            const TD *__restrict__ X, int xldc, int xldbc,
                   typename TM<TD, typename TM<TB,TC>::Res>::Res *__restrict__ res, int lda)
{
   if (Na == 0 || Nb == 0 || Nc == 0 || Nd == 0) return;
   for (int a = 0; a < Na; ++a)  {
      for (int b = 0; b < Nb; ++b)  {
         for (int d = 0; d < Nd; ++d)  {
            typename TM<TD,TC>::Res S = 0.;
            for (int c = 0; c < Nc; ++c)  {
               S += C[d + ldc * c] * X[c + xldc * b + xldbc * a];
            }
            res[d + lda * a] += B[d + ldb * b] * S;
         }
      }
   }
}


template<int ba, int bb, int bd, class TB, class TC, class TD>
inline void gemm3m_kernel(int Nc,
            const TB *__restrict__ B, int ldb,
            const TC *__restrict__ C, int ldc,
            const TD *__restrict__ X, int xldc, int xldbc,
                   typename TM<TD, typename TM<TB,TC>::Res>::Res *__restrict__ res, int lda)
{
   for (int a = 0; a < ba; ++a)  {
      for (int d = 0; d < bd; ++d)  {
         for (int b = 0; b < bb; ++b)  {
            typename TM<TD,TC>::Res S = 0.;
            for (int c = 0; c < Nc; ++c)  {
               S += C[d + ldc * c] * X[c + xldc * b + xldbc * a];
            }
            res[d + lda * a] += B[d + ldb * b] * S;
         }
      }
   }
}

#if (USE_AVX == 1)
typedef double v2df __attribute__((vector_size(16)));
typedef double v4df __attribute__((vector_size(32)));
typedef long unsigned v4li __attribute__((vector_size(32)));

// --- load 2 real values and duplicate into a 256bit register
//     sort order: {c0, c1, c0, c1}
inline v4df loadcdup (const double *C)
{
   return _mm256_broadcast_pd ((const v2df*)C);
}

/// resMem += SxComplex16(re1+re2,im1+im2)
inline void operator+= (SxComplex16 &__restrict__ resMem, // result in memory
                        v4df                      resReg) // result in register
{
   // resReg sortOrder: {re1, re2, im1, im2}

   // get imaginary parts
   v2df resIm = _mm256_extractf128_pd(resReg,1);
   // get real parts
   v2df resRe = _mm256_castpd256_pd128(resReg);
   // sum real and imaginary parts
   v2df reIm = _mm_hadd_pd(resRe, resIm);
   // get current res value
   v2df old = _mm_loadu_pd((double*)(&resMem));
   // add new part and store in res
   _mm_storeu_pd((double*)(&resMem), old + reIm);
}

inline void hadd_store(v2df *mem, const v4df &resReg)
{
   // get imaginary parts
   v2df resIm = _mm256_extractf128_pd(resReg,1);
   // get real parts
   v2df resRe = _mm256_castpd256_pd128(resReg);
   // sum real and imaginary parts
   *mem = _mm_hadd_pd (resRe, resIm);
}

// --- reorder data into contiguous array
template<int ba, int bd>
inline void reorder_ri(const SxComplex16 *__restrict__ A, int lda, SxComplex16 *__restrict__ res)
{
   for (int ia = 0; ia < ba; ++ia)
      for (int id = 0; id < bd; id+=4)  {
         //res[id + bd * ia] = A[id + lda * ia];
         v4df a0 = loadcdup((const double*)(A + id + lda * ia));
         v4df a1 = loadcdup((const double*)(A + id + lda * ia + 1));
         v4df a2 = loadcdup((const double*)(A + id + lda * ia + 2));
         v4df a3 = loadcdup((const double*)(A + id + lda * ia + 3));

         // reorder data: from (re1, i1, re2, i2) to (re1, re2, i1, i2)
         v4df a01ri = _mm256_shuffle_pd(a0, a1, 0xc);
         v4df a23ri = _mm256_shuffle_pd(a2, a3, 0xc);
         // save
         _mm256_storeu_pd((double*)(res + id + bd * ia), a01ri);
         _mm256_storeu_pd((double*)(res + id + bd * ia + 2), a23ri);
         // prefetch next bd block into level 3 cache
         //_mm_prefetch(A + id + lda * ia + bd, _MM_HINT_T2);
      }
}
// --- reorder data into contiguous array
// output (re1, re2, re3, re4, im1, im2, im3, im4)
template<int ba, int bd>
inline void reorder4_ri(const SxComplex16 *__restrict__ A, int lda, SxComplex16 *__restrict__ res)
{
   for (int ia = 0; ia < ba; ++ia)  {
      for (int id = 0; id < bd; id+=4)  {
         v2df a0 = _mm_loadu_pd((const double*)(A + id + lda * ia));
         v2df a2 = _mm_loadu_pd((const double*)(A + id + lda * ia + 2));
         v2df a1 = _mm_loadu_pd((const double*)(A + id + lda * ia + 1));
         v2df a3 = _mm_loadu_pd((const double*)(A + id + lda * ia + 3));

         v4df a02 = _mm256_insertf128_pd(_mm256_castpd128_pd256(a0), a2, 1);
         v4df a13 = _mm256_insertf128_pd(_mm256_castpd128_pd256(a1), a3, 1);
         v4df a0123r = _mm256_shuffle_pd(a02, a13, 0x0);
         v4df a0123i = _mm256_shuffle_pd(a02, a13, 0xf);
         _mm256_storeu_pd((double*)(res + id + bd * ia), a0123r);
         _mm256_storeu_pd((double*)(res + id + bd * ia + 2), a0123i);

         // prefetch next bd block into level 3 cache
         _mm_prefetch(A + id + lda * ia + bd, _MM_HINT_T2);
      }
   }
}



// aux function: get 2 consecutive AbPairs
//     AbPair: multiply 2 pairs of complex numbers and
//             put them into one 256 bit register
//             sort order: {re1, re2, im1, im2}
template<int bd>
inline void getAbPairs2(const SxComplex16 *__restrict__ Ab,
                        v4df &Breim, v4df &Bimre, v4df *AB)
{
   // load A0 & A1
   v4df A0 = _mm256_loadu_pd((const double*)Ab);
   // re * re, im * im
   v4df AB0rrii =  Breim * A0;
   v4df A1 = _mm256_loadu_pd((const double*)(Ab + bd));
   // re * re, im * im
   v4df AB1rrii =  Breim * A1;
   // re * im, im * re
   v4df AB0riri = Bimre * A0;
   // re * im, im * re
   v4df AB1riri = Bimre * A1;

   v4df AB01r = _mm256_hsub_pd (AB0rrii, AB1rrii);
   v4df AB01i = _mm256_hadd_pd (AB0riri, AB1riri);
   v4df AB01d0 = _mm256_permute2f128_pd(AB01r,AB01i,0x20);
   v4df AB01d1 = _mm256_permute2f128_pd(AB01r,AB01i,0x31);
   AB[0] = _mm256_shuffle_pd(AB01d0, AB01d1, 0x0);
   AB[1] = _mm256_shuffle_pd(AB01d0, AB01d1, 0xf);
}

// aux function: get 2 consecutive AbPairs
//     AbPair: multiply 2 pairs of complex numbers and
//             put them into one 256 bit register
//             sort order: {re1, re2, im1, im2}
template<int bd>
inline void getAbPairs2new(const SxComplex16 *__restrict__ Ab,
                        v4df &Brrii, v4df &Bniirr, v4df *AB)
{
   //v4df A00 = loadcdup((const double*)Ab);
   //v4df A01 = loadcdup((const double*)(Ab + 1));
   //v4df A0r = _mm256_shuffle_pd(A00, A01, 0x0);
   //v4df A0i = _mm256_shuffle_pd(A00, A01, 0xf);
   v4df A0r = loadcdup((const double*)Ab);
   v4df A0i = loadcdup((const double*)(Ab + 1));
   v4df a01 = A0r * Brrii;
   v4df a02 = A0i * Bniirr;
   //v4df A10 = loadcdup((const double*)(Ab + bd));
   //v4df A11 = loadcdup((const double*)(Ab + bd + 1));
   //v4df A1r = _mm256_shuffle_pd(A10, A11, 0x0);
   //v4df A1i = _mm256_shuffle_pd(A10, A11, 0xf);
   v4df A1r = loadcdup((const double*)(Ab + bd));
   v4df A1i = loadcdup((const double*)(Ab + bd + 1));
   v4df a11 = A1r * Brrii;
   v4df a12 = A1i * Bniirr;
   AB[0] = a01 + a02;
   AB[1] = a11 + a12;
}

template<int bd>
inline void getAbPairs8(const SxComplex16 *__restrict__ Ab,
                        const SxComplex16 *__restrict__ Bb,
                        v4df *AB)
{
   // load B0 & B1
   v4df Breim = _mm256_loadu_pd((const double*)Bb);
   v4df Bimre = _mm256_permute_pd(Breim,5);
   getAbPairs2<bd>(Ab, Breim, Bimre, AB);
   getAbPairs2<bd>(Ab + 2 * bd, Breim, Bimre, AB + 2);
   getAbPairs2<bd>(Ab + 4 * bd, Breim, Bimre, AB + 4);
   getAbPairs2<bd>(Ab + 6 * bd, Breim, Bimre, AB + 6);
}

template<int bd>
inline
void gemmm_b_kernel_avx_8_1_bd (int Nc,
                     const SxComplex16 *__restrict__ Ab,
                     const SxComplex16 *__restrict__ Bb,
                     const double *__restrict__ C, int ldc,
                     SxComplex16 *__restrict__ res,
                     int rldc, int rldbc)
{
   const int ba = 8;
   const int bb = 1;
   const int bd2 = bd/2;
   const int Nab = ba * bb;
   SX_CHECK  (bd % 2 == 0, bd);
   SX_CHECK (Nab % 8 == 0, ba, bb);
   v4df AB[bd2 * Nab];
   for (int id2 = 0; id2 < bd2; id2++)
      getAbPairs8<bd> (Ab + 2 * id2, Bb + 2 * id2, AB + Nab * id2);
   for (int ab = 0; ab < Nab; ab+=8)  {
      for (int ic = 0; ic < Nc; ++ic)  {
         v4df cdup = loadcdup (C + ldc * ic);
         v4df *data =  AB + ab;
         v4df S0 = cdup * data[0];
         v4df S1 = cdup * data[1];
         v4df S2 = cdup * data[2];
         v4df S3 = cdup * data[3];
         v4df S4 = cdup * data[4];
         v4df S5 = cdup * data[5];
         v4df S6 = cdup * data[6];
         v4df S7 = cdup * data[7];
         data += Nab;
         for (int id = 1; id < bd2; ++id, data += Nab)  {
            cdup = loadcdup (C + 2 * id + ldc * ic);
            S0 += cdup * data[0];
            S1 += cdup * data[1];
            S2 += cdup * data[2];
            S3 += cdup * data[3];
            S4 += cdup * data[4];
            S5 += cdup * data[5];
            S6 += cdup * data[6];
            S7 += cdup * data[7];
         }
         res[ic + rldbc * 0] += S0;
         res[ic + rldbc * 1] += S1;
         res[ic + rldbc * 2] += S2;
         res[ic + rldbc * 3] += S3;
         res[ic + rldbc * 4] += S4;
         res[ic + rldbc * 5] += S5;
         res[ic + rldbc * 6] += S6;
         res[ic + rldbc * 7] += S7;
      }
   }
}

template<int bufc>
inline
void wb_addOut1_even(v2df*__restrict__ writebuf, int nc, SxComplex16 * __restrict__ res)
{
   // --- add writebuf to res in pairs of 2
   for (int ic = 0; ic < nc; ic +=2, writebuf += 2 * bufc)  {
      v4df old = _mm256_loadu_pd ((double*)(res + ic));
      v4df res1 = _mm256_castpd128_pd256(writebuf[0]);
      v2df res2 = writebuf[bufc];
      v4df res12 = _mm256_insertf128_pd(res1, res2, 1);
      _mm256_storeu_pd((double*)(res + ic), old + res12);
   }
}

template<int bufc>
inline
void wb_addOut1_odd(v2df*__restrict__ writebuf, int nc, SxComplex16 * __restrict__ res)
{
   // --- handle odd nc
   v2df old = _mm_loadu_pd((double*)(res));
   _mm_storeu_pd((double*)(res), old + writebuf[0]);
   writebuf += bufc;
   nc--;
   res++;
   wb_addOut1_even<bufc>(writebuf, nc, res);
}

template<int ba>
inline
void wb_addOut(v2df *writebuf, int nc, SxComplex16 * __restrict__ resEnd, int rldbc)
{
   SxComplex16 * __restrict__ res = resEnd - nc;
   if (nc & 1)  {
      for (int ia = 0; ia < ba; ++ia)
	 wb_addOut1_odd<ba> (writebuf + ia, nc, res + rldbc * ia);
   } else {
      for (int ia = 0; ia < ba; ++ia)
	 wb_addOut1_even<ba> (writebuf + ia, nc, res + rldbc * ia);
   }
}

template<int bd>
inline
void gemmm_b_kernel_avx_12_1_bd (int Nc,
                     const SxComplex16 *__restrict__ Ab,
                     const SxComplex16 *__restrict__ Bb,
                     const double *__restrict__ C, int ldc,
                     SxComplex16 *__restrict__ res,
                     int /*rldc*/, int rldbc)
{
   const int ba = 12;
   const int bb = 1;
   const int bd2 = bd/2;
   const int Nab = ba * bb;
   SX_CHECK  (bd % 2 == 0, bd);
   SX_CHECK (Nab % 12 == 0, ba, bb);
   v4df AB[bd2 * Nab];
   const int bufc = 12;
   const v4li signflip = { 0x8000000000000000L,
                           0x8000000000000000L,
                           0x0000000000000000L,
                           0x0000000000000000L };
   for (int id2 = 0; id2 < bd2; id2++)  {
      /*
      v4df Breim = _mm256_loadu_pd((const double*)(Bb + 2 * id2));
      v4df Bimre = _mm256_permute_pd(Breim,5);
      getAbPairs2<bd>(Ab + 2 * id2, Breim, Bimre, AB + Nab * id2);
      getAbPairs2<bd>(Ab + 2 * bd + 2 * id2, Breim, Bimre, AB + 2 + Nab * id2);
      getAbPairs2<bd>(Ab + 4 * bd + 2 * id2, Breim, Bimre, AB + 4 + Nab * id2);
      getAbPairs2<bd>(Ab + 6 * bd + 2 * id2, Breim, Bimre, AB + 6 + Nab * id2);
      getAbPairs2<bd>(Ab + 8 * bd + 2 * id2, Breim, Bimre, AB + 8 + Nab * id2);
      getAbPairs2<bd>(Ab +10 * bd + 2 * id2, Breim, Bimre, AB +10 + Nab * id2);
      */
      v4df B0 = loadcdup((const double*)(Bb + 2 * id2));
      v4df B1 = loadcdup((const double*)(Bb + 2 * id2 + 1));
      v4df Brrii = _mm256_shuffle_pd(B0, B1, 0xc);
      v4df Biirr = _mm256_shuffle_pd(B0, B1, 0x3);
      v4df Bniirr = _mm256_xor_pd(Biirr,(const v4df)signflip);
      getAbPairs2new<bd>(Ab + 2 * id2, Brrii, Bniirr, AB + Nab * id2);
      getAbPairs2new<bd>(Ab + 2 * bd + 2 * id2, Brrii, Bniirr, AB + 2 + Nab * id2);
      getAbPairs2new<bd>(Ab + 4 * bd + 2 * id2, Brrii, Bniirr, AB + 4 + Nab * id2);
      getAbPairs2new<bd>(Ab + 6 * bd + 2 * id2, Brrii, Bniirr, AB + 6 + Nab * id2);
      getAbPairs2new<bd>(Ab + 8 * bd + 2 * id2, Brrii, Bniirr, AB + 8 + Nab * id2);
      getAbPairs2new<bd>(Ab +10 * bd + 2 * id2, Brrii, Bniirr, AB +10 + Nab * id2);
      //for (int a = 0; a < 12; a+=2)
      //   getAbPairs2new<bd>(Ab + a * bd + 2 * id2, Brrii, Bniirr, AB + a + Nab * id2);
   }
   if (rldbc % 128 == 0)  {
      // --- use writebuffer to avoid cache eviction
      v2df writebuffer[12 * bufc];
      v2df *wb = writebuffer;
      for (int ab = 0; ab < Nab; ab+=12)  {
         for (int ic = 0; ic < Nc; ++ic, wb += 12)  {
            if (wb == writebuffer + 12 * bufc)  {
               // clear writebuffer
               wb_addOut<12>(writebuffer, bufc, res + ic, rldbc);
               wb = writebuffer;
            }
            v4df cdup = loadcdup (C + ldc * ic);
            v4df *data =  AB + ab;
            v4df S0 = cdup * data[0];
            v4df S1 = cdup * data[1];
            v4df S2 = cdup * data[2];
            v4df S3 = cdup * data[3];
            v4df S4 = cdup * data[4];
            v4df S5 = cdup * data[5];
            v4df S6 = cdup * data[6];
            v4df S7 = cdup * data[7];
            v4df S8 = cdup * data[8];
            v4df S9 = cdup * data[9];
            v4df S10 = cdup * data[10];
            v4df S11 = cdup * data[11];
            data += Nab;
            for (int id = 1; id < bd2; ++id, data += Nab)  {
               cdup = loadcdup (C + 2 * id + ldc * ic);
               S0 += cdup * data[0];
               S1 += cdup * data[1];
               S2 += cdup * data[2];
               S3 += cdup * data[3];
               S4 += cdup * data[4];
               S5 += cdup * data[5];
               S6 += cdup * data[6];
               S7 += cdup * data[7];
               S8 += cdup * data[8];
               S9 += cdup * data[9];
               S10 += cdup * data[10];
               S11 += cdup * data[11];
            }
            hadd_store (wb + 0, S0);
            hadd_store (wb + 1, S1);
            hadd_store (wb + 2, S2);
            hadd_store (wb + 3, S3);
            hadd_store (wb + 4, S4);
            hadd_store (wb + 5, S5);
            hadd_store (wb + 6, S6);
            hadd_store (wb + 7, S7);
            hadd_store (wb + 8, S8);
            hadd_store (wb + 9, S9);
            hadd_store (wb + 10, S10);
            hadd_store (wb + 11, S11);
         }
         wb_addOut<12>(writebuffer, Nc % bufc, res + Nc, rldbc);
      }
   } else {
      for (int ab = 0; ab < Nab; ab+=12)  {
         for (int ic = 0; ic < Nc; ++ic)  {
            v4df cdup = loadcdup (C + ldc * ic);
            v4df *data =  AB + ab;
            v4df S0 = cdup * data[0];
            v4df S1 = cdup * data[1];
            v4df S2 = cdup * data[2];
            v4df S3 = cdup * data[3];
            v4df S4 = cdup * data[4];
            v4df S5 = cdup * data[5];
            v4df S6 = cdup * data[6];
            v4df S7 = cdup * data[7];
            v4df S8 = cdup * data[8];
            v4df S9 = cdup * data[9];
            v4df S10 = cdup * data[10];
            v4df S11 = cdup * data[11];
            data += Nab;
            for (int id = 1; id < bd2; ++id, data += Nab)  {
               cdup = loadcdup (C + 2 * id + ldc * ic);
               S0 += cdup * data[0];
               S1 += cdup * data[1];
               S2 += cdup * data[2];
               S3 += cdup * data[3];
               S4 += cdup * data[4];
               S5 += cdup * data[5];
               S6 += cdup * data[6];
               S7 += cdup * data[7];
               S8 += cdup * data[8];
               S9 += cdup * data[9];
               S10 += cdup * data[10];
               S11 += cdup * data[11];
            }
            res[ic + rldbc * 0] += S0;
            res[ic + rldbc * 1] += S1;
            res[ic + rldbc * 2] += S2;
            res[ic + rldbc * 3] += S3;
            res[ic + rldbc * 4] += S4;
            res[ic + rldbc * 5] += S5;
            res[ic + rldbc * 6] += S6;
            res[ic + rldbc * 7] += S7;
            res[ic + rldbc * 8] += S8;
            res[ic + rldbc * 9] += S9;
            res[ic + rldbc * 10] += S10;
            res[ic + rldbc * 11] += S11;
         }
      }
   }
}

inline void getAbQuad(const SxComplex16 *__restrict__ Ab,
                      const SxComplex16 *__restrict__ Bb,
                      v4df &ABr, v4df &ABi)
{
   v4df B0reim = _mm256_loadu_pd((const double*)Bb);
   v4df B0imre = _mm256_permute_pd(B0reim,5);
   v4df B1reim = _mm256_loadu_pd((const double*)(Bb + 2));
   v4df B1imre = _mm256_permute_pd(B1reim,5);
   // load A0 & A1
   v4df A0 = _mm256_loadu_pd((const double*)Ab);
   // re * re, im * im
   v4df AB0rrii =  B0reim * A0;
   v4df A1 = _mm256_loadu_pd((const double*)(Ab + 2));
   // re * re, im * im
   v4df AB1rrii =  B1reim * A1;
   // re * im, im * re
   v4df AB0riri = B0imre * A0;
   // re * im, im * re
   v4df AB1riri = B1imre * A1;

   // get (d0 d2 d1 d3).{r/i}
   v4df AB01r = _mm256_hsub_pd (AB0rrii, AB1rrii);
   v4df AB01i = _mm256_hadd_pd (AB0riri, AB1riri);
   // reorder: get (d2.r d2.i d1.r d1.i) (all wrong lane)
   v4df ABd21 = _mm256_shuffle_pd(AB01r, AB01i, 0x3);
   // swap 128bit lanes
   v4df ABd12 = _mm256_permute2f128_pd(ABd21,ABd21,0x21);
   // get almost final result (d2 and d3 still swapped)
   v4df ABd0132r =  _mm256_shuffle_pd(AB01r, ABd12, 0x4);
   v4df ABd0132i =  _mm256_shuffle_pd(AB01i, ABd12, 0xe);
   ABr = _mm256_permute_pd(ABd0132r, 6);
   ABi = _mm256_permute_pd(ABd0132i, 6);
}

inline void getAbQuadnew(const SxComplex16 *__restrict__ Ab,
                      const SxComplex16 *__restrict__ Bb,
                      v4df &ABr, v4df &ABi)
{
   v4df Br = _mm256_loadu_pd((const double*)Bb);
   v4df Bi = _mm256_loadu_pd((const double*)(Bb + 2));
   v4df Ar = _mm256_loadu_pd((const double*)Ab);
   v4df Ai = _mm256_loadu_pd((const double*)(Ab + 2));
   ABr = Ar * Br - Ai * Bi;
   ABi = Ar * Bi + Ai * Br;
}


inline
void gemmm_b_kernel_avx1_1_16 (int Nc,
                     const SxComplex16 *__restrict__ Ab,
                     const SxComplex16 *__restrict__ Bb,
                     const double *__restrict__ C, int ldc,
                     SxComplex16 *__restrict__ res,
                     int /*rldc*/, int /*rldbc*/)
{
   v4df ab0r, ab0i, ab4r, ab4i, ab8r, ab8i, ab12r, ab12i;
   /*
   getAbQuad(Ab     , Bb     , ab0r,  ab0i);
   getAbQuad(Ab +  4, Bb +  4, ab4r,  ab4i);
   getAbQuad(Ab +  8, Bb +  8, ab8r,  ab8i);
   getAbQuad(Ab + 12, Bb + 12, ab12r, ab12i);
   */
   getAbQuadnew(Ab     , Bb     , ab0r,  ab0i);
   getAbQuadnew(Ab +  4, Bb +  4, ab4r,  ab4i);
   getAbQuadnew(Ab +  8, Bb +  8, ab8r,  ab8i);
   getAbQuadnew(Ab + 12, Bb + 12, ab12r, ab12i);

   for (int ic = 0; ic < Nc; ic++)  {
      v4df c = _mm256_loadu_pd(C + ldc * ic);
      v4df resr = ab0r * c;
      v4df resi = ab0i * c;
      c = _mm256_loadu_pd(C + 4 + ldc * ic);
      resr += ab4r * c;
      resi += ab4i * c;
      c = _mm256_loadu_pd(C + 8 + ldc * ic);
      resr += ab8r * c;
      resi += ab8i * c;
      c = _mm256_loadu_pd(C + 12 + ldc * ic);
      resr += ab12r * c;
      resi += ab12i * c;

      // horizontal sum, store
      v4df halfSum = _mm256_hadd_pd(resr, resi);
      v2df resUp  = _mm256_extractf128_pd(halfSum,1);
      v2df resLow = _mm256_castpd256_pd128(halfSum);
      v2df old = _mm_loadu_pd((double*)(res + ic));
      _mm_storeu_pd((double*)(res + ic), old + resUp + resLow);
   }
}

template<int bd>
inline
void gemmm_b_kernel_avx_1_1_bd (int Nc,
                     const SxComplex16 *__restrict__ Ab,
                     const SxComplex16 *__restrict__ Bb,
                     const double *__restrict__ C, int ldc,
                     SxComplex16 *__restrict__ res,
                     int rldc, int rldbc)
{
   SX_CHECK  (bd % 16 == 0, bd);
   for (int id = 0; id < bd; id += 16)
      gemmm_b_kernel_avx1_1_16(Nc, Ab + id, Bb + id, C + id, ldc,
                               res, rldc, rldbc);
}

inline
void CX_mulB0(const v4df CX_d1234re, const v4df CX_d1234im,
             const SxComplex16 *__restrict__ B,
             v4df &BCXd12_rrii, v4df &BCXd12_riri,
             v4df &BCXd34_rrii, v4df &BCXd34_riri)
{
   // CX_d1234 contains the c-sum of C * X for current (a,b) and (d1,d2,d3,d4)
   // resort numbers into (CX_d1.re CX_d1.im CX_d2.re CX_d3.im)
   //                 and (CX_d3.re CX_d3.im CX_d4.re CX_d4.im)
   v4df d13reim = _mm256_shuffle_pd (CX_d1234re, CX_d1234im, 0x0);
   v4df d24reim = _mm256_shuffle_pd (CX_d1234re, CX_d1234im, 0xf);
   v4df CX_d12 = _mm256_permute2f128_pd(d13reim,d24reim,0x20);
   v4df CX_d34 = _mm256_permute2f128_pd(d13reim,d24reim,0x31);

   // --- load and multiply with B
   v4df Bd12 =  _mm256_loadu_pd((const double*)(B));
   v4df Bd34 =  _mm256_loadu_pd((const double*)(B + 2));
   BCXd12_rrii = CX_d12 * Bd12;
   BCXd34_rrii = CX_d34 * Bd34;
   BCXd12_riri = CX_d12 * _mm256_permute_pd(Bd12,5);
   BCXd34_riri = CX_d34 * _mm256_permute_pd(Bd34,5);
}

inline
void CX_mulB(const v4df CX_d1234re, const v4df CX_d1234im,
             const SxComplex16 *__restrict__ B,
             v4df &BCXd1324r, v4df &BCXd1324i)
{
   v4df BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri;
   CX_mulB0(CX_d1234re, CX_d1234im, B,
            BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);
   // BCXd1324{r/i} is (BCX_d1 BCX_d3 BCX_d2 BCX_d4).{re/im}
}


inline
void CX_mulB(const v4df CX_d1234re, const v4df CX_d1234im,
             const SxComplex16 *__restrict__ B,
             v4df &BCXd12_rrii, v4df &BCXd12_riri,
             v4df &BCXd34_rrii, v4df &BCXd34_riri)
{
   // CX_d1234 contains the c-sum of C * X for current (a,b) and (d1,d2,d3,d4)
   // resort numbers into (CX_d1.re CX_d1.im CX_d2.re CX_d3.im)
   //                 and (CX_d3.re CX_d3.im CX_d4.re CX_d4.im)
   v4df d13reim = _mm256_shuffle_pd (CX_d1234re, CX_d1234im, 0x0);
   v4df d24reim = _mm256_shuffle_pd (CX_d1234re, CX_d1234im, 0xf);
   v4df CX_d12 = _mm256_permute2f128_pd(d13reim,d24reim,0x20);
   v4df CX_d34 = _mm256_permute2f128_pd(d13reim,d24reim,0x31);

   // --- load and multiply with B, add to result
   v4df Bd12 =  _mm256_loadu_pd((const double*)(B));
   v4df Bd34 =  _mm256_loadu_pd((const double*)(B + 2));
   BCXd12_rrii += CX_d12 * Bd12;
   BCXd34_rrii += CX_d34 * Bd34;
   BCXd12_riri += CX_d12 * _mm256_permute_pd(Bd12,5);
   BCXd34_riri += CX_d34 * _mm256_permute_pd(Bd34,5);
}

inline
void addToRes(v4df &BCXd1324r, v4df &BCXd1324i, SxComplex16 *__restrict__ res)
{
   // --- resort into (BCX_d1.re BCX_d1.im BCX_d2.re BCX_d2.im)
   //            and  (BCX_d3.re BCX_d3.im BCX_d4.re BCX_d4.im)
   //     and add to result
   v4df BCXd12reim =  _mm256_shuffle_pd (BCXd1324r, BCXd1324i , 0x0);
   v4df old =  _mm256_loadu_pd((double*)res);
   _mm256_storeu_pd((double*)res, old + BCXd12reim);
   v4df BCXd34reim =  _mm256_shuffle_pd (BCXd1324r, BCXd1324i , 0xf);
   old =  _mm256_loadu_pd((double*)(res + 2));
   _mm256_storeu_pd((double*)(res + 2), old + BCXd34reim);
}

inline void gemm3m_kernel_avx_1_1_4(int Nc,
            const SxComplex16 *__restrict__ B, int /* ldb */,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int /*xldc*/, int /*xldbc*/,
                  SxComplex16 *__restrict__ res, int /*lda*/)
{
   /*
   for (int a = 0; a < ba; ++a)  {
         for (int d = 0; d < bd; ++d)  {
      for (int b = 0; b < bb; ++b)  {
            TM<CRTYPE_D,CRTYPE_C>::Res S = 0.;
            for (int c = 0; c < Nc; ++c)  {
               S += C[d + ldc * c] * X[c + xldc * b + xldbc * a];
            }
            res[d + lda * a] += B[d + ldb * b] * S;
         }
      }
   }
   */
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));
   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;
   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;
   }
   v4df BCXd1324r, BCXd1324i;
   CX_mulB(CX_d1234re, CX_d1234im, B, BCXd1324r, BCXd1324i);

   addToRes (BCXd1324r, BCXd1324i, res);
}

inline void gemm3m_kernel_avx_1_6_4(int Nc,
            const SxComplex16 *__restrict__ B, int ldb,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int xldc, int /*xldbc*/,
                  SxComplex16 *__restrict__ res, int /* lda */)
{
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));

   // b
   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;

   // b + 1
   X_re = _mm256_broadcast_sd ((const double*)(X + xldc));
   v4df CX_d1234re_b1 =  C_d1234 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + xldc) + 1);
   v4df CX_d1234im_b1 =  C_d1234 * X_im;

   // b + 2
   X_re = _mm256_broadcast_sd ((const double*)(X + 2 * xldc));
   v4df CX_d1234re_b2 =  C_d1234 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + 2 * xldc) + 1);
   v4df CX_d1234im_b2 =  C_d1234 * X_im;

   // b + 3
   X_re = _mm256_broadcast_sd ((const double*)(X + 3 * xldc));
   v4df CX_d1234re_b3 =  C_d1234 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + 3 * xldc) + 1);
   v4df CX_d1234im_b3 =  C_d1234 * X_im;

   // b + 4
   X_re = _mm256_broadcast_sd ((const double*)(X + 4 * xldc));
   v4df CX_d1234re_b4 =  C_d1234 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + 4 * xldc) + 1);
   v4df CX_d1234im_b4 =  C_d1234 * X_im;

   // b + 5
   X_re = _mm256_broadcast_sd ((const double*)(X + 5 * xldc));
   v4df CX_d1234re_b5 =  C_d1234 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + 5 * xldc) + 1);
   v4df CX_d1234im_b5 =  C_d1234 * X_im;

   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      //C_d1234 = *((v4df*)(C + ldc * c));

      // b
      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;

      // b + 1
      X_re = _mm256_broadcast_sd ((const double*)(X + c + xldc));
      CX_d1234re_b1 +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + xldc) + 1);
      CX_d1234im_b1 +=  C_d1234 * X_im;

      // b + 2
      X_re = _mm256_broadcast_sd ((const double*)(X + c + 2 * xldc));
      CX_d1234re_b2 +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + 2 * xldc) + 1);
      CX_d1234im_b2 +=  C_d1234 * X_im;

      // b + 3
      X_re = _mm256_broadcast_sd ((const double*)(X + c + 3 * xldc));
      CX_d1234re_b3 +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + 3 * xldc) + 1);
      CX_d1234im_b3 +=  C_d1234 * X_im;

      // b + 4
      X_re = _mm256_broadcast_sd ((const double*)(X + c + 4 * xldc));
      CX_d1234re_b4 +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + 4 * xldc) + 1);
      CX_d1234im_b4 +=  C_d1234 * X_im;

      // b + 5
      X_re = _mm256_broadcast_sd ((const double*)(X + c + 5 * xldc));
      CX_d1234re_b5 +=  C_d1234 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + 5 * xldc) + 1);
      CX_d1234im_b5 +=  C_d1234 * X_im;
   }
   // b + 0
   v4df BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri;
   CX_mulB0(CX_d1234re, CX_d1234im, B,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d1234re_b1, CX_d1234im_b1, B + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 2
   CX_mulB(CX_d1234re_b2, CX_d1234im_b2, B + 2*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 3
   CX_mulB(CX_d1234re_b3, CX_d1234im_b3, B + 3*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 4
   CX_mulB(CX_d1234re_b4, CX_d1234im_b4, B + 4*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 5
   CX_mulB(CX_d1234re_b5, CX_d1234im_b5, B + 5*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);
   v4df BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   v4df BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res);
}

inline void gemm3m_kernel_avx_1_3_8(int Nc,
            const SxComplex16 *__restrict__ B, int ldb,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int xldc, int /* xldbc */,
                  SxComplex16 *__restrict__ res, int /* lda */)
{
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));
   v4df C_d5678 = _mm256_loadu_pd((const double*)(C + 4));

   // b
   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df CX_d5678re =  C_d5678 * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;
   v4df CX_d5678im =  C_d5678 * X_im;

   // b + 1
   X_re = _mm256_broadcast_sd ((const double*)(X + xldc));
   v4df CX_d1234re_b1 =  C_d1234 * X_re;
   v4df CX_d5678re_b1 =  C_d5678 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + xldc) + 1);
   v4df CX_d1234im_b1 =  C_d1234 * X_im;
   v4df CX_d5678im_b1 =  C_d5678 * X_im;

   // b + 2
   X_re = _mm256_broadcast_sd ((const double*)(X + 2 * xldc));
   v4df CX_d1234re_b2 =  C_d1234 * X_re;
   v4df CX_d5678re_b2 =  C_d5678 * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + 2 * xldc) + 1);
   v4df CX_d1234im_b2 =  C_d1234 * X_im;
   v4df CX_d5678im_b2 =  C_d5678 * X_im;

   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      C_d5678 = _mm256_loadu_pd((const double*)(C + 4 + ldc * c));
      //C_d1234 = *((v4df*)(C + ldc * c));

      // b
      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      CX_d5678re +=  C_d5678 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;
      CX_d5678im +=  C_d5678 * X_im;

      // b + 1
      X_re = _mm256_broadcast_sd ((const double*)(X + c + xldc));
      CX_d1234re_b1 +=  C_d1234 * X_re;
      CX_d5678re_b1 +=  C_d5678 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + xldc) + 1);
      CX_d1234im_b1 +=  C_d1234 * X_im;
      CX_d5678im_b1 +=  C_d5678 * X_im;

      // b + 2
      X_re = _mm256_broadcast_sd ((const double*)(X + c + 2 * xldc));
      CX_d1234re_b2 +=  C_d1234 * X_re;
      CX_d5678re_b2 +=  C_d5678 * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + 2 * xldc) + 1);
      CX_d1234im_b2 +=  C_d1234 * X_im;
      CX_d5678im_b2 +=  C_d5678 * X_im;

   }
   // b + 0
   v4df BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri;
   CX_mulB0(CX_d1234re, CX_d1234im, B,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d1234re_b1, CX_d1234im_b1, B + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 2
   CX_mulB(CX_d1234re_b2, CX_d1234im_b2, B + 2*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   v4df BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   v4df BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res);


   // b + 0
   CX_mulB0(CX_d5678re, CX_d5678im, B + 4,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d5678re_b1, CX_d5678im_b1, B + 4 + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 2
   CX_mulB(CX_d5678re_b2, CX_d5678im_b2, B + 4 + 2*ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res + 4);

}

inline void gemm3m_kernel_avx_1_2_12(int Nc,
            const SxComplex16 *__restrict__ B, int ldb,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int xldc, int /*xldbc*/,
                  SxComplex16 *__restrict__ res, int /*lda*/)
{
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));
   v4df C_d5678 = _mm256_loadu_pd((const double*)(C + 4));
   v4df C_d9abc = _mm256_loadu_pd((const double*)(C + 8));

   // b
   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df CX_d5678re =  C_d5678 * X_re;
   v4df CX_d9abcre =  C_d9abc * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;
   v4df CX_d5678im =  C_d5678 * X_im;
   v4df CX_d9abcim =  C_d9abc * X_im;

   // b + 1
   X_re = _mm256_broadcast_sd ((const double*)(X + xldc));
   v4df CX_d1234re_b1 =  C_d1234 * X_re;
   v4df CX_d5678re_b1 =  C_d5678 * X_re;
   v4df CX_d9abcre_b1 =  C_d9abc * X_re;
   X_im = _mm256_broadcast_sd ((const double*)(X + xldc) + 1);
   v4df CX_d1234im_b1 =  C_d1234 * X_im;
   v4df CX_d5678im_b1 =  C_d5678 * X_im;
   v4df CX_d9abcim_b1 =  C_d9abc * X_im;

   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      C_d5678 = _mm256_loadu_pd((const double*)(C + 4 + ldc * c));
      C_d9abc = _mm256_loadu_pd((const double*)(C + 8 + ldc * c));

      // b
      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      CX_d5678re +=  C_d5678 * X_re;
      CX_d9abcre +=  C_d9abc * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;
      CX_d5678im +=  C_d5678 * X_im;
      CX_d9abcim +=  C_d9abc * X_im;

      // b + 1
      X_re = _mm256_broadcast_sd ((const double*)(X + c + xldc));
      CX_d1234re_b1 +=  C_d1234 * X_re;
      CX_d5678re_b1 +=  C_d5678 * X_re;
      CX_d9abcre_b1 +=  C_d9abc * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c + xldc) + 1);
      CX_d1234im_b1 +=  C_d1234 * X_im;
      CX_d5678im_b1 +=  C_d5678 * X_im;
      CX_d9abcim_b1 +=  C_d9abc * X_im;

   }
   // b + 0
   v4df BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri;
   CX_mulB0(CX_d1234re, CX_d1234im, B,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d1234re_b1, CX_d1234im_b1, B + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   v4df BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   v4df BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res);

   // b + 0
   CX_mulB0(CX_d5678re, CX_d5678im, B + 4,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d5678re_b1, CX_d5678im_b1, B + 4 + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res + 4);

   // b + 0
   CX_mulB0(CX_d9abcre, CX_d9abcim, B + 8,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   // b + 1
   CX_mulB(CX_d9abcre_b1, CX_d9abcim_b1, B + 8 + ldb,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res + 8);
}

inline
void mulAddToRes(v4df &CX_d1234re, v4df &CX_d1234im,
              const SxComplex16 *__restrict__ B,
              SxComplex16 *__restrict__ res)
{
   v4df BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri;
   CX_mulB0(CX_d1234re, CX_d1234im, B,
           BCXd12_rrii, BCXd12_riri, BCXd34_rrii, BCXd34_riri);

   v4df BCXd1324r = _mm256_hsub_pd (BCXd12_rrii, BCXd34_rrii);
   v4df BCXd1324i = _mm256_hadd_pd (BCXd12_riri, BCXd34_riri);

   addToRes (BCXd1324r, BCXd1324i, res);
}

inline void gemm3m_kernel_avx_1_1_24(int Nc,
            const SxComplex16 *__restrict__ B, int /* ldb */,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int /* xldc */, int /*xldbc*/,
                  SxComplex16 *__restrict__ res, int /*lda*/)
{
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));
   v4df C_d5678 = _mm256_loadu_pd((const double*)(C + 4));
   v4df C_d9abc = _mm256_loadu_pd((const double*)(C + 8));
   v4df C_ddefg = _mm256_loadu_pd((const double*)(C + 12));
   v4df C_dhijk = _mm256_loadu_pd((const double*)(C + 16));
   v4df C_dlmno = _mm256_loadu_pd((const double*)(C + 20));

   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df CX_d5678re =  C_d5678 * X_re;
   v4df CX_d9abcre =  C_d9abc * X_re;
   v4df CX_ddefgre =  C_ddefg * X_re;
   v4df CX_dhijkre =  C_dhijk * X_re;
   v4df CX_dlmnore =  C_dlmno * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;
   v4df CX_d5678im =  C_d5678 * X_im;
   v4df CX_d9abcim =  C_d9abc * X_im;
   v4df CX_ddefgim =  C_ddefg * X_im;
   v4df CX_dhijkim =  C_dhijk * X_im;
   v4df CX_dlmnoim =  C_dlmno * X_im;

   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      C_d5678 = _mm256_loadu_pd((const double*)(C + 4 + ldc * c));
      C_d9abc = _mm256_loadu_pd((const double*)(C + 8 + ldc * c));
      C_ddefg = _mm256_loadu_pd((const double*)(C + 12 + ldc * c));
      C_dhijk = _mm256_loadu_pd((const double*)(C + 16 + ldc * c));
      C_dlmno = _mm256_loadu_pd((const double*)(C + 20 + ldc * c));

      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      CX_d5678re +=  C_d5678 * X_re;
      CX_d9abcre +=  C_d9abc * X_re;
      CX_ddefgre +=  C_ddefg * X_re;
      CX_dhijkre +=  C_dhijk * X_re;
      CX_dlmnore +=  C_dlmno * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;
      CX_d5678im +=  C_d5678 * X_im;
      CX_d9abcim +=  C_d9abc * X_im;
      CX_ddefgim +=  C_ddefg * X_im;
      CX_dhijkim +=  C_dhijk * X_im;
      CX_dlmnoim +=  C_dlmno * X_im;
   }
   mulAddToRes (CX_d1234re, CX_d1234im, B + 0, res + 0);
   mulAddToRes (CX_d5678re, CX_d5678im, B + 4, res + 4);
   mulAddToRes (CX_d9abcre, CX_d9abcim, B + 8, res + 8);
   mulAddToRes (CX_ddefgre, CX_ddefgim, B + 12, res + 12);
   mulAddToRes (CX_dhijkre, CX_dhijkim, B + 16, res + 16);
   mulAddToRes (CX_dlmnore, CX_dlmnoim, B + 20, res + 20);
}

inline void gemm3m_kernel_avx_1_1_12(int Nc,
            const SxComplex16 *__restrict__ B, int /*ldb*/,
            const double *__restrict__ C, int ldc,
            const SxComplex16 *__restrict__ X, int /*xldc*/, int /*xldbc*/,
                  SxComplex16 *__restrict__ res, int /*lda*/)
{
   v4df C_d1234 = _mm256_loadu_pd((const double*)(C));
   v4df C_d5678 = _mm256_loadu_pd((const double*)(C + 4));
   v4df C_d9abc = _mm256_loadu_pd((const double*)(C + 8));

   v4df X_re = _mm256_broadcast_sd ((const double*)(X));
   v4df CX_d1234re =  C_d1234 * X_re;
   v4df CX_d5678re =  C_d5678 * X_re;
   v4df CX_d9abcre =  C_d9abc * X_re;
   v4df X_im = _mm256_broadcast_sd ((const double*)(X) + 1);
   v4df CX_d1234im =  C_d1234 * X_im;
   v4df CX_d5678im =  C_d5678 * X_im;
   v4df CX_d9abcim =  C_d9abc * X_im;

   for (int c = 1; c < Nc; ++c)  {
      C_d1234 = _mm256_loadu_pd((const double*)(C + ldc * c));
      C_d5678 = _mm256_loadu_pd((const double*)(C + 4 + ldc * c));
      C_d9abc = _mm256_loadu_pd((const double*)(C + 8 + ldc * c));

      X_re = _mm256_broadcast_sd ((const double*)(X + c));
      CX_d1234re +=  C_d1234 * X_re;
      CX_d5678re +=  C_d5678 * X_re;
      CX_d9abcre +=  C_d9abc * X_re;
      X_im = _mm256_broadcast_sd ((const double*)(X + c) + 1);
      CX_d1234im +=  C_d1234 * X_im;
      CX_d5678im +=  C_d5678 * X_im;
      CX_d9abcim +=  C_d9abc * X_im;
   }
   mulAddToRes (CX_d1234re, CX_d1234im, B + 0, res + 0);
   mulAddToRes (CX_d5678re, CX_d5678im, B + 4, res + 4);
   mulAddToRes (CX_d9abcre, CX_d9abcim, B + 8, res + 8);
}

void gemm3mPrepareCdup (const double *C, int ldc, v4df *Cdup, int Nc)
{
   for (int ic = 0; ic < Nc; ++ic, Cdup+=4, C+= ldc) {
      double *C1324_d = reinterpret_cast<double *>(Cdup);
      C1324_d[0] = C[0];
      C1324_d[1] = C[0];
      C1324_d[2] = C[1];
      C1324_d[3] = C[1];
      C1324_d[4] = C[2];
      C1324_d[5] = C[2];
      C1324_d[6] = C[3];
      C1324_d[7] = C[3];
      C1324_d[8] = C[4];
      C1324_d[9] = C[4];
      C1324_d[10] = C[5];
      C1324_d[11] = C[5];
      C1324_d[12] = C[6];
      C1324_d[13] = C[6];
      C1324_d[14] = C[7];
      C1324_d[15] = C[7];
   }
}

inline void gemm3mPrepareBdup (const SxComplex16 *B, int ldb, v4df *Bdup, int Nb, int bd2)
{
   const v4li signflip = { 0x8000000000000000L,
                           0x0000000000000000L,
                           0x8000000000000000L,
                           0x0000000000000000L };
   for (int ib = 0; ib < Nb; ++ib, Bdup+=4, B+= ldb) {
      for (int d = 0; d < bd2; d+=4)  {
         _mm_prefetch(B + d + 32, _MM_HINT_T1); // L2 cache (read)
         v4df Breim = _mm256_loadu_pd((const double*)(B + d));
         Bdup[0 + d * Nb] = _mm256_permute_pd(Breim,0); // rere
         v4df Bimim =  _mm256_permute_pd(Breim,0xf);
         Bdup[1 + d * Nb] = _mm256_xor_pd(Bimim,(const v4df)signflip);

         Breim = _mm256_loadu_pd((const double*)(B+d+2));
         Bdup[2 + d * Nb] = _mm256_permute_pd(Breim,0); // rere
         Bimim =  _mm256_permute_pd(Breim,0xf);
         Bdup[3 + d * Nb] = _mm256_xor_pd(Bimim,(const v4df)signflip);
      }
   }
}

inline void
gemm3m_avx_CX_3_1_8 (int Nc,
                     const v4df* Cdup,
                     const SxComplex16 *__restrict__ X, int xldbc,
                           v4df *__restrict__ CX)
{
   using namespace std;

   v4df C1122 = Cdup[0];
   v4df X_a0 = loadcdup((const double*)X            );
   v4df X_a1 = loadcdup((const double *)(X + 1 * xldbc));
   v4df X_a2 = loadcdup((const double *)(X + 2 * xldbc));
   v4df CX12_a0 = C1122 * X_a0;
   v4df CX12_a1 = C1122 * X_a1;
   v4df CX12_a2 = C1122 * X_a2;

   //asm volatile("" ::: "memory");
   v4df C3344 = Cdup[1];

   v4df CX34_a0 = C3344 * X_a0;
   v4df CX34_a1 = C3344 * X_a1;
   v4df CX34_a2 = C3344 * X_a2;

   v4df C5566 = Cdup[2];
   v4df CX56_a0 = C5566 * X_a0;
   v4df CX56_a1 = C5566 * X_a1;
   v4df CX56_a2 = C5566 * X_a2;

   //asm volatile("" ::: "memory");
   v4df C7788 = Cdup[3];
   v4df CX78_a0 = C7788 * X_a0;
   v4df CX78_a1 = C7788 * X_a1;
   v4df CX78_a2 = C7788 * X_a2;
   Cdup+=4; X++;
   for (int c = 1; c < Nc; c++, Cdup+=4, X++)  {

      C1122 = Cdup[0];
      X_a0 = loadcdup((const double*)X            );
      X_a1 = loadcdup((const double *)(X + 1 * xldbc));
      X_a2 = loadcdup((const double *)(X + 2 * xldbc));
      CX12_a0 += C1122 * X_a0;
      CX12_a1 += C1122 * X_a1;
      CX12_a2 += C1122 * X_a2;

      C3344 = Cdup[1];

      CX34_a0 += C3344 * X_a0;
      CX34_a1 += C3344 * X_a1;
      CX34_a2 += C3344 * X_a2;

      C5566 = Cdup[2];
      CX56_a0 += C5566 * X_a0;
      CX56_a1 += C5566 * X_a1;
      CX56_a2 += C5566 * X_a2;

      C7788 = Cdup[3];
      CX78_a0 += C7788 * X_a0;
      CX78_a1 += C7788 * X_a1;
      CX78_a2 += C7788 * X_a2;
   }
   CX[0] = CX12_a0;
   CX[1] = CX12_a1;
   CX[2] = CX12_a2;
   CX[3] = CX34_a0;
   CX[4] = CX34_a1;
   CX[5] = CX34_a2;
   CX[6] = CX56_a0;
   CX[7] = CX56_a1;
   CX[8] = CX56_a2;
   CX[9] = CX78_a0;
   CX[10] = CX78_a1;
   CX[11] = CX78_a2;
}

inline void
gemm3m_avx_CX_1_1_8 (int Nc,
                     const v4df* Cdup,
                     const SxComplex16 *__restrict__ X,
                           v4df *__restrict__ CX)
{
   using namespace std;

   v4df C1122 = Cdup[0];
   v4df X_a0 = loadcdup((const double*)X            );
   v4df CX12_a0 = C1122 * X_a0;

   v4df C3344 = Cdup[1];

   v4df CX34_a0 = C3344 * X_a0;

   v4df C5566 = Cdup[2];
   v4df CX56_a0 = C5566 * X_a0;

   //asm volatile("" ::: "memory");
   v4df C7788 = Cdup[3];
   v4df CX78_a0 = C7788 * X_a0;
   Cdup+=4; X++;
   for (int c = 1; c < Nc; c++, Cdup+=4, X++)  {

      C1122 = Cdup[0];
      X_a0 = loadcdup((const double*)X            );
      CX12_a0 += C1122 * X_a0;

      C3344 = Cdup[1];

      CX34_a0 += C3344 * X_a0;

      C5566 = Cdup[2];
      CX56_a0 += C5566 * X_a0;

      C7788 = Cdup[3];
      CX78_a0 += C7788 * X_a0;
   }
   CX[0] = CX12_a0;
   CX[1] = CX34_a0;
   CX[2] = CX56_a0;
   CX[3] = CX78_a0;
}

inline void gemm3m_avx_CXmulBd_3_1_4 (int Nb,
                        const v4df *__restrict__ CX,
                        const v4df *__restrict__ Bd,
                        SxComplex16 *__restrict__ res, int lda)
{
   using namespace std;
   v4df res_a0 = _mm256_loadu_pd((const double*)(res          ));
   v4df res_a1 = _mm256_loadu_pd((const double*)(res +     lda));
   v4df res_a2 = _mm256_loadu_pd((const double*)(res + 2 * lda));
   v4df res_a0_d34 = _mm256_loadu_pd((const double*)(res + 2         ));
   v4df res_a1_d34 = _mm256_loadu_pd((const double*)(res + 2 +     lda));
   v4df res_a2_d34 = _mm256_loadu_pd((const double*)(res + 2 + 2 * lda));

   for (int b = 0; b < Nb; ++b, Bd += 4, CX += 12)  {

      v4df Brere = Bd[0];
      v4df CX_a0 = CX[0];
      v4df CX_a1 = CX[1];
      v4df CX_a2 = CX[2];

      res_a0 += CX_a0 * Brere;
      res_a1 += CX_a1 * Brere;
      res_a2 += CX_a2 * Brere;

      v4df Brere_d34 = Bd[2];
      v4df CX_a0_d34 = CX[3];
      v4df CX_a1_d34 = CX[4];
      v4df CX_a2_d34 = CX[5];

      res_a0_d34 += CX_a0_d34 * Brere_d34;
      res_a1_d34 += CX_a1_d34 * Brere_d34;
      res_a2_d34 += CX_a2_d34 * Brere_d34;

      v4df BnImim = Bd[1];
      v4df CX_a0_imre, CX_a1_imre, CX_a2_imre;
      CX_a0_imre = _mm256_permute_pd(CX_a0,5);
      CX_a1_imre = _mm256_permute_pd(CX_a1,5);
      CX_a2_imre = _mm256_permute_pd(CX_a2,5);
      res_a0 += CX_a0_imre * BnImim;
      res_a1 += CX_a1_imre * BnImim;
      res_a2 += CX_a2_imre * BnImim;

      v4df BnImim_d34 = Bd[3];
      v4df CX_a0_imre_d34, CX_a1_imre_d34, CX_a2_imre_d34;
      CX_a0_imre_d34 = _mm256_permute_pd(CX_a0_d34,5);
      CX_a1_imre_d34 = _mm256_permute_pd(CX_a1_d34,5);
      CX_a2_imre_d34 = _mm256_permute_pd(CX_a2_d34,5);
      res_a0_d34 += CX_a0_imre_d34 * BnImim_d34;
      res_a1_d34 += CX_a1_imre_d34 * BnImim_d34;
      res_a2_d34 += CX_a2_imre_d34 * BnImim_d34;
   }
   //TODO: widen to a=8
   _mm256_storeu_pd((double *)(res          ), res_a0);
   _mm256_storeu_pd((double *)(res +     lda), res_a1);
   _mm256_storeu_pd((double *)(res + 2 * lda), res_a2);
   _mm256_storeu_pd((double *)(res + 2          ), res_a0_d34);
   _mm256_storeu_pd((double *)(res + 2 +     lda), res_a1_d34);
   _mm256_storeu_pd((double *)(res + 2 + 2 * lda), res_a2_d34);
}

inline void gemm3m_avx_CXmulBd_1_3_4 (int Nb,
                        const v4df *__restrict__ CX,
                        const v4df *__restrict__ Bd,
                        SxComplex16 *__restrict__ res)
{
   using namespace std;
   v4df res_b0 = _mm256_loadu_pd((const double*)(res));
   v4df res_b1 = {0.,0.,0.,0.}, res_b2 = {0., 0., 0., 0.};
   v4df res_b0_d34 = _mm256_loadu_pd((const double*)(res + 2));
   v4df res_b1_d34 = {0., 0., 0., 0.}, res_b2_d34 = {0., 0., 0., 0.};

   for (int b = 0; b < Nb; b+=3, Bd += 12, CX += 12)  {

      v4df CX_b0 = CX[0];
      v4df CX_b1 = CX[1];
      v4df CX_b2 = CX[2];

      res_b0 += CX_b0 * Bd[0];
      res_b1 += CX_b1 * Bd[4];
      res_b2 += CX_b2 * Bd[8];

      v4df CX_b0_d34 = CX[3];
      v4df CX_b1_d34 = CX[4];
      v4df CX_b2_d34 = CX[5];

      res_b0_d34 += CX_b0_d34 * Bd[2];
      res_b1_d34 += CX_b1_d34 * Bd[6];
      res_b2_d34 += CX_b2_d34 * Bd[10];

      v4df CX_b0_imre = _mm256_permute_pd(CX_b0,5);
      v4df CX_b1_imre = _mm256_permute_pd(CX_b1,5);
      v4df CX_b2_imre = _mm256_permute_pd(CX_b2,5);
      res_b0 += CX_b0_imre * Bd[1];
      res_b1 += CX_b1_imre * Bd[5];
      res_b2 += CX_b2_imre * Bd[9];

      v4df CX_b0_imre_d34 = _mm256_permute_pd(CX_b0_d34,5);
      v4df CX_b1_imre_d34 = _mm256_permute_pd(CX_b1_d34,5);
      v4df CX_b2_imre_d34 = _mm256_permute_pd(CX_b2_d34,5);
      res_b0_d34 += CX_b0_imre_d34 * Bd[3];
      res_b1_d34 += CX_b1_imre_d34 * Bd[7];
      res_b2_d34 += CX_b2_imre_d34 * Bd[11];
   }
   res_b0 += res_b1;
   res_b0_d34 += res_b1_d34;
   res_b0 += res_b2;
   res_b0_d34 += res_b2_d34;
   _mm256_storeu_pd((double *)(res    ), res_b0);
   _mm256_storeu_pd((double *)(res + 2), res_b0_d34);
}


inline void gemm3m_avx_CXmulBd_1_1_4 (int Nb,
                        const v4df *__restrict__ CX,
                        const v4df *__restrict__ Bd,
                        SxComplex16 *__restrict__ res, int /* lda */)
{
   using namespace std;
   v4df res_a0 = _mm256_loadu_pd((const double*)(res          ));
   v4df res_a0_d34 = _mm256_loadu_pd((const double*)(res + 2         ));

   for (int b = 0; b < Nb; ++b, Bd += 4, CX += 4)  {

      v4df Brere = Bd[0];
      v4df CX_a0 = CX[0];

      res_a0 += CX_a0 * Brere;

      v4df Brere_d34 = Bd[2];
      v4df CX_a0_d34 = CX[1];

      res_a0_d34 += CX_a0_d34 * Brere_d34;

      v4df BnImim = Bd[1];
      v4df CX_a0_imre = _mm256_permute_pd(CX_a0,5);
      res_a0 += CX_a0_imre * BnImim;

      v4df BnImim_d34 = Bd[3];
      v4df CX_a0_imre_d34;
      CX_a0_imre_d34 = _mm256_permute_pd(CX_a0_d34,5);
      res_a0_d34 += CX_a0_imre_d34 * BnImim_d34;
   }
   _mm256_storeu_pd((double *)(res          ), res_a0);
   _mm256_storeu_pd((double *)(res + 2          ), res_a0_d34);
}

inline void gemm3m_avx_CXmulB_1_1_4 (int Nb,
                        const v4df *__restrict__ CX,
                        const SxComplex16 *__restrict__ B, int ldb,
                        SxComplex16 *__restrict__ res, int /* lda */)
{
   using namespace std;
   v4df res_a0 = _mm256_loadu_pd((const double*)(res          ));
   v4df res_a0_d34 = _mm256_loadu_pd((const double*)(res + 2         ));

   const v4li signflip = { 0x8000000000000000L,
                           0x0000000000000000L,
                           0x8000000000000000L,
                           0x0000000000000000L };
   for (int b = 0; b < Nb; ++b, B += ldb, CX += 4)  {

      v4df Breim = _mm256_loadu_pd((const double*)(B));
      v4df Brere = _mm256_permute_pd(Breim,0);
      v4df CX_a0 = CX[0];

      res_a0 += CX_a0 * Brere;

      v4df Breim_d34 = _mm256_loadu_pd((const double*)(B + 2));
      v4df Brere_d34 = _mm256_permute_pd(Breim_d34,0);
      v4df CX_a0_d34 = CX[1];

      res_a0_d34 += CX_a0_d34 * Brere_d34;

      v4df Bimim =  _mm256_permute_pd(Breim,0xf);
      v4df BnImim = _mm256_xor_pd(Bimim,(const v4df)signflip);
      v4df CX_a0_imre = _mm256_permute_pd(CX_a0,5);
      res_a0 += CX_a0_imre * BnImim;

      v4df Bimim_d34 =  _mm256_permute_pd(Breim_d34,0xf);
      v4df BnImim_d34 = _mm256_xor_pd(Bimim_d34,(const v4df)signflip);
      v4df CX_a0_imre_d34;
      CX_a0_imre_d34 = _mm256_permute_pd(CX_a0_d34,5);
      res_a0_d34 += CX_a0_imre_d34 * BnImim_d34;
   }
   _mm256_storeu_pd((double *)(res    ), res_a0);
   _mm256_storeu_pd((double *)(res + 2), res_a0_d34);
}

namespace {
   inline int min (int a, int b) { return (a < b) ? a : b; }
}

template<int ba, int bb, int bd>
void pgemm3m_avx(int Na, int Nb, int Nc, int Nd,
            const SxComplex16 *B, int ldb,
            const double *C, int ldc,
            const SxComplex16 *X, int xldc, int xldbc,
                  SxComplex16 *__restrict__ res, int lda)
{
   using namespace std;
   if (Na <= 0 || Nb <= 0 || Nc <= 0 || Nd <= 0) return;
   const int ba2 = 12;
   const int bb2 = 16;
   const int bb3 = 15;
   const int bd2 = (Nb <= 24) ? 128 : 24;
   int Nab = ba * (Na / ba);
   int Nbb = bb * (Nb / bb);
   int Ndb = bd * (Nd / bd);
   int Nb3 = 3 * (Nb / 3);
   void *mem = NULL;
   int err = posix_memalign(&mem, sizeof(v4df), sizeof(v4df) * Nc * bd2/2);
   if (err) { SX_EXIT; }
   v4df *C1324 = (v4df*)mem;
   err = posix_memalign(&mem, sizeof(v4df), sizeof(v4df) * Nb * bd2);
   if (err) { SX_EXIT; }
   v4df *Bdup = (v4df*)mem;
   for (int d2 = 0; d2 < Ndb; d2+=bd2)  {
      for (int d = 0; d < bd2 && d2 + d < Ndb; d+=8)
         gemm3mPrepareCdup(C + d2 + d, ldc, C1324 + d/2 * Nc, Nc);
      //if (Nab>0)
      gemm3mPrepareBdup(B + d2, ldb, Bdup, Nb, min(bd2, Ndb - d2));
      for (int a2 = 0; a2 < Nab; a2+=ba2)  {

         for (int b2 = 0; b2 < Nbb; b2+=bb2)  {
            for (int a = 0; a < ba2 && a2 + a < Nab; a+=ba)  {
               v4df CX[12 * bb2];
               for (int d = 0; d < bd2 && d2 + d < Ndb; d+=bd)  {
                  /* prefetching seems not a good idea
                  // --- L1 cache (write)
                  _mm_prefetch(res + d2 + d + lda * (a2 + a), _MM_HINT_ET0);
                  _mm_prefetch(res + d2 + d + 4 + lda * (a2 + a), _MM_HINT_ET0);
                  _mm_prefetch(res + d2 + d + lda * (a2 + a + 1), _MM_HINT_ET0);
                  _mm_prefetch(res + d2 + d + 4 + lda * (a2 + a + 1), _MM_HINT_ET0);
                  _mm_prefetch(res + d2 + d + lda * (a2 + a + 2), _MM_HINT_ET0);
                  _mm_prefetch(res + d2 + d + 4 + lda * (a2 + a + 2), _MM_HINT_ET0);
                  */
                  for (int b = 0; b < bb2 && b2 + b < Nbb; b+=bb)  {
                     gemm3m_avx_CX_3_1_8(Nc, C1324 + d/2 * Nc,
                        X + xldc * (b2 + b) + xldbc * (a2 + a), xldbc,
                        CX + 12 * b);
                     // fetch into L1 cache (read) -- seems not to be a good idea
                     // _mm_prefetch(B + d2 + d + ldb * (b2 + b), _MM_HINT_T0);
                     // _mm_prefetch(B + d2 + d + 4 + ldb * (b2 + b), _MM_HINT_T0);
                  }
                  gemm3m_avx_CXmulBd_3_1_4(min(bb2, Nbb-b2), CX, Bdup + 4 * b2 + Nb * d,
                                    res + d2 + d + lda * (a2 + a), lda);
                  gemm3m_avx_CXmulBd_3_1_4(min(bb2, Nbb-b2), CX + 6, Bdup + 4 * b2 + Nb*(d + 4),
                                    res + d2 + d + 4 + lda * (a2 + a), lda);
               }
            }
         }
      }
      // --- clean-up a
      for (int a = Nab; a < Na; a++)  {
         v4df CX[4 * bb3];
         for (int b3 = 0; b3 < Nb3; b3+=bb3)  {
            for (int d = 0; d < bd2 && d2 + d < Ndb; d+=bd)  {
               for (int b = 0; b < bb3 && b3 + b < Nb3; b+=3)  {
		  /*
		  if (Nab == 0)  {
		     gemm3mPrepareBdup(B + d2 + d + ldb * (b3 + b), ldb, Bdup + 4 * (b3 + b) + Nb * d, 3, 4);
		     gemm3mPrepareBdup(B + d2 + d + 4 + ldb * (b3 + b), ldb, Bdup + 4 * (b3 + b) + Nb * (d + 4), 3, 4);
		  }
		  */
                  gemm3m_avx_CX_3_1_8(Nc, C1324 + d/2 * Nc,
                        X + xldc * (b3 + b) + xldbc * a, xldc,
                        CX + 4 * b);
               }
               gemm3m_avx_CXmulBd_1_3_4(min(bb3, Nb3-b3), CX, Bdup + 4 * b3 + Nb * d,
                     res + d2 + d + lda * a);
               gemm3m_avx_CXmulBd_1_3_4(min(bb3, Nb3-b3), CX + 6, Bdup + 4 * b3 + Nb * (d + 4),
                     res + d2 + (d + 4) + lda * a);
            }
         }
         // clean-up b within clean-up a
         if (Nb3 < Nb)  {
            for (int d = 0; d < bd2 && d2 + d < Ndb; d+=bd)  {
               for (int b = Nb3; b < Nb; b++)  {
                  gemm3m_avx_CX_1_1_8(Nc, C1324 + d/2 * Nc,
                        X + xldc * b + xldbc * a,
                        CX + 4 * (b-Nb3));
               }
               /*
               gemm3m_avx_CXmulBd_1_1_4(Nb-Nb3, CX, Bdup + 4 * Nb3 + Nb * d,
                     res + d2 + d + lda * a, lda);
               gemm3m_avx_CXmulBd_1_1_4(Nb-Nb3, CX + 2, Bdup + 4 * Nb3 + Nb*(d + 4),
                     res + d2 + d + 4 + lda * a, lda);
                     */
               gemm3m_avx_CXmulB_1_1_4(Nb-Nb3, CX, B + d2 + d + Nb3 * ldb, ldb,
                     res + d2 + d + lda * a, lda);
               gemm3m_avx_CXmulB_1_1_4(Nb-Nb3, CX + 2, B + d2 + d + 4 + Nb3 * ldb, ldb,
                     res + d2 + d + 4 + lda * a, lda);
            }
         }
      }
   }
   // clean-up d
   pgemm3m<SxComplex16, double, SxComplex16> (Na, Nb, Nc, Nd - Ndb,
            B + Ndb, ldb,
            C + Ndb, ldc,
            X, xldc, xldbc,
            res + Ndb, lda);
   free (Bdup);
   free (C1324);
}

#endif

template<class T, int bb, int bd>
class Block
{
   public:
      // Block (on stack)
      T Bb[bb * bd];
      // cast to pointer
      operator const T* () { return Bb; }
      // initializer
      Block (const T *__restrict__ Bbegin, int ldb)
      {
         copy<bb,bd>(Bbegin, ldb, Bb);
      }
};

// special case bb=1: do not copy
template<class T, int bd>
class Block<T, 1, bd>
{
   public:
      const T* Bb;
      operator const T* () { return Bb; }
      Block (const T*__restrict BbIn, int) : Bb(BbIn)
      {
      }
};

// --- block driver routine
#if (USE_AVX == 1)
#define GEMMM_B_KERNEL gemmm_b_kernel_avx_12_1_bd<bd>
#define GEMMM_B_KERNEL1 gemmm_b_kernel_avx_1_1_bd<bd>
#define PREPARE_A reorder_ri
#define PREPARE1 reorder4_ri
#else
#define GEMMM_B_KERNEL gemmm_b_kernel2<ba,bb,bd,TA,TB,TC>
#define GEMMM_B_KERNEL1 gemmm_b_kernel2<1,bb,bd,TA,TB,TC>
#define PREPARE_A copy
#define PREPARE1 copy
#endif
template<int ba, int bb, int bd, class TA, class TB, class TC>
void gemmm_b3 (int Na, int Nb, int Nc, int Nd,
             const TA *__restrict__ A, int lda,
             const TB *__restrict__ B, int ldb,
             const TC *__restrict__ C, int ldc,
             typename TM<TA, typename TM<TB,TC>::Res>::Res *__restrict__ res,
             int rldc, int rldbc)
{
   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         for (int ic = 0; ic < Nc; ++ic)  {
            res[ic + rldc * ib + rldbc * ia] = 0.;
         }
      }
   }
   const int ba2 = 4*ba;
   int Nab = ba * (Na / ba);
   int Nab2 = ba2 * (Na / ba2);
   int Nbb = bb * (Nb / bb);
   //int Ndb = bd * (Nd / bd);
   const int bd2 = 1*bd;
   int Ndb2 = bd2 * (Nd / bd2);
   for (int id = 0; id < Ndb2; id += bd2)  {
      for (int ia = 0; ia <= Nab2; ia += ba2)  {
         TA Ab[ba2 * bd2];
         for (int ja = 0; ja < ba2 && ia + ja < Nab; ja += ba)
            for (int jd = 0; jd < bd2; jd += bd)
               PREPARE_A<ba,bd>(A + id + jd + (ia + ja)*lda, lda,
                           Ab + jd*ba + ja * bd2);
         for (int ib = 0; ib < Nbb; ib += bb)  {
            for (int ja = 0; ja < ba2 && ia + ja < Nab; ja += ba)  {
               for (int jd = 0; jd < bd2; jd += bd)  {
                  // copy block if bb>1
                  Block<TB,bb,bd> Bb(B + id + jd + ib*ldb, ldb);
                  GEMMM_B_KERNEL(Nc, Ab + ba * jd + bd2 * ja, Bb,
                                 C + id + jd, ldc,
                                 res + rldc * ib + rldbc * (ia + ja),
                                 rldc, rldbc);
               }
            }
         }
      }
      // cleanup a
      if (Nab < Na)  {
         TA Ab[ba * bd2];
         for (int ja = 0; ja < Na-Nab; ++ja)
            for (int jd = 0; jd < bd2; jd += bd)
               PREPARE1<1,bd>(A + id + jd + (Nab + ja)*lda, lda,
                              Ab + jd + bd2 * ja);
         for (int ib = 0; ib < Nbb; ib += bb)  {
            TB Bb[bb * bd2];
            for (int jd = 0; jd < bd2; jd += bd)
               PREPARE1<bb,bd>(B + id + jd + ib*ldb, ldb, Bb + jd*bb);
            for (int ia = Nab; ia < Na; ia++)  {
               for (int jd = 0; jd < bd2; jd += bd)  {
                  GEMMM_B_KERNEL1(Nc, Ab + jd + (ia - Nab) * bd2, Bb + jd * bb,
                                  C + id + jd, ldc,
                                  res + rldc * ib + rldbc * ia, rldc, rldbc);
               }
            }
         }
      }
   }
   // cleanup b
   pgemmm(Na, Nb-Nbb, Nc, Ndb2,
         A, lda, B + Nd * Nbb, ldb, C, ldc,
         res + rldc * Nbb       , rldc, rldbc);
   // cleanup d
   pgemmm(Na, Nb, Nc, Nd-Ndb2,
         A + Ndb2, lda, B + Ndb2, ldb, C + Ndb2, ldc,
         res                    , rldc, rldbc);
}

/*
// KernelA       : large ba, small bb, blocking on d (bd)
// KernelB       : small ba, large bb, blocking on d (bd)
// CleanupKernel : ba=1, bb=1, blocking on d (bd)
template<class KernelA, class KernelB, class CleanupKernel>
void gemmm_driver (int Na, int Nb, int Nc, int Nd,
             const SxComplex16 *__restrict__ A, int lda,
             const SxComplex16 *__restrict__ B, int ldb,
             const double      *__restrict__ C, int ldc,
                   SxComplex16 *__restrict__ res,
             int rldc, int rldbc)
{
   SX_CHECK_VARS (int(KernelA::bd) == int(KernelB::bd),
	          KernelA::bd, KernelB::bd);
   SX_CHECK_VARS (int(KernelA::bd) == int(CleanupKernel::bd),
	          KernelA::bd, CleanupKernel::bd);

   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         for (int ic = 0; ic < Nc; ++ic)  {
            res[ic + rldc * ib + rldbc * ia] = 0.;
         }
      }
   }

   const int bd = KernelA::bd;
   // enhance bd block size to 64 for outermost d-loop
   const int bd2 = max(64 / bd, 1)*bd;
   // use up to 48 KB for prepared A data
   const int ba2 = max(48 * 1024 / 16 / bd2 / KernelA::ba,1) * KernelA::ba;

   int Nab = KernelA::ba * (Na / KernelA::ba);
   int Nab2 = ba2 * (Na / ba2);
   int Nbb = KernelA::bb * (Nb / KernelA::bb);

   int Ndb2 = bd2 * (Nd / bd2);

   int Nab_B = Nab + KernelB::ba * ((Na-Nab) / KernelB::ba);
   int Nbb_B = KernelB::bb * (Nb / KernelB::bb);
   for (int id = 0; id < Ndb2; id += bd2)  {
      for (int ia = 0; ia <= Nab2; ia += ba2)  {
         TA Ab[ba2 * bd2];
         for (int ja = 0; ja < ba2 && ia + ja < Nab; ja += ba)
            for (int jd = 0; jd < bd2; jd += bd)
               KernelA::prepareA (A + id + jd + (ia + ja)*lda, lda,
                           Ab + jd*ba + ja * bd2);
         for (int ib = 0; ib < Nbb; ib += bb)  {
	    for (int jd = 0; jd < bd2; jd += bd)  {
	       // copy block if bb>1
	       Block<SxComplex16,bb,bd> Bb(B + id + jd + ib*ldb, ldb);
	       for (int ja = 0; ja < ba2 && ia + ja < Nab; ja += ba)  {
		  KernelA::run(Nc, Ab + ba * jd + bd2 * ja, Bb,
                               C + id + jd, ldc,
                               res + rldc * ib + rldbc * (ia + ja),
                               rldc, rldbc);
               }
            }
         }
	 // cleanup b for Kernel A
         for (int ib = Nbb; ib < Nb; ib++)  {
	    for (int jd = 0; jd < bd2; jd += bd)  {
	       for (int ja = 0; ja < ba2 && ia + ja < Nab; ja += ba)  {
		  KernelA::cleanupB(Nc, Ab + ba * jd + bd2 * ja,
			            B + id + jd + ib*ldb,
                                    C + id + jd, ldc,
                                    res + rldc * ib + rldbc * (ia + ja),
                                    rldc, rldbc);
               }
            }
	 }
      }
      // use B kernel for partial A cleanup (small ba, large bb) 
      if (Nab < Nab_B) {
	 SxComplex16 Ab[KernelA::ba * bd2];
	 int Na_B = Nab_B - Nab;
	 for (int ia = Nab; ia < Nab_B; ia += KernelB::ba)  {
	    for (int jd = 0; jd < bd2; jd += bd)
	       KernelB::prepareA (A + id + jd + ia*lda, lda,
		                  Ab + jd*KernelB::ba + (ia-Nab)*bd2);
	 }
      	 for (int ib = 0; ib < Nbb_B; ib += KernelB::bb)  {
            TB Bb[KernelB::bb * bd2];
            for (int jd = 0; jd < bd2; jd += bd)
               KernelB::prepareB(B + id + jd + ib*ldb, ldb, Bb + jd*KernelB::bb);
	    for (int jd = 0; jd < bd2; jd += bd)  {
	       for (int ja = 0; ja < Na_B; ja += KernelB::ba)  {
		  KernelB::run(Nc, Ab + KernelB::ba * jd + bd2 * ja, Bb + jd*KernelB::bb,
                               C + id + jd, ldc,
                               res + rldc * ib + rldbc * (Nab + ja),
                               rldc, rldbc);
               }
            }
         }
	 if (Nbb_B < Nb)  {
	    // --- cleanup b for Kernel B using CleanupKernel
	    for (int ja = 0; ja < Nab_B-Nab; ++ja)
	       for (int jd = 0; jd < bd2; jd += bd)
		  CleanupKernel::prepareA(A + id + jd + (Nab + ja)*lda, lda,
                                          Ab + jd + bd2 * ja);
	    for (int ib = Nbb_B; ib < Nbb; ib++)  {
	       TB Bb[bd2];
	       for (int jd = 0; jd < bd2; jd += bd)
		  CleanupKernel::prepareB(B + id + jd + ib*ldb, ldb, Bb + jd*bb);
	       for (int ia = Nab; ia < Nab_B; ia++)  {
		  for (int jd = 0; jd < bd2; jd += bd)  {
		     CleanupKernel::run(Nc, Ab + jd + (ia - Nab) * bd2, Bb + jd * KernelB::bb,
                                        C + id + jd, ldc,
                                        res + rldc * ib + rldbc * ia);
               }
            }
         }
      }
      // --- cleanup a using CleanupKernel
      if (Nab_B < Na)  {
         TA Ab[KernelB::ba * bd2];
         for (int ja = 0; ja < Na-Nab_B; ++ja)
            for (int jd = 0; jd < bd2; jd += bd)
               CleanupKernel::prepareA(A + id + jd + (Nab_B + ja)*lda, lda,
                                       Ab + jd + bd2 * ja);
         for (int ib = 0; ib < Nb; ib++)  {
            TB Bb[bd2];
            for (int jd = 0; jd < bd2; jd += bd)
               CleanupKernel::prepareB(B + id + jd + ib*ldb, ldb, Bb + jd*bb);
            for (int ia = Nab_B; ia < Na; ia++)  {
               for (int jd = 0; jd < bd2; jd += bd)  {
		  CleanupKernel::run(Nc, Ab + jd + (ia - Nab_B) * bd2, Bb + jd * KernelB::bb,
                                  C + id + jd, ldc,
                                  res + rldc * ib + rldbc * ia);
               }
            }
         }
      }
   }
   // cleanup d
   pgemmm(Na, Nb, Nc, Nd-Ndb2,
         A + Ndb2, lda, B + Ndb2, ldb, C + Ndb2, ldc,
         res                    , rldc, rldbc);
}
*/


#if (USE_AVX == 1)
#ifndef GEMM3MKERNEL
#  define GEMM3MKERNEL gemm3m_kernel_avx_1_2_12
#endif
#ifndef GEMM3MKERNEL1
#  define GEMM3MKERNEL1 gemm3m_kernel_avx_1_1_12
#endif
#else
#define GEMM3MKERNEL gemm3m_kernel<ba,bb,bd,TB,TC,TD>
#define GEMM3MKERNEL1 gemm3m_kernel<1,1,bd,TB,TC,TD>
#endif

template<int ba, int bb, int bd, class TB, class TC, class TD>
void pgemm3m_b(int Na, int Nb, int Nc, int Nd,
            const TB *B, int ldb,
            const TC *C, int ldc,
            const TD *X, int xldc, int xldbc,
                  typename TM<TD, typename TM<TB,TC>::Res>::Res *__restrict__ res, int lda)
{
   if (Na <= 0 || Nb <= 0 || Nc <= 0 || Nd <= 0) return;
   int Nab = ba * (Na / ba);
   int Nbb = bb * (Nb / bb);
   int Ndb = bd * (Nd / bd);
   for (int d = 0; d < Ndb; d+=bd)  {
      for (int a = 0; a < Nab; a+=ba)  {
         for (int b = 0; b < Nbb; b+=bb)  {
            GEMM3MKERNEL
               (Nc, B + d + ldb * b, ldb, C + d, ldc,
                X + xldc * b + xldbc * a, xldc, xldbc,
                res + d + lda * a, lda);
         }
         // clean-up b
         for (int b = Nbb; b < Nb; ++b) {
            GEMM3MKERNEL1
               (Nc, B + d + ldb * b, ldb, C + d, ldc,
                X + xldc * b + xldbc * a, xldc, xldbc,
                res + d + lda * a, lda);
         }
      }
      // clean-up a
      pgemm3m<TB,TC,TD> (Na - Nab, Nb, Nc, bd,
               B + d, ldb,
               C + d, ldc,
               X + xldbc * Nab, xldc, xldbc,
               res + d + lda * Nab, lda);
   }
   // clean-up d
   pgemm3m<TB,TC,TD> (Na, Nb, Nc, Nd - Ndb,
            B + Ndb, ldb,
            C + Ndb, ldc,
            X, xldc, xldbc,
            res + Ndb, lda);
}

// --- omp driver (parallelize over a)
#ifdef USE_OPENMP
template<int ba, int bb, int bd, class TA, class TB, class TC>
void gemmm_b2_omp (int Na, int Nb, int Nc, int Nd,
             const TA *__restrict__ A, int lda,
             const TB *__restrict__ B, int ldb,
             const TC *__restrict__ C, int ldc,
             typename TM<TA, typename TM<TB,TC>::Res>::Res *__restrict__ res,
             int rldc, int rldbc)
{
#pragma omp parallel
   {
      int nt = omp_get_num_threads ();
      int nta = nt, ntb = 1;
      int tba, tbb;
      int rateBlock = 15, rateClean = 10; // approximate numbers
      do {
         tba = Na / nta;
         if (tba * nta < Na) tba++;
         tbb = Nb / ntb;
         if (tbb * ntb < Nb) tbb++;

         // cannot exploit a-blocking anyway, parallelize a
         if (Na < ba && Na >= nta) break;

         // round tba up to next a-block size, if worthwhile
         if ((tba % ba) * rateBlock > ba * rateClean && tbb >= bb)
            tba = (tba / ba + 1) * ba;
         // round tbb up to next b-block size, if worthwhile
         if ((tbb % bb) * rateBlock > bb * rateClean && tba >= ba)
            tbb = (tbb / bb + 1) * bb;
         if (nta == 1) break; // nothing to parallelize over a


         // no significant loss in efficiency due to a-parallelism
         if ((tba % ba ) * 4 < (tba / ba) * ba) break;

         // parallelize more aggressively over b
         for (ntb++; ntb < nt; ++ntb)
            if (nt % ntb == 0) break;
         nta = nt / ntb;
      } while (true);

      int iThread = omp_get_thread_num ();
      int iba = iThread / ntb;
      int ibb = iThread % ntb;

      /*
      if (iThread == 0) {
         std::cout << "nta = " << nta << " ntb = " << ntb << std::endl;
         std::cout << "tba = " << tba << " tbb = " << tbb << std::endl;
      }
      */

      int Nat = tba;
      if (iba * tba + Nat > Na) Nat = Na - iba * tba;
      if (iba * tba >= Na) Nat = 0;

      int Nbt = tbb;
      if (ibb * tbb + Nbt > Nb) Nbt = Nb - ibb * tbb;
      if (ibb * tbb >= Nb) Nbt = 0;

      gemmm_b3<ba,bb,bd,TA,TB,TC>(Nat, Nbt, Nc, Nd,
                                 A + lda * iba * tba, lda,
                                 B + ldb * ibb * tbb, ldb,
                                 C, ldc,
                                 res + rldc * ibb * tbb + rldbc * iba * tba,
                                 rldc, rldbc);
   }
}

// --- omp driver (parallelize over d)
template<int ba, int bb, int bd, class TB, class TC, class TD>
void pgemm3m_b_omp(int Na, int Nb, int Nc, int Nd,
            const TB *B, int ldb,
            const TC *C, int ldc,
            const TD *X, int xldc, int xldbc,
            typename TM<TD, typename TM<TB,TC>::Res>::Res *res, int lda)
{
#pragma omp parallel
   {
      int nt = omp_get_num_threads ();
      int tbd = Nd / nt;
      if (tbd * nt < Nd) tbd++;
      int blockSpeed = 30; // approximate speed for blocked algorithm
      int cleanSpeed = 1; // approximate speed for cleanup algorithm
      if ((tbd % bd) * blockSpeed >= bd * cleanSpeed) {
         // increase tbd to next multiple of block size
         // because this is faster than doing the cleanup
         int tbdNew = bd * (tbd / bd + 1);
         // check workload for last thread (needs to do cleanup)
         int lastThread = tbd % bd - (tbdNew - tbd) * nt;
         int lastClean = lastThread % bd;
         if (lastClean < 0) lastClean += bd;
         int lastBlock = lastThread - lastClean;
         // check that we do not increase time on final thread beyond savings
         if ( (tbd % bd - lastClean) * blockSpeed > lastBlock * cleanSpeed)
            tbd = tbdNew;
      }

      int ibd = omp_get_thread_num ();
      int Ndt = tbd;
      if (ibd * tbd + Ndt > Nd) Ndt = Nd - ibd * tbd;
      if (ibd * tbd >= Nd) Ndt = 0;
#if (USE_AVX == 1)
      pgemm3m_avx<ba,bb,bd>
#else
      pgemm3m_b<ba,bb,bd>
#endif
         (Na, Nb, Nc, Ndt,
                         B + ibd * tbd, ldb, C + ibd * tbd, ldc,
                         X, xldc, xldbc,
                         res + ibd * tbd, lda);
   }
}
#endif

} // anonymous namespace

#ifdef SXGEMMM_TUNE

#if (USE_AVX == 1)
#  define BA 12
#  define BB 1
#  ifndef BD
#    define BD 64
#  endif
# if (BD % 16 != 0)
#error BD block size must be multiple of 16
#endif
#endif

// --- the simple algorithm
void gemmm (int Na, int Nb, int Nc, int Nd,
           CRTYPE_A *__restrict__ A, int lda,
           CRTYPE_B *__restrict__ B, int ldb,
           CRTYPE_C *__restrict__ C, int ldc,
           CRTYPE_D *__restrict__ res, int rldc, int rldbc)
{
   for (int ia = 0; ia < Na; ++ia)  {
      for (int ib = 0; ib < Nb; ++ib)  {
         for (int ic = 0; ic < Nc; ++ic)  {
            CRTYPE_D S = 0.;
            for (int id = 0; id < Nd; ++id)  {
               S += A[id + lda * ia]
                  * B[id + ldb * ib]
                  * C[id + ldc * ic];
            }
            res[ic + rldc * ib + rldbc * ia] = S;
         }
      }
   }
}

#elif defined SXGEMM3M_TUNE
void gemm3m(int Na, int Nb, int Nc, int Nd,
            const CRTYPE_B *B, int ldb,
            const CRTYPE_C *C, int ldc,
            const CRTYPE_D *X, int xldc, int xldbc,
                  CRTYPE_A *res, int lda)
{
   for (int a = 0; a < Na; ++a)
      for (int d = 0; d < Nd; ++d)
         res[d + lda * a] = 0.;
   for (int a = 0; a < Na; ++a)  {
      for (int b = 0; b < Nb; ++b)  {
         for (int d = 0; d < Nd; ++d)  {
            TM<CRTYPE_D,CRTYPE_C>::Res S = 0.;
            for (int c = 0; c < Nc; ++c)  {
               S += C[d + ldc * c] * X[c + xldc * b + xldbc * a];
            }
            res[d + lda * a] += B[d + ldb * b] * S;
         }
      }
   }
}

#else
// --- here comes the externally visible interface

#if (USE_AVX == 1)
#  define SXGEMMM_BA 12
#  define SXGEMMM_BB 1
#  ifndef SXGEMMM_BD
#    define SXGEMMM_BD 64
#  else
#    if (SXGEMMM_BD % 16 != 0)
#       error "BD block size must be multiple of 16"
//    reason: use of gemmm_b_kernel_avx1_1_16 cleanup kernel
#    endif
#  endif
#  define SXGEMM3M_BA 3
#  define SXGEMM3M_BB 1
#  define SXGEMM3M_BD 8
#else
// --- blocking parameters optimal for my Intel Xeon CPU E5-2667 0 @ 2.90GHz
#  if ! defined SXGEMMM_BA
#    define SXGEMMM_BA 2
#  endif
#  if ! defined SXGEMMM_BB
#    define SXGEMMM_BB 1
#  endif
#  if ! defined SXGEMMM_BD
#    define SXGEMMM_BD 4
#  endif
#  if ! defined SXGEMM3M_BA
#    define SXGEMM3M_BA 1
#  endif
#  if ! defined SXGEMM3M_BB
#    define SXGEMM3M_BB 1
#  endif
#  if ! defined SXGEMM3M_BD
#    define SXGEMM3M_BD 1
#  endif
#endif

#ifdef USE_OPENMP
#define GEMMM gemmm_b2_omp
#define PGEMM3M pgemm3m_b_omp <SXGEMM3M_BA,SXGEMM3M_BB,SXGEMM3M_BD, SxComplex16, double, SxComplex16>
#else
#define GEMMM gemmm_b3
#if (USE_AVX==1)
#define PGEMM3M pgemm3m_avx <SXGEMM3M_BA,SXGEMM3M_BB,SXGEMM3M_BD>
#else
#define PGEMM3M pgemm3m_b <SXGEMM3M_BA,SXGEMM3M_BB,SXGEMM3M_BD, SxComplex16, double, SxComplex16>
#endif
#endif
void sxgemmm (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *__restrict__ A, ssize_t lda,
              const SxComplex16 *__restrict__ B,
              const double      *__restrict__ C,
              SxComplex16 *__restrict__ res, ssize_t rldc, ssize_t rldbc)
{
#ifndef NDEBUG
//#pragma GCC reset_options
   // check validity of input as we disable nan-handling by optimization
   for (ssize_t i = 0; i < Na * lda; ++i)
      SX_CHECK_NUM(A[i]);
   for (ssize_t i = 0; i < Nb * Nd; ++i)
      SX_CHECK_NUM(B[i]);
   for (ssize_t i = 0; i < Nc * Nd; ++i)
      SX_CHECK_NUM(C[i]);
#endif
   GEMMM <SXGEMMM_BA,SXGEMMM_BB,SXGEMMM_BD,SxComplex16,SxComplex16,double>
      (int(Na), int(Nb), int(Nc), int(Nd), A, (int)lda, B, int(Nd), C, int(Nd),
       res, int(rldc), int(rldbc));
#ifndef NDEBUG
   for (ssize_t a = 0; a < Na; ++a)
      for (ssize_t b = 0; b < Nb; ++b)
         for (ssize_t c = 0; c < Nc; ++c)
            SX_CHECK_NUM(res[c + rldc * b + rldbc * a]);
#endif
}

void sxpgemm3m(ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *B,
              const double      *C,
              const SxComplex16 *X, ssize_t xldc, ssize_t xldbc,
                    SxComplex16 *res)
{
#ifndef NDEBUG
//#pragma GCC reset_options
   // check validity of input as we disable nan-handling by optimization
   for (ssize_t i = 0; i < Nb * Nd; ++i)
      SX_CHECK_NUM(B[i]);
   for (ssize_t i = 0; i < Nc * Nd; ++i)
      SX_CHECK_NUM(C[i]);
   for (ssize_t a = 0; a < Na; ++a)
      for (ssize_t b = 0; b < Nb; ++b)
         for (ssize_t c = 0; c < Nc; ++c)
            SX_CHECK_NUM(X[c + xldc * b + xldbc * a]);
#endif
   PGEMM3M
      (int(Na), int(Nb), int(Nc), int(Nd), B, int(Nd), C, int(Nd),
       X, int(xldc), int(xldbc), res, int(Nd));
#ifndef NDEBUG
   for (ssize_t i = 0; i < Na * Nd; ++i)
      SX_CHECK_NUM(res[i]);
#endif
}
#endif /* SXGEMMM_TUNE */
#endif /* GNUG */

#endif /* WIN32 */
