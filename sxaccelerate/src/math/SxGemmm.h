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

#ifndef _SX_GEMMM_H_
#define _SX_GEMMM_H_

#include <SxMath.h>
#include <SxComplex.h>
#include <SxConfig.h>

#ifndef WIN32
SX_EXPORT_MATH void sxgemmm (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *A, ssize_t lda,
              const SxComplex16 *B,
              const double      *C,
              SxComplex16 *res, ssize_t rldc, ssize_t rldbc);

// wrapper for lda = Nd
inline
void sxgemmm (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *__restrict__ A,
              const SxComplex16 *__restrict__ B,
              const double      *__restrict__ C,
              SxComplex16 *__restrict__ res, ssize_t rldc, ssize_t rldbc)
{
   sxgemmm(Na, Nb, Nc, Nd, A, Nd, B, C, res, rldc, rldbc);
}

SX_EXPORT_MATH void sxpgemm3m(ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *B,
              const double      *C,
              const SxComplex16 *X, ssize_t xldc, ssize_t xldbc,
                    SxComplex16 *res);
#endif

#endif /* _SX_GEMMM_H_ */
