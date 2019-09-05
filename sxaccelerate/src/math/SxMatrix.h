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

#include <SxVector.h>

#ifndef _SX_MATRIX_H_
#define _SX_MATRIX_H_

#ifdef SXMATBASIS
#   undef SXMATBASIS
#endif
#define SXMATBASIS  SxVector

#ifdef SXMAT
#   undef SXMAT
#endif
#define SXMAT       SxMatrix


#ifdef SXSYMMAT
#   undef SXSYMMAT
#endif
#define SXSYMMAT    SxSymMatrix

#ifdef MEMHANDLE
#   undef MEMHANDLE
#endif
#define MEMHANDLE    SxVectorHandle

#include <SxMat.h>

template <class T>
class SxVecTraits<SxVector<T> > 
{
   public:
   typedef T ScalarType;
   typedef SxMatrix<T> MatType;
};

#endif /* _SX_MATRIX_H_ */
