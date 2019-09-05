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

#ifndef _SX_SYM_MATRIX_H_
#define _SX_SYM_MATRIX_H_

#include <SxVector.h>
#include <SxMatrix.h>
#include <SxPrecision.h>

#ifdef SXMATBASIS
#   undef SXMATBASIS
#endif
#define SXMATBASIS  SxVector

#ifdef SXSYMMAT
#   undef SXSYMMAT
#endif
#define SXSYMMAT    SxMatrix


#ifdef SXSYMMAT
#   undef SXSYMMAT
#endif
#define SXSYMMAT    SxSymMatrix

#ifdef MEMHANDLE
#   undef MEMHANDLE
#endif
#define MEMHANDLE    SxVectorHandle

#include <SxSymMat.h>

#endif /* _SX_SYM_MATRIX_H_ */
