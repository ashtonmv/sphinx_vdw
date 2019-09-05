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

#ifndef _SX_AVX_H_
#define _SX_AVX_H_

#ifdef __GNUG__
#  if (! defined USE_AVX) && defined __AVX__
#    define USE_AVX 1
#  endif
#endif

#endif /* _SX_AVX_H_ */
