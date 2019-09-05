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

#ifndef _SX_MATH_H_
#define _SX_MATH_H_


#ifdef WIN32
#  if defined(_EXPORT_sxmath)
#     define SX_EXPORT_MATH __declspec(dllexport)
#  else
#     define SX_EXPORT_MATH __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_MATH  
#endif


#endif /* _SX_MATH_H_ */
