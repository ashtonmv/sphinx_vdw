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

#ifndef _SX_NON_STD_H_
#define _SX_NON_STD_H_


#ifdef WIN32
#  if defined(_EXPORT_sxnonstd)
#     define SX_EXPORT_NONSTD __declspec(dllexport)
#  else
#     define SX_EXPORT_NONSTD __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_NONSTD
#endif

#endif /* _SX_NON_STD_H_ */
