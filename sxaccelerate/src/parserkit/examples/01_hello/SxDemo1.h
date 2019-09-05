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

#ifndef _SX_DEMO_1_H_
#define _SX_DEMO_1_H_


#ifdef WIN32
#  if defined(_EXPORT_sxdemo1)
#     define SX_EXPORT_DEMO_1 __declspec(dllexport)
#  else
#     define SX_EXPORT_DEMO_1 __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_DEMO_1
#endif

#endif /* _SX_DEMO_1_H_ */
