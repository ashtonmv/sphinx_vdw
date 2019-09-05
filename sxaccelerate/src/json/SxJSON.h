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

#ifndef _SX_JSON_H_
#define _SX_JSON_H_


#ifdef WIN32
#  if defined(_EXPORT_sxjson)
#     define SX_EXPORT_JSON __declspec(dllexport)
#  else
#     define SX_EXPORT_JSON __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_JSON
#endif

#endif /* _SX_JSON_H_ */
