// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_STRUCT_H_
#define _SX_STRUCT_H_


#ifdef WIN32
#  if defined(_EXPORT_sxstruct)
#     define SX_EXPORT_STRUCT __declspec(dllexport)
#  else
#     define SX_EXPORT_STRUCT __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_STRUCT
#endif

#endif /* _SX_STRUCT_H_ */
