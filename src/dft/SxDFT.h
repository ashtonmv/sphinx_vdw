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

#ifndef _SX_DFT_H_
#define _SX_DFT_H_


#ifdef WIN32
#  if defined(_EXPORT_sxdft)
#     define SX_EXPORT_DFT __declspec(dllexport)
#  else
#     define SX_EXPORT_DFT __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_DFT
#endif

#endif /* _SX_DFT_H_ */
