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

#ifndef SX_EXPORT_DIRAC
#ifdef WIN32
#  if defined(_EXPORT_sxdirac)
#     define SX_EXPORT_DIRAC __declspec(dllexport)
#  else
#     define SX_EXPORT_DIRAC __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_DIRAC
#endif
#endif


