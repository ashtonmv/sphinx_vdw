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

#include <SxInfNan.h>
#include <SxConfig.h>
#include <math.h>
#include <stdio.h>

#ifdef WIN32
   bool sxIsNan (int)       { return false; }
   bool sxIsNan (float)     { return false; }
   bool sxIsNan (double)    { return false; }
   bool sxIsInf (int)       { return false; }
   bool sxIsInf (float)     { return false; }
   bool sxIsInf (double)    { return false; }
#elif defined (MACOSX)
   bool sxIsNan (int i)     { return isnan(static_cast<float>(i)); }
   bool sxIsNan (float i)   { return isnan(i); }
   bool sxIsNan (double i)  { return isnan(i); }
   bool sxIsInf (int i)     { return isinf(static_cast<float>(i)); }
   bool sxIsInf (float i)   { return isinf(i); }
   bool sxIsInf (double i)  { return isinf(i); }
#else
   // --- isnan and isinf macros produce -Wconversion warnings
   //     missing cast to float for __isnanf() in math.h on Debian, CentOS
   bool sxIsNan (int i)     { return __isnanf(static_cast<float>(i)); }
   bool sxIsNan (float i)   { return __isnanf(i); }
   bool sxIsNan (double i)  { return __isnan(i); }
   bool sxIsInf (int i)     { return __isinff(static_cast<float>(i)); }
   bool sxIsInf (float i)   { return __isinff(i); }
   bool sxIsInf (double i)  { return __isinf(i); }
#endif

