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
// This is a random namber generator ("universal")
// Source: R. Wieczorkowski, R. Zielincki "Komputerowe generatory
// licyb losowych" WNT, Warszawa 1997, page 39-40
//----------------------------------------------------------------------------

#ifndef _SX_RAND_H_
#define _SX_RAND_H_

#include <SxMath.h>

void initRandom (int i, int j, int k, int l);

class SX_EXPORT_MATH SxRandom
{
   public:
      SxRandom () {
         ip = 1;
      }
      static double uu[97];
      static long int ip;
      static long int jp;
      static double cc;
      static double cd;
      static double cm;

      static double get ();

};

#endif // _SX_RAND_H_
