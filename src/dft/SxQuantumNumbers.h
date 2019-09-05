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

#ifndef _SX_QNUMBERS_H_
#define _SX_QNUMBERS_H_

#include <SxError.h>
#include <stdio.h>
#include <iostream>

class SxQuantumNumbers;
std::ostream &operator<< (std::ostream &, const SxQuantumNumbers &);



/** The SxQuantumNumbers class is a container class for quantum numbers.
    @ingroup DFT
    @author Sixten Boeck
    */
class SxQuantumNumbers
{
   public:
      int i;
      int iSpin;
      int ik;
      int iSpecies, iAtom, n, l, m;

      enum AngularMomentum { none=-1, s=0, p=1, d=2, f=3 };

      static char lToChar (int l)  {
         SX_CHECK (l >= 0, l);
         switch (l)  {
            case 0  : return 's'; break;
            case 1  : return 'p'; break;
            case 2  : return 'd'; break;
            case 3  : return 'f'; break;
            case 4  : return 'g'; break;
            case 5  : return 'h'; break;
            case 6  : return 'i'; break;
            case 7  : return 'j'; break;
            case 8  : return 'k'; break;
            default : SX_EXIT;
         }
         return '?';
      }

      SxQuantumNumbers (int iIn=0, int iSpinIn=0, int ikIn=0)  {
         SX_CHECK (iIn >= 0, iIn);
         SX_CHECK (iSpinIn >= 0 && iSpinIn < 2, iSpinIn);
         SX_CHECK (ikIn >= 0, ikIn);

         i     = iIn;
         iSpin = iSpinIn;
         ik    = ikIn;

         iSpecies = iAtom = n = l = m = -1;
      }

      SxQuantumNumbers (int isIn, int nIn, int lIn, int mIn)  {
         SX_CHECK  (isIn >= 0, isIn);
         SX_CHECK  (nIn  >= 0, nIn);
         SX_CHECK  (lIn  >= 0 && lIn < 9, lIn);
         SX_CHECK (mIn  >= -lIn && mIn <= lIn, mIn, lIn);
         iSpecies = isIn;
         n        = nIn;
         l        = lIn;
         m        = mIn;
         i = iSpin = ik = iAtom = -1;
      }

      SxQuantumNumbers (int isIn, int iaIn, int nIn, int lIn, int mIn)  {
         SX_CHECK  (isIn >= 0, isIn);
         SX_CHECK  (iaIn >= 0, iaIn);
         SX_CHECK  (nIn  >= 0, nIn);
         SX_CHECK  (lIn  >= 0 && lIn < 9, lIn);
         SX_CHECK (mIn  >= -lIn && mIn <= lIn, mIn, lIn);
         iSpecies = isIn;
         iAtom    = iaIn;
         n        = nIn;
         l        = lIn;
         m        = mIn;
         i = ik = iSpin = -1;
      }

      SxQuantumNumbers (int isIn, int iaIn, int nIn, int lIn, int mIn,
                        int iSpinIn, int ikIn)  {
         SX_CHECK  (iSpinIn >= -1 && iSpinIn < 2,  isIn);
         SX_CHECK  (ikIn >= -1, ikIn);
         SX_CHECK  (isIn >= -1, isIn);
         SX_CHECK  (iaIn >= -1, iaIn);
         SX_CHECK  (nIn  >= -1, nIn);
         SX_CHECK  (lIn  >= -1 && lIn < 9, lIn);
         SX_CHECK (mIn  >= -lIn && mIn <= lIn, mIn, lIn);
         iSpecies = isIn;
         iAtom    = iaIn;
         n        = nIn;
         l        = lIn;
         m        = mIn;
         iSpin    = iSpinIn;
         ik       = ikIn;
      }

      SxQuantumNumbers (const SxQuantumNumbers &in)  {
         i        = in.i;
         iSpin    = in.iSpin;
         ik       = in.ik;
         iSpecies = in.iSpecies;
         iAtom    = in.iAtom;
         n        = in.n;
         l        = in.l;
         m        = in.m;
      }

      bool operator== (const SxQuantumNumbers &o) const
      {
         return (   iSpecies == o.iSpecies
                 && iAtom    == o.iAtom
                 && n        == o.n
                 && l        == o.l
                 && m        == o.m
                 && i        == o.i
                 && iSpin    == o.iSpin
                 && ik       == o.ik);
      }

};


inline std::ostream &operator<< (std::ostream &s, const SxQuantumNumbers &in)
{
   s << "QN:";
   if (in.iSpecies != -1)  s << " is="    << in.iSpecies;
   if (in.iAtom    != -1)  s << " ia="    << in.iAtom;
   if (in.n        != -1)  s << " n="     << in.n;
   if (in.l        != -1)  s << " l="     << in.l;
   if (in.m        != -1)  s << " m="     << in.m;
   if (in.i        != -1)  s << " i="     << in.i;
   if (in.iSpin    != -1)  s << " iSpin=" << in.iSpin;
   if (in.ik       != -1)  s << " ik="    << in.ik;
   s << " ";
   return s;
}

#endif /* _SX_Q_NUMBERS_H_ */
