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

#ifndef _SX_ORBITALS_H_
#define _SX_ORBITALS_H_

class SxOrbital
{
   public:
      SxOrbital (int is_=0, int ia_=0, int n_=0, int l_=0, int m_=0)  
      {
         iSpecies = is_; iAtom = ia_; n = n_;  l = l_; m  = m_;
      }
      int iSpecies, iAtom, n, l, m;
      void print () {
         printf ("SxOrbital: is=%d, ia=%d, n=%d, l=%d, m=%d\n",
                 iSpecies, iAtom, n, l, m);
      }
      bool operator== (const SxOrbital &o) const
      {
         return (   iSpecies == o.iSpecies
                 && iAtom    == o.iAtom
                 && n        == o.n
                 && l        == o.l
                 && m        == o.m);
      }
      bool operator< (const SxOrbital &o) const
      {
         return (   iSpecies < o.iSpecies
                 && iAtom    < o.iAtom
                 && n        < o.n
                 && l        < o.l
                 && m        < o.m);
      }
      bool operator> (const SxOrbital &o) const
      {
         return ( !(*this < o) );
      }
};

#endif /* _SX_SX_ORBITALS_H_ */
