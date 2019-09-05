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

#include <SxPartialWaveBasis.h>
#include <SxPAWPot.h>

SxPartialWaveBasis::SxPartialWaveBasis (const SxConstPtr<SxPAWPot> &pot,
                                        const SxAtomicStructure &str)
   : nElements(0)
{
   SX_CHECK (pot);
   potPtr = pot;

   SX_CHECK (potPtr->getNSpecies () == str.getNSpecies (),
             potPtr->getNSpecies (), str.getNSpecies ());
               
   for (int is = 0; is < str.getNSpecies (); ++is)
      nElements += potPtr->getNProj (is) * str.getNAtoms (is);

   SX_CHECK (nElements > 0);
}

void SxPartialWaveBasis::print () const
{
   cout << SX_SEPARATOR;
   cout << "| Partial wave basis |p>" << endl;
   cout << "| Number of partials:      " << nElements << endl;
   cout << "| Projections to Gk-basis: " 
        << (projectors ? "enabled" : "disabled") << endl;
   cout << SX_SEPARATOR;
}

SxDiracVec<TGBasisType>
SxPartialWaveBasis::toPWBasis (const SxGBasis *G,
                               const SxDiracVec<Complex16> &pCoeff) const
{
   SX_CHECK (projectors);
   SxDiracVec<TGBasisType> res = projectors->toPWBasis (G, pCoeff);
   SX_CHECK (res.getBasisPtr () == G);
   return res;
}

SxPtr<SxAOBasis> 
SxPartialWaveBasis::createProjBasis (const SxGkBasis &gk,
                                     const SxPAWPot &pawPot)
{
   SxArray<SxDiracMat<Double> > projRad(pawPot.getNSpecies ());
   for (int is = 0; is < pawPot.getNSpecies (); ++is)  {
      if (    pawPot.pPsFine.getSize () > 0
           && pawPot.pPsFine(is).getSize () > 0)
      {
         projRad(is) = pawPot.pPsFine(is);
      } else {
         projRad(is) = pawPot.pPS(is);
      }
   }
   return SxPtr<SxAOBasis>::create (gk, projRad, pawPot.lPhi);
}


