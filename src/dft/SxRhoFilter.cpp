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

#include <SxRhoFilter.h>
#include <SxPAWHamiltonian.h>
#include <SxProjector.h>

SxDensity SxRhoFilter::operator| (const SxDensity &in)
{
   if (mode == NoFilter) return in.getCopy ();
   if (in.checkType<SxRho> ()) return in.getCopy ();

   // --- PAW rho filtering
   SX_CHECK (mode == PseudoRho);
   SX_CHECK (in.checkType<SxPAWRho> ());
   const SxPAWRho &rho = in.getRef<SxPAWRho> ();
   SxPtr<SxRho> res = SxPtr<SxRho>::create ();
   int nSpin = rho.getNSpin ();
   res->rhoR.resize (nSpin);
   SX_CHECK (rho.pwRho.rBasisPtr);
   res->rBasisPtr = rho.pwRho.rBasisPtr;

   if (!gBasisPtr)  {
      // --- set up filtering G basis
      const SxAtomicStructure* structPtr
         = res->rBasisPtr->getGBasis ().structPtr;
      SX_CHECK (structPtr);
      gBasisPtr = SxPtr<SxGBasis>::create (res->rBasisPtr->fft3d.mesh,
                                           *structPtr, gCut);
   } 
   const SxGBasis &gBasis = *gBasisPtr;
   
   // --- prepare aux data for compensation charges
   int lMaxRho = rho.potPtr->lMaxRho.maxval ();
   SxDiracMat<Double> YlmGl(gBasis.ng, sqr(lMaxRho + 1));
   YlmGl.setBasis (gBasis);
   for (int lm = 0, l = 0; l <= lMaxRho; ++l) {
      SxDiracVec<Double> gl = pow (gBasis.g2, 0.5 * l);
      for (int m = -l; m <= l; ++m,++lm)
         YlmGl.colRef(lm) <<= gBasis.getYlm (l,m) * gl;
   }
   double rNorm = FOUR_PI / sqrt(res->rBasisPtr->cell.volume);

   const SxPAWPot &pawpot = *rho.potPtr;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // Fourier filter the pseudo density
      PsiG rhoG = gBasis | rho.pwRho.rhoR(iSpin);

      // --- now add the compensation charges
      for(int is = 0; is < rho.Dij.getNSpecies (); is++)  {
         //double rc2 = sqr(SxPAWHamiltonian::rSoft);
         double rc2 = sqr(pawpot.rc(is));
         SxDiracVec<Double> gauss = exp ( -(0.25 * rc2) * gBasis.g2);
         for(int ia = 0; ia < rho.Dij.getNAtoms (is); ia++)  {
            double dFac = 1.;
            PsiG T = gBasis.getPhaseFactors (is,ia);
            for (int l = 0; l <= pawpot.lMaxRho(is); ++l, dFac *= 2 * l + 1) {
               // i^l factor
               SxComplex16 il = 1.;
               if      ((l & 3) == 1) il = I; 
               else if ((l & 3) == 2) il = -1.;
               else if ((l & 3) == 3) il = -I;
               for (int m = -l; m <= l; ++m)  {
                  int lm = SxYlm::combineLm(l,m);
                  double Qlm = 0.;
                  SX_LOOP2(ipt,jpt)  {
                     int l1 = pawpot.lPhi(is)(ipt), 
                         off1 = pawpot.offset(is)(ipt) + l1;
                     int l2 = pawpot.lPhi(is)(jpt),
                         off2 = pawpot.offset(is)(jpt) + l2;
                     // now the angular momentum sums
                     for (int m1 = -l1; m1 <= l1; ++m1)  {
                        int lm1 = SxYlm::combineLm(l1,m1);
                        int ipl = off1 + m1;
                        for (int m2 = -l2; m2 <= l2; ++m2)  {
                           int lm2 = SxYlm::combineLm(l2,m2);
                           int jpl = off2 + m2;
                           double cg = pawpot.clebschGordan(lm1, lm2, lm);
                           // no i^l factors for PAW radial ylm
                           // but is l1+l2-l (from cg definition)
                           if ((l1 + l2 - l) & 2) cg = -cg;
                           Qlm += rho.Dij(iSpin,is,ia)(ipl,jpl)
                                * cg
                                * pawpot.QijL(is)(ipt,jpt)(l);
                        }
                     }
                  }
                  double N = rNorm * SxYlm::getYlmNormFactor(l,m) / dFac;
                  double weight = N * Qlm;
                  PsiG grl = (YlmGl.colRef(lm) * gauss) * T;
                  rhoG.plus_assign_ax(il * weight,grl);
               }
            }
         }
      }
      // --- get in R space (TODO: leave it in G)
      res->rhoR(iSpin) = *res->rBasisPtr | rhoG;
   }
   return SxDensity (res);
}

SxRhoFilter::SxRhoFilter (const SxSymbolTable *table)
   : mode (PseudoRho), gCut (10.)
{
   SX_CHECK (table);
   try {
      if (table->contains ("gCut"))
         gCut = table->get ("gCut")->toReal ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
}
