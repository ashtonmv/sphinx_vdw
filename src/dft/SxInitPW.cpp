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

#include <SxInitPW.h>
#include <SxPWHamiltonian.h>

SxInitPW::SxInitPW (const SxRBasis &R, int nSpin, double nElectrons,
                    const SxPtr<SxPseudoPot> &potPtrIn)
   : potPtr(potPtrIn),
     rho(R, nSpin, nElectrons)
{
   // empty
}

void SxInitPW::randomRho ()
{
   rho.randomize ();
}

void SxInitPW::atomicRho (const SxSymbolTable *rhoGroup,
                          const SxAtomicStructure &structure)
{
   SxPseudoPot &psPot = dynamic_cast<SxPseudoPot&> (*potPtr);
   int nSpin = rho.getNSpin ();
   
   int iSpecies, nSpecies = psPot.nSpecies;
   int iSpin;
   int l, lMax;
   double spinMoment, s = 1., sumFocc;
   SxVector<TReal8>  foccSpin(nSpecies);
   // --- create local LCAO occupations
   SxArray<SxArray<SxVector<TReal8> > >  foccAtomicRho(nSpecies); 
   // --- resize local LCAO occupations
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lMax = psPot.lMax(iSpecies);
      foccAtomicRho(iSpecies).resize(lMax+1);
      for (l=0; l <= lMax; l++)  {
         foccAtomicRho(iSpecies)(l).resize(nSpin);
      }
      foccSpin(iSpecies) = (TReal8::Type)
                           psPot.foccAtomicRho(iSpecies)(l=0).getSize();
   }
   if (rhoGroup && rhoGroup->contains("spinMoment"))  {
      double nElectrons = rho.nElectrons;
      spinMoment = rhoGroup->get("spinMoment")->toReal();
      s          = (nElectrons - spinMoment)/(nElectrons + spinMoment);
   }
   
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lMax = psPot.lMax(iSpecies);
      if (fabs(s-1.) > 1e-12 && (nSpin == 2))  {
         for (l=0; l <= lMax; l++)  {
            sumFocc = psPot.foccAtomicRho(iSpecies)(l).sum();
            foccAtomicRho(iSpecies)(l)( SPIN_UP ) = sumFocc * 1./(1.+s);
            foccAtomicRho(iSpecies)(l)(SPIN_DOWN) = sumFocc *  s/(1.+s);
         }      
      }  else  {
         if (foccSpin(iSpecies) == nSpin)  {
            foccAtomicRho(iSpecies) = psPot.foccAtomicRho(iSpecies);
         }  else  {
            for (l=0; l <= lMax; l++)  {
               sumFocc = psPot.foccAtomicRho(iSpecies)(l).sum();
               for (iSpin=0; iSpin < nSpin; iSpin++)  {
                  foccAtomicRho(iSpecies)(l)(iSpin) = 
                     sumFocc / (Real8)nSpin;
               }
            }
         }
      }
   }
   
   // --- copy local occupations to PP object 
   psPot.foccAtomicRho = foccAtomicRho;
   
   // --- calculate atomic charge density
   rho.atomicChargeDensity (structure, psPot);
   
}

void SxInitPW::readRho (const SxString &file)
{
   rho.readRho (file);
}

void SxInitPW::addCharge (const SxMeshR &chargeR)
{
   int nSpin = rho.getNSpin ();
   double totalExtra = chargeR.sum () * rho.rBasisPtr->dOmega;
   rho.nElectrons += totalExtra;
   rho.normalizeRho ();
   rho.nElectrons -= totalExtra;
   for (int iSpin = 0.; iSpin < nSpin; ++iSpin)
      rho(iSpin).plus_assign_ax(-1./nSpin, chargeR);
}


SxPtr<SxPWSet> SxInitPW::setupWaves (const SxSymbolTable *wavesGroup,
                                     const SxPtr<SxGkBasis> &gkPtr,
                                     int nStates, int nSpin)
{
   return setupPW (wavesGroup, gkPtr, nStates, nSpin);
}

void SxInitPW::randomize (const SxPtr<SxPWSet> &wavesPtr)
{
   SxPtr<SxPW> (wavesPtr)->randomize ();
}

void SxInitPW::rhoFromWaves (const SxPtr<SxPWSet> &wavesPtr, const Focc &focc)
{
   rho.computeRho (focc, *SxPtr<SxPW> (wavesPtr));
}

SxPtr<SxHamiltonian>
SxInitPW::setupHam (const SxPtr<SxPWSet> &wavesPtr,
                    const SxPtr<SxSpeciesData> &,
                    const SxAtomicStructure &structure)
{
   SX_CHECK (dynamic_cast<SxPseudoPot*> (potPtr.getPtr ()));
   SxPseudoPot &psPot = dynamic_cast<SxPseudoPot&> (*potPtr);
   SxPW &waves = *SxPtr<SxPW>(wavesPtr);
   return SxPtr<SxPWHamiltonian>::create (waves, rho.rhoR, psPot, structure);
}
