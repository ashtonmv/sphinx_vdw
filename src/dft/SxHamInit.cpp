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

#include <SxHamInit.h>
#include <SxError.h>
#include <typeinfo>
#include <SxPW.h>
#include <SxPAWPot.h>

void SxHamInit::notImplemented (const char* func) const
{
   cout << func << " not implemented (in "
        << typeid(*this).name () << ")" << endl;
}

void SxHamInit::randomRho ()
{
   notImplemented (SX_FUNC);
}


void SxHamInit::atomicRho (const SxSymbolTable *,
                                const SxAtomicStructure &)
{
   notImplemented (SX_FUNC);
}

void SxHamInit::readRho (const SxString &)
{
   notImplemented (SX_FUNC);
}


void SxHamInit::rhoFromWaves (const SxPtr<SxPWSet> &,
                                 const Focc &)
{
   notImplemented (SX_FUNC);
}

void SxHamInit::addCharge (const SxMeshR &)
{
   notImplemented (SX_FUNC);
}


void SxHamInit::randomize (const SxPtr<SxPWSet> &)
{
   notImplemented (SX_FUNC);
}

SxPtr<SxPWSet> SxHamInit::setupPW (const SxSymbolTable *wavesGroup,
                                        const SxPtr<SxGkBasis> &gkPtr,
                                        int nStates, int nSpin)
{
   SxString tmpDir = "";
   if (wavesGroup
       && wavesGroup->contains("keepWavesOnDisk")
       && wavesGroup->get("keepWavesOnDisk")->toAttribute ())
   {
      tmpDir = ".";
   }

   return SxPtr<SxPW>::create (nStates, nSpin, gkPtr, tmpDir);
}

SxPtr<SxSpeciesData> SxHamInit::setupPotential (const SxSymbolTable *table)
{
   const SxSymbolTable *top = table->topLevel();
   // --- normconserving pp
   if (   top->containsGroup ("PWHamiltonian")
       || top->containsGroup ("Hamiltonian")
       || top->containsGroup ("pseudoPot"))
   {
     return SxPtr<SxPseudoPot>::create (top);
   }
   // --- PAW
   if (top->containsGroup ("PAWHamiltonian") || top->containsGroup ("pawPot")) {
      return SxPtr<SxPAWPot>::create (top);
   }
   // --- generic
   return SxPtr<SxSpeciesData>::create (top);
}


int SxHamInit::getNComp (const SxSymbolTable *top)
{
   if (top->containsGroup("KP8x8Hamiltonian"))  {
      cout << "initialize basis for 8-Band k.p" << endl;
      return 8;
   }
   if (top->containsGroup("kpHamiltonian"))  {
      const SxSymbolTable *ham = top->getGroup("kpHamiltonian");
      int nBands = ham->get("nBands")->toInt();
      cout << "initialize basis for " << nBands << "-Band k.p" << endl;
      return nBands;
   }
   if (top->containsGroup("StrainField"))  {
      cout << "initialize basis for strain field calculation" << endl;
      return 3;
   }
   return 1;
}
