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

#include <SxInitKP.h>
#include <SxPW.h>
#include <SxKdotP8x8.h>
#include <SxKdotP.h>
#include <SxStrain.h>

SxInitKP::~SxInitKP () {}

SxPtr<SxHamiltonian>
SxInitKP::setupHam (const SxPtr<SxPWSet> &wavesPtr,
                    const SxPtr<SxSpeciesData> &,
                    const SxAtomicStructure &)
{
   SxPW &waves = *SxPtr<SxPW>(wavesPtr);
   RhoR rhoR(1);
   rhoR(0) = SxMeshR (rPtr);
   rhoR(0).set (0.);
   if (type == KdotP8x8)  {
      cout << "8-band k.p Hamiltonian used here...";
      return SxPtr<SxKdotP8x8>::create (waves, rhoR);
   } else if (type == KdotPgeneral)  {
      cout << "n-band k.p Hamiltonian used here...\n";
      return SxPtr<SxKdotP>::create (waves, rhoR);
   } else if (type == StrainField)  {
      cout << "Strain field calculation...";
      return SxPtr<SxStrain>::create (waves, rPtr);
   } else {
      SX_EXIT;
   }
   return SxPtr<SxHamiltonian> ();
}

void SxInitKP::randomize (const SxPtr<SxPWSet> &wavesPtr)
{
   SxPtr<SxPW> (wavesPtr)->randomize ();
}
