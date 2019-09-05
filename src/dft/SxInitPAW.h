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

#ifndef _SX_INIT_PAW_H_
#define _SX_INIT_PAW_H_

#include <SxDFT.h>
#include <SxHamInit.h>
#include <SxPAWHamiltonian.h>
class SxHamSolver;

/** \brief Initializer class for PAW


    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxInitPAW : public SxHamInit
{
   public:
      /// Pointer to species (i.e. potential) data
      SxPtr<SxSpeciesData> potPtr;
      /// PAW Hamiltonian
      SxPtr<SxPAWHamiltonian> pawHam;
   public:
      SxInitPAW (const SxSymbolTable *table,
                 const SxPtr<SxSpeciesData> &potPtrIn,
                 const SxPtr<SxGkBasis> &gkPtr,
                 const SxHamSolver *hamSolver);
      /// Randomize density
      virtual void randomRho ();
      /// Atomic density
      virtual void atomicRho (const SxSymbolTable *rhoGroup,
                              const SxAtomicStructure &structure);
      /// Density from file
      virtual void readRho (const SxString &file);
      /// Density from wave functions
      virtual void rhoFromWaves (const SxPtr<SxPWSet> &wavesPtr,
                                 const Focc &focc);
      /// Add extra charge to density
      virtual void addCharge (const SxMeshR &chargeR);
      /// Initialize waves
      virtual SxPtr<SxPWSet> setupWaves (const SxSymbolTable *wavesGroup,
                                         const SxPtr<SxGkBasis> &gkPtr,
                                         int nStates, int nSpin);
      /// Randomize waves
      virtual void randomize (const SxPtr<SxPWSet> &wavesPtr);
      /// Set up the Hamiltonian
      virtual SxPtr<SxHamiltonian> setupHam (const SxPtr<SxPWSet> &,
                                          const SxPtr<SxSpeciesData> &,
                                          const SxAtomicStructure &);
      /// Virtual destructor
      virtual ~SxInitPAW () {}

};

#endif /* _SX_INIT_PAW_H_ */

