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

#ifndef _SX_HAM_INIT_H_
#define _SX_HAM_INIT_H_

#include <SxDFT.h>
#include <SxSpeciesData.h>
#include <SxPWSet.h>
#include <SxAtomicStructure.h>
#include <SxHamiltonian.h>

/** \brief Generic initial guess class

    \author C. Freysoldt, freysoldt@mpie.de */
class SxHamInit
{
   protected:
      /// Warn about missing implementation
      void notImplemented (const char* func) const;
      /// Default implementation of setupWaves for SxPW
      static
      SxPtr<SxPWSet> setupPW (const SxSymbolTable *wavesGroup,
                              const SxPtr<SxGkBasis> &gkPtr,
                              int nStates, int nSpin);
   public:
      // --- virtual interface for initialGuess

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
                                         int nStates, int nSpin) = 0;
      /// Randomize waves
      virtual void randomize (const SxPtr<SxPWSet> &wavesPtr);
      /// Set up the Hamiltonian
      virtual SxPtr<SxHamiltonian> setupHam (const SxPtr<SxPWSet> &,
                                          const SxPtr<SxSpeciesData> &,
                                          const SxAtomicStructure &) = 0;

      /// Virtual destructor
      virtual ~SxHamInit () {}


      // --- static setup functions -------------------------

      /** \brief Set a potential from input file
        @param  table   symbol table
        @return         pointer to the potential

        This routine allows to set up the potential without
        specifying in the code what potential will be used.
        Instead, the symbol table (i.e. the input file) will
        be scanned for known potential types.

        The main use of this routine is to set up a species data
        as seen by the main code when part of the initialization
        is potential-dependent. For instance, the valenceCharge for
        PAW potentials may be read from the potential file.

        */
      static
      SxPtr<SxSpeciesData> setupPotential (const SxSymbolTable *table);
      /** \brief Determine number of components in G/R basis

         Some Hamiltonians work with multi-component wavefunctions,
         e.g. SxKdotP, SxStrain, etc.

         This static routine defines the number of components from
         the input file's symbol table.
      */
      static int getNComp (const SxSymbolTable *top);
};

#endif /* _SX_HAM_INIT_H_ */
