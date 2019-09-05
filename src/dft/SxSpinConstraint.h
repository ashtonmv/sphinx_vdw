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

#ifndef _SX_SPINCONSTRAINT_H_
#define _SX_SPINCONSTRAINT_H_

#include <SxSpinConstraint.h>
#include <SxFermi.h>
#include <SxPAWSet.h>
#include <SxPAWRho.h>
#include <SxRho.h>
#include<SxHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxSpectrum.h>
#include <SxPartialWaveBasis.h>
#include <SxAtomicStructure.h>
#include <SxProjector.h>
#include <SxRBasis.h>
#include <SxGrid.h>


class SX_EXPORT_DFT SxSpinConstraint
{
   public:

      SxPtr<SxPAWPot> pawPotPtr;
      SxVector<Double> nuA;
      SxVector<Double> targetSpin;
      /// whether or not a specific atom is spin-constrained (:iTlAtom)
      SxArray<bool> constrainedAtom;
      SxSpinConstraint( const SxPtr<SxPAWPot> &pawPotPtrIn);
      /// Read constraints from symbol table
      void read (const SxSymbolTable *topLvl,
                 const SxAtomicStructure &structure);
      
      void setNuA (double value);
      /** \brief Set the constraints
          @param targetSpinIn - the desired target spins
          @param str - structure for symmetry checking
        */
      void setConstraints ( const SxVector<Double> &targetSpinIn,
                           const SxAtomicStructure &str);
      /// Check compatibility of spin constraint with symmetry
      void checkSym (const SxAtomicStructure &str) const;
      double calcKappaOpt (const SxVector<Double> &MSpinIn, const SxVector<Double> &MSpinPlusIn, double kappa);
      /// Compute new nu Lagrangian parameters for given waves
      SxVector<Double> computeNu (const SxPtr<SxPAWSet> &wavesPtr, SxFermi &fermi, const Real8 ekt, double epsilon = 1e-20);

   private:
      /// Get partial volume in state-space for all atoms
      SxArray<SxMatrix<Complex16> > getOmegaNN (int iSpin, int ik) const;
      /// Get spins for current nuA with subspace diagonalization
      SxVector<Double> getSpinsDiag ( const Eps &eps0In, SxFermi &fermiIn, const double &ektIn , SxPtr<SxPAWSet> wavesPtr = SxPtr<SxPAWSet>() );
      /// Pointer to the currently used waves
      const SxPAWSet *wvPtr;

};
#endif  /*_SX_SPINCONSTRAINT_H_*/
