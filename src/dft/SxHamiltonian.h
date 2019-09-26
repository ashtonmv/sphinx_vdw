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

#ifndef _SX_HAMILTONIAN_H_
#define _SX_HAMILTONIAN_H_

#include <SxPrecision.h>
#include <SxFermi.h>
#include <SxSymbolTable.h>
#include <SxDFT.h>

#include <SxPsiSet.h>
#include <SxPWSet.h>
#include <SxAtomicStructure.h>
#include <SxDensity.h>
#include <SxPWOverlap.h>
#include <SxXC.h>
#include <SxVDW.h>

/** \brief Abstract Hamiltonian \f$\hat{H}\f$ class

  \sa     \ref page_dirac
  \author Sixten Boeck
  */
class SX_EXPORT_DFT SxHamiltonian
{
   public:
      SxHamiltonian () {  }
      virtual ~SxHamiltonian ()     {  }

      virtual void compute (const SxFermi &,
                            bool /*tauChanged*/=true, 
                            bool /*rhoChanged*/=true)  
      {
        SX_EXIT; 
      };

      /** \brief Compute new energy*/
      virtual PrecEnergy getEnergy (const SxPsiSet &, const SxFermi &)
      {
         SX_EXIT; return 0.;
      }

      /** Get last computed energy */
      virtual PrecEnergy getEnergy () const
      {
         SX_EXIT; return 0.;
      }
       
      /** Get last computed double counting correction */
      virtual PrecEnergy getDoubleCounting () const
      {
         SX_EXIT; return 0.;
      }
      
      virtual void printEnergies () const = 0;
      virtual void writePotentials () const { }
      virtual void read (const SxSymbolTable *) = 0;
      /** \brief Update internal settings from some input group
          @param table The input group to read from
          @return true if anything changed

          @note: Must be overloaded by derived Hamiltonian if there
                 is any changeable feature. Must also overload
                 backToDefault then.
       */
      virtual bool rereadTable (const SxSymbolTable *table);
      /** \brief Set internal settings to default group
          @param table some group in the symbol tree
       */
      virtual void backToDefault (const SxSymbolTable *table);

      virtual void setXCMeshDensity (int) { }
      /** \brief Compute new density from waves and occupations */
      virtual void computeRho (const Focc &, const SxPsiSet &) { }//= 0;
      /** \brief Preconditioner
        @param psi wave function to be preconditioned
        */
      virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &) const
      { SX_EXIT; return SxDiracVec<TPrecCoeffG::TReal> (); }
      /** \brief Get access to density */
      virtual SxDensity& getRho () { SX_EXIT; return *(SxDensity *)NULL; };
      
      // VDW Correction
      bool applyVDWCorrection;
      SxVDW VDW;

      static const SxSymbolTable *getHamiltonianGroup (const SxSymbolTable *);
      static int getNSpin (const SxSymbolTable *);
      static int getNEmptyStates (const SxSymbolTable *);
      static int getNEmptyStatesChi (const SxSymbolTable *);
      /** \brief Compute the number of states in the calculation
          \todo: find better formula for nSpin = 2
        */
      static int getNStates (const SxSymbolTable *, double nElectrons, 
                             double spinMoment = 0.);
      static double getNExcessElectrons (const SxSymbolTable *);
      /** \brief Compute the number of electrons in the calculation
          @param table symbol table
          @param valCharge valence charges for each species
          @param str the atomic structure (for getting number of atoms)
        */
      static double getNElectrons (const SxSymbolTable *table,
                                   const SxArray<double> &valCharge,
                                   const SxAtomicStructure &str);
      static double getEkt (const SxSymbolTable *);

      static bool withDipoleCorrection (const SxSymbolTable *);
      /** \brief Current structure */
      SxAtomicStructure structure;
      virtual PrecEnergy getETrial()   {
         return 0.;
      }
      
      /*
      enum WhatToCompute {
         RhoParts       = 0x0001,
         WaveParts      = 0x0002,
         StructureParts = 0x0004,
         ComputeEnergy  = 0x0100,
         ComputeForces  = 0x0200
      };
      */

      /// Provide overlap operator
      virtual SxConstPtr<SxOverlapBase> getS () const;

      virtual PsiG apply (const PsiG &) const
      {
         SX_EXIT; return PsiG ();
      }
      
      PsiG operator* (const PsiG &psi) const
      {
         return apply (psi);
      }

      PsiG operator| (const PsiG &psi) const
      {
         return apply (psi);
      }

};

#endif /* _SX_HAMILTONIAN_H_ */
