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

#ifndef _SX_TAYLOR_EXP_POTENTIAL_H_
#define _SX_TAYLOR_EXP_POTENTIAL_H_

#include <SxClassic.h>
#include <SxPotential.h>
#include <SxArtifactFilter.h>

/** \brief Potential based on Taylor Expansion around equilibrium structure 

    \b contains hessian and shall contain also higher order derivatives

    \author Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_CLASSIC SxTaylorExpPotential : public SxPotential
{
   public:

      SxTaylorExpPotential ();
      SxTaylorExpPotential (const SxSymbolTable *);
      SxTaylorExpPotential (const SxAtomicStructure &, const  SxSymbolTable *);

      virtual ~SxTaylorExpPotential ();
      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);
      
      virtual SxAtomicStructure getForces (const SxAtomicStructure &,
                                           const SxSymbolTable * = NULL);

      virtual SxSpeciesData getSpeciesData () const;
      virtual PrecEnergy getEnergy () const;
 
      SxAtomicStructure loadStructureFHI98 (const SxString &);

      SxArray<SxAtomicStructure> getNumericalHesseMatrix 
         (const SxAtomicStructure &, const double &);
      SxArray<SxArray<SxAtomicStructure> > getNumericalTTensor (const SxAtomicStructure &,
            const double &);      
      void exportHesse (const SxArray<SxAtomicStructure> &, const SxString &);
      void exportTTensor (const SxArray<SxArray<SxAtomicStructure> > &, const SxString &);      
      double e(const SxVector<Int> &);
      int factorial(int);
      
      void symmetrizeHessian ();
      bool containsContribution (const SxString &);
      SxMatrix<Double> hessian;
      SxArray<SxArray<SxVector<Double> > >  ttensor;
      int contributions;
      int nAtoms, nDoF;
      SxAtomicStructure equTau;
      SxSpeciesData speciesData;
      double totalEnergy;
      SxSymbolTable *hessianGroup;
};

#endif /* _SX_TAYLOR_EXP_POTENTIAL */
