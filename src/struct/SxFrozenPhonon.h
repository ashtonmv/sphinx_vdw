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

#ifndef _SX_FROZEN_PHONON_H_
#define _SX_FROZEN_PHONON_H_

#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxArray.h>
#include <SxPotential.h>
#include <SxSpeciesData.h>
#include <SxArtifactFilter.h>
#include <SxHamSolver.h>
#include <SxHessianOps.h>
#include <SxTypes.h>
#include <SxSecondaryStructure.h>
#include <SxThermodynamics.h>
#include <SxSynchronousTransit.h>
#include <SxStruct.h>

/** \brief Frozen Phonon 

    \b SxFrozenPhonon = S/PHI/nX Frozen Phonon
    
    Performs evaluation of Hessian matrix for a given structure 
    on basis of a finite-differences scheme. 
    Desired features:
      - optional outout of phonon dispersion relation 
      - optional output of vibrational thermodynamic properties
      - full adaption of symmetry-properties to reduce numerical effort 
    ... works for helices yet (03/17/05)

    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_STRUCT SxFrozenPhonon
{
   public:

      //---------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //---------------------------------------------------------------------
      //@{
      SxFrozenPhonon (const SxAtomicStructure &, SxPotential *);
      SxFrozenPhonon ();
      ~SxFrozenPhonon ();
      //@}
      //---------------------------------------------------------------------
      /**@name Interface to input file
         Controlling the class by the \ref tutor_parser. */
      //---------------------------------------------------------------------
      //@{
      void print (const SxSymbolTable *);
		void execute (const SxSymbolTable *, bool calc=true);
		//@}
      static bool isRegistered (const SxSymbolTable *);
   
	protected:
      /**\brief species data*/
      SxSpeciesData speciesData;
      /**\brief initial structure*/
      SxAtomicStructure inputTau;
      /**\brief forces acting on initial structure*/
      SxAtomicStructure inputForces;
      /**\brief potential to be analysed*/
     	SxPotential *potential;  // TODO: should be const
      /**\brief number of degrees of freedom*/
      int nDoF;
      /**\brief masses stored in a vector*/
      SxVector<Double> massVec;
      /**\brief hessian matrix serves as basis 
                for refinement calculations */
      SxHessianOps basisHOps;
      /**\brief electronic loops*/
		SxArray<const SxSymbolTable *> elMinimCmds;
      
      /**\brief is the system a helix ?*/
      bool secondaryStructure;
      /**\brief object to deal with the structure and dynamics of the helix*/
      SxSecondaryStructure peptideChain;
      /**\brief contains the phonon basis in case of refinement calculations*/
      SxSecondaryStructure basisPeptideChain;

		double getEPot ();
      /**\brief is getting the forces acting on the given structure*/
      SxAtomicStructure getForces (const SxAtomicStructure &);
      /**\brief gets the force response for a given displacement direction*/
		SxAtomicStructure getResponse (const SxAtomicStructure &, double, bool);
      /**\brief gets the i'th displacement in a given basis*/
      SxAtomicStructure getDisplacement (int , const SxString &);
      /**\brief composes the hessian in cartesian coordinates 
        for a given response matrix in a given basis*/
      SxMatrix<Double> getHessian (const SxMatrix<Double> &, const SxString &);
};


#endif /* _SX_FROZEN_PHONON_H_ */
