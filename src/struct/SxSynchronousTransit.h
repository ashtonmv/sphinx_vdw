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

#ifndef _SX_SYNCHRONOUS_TRANSIT_H_
#define _SX_SYNCHRONOUS_TRANSIT_H_

#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxArray.h>
#include <SxPotential.h>
#include <SxSpeciesData.h>
#include <SxHamSolver.h>
#include <SxHessianOps.h>
#include <SxTypes.h>
#include <SxStruct.h>

/** \brief Synchronous Transit 

  Samples the PES along trajectories corresponding to 
  synchronous (linear) displacements of the atoms
  Is helpful for transition state searches and
  refinement calculations of vibrational frequencies 

    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_STRUCT SxSynchronousTransit
{
   public:

      //---------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //---------------------------------------------------------------------
      //@{
      SxSynchronousTransit (const SxAtomicStructure &, SxPotential *);
      SxSynchronousTransit ();
      ~SxSynchronousTransit ();
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
      /**\brief sets the electronic minimizer*/
      void setElMinimCmds (const SxArray<const SxSymbolTable *> &in);

      /**\brief performs a refinement calculations for a dynamical matrix 
                by sampling the PES along the eigenmodes and replacing 
                the eigenvalues by the obtained curvatures
                refineUpTo = upper bound for frequencies to correct
                dX = displacement length
       */
      SxHessianOps getRefinedDynamical 
         (const SxHessianOps &in, double refineUpTo, double dX);
      /**\brief returns the potential energies for a given input array 
                of structures 
       */
      SxArray<double> getEPots  (const SxArray<SxAtomicStructure> &traj);
   
	protected:
      /**\brief species data*/
      SxSpeciesData speciesData;
      /**\brief initial structure*/
      SxAtomicStructure inputTau;
      /**\brief potential to be analysed*/
     	SxPotential *potential;  // TODO: should be const
      /**\brief number of degrees of freedom*/
      int nDoF;
      /**\brief masses stored in a vector*/
      SxVector<Double> massVec;
		SxArray<const SxSymbolTable *> elMinimCmds;

		double getEPot ();
      /**\brief is getting the forces acting on the given structure*/
      SxAtomicStructure getForces (const SxAtomicStructure &);
};


#endif /* _SX_SYNCHRONOUS_TRANSIT_H_ */
