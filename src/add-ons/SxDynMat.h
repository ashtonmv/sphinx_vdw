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

#ifndef _SX_DYN_MAT_H_
#define _SX_DYN_MAT_H_

#include <SxString.h>
#include <SxAtomicStructure.h>
#include <SxMath.h> 
#include <SxConstants.h>
#include <SxVector.h>
#include <SxMatrix.h> 
#include <SxSymMatrix.h> 
#include <SxGrid.h>
#include <SxKPoints.h>
#include <SxExt.h>

/** \brief Compute dynamical matrix from a set of displacements/forces calculations. 
    
  The single input file contains the information on the species and all sets
  (structures) of coordinates and forces (and possibly spins). The first set
  needs to be of the undisplaced structure. Arbitrary structures and
  displacements are allowed as long as these provide sufficient data to
  construct the dynamical matrix. All symmetries are determined and applied.
  The dynamical matrix is written to a S/PHI/nX binary file in the Eigenbasis
  format (exact q-points and their Eigenvalues and Eigenvectors) together with
  the primitive structure.
  TODO add some equations

    \b SxDynMat = S/PHI/nX Compute Dynamical Matrix

    \author Matth\'e Uijttewaal, uijttewaal@mpie.de 
 */
class SX_EXPORT_EXT SxDynMat 
{
    public:
      // constructor
      SxDynMat () { /*empty*/ }

      // destructor
     ~SxDynMat ()  { /*empty*/ }

    protected:
      // vectorial phase factor (move to SxMath?)
      inline SxComplex<double> expI (const double &phi) const
       { return cos (phi) + I * sin (phi); }
      inline SxVector<Complex16> expI (const SxVector<Double> &phi) const {
         return cos (phi) + I * sin (phi);
      }
   
    public:
      /*undisplaced structure; Note: the no. species (nSpec) could differ from 
        the no. atoms in the primitive cell (nPrimAt).*/
      SxAtomicStructure str; 

   protected:
      //reciprocal masses for all nPrimAt (used in hermit)
      SxVector<Double> recMass;
      
      /*all displacement structures; the size increases from (larger than) the 
        no. degrees of freedom (nDof) to (larger than) the number of modes
        (nMode) of the system*/
      SxArray<SxAtomicStructure> dTauStr;

      /*Fourier transformed (complex!) displacements in matrix form. It has dimensions: 
        no. (exact) q-points, nQ * (3 * nPrimAt, nMode)*/
      SxArray<SxMatrix<Complex16> > dTau; 
      
      //all force structures corresponding to dTauStr
      SxArray<SxAtomicStructure> forcStr;

      SxAtomicStructure backgroundForces;

      //Fourier transformed (complex) forces corresponding to dTau
      SxArray<SxMatrix<Complex16> > forces; 

          /**\brief edit the input data: split the species, correct for
           * drift & initial forces, perform mass scaling, determine the
           * displacements and set the primitive cell with the spin allowed symmetries

           $F_i = F_i - F_0$
           $F_i = F_i - sum_a F_ia / nA$
           $F_i = F_i * sqrt (M)$
           $dTau_i = str_i - str_0$
           $dTau_i = dTau_i / sqrt(M)$ Add more indices?
           @param initForc the forces for the initial structure
           @param initSpin the spins of the initial structure
           @recMass the reciprocal masses
           @return the primitive cell with the symmetries of the structure
        */
      SxCell editInput (const SxAtomicStructure &initForc, 
                        const SxAtomicStructure &initSpin,
                        const SxVector<Double> &recMass);

   public:
      /**\brief read file with coords and forces ( and possibly spins) for all 
        structures. 
  
         Forces for the initial structure are optional.
         The displacement structures, the reciprocal masses, the primitive cell
         and the symmetry operations allowed by the spins are set. 
         The forces are corrected for drift and scaled with the masses.
         @param input the input file
         @param output the output file
         @param epsE the accuracy of the structure
         @return the primitive cell with the symmetries
       */
      SxCell readInput (const SxString &input, 
                              ofstream &output, 
                        const double &epsE);
   
      /**\brief get the exact q-points of the supercell as the multiples of
       * reciprocal cell vectors in the reciprocal of the primitive cell

        $xQ = n1 * R1 + n2 * R2 + n3 * R3$, where the sum of the vectors is
        within the reciprocal primitive cell 
        @param primCell the primitive cell
        @param output the output file
        @return the array of exact q-points
        */
      SxArray<Coord> getExactQ (const SxCell &primCell, ofstream &output) const;

   protected:      
      /**\brief symmetrise the forces of mode iDof if iSym is a symmetry of the
       * corresponding displacement, internal routine

         $S ^ str_Pi <reorder> sgn R ^ dTau_Pi$ if equal -> 
         $F_i/2 + sgn R ^ F_Pi/2$ &
         $dTau_i/2 + sgn R ^ dTau_Pi/2$ 
         @param iSym the index of the symmetry operation to apply
         @param iDof the mode number to apply it to
         @return symmetrised?
       */
      bool symMod (const int iSym, const int iDof); 

      /**\brief generate new force (+displacement) modes, internal routine
          N.B. multiple equivalent dTau could be constructed!

          $R ^ dTau_Pi, R ^ F_Pi$ for all displ. changing symm!
         @param genSyms symmetry operations to use
         @return total no. modes
       */  
      int genNewMods (const SxArray<SxArray<SxSymOp> > &genSyms, //:nDof,nSym
                      const bool debug);

   public:
      /**\brief apply all symmetry operations (within the primitive cell) to
       * symmetrise the forces and to generate new, symmetry mapped forces and displacements

         ${S} / kernel dTau_i$
         @param output file to print extra information
         @return no. new modes
       */
      int applySymOps (ofstream &output, const bool debug);

   protected:
      /**\brief transform dTauStr and forcStr to q-space (they become matrices)

        $dTauQ = sum_cells dTau_i * exp(I Q ^ str)$
        $FQ = sum_cells F_i * exp(I Q ^ str)$
        @param exactQ the exact q-points of the supercell
        */
      void fourier (const SxArray<Coord> &exactQ);
  
      /* make the dynamical matrix hermitian (real eigenvalues) without breaking the
       * translation symmetry (3 frequencies at gamma are zero) 
        TODO use SxArtifactFilter instead?

        $dynMat/2 + dynMat^+/2$
        $dynMat - (sum dynMat / sqrt (M))*3/nMode*sqrt (M)$, repeat
       @param dynMatQ the dynamical matrix to be corrected
       @param driftCor correct drift as well?
       @param debug are we in debug mode?
       @return dynMatQ hermitian?
       */
      bool hermit (SxMatrix<Complex16> &dynMatQ, 
                   bool driftCor, 
                   const bool debug);

   public:
      /**\brief get the dynamical matrix for the exact q-vectors from a
       * transformation of the forces with dTau
  
         $dynMatQ = - forcesQ ^ dTauQ^-1 * scaling$ -> eigensystem
         @param exactQ q-points for which to determine dynamical matrix
         @param output file to print some information
         @return array of complex Hermitian dynamical matrices
       */   
      SxArray<SxMatrix<Complex16>::Eigensystem> getDynMat (
            const SxArray<Coord> &exactQ, 
            const bool debug);

      /**\brief print the dynamical matrices to a S/PHI/nX binary file
        
         This includes the primitive structure (Bohr), the q-points (1/Bohr),
         the eigenvals (meV^2) and eigenvecs (sqrt(u)*Bohr).
         @param eig the eigensystems of the dynamical matrices
         @param qPoints the q-points to which the dynMats correspond 
         @param primCell the primitive cell
       */
      void print (const SxArray<SxMatrix<Complex16>::Eigensystem> &eig,
                  const SxArray<Coord> &qPoints,
                  const SxCell &primCell) const;

      void printHesse (const SxArray<SxMatrix<Complex16>::Eigensystem> &eig,
                  const SxArray<Coord> &qPoints, const SxCell &primCell, const bool hesse,
                  const bool sxhesse,  const bool bgForces, const double disp) const;

};
#endif /* _SX_DYN_MAT_H_ */
