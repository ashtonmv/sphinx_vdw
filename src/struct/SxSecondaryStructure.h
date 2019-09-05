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

#ifndef _SX_LINEAR_CHAIN_H_
#define _SX_LINEAR_CHAIN_H_

#include <SxVector.h>
#include <SxAtomicStructure.h>
#include <SxArtifactFilter.h>
#include <SxSpeciesData.h>
#include <SxStruct.h>

/**
  @brief A Class to deal with structure and vibrations of linear chains,
         helices etc.
  */
class SX_EXPORT_STRUCT SxSecondaryStructure 
{
      
   public:
      SxSecondaryStructure ();
      SxSecondaryStructure (const SxString &, const SxString &,  
                     const SxString &, int, bool);
      
      SxSecondaryStructure (const SxAtomicStructure &, const SxSpeciesData &,
                     const SxString &,  
                     const SxString &, int, bool);
      
      ~SxSecondaryStructure ();
      /**\brief sets the hessian matrix (input hessian should 
                contain column-wise the 
                force response for a reference peptide, 
                but not necessarily must be 
                a full hessian*/
      void setHessian (const SxMatrix<Double> &, int, bool);
      /**\brief sets the hessian matrix in case 
                it is output of a refinement calc. */
      void setHessianFromRefinement 
         (const SxMatrix<Double> &, const SxMatrix<Double> &, int, bool);
      /**\brief gets symmetry center*/
      SxArray<double> getSymCenter ();
      /**\brief gets (symmetrized) Coords*/
      SxAtomicStructure getCoords ();

      /**\brief gets full hessian in cartesian coordinates*/
      SxMatrix<Double> getFullCartHessian ();
      /**\brief gets full hessian in cylindrical coordinates*/
      SxMatrix<Double> getFullCylHessian ();
      /**\brief get the pitch of the helix/FES */
      double getZPitch ();
      /**\brief gets the number of peptides in the helix/FES */
      int getNPeptides ();
      /**\brief gets the number of atoms per petide */
      int getNAtomsPP ();
      /**\brief gets the number of degrees of freedom per supercell*/
      int getNDoF ();
      /**\brief gets the number of degrees of freedom per irreducible unit*/
      int getNDoFPIU ();
      /**\brief gets frequency for given phase angle (point in BZ)*/
      SxComplex16 getFrequency (int, double, const SxString &);
      /**\brief gets a spline interpolation for the lowest frequency 
                branch of helices by employing symmetry constraints 
                at twist angle and 180 degs and taking the spline points 
                at numerically stable points from the fourier interpolation */
      SxArray<double> splineItLowestBranch (const SxArray<double> &);
      /**\brief resizes the kMesh and computes the phonon-dispersion relation 
                on this mesh, if kMesh size is larger then number of peptides,
                the kMesh contains also (fourier)-interpolated values*/
      void resizeKMesh (int);
      /**\brief prints phonon dispersion to xmgrace readable file*/
      void printPhononDispersionRelation 
         (const SxString&, const SxString&, double);
      /**\brief prints out specific force constants of the helix*/
      void printFC (const SxSymbolTable *);
      /**\brief */
      void rescaleForceConstants (int, double, int);
      /**\brief gets the whole freqenucy spectrum*/
      SxVector<Complex16> getFrequencies 
         (const SxString&, double);
      /**\brief gets the symmetry reduced force constant matrix
                for a given point in the vibrational BZ*/
      SxMatrix<Complex16> getK (double); 
         /**\brief transforms a DoF in peptide-order to a Dof in input-order*/
      int getDoFInputOrder (int);
      /**\brief gets the cartesian block-matrix corresponding to the 
                interaction of a peptide with it's n'th nearest neighbor*/
      SxMatrix<Double> getKcart (int);
      /**\brief gets interaction coefficient for the n'th nearest neigbor
        and a given phase angle  */
      SxComplex16 getIaCoeff (double, const SxVector<Complex16> &, int);
      /**\brief gets expectation value of dynamical matrix for a 
        given input vector
        and a given phase angle  */
      SxComplex16 getEValue (double, const SxVector<Complex16> &, 
                             const SxArray<SxArray<SxMatrix<Double> > > &);
      /**\brief gets first derivative of interaction coefficient for 
        the n'th nearest neigbor and a given phase angle  */
      SxComplex16 getIaCoeffdWdPhi (double, const SxVector<Complex16> &, int, double);

      /**\brief gets first derivative of interaction coefficient for 
        the n'th nearest neigbor and a given phase angle erases matrix 
        elements according to low importance  */
      SxComplex16 getIaCoeffdWdPhiThinned (double, const SxVector<Complex16> &, 
           const SxMatrix<Double> &, int, double);
      /**\brief gets importance measure for force constants  */
      SxMatrix<Double> getImportanceMatrix (double, const SxVector<Complex16> &, int);

      /**\brief gets importance measure for force constants reg. freq. */
      SxMatrix<Double> getImportanceMatrixFreq (int, int);
      /**\brief permutates an vector from input to peptide order*/
      SxVector<Double> getPeptideOrderedVec (const SxVector<Double> &);
      /**\brief permutates an vector from peptide to input order*/
      SxVector<Double> getInputOrderedVec (const SxVector<Double> &);
      /**\brief returns a displacement vector for a refinement calc.*/
      SxVector<Double> getDisplacementVector (int);
      /**\brief returns a normed displacement vector*/
      SxVector<Double> getNormedDisplacementVector (int);
      /**\brief gets a suitable resolution for the phonon dispersion relation*/
      double getPDRes ();
      /**\brief transforms a displacement vector from cartesian 
                to cylindrical coordinates*/
      SxVector<Double> displCartToCyl (const SxVector<Double> &);
      /**\brief transforms a displacement vector from cylindrical
                to cartesian coordinates*/
      SxVector<Double> displCylToCart (const SxVector<Double> &);
      /**\brief prints out the interaction blocks in terms 
                of cartesian coordinates*/
      void printBlocks (const SxString &);
      /**\brief prints a matrix*/
      void printBlock (const SxMatrix<Double> &, const SxString &, bool);
      /**\brief gets the number of atoms in between two atoms along the chain (to define nearest neighbor interactions)*/
      int getInChainDistance (int, int, int, int);
      /**\brief determines whether the interaction in between a given pair of atoms is inside a given cutoff shell or not*/
      bool inShell (int, int, int, int, int);


   protected:
      /**\brief updates cylindrical coordinates 
                attention: cartesian structure is shifted to 
                symmetry center before*/
      void updateCylCoords ();
      /**\brief updates cartesian coordinates*/
      void updateCartCoords ();
      /**\brief updates the B-matrix for transf. in between helix symmetry
        and cartesian coordinates (attention: cylindrical coordinates must be 
       generated first*/
      void updateBMatrix ();
      /**\brief switches to peptide order*/
      void switchToPeptideOrder ();
      /**\brief mapps input to peptide order for structure and masses*/
      void peptideOrderTau ();
      /**\brief switches to input order*/
      void switchToInputOrder ();
      /**\brief mapps peptide to input order for structure and masses*/
      void inputOrderTau ();
      /**\brief loads periodicity information from a file*/
      void setPeriodicity(const SxString &);
      /**\brief loads atomic structure + species information (fhi98md format)*/
      void setStructure (const SxString &);
      /**\brief symmetrizes coords (according to helical symmetry)
                the implemented symmetrization is a bit ugly 
                (but it works*/
      void symmetrizeCoords ();
      /**\brief cuts the interaction beyond the given threshhold value*/
      SxArray<SxArray<SxMatrix<Double> > > 
         applySupercellCutoff 
         (const SxArray<SxArray<SxMatrix<Double> > > &, int);
      /**\brief applies symmetry on blockwise stored hessian*/
      SxArray<SxArray<SxMatrix<Double> > > 
         applySymmetry (const SxArray<SxArray<SxMatrix<Double> > > &, bool);
      /**\brief applies point-group symmetry for beta-sheets*/
      SxMatrix<Double>  
         applyBSheetSymmetry (const SxMatrix<Double> &);
     /**\brief gets periodic indices*/
      SxArray<int> getPeriodicIndices (int); 

      /**\brief sets the chain type (helix, fes or linear*/
      void setChainType (const SxString &);

      /**\brief performs memory allocation for blockwise storage of hessian*/
       SxArray<SxArray<SxMatrix<Double> > > getInitialBlockArray ();
      /**\brief gets a blocked nPeptides*nPeptide peptide wise stored 
        array from an input ordered Hessian  */
      SxArray<SxArray<SxMatrix<Double> > > getPeptideOrderedBlocks 
         (const SxMatrix<Double> &);
       /**\brief gets a blocked nPeptides*nPeptide  array from a Hessian  */
            SxArray<SxArray<SxMatrix<Double> > >  getBlocks 
         (const SxMatrix<Double> &);
      /**\brief writes blockwise storage to matrix*/
      SxMatrix<Double> getHessianFromBlocks 
         (const SxArray<SxArray<SxMatrix<Double> > > &);
      /**\brief writes blockwise peptide-ordered storage to input ordered hessian*/  
        SxMatrix<Double> getHessianFromBlocksInputOrder
        (const SxArray<SxArray<SxMatrix<Double> > > &);

      /**\brief transforms hessian from cylindrical to cartesian coordinates*/
      SxMatrix<Double> getCartHessian 
         (const SxMatrix<Double> &, const SxAtomicStructure &);
      /**\brief transforms hessian from cartesian to cylindrical coordinates*/
      SxMatrix<Double> getCylHessian 
         (const SxMatrix<Double> &, const SxAtomicStructure &);
      /**\brief transforms hessian from cartesian to FES coordinates*/
      SxMatrix<Double> getFESHessian 
         (const SxMatrix<Double> &);
      /**\brief gets partial derivatives, which are needed by the above transformations*/
      double getT(const SxString &, int, const SxAtomicStructure &); 
      /**\brief computes the exactly known points in the BZ*/
      void computeKMesh ();
      /**\brief computes normalDisplacements and curvatures (see below)*/
      void computeNormalDisplacements ();
       

      /**\brief symmetrizes a matrix (should be elsewhere)*/
      SxMatrix<Double>  getSymmetrizedMatrix (const SxMatrix<Double> &); 

      /**\brief utility functions for spline interpolation*/
      void spline 
         (const SxArray<double> &, const SxArray<double> &, int , 
          double, double,  SxArray<double> *);
      
      double splint 
         (const SxArray<double> &, const SxArray<double> &, 
          const SxArray<double> &, int n, double x);

      
      /**\brief contains periodicity information*/
      SxList<SxVector<Double> > periodicity;
      
      
      
      /**\brief cartesian coordinates*/
      SxAtomicStructure tauCart;
      /**\brief cartesian coordinates shifted to center of symmetry in xy-plane*/
      SxAtomicStructure tauCartShifted;
      /**\brief cylindrical coordinates*/
      SxAtomicStructure tauCyl;
      /**\brief B-Matrix (see above: updateBMatrix () )*/
      SxMatrix<Double> B;
      /**\brief symmetry center (x, y)*/
      SxArray<double> symCenter;
      
      /**\brief hessian in cartesian coordinates*/
      SxMatrix<Double> hessianCart;
      /**\brief hessian in cylindrical coordinates*/
      SxMatrix<Double> hessianCyl;

      /**\brief blockwise, peptide ordered storage of the hessian 
                in cartesian coordinates*/
      SxArray<SxArray<SxMatrix<Double> > > blockArrayCart;
      /**\brief blockwise, peptide ordered storage of the hessian 
                in cylindrical coordinates*/
      SxArray<SxArray<SxMatrix<Double> > > blockArrayCyl;

      /**\brief contains a suitable set of displacements to perform 
                refinement calculations of the harmonic frequencies
                within the frozen Phonon approach */
      SxMatrix<Double> normalDisplacements;
      
      /**\brief contains the curvatures, which determine the 
                length of finite displacement*/
      SxVector<Double> normalCurvatures;
      
      /**\brief contains the points of 
                in the phonon dispersion relation; 
                by default the mesh size equals the number of 
                peptides in the supercell, meshsize can be changed 
                by calling the resizeKMesh() function */
      SxArray<SxVector<Complex16> > kMeshFreqs;

      


      /**\brief ionic masses */
      SxAtomicStructure masses;

      /**\brief number of peptides  
                contained in the supercell */
      int nPeptides;
      /**\brief number of irreducible units 
                contained in the supercell */
      int nUnits;
      /**\brief number of kpoints */
      int nKPoints;  
      /**\brief number of atoms per Peptide*/
      int nAtomsPP;
      /**\brief number of atoms per irreducible Unit*/
      int nAtomsPU;
      /**\brief twist angle in rad */
      double twistRad;
      /**\brief twist angle in degs*/
      double twistDeg;
      /**\brief zPitch */
      double zPitch;
      /**\brief number of turns*/
      double nTurns;
      /**\brief total number of degrees of freedom*/
      int nDoF;
      /**\brief degrees of freedom per Peptide*/
      int nDoFPP;   
      /**\brief number of degrees of freedom 
                per irreducible unit  */
      int nDoFPU;
      /**\brief number of rotations of the atomic basis, which 
        are vibrations of 
                the chain (i.e are not free rotations)*/
      int vibRots;

      /**\brief chainType */
      SxString chainType;

      /**\brief flag which indicates, whether data is actually in peptide 
                or in input order */

      bool isInInputOrder;
      


};      
#endif /* _SX_LINEAR_CHAIN_H_ */
