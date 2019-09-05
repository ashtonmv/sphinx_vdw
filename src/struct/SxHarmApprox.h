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
#ifndef _SX_HARMAPPROX_H_
#define _SX_HARMAPPROX_H_

#include <SxMatrix.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxPotential.h>
#include <SxHamSolver.h>
#include <SxStruct.h>

/** \brief Calculates phonons in the harmonic approximation

    \b SxHarmApprox = SPHInX Harmonic Approximation

    The dynamical matrix, which is the crucial item in the harmonic 
    approximation, can be computed in the direct space using symmetries 
    of the according structure. Also the dynamical matrix in momentum 
    space can be calculated via different methods. Futhermore the eigenvalues
    and eigenvectors of that dynMatK can be calculated.

    Names convention:


    \author Blazej Grabowski, blazejgrabowski@gmx.de */

class SX_EXPORT_STRUCT SxHarmApprox
{
   public:
      SxHarmApprox ();
      SxHarmApprox (bool forceNoSym);
      ~SxHarmApprox ();
      
      /** brief Forces the cutoff to have a fixed value given by the input. 
         The cutoff indicates the sphere around the moved atom in which
         the forces are taken into account. */
//      void forceRCut(double rCut);

      /** brief Calculates the cutoff which indicates the sphere around 
        the moved atom in which the forces are taken into account. */
//      double getRCut(const SxAtomicStructure &str, const SxCell &cell);

      /** brief Checks if the difference vector between the two input vectors
        lies in the rCut sphere. Takes periodic images also into account! 
       It returns true in the fourth element of the output vector if this
       is the case. The first three elements include the relative vector.
       Otherwise it returns false.*/
//      SxVector<Double> insideSphere(const SxVector3<Double> &tau1, const 
//            SxVector3<Double> &tau2);

      /** brief Compares two vectors */
      bool equal(const Coord &a, const Coord &b);

      /** brief Returns the index of the vector if it is included in the
        box, otherwise it returns -1. 

        /param a vector which is checked 
        /param box */
      int getBoxIndex(const Coord &a, const SxArray<Coord> &box);
  
      /** brief Computes the lattice vector in units of the small cell
        belonging to vector a 

        /param a 
        /param cell */
      Coord getRLat(const Coord &a, const SxCell &cell);

      /** brief Structure check

        Checks if the input structure can be recreated by moving
        along the small cell vetors */
      void checkStr(const SxAtomicStructure &str, const SxCell &smallCell);

      /** brief Sets additional properties of the structure 
        
        */
      void setIndexesAndMoved();

      
      /** brief Center both structures around (0,0,0)

        For calculation of the dynamical matrix it is usefull to
        have the small structure which specifies the moved atoms
        in the center of the big structure. */
      void centerStr();
      
      /** brief Sets class variables
      
      Moves also the structure so that the first atom lies on 0,0,0. */
      void set(/*SxPotential *potential_, */
            const SxSymbolTable *symbolTable, const SxCell &smallCell_, 
            double deviation_);

      /** brief Prints structure, rLat-, tauIndex, moved, rLat-, tauBox */
      void printStr();

      /** brief Prints the dynamical matix */
      void printDynMatR(bool detail);
      
      /** brief Saves dynamical matrix */
      void saveDynMatR(const SxString &outFile);

      /** brief Loads dynamical matrix from file */
      void loadDynMatR(const SxString &dynMatRFile);
      
      /** brief Saves the big structure */
      void saveStr(const SxString &file, const SxAtomicStructure &s);
      
      SxAtomicStructure getForcesDFT(const SxAtomicStructure &str,
      	const SxSymbolTable *table);
         
      /** brief Computes the dynamical matrix 
        
        The dynamical matrix is computed in direct space only for the
        displacment in one direction indicated by input parameter alpha
        but for al atoms of all species in the small structure.        

        /param alpha displacment direction (0 = x, 1 = y, 2 = z) */     
      void computeDynMatR1d(int alpha, const SxSymbolTable *symbolTable);

      /** brief Computes dynamical matrixsymmetry operations

        Checks if the deviation in y direction or z, respective, can be 
        constructed from the x deviation. For this purpose the isolated
        structure (vacuum arround the supercell) has to have symmetry
        operations that convert x into y axis and also for x - z. It 
        returns the operations. */
      SxArray<SxArray<SymMat> > getDynMatRSymOp();

      /** brief Fills dynamical matrix through symmetries */
      void fillDynMatR(int beta, const SymMat &S);
      
      void cut();
      
      /** brief Computes whole dynamical matrix */
      void computeDynMatR(const SxSymbolTable *symbolTable);

      void symmetrizeDynMatR();
      
      void extendDynMatR();
      
      double getRCut();

      void cutSphereDynMatR();
      
      /** No Computation, just filling an existing dynamical matrix */
      bool setLoadFillSave(const SxSymbolTable *symbolTable, 
         const SxCell &smallCell_, const SxString &inputDynMatR,
         const SxString &outputDynMatR);

      /** brief Computes eigen frequencies */
      SxArray<SxVector<Complex16> > computeEigFreq(const SxArray<Coord> &kPoints);


   protected:
 
       /** brief Returns a row of the sub dynamical matrix */
       Coord getSubDynMatRrow(int n, int r, int c, int beta);

       /** brief Sets a row of the sub dynamical matrix*/
       void setSubDynMatRrow(int n, int r, int c, int beta, const Coord row);
      
       // SxPotential *potential; TODO 
       SxSpeciesData     species;
       double deviation, rCut;
       bool forceRCut, forceNoSym;
       SxArray<SxMatrix<Double> > dynMatR, dynMatK;
       SxAtomicStructure str;
       SxArray<SxVector<Int> > rLatIndex, tauIndex, moved;
       SxCell smallCell;
       SxArray<Coord> rLatBox, tauBox;
       int nRLat, nTau;
};

#endif /* _SX_HARMAPPROX_H_ */
