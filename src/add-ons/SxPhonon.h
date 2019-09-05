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

#ifndef _SX_PHONON_H_
#define _SX_PHONON_H_

#include <SxString.h>
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxMath.h> 
#include <SxArray.h>
#include <SxMurn.h>
#include <SxAtomicStructure.h>
#include <SxDynMat.h>
#include <SxKPoints.h>
#include <SxExt.h>


/** \brief Compute the phonon frequencies (meV) for the provided q-points from
 * the dynamical matrix in the eigenvector basis. 
 
 Frequency refinement is also supported. This uses the eigenvectors from one set
 of calculations and determines the corresponding eigenvalues from a single
 calculation in which all phonons have a certain amplitude. This could be used
 for efficient calculation of the volume dependence of phonon frequencies.
 TODO add some equations

    \b SxPhonon = S/PHI/nX Compute Phonon frequencies

    \author Matth\'e Uijttewaal, uijttewaal@mpie.de 
 */
class SX_EXPORT_EXT SxPhonon 
{
   protected:
      // --phase factors (move to SxMath?)
      inline SxComplex<double> expI (const double &phi) const 
      { return cos (phi) + I * sin (phi); }
      
      inline SxVector<Complex16> expI (const SxVector<Double> &phi) const 
      { return cos (phi) + I * sin (phi); }
      
      // some coordinate functions (move to SxMath?)
      inline Coord cosC (const Coord &c) const
      { return Coord (cos (c(0)), cos (c(1)), cos (c(2))); }

      inline Coord absSqrt (const Coord &coord) const { 
         Coord absC = coord.absSqr (); 
         return Coord (sqrt (absC(0)), sqrt (absC(1)), sqrt (absC(2))); 
      }

      // --take a sqrt of a matrix, retaining its signs (move to SxMath?)
      SxMatrix<Double> sgnSqrt (const SxMatrix<Double> &x) const { 
         ssize_t nElt = x.getSize (), nCol = x.nCols ();

         //determine the correct signs
         SxMatrix<Double> sgnMat (nElt / nCol, nCol);
         for (ssize_t iElt = 0; iElt < nElt; iElt++)  { 
            if (x(iElt) > 1.e-11)  {sgnMat(iElt) = 1.;}
            else if (x(iElt) < -1.e-11)  {sgnMat(iElt) = -1.;}
            else  {sgnMat(iElt) = 0.;}
         }

         //calculate the result
         SxMatrix<Double> res (sqrt (x.abs ()) * sgnMat);
         res.reshape (nElt / nCol, nCol);
         return res;
      }

      // -squares of a number (move to SxMath?)
      inline double sqr (const double &x) const {return x*x;}
      inline int sqr (const int &x) const {return x*x;}

      // cube of a number (move to SxMath?)
      inline double cube (const double &x) const {return x*x*x;}

   public:
      // (almost) empty constructor
      SxPhonon ()  { scal = 1.; }
      
      // destructor
     ~SxPhonon ()  { /*empty*/ }

     /**\brief read the input data

       @param setInput file with the q-point settings
       @param dynMatInput file with the dynamical matrices
       @param epsE the epsEqual value of the structure
       @return the eigenvalues of the dynamical matrix
       */
      SxArray<SxMatrix<Complex16>::Eigensystem> readInput (
            const SxString &setInput,
            const SxString &dynMatInput,
            const double &epsE);

      /**\brief get the lattice shifts of the primitive cell
         
        uses DynMat::getExactQ; last three are the superCell lattice vectors
       \note is there a direct relation exactQ <-> lattice shifts?
       @param output the output file
       @return the lattice shifts + supCell vectors
       */
     SxArray<Coord> getLatShifts (ofstream &output) const; 

   protected:
      int sQ, sMod; //special qPoint & mode index
      SxKPoints qPoints; //:nQ; q-points of qPath/qMesh
      SxArray<Coord> exactQ; //:nCell
      SxAtomicStructure primStr; //:3*nAtom; primitive structure
      double scal; //scaling factor for phonon refinement

     /**\brief Fourier reduction factor at supercell borders

       General 'reduce' function: in principle, near the borders and for
       inexact qPoints the fourier factors need to be corrected for image
       contributions. However, this can only be done if there is an inversion
       for the displaced atom, but this also means that there must be an atom
       at the other side of the border and so that any correction will be
       small (atoms can't be too close together). The only correction that
       remains is for atoms precisely at the border.
       TODO: simplify?
      @param iQ qPoints index
      @param supCell the supercell (for the borders)
      @return fourier reduction matrices for superstructure
      */
      SxMatrix<Complex16> fourierReduce (
            const int iQ, 
            const Coord &latShift,
            const SxCell &supCell,
            const bool debug) const;
      
   public:
      /**\brief fourier transform the dynamical matrix from the exact qPoints
       * to the eigensystem of the provided set of qPoints

        \note the fourier phases are reduced at the borders
       @param inEigs input eigensystems
       @param latShifts the lattice shifts + the supercell
       @param debug are we in debug mode?
       @return the transformed dynMats
       */
      SxArray<SxMatrix<Complex16>::Eigensystem> transformDynMat (
            const SxArray<SxMatrix<Complex16>::Eigensystem> &inEigs, 
            const SxArray<Coord> &latShifts,
            const bool debug);
      
     /**\brief print frequencies (meV) of dynMat and frequency and eigenvector
      *  of q-point sQ and mode sMod. If no q-point and/or mode is provided
      *  the minima are printed. Count also the number of negative frequencies. 
      
      @param eigs the resulting eigensytems
      @param sxbOut the binary output file
      @param output the ASCII output file
      */
      void printFreq (const SxArray<SxMatrix<Complex16>::Eigensystem> &eigs,
                      const SxString &sxbOut,
                            ofstream &output) const;

      void printFull (const SxArray<SxMatrix<Complex16>::Eigensystem> &eigs) const;


/** \brief refine already calculated frequencies using known eigenvectors and
  one calculation with a displacement pattern inlcuding all eigenvectors */
 
      /**\brief read refineInput file (displacements, forces)
        
        The drift is corrected and the species are splitted. The primStr and
        qPoints are not reset, only a scaling factor is saved.
        @param input input file
        @param latShifts lattice shifts to check the superstructure (+supercell)
        @return atomic structures of displacements and forces
        */
     SxArray<SxAtomicStructure> readRefine (
           const SxString &input,
           const SxArray<Coord> &latShifts);
 
     /**\brief decompose force and displacement structures in phonon modes

       \note experimental feature
       @param strs atomic structures of displacements and forces
       @param eigVecs eigenvectors
       @param latShifts lattice shifts, last three are supCell
       @return relative phonon components of forces and displacements: values
       */
     SxArray<SxVector<Double> > decompose (
           const SxArray<SxAtomicStructure> &strs, 
           const SxArray<SxMatrix<Complex16> > &eigVecs,
           const SxArray<Coord> &latShifts) const;
};
  
#endif /* _SX_PHONON_H_ */
