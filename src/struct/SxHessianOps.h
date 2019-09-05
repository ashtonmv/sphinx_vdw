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

#ifndef _SX_HESSIANOPS_H_
#define _SX_HESSIANOPS_H_

#include <SxMatrix.h>
#include <SxVector.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxElemDB.h>
#include <SxStruct.h>


/** \brief sxhessian ops 

    \b SxHessianOps = a tool to deal with hessian/dynamical matrices

    \ingroup   group_structure
    \author    Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_STRUCT SxHessianOps
{
   public:
      /**\brief dynamical matrix (hessian in mass-weighted coordinates) 
                MUST BE of Complex type in order to be applicable 
                when calculating e.g. phonon dispersion relations*/
      SxMatrix<Complex16> dynamical;
      /**\brief hessian matrix; MUST BE Complex*/
      SxMatrix<Complex16> hessian;
      /**\brief masses*/
      SxVector<Double> massVec;
      /**\brief eigensystems of dynamical and hessian*/
      SxMatrix<Complex16>::Eigensystem eigHessian, eigDyn;
      /**\brief number of degrees of freedom in the system*/
      int nDoF;

      //---------------------------------------------------------------------
      /**\brief Constructors and Destructors */
      //---------------------------------------------------------------------
      SxHessianOps ();
      SxHessianOps (const SxHessianOps &);
      SxHessianOps (const SxMatrix<Complex16> &, const SxVector<Double> &);
      ~SxHessianOps ();

      /** = operator copies the content*/
      SxHessianOps &operator= (const SxHessianOps &in);
      
      void set (const SxMatrix<Complex16> &, const SxVector<Double> &);
      void set (const SxHessianOps &);
      /**\brief sets the hessian from refinment calculations*/
      void setFromRefinement (const SxHessianOps &, const SxMatrix<Double> &);
      
      /**\brief   gets eigenvelocities (eigenvectors of d. m.)*/
      SxMatrix<Double> getEigenVelocities ();
      /**\brief   gets eigenfrequencies in 1/cm (sqrt of eigenvalues d. m.)*/
      SxVector<Complex16> getEigenFrequencies ();
      /**\brief   get reduced masses */
      SxVector<Double> getReducedMasses ();
      /**\brief  the curvatures along eigenmodes can be changed */
      void setCurvaturesDynBasis (SxList<SxList<double> > &);
      /**\brief  the frequencies can be changed */
      void setFrequencies (SxList<SxList<SxComplex16> > &);
      /**\brief  gets a suitable displacement vector for refinement 
                 calculations*/
      SxVector<Double> getDisplacementVector (int );
      /**\brief  gets the hessian matrix*/
      SxMatrix<Complex16> getHessian ();
      /**\brief  gets the dynamical matrix*/
      SxMatrix<Complex16> getDynamical ();
          /**\brief  returns the number of degrees of freedom*/
      int getNDoF () const;
      /**\brief writes hessian matrix and eigenfrequencies*/
      void write (const SxString &);
      /**\brief prints molden-format file (contains eigenmodes)*/
      void printMolden (const SxString &, const SxAtomicStructure &, int);
      /**\brief returns symmetrized matrix for an non-symmetric input 
                matrix (utility routine)*/
      SxMatrix<Double> getSymmetrizedMatrix 
         (const SxMatrix<Double> &);
      /**\brief returns a matrix with orthonormalized columns 
                (utility routine)*/
      SxMatrix<Complex16> getOrthonormalizedMatrix 
         (const SxMatrix<Complex16> &);

   protected: 
      /**\brief is switched to true, whenever one of the set-procedures 
                is called; is checked whenever output-routines are 
                called; if false, eigenvalues -modes and -frequencies 
                are recalculated */
      bool needsUpdate; 
      /**\brief   performs necessary diagonalization of dynamical matrix 
                  and hessian*/
      void update ();
      /**\brief   sets the hessian matrix */
      void setHessian (const SxMatrix<Complex16> &);
      /**\brief   sets the masses */
      void setMasses (const SxVector<Double> &);

};


#endif /* _SX_TIME_CORR_H_ */
