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
// AUTHOR: Lars Ismer, ismer@fhi-berlin.mpg.de

#ifndef _SX_ARTIFACT_FILTER_H_
#define _SX_ARTIFACT_FILTER_H_

#include <SxOperator.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxSpeciesData.h>
#include <SxGeom.h>


/** \brief Apply constraints to forces

    \b SxArtifactFilter = S/PHI/nX Artifact Filter

    This filter is used to project out artifacts of conifgurational vectors
    (like forces, velocities etc.) or matrices (hessian, dynamical).
    In the current implementation it can (only) 
    eliminate artificial translational and rotational parts, i.e. performs 
    a transformation to internal coordinates.
    
    \author Lars Ismer, ismer@fhi-berlin.mpg.de */

class SX_EXPORT_GEOM SxArtifactFilter 
   : public SxOperatorBase<SxVector<TPrecTauR> > SXOP_LINKFIX
{
   public:
 
      SxArtifactFilter ();
      virtual ~SxArtifactFilter ();

      SXOPERATOR_GETCOPY(SxArtifactFilter, SxVector<TPrecTauR>);

      /** \brief artifact filter is setted up with an input structure 
                 and an input-freezeMode (determines which degrees of 
                 freedom shall be considered as external)*/
      void set (const SxAtomicStructure &, const SxString &, bool);
      /** \brief changes/sets up freezeMode*/
      void setFreezeMode (const SxString &);
      /** \brief changes/sets up structure*/
      void setStructure (const SxAtomicStructure &);
      /** \brief determines wether projections should be done in 
                 mass-weighted coordinates */
      void setMassWeighting (bool);
      /** \brief masses are setted up by a nDoF-dimensional vector*/
      void setMasses (const SxVector<Double> &);
      /** \brief masses are setted up by SxSpeciesData*/
      void setMasses (const SxSpeciesData &);
      /** \brief not under usage yet*/
      void validate (const SxVector<Int> &equivalentIdx);
      /** \brief optimizes D-Matrix with respect to maximum overlap
        of eigenvectors of manipulated matrix with eigenvectors of 
        unmanipulated matrux
           */
      void optimizeDMatrix (const SxMatrix<Double> &);

   protected:
      /** \brief D-Matrix performs projection to internal coordinates 
           Reference: Vibrational Analysis in Gaussian, J.W. Ochterskii,
           http://www.lct.jussieu.fr/manuals/Programmes/Gaussian98/vib/vib.pdf
           */
      SxMatrix<Double> getDMatrix ();
      /** \brief translates the structure to the center of mass, returnvalue 
                 is the 3-d coordinate of the center of mass*/
      SxVector3<Double> translateToCOM (SxVector<Double> *);
      /** \brief gets the inertia tensor (tensor of 2nd order) of the system*/ 
      SxMatrix<Double> getInertiaTensor (const SxVector<Double> &);
      /** \brief gets the matrix*/
      SxVector<Double> getMasses ();

      /** \brief performs the projection/filtering on a given vector*/
      virtual SxVector<TPrecTauR> operator*(const SxVector<TPrecTauR> &) const; 
      virtual void applyInPlace (SxVector<TPrecTauR> &) const { SX_EXIT; }
      
      SxString freezeMode;
      SxAtomicStructure tau;
      bool massWeighting;
      bool needsUpdateDM;
      SxVector<Double> massVec;
      SxMatrix<Double> DMatrix;
};

#endif /* _SX_ARTIFACT_FILTER_H_ */
