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

#ifndef _SX_PAW_OVERLAP_H_
#define _SX_PAW_OVERLAP_H_

#include <SxDFT.h>
#include <SxOverlap.h>
#include <SxPartialWaveBasis.h>

/** \brief Projector augmented wave overlap operator

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWOverlap : public SxOverlapBase
{
   protected:
      /// PW projectors
      SxConstPtr<SxPartialWaveBasis> pBasis;
      /// PAW potential
      SxConstPtr<SxPAWPot>       pawPot;
   public:
      /** \brief Constructor
        */
      SxPAWOverlap (const SxConstPtr<SxPartialWaveBasis> &pIn,
                    const SxConstPtr<SxPAWPot>           &pot)
         : pBasis(pIn), pawPot(pot)
      { /* empty */ }

      /// Norm square
      virtual double normSqr (const SxDiracVec<Complex16> &psi) const;

      /// Scalar product
      virtual SxComplex16 dot (const SxDiracVec<Complex16> &x,
                               const SxDiracVec<Complex16> &y) const;

   protected:
      /// Partial wave contribution to scalar product
      SxComplex16 dotPartial (const SxDiracVec<Complex16> &px,
                              const SxDiracVec<Complex16> &py) const;

   public:
      /// Apply overlap operator
      virtual SxDiracVec<Complex16>
      apply (const SxDiracVec<Complex16> &psi) const;

      /// Set states x orthogonal to states y
      virtual void setOrthogonal (SxDiracVec<Complex16> *xPtr,
                                  const SxDiracVec<Complex16> &y) const;

      /// Get PAW correction to overlap matrix from projections <p|psi>
      SxDiracMat<TPrecCoeffG>
      getDeltaS (const SxDiracVec<TPrecCoeffG> &px, 
                 const SxDiracVec<TPrecCoeffG> &py) const;

      /// Get PAW overlap matrix elements
      virtual SxDiracMat<Complex16>
      getMatrix (const SxDiracVec<Complex16> &x,
                 const SxDiracVec<Complex16> &y) const;

      /// Orthonormalize states x
      virtual void orthonormalize (SxDiracVec<Complex16> *xPtr,
                                   const SxOrthoMethod how=GramSchmidt) const;
      /// Normalize states x
      virtual void normalize (SxDiracVec<Complex16> *xPtr) const;

};

#endif /* _SX_PAW_OVERLAP_H_ */

