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

#ifndef _SX_PW_OVERLAP_H_
#define _SX_PW_OVERLAP_H_

#include <SxOverlap.h>
#include <SxDFT.h>

/** \brief Plane-wave overlap operator interface

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPWOverlap : public SxOverlapBase
{
   public:
      /// Norm square
      virtual double normSqr (const SxDiracVec<Complex16> &psi) const;

      /// Scalar product
      virtual SxComplex16 dot (const SxDiracVec<Complex16> &x,
                               const SxDiracVec<Complex16> &y) const;

      /// Apply overlap operator
      virtual SxDiracVec<Complex16> 
      apply (const SxDiracVec<Complex16> &psi) const;
      
      /// Set states x orthogonal to states y
      virtual void setOrthogonal (SxDiracVec<Complex16> *xPtr,
                                  const SxDiracVec<Complex16> &y) const;

      /// Orthonormalize states x
      virtual void orthonormalize (SxDiracVec<Complex16> *xPtr,
                                   const SxOrthoMethod how=GramSchmidt) const;
             
};

#endif /* _SX_PW_OVERLAP_H_ */
