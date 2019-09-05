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

#ifndef _SX_DRIFT_FILTER_H_
#define _SX_DRIFT_FILTER_H_

#include <SxOperator.h>
#include <SxAtomicStructure.h>
#include <SxStruct.h>


/** \brief Project out center-of-mass drift

    \b SxDriftFilter = S/PHI/nX center-of-mass drift filter 

    This class provides a drift filter for forces or displacements
    that should have a vanishing sum over all atoms.

    The filter subtracts the average drift.

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_STRUCT SxDriftFilter : public SxOperatorBase<SxAtomicStructure>
{
   public:
      /// Empty constructor
      SxDriftFilter () {/* empty */}

      /// Empty destructor
      virtual ~SxDriftFilter () { /* empty */ }

      /// Apply the center-of-mass-filter out-of-place
      virtual SxAtomicStructure 
      operator*(const SxAtomicStructure &forces) const
      {
         SxVector3<Double> comForce = forces.sum () / forces.getNAtoms ();
         //cout << "center-of-mass drift (projected out): " << comForce << endl;
         return forces - comForce;
      }

      /// Apply the center-of-mass-filter in-place
      virtual void applyInPlace (SxAtomicStructure &forces) const
      {
         SxVector3<Double> comForce = forces.sum () / forces.getNAtoms ();
         //cout << "center-of-mass drift (projected out): " << comForce << endl;
         forces -= comForce;

      }

      /** \brief Get SxPtr copy
          \note Standard copy constructor is OK, no dangerous members
      */
      SXOPERATOR_GETCOPY(SxDriftFilter,SxAtomicStructure);

};

#endif /* _SX_DRIFT_FILTER_H */
