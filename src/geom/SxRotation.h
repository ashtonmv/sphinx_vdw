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

#ifndef _SX_ROTATION_H_
#define _SX_ROTATION_H_

#include <SxGeom.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxSymType.h>
#include <SxString.h>

/// A real space coordinate (or vector)
typedef SxVector3<Double>  Coord;
/// A real space rotation matrix
typedef SxMatrix3<Double>  SymMat;

/** @brief (3x3) matrix with approximate equality operator

    @ingroup group_structure
    @author  C. Freysoldt, freyso@fhi-berlin.mpg.de */
class SX_EXPORT_GEOM SxRotation : public SymMat 
{
   public:
      SxRotation () : SymMat () {}
      SxRotation (const SymMat &in) : SymMat(in) {}
      /** \brief Constructor from rotation axis and rotation angle

        \param axis  The rotation axis (non-zero!)
        \param angle Rotation angle in rad (full rotation \f$2\pi$)
        */
      SxRotation (const Coord &axis, double angle);
      bool operator== (const SymMat &in);
      /** \brief Get symmetry element name of a matrix
        */
      static SxString getName (const SymMat& in);
      /** \brief Get symmetry element name
        */
      SxString getName ()   { return getName (*this); }
      /// Get symmetry type
      static SxSymType getType (const SymMat &in);
      /// Get symmetry type
      SxSymType getType () { return getType (*this); }
      /// Get hash
      static size_t hash (const SxRotation &);
};

#endif /* _SX_ROTATION_H */
