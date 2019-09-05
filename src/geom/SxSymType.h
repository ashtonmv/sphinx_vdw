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

#ifndef _SX_SYM_TYPE_H_
#define _SX_SYM_TYPE_H_

#include <SxVector3.h>
#include <SxGeom.h>
#include <SxString.h>

/** \brief Symmetry type classification

    \b SxSymType = SPHInX Symmetry Type Classification

    This container class is used for classification of a symmetry element.
    It is mainly used as return type of the SxSymFinder class.

    \ingroup  group_structure
    \author   Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_GEOM SxSymType  
{
   public:
      enum Classification { None, Identity, Inversion, Mirror, 
                            Rotation, MirrorRotation };

      Classification     classification;

      /** \brief Rotational symmetry axis count

          A rotational symmetry operator is defined for only certain degrees
          (\f$ \alpha_{S} = 360, 180, 120, 90, 60 \f$ degrees). 
          The axis count is defined as number
          of applicable rotations to reach a full circle
          \f[
             n_{\mathrm{axis count}} = 360 / \alpha_{S}
          \f]
          Hence, only 1, 2, 3, 4, and 6 are possible axis count numbers */
      int                axisCount;  // 0, 1, 2, 3, 4, 6
      /** \brief Coordinates of a symmetry element

          Often an additional vector is required in order to specify a symmetry
          operator uniquely. In case of a mirror plane (SxSymType::Mirror)
          the opCoord is the normal vector to the mirror plane.
          All rotational symmetry operators (SxSymType::Rotation and
          SxSymType::MirrorRotation) define the rotational axis using the
          opCoord variable.
          */
      SxVector3<Double>  opCoord;
      /** \brief Symmetry operator identification string

          The identifier is a human readable string defining the symmetry
          element according to the Schoenfliess nomenclature:
          - Identity: E  (german, [E]inheitsoperator)
          - Inversion: -1
          - Rotations: C2, C3, C4, C6
          - Mirror rotations: S2, S3, S4, S6 */
      SxString           identifier;

      /** \brief Standard constructor

          \par Example:
\code
// 4 fold axis parallel to the z axis
SxSymType C4 (SxSymType::Rotation, 4, SxVector3<Double>(0,0,1));

// xz mirror plane
SxSymType M (SxSymType::Mirror, 1, SxVector3<Double> (0,1,0));
\endcode */
      SxSymType (Classification c_=None, 
                 int axisCount_=0, 
                 const SxVector3<Double> &opCoord_=SxVector3<Double>());
      /** \brief copy constructor
      
          The copy constructor copies all elements of the input except
          the identification member. It is recomputed. Therefore the copy
          constructor calls SxSymType::updateIdentifier.  */
      SxSymType (const SxSymType &);

      /** \brief Comparison operator

          This operator compares the classification only
          \par Example:
\code
SxSymType S (SxSymType::Inversion);
if (S == SxSymType::Rotation)  {
   ...
}
\endcode
      */
      bool operator== (Classification)    const;
      /** \brief Comparison operator

          This operator compares two SxSymType objects. Note, that the
          SxSymType::identifier member is not used for comparison. */
      bool operator== (const SxSymType &) const;

      /** @{
        
          \brief Inequality comparison operator

          This operator has been implemented in order to allow SxSymType
          to be used in SxList and SxSortedList classes.
          The sort order is
          -# Mirror Rotation, axis count = 6
          -# Rotation, axis count = 6
          -# Mirror Rotation, axis count = 4
          -# Rotation, axis count = 4
          -# Mirror Rotation, axis count = 3
          -# Rotation, axis count = 3
          -# Mirror Rotation, axis count = 2
          -# Rotation, axis count = 2
          -# Mirror
          -# Inversion
          -# Identity */
      inline bool operator<  (const SxSymType &in) const
      { 
         return getRank () < in.getRank (); 
      }
      bool operator>  (const SxSymType &in) const
      { 
         return getRank () > in.getRank (); 
      }

      bool isSelfInverse () const
      {
         switch (classification)  {
            case Identity:   return true;
            case Mirror:     return true;
            case Inversion:  return true;
            case MirrorRotation: // same as Rotation
            case Rotation: {
                              SX_CHECK (axisCount >= 2, axisCount);
                              return axisCount == 2;
                           }
            default:         SX_EXIT;
         }
         return false;
      }
      
   protected:
      /// Map symmetry type to number (for < and >)
      int getRank () const;
      /** @} */

   protected:

      /** \brief Recompute human readable symmetry identifier

          This routine sets up the identifier member according to the 
          Schoenfliess nomenclature.
          \sa SxSymType::identifier  */
      void updateIdentifier ();
};

#endif /* _SX_SYM_TYPE_H_ */
