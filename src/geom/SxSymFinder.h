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

#ifndef _SX_SYM_FINDER_H_
#define _SX_SYM_FINDER_H_

#include <SxAtomicStructure.h>
#include <SxUniqueList.h>
#include <SxLinEquation.h>
#include <SxRotation.h>
#include <SxGeom.h>

/** @brief Symmetrizer class

    @ingroup group_structure
    @author  C. Freysoldt, freyso@fhi-berlin.mpg.de */
class SX_EXPORT_GEOM SxSymFinder  {
   public:
      /// Constructor
      SxSymFinder () {/* empty */}

      /** Constructor with implicit compute (structure)
        \note Do not call compute afterwards.
        \example
        \code
// move to high-symmetry position
structure += SxSymFinder(structure).getHighSymShift ();
// wrap atoms around
structure %= structure.cell;
        \endcode
      */
      inline SxSymFinder (const SxAtomicStructure &structure)
      {
         compute (structure);
      }
      
      /// key: Set of equations that characterize shifts
      /// value: Set of rotation matrices that become symmorphic upon shifts
      typedef SxUniqueList<SxRotation,SxRotation> RotationList;
      SxMap<SxLinEquation, RotationList, SxNull> equations;

      /// This routine runs the symmetry finder
      void compute (const SxAtomicStructure &structure);

      /// Returns a shift vector to a maximum symmetry position
      Coord getHighSymShift () const;

};

#endif /* _SX_SYM_FINDER_H_ */
