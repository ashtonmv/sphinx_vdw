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

#ifndef _SX_SPECIES_REF_H_
#define _SX_SPECIES_REF_H_
#include <SxOperator.h>
#include <SxAtomicStructure.h>
#include <SxGeom.h>

/** \brief Return on-the-fly reference to atoms of one species

    \b SxSpeciesRef = S/PHI/nX species on-the-fly reference

    In some algorithms, all atoms from one species must be
    selected. This filter provides a means to do this on-the-fly
    efficiently.
    
    \note
    The resulting atomic structure must only be used as a
    temporary object in computational expressions.
    No direct assignment, no initialization.
    \example
    \code
deltaSqr = delta.absSqr ();
for (is = 0; is < delta.getNSpecies (); ++is)
   forces += z(is) * (SxSpeciesRef(is) | delta) 
                   / deltaSqr (delta.getRange(is));
    \endcode

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxSpeciesRef 
	: public SxOperatorBase<SxAtomicStructure> SXOP_LINKFIX
{
   public:
      /// Species to select
      int iSpec;
      /// Which structure to derive from:
      enum ParentSelection { 
         /// Derive result from the filtered structure (default)
         MapToFiltered, 
         /// Derive result from the parent of the filtered structure
         MapToParentOfFiltered 
      } mode;
      
      /// Constructor
      SxSpeciesRef (int i, enum ParentSelection mapping = MapToFiltered) 
         : iSpec(i), mode(mapping) { /* empty */}
      
      /// Filter routine
      virtual SxAtomicStructure operator* (const SxAtomicStructure &in) const;
      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT; }
      /// Standard getCopy
      SXOPERATOR_GETCOPY (SxSpeciesRef, SxAtomicStructure);
};
#endif /* _SX_SPECIES_REF_H_ */
