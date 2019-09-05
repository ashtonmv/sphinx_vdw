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

#ifndef _SX_EQUIVALENCE_H_
#define _SX_EQUIVALENCE_H_

#include <SxGeom.h>
#include <SxAtomicStructure.h>

/** \brief Manage symmetry-equivalent atoms in structures


    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxEquivalence
{
   protected:
      /// Primitive structure
      SxAtomicStructure primStr;

   public:
      /// List of matrix symmetries (iEqClass:iSym)
      SxArray<SxArray<int> > matSym;

      /// Map atoms to equivalent class
      SxArray<int> equivId;

      /// Which rotation brings representative atom to this atom?
      SxArray<int> mappingRot;

      /// Get number of equivalence classes
      inline ssize_t getSize () const { return matSym.getSize (); }

      /// Setup from a primitive structure
      void setup (const SxAtomicStructure &primStrIn);

      /// Setup from the primitive structure
      void reset () { setup (primStr); }

      /// Print in sx format
      void fprintsx (FILE *fp) const;

      /// Read from symbol table
      void read (const SxSymbolTable *table);

      /// Empty constructor
      SxEquivalence () {}

      /// Constructor from a primitive structure
      SxEquivalence (const SxAtomicStructure &primStrIn) { setup (primStrIn); }

      /** \brief Transfer primitive structure's equivalence relations to a more
          complex structure

          This routine maps the complex structure's atoms to the primitive
          structure, assigning the proper equivalence relations.
          */
      SxEquivalence mapToComplex (const SxAtomicStructure &supercell,
                                  double dist = 1e-3) const;

      /// Map symmetries to an alternative list of symmetries
      void mapSyms (const SxArray<SymMat> &newSym);
};

#endif /* _SX_EQUIVALENCE_H_ */
