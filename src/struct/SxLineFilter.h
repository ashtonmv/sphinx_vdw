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

#ifndef _SX_LINE_FILTER_H_
#define _SX_LINE_FILTER_H_

#include <SxOperator.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxStruct.h>

/** \brief Apply constraints to forces

    \b SxLineFilter = S/PHI/nX Line constraint filter

    This filter is used to project selected coordinates onto certain
    lines.
    Usually it will be applied to forces in order to apply line constraints
    in structure optimization routines.

    \ingroup group_structure;
    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_STRUCT SxLineFilter
: public SxOperatorBase<SxAtomicStructure> SXOP_LINKFIX
{
   public:

      SxLineFilter ();
      SxLineFilter (const SxSymbolTable *, const SxAtomicStructure &);
      virtual ~SxLineFilter ();

      virtual SxAtomicStructure operator* (const SxAtomicStructure &) const;
      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT; }


   protected:
      SxAtomicStructure lines;

      SXOPERATOR_GETCOPY (SxLineFilter, SxAtomicStructure);
   public:
      /// True if filter doesn't do anything
      bool isEmpty () const { return lines.getSize () == 0; }
};

#endif /* _SX_LINE_FILTER_H_ */
