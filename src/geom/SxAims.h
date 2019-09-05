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

#ifndef _SX_AIMS_H_
#define _SX_AIMS_H_

#include <SxGeom.h>
#include <SxList.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxPrecision.h>
#include <SxString.h>
#include <SxAtomicStructure.h>

#include <stdlib.h>
#include <stdio.h>

/** \brief AIMS geometry Converter

    \b SxAims = S/PHI/nX Aims File Format Converter

    \author Björn Lange, bjoern.lange@duke.edu */
class SX_EXPORT_GEOM SxAims
{
   public:

      SxAims ();
      SxAims (const SxString &file_);

      ~SxAims ();

      /** Specifies the Aims filename.
          \param file_ Specifies the filename for read(). */
      void setFilename (const SxString &file_);

      /** Returns the atomic structure. */
      SxAtomicStructure &getStructure ();

      /** Returns unique species. */
      SxList<SxString> &getUniqueSpecies ();

      /** Read Aims file. */
      void read ();

      /** Write data to Aims file 
        */
      void write (const SxAtomicStructure &structure_,
                  const SxArray<SxString> &chemNames_);

   protected:

      SxString filename;

      SxAtomicStructure structure;
      SxList<SxString> uniqueSpecies;

};

#endif /* _SX_AIMS_H_ */
