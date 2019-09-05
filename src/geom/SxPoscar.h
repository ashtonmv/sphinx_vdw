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

#ifndef _SX_POSCAR_H_
#define _SX_POSCAR_H_

#include <SxGeom.h>
#include <SxList.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxPrecision.h>
#include <SxString.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

#include <stdlib.h>
#include <stdio.h>

/** \brief S/PHI/nX to VASP Structure Converter

    \b SxPoscar = S/PHI/nX to VASP-POSCAR Atomic Structure converter

    \author Sixten Boeck, boeck@mpie.de
    */
class SX_EXPORT_GEOM SxPoscar
{
   public:

      SxPoscar ();
      SxPoscar (const SxString &file_);

      ~SxPoscar ();

      /** Specifies the POSCAR filename.
          \param file_ Specifies the filename for read(). */
      void setFilename (const SxString &file_);

      /** Returns the atomic structure.
        */
      SxAtomicStructure getStructure () const;

      SxArray<SxString> getUniqueSpecies() const;

      /** Read from a POSCAR file. */
      void read ();

      /** Write data to POSCAR file 
        */
      void write (const SxAtomicStructure &structure_,
                  const SxArray<SxString> &species);

   protected:

      SxString filename;

      // --- frames
      SxAtomicStructure structure;

      SxArray<SxString> chemNames;

};

#endif /* _SX_POSCAR_H_ */
