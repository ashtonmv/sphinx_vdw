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

#ifndef _SX_PDB_FAST_H_
#define _SX_PDB_FAST_H_

#include <SxList.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxPrecision.h>
#include <SxString.h>
#include <SxAtomicStructure.h>
#include <SxGeom.h>

#include <stdlib.h>
#include <stdio.h>

/** \brief PDB Converter (Faster SxPDB)

    \b SxPDBFast = S/PHI/nX PDB File Format Converter

    \par Filter:
    \li SxAtomicStructure+chemical symbols ->\b SxPDBFast -> file
    \li file ->\b SxPDBFast -> SxAtomicStructure

    \par Large PDB files with SxPDBFast:
    \li 1e6t_pr0039.pdb (18MB, 200.000 atoms)
    \li time to load: 4.85 s

    \par Large PDB files with SxPDB:
    \li 1e6t_pr0039.pdb (18MB, 200.000 atoms)
    \li time to load: 30 min and still waiting (presumption few days)

    \par Large PDB files:
    http://ndbserver.rutgers.edu/ftp/NDB/coordinates/na-replaced/

    \author Sixten Boeck, boeck@mpie.de
    \author Vaclav Bubnik, bubnik@mpie.de */
class SX_EXPORT_GEOM SxPDBFast
{
   public:

      SxPDBFast ();
      SxPDBFast (const SxString &file_);

      ~SxPDBFast ();

      /** Specifies the PDB filename */
      void setFilename (const SxString &file_);

      SxAtomicStructure &getStructure ();
      SxList<SxString> &getUniqueSpecies ();

      /** Read from a PDB file */
      void read ();

      /** Write the data to a PDB file */
      void write (const SxAtomicStructure &structure_,
                  const SxArray<SxString> &chemName);

   private:
      SxString filename;
      SxAtomicStructure  structure;
      SxList<SxString>   uniqueSpecies;
};

#endif /* _SX_PDB_FAST_H_ */
