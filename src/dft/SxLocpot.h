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

#ifndef _SX_LOCPOT_H_
#define _SX_LOCPOT_H_

#include <SxDFT.h>
#include <SxString.h>
#include <SxRho.h>
#include <SxAtomicStructure.h>

/** \brief S/PHI/nX to VASP LOCPOT Converter

    \b SxPoscar = S/PHI/nX to VASP-LOCPOT converter

    \author Bjoern Lange, b.lange@mpie.de
    */
class SX_EXPORT_DFT SxLocpot
{
   public:

      SxLocpot ();
      SxLocpot (const SxString &file_);

      ~SxLocpot ();

      /** Specifies the LOCPOT filename.
          \param file_ Specifies the filename for read(). */
      void setFilename (const SxString &file_);

      SxCell getCell () const { return cell;}
      SxMesh3D getMesh () const { return mesh;}
      SxMeshR getPotential () const { return potential;}
      SxAtomicStructure getStructure () const { return structure;}


      /** Read from a LOCPOT file. */
      void read ();

   protected:

      SxCell readVASPCell (FILE *fp, const SxString &file);
      
      SxString filename;

      SxCell cell;

      SxMesh3D mesh;

      SxAtomicStructure structure;

      SxMeshR potential;

};

#endif /* _SX_LOCPOT_H_ */
