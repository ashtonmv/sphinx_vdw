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

#ifndef _SX_XYZ_H_
#define _SX_XYZ_H_

#include <SxGeom.h>
#include <SxList.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxPrecision.h>
#include <SxString.h>
#include <SxAtomicStructure.h>

#include <stdlib.h>
#include <stdio.h>

/** \brief XYZ Converter

    \b SxXYZ = S/PHI/nX XYZ File Format Converter

    \par References:
    \li http://hackberry.trinity.edu/IJC/Text/xmolxyz.html
    \li http://www.ch.ic.ac.uk/chemime/

    \todo SxXYZ::write, load n-th frame

    \author Sixten Boeck, boeck@mpie.de
    \author Vaclav Bubnik, bubnik@mpie.de */
class SX_EXPORT_GEOM SxXYZ
{
   public:

      SxXYZ ();
      SxXYZ (const SxString &file_);

      ~SxXYZ ();

      /** Specifies the XYZ filename.
          \param file_ Specifies the filename for read(). */
      void setFilename (const SxString &file_);

      /** Returns the number of loaded frames.
          \return the number of loaded frames in read(). */
      int getNFrames ();

      /** Returns the atomic structure.
          \param idx_ Specifies the frame [0...getNFrames)
          \return the atomic structure from the specified frame. */
      SxAtomicStructure &getStructure (int idx_);

      /** Returns unique species.
          \param idx_ Specifies the frame [0...getNFrames)
          \return the species from the specified frame. */
      SxList<SxString> &getUniqueSpecies (int idx_);

      /** Read all frames from a XYZ file. */
      void read ();

      /** Read specified frames from a XYZ file.
          \param startFrame_  read from the frame#
          \param endFrame_   read to the frame# (-1 all)
          \param step_ */
      void read (int startFrame_, int endFrame_, int step_);

      /** Write data to XYZ file 
        */
      void write (const SxAtomicStructure &structure_,
                  const SxArray<SxString> &chemName,
                  bool writeAMat=false);

   protected:

      SxString filename;

      // --- frames
      SxArray<SxAtomicStructure> structures;
      SxArray<SxList<SxString> > uniqueSpecies;

};

#endif /* _SX_XYZ_H_ */
