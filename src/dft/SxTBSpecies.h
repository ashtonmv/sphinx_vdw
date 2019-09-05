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


#ifndef _SX_TB_SPECIES_H_
#define _SX_TB_SPECIES_H_

#include <SxSymbolTable.h>
#include <SxSpeciesData.h>
#include <SxArray.h>

/** \brief Species Data for Tight-Binding (DFTB)

    \b SxTBSpecies = SPHInX Tight-Binding Species 

    \author Hazem Abu-Farsakh h.farsakh@mpie.de */

class SxTBSpecies : public SxSpeciesData
{
   public:
      SxTBSpecies ();
      SxTBSpecies (const SxSymbolTable *table);
      ~SxTBSpecies ();

      /** brief User-provided name for the tight-binding species. They
        are used to construct the slater-koster file names */
      SxArray<SxString> skElementName;

      /** brief Maximum angular momentum. For the moment this has no effect.
                The value will be taken from the SK (potential) files.  */
      SxVector<Int>     lMax;

      /** brief Path of the tight-binding slater-koster files */
      SxString skFilesPath;

      /** brief format of the tight-binding slater-koster files 
                1: original format of the dftb code (Frauenheim's group) 
                2: S/PHI/nX format */
      int      skFilesFormat;
};

#endif /* _SX_TB_SPECIES_H_ */
