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

#include <SxError.h>
#include <stdio.h>
#include <iostream>
#include <SxIO.h>
#include <SxTypes.h>
#include <SxUtil.h>

/** \brief Class for reading special file formates used in DFTB Slater-Koster files  

    \b SxDFTBFilesReader = S/PHI/nX DFTB Slater-Koster files reading assitant 

    This class allows to read files that contain numbers in the following format: 
    n*number, where n is the repettion. This format is used in the
    slater-koster files used in the DFTB code (see www.dftb.org)

    \author Christoph Freysoldt, freysoldt@mpie.de and Hazem Abu-Farsakh, h.farsakh@mpie.de
*/

class SxDFTBFilesReader
{
   public:
      FILE *fpIn;

      SxDFTBFilesReader (FILE *inFile = NULL)
         : fpIn(inFile), occurence (0)
      {
         // empty
      }

      bool   hasMoreReals () { return occurence != 0; }
      double getNextReal ();

   protected:
      int occurence;
      double val;
};

