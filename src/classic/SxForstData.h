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

#ifndef _SX_FORST_DATA_H_
#define _SX_FORST_DATA_H_

#include <SxClassic.h>
#include <SxSpeciesData.h>
#include <SxSymbolTable.h>

/** \brief ...

    \b SxForstData = S/PHI/nX ...

    ....

    \author Timo Jacob, jacob@fhi-berlin.mpg.de */
class SX_EXPORT_CLASSIC SxForstData : public SxSpeciesData
{
   public:

      /** \brief ...

          ...
          */
      SxForstData ();

      SxForstData (const SxSymbolTable *);

      /** \brief Destructor */
      ~SxForstData ();

   protected:

      // here protected members 

};

#endif /* _SX_FORST_DATA_H_ */
