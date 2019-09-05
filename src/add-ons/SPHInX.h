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



#ifndef _SPHINX_H_
#define _SPHINX_H_

#include <SxString.h>
#include <SxExt.h>

/** \brief SPHInX main executable

    This "add-on" is the SPHInX program. It contains the main entrance
    point main().

    \ingroup  group_addons
    \author   All SPHInX developers, see src/AUTHORS */
class SX_EXPORT_EXT SPHInX
{
   public:
      SPHInX ();
      SPHInX (const SxString &inFile);
      ~SPHInX ();

      int compute ();

   protected:
      SxString file;
};


#endif /* _SPHINX_H_ */
