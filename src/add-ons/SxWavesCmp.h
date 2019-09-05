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

#ifndef _SX_WAVES_CMP_H_
#define _SX_WAVES_CMP_H_

#include <SxArray.h>
#include <SxString.h>
#include <SxMatrix.h>
#include <SxGkBasis.h>
#include <SxPW.h>
#include <SxExt.h>

/** \brief Compares

    \b SxWaves = SPHInX waves.sxb file compare

    ....

    \ingroup  group_addons
    \author   Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_EXT SxWavesCmp
{
   public:

      SxWavesCmp ();
     ~SxWavesCmp ();

      void readWaves1 (const SxString &);
      void readWaves2 (const SxString &);

      SxMatrix<Double> 
         getOverlap (const SxArray<int> &states, int iSpin, int ik) const;

      SxPW getDifferences ();

      SxDiracVec<Double> absSqr (int i, int iSpin, int ik) const;

   protected:

      SxGkBasis Gk;
      SxPW      waves1, waves2;

};

#endif /* _SX_WAVES_CMP_H_ */
