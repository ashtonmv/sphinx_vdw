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

#ifndef _SX_DIELEC_H_
#define _SX_DIELEC_H_

#include <SxDirac.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxPerturbK.h>
#include <SxExt.h>

/** \brief Dielectric properties

    \b SxClass = S/PHI/nX ...

    ....

    \ingroup  group_addons
    \author   Christoph Freysoldt freyso@fhi-berlin.mpg.de */
class SX_EXPORT_EXT SxDielec
{
   public:
      /// Constructor
      SxDielec () : beta(0.) { /* empty */ }

      /// Destructor
      ~SxDielec () { /* empty */ }

      /// Calculate dielectric tensor in cartesian coordinates
      SxMatrix3<TPrecCoeffG> compute (const SxPW    &waves, 
                                      const SxFermi &fermi,
                                      const SxPerturbK &kp,
                                      double imagFreq,
                                      const SxCell &cell,
                                      bool   computeWing);

      /// Calculate dielectric tensor in cartesian coordinates
      void compute (const SxPW    &waves, 
                                  const SxFermi &fermi,
                                  const SxPerturbK &kp,
                                  SxArray<double> imagFreq,
                                  const SxCell &cell,
                                  bool   computeWing);

      /// The head of the dielectric matrix (nOmega, 3, 3)
      SxArray<SxMatrix3<Double> >  head;
      /** The wings of the dielectric matrix (nOmega, meshSize, 3)
        The wings are purely imaginary, therefore only the imaginary
        part is stored.
        */
      SxArray<SxDiracMat<Double> > wings;

      /// Inverse electronic temperature for finite-T calculations
      double beta;

};

#endif /* _SX_DIELEC_H_ */
