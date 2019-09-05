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

#ifndef _SX_RHO_FILTER_H_
#define _SX_RHO_FILTER_H_

#include <SxDFT.h>
#include <SxGBasis.h>
#include <SxDensity.h>

/** \brief Density filters

    When PAW densities are mixed with preconditioners, the problem arises
    how to include the on-site charges in the preconditioning. We solve this
    problem by separating the density in a well-behaved "mixing" part and a
    possibly problematic rest, that should not be preconditioned or
    Pulay-mixed.

    This class filters out the mixable part, which is a pure pseudo-density
    with a low cut-off. Alternatively, the filtering can be switched off
    (mode = NoFilter).

    TODO: There are a lot of optimization opportunities, e.g. exploiting
    the low cutoff, working in G-space, etc.

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxRhoFilter
{
   public:
      /// Possible Filter modes
      enum Mode { NoFilter, PseudoRho };
      /// Current filter mode
      Mode mode;

      /// G cutoff
      double gCut;

      /// The G-basis used for Fourier-filtering
      SxPtr<SxGBasis> gBasisPtr;

      /// Constructor
      SxRhoFilter (const SxSymbolTable *table);

      /// Constructor
      SxRhoFilter (Mode modeIn = NoFilter) : mode(modeIn), gCut (10.) {}

      /// Constructor
      SxRhoFilter (double gCutIn) : mode(PseudoRho), gCut (gCutIn) {}

      /// Apply the filter
      SxDensity operator| (const SxDensity &in);
};

#endif /* _SX_RHO_FILTER_H_ */
