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

#ifndef _SX_FOURIER_INTERPOL_H_
#define _SX_FOURIER_INTERPOL_H_

#include <SxCell.h>
#include <SxFFT3d.h>
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxAtomicStructure.h>
#include <SxDFT.h>

/** \brief Fourier interpolation



    \author Christoph Freysoldt freyso@fhi-berlin.mpg.de */
class SX_EXPORT_DFT SxFourierInterpol
{
   public:
      /// The mesh in reciprocal space
      SxVector<TPrecCoeffG> meshInG;
      /// The G-vectors in relative coordinates
      SxMatrix<Int> gVec;
      /// The real space cell
      SxCell cell;

      /// \name Constructors
      //@{
      /// Empty constructor
      SxFourierInterpol () {/* empty */}
      /// From mesh file
      SxFourierInterpol (SxBinIO &);
      //@}

      /// Print
      void print (ostream & = cout);
      
      /// Interpolate with discrete Fourier transform
      SxVector<Double> interpolate (const SxAtomicStructure &rVec);
      
      /// Remove coefficients with absSqr < threshold
      void condense (double threshold = 1e-20);
};

#endif /* _SX_FOURIER_INTERPOL_H_ */
