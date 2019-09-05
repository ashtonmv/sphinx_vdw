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

#ifndef _SX_RS_PROJ_H_
#define _SX_RS_PROJ_H_

#include <SxDFT.h>
#include <SxNaturalCubicSpline.h>
#include <SxRBasis.h>
#include <SxRadGBasis.h>

/** \brief Real-space projectors (using spherical harmonics)

    \b SxClass = S/PHI/nX real-space projectors

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxRSProj
{
   public:
      /// The real-space functions (:is)(:ipt)
      SxArray<SxArray<SxNaturalCubicSpline> > rProj;

      /// Local projector quantum numbers
      class ProjQuantumNumbers {
         public:
            int n, l, m;
            ProjQuantumNumbers (int nIn = -1, int lIn = -1, int mIn = -1)
               : n(nIn), l(lIn), m(mIn) 
            { }
      };

      /// Map of projectors (:is)(:ipl)
      SxArray<SxArray<ProjQuantumNumbers> > projId;

      /// Number of projectors
      int npg;

      /// Real-space cutoffs (:is)
      SxArray<double> rCut;

      /// Real-space basis
      SxPtr<SxRBasis> rBasisPtr;

      /// Max. l value (:is)
      SxArray<int> lmax;

      /// Radial mesh resolution
      double dr;

      /// Number of negative radial points in spline
      int nrn;

      /// Smoothening beta
      double beta;

      /// Project
      SxDiracVec<Complex16> project (const SxDiracVec<TPrecCoeffG> &psi) const;

      /// Gradient
      SxDiracVec<TPrecCoeffG> gradient (const SxDiracVec<Complex16> &pPsi,
                                        const SxGBasis &gk) const;

      /// Max. number of states
      int nPsiMax;
   protected:
      /// Project
      SxDiracVec<Complex16> project (const SxAtomicStructure &str,
                                     const SxDiracVec<Complex16> &uk,
                                     const Coord &kVec) const;
      /// Multi-state project
      SxDiracVec<Complex16> project (const SxAtomicStructure &str,
                                     const SxArray<SxDiracVec<Complex16> > &uk,
                                     const Coord &kVec) const;
      /// Gradient
      SxDiracVec<Complex16> gradient (const SxAtomicStructure &str,
                                      const SxDiracVec<Complex16> &pPsi,
                                      const Coord &kVec) const;
      /// Gradient
      SxArray<SxDiracVec<Complex16> >
      gradientN (const SxAtomicStructure &str,
                 const SxDiracVec<Complex16> &pPsi,
                 const Coord &kVec) const;

      /// Hardening factor (cached)
      mutable SxDiracVec<Double> antiGauss;
   public:

      /// Constructor
      SxRSProj (const SxAtomicStructure &str,
                const SxArray<SxArray<SxDiracVec<Double> > > pIn,
                const SxRadGBasis &radG,
                const SxPtr<SxRBasis> &rBasisIn,
                double betaIn,
                double drIn, double rMax, double dPhi,
                int nrnIn);

      /// Constructor
      SxRSProj (const SxSymbolTable *table,
                const SxAtomicStructure &str,
                const SxArray<SxArray<SxDiracVec<Double> > > pIn);

      /// Internal setup routine
      void setup (const SxAtomicStructure &str,
                  const SxArray<SxArray<SxDiracVec<Double> > > pIn,
                  const SxRadGBasis &radG,
                  double rMax,
                  double dPhi,
                  double eCut = -1.);

};

#endif /* _SX_RS_PROJ_H_ */
