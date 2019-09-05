//     
//      Contact:    Sixten Boeck, boeck@mpie.de
//                  Algorithm Design and Modeling Group
//                  Computational Materials Design
//                  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//     
//      Authors:    see src/AUTHORS
//     
// ---------------------------------------------------------------------------

#ifndef _SX_RADIALG_BASIS_H_
#define _SX_RADIALG_BASIS_H_

#include <SxPrecision.h>
#include <SxQuantumNumbers.h>
#include <SxBasis.h>
#include <SxGBasis.h>
#include <SxDirac.h>
#include <SxArray.h>
#include <SxTimer.h>
#include <SxDFT.h>
#include <SxCubicSpline.h>

typedef SxArray<SxDiracMat<Complex16> > SxGkVec; //(ik,iOrbital:ig)

class SX_EXPORT_DFT SxRadGBasis : public SxBasis 
{
   public:

      enum gridType {Linear, Quadratic, Mixed, Hyperbel, Logarithmic};

      typedef TRadGBasisType               TBasisType;
      typedef SxDiracVec<TRadGBasisType>   TPsi;

      SxRadGBasis ();
      ~SxRadGBasis ();
      SxRadGBasis (double gMin, double gMax, int nPoints, enum gridType mode = Linear);
      SxRadGBasis (const SxDiracVec<Double> &basis);
      SxRadGBasis (const SxDiracVec<Double> &dgIn, double factor, double g0);
      void set (double gMin, double gMax, int nPoints, enum gridType mode);
      void set (const SxDiracVec<Double> &basis);
      void set (const SxDiracVec<Double> &dgIn, double factor, double g0);
      const SxDiracVec<Double> &getRadGFunc () const;
      const SxDiracVec<Double> &getRadDG () const;

      virtual SxString getType () const { return "|g>"; };
      virtual ssize_t getNElements () const { return g.getSize (); };

      REGISTER_PROJECTOR (SxRadGBasis, SxRadGBasis, identity);
      REGISTER_PROJECTOR (SxRadGBasis, SxRadBasis, toRadBasis);
      REGISTER_PROJECTOR (SxRadGBasis, SxRadRBasis, toRadRBasis);
      REGISTER_PROJECTOR (SxRadGBasis, SxGBasis, toGBasis);

      SxDiracVec<TBasisType> identity (const SxRadGBasis *basisPtr,
                                       const TPsi &vec) const;
      SxDiracVec<TRadBasisType> toRadBasis (const SxRadBasis *radBasisPtr,
                                            const TPsi &vec) const;
      SxDiracVec<TRadBasisType> toRadRBasis (const SxRadRBasis *radRBasisPtr,
                                            const TPsi &vec) const;
      SxDiracVec<TGBasisType> toGBasis (const SxGBasis *gBasisPtr,
                                        const TPsi &vec) const;
      // int(... G^2dG) 
      double integrate (const TPsi &integrand, bool useSpline = true) const;
      inline double simpsonTerm (double x1, double x2, double x3, double x4, double f1, double f2, double f3, double f4) const;
      inline double simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const;
      inline double simpsonTerm (double x1, double x2, double f1, double f2) const;

      /// Trace implementation: integrate over radial space
      virtual Real8 tr(const SxDiracVec<Double> &x) const
      {
         return integrate (x, false);
      }
      

      SxDiracVec<Double> toSpline (const SxDiracVec<Double> &vec) const;
      SxDiracVec<Double> toVec (const SxDiracVec<Double> &vec) const;

      static SxDiracVec<Double> jsb (int l, const SxDiracVec<Double> &z);

   protected:

      SxDiracVec<Double> g;
      SxDiracVec<Double> dg;
};

namespace Timer {
   enum RadGBasisTimer {
      radG2Rad,
      radG2RadR,
      radG2Gk,
      integrateG,
      VecToSplineG,
      SplineToVecG,
      phaseFactors,
      splineGetY
   };
}

SX_REGISTER_TIMERS (Timer::RadGBasisTimer)
{
   using namespace Timer;
   regTimer (radG2Rad,   "RadG to Rad Basis");
   regTimer (radG2RadR,  "RadG to RadR Basis");
   regTimer (radG2Gk,    "RadG to Gk Basis");
   regTimer (integrateG, "RadG Integration");
   regTimer (VecToSplineG,"Vector to Spline G");
   regTimer (SplineToVecG,"Spline to Vector G");
   regTimer (phaseFactors,"calc Phasefactors");
   regTimer (splineGetY,  "Spline to |G+k|");
}

#endif /* _SX_RADIALG_BASIS_H_ */
