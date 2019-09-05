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

#ifndef _SX_RADIALR_BASIS_H_
#define _SX_RADIALR_BASIS_H_

#include <SxPrecision.h>
#include <SxQuantumNumbers.h>
#include <SxBasis.h>
#include <SxDirac.h>
#include <SxArray.h>
#include <SxTimer.h>
#include <SxDFT.h>
#include <SxCubicSpline.h>


class SX_EXPORT_DFT SxRadRBasis : public SxBasis 
{
   public:

      enum gridType {Linear, Quadratic, Mixed, Hyperbel, Logarithmic};

      typedef TRadRBasisType               TBasisType;
      typedef SxDiracVec<TRadRBasisType>   TPsi;

      SxRadRBasis ();
      ~SxRadRBasis ();
      SxRadRBasis (double rMin, double rMax, int nPoints, enum gridType mode = Linear);
      SxRadRBasis (const SxDiracVec<TBasisType> &basis);
      SxRadRBasis (const SxDiracVec<TBasisType> &drIn, double factor, double r0);
      void set (double rMin, double rMax, int nPoints, enum gridType mode);
      void set (const SxDiracVec<TBasisType> &basis);
      void set (const SxDiracVec<TBasisType> &drIn, double factor, double r0);
      const SxDiracVec<TBasisType> & getRadRFunc () const;
      const SxDiracVec<TBasisType> & getRadDR () const;

      virtual SxString getType () const { return "|r>"; };
      virtual ssize_t getNElements () const { return r.getSize (); };

      REGISTER_PROJECTOR (SxRadRBasis, SxRadRBasis, identity);
      REGISTER_PROJECTOR (SxRadRBasis, SxRadBasis, toRadBasis);
      REGISTER_PROJECTOR (SxRadRBasis, SxRadGBasis, toRadGBasis);

      SxDiracVec<TBasisType> identity (const SxRadRBasis *basisPtr,
            const TPsi &vec) const;
      SxDiracVec<TRadBasisType> toRadBasis (const SxRadBasis *radBasisPtr,
            const TPsi &vec) const;
      SxDiracVec<TRadBasisType> toRadGBasis (const SxRadGBasis *radBasisPtr,
            const TPsi &vec) const;
      // int(... R^2dR) 
      double integrate (const TPsi &integrand, bool useSpline = true) const;
      inline double simpsonTerm (double x1, double x2, double x3, double x4, double f1, double f2, double f3, double f4) const;
      inline double simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const;
      inline double simpsonTerm (double x1, double x2, double f1, double f2) const;
      

      SxDiracVec<Double> toSpline (const SxDiracVec<Double> &vec) const;
      SxDiracVec<Double> toVec (const SxDiracVec<Double> &vec) const;
      
      static SxDiracVec<Double> jsb (int l, const SxDiracVec<Double> &z);

   protected:

      SxDiracVec<Double> r;
      SxDiracVec<Double> dr;
};

namespace Timer {
   enum RadRBasisTimer {
      radR2RadG,
      radR2Rad,
      integrateR,
      VecToSplineR,
      SplineToVecR
   };
}

SX_REGISTER_TIMERS (Timer::RadRBasisTimer)
{
   using namespace Timer;
   regTimer (radR2RadG,   "RadR to RadG Basis");
   regTimer (radR2Rad,    "RadR to Rad Basis");
   regTimer (integrateR,  "radR Integration");
   regTimer (VecToSplineR,"Vector to Spline R");
   regTimer (SplineToVecR,"Spline to Vector R");
}

#endif /* _SX_RADIALR_BASIS_H_ */
