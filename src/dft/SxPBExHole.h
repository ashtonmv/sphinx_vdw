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

#ifndef _SX_PBE_X_HOLE_H_
#define _SX_PBE_X_HOLE_H_

#include <SxDFT.h>
#include <SxGaussQuad.h>

// References
// (1) J. Chem. Phys. 109,3313 (1998)
// (2) HSE implementation notes

/** \brief PBE exchange hole model

    This class yields the model for the PBE exchange hole as described in
    J. Chem. Phys. 109,3313 (1998).

    It is needed also for the HSE functional.

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPBExHole
{
   public:
      /// Constant model parameters
      static const double A;
	  static const double B;
	  static const double C;
	  static const double D;
	  static const double E;
      /// Coefficients for H(s)
      static const double a1;
	  static const double a2;
	  static const double a3;
	  static const double a4;
	  static const double a5;

      /// s-dependent model parameters
      double F, G, H;

      /// s-dependent model parameters: s-derivative
      double Fs, Gs, Hs;

      /// Compute H from interpolating analytic fit
      inline void computeH (double s, double s2, double s4)  {
         // Eq. A5
         H = (a1 + a2 * s2) * s2
           / (1. + s4 * (a3 + s * a4 + s2 * a5));
      }

      /// Compute H from interpolating analytic fit
      inline double computeH (double s)  {
         double s2 = s * s;
         computeH (s, s2, s2 * s2);
         return H;
      }
      
      /// Compute dH/ds from interpolating analytic fit
      void computeHs (double s, double s2, double s4);

      /// Compute s-dependent model parameters F, G, and H
      void computeFGH (double s, bool deriv = false);
      
      /// Compute PBE exchange hole
      double exchangeHolePBE(double s, double y);

      /// Testing routine
      void test (int argc, char **argv);

      /// Screened exchange enhancement factor
      double Fx,
      /// Screened exchange enhancement factor s-derivative
             Fxs,
      /// Screened exchange enhancement factor nu-derivative
             FxNu;
      
      /// Screened exchange enhancement factor
      double screenedFx (double s, double nu);

      /// Numerical integral from Henderson paper
      double hendersonK (double alpha, double beta);
   private:
      /// Gauss-Legendre quadrature
      SxGaussQuad gaussLegendre, 
      /// Gauss-Laguerre quadrature
                  gaussLaguerre;

      /// Auxiliary integral for Fx (4-th order approximation to exp(-x)
      static double integralP4 (double alpha, double beta, double z);

      static double P4 (double x)
      {
         return 1. - x * (1. - 0.5 * x * (1. - x/3. * (1. - 0.25 * x)));
      }
};

#endif /* _SX_PBE_X_HOLE_H_ */
