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

#ifndef _SX_RSPROJ_H_
#define _SX_RSPROJ_H_

#include <SxFFT1d.h>
#include <SxPtr.h>
#include <SxRBasis.h>
#include <SxTimer.h>
#include <SxAtomicStructure.h>
#include <SxDFT.h>

/** \brief Real-space projectors

    \b SxClass = S/PHI/nX real-space projectors

    This class provides real-space projector functionality via
    the Euler Exponential spline in G-space method (NLEES-G).

      --- ref1: H.-S. Lee, M. E. Tuckerman, G. J. Martyna
                Efficient Evaluation of Nonlocal Pseudopotentials via
                Euler Exponential spline interpolation
                ChemPhysChem 2005, 6, 1827-1835


    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxEESGProj
{
   protected:
      /** The order of the cardinal B splines */
      int splineOrder;

      /** The real-space basis for interpolation (s-space) */
      SxPtr<SxRBasis> sBasisPtr;

      /** Phase-like reciprocal-space co-factors in 1D */
      SxVector<Complex16> dpx, dpy, dpz;

      /// cached G-basis
      mutable const SxGBasis *cachedG;
      /// cached G vectors in packed form
      mutable SxArray<ssize_t> packedGrel;
      /// Update packedGrel
      void updateG (const SxGBasis *gPtr) const;
      
   public:
      /** \brief Constructor
          @param cell  Real-space cell
          @param mesh  s-space mesh
          @param order spline order
        */
      SxEESGProj (const SxCell &cell, const SxVector3<Int> &mesh, int order);

      /** Get the s-space basis
          \note The |G+k>-bases need to register this basis.
       */
      const SxRBasis &getSBasis () const { return *sBasisPtr; }

      /** \brief Perform projections via EES-G

          @param refPhi  right side function
          @param psi     left side function
          @param centers centers for phi
          @return the projections Z

          \f[
          Z(\mathbf x) =
          \langle \phi| e^{i \mathbf G \cdot \mathbf x} | \psi \rangle
          \f]

          \note The |G+k>-basis needs to have registered the s-basis.
        */
      SxVector<Complex16> project (const SxDiracVec<Complex16> &refPhi,
                                   const SxDiracVec<Complex16> &psi,
                                   const SxAtomicStructure     &centers) const;
      /** \brief Perform projections + their position gradients via EES-G

          @param refPhi  right side function
          @param psi     left side function
          @param centers centers for phi
          @return 4 x nCenters matrix, containing the projections Z (row 0)
                  and the gradient with respect to center positions (row 1-3)

          \f[
          Z(\mathbf x) =
          \langle \phi| e^{i \mathbf G \cdot \mathbf x} | \psi \rangle
          \f]

          \note The |G+k>-basis needs to have registered the s-basis.
        */
      SxVector<Complex16>
      projectGrad (const SxDiracVec<Complex16> &refPhi,
                   const SxDiracVec<Complex16> &psi,
                   const SxAtomicStructure     &centers) const;
      /** \brief Get projection gradient via EES-G

          @param proj    projections Z(x)
          @param refPhi  right side function
          @param centers centers for phi
          @return Gradient of projections wrt conjugate plane-wave coefficients

          \f[
          \sum_{\mathbf x} 
          \frac{\partial Z^*(\mathbf x)}{\partial \langle\psi|\mathbf G>}
          Z(\mathbf x)
          \f]

          \note The |G+k>-basis needs to have registered the s-basis.
        */
      SxDiracVec<Complex16> gradient (const SxVector<Complex16>   &proj,
                                      const SxDiracVec<Complex16> &refPhi,
                                      const SxAtomicStructure     &centers)
                                      const;

      /** \brief Compute cardinal B spline values for offset xx
           @param xx  offset, between 0 and 1
           @param p   spline order
           @param res preallocated return array

           Returns res(i) = cardBSpline(xx+i,p) (i < p)
        */
      static void cardBSplineAll (double xx, int p, SxVector<Double> &res);

      /** \brief Compute cardinal B spline values for some x
           @param x    argument, between 0 and p
           @param p    spline order
           @param work preallocated work space
           @return Cardinal B spline of order p at x

           \note This is a wrapper to cardBSplineAll. If you need
           the values for several x, it is more efficient to get
           all the values at once.
        */
      static double cardBSpline (double x, int p, SxVector<Double> &work);
      
      /** \brief This computes exp(i kr) */
      static inline SxComplex16 getPhase (double kr)
      {
         SxComplex16 res;
         sincos(kr, &res.im, &res.re);
         return res;
      }

      /** \brief Compute the reciprocal space phase-like co-factor \f$d_p\f$
                 in 1D
        @param p    spline order
        @param fft  a 1D FFT of size N
        @return     d_p(g) for g = 0 .. N-1

        \f[
        d_p (g) = e^{2\pi i g p / N} 
                  [ \sum_s M_p (s) e^{2 \pi i g s / N} ]^{-1}
        \f]
        cf. Eq. 11 of ref. 1 
        */
      static SxVector<Complex16> getDp (int p, SxFFT1d &fft);

      enum Timer { 
         /// Computation of cardinal B splines
         CardBSpline,
         /// Setup of \f$d_p\f$
         DpSetup,
         /// Projection
         Projection,
         /// Gradient
         Gradient,
         /// phi * psi * dp
         PhiPsiDp,
         /// Q*phi*dp
         GradientMul,
         /// Mp setup in projection
         ProjMp,
         /// FFT in projection
         GtoS,
         /// FFT in projection
         StoG,
         /// Summation
         ProjSumS,
         /// Q(s) setup
         GradientQ
      };
};

SX_REGISTER_TIMERS (SxEESGProj::Timer)
{
   regTimer (SxEESGProj::CardBSpline, "Mp splines");
   regTimer (SxEESGProj::DpSetup,     "dp setup");
   regTimer (SxEESGProj::Projection,  "<phi|psi>");
   regTimer (SxEESGProj::Gradient,    "gradient <phi|psi>");
   regTimer (SxEESGProj::PhiPsiDp,    "phi*psi*dp");
   regTimer (SxEESGProj::GradientMul, "Q*phi*dp");
   regTimer (SxEESGProj::ProjMp,      "<phi|psi> Mp");
   regTimer (SxEESGProj::GtoS,        "<phi|psi> G->s");
   regTimer (SxEESGProj::StoG,        "<phi|psi> s->G");
   regTimer (SxEESGProj::ProjSumS,    "<phi|psi> sum");
   regTimer (SxEESGProj::GradientQ,   "Q(s)");
}

#endif /* _SX__H_ */
