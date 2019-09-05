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

#ifndef _SX_RADIAL_ATOM_H_
#define _SX_RADIAL_ATOM_H_

#include <SxDFT.h>
#include <SxRadialMesh.h>
#include <SxSphereGrid.h>
#include <SxXC.h>
#include <SxXCFunctional.h>

/** \brief Functionality for DFT atom (in radial coordinates)


    \author C. Freysoldt */
class SX_EXPORT_DFT SxRadialAtom
{
   public:

      /// Constructor
      SxRadialAtom (SxXC::XCFunctional xcIn, 
                    SxSphereGrid::GridType gridTypeIn = SxSphereGrid::Grid_110);

      /** \brief Compute radial Hartree potential component from
                 rho component.

          The expression used is
          \f[
          V_{lm}(r) = \frac{4\pi}{(2l+1)r^{l+1}} 
                    \int_0^r (r')^2 dr' ~(r')^l \rho_{lm}(r')
                    + \frac{4\pi r^l}{(2l+1)}
                    \int_r^{\infty} (r')^2 d r' ~ (r')^{-l} \rho_{lm}(r')
          \f]
          In practice, the first integral is split into the numerical
          part on the radial grid starting from r0, and a correction term
          \f[
          \int_0^r0 (r')^2 dr' ~(r')^l \rho_{lm}(r')
          \approx \rho_{lm}(r0) \frac{r0^{l+3}}{l+3} 
                  - \frac{d\rho_{lm}}{dr} \frac{r0^{l+4}}{(l+4)(l+3)}
          \f] 
          which results from
          \f[
          \rho_{lm}(r) \approx \rho(r0) + (r-r0) \frac{d\rho_{lm}}{dr}
          \f]
          The derivative is estimated as finite difference quotient of the
          first two points on the radial grid.

          The actual integrations are performed by advancing r and r'
          simultaneously, using Simpson integration.
          
        */
      static SxDiracVec<Double> 
      getHartreePotential (const SxDiracVec<Double> &rho);

      static SxDiracVec<Double> 
      getHartreeSmooth (const SxDiracVec<Double> &rho);

      /** \brief Compute screened Hartree energy
        @param rho1 (l,m) expansion coefficient of rho
        @param rho2 (l,m) expansion coefficient of rho'
        @param omega screening parameter
        @param dg    reciprocal space integration
        @return The screened Hartree energy
          This computes
          \f[
          \int d^3\mathbf r\, d^3\mathbf r' ~ \rho(\mathbf r)
           V^\textrm{SR}_\omega(|\mathbf r - \mathbf r'|) \rho(\mathbf r')
          \f]
          where
          \f[
          V^\textrm{SR}_\omega(|\mathbf r - \mathbf r'|)
          \f]
        */
      static double computeScreenedHartree (const SxDiracVec<Double> &rho1,
                                            const SxDiracVec<Double> &rho2,
                                            double omega, double dg);

      /// Exchange-correlation energy of last computeXC call
      double eXc;

      /// Grid for angular integrations
      SxSphereGrid::GridType gridType;

      /// Exchange-correlation functional
      SxXC::XCFunctional xcFunctional;

      /** Compute exchange-correlation energy & potential
          @param   rho density on radial mesh
          @param   pointer to eXcPtr exchange-correlation energy;
                   exchange-correlation energy will be added to this
          @return  xc potential

          The exchange-correlation energy is computed as
          \f[
          eXC = \int r^2 dr \rho_{00}(r) \varepsilon_{xc}(\rho_{00}(r))
              + \frac 12 \sum_{l\ne0,m} \rho_{lm}(r)^2 f_{xc}(\rho_{00}(r))
          \f]
          See Ref. 1, Eq. (29)

          \note At present, only LDA in Perdew-Zunger parametrization
                is available.

       */
      static
      SxArray<SxRadialMesh> computeXC2 (const SxArray<SxRadialMesh> &rho,
                                        double *eXcPtr);

      /** Compute exchange-correlation energy & potential
          @param   rho density on radial mesh
          @param   pointer to eXcPtr exchange-correlation energy;
                   exchange-correlation energy will be added to this
          @return  xc potential

          This routines uses an angular grid integration.

       */
      SxArray<SxRadialMesh> computeXC (const SxArray<SxRadialMesh> &rho);
      /** \brief Compute xc potential for spherical density
        @param rad          radial mesh (logarithmic grid)
        @param logDr        logarithmic grid increment
        @param rho          density (not its L=0 expansion coefficient!)
        @param xcFunctional xc functional (LDA or PBE)
        @param eXcPtr       if non-NULL, xc energy will be stored here
        @param simpsonWeights if true, vXc is the analytic derivative of
                              eXc including the weights from Simpson
                              integration.
        */
      static SxDiracVec<Double> 
      computeXC (const SxDiracVec<Double> &rad,
                 double logDr,
                 const SxDiracVec<Double> &rho,
                 SxXC::XCFunctional xcFunctional,
                 double *eXcPtr = NULL,
              SxXCFunctional::WhatToCompute xcMode = SxXCFunctional::ComputeXC,
                 bool simpsonWeights = false);

      SxRadialMesh computeXC (const SxRadialMesh &rho)
      {
         SxArray<SxRadialMesh> rhoSpin(1);
         rhoSpin(0)= rho;
         //vxcSpin = computeXC2(rhoSpin, &eXc);
         return computeXC(rhoSpin)(0);
      }

      //\name Auxiliaries
      //@{
      /// Integrate from 0 to r0
      static
      double integrate0 (const SxDiracVec<Double> &x);

      /** \brief Apply laplacian
          @param f     function
          @param r     radial vectors (logarithmic grid)
          @param logDr logarithmic grid step
          @param l     angular momentum
          @return      the Laplacian applied

          \f[
          \nabla^2 f(r) = \frac{1}{r^2} \frac{\partial}{\partial r}
                          [r^2 \frac{\partial}{\partial r} f(r)]
                        - \frac{l(l+1)}{r^2}
          \f]

          Note for \f$r_i = r_0 e^{\lambda i}\f$, we have
          \f[
          \frac{\partial f(r)}{\partial r} = 
          \frac{\partial f(r_i)}{\partial i}/(\lambda r_i)
          \f]

        */
      static SxDiracVec<Double> laplace (const SxDiracVec<Double> &f,
                                         const SxDiracVec<Double> &r,
                                         double logDr, int l);

      /** \brief Apply laplacian
          @param f     function
          @param l     angular momentum (or take from handle)
          @return      the Laplacian applied

          @note The radial basis is taken from the basis pointer.
          \f[
          \nabla^2 f(r) = \frac{1}{r^2} \frac{\partial}{\partial r}
                          [r^2 \frac{\partial}{\partial r} f(r)]
                        - \frac{l(l+1)}{r^2}
          \f]

          Note for \f$r_i = r_0 e^{\lambda i}\f$, we have
          \f[
          \frac{\partial f(r)}{\partial r} = 
          \frac{\partial f(r_i)}{\partial i}/(\lambda r_i)
          \f]

        */
      static
      SxDiracVec<Double> laplace (const SxDiracVec<Double> &f, int l = -1)
      {
         SX_CHECK (f.getBasisPtr ());
         const SxRadBasis &rad = f.getBasis<SxRadBasis> ();
         int is = f.handle->auxData.is;
         if (l == -1) l = f.handle->auxData.l;
         SX_CHECK (is < rad.getNSpecies (), is, rad.getNSpecies ());
         return laplace (f, rad.radFunc(is), rad.logDr(is), l);
      }

      //@}

};

inline double weightSimpson (int i, int n)
{
   SX_CHECK (i >= 0 && i < n, i, n);
   if (i == 0) return 1. / 3.;
   int d = n - i;
   if (d > 4) return ((i & 1) ? 4. : 2. ) / 3.;
   if (n & 1)  {
      if (d == 1) return 1./3.;
      if (d & 1) return 2./3.;
      return 4. / 3.;
   } // else
   if (d == 1) return 3./8.;
   if (d == 4) return 17./24.;
   return 9. / 8.;
}

#endif /* _SX_RADIAL_ATOM_H_ */
