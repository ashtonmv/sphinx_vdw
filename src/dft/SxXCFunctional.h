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

#ifndef _SX_XC_FUNCTIONAL_H_
#define _SX_XC_FUNCTIONAL_H_

#include <SxArray.h>
#include <SxConstants.h>
#include <SxDFT.h>
#include <SxPBExHole.h>

/** \brief   Container of Exchange-correlation functionals.
  Here the following exchange-correlation functionals are defined:
       - LDA functionals:
          -# Perdew/Zunger '81
             <a href="http://prola.aps.org/abstract/PRB/v23/i10/p5048_1">
             PRB 23,5048(1981)</a>
       - GGA functionals:
          -# Perdew/Wang '91 (prepared)
             <a href="http://prola.aps.org/abstract/PRB/v45/p13244_1">
             PRB 45, 13244 (1992)</a>.
          -# Perdew/Burke/Enzerhof
             <a href="http://prola.aps.org/abstract/PRL/v77/p3865_1">
             PRL 77, 3865 (1996)</a>,
             <a href="http://prola.aps.org/abstract/PRB/v54/p16533_1">
             PRB 54, 16533 (1996)</a>
*/
class SX_EXPORT_DFT SxXCFunctional {
   public:
      enum WhatToCompute {
         /// Compute xc energy
         ComputeEnergy    = 0x0001,
         /// Compute potential
         ComputePotential = 0x0002,
         /// Compute kernel
         ComputeKernel    = 0x0004,
         /// Compute energy and potential
         ComputeEandV     = 0x0003,
         /// Compute exchange part
         ComputeX         = 0x0010,
         /// Compute correlation part
         ComputeC         = 0x0020,
         /// Compute exchange & correlation
         ComputeXC        = 0x0030,
         /// Compute default
         ComputeExcVxc    = 0x0033
      };

   protected:
      /// Compute mode (flags from WhatToCompute)
      int mode;

      /// Number of spin channels
      int nSpin;

      /// Pointer to PBE exchange hole class
      SxPtr<SxPBExHole> xHolePtr;

   public:
      /// Exchange-correlation energy density
      double eXc;

      /// Exchange-correlation energy density per electron
      double epsXc;

      /// Exchange-correlation potential (iSpin)
      SxArray<double> vXc;

      /// GGA only: derivative of eXc wrt \f$|\nabla \varrho_\sigma|\f$
      SxArray<double> eXc_grad;

      /// GGA only: derivative of eXc wrt \f$|\nabla \varrho|\f$
      double eXc_gradRhoTl;

      /// Exchange-correlation kernel (iSpin:jSpin)
      SxArray<SxArray<double> > fXc;

      /// HSE screening constant
      static double omegaHSE;

   protected:
      /// for PBE/LSDA_PW: derivative of correlation energy wrt rs
      double ec_rs, 
      /// for PBE/LSDA_PW: derivative of correlation energy wrt zeta
             ec_z;
      /// rs for total density
      double rsTot;
   public:
      /// Change number of spins
      void resize (int nSpinIn);

      /// Initialize
      void init (int nSpinIn, int modeIn)  {
         mode = modeIn;
         resize (nSpinIn);
      }

      /// Get number of spins
      int getNSpin () const { return nSpin; }

      /// Constructor
      SxXCFunctional (int nSpinIn = 1, int modeIn = ComputeExcVxc)
      {
         init (nSpinIn, modeIn);
      }

      /** \brief Compute LDA (Perdew/Zunger '81)
        \note No matter if ComputeEnergy is set, the xc energy is
              always computed.
        */
      void computeLDA (double rho);

      /** This defines the LDA functional according to
          <a href="http://prola.aps.org/abstract/PRB/v23/i10/p5048_1">
          PRB 23,5048(1981)</a>. The exchange-correlation energy reads
          \f[
            \epsilon_xc(r_s,\zeta)
               = (a_x^U + a_x^P) \frac{1}{r_s}
               + \epsilon^U(r_s)
               + f(\zeta)[\epsilon_c^P(r_s) - \epsilon_c^U(r_s)]
          \f] with f being
          \f[
            f(\zeta) = \frac{(1+\zeta)^{4/3} + (1-\zeta)^{4/3}-2}
                            {2^{4/3} - 2}.
          \f]
          @brief The spin-polarized LDA functional. */
      void computeLSDA (double rhoUp, double rhoDown);

      /** \brief Perdew-Wang parametrization of the LSDA functional
       */
      void computeLSDA_PW (double rhoUp, double rhoDown);

      void computePBE (double rhoUp, double gradRhoNormUp,
                       double rhoDown, double gradRhoNormDown,
                       double gradNorm);

      /** \brief Compute short-ranged PBE-exchange
          @param rho      density for one spin channel
          @param gradRho  density for this spin channel
          @param iSpin    the spin channel
          @param scaling  scaling factor

          \note The scaled contribution is added to eXc/vXc.

          @note The screening length is defined by omegaHSE
        */
      void addScreenedExPBE (double rho, double gradRho, int iSpin,
                                 double scaling = 1.);

      /** Compute local HSE functional

         This computes
         \f[
         \frac 34 E_x^{SR-PBE} + E_x^{LR-PBE} + E_c^{PBE}
         = E_{xc}^{PBE} - \frac 14 E_{xc}^{SR-PBE}
         \f]
        */
      void computeHSE (double rhoUp, double gradRhoNormUp,
                       double rhoDown, double gradRhoNormDown,
                       double gradNorm);

      /** Compute local PBE hybrid functional

         This computes
         \f[
         (1-x) E_x^{PBE} + E_c^{PBE}
         \f]
        */
      void computePBEHybrid (double rhoUp, double gradRhoNormUp,
                             double rhoDown, double gradRhoNormDown,
                             double gradNorm, double alphaExact);


      ///@{
      /// @name Auxiliaries
      /// Compute Wigner-Seitz radius from density\f$r_s\f$
      static inline double getRs (double rho)
      {
         return rho > 1e-50 ? cbrt (3. / (FOUR_PI * rho))
                            : 2.8794119e+16 /* value for 1e-50 */;
      }

      /// Compute spin polarization
      static inline double getZeta (double rhoUp, double rhoDown)
      {
         double rhoTot = fabs(rhoUp) + fabs(rhoDown);
         return (rhoTot > 0.) ? (fabs(rhoUp) - fabs(rhoDown)) / rhoTot : 0.;
      }

      static inline double getKf (double rS)  {
         static const double kfFac = pow  (9./4. * PI, 1./3.);
         return kfFac / rS;
      }
      ///@}
      
      /// Test LDA implementation
      void test (double rhoMin, double rhoMax, double dRho);

      //---------------------------------------------------------------------
      /**@name Perdew/Wang '91 */
      //@{
      /** Computes according to 
        <a href="http://prola.aps.org/abstract/PRB/v45/p13244_1">
        PRB 45, 13244 (1992), [REF]</a> (Eqn. A5). 
        \f[
           \frac{\delta G}{\delta r_s}
              = -2 A \alpha_1 \ln (1+\frac{1}{Q_1})
              - \frac{Q_0 Q_1'}{Q_1^2+Q_1}
        \f]
         @note  This routine is used by PBE (SxXC::correlationPBE) as well.
         @param A   A, as defined in REF, table 1
         @param a1  \f$ \alpha_1 \f$ as in REF, table 1
         @param b1  \f$ \beta_1  \f$ as in REF, table 1
         @param b2  \f$ \beta_2  \f$ as in REF, table 1
         @param b3  \f$ \beta_3  \f$ as in REF, table 1
         @param b4  \f$ \beta_4  \f$ as in REF, table 1
         @param rs  Wigner-Seitz radius
         @param grs Returned value: \f$ \frac{\delta G}{\delta r_s} \f$
         @return    \f$ g = Q_0 Q_1' \f$
       */
      static double G_PW91 (double A, double a1, double b1, double b2, 
                            double b3, double b4, double rs, double *grs);
      void testPW91 ();
      //@}

      /// The hybrid functional mixing factor
      static double alphaHybrid;
};


#endif /* _SX_XC_FUNCTIONAL_H_ */
