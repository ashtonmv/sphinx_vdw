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

#ifndef _SX_PRECONDITIONER_H_
#define _SX_PRECONDITIONER_H_
#include <SxSymbolTable.h>
#include <SxDirac.h>
#include <SxTypes.h>
#include <SxDensity.h>
#include <SxDFT.h>
#include <SxDipoleCorrZ.h>

/**
 \brief Density preconditioner library
 \author C. Freysoldt, freysoldt@mpie.de

 \b SxPreconditioner = S/PHI/nX density preconditioners
 
 \par General
 It can be shown that the ideal preconditioner is given by the inverse
 dielectric function when the system is in the linear regime, i.e. the
 density is so close to the true ground state density \f$\rho_0\f$ that
 deviations in the effective potential are approximately proportional
 (yet nonlocal) to the deviations in the density:
 \f[
 v_{eff}[\rho + \delta \rho] = v_{eff}[\rho] 
                               + \hat V[\rho] \cdot \delta\rho
 \f]
 with a linear operator \f$\hat V\f$ which includes classical (Hartree)
 as well as quantum (exchange-correlation) fields.

 With a test density \f$\rho_{in} = \rho_0 + \delta\rho_i \f$ one
 obtains a potential \f$v = v_0 + \delta v_i\f$. This potential
 induces an extra density \f$\delta\rho_o\f$ which in the linear regime
 is given by the polarizability:
 \f[
 \delta \rho_o = \hat P \cdot \delta v_i 
               = \hat P \hat V \delta \cdot \rho_i
 \f]
 On the other hand, we can obtain the difference between the
 input and output deviations as the residuum R:
 \f[
 R = \delta \rho_o - \delta \rho_i 
   = (\hat P \hat V - 1) \delta \rho_i
 \f]
 We can obtain the groundstate density from 
 \f$\rho_{in} = \rho_0 + \delta\rho_i\f$ by
 \f[
 \rho_0 = \rho_{in} - \delta\rho_i
        = \rho_{in} + (1 - \hat P \hat V)^{-1} R
 \f]
 The term in the brackets is nothing else than the 
 dielectric function (well, its adjoint, but we only
 have Hermitean approximations where that doesn't matter), hence
 \f[
 \rho_o = \rho_{in} + \epsilon^{-1} R
 \f]
 which shows that the ideal preconditioner is the inverse
 dielectric matrix when the system is in the linear regime, and
 that the ground state density is obtained in a single (!) step.
 Of course, in practice one may not be close to the linear regime,
 and more importantly, the inverse dielectric function is difficult
 to obtain.

 The Kerker mixer (with a scaling of 1) corresponds to
 the inverse dielectric function of the homogenous electron gas
 (HEG) in the Thomas-Fermi approximation.
 The inverse dielectric function is given by
 \f[
 \epsilon^{-1}_{HEG}[\rho](G) = \frac{G^2}{G^2 + q_{TF}^2[\rho]}
 \f]
 with the Thomas-Fermi screening length \f$q_{TF}\f$ given by
 \f[
    q_{TF}^2[\rho] = \frac{4}{\pi} \sqrt[3]{3\pi^2 \rho}
 \f]
 For the Kerker mixer, $q_{TF}$ is treated as a free parameter
 (Kerker damping).

 A more complex dielectric function is given by the Lindhard
 function, which is the analytic solution for the HEG in the 
 random-phase approximation. It reads
 \f[
 \epsilon^{-1}_{HEG}[\rho](G) = \frac{G^2}{G^2 + q_{TF}^2[\rho] \cdot f}
 \f]
 with a correction factor
 \f[
 f = \frac12 + \frac{1-x^2}{4x} \ln \left| \frac{1+x}{1-x}\right|
 \qquad x = \frac{2G}{\pi}{q_{TF}^2[\rho]}
 \f]

 For semiconductors, we employ the model dielectric function of
 Cappellini, Del Sole, Reining, and Bechstedt
 \f[
 \varepsilon(G,\varepsilon_0,q_{TF})
 = 1 + \frac{1}{(1 + \frac{1}{\varepsilon_0 - 1}
                   + \alpha (G/q_{TF})^2
                   + \frac{G^4}{4 \omega_p^2}}
 \f]
 where \f$\omega_p\f$ is the plasma frequency
 \f[
 \omega_p^2 = \frac{\pi^2}{48} q_{TF}^6
 \f]
 
*/
class SX_EXPORT_DFT SxPreconditioner
{
   public:
      SxPreconditioner ();
      ~SxPreconditioner () { };

      /** \brief Initialize from symbol table */
      void readTable (const SxSymbolTable *);

      /// Possible preconditioner types
      enum Type {
         /// No preconditioner (only scaling)
         Identity = 0,
         /// Kerker preconditioner
         Kerker, 
         /// Preconditioner based on Lindhard dielectric function
         Lindhard, 
         /// Preconditioner based on CSRB dielectric function
         CSRB,
         /// Kerker-type preconditioner where rho(G) is small
         Freysoldt,
         /// Elliptic preconditioner
         Elliptic,
         /// Number of types
         nTypes // This must always be last
      };

      
      /// Apply preconditioner
      SxMeshR operator* (const SxMeshR &) const; 

      /// Apply preconditioner
      SxDensity operator* (const SxDensity &) const; 
      
      /// Update preconditioner with latest density
      void setRho (const SxMeshR &);

      /// Update preconditioner with latest density
      void setRho (const SxDensity &);

      /// Print preconditioner info
      void print () const;

      /** \name Preconditioner parameters */
      //@{
      /// Global scaling factor
      double scaling;

      /// Spin scaling factor
      double spinScaling;

      /// Type of screening model
      enum Type type;
         
      /// Kerker damping (for Kerker)
      double kerkerDamping;

      /// Dielectric constant (for CSRB)
      double dielecConstant;
      
      ///@}

      void useKerker (double scalingIn, double dampingIn);

   protected:
      /** \name Internal parameters */
      //@{
      /// Reciprocal space preconditioner
      SxDiracVec<TReal8> K;
      ///@}

      /** \name Internal functions */
      ///@{
       
      /// Setup preconditioner elements in G space
      void setupK (const SxGBasis &G, double q2tf);

      /// Elliptic parameters
      SxDiracVec<Double> ellipticA, ellipticB;

      /// Dipole correction for elliptic preconditioner
      SxPtr<SxDipoleCorrZ> dipoleCorr;

      /// Elliptic solver
      SxMeshR solveElliptic (const SxMeshR &meshR) const;

      /// Elliptic operator (in G space)
      PsiG applyElliptic (const PsiG &resG) const;

   public:
      /** \brief Inverse of the Cappellini/Del Sole/Reining/Bechstedt 
                 model dielectric function
                 
        ref 4: G. Cappellini, R. Del Sole, L. Reining, F. Bechstedt,
               "Model dielectric function for semiconductors"
               Phys. Rev. B 47, 9892 (1993)
      */
      static double invCSRB (double g2, double q2tf, double eps0);
      
      /** \brief Inverse of the Lindhard model dielectric function */
      static double invLindhard(double g2, double q2tf);
      ///@}
   public:
      /// Get the type of the preconditioner
      Type getType () const { return type; }
};

#endif /* _SX_PRECONDITIONER_H_ */
