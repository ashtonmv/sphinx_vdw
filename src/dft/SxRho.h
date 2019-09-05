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

#ifndef _SX_RHO_H_
#define _SX_RHO_H_

#include <SxBinIO.h>
#include <SxRBasis.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxAtomicStructure.h>
#include <SxPseudoPot.h>
#include <SxDensity.h>
#include <SxDFT.h>

/**
  \brief Electronic charge density \f$\varrho(R)\f$

  \b SxRho = S/PHI/nX Rho (\f$\varrho\f$)

  This class represents electronic charge densities in realspace.
  In the spin-compensated case this class contains 1 density, otherwise
  2 (\f$\varrho_\uparrow\f$ and \f$\varrho_\downarrow\f$.
  
  \ingroup group_dft
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxRho : public SxDensity
{
   public:

      RhoR rhoR;

      /** \brief standard constructor */
      SxRho ();

      void operator= (const SxRho &in);

      /** \brief initialization from (spin) density */
      SxRho (const RhoR &);

      /** \brief initialization from single density */
      SxRho (const SxMeshR &);

      /** \brief main constructor

          Initalize local variables, such as the pointer to the
          |R> basis as well as the number of electrons used for
          later normalization */
      SxRho (const SxRBasis &, int nSpin=1, Real8 nElectrons_=0.);
      /**
          \brief Constructor from file
          
          This constructs a density by reading it from a mesh file.
          Use this if you have no idea what kind of density you expect
          (nSpin, nElectrons unknown).
        */
      SxRho (SxBinIO &io, const SxRBasis* rBasisPtrIn = NULL);

      virtual ~SxRho ();

      /// Return number of spin channels
      inline int getNSpin () const
      {
         return int(rhoR.getSize ());
      }

      /** \brief extract a density beloning to a specified spin channel */
      SxMeshR &operator() (int iSpin);
      /** \brief extract a density beloning to a specified spin channel */
      const SxMeshR &operator() (int iSpin) const;
      /** \brief extract a density beloning to a specified spin channel */
      SxMeshR &operator() (SxAutoLoop &iSpin)
      {
         iSpin.setLimit (getNSpin ());
         return operator()(int(iSpin.i));
      }
      /** \brief extract a density beloning to a specified spin channel */
      const SxMeshR &operator() (SxAutoLoop &iSpin) const
      {
         iSpin.setLimit (getNSpin ());
         return operator()(int(iSpin.i));
      }
      /** \brief Returns the used |R> basis pointer */
      const SxBasis *getBasisPtr() const;
      const SxRBasis *getRBasisPtr() const;

      /** \brief compute atomic charge density

          This routine computes an atomic charge density according to
          \f[
             \varrho(G) 
                = 
             \sum_r
                \langle G | r \rangle
                \langle r | R_{i_s,i_a,n,l,m ,Y_{l=s,m}} \rangle
          \f] 
          and transforms the result to realspace via FFT. The radial
          function is computed as
          \f[
             R_{i_s,i_a,n,l,m}
             =
             f^{\mathrm{init}}_{i_s,i_a,n,l,m}
             \Psi^{\mathrm{ps} 2}_{i_s}
          \f]
          with \f$ \Psi^{\mathrm{ps}} \f$ being the pseudopotential
          wavefunctions given on a radial grid.
       */

//      RhoR &atomicChargeDensity (const SxAtomicStructure &, 
//                                 const SxPseudoPot &pot,
//                                 Real8 spinMoment=0.);
      // --- *** UGLY QUICK&DIRTY TO COPY FROM REL-1-0 ***
      //     *** TOBE REWRITEN IN A CLEAR WAY LATER *** {
      RhoR &atomicChargeDensity (const SxAtomicStructure &, 
                                 const SxPseudoPot &pot);
      // --- *** UGLY QUICK&DIRTY TO COPY FROM REL-1-0 ***
      //     *** TOBE REWRITEN IN A CLEAR WAY LATER *** }

      /** \brief Displace atoms via Hirshfeld decomposition
          @param toStr    target atomic structure
          @param atomRhoG ideal atomic shape in G space (for each species)
          @allAtomRptr    overlap of atomic reference densities; this will
                          be overwritten

          @note: this should usually be done to the deformation
                 density (after subtraction of reference atomic densities)
                 Therefore, the overlap of atomic densities is usually
                 available.
          @note: the original structure is taken from the G basis from
                 atomRhoG

          Hirshfeld decomposition is given by a weighting function
          \f[
               w_i(r) = \rho_i(r) / \sum_j \rho_j(r)
          \f]
          where \f$\rho_i(r)\f$ is the reference atomic density for atom i.

          This function displaces each atomic contribution
          \f[
              \Delta \rho_i(r) = w_i(r) \Delta\rho(r)
          \f]
          to its new position.

          Theoret. Chim. Acta (Berl.) 44, 129-138 (1977)
      */
      void displaceHirshfeld (const SxAtomicStructure &toStr,
                              const SxArray<PsiG> &atomRhoG,
                              SxDiracVec<Double> *allAtomRptr = NULL);
      
      /** \brief Compute charge density

          The charge density will be computed according to
          \f[
             \varrho(R) 
                = \sum_{i,\sigma,k} \omega_k f^{\mathrm{occ}}_{i,\sigma,k}
                  | \langle R | \Psi_{i,\sigma,k} \rangle |^2
          \f]
          with \f$\omega_k\f$ being the weight factors of the k points
          and \f$ f^{\mathrm{occ}}_{i,\sigma,k}\f$ the occupation numbers.

          \param  focc     occupation numbers
          \param  psiSet   wavefunction set \f$\Psi_{i,\sigma,k}\f$
          \todo should work with SxPsiSet instead */
      RhoR &computeRho (const Focc &focc, const SxPWSet &psiSet);
      /** \brief compute induced charge density

          The induced charge density can be computed according to
          \f[
             \varrho(R) 
                = \sum_{i,\sigma,k} \omega_k f^{\mathrm{occ}}_{i,\sigma,k}
                     \langle \Psi^{(1)} | \mathbf{R} \rangle
                     \langle \mathbf{R} | \Psi^{(0)} \rangle
                  + c.c.
                = 2 \sum_{i,\sigma,k} \omega_k f^{\mathrm{occ}}_{i,\sigma,k}
                     \Re{\langle \Psi^{(1)} | \Psi^{(0)} \rangle}
          \f]
          \param focc     Fermi occupation numbers
          \param psiSet0  \f$ \langle \mathbf{G} | \Psi^{(0)} \rangle \f$
          \param psiSet1  \f$ \langle \mathbf{G} | \Psi^{(1)} \rangle \f$ 
          \author Matthias Wahn, wahn@fhi-berlin.mpg.de
       */
      static RhoR computeRho (const Focc &focc, 
                              const SxPW &psiSet0, const SxPW &psiSet1);
      Real8 getNorm () const;
      void checkNorm ();
      void checkNorm (RhoR &) const;
      /** \brief Normalize charge density to the number of electrons.
          
           This routine normalizes the charge density that the integral
           contains the total number of electrons in the system, i.e.,
           \f[
              n_{\mathrm{electrons}} = \int_\Omega \varrho(R) d\Omega
           \f]
           with \f$\Omega\f$ being the unit cell. Numerically the integral
           is represented on the regular realspace grid as
           \f[
              n_{\mathrm{electrons}} = \sum_i \varrho(R_i) \Delta \Omega
           \f]
           where \f$\Delta \Omega = \frac{\Omega}{n_\mathrm{mesh-elements}}\f$.
           \sa SxAtomicStructure
           \sa SxRBasis
           \sa SxRBasis::dOmega
           \sa SxFFT3d
           */
      void normalizeRho ();
      void normalizeRho (RhoR *, Real8) const;
      /** \brief Initialize density with random numbers

          This routine initializes the charge denisty with random numbers.
          \f[
             \varrho_{\sigma}(R) = \mathrm{rand}()
          \f]
          The computed randomized density is normalized by an internal
          call of ::normalizeRho.

          \par Example:
\code          
   SxRBasis R (mesh, cell);
   SxRho rhoR (R, nSpin, nElectrons);
   rhoR.randomize ();
\endcode
      */
      void randomize ();
//      SxRho &mixRhoLinearly (const RhoR &, Real8 mixFactor=0.05);
//
      virtual void readRho  (const SxBinIO &);
      /// Synchronize across MPI tasks
      virtual void syncMPI ();

      using SxDensity::writeRho;
      using SxDensity::readRho;
      virtual void writeRho (SxBinIO &) const;
//
//   protected:

      /** \brief Pointer to the |R> basis */
      const SxRBasis *rBasisPtr;

      /** \brief the number of electrons in the system

          This local variable is required to check the norm of the
          charge density as well as the renomalization. It is set
          in the constructor */
      Real8 nElectrons;

      /// \name SxDensity interface
      ///@{

      /// Assign from a density
      virtual void operator= (const SxDensity &in);
      /// Add a density
      virtual void operator+= (const SxDensity &x);
      /// Subtract a density
      virtual void operator-= (const SxDensity &x);
      /// axpy-like operation
      virtual void plus_assign_ax (double a, const SxDensity &x);
      /// axpy-like operation for the spin density
      virtual void plus_assign_aspin (double a, const SxDensity &x);
      /// Scalar product
      virtual double operator| (const SxDensity &x) const;
      /// Square norm
      virtual double normSqr () const;
      ///@}

      /// Get change in density (as pointer)
      virtual SxDensity operator- (const SxDensity &x) const;
      /// Get a copy (as a pointer)
      virtual SxDensity getCopy () const;
      /// Get spin density (as a pointer)
      virtual SxDensity spin () const;
      /// Check if this is a spin-polarized density
      virtual bool hasSpin () const
      {
         return (getNSpin () == 2);
      }


      /// Renormalize (additive)
      virtual void renormalize ();

};

namespace Timer {
   enum RhoTimer {
      RhoMixing,
      Rho
   };
}

SX_REGISTER_TIMERS (Timer::RhoTimer)
{
   using namespace Timer;
   regTimer (RhoMixing,  "density mixing");
   regTimer (Rho,        "Charge density");

}
      

#endif /* _SX_RHO_H_ */
