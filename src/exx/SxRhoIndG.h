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

#ifndef _SX_RHO_IND_G_H_
#define _SX_RHO_IND_G_H_

#include <SxExx.h>
#include <SxDirac.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxRBasis.h>
#include <SxGBasis.h>
#include <SxGkBasis.h>
#include <SxFockGk.h>

/** \brief Change of charge density induced by a small change of the potential

    \b SxClass = S/PHI/nX induced charge density

    Computes the change in the charge density "induced" by the exchange (Fock)
    operator \f$ V_{\rm x} \f$. In G space it reads
    \f[
      \rho_{\rm ind}^{\rm stEXX}({\mathbf G})
        := E({\mathbf G}) + E^{\ast}(-{\mathbf G})\;,
    \f]
    with
    \f[
      E({\mathbf G}) := \frac{1}{\Omega_{\rm cell}}
        \sum_{v, c, {\mathbf k}}
        \frac{\langle v{\mathbf k} | V_{\rm x} | c{\mathbf k} \rangle
              \langle c{\mathbf k} |
              {\rm e}^{\rm i}{\mathbf G}{\mathbf r} |
              v{\mathbf k} \rangle}
             {\varepsilon_{v{\mathbf k}} - \varepsilon_{c{\mathbf k}}}
      \;.
    \f]
    The prefactor of 2, stated in Eq. (17) of
    <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 10031 (1999) </a> is wrong and rather belongs to Eq. (8) of the same
    reference (where it is missing), expressing there, that one deals
    with non-spin-polarised systems only.

    \f$v\f$ and \f$c\f$ in the above equation refer to the valence and
    conduction bands, respectively, of a semiconductor or insulator.

    For detailed information see <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 10031 (1999) </a>.

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de */
class SX_EXPORT_EXX SxRhoIndG
{
   public:

      /** constructor */
      SxRhoIndG ();

      /** constructor */
      SxRhoIndG (const SxGBasis &G, const SxGkBasis &gk, int nG_);

      /** destructor */
      ~SxRhoIndG ()  { /* empty */ }

      /** computes the induced charge density \f$ \rho_{\mathrm ind} \f$ in
          reciprocal space

          \TODO Check that G and -G component succeed each other!
       */
      void compute (const SxPW &waves, const SxFermi &fermi, SxFockGk &FockOp);

      /** get the induced charge density in real space */
      SxMeshR getRhoIndR () const;

      /** get the induced charge density in reciprocal space */
      SxMeshG getRhoIndG () const;

      /** yields the contribution of the exchange (Fock) operator to the total
          energy of the system */
      PrecEnergy getXEnergy () const;

      ///@name Debugging
      //@{
      /** \brief test FFT prefactor

        When porting the code to other systems, there might occur
        some inconsistency with the FFT prefactors. This function is
        meant to have an easy check on this.

        The computation of \f$ \rho_{\mathrm ind} \f$ requires expressions
        of the form
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | c\mathbf{k} \rangle\; ,
        \f]
        where \f$ |v\mathbf{k}\rangle \f$ and \f$ |v\mathbf{k}\rangle \f$
        shall denote valence and conduction band states, respectively,
        and the operator is given in its real space form.

        This is conveniently performed using the Fourier trafo:
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | c\mathbf{k} \rangle
           = \int\limits_{\Omega} \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                   \mathbf{r}} \, \psi_{v\mathbf{k}}^{\ast}(\mathbf{r})\,
                   \psi_{c\mathbf{k}}(\mathbf{r})\, \mathrm{d}\mathbf{r}\; ,
        \f]
        with \f$\Omega\f$ being the volume of the unit cell.

        This test now just performs
        \f[ \langle v\mathbf{k} | \mathrm{e}^{-\mathrm{i}\mathbf{G}\cdot
                    \mathbf{r}} | v\mathbf{k} \rangle
        \f]
        for the lowest valence band at one \f$\mathbf{k}\f$-point. Note,
        that also the ket is a valence band here. If all prefactors are
        correct, the \f$(\mathbf{G}=\mathbf{0})\f$-component must be one,
        due to normalization.
       */
      void checkFFTfactor (const SxPW &waves);
      //@}

   protected:

      /** the induced charge density in reciprocal space */
      SxMeshG           rhoIndG;

      /** pointer to basis in reciprocal space */
      const SxGBasis   *GPtr;

      /** pointer to the |G+k> basises */
      const SxGkBasis  *gkPtr;

      /** real space basis */
      SxRBasis          R;

      /** number of G components used in \f$ \rho_{\mathrm ind} \f$ */
      int               nG;

      /** number of k-points */
      const int         nk;

      /** number of valence bands */
      int               nv;

      /** number of conduction bands */
      int               nc;

      /** FFT prefactor used for the transformation of a function with the
          dimension of a charge density

          \note Check this factor using SxRhoInd::checkFFTfactor.
       */
      SX_FFT_REAL       fftFactorRtoG;

      /** the contribution of the exchange (Fock) operator to the total energy
          of the system */
      PrecEnergy        eX;
};

#endif /* _SX_RHO_IND_G_H_ */
