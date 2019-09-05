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

#ifndef _SX_RHO_MIXER_G_H_
#define _SX_RHO_MIXER_G_H_

#include <SxExx.h>
#include <SxRho.h>
#include <SxXC.h>

/**  brief Charge densitiy mixer

    \b SxRhoMixing = S/PHI/nX Charge density (\f$\rho\f$) mixing class

    \ref Kresse96a and \ref Kresse96b.

    ....

    \author John Doe, johndoe@sfhingx.de */
class SX_EXPORT_EXX SxRhoMixerG
{
   public:

     enum MixerType       { Linear, Pulay };
     enum NormHandling    { RenormOn, RenormOff };
     
     class SxPreconditioner
     {
        public:
           SxPreconditioner () { };
           virtual ~SxPreconditioner () { };
           
           virtual SxMeshG operator* (const SxMeshG &) const { 
              SX_EXIT; 
              return SxMeshG();
           }
     };

     class SxIdentity : public SxPreconditioner
     {
        public:
           double a;
           SxIdentity () : SxPreconditioner () { a = 1.; }
           virtual ~SxIdentity () { }
           virtual SxMeshG operator* (const SxMeshG &m) const { return a*m; }
     };
     SxIdentity identity;
     
     class SxKerker : public SxPreconditioner
     {
        public:
           SxKerker (const SxGBasis &, double a=0.8, double q=1.5);
           virtual ~SxKerker () { }
           
           virtual SxMeshG operator* (const SxMeshG &) const;
           
           SxDiracVec<TReal8> K;
     };

      /**
        \param  rhoMixing   mixing parameter for linear mixer
        */
      SxRhoMixerG (MixerType type=Linear, 
                   double rhoMixingIn = 0.1, 
                   int maxStepsIn=1,
                   double spinMix=1.);
      ~SxRhoMixerG ();

      void setType (MixerType);
      void setPreconditioner (const SxPreconditioner &);
      void setMaxSteps (int);
      /** \brief set parameters for density normalization

          The mixing scheme contains numerical operations which lead to
          small diviations in the charge density norm. In order to regain
          stability in the algorithm it is necessary to renomalize the
          output density. So the routine ::getMixedRho returns an explicitly
          normalized density. This routine sets up the required normalization
          parameters
          \param rhoNormIn  norm of the density, usually number of electons
          \param dOmegaIn   infinitismal volume of the unit cell,
                            isually volume / FFT mesh size. */
      void setNorm (double rhoNormIn, double dOmegaIn);

      /** \brief set modus to handle the renormalization

          In case one does not want the output of the mixing routine to be
          renormalized to a certain value, the renormalization procedure
          (default) can be disabled. This is, e.g., necessary when the output
          is not a charge density, but a potential, as happens in the gsEXX
          scheme.
          \param modus  denotes the modus of NormHandling. Expects the values
                        'RenormOn' (default) or 'RenormOff'. */
      void setNormModus (NormHandling modus=RenormOn);

      void addRhoIn  (const RhoR &);
      void addRhoInG (const RhoG &);
      void addRhoIn  (const SxVector<Double> &);
      void addRhoOut (const RhoR &);
      void addRhoOut (const SxVector<Double> &);
      /** \brief Append residue to mixer's history

          This is an alternative way of setting up the mixer's history. Instead of
          providing previous \f$\varrho^{\mathrm{in}}\f$ and \f$\varrho^{\mathrm{out}}\f$
          it is possible to provide \f$\varrho^{\mathrm{in}}\f$ and \b R.
          Don't mix up both schemes!!!
          \sa addRhoIn
          \sa addRhoOut */
      void addResidue (const RhoR &);
      void addResidueG (const RhoG &);

      /// \brief Computes the mixed (optimal) density
      RhoR getMixedRho () const;
      RhoG getMixedRhoG () const;
      RhoG getMixedRhoLinear () const;
      RhoG getMixedRhoPulay () const;

      RhoG getRhoIn (int i) const;
      double getNormR () const  { return normR; }
      double getNormRPol () const  { return normRPol; }

   protected:

      /// \brief History buffer of previous input densities
      mutable SxList<RhoG> rhoIn;  // :it,:iSpin

      /// \brief History buffer of previous output densities
      mutable SxList<RhoG> rhoOut;  // :it

      /// \brief Residual vector
      mutable SxList<RhoG> R;    // :it

      /// \brief change of the residual vector
      mutable SxList<RhoG> dR;   // :it

      /// \brief Change of the input densities
      mutable SxList<RhoG> dRhoIn; // :it

      /** \brief contains the maximum number of densities which can
                 be held in the history buffers

          The default is 1 (linear mixer) */
      int maxSteps;


      /** \brief Mixer type

          The default is a linear mixer */
       MixerType type;

      /** \brief Handling of the normalization

          'RenormOn' renormalizes the results afterwards to the value set by
          setNorm, 'RenormOff' does not. */
      NormHandling renormModus;

      /// \brief Mixing parameter for the linear mixer
      double rhoMixing;

      /// \brief Mixing parameter for spin mixer
      double spinMixing;

      /// \brief Residual vector norm \f$ |R|= \sqrt{\langle R|R \rangle } \f$
      mutable double normR;
      mutable double normRPol;

      /// \brief infinitesimal volume of the super cell
      double dOmega;

      /// \brief norm of the density, used for renomalization
      double rhoNorm;

      const SxPreconditioner *precondPtr;

      const SxXC::Type UP, DN;


      void removeFirst () const;
};

#endif /* _SX_RHO_MIXER_G_H_ */
