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

#ifndef _SX_RHO_MIXER_H_
#define _SX_RHO_MIXER_H_

#include <SxDensity.h>
#include <SxPreconditioner.h>
#include <SxRhoFilter.h>
#include <SxXC.h>
#include <SxDFT.h>

/**  brief Charge densitiy mixer

    \b SxRhoMixing = SPHInX Charge density (\f$\rho\f$) mixing class

    \ref Kresse96a and \ref Kresse96b.

    ....

    \author John Doe, johndoe@mpie.de */
class SX_EXPORT_DFT SxRhoMixer
{
   public:

      enum MixerType       { Linear, Pulay };
      enum NormHandling    { RenormOn, RenormOff };

      /**
        \param  rhoMixing   mixing parameter for linear mixer
        */
      SxRhoMixer (MixerType type=Linear, 
                  double rhoMixingIn = 0.1, 
                  int maxStepsIn=1);
      /** \brief Initialize by reading parameters from symbol table
          \param table the symbol table
          \param spin If false, the spin mixing is set to -1.
        */
      SxRhoMixer (const SxSymbolTable *table, bool spin = true);

      ~SxRhoMixer ();

      /** \brief Read parameters from symbol table
          \param table the symbol table
          \param G G-basis for preconditioner
       */
      void readTable (const SxSymbolTable *cmd);

      void setType (MixerType);

      void setMaxSteps (int);
      int getMaxSteps ();

      /** \brief set modus to handle the renormalization

          In case one does not want the output of the mixing routine to be
          renormalized to a certain value, the renormalization procedure
          (default) can get switched off. This is, e.g., necessary if the
          output is not a charge density, but a potential, as happens in the
          gsEXX scheme.
          \param modus  denotes the modus of NormHandling. Expects the values
                        'RenormOn' (default) or 'RenormOff'. */
      void setNormModus (NormHandling modus=RenormOn);

      void addRhoIn  (const SxDensity &);
      void addRhoIn  (const SxVector<Double> &);
      void addRhoOut (const SxDensity &);
      void addRhoOut (const SxVector<Double> &);
      /** \brief Append residue to mixer's history

          This is an alternative way of setting up the mixer's history. 
          Instead of providing previous \f$\varrho^{\mathrm{in}}\f$ and 
          \f$\varrho^{\mathrm{out}}\f$ it is possible to provide
          \f$\varrho^{\mathrm{in}}\f$ and \b R.
          Don't mix up both schemes!!!
          \sa addRhoIn
          \sa addRhoOut */
      void addResidue (const SxDensity &resIn);

      /// \brief Computes the mixed (optimal) density
      SxDensity getMixedRho ();
      SxDensity getMixedRhoLinear ();
      SxDensity getMixedRhoPulay ();

      SxDensity getRhoIn (int i) const;
      double getNormR () const    { return normR; }
      double getNormRPol () const { return normRPol; }

      /** The preconditioned Pulay mixer may produce negative densities
          as well as ill-shaped densities in the vacuum region of a slab
          system. This creates massive convergence problems when the
          dipole correction is used, as this requires a well defined
          density minimum at low densities (otherwise the dipole and
          hence the potential fluctuates wildly).

          In order to overcome this problem, the density mixer can be
          forced to use linear mixing when the "output" density is small.
          This behaviour is switched on with this flag. It is set by
          SxHamSolver.

          \brief Enable linear mixing in low-density regions
          \sa maxVacDensity, vacSmoothFactor
        */
      bool linearVacMixing;
      /** In order to smoothly switch on the linear vacuum mixing, the
          following weighting function is used:
          \f[
          \sigma = \frac{1}{1 + (\rho/\rho_0)^{-1/\beta}}
          \f]

          This is a Fermi-like function for the logarithm of the density.
          
          \f$\rho_0\f$ is set to the root mean square residue. However,
          in early stages this may be too large (effectively using linear
          scaling everywhere), so \f$\rho_0\f$ is reduced to maxVacDensity
          when the residue is large.

          \note
          The scaling exponent \f$\beta\f$ is determined from vacSmoothFactor.
          If the scaling is less than 0.01 (1%) from 0 or 1, the border values
          are used, i.e. 0.009 -> 0 or 0.991 -> 1.
          \note
          Currently this value is set to a ill-tested value. If necessary,
          one could add this as a parameter to the input file.

          \sa vacSmoothFactor
          \brief Maximum vacuum density for linear mixing in the vacuum.
        */
      double maxVacDensity;
      /** This defines the range of the scaling for the linear mixing
          in the vacuum. The logarithm of this value plays the role
          of the Fermi energy in the Fermi-like scaling function for
          the logarithm of the density.
          \sa maxVacDensity
          \note This is set to 1.5, which is not well tested.
          \note With
          \f[
              f = vacSmoothFactor ^{-ln(0.01)}
          \]
          the border scaling 0.01 is reached for rho0 * f, and
          0.99 for rho0 / f.
          \todo Find a more intuitive parameter. Maybe beta, or 1/beta, or f.

          \brief Smoothening parameter for linear mixing in the vacuum
        */
      double vacSmoothFactor;
      /// Whether to use adaptive scaling in Pulay mixer
      bool adaptiveScaling;
      /// Printout function
      void print () const;

      /// Request residue profiles (for debugging)
      bool residueProfile;

   protected:

      /// \brief History buffer of previous input densities
      SxList<SxDensity> rhoIn;  // :it

      /// \brief History buffer of previous output densities
      SxList<SxDensity> rhoOut;  // :it

      /// \brief History buffer of previous residues
      SxList<double> res2;

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

   public:
      /// \brief Mixing parameter for the linear mixer
      double rhoMixing;

      /// \brief Mixing parameter for spin mixer
      double spinMixing;

      /// \brief Residual vector norm \f$ |R|= \sqrt{\langle R|R \rangle } \f$
      double normR;
      double normRPol;

   protected:
      const SxXC::Type UP, DN;

      void removeFirst ();

      /// Reset Pulay mixer (remove history except for one)
      void reset ();

   public:
      /// Preconditioner
      SxPreconditioner preconditioner;

      /// Filter
      SxPtr<SxRhoFilter> filterPtr;

};

#endif /* _SX_RHO_MIXER_H_ */
