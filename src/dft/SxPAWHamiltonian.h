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

#ifndef _SX_PAW_HAMILTONIAN_H_
#define _SX_PAW_HAMILTONIAN_H_

/*
References:
 1. P. E. Bloechl, Phys. Rev. B 50, 17953 (1994).
 2. C. Freysoldt, PAW implementation notes
*/


#include <SxDFT.h>
#include <SxHamiltonian.h>
#include <SxPAWPot.h>
#include <SxRadialMesh.h>
#include <SxRadMat.h>
#include <SxXC.h>
#include <SxPAWRho.h>
#include <SxAOBasis.h>
#include <SxPartialWaveBasis.h>
#include <SxDipoleCorrZ.h>
#include <SxPAWExchange.h>
#include <SxHubbardU.h>

/**\brief PAW Hamiltonian

  \b SxPAWHamiltonian = S/PHI/nX Projector Augmented Waves Hamilton operator


  \ingroup group_dft
  \author Christoph Freysoldt, freysoldt@mpie.de and Bjoern Lange, lange@mpie.de

  */
class SX_EXPORT_DFT SxPAWHamiltonian : public SxHamiltonian
{
   public:

      /// Empty constructor
      SxPAWHamiltonian ();
      /** \brief Constructor/initializer from density in file
          @param G         global G basis (with registered R and structure)
          @param potPtrIn  PAW potentials
          @param Gk        global G+k basis
          @param rhoFile   filename of sxb file with PAW density
          @param table     symbol table (for xc etc.)

          @note Does not yet work for hybrid functionals.

        */
      SxPAWHamiltonian (const SxGBasis  &G,
                        const SxPtr<SxPAWPot> &potPtrIn,
                        const SxGkBasis &Gk,
                        const SxString  &rhoFile,
                        const SxSymbolTable *table);
      virtual ~SxPAWHamiltonian ();

      virtual void compute (bool tauChanged=true, bool rhoChanged=true);

      virtual void printEnergies () const;
      virtual void read (const SxSymbolTable *);
      /** \brief Update internal settings from some input group
          @param table The input group to read from
          @return true if anything changed

          At present, this routine supports
          - switching off the dipole correction
          - changing the xc mesh

       */
      virtual bool rereadTable (const SxSymbolTable *table);
      /** \brief Set internal settings to default
          @param table some group in the symbol tree
       */
      virtual void backToDefault (const SxSymbolTable *table);

      /** \brief Setup internal variables relevant for HubbardU
          @param gk The |G+k> basis
          @param nSpin number of spin channels

          Customized-orbital HubbardU requires additional atomic-orbital
          projectors, density matrices and associated Hamiltonians. This
          routine sets them up. It must be called after read, and only if
          HubbardU is used.
        */
      void setupForHubbard (const SxGkBasis &gk, int nSpin);


   public:

      static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);

      /// PAW Hamiltonian components to compute
      enum Contrib {
         /// Nuclear Coulomb potential
         CalcNuc        = 0x00000001,
         /// Local vBar potential
         CalcBar        = 0x00000002,
         /// Exchange-correlation on radial grid
         CalcXcRad      = 0x00000010,
         /// Exchange-correlation in PS part
         CalcXcPS       = 0x00000020,
         /// Exchange-correlation
         CalcXc         = 0x00000030,
         /// 1-center kinetic energy correction
         CalcKinRad     = 0x00000040,
         /// Pseudo-kinetic energy
         CalcKinPS      = 0x00000080,
         /// Kinetic energy
         CalcKin        = 0x000000c0,
         /// Compensation charge density
         CalcNHatRad    = 0x00000100,
         /// Compensation charge density
         CalcNHatGHard  = 0x00000200,
         /// Compensation charge density
         CalcNHatGSoft  = 0x00000400,
         /// Compensation charge density
         CalcNHatG      = 0x00000600,
         /// Compensation charge density
         CalcNHat       = 0x00000700,
         /// Core energy
         CalcCore       = 0x00000800,
         /// Hartree potential in radial grid
         CalcHartreeRad = 0x00001000,
         /// Hartree potential from rhoPS
         CalcHartreePS  = 0x00002000,
         /// Dipole correction
         CalcDipoleC    = 0x00004000,
         /// Hard-soft compensation charge correction
         CalcU          = 0x00008000,
         /// Hartree potential
         CalcHartree    = 0x0000ff00,
         /// Local PS potential contributions (and vBar)
         CalcVPS        = 0x00006022,
         /// Core-valence exact exchange
         CalcCorValExch = 0x00010000,
         /// Valence exact exchange
         CalcValExch    = 0x00020000,
         /// Hubbard U
         CalcHubbardU   = 0x00100000,
         /// Calculate everything
         CalcAll        = 0xffffffff
      };

      /// Hamiltonian contributions (set of flags)
      long int hContrib;

      /// Dipole correction
      SxPtr<SxDipoleCorrZ> dipoleCorr;

      /// Hubbard U
      SxPtr<SxHubbardU> hubbardU;

      /// External potential
      SxMeshR vExt;

      /** \brief Symmetry group for real Ylm (:iSym,:l,(l+m x l+m) ) */
      SxYlmRotGroup ylmRot;

      /** \brief Set up projectors in reciprocal space
        */
      void setupProj (const SxGkBasis &gk);

      /// Radial multipole moments
      SxArray<SxArray<SxVector<Double> > > Qrl;

      /// Total energy
      double eTotal;
      /// exchange-correlation energy
      double eXc;
      /// Hartree energy
      double eHartree;
      /// \f$\overline v\f$ energy
      double eBar;
      /// External potential energy
      double eExt;
      /// Double counting energy
      double eDoubleCounting;
      /// Core energy
      double eCore;

      /// pseudoized potential (Ref. 1, Eq. 34)
      SxArray<SxMeshR> vPS;

      /** \brief Compute density-related terms in Hamiltonian
        */
      void computeRhoTerms ();

      /// 1-center Hamiltonian (potentials) (:iSpin)(:is,:ia)(:ip,:jp)
      SxRadMat Vij;
      /// Atomic spin bias for spin constraints
      SxVector<Double> nuA;

      /// soft PAW radius
      static double rSoft;

      /// Compute eHartreeU
      double computeU (SxArray<SxArray<SxVector<Double> > > *dEdQrl);

      /// Pointer to PAW potential container
      SxPtr<SxPAWPot> potPtr;

      /// Pointer to projectors
      SxPtr<SxAOBasis> projBasis;

      /// Pointer to partial-wave basis
      SxPtr<SxPartialWaveBasis> pBasis;

      /** Pointer to the xc functionals */
      SxPtr<SxXC> xcPtr;

      /// Pointer to exact exchange operator
      SxPtr<SxPAWExchange> exchangePtr;

      // TODO: SxConstPtr<SxRBasis>
      SxRBasis *rPtr;

      /// \f$ \overline v\f$ in G-space (:is)
      SxArray<PsiG> vBarG;

      /// rhoCorePS in G-space (:is)
      SxArray<PsiG> rhoCorePSG;

      /// Set up \f$ \overline v\f$ and rhoCore in G-space
      void setupPSRefG (const SxGBasis &gBasis);

      /// Compute 1-center density matrix
      SxRadMat computeDij (const SxPWSet &waves,
                                    const Focc &focc) const;

      /// Compute 1-center density matrix
      static
      SxRadMat computeDij (const SxPWSet &waves,
                           const Focc &focc,
                           const SxConstPtr<SxPAWPot> &potPtr,
                           const SxPartialWaveBasis &pBasis,
                           const SxYlmRotGroup &ylmRot);

      /// Compute energy & Hamiltonian
      void compute (const SxPWSet &waves, const SxFermi &fermi,
                    bool computeNewRho = true);

      /// Apply Hamiltonian
      virtual PsiG apply (const PsiG &psi) const;

      /// Return last computed energy
      virtual PrecEnergy getEnergy () const { return eTotal; }

      virtual PrecEnergy getDoubleCounting () const { return eDoubleCounting; }

      /** \brief Compute new energy*/
      virtual PrecEnergy getEnergy (const SxPsiSet &, const SxFermi &);

      /** \brief Get Harris-Foulkes energy */
      double getHarrisEnergy (const SxFermi &fermi);

      /** \brief Plane-wave preconditioner */
      virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &psi) const;

      /** Compute new density */
      virtual void computeRho (const Focc &, const SxPsiSet &);

      /** Compute pseudo-core density */
      SxMeshR computeCorePS () const;

      /** Compute subspace Hamiltonian matrix
          \note This is a memory-friendly blocked version
      */
      SxDiracMat<Complex16> getHnm (const PsiGI &waveBlock) const;

      /// Forces
      SxAtomicStructure forces;

      /// Compute forces given in Ref.1 eq(50,57,58,59)
      void computeForces (const SxPWSet &waves, SxFermi &fermi);

      /// Compute the gradient of the U-functional given in Ref.1 eq(28)
      SxAtomicStructure computeGradU ();

      /** \brief Computes the local compensation charge density given in Ref1. eq(22)

          @param is Species id
          @param ia Atom id within species
          @param YlmGl precomputed sperical harmonics
          @param rc2 squared gaussian width


      */
      PsiG computeRhoHatG (int is, int ia, SxDiracMat<Double> &YlmGl, double rc2);

      SxArray<SxRadMat> computeUnSymmGradDij (const SxFermi &fermi,
                                              const int iSpin,
                                              const int ik,
                                              SxDiracVec<Complex16> &P,
                                              SxArray<SxDiracMat<Complex16> > &gradP);

      /** \brief Scalar product (with gradient overlap operator) <i|\nabla O|j>
        @param iSpecies  which species
        @param iAtom     which atom
        @param i         state index
        @param j         other state index
        @param P         precompute  projection          npg x nStates
        @param gradP     precomputed projection gradient npg x nStates
        */
      SxComplex16
      pawGradDot (int iSpecies, int iAtom, int i, int j,
                  const SxDiracMat<Complex16> &P,
                  const SxDiracMat<Complex16> &gradP) const;

      SxArray<SxRadMat> symmetrizeGradDij (SxArray<SxRadMat> &GDij);

      /// Provide overlap operator
      virtual SxConstPtr<SxOverlapBase> getS () const;

      /// PAW density
      SxPAWRho pawRho;

      /** \brief Get access to density */
      virtual SxDensity& getRho () { return pawRho; }

};

namespace Timer {
   enum PAWHamTimer { Dij, ComputeU, ComputeHam,
                      PAW_H_psi, PAW_H_nl, PAW_H_loc, PAW_H_hubMO, RadXC,
                      ComputeHamRho, PAWdot, PAWdotP, ComputeGradU,
                      F11, F12, F13, F14, F2, F3, Force,
                      GradDij, pawGradDot, GDijSym,
                      khrMarker1, khrMarker2, khrMarker3};
}

SX_REGISTER_TIMERS(Timer::PAWHamTimer)
{
   regTimer (Timer::Dij, "PAWHam Dij");
   regTimer (Timer::ComputeU, "PAWHam: U");
   regTimer (Timer::ComputeHam, "PAWHam: compute");
   regTimer (Timer::ComputeHamRho, "PAWHam: rho terms");
   regTimer (Timer::PAW_H_psi, "PAWHam: H*psi");
   regTimer (Timer::PAW_H_nl , "PAWHam: H*psi nl");
   regTimer (Timer::PAW_H_loc, "PAWHam: H*psi loc");
   regTimer (Timer::PAW_H_hubMO, "PAWHam: Hubbard MO");
   regTimer (Timer::PAWdot, "PAW dot");
   regTimer (Timer::PAWdotP, "PAW dot: project");
   regTimer (Timer::ComputeGradU, "PAWHam: GradU");
   regTimer (Timer::F11, "PAWHam: Force 1-1");
   regTimer (Timer::F12, "PAWHam: Force 1-2");
   regTimer (Timer::F13, "PAWHam: Force 1-3");
   regTimer (Timer::F14, "PAWHam: Force 1-4");
   regTimer (Timer::F2, "PAWHam: Force 2");
   regTimer (Timer::F3, "PAWHam: Force 3");
   regTimer (Timer::Force, "PAWHam: Force tot");
   regTimer (Timer::GradDij, "PAWHam: GradDij");
   regTimer (Timer::pawGradDot, "PAWHam: PAWGradDot");
   regTimer (Timer::GDijSym, "PAWHam: GradDij Sym");
   regTimer (Timer::khrMarker1, "PAW::compute Marker1");
   regTimer (Timer::khrMarker2, "PAW::compute Marker2");
   regTimer (Timer::khrMarker3, "PAW::compute Marker3");
}



#endif /* _SX_PAW_HAMILTONIAN_H_ */
