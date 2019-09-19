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

#ifndef _SX_PW_HAMILTONIAN_H_
#define _SX_PW_HAMILTONIAN_H_

#include <SxHamiltonian.h>
#include <SxQuantumNumbers.h>
#include <SxXC.h>
#include <SxRBasis.h>
#include <SxRadBasis.h>
#include <SxPW.h>
#include <SxMuPW.h>
#include <SxRho.h>
#include <SxPseudoPot.h>
#include <SxAtomicStructure.h>
#include <SxTypes.h>
#include <SxFermi.h>
#include <SxEESGProj.h>
#include <SxRSProj.h>
#include <SxDFT.h>
#include <SxVDW.h>

/** \brief Plane-wave Hamiltonian \f$
              \langle \bf {G+k} | \hat{H} | \bf{G+k} \rangle
           \f$

    \b SxPWHamiltonian = S/PHI/nX Plane-Wave Hamiltonian
           
    This class represents the Hamiltonian in a plane-wave basis-set. The
    basis-set is a ::SxPWSet or any object derived from it (this allows to use
    actually any kind of basis-set to be projected onto plane-waves on demand).

    The class computes the total energy and its energy contributions, the
    gradient \f$ H |\Psi\rangle \f$, as well as the force contributions.

    The major interface routines are ::compute, ::update, ::getEnergy, and
    the overloaded * operator (::operator*).

    \sa      \ref page_dirac
    \sa      SxPW
    \ingroup group_dirac
    \ingroup group_dft
    \author  Sixten Boeck
   */
class SX_EXPORT_DFT SxPWHamiltonian : public SxHamiltonian
{
   public:

      /** \brief Hamiltonian contribution identifies
          \sa    ::contrib */
      enum Contrib { CALC_NONE = 0x00, CALC_KIN = 0x01, CALC_HARTREE = 0x02,
                     CALC_X    = 0x04, CALC_C   = 0x08,
                     CALC_LOC  = 0x10, CALC_NL  = 0x20, CALC_EXT= 0x40,
                     CALC_SCR  = 0x80,
                     CALC_RHO  = 0x100,
                     CALC_XC   = CALC_X | CALC_C,
                     CALC_EFF  = CALC_HARTREE | CALC_XC | CALC_LOC | CALC_RHO,
                     CALC_DEFAULT = CALC_KIN | CALC_EFF | CALC_NL | CALC_SCR,
                     CALC_ALL  = 0xffff };
      enum Preconditioner { Payne, Arias };

      /** \brief Empty standard constructor

          This constructor is just a placeholder. To instatiate the 
          PW Hamiltonian properly the wavefunctions, the charge density,
          and the atomic structure must be provided to the constructor. */
      SxPWHamiltonian ();

      /** \brief Setup the PW Hamiltonian

          This constructor initializes the necessary memory and computes
          all required formfactors by calling 
          - computePhiGauss, 
          - computeESelf,
          - computePhiLocPseudo, and
          - computePhiNonLocal
          Note, that this constructor does \b not call ::compute or
          ::update! They have to be called explicitly.
          \param waves  Bloch wave coefficients \f$ \Psi_{i \sigma k}\f$
          \param rhoR   Electronic charge density \f$ \varrho(R) \f$ */
      SxPWHamiltonian (const SxPWSet &, const RhoR &, 
                       const SxPseudoPot &, const SxAtomicStructure &,
                       bool rsProjNew = false);
      /** \brief Setup the PW Hamiltonian
          @param G           global G basis (with registered R and structure)
          @param Gk          global G+k basis
          @param potPtr      pseudopotentials
          @param rhoFile     filename of sxb file with PAW density
          @param table       symbol table (for xc etc.)
          @param wavesPtrIn  symbol table (for xc etc.)

          This constructor initializes the Hamiltonian ready to use
          from a density in a file, a symbol table, and the necessary
          basis sets.

          For computing the energy, however, the wavesPtr must be set, too,
          and the ::compute routine must be called again.

          */
      SxPWHamiltonian (const SxGBasis  &G,
                       const SxGkBasis &Gk,
                       const SxPtr<SxPseudoPot> &potPtr,
                       const SxString  &rhoFile,
                       const SxSymbolTable *table,
                       const SxPtr<SxPWSet> &wavesPtrIn = SxPtr<SxPWSet> ());

      /** \brief Destructor of the PW Hamiltonian */
      virtual ~SxPWHamiltonian ();

      /// Get the PWHamiltonian group or return NULL
      static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);

      /** \brief Read user input parameter from the input file parser

          Read the group "PWHamiltonian{}" from the input file and setup the
          local member variables correspondingly. These are
          - xFunctional, cFunctional, and xcFunctional
          - electronic smearing temperature ::ekt
         */
      virtual void read  (const SxSymbolTable *);
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

      /** \brief check the number of empty states

          When diagonalising the Hamiltonian in G space representation, the
          maximal number of states you can obtain is equal to the number of
          G vectors used to represent the Hamiltonian in matrix form.

          Trying to get more than \f$nf$ states from a
          \f$|\mathbf{G}+\mathbf{k}\rangle\f$ basis set with \f$n\f$ G vectors
          is not meaningful and may lead to convergence problems or other
          problems. This function is therefore meant to check, whether the
          number of states being tried to be accessed is within an appropriate
          range. */
      void validateNStates (int nStates);   

      /** \brief (Re)compute all contributions when atoms have moved.

          This routine updates all contributions of the Hamiltonian,
          even the formfactors. It has to be called whenever the atoms
          have moved (e.g. after each structure optimization or
          molecular dynamics step). This routine recomputes first the
          formfactors be calling
          - ::computeRhoGauss
          - ::computeScrPotential
          - SxXC::computeXC
          It also recomputes ::vLocR by FFT of
          \f[
             v_{\mathrm{loc}}^{\mathrm{ps}}(\mathbf{G}) 
                = S(\mathbf{G}) 
                  \Phi_{\mathrm{loc}}^{\mathrm{ps}}(\mathbf{G})
          \f]
          Afterwards ::update is called in order to recompute also those
          terms of the PW Hamitltonian which depend only on the charge
          density and/or the wavefunctions.

          \param fermi  Current Fermi distribution \f$ f^{\mathrm{occ}}_{i \sigma k}\f$
          \param tauChanged recompute all terms depending on \f$ \tau_{i_s i_a}\f$
          \param rhoChanged recompute all terms depending on \f$ \varrho(R) \f$
          */
      virtual void compute (const SxFermi &,
                            bool tauChanged=true, bool rhoChanged=true);

           /** \brief Compute new density
          @param focc occupation numbers
          @param psi wave functions
      */
      virtual void computeRho (const Focc &focc, const SxPsiSet &psi);
 
      /** \brief Recompute all terms when \f$ \varrho(R) \f$ and/or \f$ \Psi_{i \sigma k}\f$ have change but not \f$ \tau_{i_s i_a}\f$.

          In order to safe computer time formfactors (which depend
          on the atomic positions \f$ \tau_{i_s i_a} \f$ are computed and
          saved in member variables (see ::compute). 
          During the electronic minimization these formfactors are constant
          and need not to be recomputed. This routine updates only terms
          depending on the charge density and/or the wavefunctions.
          It computes
          - the kinetic energy contribution ::eKin
          - the Hartree contribution ::vHartreeR, ::eHartree, and ::eHartreeElec
          - the local pseudopotential energy ::eLocPs
          - the XC potential by calling SxXC::updateXC
          - and (according to ::contrib) ::vEffR */
       void update (const SxFermi &);

      /** \todo caching of psi -> <PsiH|Psi> */
      PsiG operator* (const PsiG &) const;
      PsiG H_PsiG  (const PsiG &) const;
      PsiG H_PsiGI (const PsiG &) const;
      virtual PsiG apply (const PsiG &psi) const 
      {
         return SxPWHamiltonian::operator* (psi); 
      }

      SxDiracVec<TPrecCoeffG::TReal> 
         preconditioner (const PsiG &, Preconditioner type) const;

      virtual SxDiracVec<TPrecCoeffG::TReal> preconditioner (const PsiG &psi) const {
         return preconditioner (psi, Payne);
      }

      /// Density
      SxRho rho;
      
      /** \brief External potential in R space (:iSpin,:iR)
        
          This is an additional external local potential and must
          not be confused with the "external potential" from
          electron structure theory, where the ion-electron interaction
          is called the external potential.

          \sa vExtActsOnNuclei
        */
      SxArray<SxMeshR>    vExtR;   // :iSpin
      /** Whether or not vExt acts on the nuclei */
      bool vExtActsOnNuclei;
      /** \brief Effective potential in R space (:iSpin,:iR)
        */
      SxArray<SxMeshR>    vEffR;   // :iSpin
      const SxPWSet      *wavesPtr;
      const SxFermi      *fermiPtr;
            SxPseudoPot   psPot;   // TODO: should be const
      const SxGBasis     *gBasisPtr;

      double              ekt;

      /** \brief Block size for blocked non-local pseudopotential algorithm
          
          The non-local pseudopotential includes a sum over the orbitals.
          In the standard algorithm, each orbital is constructed from the
          reference orbital and the phase factor and is then used.

          In the blocked algorithm, nlBlocksize orbitals are constructed at
          a time, and then used together which allows to use efficient
          matrix routines. This is still a naive approach as it doesn't
          exploit intermediate summations over orbitals or atoms that would
          allow even more efficient algorithms.
        */
      int nlBlockSize;
      
      /** Pointer to the xc functionals */
      SxPtr<SxXC> xcPtr;
      
      // meshDensity parameter for calculation of XC contributions
      int xcMeshDensity;

      Contrib contrib;


      // --- Energy contributions
      PrecEnergy  eKin;
      PrecEnergy  eHartreeElec;
      PrecEnergy  eHartreeGauss;
      PrecEnergy  eHartree;
      PrecEnergy  eLocPs;
      PrecEnergy  eExt;
      PrecEnergy  eNl;
      PrecEnergy  eEwald;
      PrecEnergy  eSelf;
      PrecEnergy  eTotal;
      PrecEnergy  eExternal;
      PrecEnergy  eDoubleCounting;
      PrecEnergy  eVDW;

      bool applyVDWCorrection;
      SxVDW vdwCorrection;

      // --- Hartree contribution
      SxArray<SxDiracVec<TPrecPhi> >  phiGaussG;   //:is,:ig
      SxMeshG                         rhoGaussG;   // :ig
      SxMeshR                         vHartreeR;
      void computePhiGauss (const SxGBasis &);
      void computeRhoGauss (const SxGBasis &);
      void computeESelf ();
      bool dipoleCorrection;
      // external field along z when dipole correction is used
      double zField;
      // last value of zAlign for charged slabs
      double zAlignLast;

      // --- Forces
      SxAtomicStructure  fHartree;
      SxAtomicStructure  fLocPs;
      SxAtomicStructure  fNl;
      SxAtomicStructure  fXC;
      SxAtomicStructure  fScr;
      SxAtomicStructure  fTotal;
      bool calcForces;

      // --- temp
      SxAtomicStructure structure;

      // --- local pseudopotential contribution
      SxMeshR                         vLocR;
      SxArray<SxDiracVec<TPrecPhi> >  phiLocPseudo;  //:is,:ig
      void computePhiLocPseudo (const SxGBasis &);

      // --- non-local contribution

      /** \f$\langle \mu | V_{ps} | \mu \rangle \f$ */
      SxArray<PrecEnergy>             eKB;       //:{is,n,l,m}

      // sum_r < G+k | r > < r | Vnl |  >
      SxArray<SxArray<PsiG> >         phiNl;     //:ik,:{is,n,l,m},:ig
      SxPtr<SxAOBasis>                nlProj;

      SxArray<int>                    phiOrbNl;  //:{is,n,l,m,ik}
      SxArray<SxQuantumNumbers>       psiOrbNl;  //:{is,ia,n,l,m,ik}


      void computePhiNonLocal (const SxGkBasis &, bool rsProjNew = false);
      void computeVNl (int iSpin, int ik, const Focc &focc);

      void computeScrPotential ();

      /** Returns the Hamiltonian (at a k-point \f$ \mathbf k \f$
          in matrix form, i.e.
          \f$ \langle \mathbf{G}+\mathbf{k}
                  H | \mathbf{G'}+\mathbf{k} \rangle \f$

          @param   gkBasis the set of |G+k> basises
          @param   iSpin   spin quantum number
          @param   ik      k-point
          @param   nG      How many G vectors are considered in
                           \f$H_\mathbf k\f$.
                           If not provided or set to -1, the full matrix (up to
                           the number of plane waves at \f$ \mathbf k \f$) will
                           be evaluated. Othewise the corresponding subspace
                           is used.
          @returns \f$ H_{\mathbf k}(\mathbf{G},\mathbf{G'}) :=
                   \f$ \langle \mathbf{G}+\mathbf{k}
                           H | \mathbf{G'}+\mathbf{k} \rangle \f$
       */
      SxDiracSymMat<TPrecCoeffG> getMatrix (const SxGkBasis &gkBasis,
                                            int iSpin, int ik, int nG=-1);

      // --- Utitlity functions
      virtual PrecEnergy getEnergy (const SxPsiSet &, const SxFermi &);
      /// Return last computed energy
      virtual PrecEnergy getEnergy () const
      {
         return eTotal;
      }
      virtual PrecEnergy getDoubleCounting () const
      {
         return eDoubleCounting;
      }
      virtual void writeRho (const SxString &) const;

      void set (const SxPWSet &, const SxFermi &);
      
      void setXCMeshDensity (int);
      /// Return reference to density
      virtual SxDensity& getRho ();

      inline const SxPWSet &getWavesRef () const { 
         SX_CHECK (wavesPtr);
         return *wavesPtr; 
      }
      void printEnergies () const;
      void writeVElStat (const SxString &) const;

      void printContribMask () const;

      bool printHartree;

      SxPtr<SxEESGProj> rsProjector;
      SxPtr<SxRSProj> rsProj2;
      
      SxMatrix<TPrecCoeffG> getNlMatrix (const SxMuPW &mu, int ik,
                                         const SxAtomicStructure &str) const;

};

namespace Timer {
   enum PWHamTimer {
      hUpdate,
      Kin,
      Hartree,
      XC,
      locPs,
      nonLoc,
      hPsi,
      tPsi,
      vEffPsi,
      vNlPsi,
      phiGauss,
      phiLocPs,
      phiNonLoc,
      ewald
   };
}

SX_REGISTER_TIMERS (Timer::PWHamTimer)
{
   using namespace Timer;
   regTimer (hUpdate,    "H update");
   regTimer (Kin,        "Kinetic energy");
   regTimer (Hartree,    "Hartree potential");
   regTimer (XC,         "Exchange-correlation");
   regTimer (locPs,      "Local pseudo pot.");
   regTimer (nonLoc,     "NL pseudo pot.");
   regTimer (phiGauss,   "Phi Gauss");
   regTimer (phiLocPs,   "Phi Local");
   regTimer (phiNonLoc,  "Phi non-local");
   regTimer (ewald,      "Ewald summation");
   regTimer (hPsi,       "|X> = H |Psi>");
   regTimer (tPsi,       "T |Psi>");
   regTimer (vEffPsi,    "vEff |Psi>");
   regTimer (vNlPsi,     "vNl |Psi>");
}

#endif /* _SX_PW_HAMILTONIAN_H_ */
