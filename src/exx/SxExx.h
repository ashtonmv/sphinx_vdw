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

#ifndef SX_EXPORT_EXX
#ifdef WIN32
#  if defined(_EXPORT_sxexx)
#     define SX_EXPORT_EXX __declspec(dllexport)
#  else
#     define SX_EXPORT_EXX __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_EXX
#endif
#endif

#ifndef _SX_EXX_H_
#define _SX_EXX_H_

#include <SxXC.h>
#include <SxPW.h>

class SxFockGk;
class SxSymFockGk;

/** \brief EXX computation

    \b SxClass = S/PHI/nX EXX class

    

    \author M. Wahn
    Reintegrated by C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_EXX SxExx : public SxXC
{
   public:

      /// Constructor
      SxExx ();

      /// Copy constructor
      SxExx (const SxXC &xc);

      virtual ~SxExx ();

      /** This reads parameters provided in the input file and initializes 
          (some) member variables of the SxExx class.
      */
      void read (const SxSymbolTable *table, const SxGBasis &G);

      /** exchange (Fock) operator in G space, used in connection with EXX
          calculations (non-symmetrized version) */
      SxPtr<SxFockGk>              fockPtr;

      /** array of exchange (Fock) operators, used in connection with gsEXX
          calculations (non-symmetrised version) */
      SxArray<SxPtr<SxFockGk> >    fock;

      /** exchange (Fock) operator in G space, used in connection with EXX
          calculations (symmetrized version) */
      SxPtr<SxSymFockGk>           symFockPtr;

      /** array of exchange (Fock) operators, used in connection with gsEXX
          calculations (symmetrised version) */
      SxArray<SxPtr<SxSymFockGk> > symFock;

      /// Update EXX potentials
      virtual void updateXC (const RhoR &rhoR, const SxPWSet *wavesPtr,
                             const SxFermi *fermiPtr);
      /** Computes the contribution of the exchange (Fock) operator to the
          system's total energy. Overwrites eX, eXc, and eExchange with the
          computed value. */
      void computeFockEnergy (const SxPW &waves, const SxFermi &fermi);
      void computeSymFockEnergy (const SxPW &waves, const SxFermi &fermi);

      /** When using a Monkhorst-Pack mesh for k-points, and especially if it's
          reduced, the formula for \f$ \chi_{0}({\mathbf G}, {\mathbf{G'}}) \f$
          given in <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 59, 10031 (1999) </a>
          does no longer hold. The error, indeed, is not that much, so that
          you can quiescently use it. However, if you intend to get correct
          results, you should use the flag \"useCorrectChi\". */
      bool useCorrectChi;

      /** the number of
          \f$\mathbf{G}\f$ vectors considered in the linear response function
          \f$\chi_{0}(\mathbf{G}, \mathbf{G'})\f$, the non-local exchange
          (Fock) operator \f$ V^{\mathrm Fock}(\mathbf{G}, \mathbf{G'})\f$,
          and the induced charge-density
          \f$\Delta\rho(\mathbf{G}) = \rho_{\mathrm ind}(\mathbf{G})\f$ */
      int                nGChi;

      //-----------------------------------------------------------------------
      /**@name Standard Exact Exchange Formalism

         The way the Standard EXX formalism has been implemented currently
         is not the most effective one concerning memory and CPU time. However,
         the implementation is therefore as 'modular' as possible, in the
         sense, that every mathematical/physical entity got an extra class:

         - SxChi0 for the computation of the linear response function
           \f$ \chi_{0}({\mathbf G}, {\mathbf G'}) \f$
         - SxFockGk for the computation of the exchange (Fock) operator
           \f$ V_{\rm x} \f$ in G space,
         - SxRhoIndG for the computation of the change
           \f$ \rho_{\rm ind}^{\rm stEXX} \f$ in charge density,
           caused by the influence of the exchange operator
           \f$V_{\rm x}\f$.
       */
      //-----------------------------------------------------------------------
      //@{
      /** This routine computes the EXX potential using the "standard method"
          presented in
          <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 59, 10031 (1999) </a>. The correlation part is set to 0. */
      void computeEXX_noCorr (const SxPWSet &waves, const SxFermi &fermi);

      /** This routine computes the exchange potential using the "Exact
          Excange" (EXX) formalism presented in
          <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 59, 10031 (1999) </a>. The correlation part is computed with LDA. */
      void computeEXX_LDA (const SxMeshR &rho,
                           const SxPWSet &waves, const SxFermi &fermi);

      /** This routine computes the exchange potential using the "Exact
          Excange" (EXX) formalism presented in
          <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 59, 10031 (1999) </a>. The correlation part is computed with PBE. */
      void computeEXX_PBE (/*const*/ SxMeshR &rho,
                           const SxPWSet &waves, const SxFermi &fermi);

      /** This routine computes the exchange potential using the "Exact
          Excange" (EXX) formalism presented in
          <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev. B, 59, 10031 (1999) </a>. The correlation part is computed with PBE (WB),
          i.e., a discretized version of the PBE functional. */
      void computeEXX_PBE_WB (/*const*/ SxMeshR &rho,
                              const SxPWSet &waves, const SxFermi &fermi);
      //@}

      //-----------------------------------------------------------------------
      /**@name Slater potential and KLI Formalism */
      //-----------------------------------------------------------------------
      //@{
      /** Computes an approximative exchange potential employing the so-called
          KLI method introduced by Krieger, Li, and Iafrate. See also
          <a href="http://prola.aps.org/pdf/PRA/v46/i9/p5453_1">
          Phys. Rev. A, 46, 5453 (1992) </a>. Its average valus is set to 0
          to facilitate comparison with EXX, Slater, and so on. Correlation
          is not considered.
       */
      void computeKLI_noCorr (const SxMeshR &rho,
                              const SxPWSet &waves, const SxFermi &fermi);

      /** The X potential is treated by means of an approximation proposed
          by Slater in
          <a href="http://prola.aps.org/abstract/PR/v81/i3/p385_1">
          Phys. Rev. 81, 385 (1951) </a>.
          Its average value is set to 0 in order to facilitate comparisons
          with EXX, KLI, and so on. Correlation is not considered.
       */
      void computeSlater (const SxMeshR &rho,
                          const SxPWSet &waves, const SxFermi &fermi);

      /** Computation of the Slater potential, proposed by Slater in
          <a href="http://prola.aps.org/abstract/PR/v81/i3/p385_1">
          Phys. Rev. 81, 385 (1951) </a>.
          @param  kli true, if routine is used for KLI or LHF
          @return the Slater potential
       */
      VxcR getVSlater (const SxMeshR &rho, const SxPWSet &waves,
                       const SxFermi &fermi, bool kli=false);
      //@}

      //-----------------------------------------------------------------------
      //   gsEXX  <-->  stEXX
      //-----------------------------------------------------------------------
      PsiG greensfunction (const PsiG &h1psi0,
                           const SxPW &waves, const SxFermi &fermi);

      PsiG probier (const PsiG &vExG, const SxPW &waves, const SxFermi &fermi);

      void compareRhoInd (const PsiG &rhoIndG, const PsiG &vExG,
                          const SxPW &waves, const SxFermi &fermi);

      void computeVExUnsym (const SxPWSet &waves, const SxFermi &fermi,
                            const SxGBasis &G);

      /** Unsymmetrised exchange potential */
      SxList<SxMeshG>    vExUnsym;               // :iSpin
      SxList<SxMeshG>    vSlaterUnsym;           // :iSpin
      SxList<SxMeshG>    vExLoopG;               // :iSpin

      /** the expectation values of the summands of the Slater potential */
      SxArray<SxArray<PrecRhoR> > slaterExpect;  // :ik, iv
};

#endif /* _SX_EXX_H_ */
