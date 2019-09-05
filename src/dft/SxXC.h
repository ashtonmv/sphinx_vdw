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

#ifndef _SX_XC_H_
#define _SX_XC_H_

#include <SxPrecision.h>
#include <SxGBasis.h>
#include <SxRBasis.h>
#include <SxPseudoPot.h>
#include <SxDirac.h>
#include <SxPWSet.h>
#include <SxFermi.h>

#include <SxDFT.h>
/** This class computes the exchange-correlation functional
    on FFT-grids. The actual functionals are implemented in SxXCFunctional.

    Furthermore non- linear core correction
    (<a href="http://prola.aps.org/abstract/PRB/v26/i4/p1738_1">PRB,26,1738(1982)</a>)
    is treated in this class.
    Some member variables are of list type due to spin-polarization.
    For spin-unpolarized calculations the spin number is just zero.

    \par Usage:
    The simplest way to use this class is just to derive from it. Then you've
    got direct access to the members #eXc and #vXc, the exchange-correlation
    energy (given in Hartree) and the corresponding exchange-correlation
    potential.
    \code
       class MyNewClass : public SxXC
       {
          public:
             MyNewClass (SxGBasis *g) : SxXC (g) { };
             ...
             void compute ()  {
                ...
                computeXC();
                updateXC(rhoR);
                ...
             }
             void update ()  {
                ...
                updateXC (rhoR);
                eTotal = ... + eXc;
                for (iSpin=... )  {
                   vEff(iSpin) = ... + vXc(iSpin);
                }
                ...
             }
             ...
       };
    \endcode
    Please look into SxHamiltonian as well.
    \todo    SMP parallelization
    \ingroup group_dft
    \author  Sixten Boeck */
class SX_EXPORT_DFT SxXC
{
   public:

      /** Enumeration type for addressing spin numbers for several member
          variables */
      enum Type {UP=0, DOWN=1};

      enum XCFunctional  {
         Unknown    = -1,
         LDA        = 0,
         PBE        = 1,
         PBE_LDA    = 2,
         READ_VXC   = 3,
         EXX        = 4,
         EXX_LDA    = 5,
         EXX_PBE    = 6,
         EXX_PBE_WB = 7,
         Slater     = 8,
         KLI        = 9,
         LDA_PW     = 10,
         PBE0       = 32,
         HSE06      = 33
      };

      static bool isHybrid (XCFunctional xc)
      {
         return (xc == PBE0 || xc == HSE06);
      }

      XCFunctional xFunctional, cFunctional, xcFunctional;

      /** After calling SxXC::updateXC this member contains the 
          exchange-correlation energy in Hartree */
      PrecEnergy         eX, eC, eXc;
      /** Exchange energy in Hartree */
      PrecEnergy eExchange;
      /** Correlation energy in Hartree */
      PrecEnergy eCorrelation;
      /** Exchange-correlation potential energy in Hartree */  
      PrecEnergy         eVxc;
      /** Contains the exchange-correlation potential after returning from
          ::update */
      SxArray<SxMeshR>    vXc;                    // :iSpin
      /** Used for the mixing of different functionals*/
      SxArray<SxMeshR>    vEx, vCor;              // :iSpin

      /** It is possible to read in the xc potential from a file
          rather than computing it. This is the corresponding file
          name. It is read from the initialGuess group by SxHamSolver. 
         \brief External xc potential file name
        */
      SxString vXcFile;
      /// This reads the xc potential from vXcFile
      void readVXC (const SxRBasis &R);
      

      //---------------------------------------------------------------------
      /**@name Non-linear core correction*/
      //@{ Non-linear core correction
      /// Does any species have non-linear core corrections?
      bool nlcc;

      /** For NLCC only. The core density (per species) in reciprocal space. 
          (iSpecies:ig)
       */
      SxArray<SxDiracVec<TPrecPhi> >  phiCore;
      /** For NLCC only. The total core density in realspace. */
      SxMeshR                      rhoCoreR;
      /** This member variable is used for representing the sum of
          both, the electronic and core charge density. This entity enters
          the several XC functionals. */
      //@}
      //---------------------------------------------------------------------


      RhoR                         rhoXcR;       // NLCC: rhoR + rhoNLCC
      int nSpin;
      
      /// Empty constructor
      SxXC () : nSpin(0) { };
      /// Copy constructor
      SxXC (const SxXC &); 
      /// Constructor for XC with possible non-linear core correction
      SxXC (const SxPseudoPot &, const SxGBasis *, int nSpin_=1);
      /// Constructor for XC without non-linear core correction
      SxXC (int nSpin_);
      virtual ~SxXC ();

      /// Initialize when basis and nSpin are set
      void init ();

      /** This reads parameters provided in the input file and initializes 
          (some) member variables of the SxXC class.
      */
      void read (const SxSymbolTable *hamiltonian, const SxCell &cell);

      /** Set the xc basis from the input file
        */
      void setXcMesh (const SxSymbolTable *meshGrp, const SxCell &cell);

      /** \brief Get the xc potential specified in the input file
          \note This function will set SxXCFunctional::omegaHSE for 
                HSE functional.
      */
      static XCFunctional getXCFunctional (const SxSymbolTable *top);

      /** @{

          \brief Enable/disable exchange and/or correlation contribution

          This functions are used to control whether the exchange and/or
          the correlation contributions to the total energy and the
          effective potentials should be calculated.
          The functions modify SxXC::calcX and SxXC::calcC
          Used in particular in connection with SxPWHamiltonian::contrib */
      void enableExchange ();
      void disableExchange ();
      void enableCorrelation ();
      void disableCorrelation ();
      /** @}  */

      /** \brief status of the exchange functional

          @return false :  exchange functional is switched off
          @return true  :  exchange functional is switched on */
      bool getExchangeStatus () const;

      /** \brief status of the correlation functional

          @return false :  correlation functional is switched off
          @return true  :  correlation functional is switched on */
      bool getCorrelationStatus () const;

      /** Computes the exchange potential energy. Don't mix that up with the
          exchange energy! */
      void computeXPotEnergy (const RhoR &rhoR);

      /** Computes the exchange, correlation, and exchange-correlation
          potential energy. */
      void computeXCPotEnergy (const RhoR &rhoR);

      /** Computes, in case of presence of partial core corrected species,
          the NLCC formfactor. It reads
          \f[
            \Phi_{i_s}^{\rm core} 
               = \frac{4\pi}{\sqrt{(\Omega)}}\int_0^\infty dr
                           r^2 j_0(|G|r) 
                           \tilde{\varrho}^{\rm core}_{i_s}.
          \f]
          @note The projection onto the plane-wave basis will be done externaly
                in SxRadialBasis::toBasis. 
          \sa   <a href="http://prola.aps.org/abstract/PRB/v26/i4/p1738_1">PRB,26,1738(1982)</a>
          \todo logDr shouldn't be unique for all potential components.
          */
      void computePhiCore (const SxPseudoPot &, const SxGBasis *);


      /** Call computeXC whenever the atomic structure has changed. In that
          case (and the presence of partial core corrected species) the 
          core charge density 
          \f[
             \tilde{n}^{\rm core}(r) = \sum_{G} e^{iG \cdot r} 
                                       \sum_{i_a} S_{i_a}(G)
                                                  \Phi_{i_a}^{\rm core}(G)
          \f]
          has to be recomputed for they
          depends on structure factors! Please refer SxXC::computePhiCore for
          a detailed description of formfactor \f$ \Phi^{\rm core} \f$.
          Afterwards returning from computeXC you should manually call 
          SxXC::updateXC for computing the response to the electronic charge 
          density.*/
      virtual void computeXC ();
      /** If the atomic stucture has not changed since the last call but the
          electronic charge density has recomputed, call this function.
          Otherwise call SxXC::computeXC first.  
          
          It is a driver routine for the exchange-correlation functional in use.
          Dependently on SxControl::xcFunctional it calls one of the following
          routines:
             - LDA functionals: 
                 - Perdew/Zunger '81: SxXC::computeLDA or SxXC::computeLSDA
                 - Perdew/Wang '91: SxXC::computeLDA_PW or SxXC::computeLSDA_PW
             - GGA functionals: 
                 - Perdew/Burke/Ernzerhof: 
                   SxXC::updateGGA_PBE_WB

          @param rhoR electonic charge density/ies 
                      (spin unpolarized/polarized case)
      */
      virtual void updateXC (const RhoR    &rhoR,
                             const SxPWSet *wavesPtr=NULL,
                             const SxFermi *fermiPtr=NULL);


      //---------------------------------------------------------------------
      /**@name LDA routines */
      //@{
      /** Compute LDA for non-polarized densities */
      void computeLDA ();  // :{X,C}
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
      void computeLSDA ();
      void computeLDA_PW ();
      void computeLSDA_PW ();
      //@}
      //---------------------------------------------------------------------


      //---------------------------------------------------------------------
      /**@name GGA utility routines
         The following section describes routines that are used by (almost)
         all gradient-corrected XC functionals.*/
      //@{

      /** Driver function to compute gradient-corrected PBE functional.
       */
      void updateGGA_PBE_WB();

      /** Driver function to compute some gradient-corrected functional
       */
      void updateGGA_WB();


      /** Returns the 1st derivate by computing. 
          \f[
            \nabla m(R) = \sum_G iGe^{iG \cdot r} 
                                 \int_\Omega d^3 r e^{-iG\cdot r} m(R)
          \f]
          @param  meshR The input entity in real space. 
                        See formula above. Refered as \em m.
          @return \f$ \nabla m \f$
       */
      SxArray<SxMeshR> gradient (SxMeshR &meshR);
      //@}
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      /**@name Perdew/Burge/Ernzerhof */
      //@{
      //void testPBE ();
      //@}

      /** shifts the XC potential by setting its (G=0)-component to 0 */
      void shiftXC (bool xc=true, bool x=false, bool c=false);

      /** \brief special R-basis for xc potential & energy

          The non-linearity of the xc functional can lead to oscillations
          of the xc-energy when the complete system is rigorously shifted
          against the FFT mesh. The wiggles scale with the grid.

          It is possible to use a finer FFT mesh for the calculation of
          xc energies by setting the xcBasisPtr to the corresponding
          r-Basis. Fourier interpolation is used to switch between the
          standard Fourier mesh and the xc mesh.
        */
      SxConstPtr<SxRBasis> xcBasisPtr;

      /// Stores the default xc basis
      SxConstPtr<SxRBasis> origXcBasis;

      /** \brief Update xcBasisPtr from some input group
          @param table The input group to read from
          @return true if anything changed

          \note
          If table is NULL or there is no xcMesh group in the table,
          xcBasisPtr will be set to origXcBasis.
        */
      bool rereadTable (const SxSymbolTable *table);
   protected:

      friend class SxExx;
      bool calcX;
      bool calcC;
};

#endif /* _SX_XC_H_ */
