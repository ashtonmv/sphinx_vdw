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

#ifndef _SX_PAW_POT_H_
#define _SX_PAW_POT_H_

#include <SxDFT.h>
#include <SxSymbolTable.h>
#include <SxArray.h>
#include <SxDirac.h>
#include <SxSpeciesData.h>
#include <SxRadialMesh.h>
#include <SxSphereGrid.h>

/*
References:
 1. P. E. Bloechl, Phys. Rev. B 50, 17953 (1994).
*/


/** \brief Container for PAW Potentials

    \b SxPAWPot = S/PHI/nX Projector Augmented Waves Potential Container

    ....

    \author Sixten Boeck, boeck@sfhingx.de
    \author Christoph Freysoldt, freysoldt@mpie.de
 */
class SX_EXPORT_DFT SxPAWPot : public SxSpeciesData
{
   public:


      SxPAWPot ();
      SxPAWPot (const SxSymbolTable *);
      ~SxPAWPot ();

      void print () const;

   public:
      /** \brief read PAW potentials from AbInit- code 

          @param fileName file name
          @param is species id
      */
      void readAbInit (const SxString &fileName, int is);

      /** \brief Read core waves in abinit format
        @param file    file name
        @param is      which species
        */
      void readCoreAbinit (const SxString &file, int is);

      /** \brief read PAW potentials from AtomPAW - code 

          @param fileName file name
          @param is species id
      */
      void readAtomPAW (const SxString &fileName, int is);

   protected:
      /// Printout debug information when reading potentials?
      bool verbose;
      /** \brief read PAW potentials from CP-PAW code 

          @param fileName file name
          @param is species id
       
          The data format can be found in the CP-PAW source code
             PAW/src/Tools/Atom/paw_atom.f
          in the subroutine WRITEOUT
       */
      void readCPPAW (const SxString &fileName, int is);
      /** \brief read PAW potentials POTCAR from VASP-code

          @param fileName   file name
          @param is species id
          @param useProjG   directly projector G splines from POTCAR
      */
      void readVasp (const SxString &fileName, int is, bool useProjG);
   public:
 
      /// Logarithmic step size log (r(i+1)/r(i)) (:iSpecies)
      SxArray<Real8>  logDr,
      /// Cut-off radius for short-ranged compensation densities
                      rc,
     /// outmost Cut-off radius for wavefunctions
                      rCore;

      /// l-value for each projector (:is,:iPhi)
      SxArray<SxArray<int> >  lPhi; 

      /// l-value for each core wave (:is,:iPsi)
      SxArray<SxArray<int> >  lCore; 

      /// max. l value (:is)
      SxVector<Int> lMax;

      /// max. l value for density expansion
      SxVector<Int> lMaxRho;

      /// Which angular integration grid to use
      SxArray<SxSphereGrid::GridType> aGridType;

      /// Return max l of all species
      int getLMax () const { return lMax.maxval (); }

      /// Offset for projector (:is,:ipt)
      SxArray<SxArray<int> > offset;

      /// Radial mesh (:is)
      SxArray<SxDiracVec<Double> > rad,
      /// pseudo-core density (:is) / Y00
                                   rhoCorePS, 
      /// all-electron core density (:is) / Y00
                                   rhoCoreAE,
      /// Local smoothening potential vBar (:is) / Y00
                                   vBar;

      /** \brief Initialization density (:is) / Y00
        \note The basis type of this quantity may vary!
         - VASP:    radial G basis with linear mesh
         - ATOMPAW: standard radial basis
        */
      SxArray<SxDiracVec<Double> > rhoInit;

      /// Basis set for rhoInit
      SxArray<SxPtr<SxBasis> > rhoInitBasis;

      /// Compensation charge monopole for init density
      SxArray<double> rhoInitQ0;

      /// Basis set for pPS
      SxArray<SxPtr<SxBasis> > projBasis;

      /// Partial waves (pseudo) (:is,:iphi)
      SxArray<SxDiracMat<Double> > phiPS, 
      /// Partial waves (all-electron) (:is,:iphi)
                                   phiAE,
      /// Core waves (all-electron) (:is,:iPsi)
                                   psiCoreAE,
      /// Projector function (:is,:iphi)
                                   pPS,
      /// Projector function (:is,:iphi) on a fine grid (for r->G transform)
                                   pPsFine;

      /// Get all pseudo wave functions for a species
      SxArray<SxArray<SxDiracVec<Double> > > getPhiPS ();

      /// Get a single pseudo wavefunction
      SxDiracVec<Double> getPhiPS (int is, int ip);

      /** \brief Kinetic energy correction matrix

          \note See Ref. 1, Sec. VI E, last paragraph (iii)
          This is the kinetic energy differences between AE and PS partial
          waves. The all-electron partial waves use scalar relativistic
          corrections in the kinetic-energy expression.
          \f[
          \langle \phi_i | -\frac 12 \nabla^2 | \phi_j\rangle
          -
          \langle \tilde\phi_i | -\frac 12 \nabla^2 | \tilde\phi_j\rangle
          \f]
        */
      /// (:is)(:ipt,:jpt)
      SxArray<SxDiracMat<Double> > deltaKin,
      /** \brief Overlap correction matrix (:is)(:ipt,:jpt)

          \note See Ref. 1, Sec. VI E, last paragraph (iii)
          \f[
          \langle \phi_i | \phi_j\rangle
          -
          \langle \tilde\phi_i | \tilde\phi_j\rangle
          \f]
       */
                                   deltaS,
      /** \brief Core exchange matrix (:is)(:ipt,:jpt)

          \f[
          \langle \phi_i | \Sigma_x^{\rm core} | \phi_j\rangle
          \f]
       */
                                   coreX;

      /// Get number of species
      int getNSpecies () const  {
         return int(lPhi.getSize ());
      }

      /** \brief Get number of projectors for one species
          \note This is the number of m-dependent projectors. 
          \sa getNProjType
          The number of projector types is given by 
          \code
pawPot.getNProjType ();
          \endcode
        */
      int getNProj (int iSpecies) const;

      /** \brief Get number of projector types for one species
        \sa getNProj
        */
      int getNProjType (int iSpecies) const
      {
         SX_CHECK (iSpecies >= 0 && iSpecies < lPhi.getSize (),
                   iSpecies, lPhi.getSize ());
         return int(lPhi(iSpecies).getSize ());
      }

      int getNProjType (const SxAutoLoop &iSpecies) const
      {
         iSpecies.setLimit (getNSpecies ());
         return getNProjType(int(iSpecies.i));
      }

      /// Clebsch-Gordan coefficients
      SxYlm::SxClebschTable clebschGordan;

      /// Occupation numbers for initialization
      SxArray<SxVector<Double> > foccInit;
      /// Density matrix for initialization
      SxArray<SxMatrix<Double> > DijInit;

      /// Core energy (nuc + xc + Hartree)
      SxVector<Double> coreEnergy;

   protected:
      /** \brief Pointer to radial basis
        \note Use setBasis to also register the radial basis for all members
              that may need it.
        */
      SxConstPtr<SxRadBasis> radBasisPtr;

      /** \brief Pointer to fine radial basis
        This is exclusively used for |G+k><G+k|r><r|pPsFine>.
        It allows to change the resolution of the standard radial grid
        without affecting the projector quality.
        */
      SxPtr<SxRadBasis> fineRadBasisPtr;
      
      void createFineBasis (int is, double r0, double rMax, int nPts);
   public:

      /// Get reference to the radial basis
      const SxRadBasis &getRadBasis () const
      {
         SX_CHECK (radBasisPtr);
         return *radBasisPtr;
      }

      /** \brief Get SxPtr to radial basis
        \note In most cases, you'll prefer getBasis
        */
      const SxConstPtr<SxRadBasis> &getBasisPtr () const
      {
         return radBasisPtr;
      }

      /// Set radial basis
      void setBasis (const SxConstPtr<SxRadBasis> &radPtr);

      // \name Utility functions
      //@{
      /** \brief Compute the radial density from the density matrix
          
        */
      SxRadialMesh computeRho (const SxMatrix<Double> &Dij, int is,
                               const SxDiracMat<Double> &phi) const;

      /** Compute radial pseudo-density from the density matrix
          \note This adds the core contribution
        */
      SxRadialMesh computeRhoPS (const SxMatrix<Double> &Dij, 
                                 int is, int nSpin) const;

      /** Compute radial all-electron density from the density matrix
          \note This adds the core contribution
        */
      SxRadialMesh computeRhoAE (const SxMatrix<Double> &Dij,
                                 int is, int nSpin) const;

      /** \brief Get initial density for single atom in G-basis
          @param G the G basis
          @param iSpecies which species
        */
      SxDiracVec<TPrecCoeffG>
      getAtomRhoG(const SxGBasis &G, int iSpecies) const;


      /** \brief Compute \f$phi_i(r) * \phi_j(r)\f$.
        @param phi which phi's to take (should be phiAE or phiPS)
        @param ipt 1st partial wave type
        @param jpt 2nd partial wave type
      */
      inline static 
      SxDiracVec<Double> nij(const SxDiracMat<Double> &phi,
                             ssize_t ipt, ssize_t jpt)
      {
         SxDiracVec<Double> res = phi.colRef(ipt) * phi.colRef(jpt);
         res.handle->auxData.is = phi.handle->auxData.is;
         res.setBasis (phi.getBasisPtr ());
         return res;
      }

      /** Get normalized generalized Gaussians

        These are given by
        \f[
        grl = \frac{2} / {\Gamma(l+3/2)} r^l e^{-(r/r_c)^2}
        \f]
        In reciprocal space, they become

        \f[
        grl(G) = \frac{1}{(2l+1)!!} |G|^l e^{-\frac 14 r_c^2 |G|^2}
        \f]

        where (2l-1)!! = (2l+1)(2l-1)(2l-3)...1

        */
      SxDiracVec<Double> getGrl (int is, int l) const;

      /** \brief Get matrix elements of potential V
        */
      SxMatrix<Double> getVMatEl (const SxRadialMesh &V,
                                  const SxDiracMat<Double> &phi) const;
      //@}
      
      /** \brief Recompute overlap correction matrix

          \note Reasons to call this routine
          -# Enforce consistency between waves and deltaS
          -# Change in radial grid (interpolation)
          -# No deltaS in PAW file
        */
      void recomputeOverlap ();

     /** \brief  Check if PAW overlap operator is positive definite
         \param iSpecies species to check
         \param crashOnError if a failure should be ignored

         the overlap operator is 1 + \sum_ij |p_i>S_ij<p_j|
         => non-trivial eigenfunctions must be linear combination of |p_k>
         i.e.
         \f[\sum_j S_{ij} <p_j|p_k> c_k = (\lambda - 1) c_i\f]
         => all eigenvalues of S_ij <p_j|p_k> must be > -1
      */
      void checkOverlap (int is, bool crashOnError = true);

      /** \brief  extend Potential to new RMax value
       */
      void extendRad(double newRMax);

      /// Multipole moments (is)(ipt,jpt)(l)
      SxArray<SxArray2<SxVector<Double> > > QijL;

      /** \brief AE norm within PAW radius (is)(ipt,jpt)

        \f[ \omega_{ij} = \int_0^{r_{PAW}} r^2 dr \phi^{AE}_i(r)\phi^{AE}_j(r)
                        * \delta_{l_i l_j} \f]

        This quantity is needed e.g. for spin constraints. It is computed in
        setupQijL
        */
      SxArray<SxMatrix<Double> > omegaPAW;

      /// Set up QijL and omegaPAW
      void setupQijL ();

      /** 4-index exchange matrices (radial part) (:is)(:L)(:ijkl)

        \note the 4-index is obtained by get4Idx

        */
      SxArray<SxArray2<double> > xKernel;

      /** \brief Moments of the PAW corrections (:is)(:m,:L,:ij)
        */
      SxArray<SxArray3<double> > MijL;
      
      /// Compute core-valence exchange
      void computeCoreX ()  
      {
         for (int is = 0; is < getNSpecies (); ++is)
            if (psiCoreAE(is).getSize () > 0)
               computeCoreX (is);
      }

      /// Compute the 4-index exchange matrices 
      void computeXKernel ();

      /// Compute core-valence exchange integrals
      void computeCoreX (int is);

      /// Get condensed 2-index
      static inline int get2Idx (int i, int j, int N)  {
         SX_CHECK (i >= 0 && i < N, i, N);
         SX_CHECK (j >= 0 && j < N, j, N);
         if (i < j) return (i * (2 * N - i - 1)) / 2 + j;
         else       return (j * (2 * N - j - 1)) / 2 + i;
      }

      /// Condense the 4-index 
      static inline int get4Idx (int i, int j, int k, int l, int N)
      {
         return get2Idx (get2Idx (i, j, N), get2Idx(k, l, N), (N*(N+1))/2);
      }

      /// Compensation charges for Kresse-Joubert xc treatment
      SxArray<SxDiracVec<Double> > kjXcShape;

};

#endif /* _SX_PAW_POT_H_ */
