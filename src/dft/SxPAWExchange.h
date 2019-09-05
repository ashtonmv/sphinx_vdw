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

#ifndef _SX_PAW_EXCHANGE_H_
#define _SX_PAW_EXCHANGE_H_

#include <SxDFT.h>
#include <SxPAWPot.h>
#include <SxRadMat.h>
#include <SxPWSet.h>
#include <SxRBasis.h>
#include <SxGBasis.h>
#include <SxPartialWaveBasis.h>

/** \brief PAW exchange operator

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWExchange
{
   public:
   
      /// Compensation moments (complex)
      typedef SxArray<SxArray<SxVector<Complex16> > > SxQlm;

      /// Pointer to the PAW potential
      SxConstPtr<SxPAWPot> potPtr;

      /// Pointer to partial wave basis for exchange waves
      SxConstPtr<SxPartialWaveBasis> pBasisX;

      /// Pointer to partial wave basis
      SxConstPtr<SxPartialWaveBasis> pBasis;

      /// Pointer to the exchange waves
      SxConstPtr<SxPWSet> wavesPtr;

      /// Pointer to SxRBasis
      SxConstPtr<SxRBasis> rPtr;

      /// Pointer to SxGBasis
      SxConstPtr<SxGBasis> gPtr;

   protected:
      /// Occupations
      Focc focc;
   public:
      void setFocc (const Focc &foccIn);

      /** \brief Symmetry group for real Ylm (:iSym,:l,(l+m x l+m) )
        */
      SxYlmRotGroup ylmRot;

      // -------------------------------------------------------

      /// Constructor
      SxPAWExchange ();

      /// Compute exchange energy
      double compute (const SxConstPtr<SxPWSet> &wavesIn,
                      const Focc &foccIn,
                      const SxRadMat &Dij);

      /// Apply exchange operator
      PsiG apply (const PsiG &psi) const;

      /// Compute on-center exchange operator
      void computeXij (const SxRadMat &Dij);

   protected:
      /// on-center exchange operator matrices (:iSpin)
      SxRadMat Xij;

      double computeXijSum (int is,
                            const SxMatrix<Double> &DijR,
                            const SxNArray<double, 2> &K,
                            int ipt, int l1, int lm1,
                            int lpt, int l4, int lm4);

      ///@name symmetry-specific data
      ///@{
      /// v(G+k-q)
      mutable SxArray<SxDiracVec<Double> > vkq;
      /// \f$|\mathbf{G+k-q}|^l Ylm(\mathbf{G+k-q})\f$
      mutable SxArray<SxDiracMat<Complex16> > YlmGl;
      /// map ir -> S^-1(ir)
      mutable SxArray<SxArray<ssize_t> > symRot;
      /// Pointer to projectors for soft compensation charges
      mutable SxPtr<SxAOBasis> pwProjComp;
      /// Index of G where G+k-q=0 (-1 if there is none) (:iSym)
      mutable SxArray<int> idxGkq0;
      /// Hard Gaussian shapes
      mutable SxArray<SxDiracMat<Double> > gHard; 
      /// Hard Gaussian shapes
      mutable SxArray<SxDiracVec<Double> > gSoft;
      /// Small G vectors
      mutable SxArray<SxArray<ssize_t> > smallG;
      /// Small G values
      mutable SxArray<SxDiracVec<Double> > smallG2;
      /// vLR for small G values
      mutable SxArray<SxDiracVec<Double> > smallV;

      ///@}

      /// U kernel (:iTlAtom, :iTlNeighbor) (lm x lm)
      SxArray<SxArray<SxMatrix<Double> > > kernelU;
      /// Neighbor list
      SxArray<SxAtomicStructure> neighborDist;
      /// Compute the U kernel
      void computeUKernel ();

      /// Determine V=0 component
      double computeV0 (const Coord &kVec) const;

      /** \brief Compute data for all rotations of q
           - vkq         v(k-q)
           - symRot      real-space rotation map
           - YlmGl       |G+k-q|^l Ylm(G+k-q)
           - pwProjComp  hard compensation projectors
        */
      void prepareRotData (int iq, Coord &kVec) const;

      /// Clear the rotated-q data
      void rotCleanup () const;

      /** \brief Compute compensation part of pseudo-exchange

        @param vG 2-orbital exchange potential with soft(!) compensation charges
        @param dkq2 \f$\mathbf k-\mathbf q\f$
        @param qlm 2-orbital compensation moments
        @return The exchange with normalized (hard) compensation charges.

        This computes
        \f[
        \int d^3 \mathbf r g_{RL}(\mathbf r) V_{ab}
        \f]

        */
      SxQlm gRlExchange (const PsiG &vG,
                         const Coord &dkq,
                         const SxQlm &qlm) const;

      /** \brief Compute hard/soft correction in compensation exchange

        @param dkq2 \f$\mathbf k-\mathbf q\f$
        @param qlm 2-orbital compensation moments
        @param grl if provided, compute exchange with single compensation
               moments instead

        This computes
        \f[
        \int d^3 \mathbf r [\hat n_{ab}^*(\mathbf r) g_{RL}(\mathbf r) 
                           V[\hat n_{ab} - \hat n_{ab}']
        \f]

        */
      SxComplex16 qRlExchange (const Coord &dkq,
                               const SxQlm &qlm,
                               SxQlm *grl = NULL,
                               const SxDiracVec<TPrecFocc> &foccJ 
                                   = SxDiracVec<TPrecFocc> ()    ) const;

      /** \brief Compute compensation part of pseudo-exchange
        @param vG 2-orbital exchange potential with soft(!) compensation charges
        @param dkq2 \f$\mathbf k-\mathbf q\f$
        @param qlm 2-orbital compensation moments
        @param pX projection onto exchange-state
        @return pseudo-exchange in partial-wave basis
        */
      SxDiracVec<TPrecCoeffG> 
      compensationExchange (const PsiG &vG,
                            const Coord &dkq,
                            const SxQlm &qlm,
                            const SxDiracVec<TPrecCoeffG> &pX,
                            const SxArray<SxComplex16> &n0) const;

      /** \brief Get multipole moments of wave function product
        */
      SxQlm getQlm (const SxDiracVec<TPrecCoeffG> pA,
                    const SxDiracVec<TPrecCoeffG> pB) const;

      /** \brief Compute PAW corrections in reciprocal space for 
                 \f$|G|<13 \omega_{HSE}\f$
          @param pa projections for state a
          @param pb projections for state b
          @param iSym current symmetry index
          @param dkq current k-q
          @param Mab0 here, the 2nd-order expansion coefficients at G+k-q=0 is
                 stored (needed for 2nd-order correction of potential)

          This is needed for calculating the screening of the PAW corrections.
        */
      PsiG computePAWCorrG (const SxDiracVec<TPrecCoeffG> &pa,
                            const SxDiracVec<TPrecCoeffG> &pb,
                            int iSym,
                            const Coord &dkq,
                            SxComplex16 *Mab0);
      /// Get |G+k-q|^2
      SxDiracVec<Double> getGkq2(const Coord &dkq) const;

      /** \brief Get compensation charge phases from multipole moments

        @param qlm    multipole moments
        @param iSym   symmetry id
        @param dkq    \f$\mathbf{k-q}\f$ for translational phase factors 

          \f[
          \phi_{is}(\mathbf G) 
          = \frac{4\pi}{\sqrt{\Omega}}
            \sum_{ia} e^{i (\mathbf{G+k-q})r_{is,ia}}
            \sum_{l,m} Q^{lm}_{is,ia} Y_{lm}(\mathbf{G+k-q})
            \frac{1}{2l+1} |\mathbf{G+k-q}|^l
          \f]

          By multiplying these phases with the reciprocal space
          shape function (Gaussians in our case), the compensation
          charges can be computed. By precomputing the phase sums,
          hard and soft compensation charges can be easily computed.
          \sa getNHat, getNHatSoft
        */
      SxArray<PsiG> getNHatUnshaped (const SxQlm &qlm,
                                     int iSym,
                                     const Coord &dkq) const;
      /** \brief Get hard compensation charge density from phases
          @param j which column to use from nHatUnshaped
        */
      PsiG getNHat (const SxArray<PsiG> &nHatUnshaped, ssize_t j = 0) const;

      /** \brief Get soft compensation charge density from phases
          @param j which column to use from nHatUnshaped
        */
      PsiG getNHatSoft (const SxArray<PsiG> &nHatUnshaped, ssize_t j = 0) const;

      /** \brief Get rotated set of projections
          @param proj projections
          @param iSym  id of symmetry
        
          This returns projections
          \f[
          \langle p_{S(i)} | \psi_{n\mathbf k}\rangle 
          = \langle p_i | \psi_{n S(\mathbf k)}\rangle
          \f]
        */
      SxDiracVec<TPrecCoeffG> getRotated (const SxDiracVec<TPrecCoeffG> &proj,
                                          const Coord &kVec, int iSym) const;
   public:
      void test (int argc, char **argv);
};

#endif /* _SX_PAW_EXCHANGE_H_ */
