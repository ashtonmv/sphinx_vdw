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

#ifndef _SX_WANNIER_H_
#define _SX_WANNIER_H_


#include <SxConfig.h>
#include <SxPW.h>
#include <SxString.h>
#include <SxGkBasis.h>
#include <SxGauss.h>
#include <SxKFinDiffs.h>
#include <SxDirTensor3.h>  // TODO: write an SxTensor4 class for including spin
#include <SxExt.h>


/** \brief Creates "maximally localized Wannier functions" from a set of Bloch
           states being represented in plane waves.

    \b SxWannier = S/PHI/nX Maximally Localized Wannier Functions

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de
            Hazem Abu-Farsakh, hazem@phys.upb.de   */

class SX_EXPORT_EXT SxWannier
{
   public:

      typedef SxArray<SxArray<SxVector3<Double> > >                WFcenters;
      typedef SxArray<double>                                      WFspreads;
      typedef SxArray<SxArray<SxArray<SxDiracMat<Complex16> > > >  WFTensor5;
      typedef SxArray<SxDirTensor3<Complex16> >                    WFTensor4;
      typedef SxArray<SxArray<SxDiracMat<Complex16> > >            WFgradients;
      typedef SxArray<SxArray<SxDiracVec<Complex16> > >            WFwannier;

      ///@name Constructors and Destructor
      //@{
      /** constructor */
      SxWannier ();

      /** constructor */
      SxWannier (const SxSymbolTable *table);

      /** destructor */
      ~SxWannier ()  { /* empty */ }
      //@}

      void print () const;

      /** create initial Wannier functions from Gauss orbitals */
      void initOrbitals ();

      /** create finite difference vectors in k space */
      void initKDiffs ();

      /** compute initial ovelap tensor elements
          \f$ M_{mn}^{(0)(\sigma)(\mathbf{k},\mathbf{b})} \f$, defined by
          Eq. (58) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>
       */
      void initOverlapMatrices ();

      /** write initial overlap matrices to output file */
      void writeMkb (const WFTensor5 &Mkb, const SxString &filename);

      /** find the transformation tensor \f$ U_{mn}^{\mathbf{k}} \f$,
          defined by the Eqs. (59) and (60) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>,
          which transformes the Bloch wavefunctions
          \f$ |u_{n\mathbf{k}} \rangle \f$ into the Maximally localized
          Wannier functions \f$ |w_{n\mathbf{R}} \rangle \f$. */
      void computeMLWFs ();

      /** compute the Wannier functions in real space */
      WFwannier computeMLWFsRealSpace ();

      /** write Wannier functions to sxb-file */
      void writeMLWFs ();

      /** logarithm of a complex number */
      // TODO get this into SxComplex.h
      SxComplex16 lnC (SxComplex16 &arg) const;

      /** exponential function of a complex number */
      // TODO get this into SxComplex.h
      SxComplex16 expC (const SxComplex16 &z) const;

      /** exponential function of anti-hermitian matrix */
      // TODO get this into Matrix class
      SxDiracMat<Complex16> mExp (const SxDiracMat<Complex16> &B) const;

      /** operator A[ B ] := (B - B ) / 2 */
      SxDiracMat<Complex16> sAnti (const SxDiracMat<Complex16> &B) const;

      /** operator S[ B ] := (B + B ) / 2i */
      SxDiracMat<Complex16> sSymm (const SxDiracMat<Complex16> &B) const;

      /** check matrix upon unitarity: $U \cdot U^{+} = E$ (!) */
      void uniCheck (const SxDiracMat<Complex16> &U) const;

      /** check matrix upon anti-hermitecity */
      void antiCheck (const SxDiracMat<Complex16> &W) const;

      /** dump a pair of 3-vectors in relative coordinates */
      void dumpPair (const Coord  &vec1,
                     const Coord  &vec2,
                     const double &r0);

      /** repeats a function periodically */
      // TODO: write an extra-add-on (maybe, as template)
      SxDiracVec<Complex16>
         repeatFunction (const SxDiracVec<Complex16> &function,
                         const RelVec &dim, const RelVec &repeat);

      /// for developing only
      void testScalarProduct ();
      void testKShift1 ();
      void testKShift2 ();
      void testComputeFunctions ();
      

   protected:
      /** file containing the underlying Bloch functions */
      SxString wavesfile;

      /** pointer to the |G+k> basises */
      SxGkBasis gkSet;  // TODO: rename into gkSet

      /** |G+0> basis */
      SxGBasis G;

      /** pointer to the unit cell */
      SxCell cell;

      /** pointer to the Bloch waves in PW representation */
      SxPW waves;

      /** pointer to the wavefunctions used as a starting point for the
          minimisation: could be the Bloch wavefunctions or their projections
          onto Gaussians */
      const SxPW *uInitPtr;

      /** whether or not to use initial guesses at the Wannier functions */
      bool useGuess;

      /** number of k-points */
      int nk;

      /** number of finite difference vectors */
      int nb;

      /** number of spin channels */
      int nSpin;

      /** maximum number of iteration steps */
      int maxSteps;

      /** convergence criterion for the minimisation */
      double dSpread;

      /** increment used in the steepest descent scheme */
      double alpha;

      /** band index of the lowest band considered */
      int iBottom;

      /** band index of the highest band considered */
      int iTop;

      /** number of bands considered */
      int nBands;

      /** number of meshpoints as 3-dim. vector */
      RelVec meshDim;

      /** cutoff given in the wavesfile */
      double eCut;

      /** repetition of Wannier functions in real space */
      SxVector3<Int> repetition;

      /** translation of Wannier functions in real space, just for
          the visualization */
      RelVec translation;

      /** the reciprocal unit cell, i.e. a matrix containing
          the \f$\mathbf{b}_i\f$-vectors */
      SxCell BCell;

      /** initial ovelap tensor \f$ M \f$, given through the elements
          \f$ M_{mn}^{(0)(\sigma)(\mathbf{k},\mathbf{b})} \f$, defined in
          Eq. (58) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>
       */
      WFTensor5 MkbInit;

      /** tranformation tensor \f$ U \f$, given by unitary trafo matrices
          with the elements \f$ U_{mn}^{(\sigma)(\mathbf{k})} \f$, defined
          via the Eqs. (59) and (60) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>
       */
      WFTensor4 Uk;

      /** transformation tensor \f$ P \f$ (preparation tensor), given by the
          products of the projection matrices of Eq. (62) and the Loewdin trafo
          \f$ \frac{1}{\sqrt{S}} \f$ in Eq. (63) of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>
       */
      SxDirTensor3<TPrecCoeffG> P;

      /** altogether transformation tensor \f$ D = P \cdot U \f$
          (\f$\mathbf{k}\f$-pointwise) turning the Bloch states into the set
          of localised Wannier functions
       */
      WFTensor4 D;

      /** compute centers of Wannier functions */
      WFcenters computeCenters (const WFTensor5 &Mkb) const;

      /** compute spread of Wannier functions */
      WFspreads computeSpreads (const WFTensor5 &Mkb,
                                const WFcenters &wrn) const;

      /** compute gradient of the spreads */
      WFgradients computeGradients (const WFTensor5 &Mkb,
                                    const WFcenters &wrn) const;

      void printCenters (const WFcenters &wrn) const;

   public:
      /** initial guesses at the Wannier functions */
      SxGauss initialGuess;

      /** finite difference vectors in k space */
      SxKFinDiffs kDiffs;
};


#endif /* _SX_WANNIER_H_ */
