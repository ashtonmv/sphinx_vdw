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

#ifndef _SX_RADIAL_BASIS_H_
#define _SX_RADIAL_BASIS_H_

#include <SxPrecision.h>
#include <SxQuantumNumbers.h>
#include <SxBasis.h>
#include <SxGBasis.h>
#include <SxRadGBasis.h>
#include <SxDirac.h>
#include <SxArray.h>
#include <SxTimer.h>
#include <SxDFT.h>

/** This class describes a radial mesh basis. It will be used e.g. for
    atomic wavefunctions and in PAW part. In contrast to most of the
    other basis modules it contains a intermediate function to convert
    a wavefunction given on a radial grid to a plane wave grid.

    \b SxRadBasis = S/PHI/nX Radial Basis

    \ingroup group_dft
    \author  Sixten Boeck, Christoph Freysoldt
    */
class SX_EXPORT_DFT SxRadBasis : public SxBasis
{
   public:

      typedef TRadBasisType               TBasisType;
      typedef SxDiracVec<TRadBasisType>   TPsi;

      SxRadBasis ();
      SxRadBasis (const SxArray<SxDiracVec<TReal8> > &radFunc,
                  const SxArray<Real8>               &logDr);
      SxRadBasis (double rMin, double rMax, int nPoints);
      SxRadBasis (const SxVector<Double> &rMin, const SxVector<Double> &rMax, 
            const SxVector<Int> &nPoints);
      SxRadBasis (const SxString &file);
      SxRadBasis (const SxBinIO &io);
      SxRadBasis (const SxSymbolTable *table);
      void set   (const SxArray<SxDiracVec<TReal8> > &radFunc,
                  const SxArray<Real8>               &logDr);
      void set   (double rMin, double rMax, int nPoints);
      void set (const SxVector<Double> &rMin, const SxVector<Double> &rMax, 
            const SxVector<Int> &nPoints);
      void setup (const SxSymbolTable *table);

      virtual ~SxRadBasis ();
      void read (const SxString &file);
      void read (const SxBinIO &io);
      SxRadBasis readMesh (const SxString &file);
      SxRadBasis readMesh (const SxBinIO &io);

      /** \brief Add radial mesh
        \param radFuncIn radial functions for new species for different l
        \param logDrIn   \f$ \ln \frac{r(i+1)}{r(i)}\f$ for different l
        \return new species id

        The different radial meshes in Dirac projections are given by the
        handle->auxData.is information of the Dirac vector. The return
        value of this function gives you the 'is' to access the radial mesh
        that is added by this function.

        Note that adding meshes in this way is less efficient than
        giving all in the constructor, and that all internal cashes
        are cleaned. Avoid calling this routine very often.

        \todo Single-shot routine for toPWBasis ()
        */
      int addMesh (const SxDiracVec<TReal8> &radFuncIn,
                         double             logDrIn);

      /// Dummy function, since number of elements may vary
      virtual ssize_t getNElements () const {
         SX_EXIT;
         return -1;
      }

   protected:
      void registerMemoryObservers ();

   public:
      int getNSpecies () { return int(radFunc.getSize()); }
      /// Return basis type
      virtual SxString getType () const { return "|r>"; }

      REGISTER_PROJECTOR (SxRadBasis, SxRadBasis,  changeRadBasis);
      REGISTER_PROJECTOR (SxRadBasis, SxGBasis,    toPWBasis);
      REGISTER_PROJECTOR (SxRadBasis, SxRadRBasis, toRadRBasis);
      REGISTER_PROJECTOR (SxRadBasis, SxRadGBasis, toRadGBasis);

      SxDiracVec<TBasisType> identity  (const SxRadBasis *,
                                        const SxDiracVec<TBasisType> &) const;
      SxDiracVec<TBasisType> changeRadBasis (const SxRadBasis *basis,
                                             const SxDiracVec<TBasisType> &vec) const;
      /**
        \brief Project to plane-wave basis

        \Note is the iAtom field in the input vector set, i.e., other
                 than -1, the translation operator is applied,
                 otherwise not...
       */
      SxDiracVec<TGBasisType> toPWBasis (const SxGBasis *,
                                         const SxDiracVec<TBasisType> &) const;
      SxDiracVec<TRadGBasisType> toRadGBasis (const SxRadGBasis *,
                                              const SxDiracVec<TBasisType> &) const;
      SxDiracVec<TRadGBasisType> toRadRBasis (const SxRadRBasis *,
                                              const SxDiracVec<TBasisType> &) const;

      /// Integrate over radial space with r^2 dr
      virtual Real8 tr (const SxDiracVec<Double> &) const;

      /** \brief Bessel functions
        */
      static  SxDiracVec<TReal8> jsb (int l, const SxDiracVec<TReal8> &);

      Real8 cosTheta (int ig, int jg, const SxGBasis &g);
      Real8 pl (int l, Real8 x);

      double integrate (const SxDiracVec<TReal8> &) const;

      /// Radial mesh, iSpecies:l:ir
      SxArray<SxDiracVec<TReal8> >  radFunc;  // :is,:r
      SxArray<Real8>                logDr;    // :is
      double getRMax (int is) const { return radFunc(is)(radFunc(is).getSize() - 1);};
      double getRMax () const;
      double getRMin (int is) const { return radFunc(is)(0);};
      ssize_t getNPoints (int is) const { return radFunc(is).getSize();};

      int getNSpecies () const {
         SX_CHECK (logDr.getSize () == radFunc.getSize (),
                   logDr.getSize (), radFunc.getSize ());
         return int(radFunc.getSize ());
      }

      /// \name Cashing
   protected:
      //@{
      /** The most expensive part in a <G+k|RY> projection are the integrals
        \f[ \int_0^\infty r^2 dr j_{l}(|G+k|\cdot r) \cdot R(r) \f]
        because they have to be calculated for each G vector (well, at least
        for each possible length |G+k|). Because this doesn't depend on m, we
        should cache these integrals for the current l and vector and all
        known G+k bases.
        \brief Cashed spherical Bessel integrals ik:ig
        */
      mutable SxArray<SxDiracVec<TRadBasisType> > cashedRl;
      /// The vector the integrals of which are cashed
      mutable SxDiracVec<Double> cashedVec;
      /// The known G bases, :ik:
      // no SxArray, because it's just for finding ik, where lists are
      // perfect
      mutable SxList<const SxGBasis *> gBases;

      //@}

   public:

      enum { S=0, IPY, IPZ, IPX, DXY, DYZ, DZ2, DXZ, DX2_Y2 };

      /** \brief Get real, unnormalized Ylm for all l<=lmax and all m
           for a shifted G basis
           \todo Move to SxGBasis (as member)
        */
      static SxDiracMat<Double> realYlm (int lmax, const SxGBasis &G,
                                         const Coord &dG);

      /**
        \brief Calculates the projection of the radial part onto plane waves

        This routine should not be used from the outside. Set up a proper
        radial basis containing all the radial meshes you need and use the
        Dirac projector notation.

        \Note The normalization factor of \f$\frac{4\pi}{\sqrt\Omega}\f$
              is NOT included here, because it is combined with the Ylm
              normalization factors later on.
        \todo Shouldn't be called from outside, so make it protected
        */
      SxDiracVec<TRadBasisType>
      toPWBasis (const SxDiracVec<TReal8> &rad,
                       const SxDiracVec<TReal8> &psi,
                       const SxGBasis &pwBasis,
                       int l, Real8 logDr) const;

      /**
        \brief cutoffvalue for (Rad|G)(G|Psi) Transformation
        */
      mutable double cutoff;

   public:
      /** \brief Specific radial basis for projections

          SxRadBasis is a container class for the radial meshes of all
          species. Sometimes you'd like to pick a particular mesh.
          This is what special basis is doing for you - it combines the
          multi-species container with meta-data (species, l, m).

          Don't treat this as a real basis.
        */
      class SpecialBasis
      {
         friend class SxRadBasis;
         /// The full radial basis
         const SxRadBasis &radBasis;
         public:
            /// Specific data for radial basis
            int is, l, m;
         protected:
            /// Constructor
            inline
            SpecialBasis (const SxRadBasis &fullBasis, int isIn,
                          int lIn = -1, int mIn = 0)
            : radBasis(fullBasis), is(isIn), l(lIn), m(mIn)
            {
               SX_CHECK (is >= 0 && is < radBasis.radFunc.getSize (),
                         is, radBasis.radFunc.getSize ());
               SX_CHECK (l < 0 || abs(m) <= l, l, m);
            }
         public:
            /// Get number of elements
            int getSize () const  {
               return (int)radBasis.radFunc(is).getSize ();
            }

            /// SxDiracVec constructor
            inline
            operator SxDiracVec<Double> () const
            {
               SxDiracVec<Double> res(getSize ());
               res.setBasis (&radBasis);
               res.handle->auxData.is = is;
               res.handle->auxData.l  = l;
               res.handle->auxData.m  = m;
               return res;
            }

      };
      /// Return species-specific radial basis
      inline SpecialBasis operator() (int is) const {
         SX_CHECK(is < radFunc.getSize (), is, radFunc.getSize ());
         return SpecialBasis (*this, is);
      }
};


namespace Timer {
   enum RadBasisTimer {
      JsbInt,
      JsbCalc,
      Phase,
      RadTotal,
      rad2radG,
      rad2radR
   };
}

SX_REGISTER_TIMERS (Timer::RadBasisTimer)
{
   using namespace Timer;
   regTimer (JsbInt,   "jsb Integrals");
   regTimer (JsbCalc,  "jsb Setup");
   regTimer (RadTotal, "radBasis");
   regTimer (rad2radG, "rad to radG Basis");
   regTimer (rad2radR, "rad to radR Basis");
}

/** \brief Stand alone radial interpolation function
     @param psi function to be interpolated
     @param r0     start of original logarithmic grid
     @param logDr  logDr of original logarithmic grid
     @param newRad new grid
     @return interpolated function

     Interpolate psi from a logarithmic grid (defined by r0 and logDr)
     to new grid given in newRad.
     \author C. Freysoldt
  */
SxDiracVec<Double> SX_EXPORT_DFT
interpolateRad(const SxDiracVec<Double> &psi,
               double r0, double logDr,
               const SxDiracVec<Double> &newRad);
#endif /* _SX_RADIAL_BASIS_H_ */
