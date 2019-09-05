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

#ifndef _SX_PAW_BASIS_H_
#define _SX_PAW_BASIS_H_

#include <SxDFT.h>
#include <SxGBasis.h>
#include <SxPartialWaveBasis.h>

/** \brief Full PAW basis

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWBasis
: public SxBasis
{
   public:
      /// Pointer to G+k basis
      SxConstPtr<SxGBasis> gBasis;

      /// Pointer to partial wave basis
      SxConstPtr<SxPartialWaveBasis> pBasis;

      /// Constructor
      SxPAWBasis (const SxConstPtr<SxGBasis> &gIn,
                  const SxConstPtr<SxPartialWaveBasis> &pIn)
         : gBasis (gIn),
           pBasis (pIn)
      {
         SX_CHECK (gIn);
         SX_CHECK (pIn);
      }  

      // --- SxBasis stuff
      typedef Complex16 TBasisType;

      virtual ~SxPAWBasis () {
         // empty
      }

      // --- overloaded virtual SxBasis functions ----------
      virtual ssize_t getNElements () const
      {
         return gBasis->getNElements () + pBasis->getNElements ();
      }

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const  {
         return "|G+k,p>";
      }

      virtual 
      SxComplex16 scalarProduct (const SxDiracVec<Complex16> &x,
                                 const SxDiracVec<Complex16> &y) const;

      REGISTER_PROJECTOR (SxPAWBasis, SxPAWBasis, identity);
      REGISTER_PROJECTOR (SxPAWBasis, SxGBasis, toPWBasis);
      REGISTER_PROJECTOR (SxPAWBasis, SxRBasis, toRBasis);
      REGISTER_PROJECTOR (SxPAWBasis, SxAOBasis, toAO);
      REGISTER_PROJECTOR (SxPAWBasis, SxPartialWaveBasis, toPartials);
      REGISTER_PROJECTOR (SxPAWBasis, SxBasis, toAny);

      /** \brief Project to G-basis */
      SxDiracVec<TGBasisType> toPWBasis (const SxGBasis *,
                                         const SxDiracVec<Complex16> &) const;
      /** \brief Project to p-basis */
      SxDiracVec<TGBasisType> toPartials (const SxPartialWaveBasis *,
                                          const SxDiracVec<Complex16> &) const;
      /** \brief Project to R-basis */
      SxDiracVec<TGBasisType> toRBasis (const SxRBasis *,
                                        const SxDiracVec<Complex16> &) const;
      /** \brief Project to R-basis */
      SxDiracVec<TGBasisType> toAO     (const SxAOBasis *,
                                        const SxDiracVec<Complex16> &) const;
      /// \brief Identity projection
      SxDiracVec<TGBasisType> identity (const SxPAWBasis *basis,
                                        const SxDiracVec<Complex16> &in) const
      {
         SX_CHECK (in.getBasisPtr () == this);
         SX_CHECK (basis == this);
         return in;
      }
      /** \brief Projector from \b PAW space to some anonymous basis

          This is the versatile interface: we try to dynamically determine
          the basis to project to.

          \sa \ref page_dirac */
      SxDiracVec<Complex16> toAny (const SxBasis *,
                                   const SxDiracVec<Complex16> &) const;

};

#endif /* _SX_PAW_BASIS_H_ */
