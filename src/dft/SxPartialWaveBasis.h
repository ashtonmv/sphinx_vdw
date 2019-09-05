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

#ifndef _SX_PARTIAL_WAVE_BASIS_H_
#define _SX_PARTIAL_WAVE_BASIS_H_

#include <SxDFT.h>
#include <SxBasis.h>
#include <SxPAWPot.h>
#include <SxAOBasis.h>

/** \brief Partial wave basis (PAW method)

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPartialWaveBasis : public SxBasis
{
   public:
      typedef Complex16 TBasisType;
   protected:
      /// Pointer to PAW potentials (access via getPotPtr or getPot)
      SxConstPtr<SxPAWPot> potPtr;

      /// Number of elements
      ssize_t nElements;

   public:
      /// Get pointer to PAW potential
      SxConstPtr<SxPAWPot> getPotPtr () const
      {
         return potPtr;
      }

      /// Return reference to PAW potential
      const SxPAWPot& getPot () const
      {
         SX_CHECK (potPtr);
         return *potPtr;
      }

      /// Pointer to projector basis
      SxConstPtr<SxAOBasis> projectors;

      /// Constructor
      SxPartialWaveBasis (const SxConstPtr<SxPAWPot> &pot,
                          const SxAtomicStructure &str);

      virtual ~SxPartialWaveBasis ()
      {
         // empty
      }

      // --- overloaded virtual SxBasis functions ----------

      /// \brief Number of sampling points of the corresponding Dirac vector.
      virtual ssize_t getNElements () const
      {
         return nElements;
      }

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const  {
         return "|p>";
      }

      REGISTER_PROJECTOR (SxPartialWaveBasis, SxPartialWaveBasis, identity);
      REGISTER_PROJECTOR (SxPartialWaveBasis, SxGBasis, toPWBasis);

      /** \brief Project to G-basis */
      SxDiracVec<TGBasisType> toPWBasis (const SxGBasis *,
                                         const SxDiracVec<Complex16> &) const;
      /// \brief Identity projection
      SxDiracVec<TGBasisType> identity (const SxPartialWaveBasis *basis,
                                        const SxDiracVec<Complex16> &in) const
      {
         SX_CHECK (in.getBasisPtr () == this);
         SX_CHECK (basis == this);
         return in;
      }

      /** \brief Create a projector basis
        @param gk G+k basis
        @param pawPot partial wave basis
        */
      static SxPtr<SxAOBasis> createProjBasis (const SxGkBasis &gk,
                                               const SxPAWPot  &pawPot);

      /** \brief Create a projector basis
        @param gk G+k basis
        */
      inline void createProjBasis (const SxGkBasis &gk)
      {
         SX_CHECK (potPtr);
         projectors = createProjBasis (gk, *potPtr);
      }
};

#endif /* _SX_PARTIAL_WAVE_BASIS_H_ */
