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

#ifndef _SX_PSI_SET_H_
#define _SX_PSI_SET_H_

#include <SxMemConsumer.h>
#include <SxDFT.h>

/** \brief Abstract class to represent any set of wave function \f$ 
              | \Psi \rangle
           \f$

  \sa      \ref page_dirac
  \ingroup group_dirac
  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxPsiSet : public SxMemConsumer
{
   public:

      /** \brief Possible wavefunction types */
      enum PsiType { NONE, PW, ATOMORB, WANNIER };

      SxPsiSet (PsiType t) : SxMemConsumer ()  { type = t; }
      virtual ~SxPsiSet ()                     { }

      /** \brief Returns the wavefunction type */
      PsiType getType () const { return type; }

   protected:

      /** \brief The type of the wavefunction. */
      PsiType type;

      SxPsiSet ()  { type = NONE; }

};

#endif /* _SX_PSI_SET_H_ */
