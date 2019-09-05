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

#ifndef _SX_FORCE_SYM_H_
#define _SX_FORCE_SYM_H_

#include <SxGeom.h>
#include <SxOperator.h>
#include <SxAtomicStructure.h>

/** \brief Symmetrize forces (or similar)

    \b SxForceSym = S/PHI/nX force (or similar) symmetrizer

    This class provides a symmetrization filter for forces.
    Applying a symmetry to the atoms corresponds to a permutation
    of atom indices:
    \f[
    S \tau_i \rightarrow \tau_{P_S(i)} + \mathbf R
    \f]
    where \f$\mathbf R\f$ is a lattice vector.

    Symmetrizing forces (or displacements, or similar) is defined
    as
    \f[
    \mathbf f_i^{\rm sym} = \sum_S  S \mathbf f_{P_S^{-1}(i)}
    \f]

    This becomes clear when you imagine
    the forces (or atom shifts) as being attached to the
    corresponding atoms (as little vectors). When you rotate the
    structure, you need to know not only how the force rotates, but
    also where the attachment point moves.


    We use SxAtomInfo to store and apply the permutations.

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_GEOM SxForceSym
	: public SxOperator<SxAtomicStructure> SXOP_LINKFIX
{
   protected:
      /// The permutations for each symmetry
      SxArray<SxAtomInfo::ConstPtr> permutations;

   public:
      /// Constructor
      SxForceSym (const SxAtomicStructure &str) { setup (str); }

      /// Empty constructor
      SxForceSym () {/* empty */}

      void setup (const SxAtomicStructure &str);

      /// How to apply the symmetrizer
      virtual SxAtomicStructure
      operator*(const SxAtomicStructure &forces) const;

      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT; }

      /** \brief Get SxPtr copy
          \note Standard copy constructor is OK, no dangerous members
      */
      SXOPERATOR_GETCOPY(SxForceSym,SxAtomicStructure);

      /// Get number of symmetries
      inline int getNSymmetries () const {
         return int(permutations.getSize ());
      }

      bool checkStr (const SxAtomicStructure &str) const;

};

#endif /* _SX_FORCE_SYM_H_ */
