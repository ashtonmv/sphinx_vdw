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

#ifndef _SX_STICKY_FILTER_H_
#define _SX_STICKY_FILTER_H_

#include <SxOperator.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>
#include <SxStruct.h>

/** \brief Apply constraints to forces

    \b SxStickyFilter = SPHInX Sticky Atoms Filter

    This filter is used to eliminate all coordinates belonging to atoms
    marked as \em sticky or \em not \em movable. Usually it will be applied 
    to forces in order to make sure that only \e movable atoms are treated
    in structure optimization routines.

    The Filter \f$\phi\f$ is 1 for all movable atoms and 0 for all 
    sticky atoms, i.e.,
    \f[
       \phi_{i_s,i_a} 
          = \left\{ 
                \begin{array}{cl}
                   1 & (i_s,i_a) \quad \mathrm{movable} \\
                   0 & \mathrm{sticky}
                \end{array}
            \right.
    \f]
    The application of the filter is
    \f[
       \mathbf{F}'_{i_s,i_a}  = \phi_{i_s,i_a} \mathbf{F}_{i_s,i_a}
    \f]
    \f[
       \mathbf{F}_{i_s,i_a}' 
           = \left\{ 
               \begin{array}{ll}
                  F_{i_s,i_a} & (i_s,i_a) \quad \mathrm{movable} \\
                  0           &  \mathrm{otherwise}
               \end{array}
             \right.
    \f]


    \ingroup group_structure;
    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_STRUCT SxStickyFilter 
: public SxOperatorBase<SxAtomicStructure>, 
  public SxOperatorBase<SxVector<TPrecTauR> >
{
   public:

      SxStickyFilter ();
      SxStickyFilter (const SxSymbolTable *);
      virtual ~SxStickyFilter ();

      /** \brief validate that symmetry equivalent atoms have same stickyness
          \param S             symmetry matrix
          \param equivalentIdx indices of equivalent atoms by symmetry */
      void validate (const SymMat &S, const SxVector<Int> &equivalentIdx);
      
      // enable types explicitly because this is a multi-type operator
      SXOPERATOR_TYPE (SxVector<TPrecTauR>);
      SXOPERATOR_TYPE (SxAtomicStructure);

      virtual
      SxVector<TPrecTauR> operator* (const SxVector<TPrecTauR> &) const;
      virtual void applyInPlace (SxVector<TPrecTauR> &) const
      {
         SX_EXIT; 
      }
      virtual SxAtomicStructure operator* (const SxAtomicStructure &) const;
      virtual void applyInPlace (SxAtomicStructure &) const { SX_EXIT; }

   protected:
      SxArray<SxVector3<Int> > sticky;  // :ia:iDoF

      
      SXOPERATOR_GETCOPY (SxStickyFilter, SxVector<TPrecTauR>);
      SXOPERATOR_GETCOPY (SxStickyFilter, SxAtomicStructure);
   public:
      /// Read access to the underlying sticky data
      const SxArray<SxVector3<Int> >& getStickyArray () const
      { 
         return sticky; 
      } 

      /// Check whether the filter is an expensive 1
      bool isUnity () const {
         SX_LOOP2(i,d) if (sticky(i)(d) != 1) return false;
         return true;
      }
};

#endif /* _SX_STICKY_FILTER_H_ */
