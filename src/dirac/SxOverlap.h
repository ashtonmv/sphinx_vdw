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

#ifndef _SX_OVERLAP_H_
#define _SX_OVERLAP_H_

#include <SxDirac.h>
#include <SxPtr.h>

/** Determines what method in orthogonalization/orthonormalization 
    will be applied.
 */
enum SxOrthoMethod { GramSchmidt, Loewdin};

/** \brief Generic overlap operator interface

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DIRAC SxOverlapBase
{
   public:
      /// Virtual destructor
      virtual ~SxOverlapBase () {}

      /// Norm square
      virtual double normSqr (const SxDiracVec<Complex16> &psi) const = 0;

      /// Scalar product
      virtual SxComplex16 dot (const SxDiracVec<Complex16> &x,
                               const SxDiracVec<Complex16> &y) const = 0;

      /// Apply overlap operator
      virtual SxDiracVec<Complex16> 
      apply (const SxDiracVec<Complex16> &psi) const = 0;
      
      /// Set states x orthogonal to states y
      virtual void setOrthogonal (SxDiracVec<Complex16> *xPtr,
                                  const SxDiracVec<Complex16> &y) const = 0;

      /// Orthonormalize states x
      virtual void orthonormalize (SxDiracVec<Complex16> *xPtr,
                                   const SxOrthoMethod how = GramSchmidt) const = 0;
             
      /// Normalize states x
      virtual void normalize (SxDiracVec<Complex16> *xPtr) const;
                               
      /// Get overlap matrix
      virtual SxDiracMat<Complex16>
      getMatrix (const SxDiracVec<Complex16> &x,
                 const SxDiracVec<Complex16> &y) const
      {
         SX_CHECK (x.getSize () > 0);
         SX_CHECK (y.getSize () > 0);
         SX_CHECK (x.getBasisPtr () == y.getBasisPtr ());
         if (x.nCols () > y.nCols ())
            return x.overlap (apply (y));
         else
            return apply (x).overlap (y);
      }

};
/** \brief Generic overlap operator

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DIRAC SxOverlap
{
   protected:
      SxConstPtr<SxOverlapBase> S;
   public:
      SxOverlap (const SxConstPtr<SxOverlapBase> &in) 
         : S(in)
      { /* empty */ }

      /// Set states x orthogonal to states y
      void setOrthogonal (SxDiracVec<Complex16> *xPtr,
                          const SxDiracVec<Complex16> &y) const
      {
         S->setOrthogonal (xPtr, y);
      }

      /// Orthonormalize states x
      void orthonormalize (SxDiracVec<Complex16> *xPtr,
                   const SxOrthoMethod how = GramSchmidt) const
      {
         S->orthonormalize (xPtr, how);
      }
             
      /// Normalize states x
      void normalize (SxDiracVec<Complex16> *xPtr) const
      {
         S->normalize (xPtr);
      }

      SxDiracMat<Complex16> 
      getMatrix (const SxDiracVec<Complex16> &x,
                 const SxDiracVec<Complex16> &y) const
      { 
         return S->getMatrix (x, y);
      }
                               
   protected:
      /// Overlap operator expression
      class SxOvlpExpr  {
         protected:
            /// Type of expression
            const enum { SPsi, PsiS } type;

            /// Overlap operator reference
            const SxOverlapBase &S;

            /// psi reference
            const SxDiracVec<Complex16> &psi;
         public:
            SxOvlpExpr (const SxOverlapBase &SIn,
                        const SxDiracVec<Complex16> &psiIn)
               : type(SPsi), S(SIn), psi(psiIn)
            { /* empty */ }

            SxOvlpExpr (const SxDiracVec<Complex16> &psiIn,
                        const SxOverlapBase &SIn)
               : type(PsiS), S(SIn), psi(psiIn)
            { /* empty */ }

            /// Perform S | psi
            operator SxDiracVec<Complex16> () const
            {
               SX_CHECK (type == SPsi);
               return S.apply (psi);
            }

            /// Perform (psi | S | psi2)
            SxComplex16 operator| (const SxDiracVec<Complex16> &psi2)
            {
               SX_CHECK (type == PsiS);
               if (&psi == &psi2) 
                  return S.normSqr (psi);
               else
                  return S.dot (psi, psi2);
            }
      };

   public:
      /// Get S|psi expression
      SxOvlpExpr operator| (const SxDiracVec<Complex16> &psi)
      {
         return SxOvlpExpr (*S, psi);
      }

      /// Apply S
      SxDiracVec<Complex16> operator* (const SxDiracVec<Complex16> &psi) const
      {
         return S->apply (psi);
      }
      
      friend SxOvlpExpr operator| (const SxDiracVec<Complex16> &,
                                   const SxOverlap &S);
};

/// Get psi|S expression
inline 
SxOverlap::SxOvlpExpr operator| (const SxDiracVec<Complex16> &psi,
                                 const SxOverlap &S)
{
   return SxOverlap::SxOvlpExpr (psi, *S.S);
}

#endif /* _SX_OVERLAP_H_ */
