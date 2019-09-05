// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_PROJ_PAIR_MATRIX_H_
#define _SX_PROJ_PAIR_MATRIX_H_

#include <SxList.h>
#include <SxVector.h>
#include <SxMath.h>

/**  \brief A projector pair matrix

  While a projector matrix is expressed as
  \f[
  \sum_k \vec p_k^T a_k \vec p_k
  \f]
  we have here the same with projector pairs \f$pk = (s_k t_k)\f$, where
  the projector weight \f$a_k\f$ transforms into a 2x2 matrix
  \f[
  \left(
  \begin{array}{cc} a_{ss} & a_{st} \\ 
                    a_{ts} & a_{tt}
  \end{array}
  \right)_k
  \f]

  \author C. Freysoldt freysoldt@mpie.de

  */
class SX_EXPORT_MATH SxProjPairMatrix {
   public:
      typedef SxVector<TPrecTauR> Vector;
      /// Initial diagonal part
      double diag;
      /// Number of projector pairs
      int nProj;
      /// The projector pairs
      SxList<Vector> sk,tk;
      /// The projector mixing matrix
      SxList<double> akss,akst,akts,aktt;
      
      SxProjPairMatrix (double diagIn, int nProjIn)
         : diag(diagIn), nProj(nProjIn) 
      {
         SX_CHECK (diag > 0., diagIn);
         SX_CHECK (nProj > 0, nProj);
      }

      /// Add a new projector pair
      void append (const Vector& s, const Vector& t,
                   double prefactor,
                   double ass, double ast, double ats, double att);

      /// Apply matrix to a vector
      Vector operator^(const Vector& in) const;
};

#endif /* _SX_PROJ_PAIR_MATRIX_H_ */
