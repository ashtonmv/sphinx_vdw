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

#ifndef _SX_MESH3D_H_
#define _SX_MESH3D_H_
#include <SxVector3.h>

/** \brief Simple 3D mesh class

    \b SxClass = S/PHI/nX 3D mesh

    The main purpose of this class is to map between linear (storage)
    and 3D indices.

    \author Christoph Freysoldt, freysoldt@mpie.de
 */ 
class SxMesh3D : public SxVector3<Int>
{
   public:
      /// Possible index vector ranges
      enum VecRange {
         /// each vector coordinate 0 <= n(i) < mesh(i)
         Positive,
         /// each vector coordinate -mesh(i)/2 < n(i) <= mesh(i)/2
         Origin,
         /// no restrictions
         Unknown
      };

      /// Standard constructor
      SxMesh3D (int l = 0, int m = 0, int n = 0) 
         : SxVector3<Int> (l, m, n) { /* empty */}

      /// Constructor from vector
      SxMesh3D (const SxVector3<Int> &mesh) : SxVector3<Int>(mesh) { }

      /// Return size
      inline ssize_t getSize () const { return product (); }
      
      /// Get linear mesh index from mesh vector (c running fastest)
      inline ssize_t getMeshIdx (int a, int b, int c, VecRange range) const
      {
         SX_CHECK (getSize () > 0);
         const SxVector3<Int> &m = *this;
         SX_CHECK (m(0) != 0 && m(1) != 0 && m(2) != 0,
                   m(0), m(1), m(2));
         switch (range)  {
            case (Unknown)  : if (abs(a) >= m(0)) a %= m(0);
                              if (abs(b) >= m(1)) b %= m(1);
                              if (abs(c) >= m(2)) c %= m(2);
                              // fall through
            case (Origin)   : if (a < 0) a += (*this)(0);
                              if (b < 0) b += (*this)(1);
                              if (c < 0) c += (*this)(2);
                              break;
            case (Positive) : break;
            default         : SX_EXIT;
         } 
         SX_CHECK (a>=0 && a < m(0), a, m(0));
         SX_CHECK (b>=0 && b < m(1), b, m(1));
         SX_CHECK (c>=0 && c < m(2), c, m(2));
         return c + m(2) * (b + m(1) * a);
      }

      /// Get linear mesh index from mesh vector
      inline ssize_t getMeshIdx (const SxVector3<Int> &n, VecRange range) const
      {
         return getMeshIdx(n(0), n(1), n(2), range);
      }
         
      /// Turn linear mesh index into vector
      inline SxVector3<Int> getMeshVec(ssize_t idx, VecRange range) const
      {
         SX_CHECK (getSize () > 0);
         SX_CHECK(idx >= 0 && idx < getSize (), idx, getSize ());
         SX_CHECK(range == Positive || range == Origin);
         const SxVector3<Int> &m = *this;
         SxVector3<Int> res;
         ssize_t a = idx / m(2);
         res(0) = int(a / m(1));
         // res(1) = a % m(1);
         res(1) = int(a - res(0) * m(1));
         // v(2) = idx % m(2)
         res(2) = int(idx - a * m(2));
         if (range == Origin)  {
            if (2 * res(0) > m(0)) res(0) -= m(0);
            if (2 * res(1) > m(1)) res(1) -= m(1);
            if (2 * res(2) > m(2)) res(2) -= m(2);
         }
         return res;
      }
};

inline std::ostream& operator<< (std::ostream &out, const SxMesh3D &mesh)
{
   out << mesh(0) << " x " << mesh(1) << " x " << mesh(2);
   return out;
}

#endif /* _SX_MESH3D_H_ */
