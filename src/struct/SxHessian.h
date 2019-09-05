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

#ifndef _SX_HESSIAN_H_
#define _SX_HESSIAN_H_

#include <SxStruct.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>

/** \brief Hessian matrix

    \b SxClass = S/PHI/nX Hessian matrix

    \author C. Freysoldt */
class SX_EXPORT_STRUCT SxHessian
{
   protected:

      /// The Hessian matrix elements (iAtom:iNeighbor)
      SxArray<SxArray<SxMatrix3<Double> > > hessian;

      /// The sparsity
      SxArray<SxVector<Int> > neighbors;

      /// The cell mapping (optional)
      SxArray<SxArray<Int> > offset;

   public:
      /// Empty constructor
      SxHessian () { }
      /// Constructor from full matrix
      SxHessian (const SxMatrix<Double> &full);

      /// Construct full matrix
      SxMatrix<Double> getFull () const;

      /// Write full matrix to binary file
      static void writeFull (SxBinIO &io,
                             const SxMatrix<Double> &full);

      static void write (const SxString &fileName,
                         const SxMatrix<Double> &hessianFull,
                         const SxAtomicStructure &str);

      /// Read full matrix from binary file
      static SxMatrix<Double> readFull (const SxBinIO &io);

};

#endif /* _SX_HESSIAN_H_ */
