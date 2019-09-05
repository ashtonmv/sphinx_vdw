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

#ifndef _SX_RAD_MAT_H_
#define _SX_RAD_MAT_H_
#include <SxDFT.h>
#include <SxMatrix.h>
#include <SxAtomicStructure.h>
#include <SxYlm.h>
#include <SxBlockDensityMatrix.h>

class SxPAWPot;
//template<class T> class SxDiracVec<T>;

/** \brief PAW 1-center matrices

    \b SxClass = S/PHI/nX PAW 1-center matrices


    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxRadMat : public SxBlockDensityMatrix
{
   protected:
      SxConstPtr<SxAtomInfo> atomInfo;
      friend class SxPAWRho;
      friend class SxPreconditioner;
   public:

      /// Assignment operator
      void operator= (const SxRadMat &in)
      {
         if (&in == this) return;
         atomInfo = in.atomInfo;
         rho = in.rho;
      }

      /// Get a 1-center matrix
      SxMatrix<Double> &operator() (int iSpin, int is, int ia)  {
         return rho(iSpin, atomInfo->offset (is) + ia);
      }

      /// Get a 1-center matrix
      const SxMatrix<Double> &operator() (int iSpin, int is, int ia) const  {
         return rho(iSpin, atomInfo->offset (is) + ia);
      }

      // autoloop wrapper
      template<class T1, class T2, class T3>
      SxMatrix<Double> &operator() (const T1 &iSpin,
                                    const T2 &is,
                                    const T3 &ia)
      {
         SxAutoLoop::setLimit(iSpin, getNSpin ());
         SxAutoLoop::setLimit(is, atomInfo->nSpecies);
         SxAutoLoop::setLimit(ia, atomInfo->nAtoms ((int)is));
         return rho((int)iSpin, atomInfo->offset((int)is) + ia);
      }

      // autoloop wrapper
      template<class T1, class T2, class T3>
      const SxMatrix<Double> &operator() (const T1 &iSpin,
                                          const T2 &is,
                                          const T3 &ia) const
      {
         SxAutoLoop::setLimit(iSpin, getNSpin ());
         SxAutoLoop::setLimit(is, atomInfo->nSpecies);
         SxAutoLoop::setLimit(ia, atomInfo->nAtoms ((int)is));
         return rho((int)iSpin, atomInfo->offset((int)is) + ia);
      }

      /// Resize to nSpin and number of atoms, don't initialize matrices
      void resize (const SxConstPtr<SxAtomInfo> &info,
                   int nSpin);
      /// Resize
      void resize (const SxConstPtr<SxAtomInfo> &info,
                   const SxPAWPot &pawPot,
                   int nSpin);
      /// Resize according to example
      void resize (const SxRadMat &in);

      void set (double x);

      /*
      /// In-place scaling
      void operator *= (double s)
      {
         for (int iSpin = 0; iSpin < getNSpin (); ++iSpin)
            for (int iTl = 0; iTl < getNSite (); ++iTl)
               rho(iSpin, iTl) *= s;
      }
      */

      void symmetrize (const SxAtomicStructure &str,
                       const SxPAWPot &pot,
                       const SxYlmRotGroup &ylmRot);

      inline int getNSpecies () const
      {
         return atomInfo->nSpecies;
      }

      inline int getNAtoms (int is) const
      {
         return atomInfo->nAtoms (is);
      }

      /** \brief Multiply matrices with projections
        @param p   \f$\langle p_i|\psi\rangle\f$

        \f[
        \sum_j \mathrm i^{l_i-l_j} A_{ij} <p_j | \psi>
        \f]

        Aij are the matrices stored in rho.
        */
      SxDiracVec<Complex16> operator^ (const SxDiracVec<Complex16> &p) const;

      // --- overloaded functions from SxDensity
      /// Write to binary file
      virtual void readRho (const SxBinIO &io);
      /// Read from binary file
      virtual void writeRho (SxBinIO &io) const;
};

/** \brief Compute trace
  \f[
  \sum_{R} \sum_{i,j\,@R} D_{ij} A_{ij}
  \f]
  */
double tr (const SxRadMat &D, const SxRadMat &A);

#endif /* _SX_RAD_MAT_H_ */
