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

#ifndef _SX_GAUSS_IJ_H_
#define _SX_GAUSS_IJ_H_

#include <SxMath.h>
#include <SxMatrix.h>
#include <SxYlm.h>

/** \brief Analytic integration of Gaussian 2-center integrals

    \note Currently, only the Hartree energy is implemented.
    Overlap matrix elements would be obtained by using n+2 instead of n from
    the radial Kernel (computeK).
    Kinetic energy matrix elements would be obtained by using n+4 instead of
    n in the radial kernel.

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_MATH SxGaussIJ
{
   public:
      // precomputed double factorial (2l+1)!! = 1 * 3 * 5 * ... (2l+1)
      SxArray<double> doubleFactorial;

      // Compute double factorial (2l+1)!! = 1 * 3 * 5 * ... (2l+1)
      void computeDoubleFactorial (int n);

   protected:
      /** \brief One-dimensional integral  K (PAW notes, Eq. 10) */
      SxMatrix<Double> K;

   public:
      /** \brief One-dimensional integral  K (PAW notes, Eq. 10)
        \f[
            K\left(n,l,\alpha^2\right) :=
            \int dk~ k^{n} j_{l}(k) e^{-\frac 14k^2/\alpha^2}
        \f]
        */
      void computeK (double a2, int lmax);

      /** \brief One-dimensional integral K with (r)^{-(n+1)} prefactor
        @param r2  \f$|\mathbf r|^2\f$
        @param rc2 \f$r_c^2\f$
        @param lMax maximum l and n value to compute

        The result is stored in K.
        */
      void computeK (double r2, double rc2, int lMax);

      /** \brief Set distance vector, Gaussian width, and maximum l
        @param r    distance vector \f$\tau_2 - \tau_1\f$. 1 will be the left
                    side (first index), 2 will be the right side (2nd index)
                    of the interaction matrix computed by computeWork.
        @param rc2  sum of the Gaussian widths
        @param lSum sum of the maximum l values
        
        This routine sets up auxiliary data for successive calls to
        #computeWork.
      */
      void set (const SxVector3<Double> &r, double rc2, int lSum);

      /** \brief Set distance vector, two Gaussian widths, and maximum l
        @param r    distance vector \f$\tau_2 - \tau_1\f$. 1 will be the left
                    side (first index), 2 will be the right side (2nd index)
                    of the interaction matrix computed by computeWork.
        @param rc2A  sum of the Gaussian widths
        @param rc2B  sum of the Gaussian widths to be subtracted
        @param lSum  sum of the maximum l values
        
        This routine sets up auxiliary data for successive calls to
        #computeWork, which will then compute E(rc2A) - E(rc2B).
        Using the delta routine is slightly more efficient than computing
        the difference from the final results, since the angular loops
        are independent of rc2 and hence executed only once.
      */
      void setDelta (const SxVector3<Double> &r,
                     double rc2A, double rc2B, int lSum);

      /** \brief Compute interaction matrix
          @param lmax1 maximum l value left side
          @param lmax2 maximum l value right side  
          @param clebschGordan Clebsch-Gordan coefficients for
                 with limits of (lmax1, lmax1 + lmax2 + order, lmax2)
          @param order  If 0, get the interaction energy. If 1, prepare
                        for forces. If 2, prepare for Hessian. Forces and
                        Hessian must then be computed using the corresponding
                        routines getForce and getHessian.
          @param resPtr If given, put results into this matrix. Otherwise,
                        the internal workspace will be used. The returned
                        matrix should never be assigned, but directly used.
          @return The interaction matrix.
          */
      SxMatrix<Double> compute (int lmax1, int lmax2, 
                                const SxYlm::SxClebschTable &clebschGordan,
                                int order = 0,
                                SxMatrix<Double> *resPtr = NULL);

      /** \brief Compute interaction matrix */
      void compute (int lmax1, int lmax2, 
                    const SxYlm::SxClebschTable &clebschGordan,
                    SxMatrix<Double> *resPtr)
      {
         SX_CHECK (resPtr);
         compute (lmax1, lmax2, clebschGordan, 0, resPtr);
      }

      /** \brief Calculate force for an interaction
        @param lm1 combined (l,m) index for left side
        @param lm2 combined (l,m) index for right side
        @param clebschGordan Clebsch-Gordan coefficients

        \note Before this routine is called, the internal workspace must
              be set up with computeWork (...., 1)
        */
      SxVector3<Double> getForce (int lm1, int lm2, 
                                  const SxYlm::SxClebschTable &clebschGordan);
      /** \brief Calculate Hessian for an interaction
        @param lm1 combined (l,m) index for left side
        @param lm2 combined (l,m) index for right side
        @param clebschGordan Clebsch-Gordan coefficients

        \note Before this routine is called, the internal workspace must
              be set up with computeWork (...., 2)
        */
      SxMatrix3<Double> getHesse (int lm1, int lm2, 
                                  const SxYlm::SxClebschTable &clebschGordan);
   protected:
      /// Distance vector squared
      double dist2;

      /// Y_lm (r)
      SxVector<Double> ylmR;

      /// Set up ylmR
      void setupYlm (const SxVector3<Double> &r, int lMax);

      /// work space (to avoid frequent memory allocation)
      SxVector<Double> workspace;

      /// What the workspace is prepared for (0=energy, 1=forces, 2=Hessian)
      int workOrder;

      /// Number of rows in workspace matrix
      int workNlm1;

      /// Number of columns in workspace matrix
      int workNlm2;

   public:
      /// Constructor
      SxGaussIJ ();

};

#endif /* _SX_GAUSS_IJ_H_ */
