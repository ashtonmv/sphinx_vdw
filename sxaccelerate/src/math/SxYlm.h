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

#ifndef _SX_YLM_H
#define _SX_YLM_H
#include <SxMath.h>
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxMatrix3.h>
#include <SxNArray.h>

typedef SxArray<SxArray<SxMatrix<Double> > > SxYlmRotGroup;

class SX_EXPORT_MATH SxYlm  
{

   public:

      /// Type of spherical harmonics for Clebsch-Gordan coefficients
      enum YlmType { RealYlm, ComplexYlm };



      typedef SxNArray<double, 3> SxClebschTable;

      /** \brief Calculate unnormalized Plm for certain m and all l up to lmax 
       */
      static void getPlmArray(int lmax, int m, double x, SxVector<Double> *out);

      /// Get normalization factor for realYlm
      static double getYlmNormFactor (int l, int m);

      /// Get unnormalized, real Ylm
      static void getYlmArray(int lmax, double x, double y, double z,
                       SxVector<Double> *resultPtr);

      /// Get unnormalized, real Ylm
      inline static
      void getYlmArray(int lmax, const SxVector3<Double> &xyz,
                       SxVector<Double> *resultPtr)
      {
         getYlmArray(lmax, xyz(0), xyz(1), xyz(2), resultPtr);
      }

         /** \brief Combine l and m into one index

             The l's and m's are ordered according to
             -# increasing l
             -# increasing m

             This leads to the formula
             idx(l,m) = l*(l+1) + m
           */
      static inline int combineLm (int l, int m)  
      { 
         SX_CHECK (l >= 0, l);
         SX_CHECK (abs(m) <= l, m, l);
         return l*(l+1) + m; 
      }

      /** \brief Wigner 3j symbols

        These are calculated according to
        
        \f[ (-1)^{j1+j2+m3}\cdot\sqrt{\frac{
        (j1+m1)!(j2+m2)!(j3+m3)!(j1-m1)!(j2-m2)!(j3-m3)!}{
        (-j1+j2+j3)!(j1-j2+j3)!(j1+j2-j3)!(1+j1+j2+j3)!}}
        \cdot \sum_k (-1)^k \frac{
        (-j1+j2+j3)!(j1-j2+j3)!(j1+j2-j3)!}{
        (j3-j1-m2+k)!(j3-j2+m1+k)!(j1+j2-j3-k)!k!(j1-m1-k!)(j2-m2-k)!}
        \f]
        where the sum runs from \f$max(j2-j3-m1, j1-j3+m2, 0)\f$ to 
        \f$min(j1+j2-j3,j1-m1,j2-m2)\f$

        The derivation of this formula is left as an exercise to
        the curious reader ;-)

        */
      static double wigner3j (int j1, int j2, int j3, int m1, int m2, int m3);

      /**
        \brief Clebsch-Gordan coefficients for spherical harmonics
        \param l1max maximum l value for l1
        \param l2max maximum l value for l2
        \param l3max maximum l value for l3
        \param type (RealYlm/ComplexYlm)

        The Clebsch-Gordan coefficients are the coupling coefficients
        for direct products of basis functions, i.e.
        
        < X(l1,m1) * X(l2,m2) | X(l3,m3) >

        The Clebsch-Gordan coefficients by this routine include a factor
        \f$i^(l1+l2-l3)\f$ which takes into account that the radial space
        spherical harmonics have a factor \f$i^l\f$ so that the reciprocal
        space ones do not have it. 

        \return The coefficients as (l1,m1) : (l2,m2) : (l3,m3),
        where (l,m) stands for the combined index #combineLm (l,m)
        */
      static SxClebschTable getClebschGordan (int l1max, int l2max, int l3max, 
                        enum SxYlm::YlmType type = SxYlm::RealYlm);

      /** \brief Compute real Ylm representation of a cartesian rotation

          @note Usually, the complete rotation group is required. Use the
          corresponding function. This function should never be called 
          directly.

          @param lmax compute matrices up to l=lmax
          @param S cartesian rotation (must be unitary!)
          @param complexClebschGordan precomputed Clebsch-Gordan coefficients
                 for complex spherical harmonics (the coefficients are 
                 real, though)
          
          Formally, the representation can be written as
          \f[
          D_{MM'}^L = \int d^3 \Omega Y^*_{LM}(\Omega) Y_{LM'}(S\Omega)
          \f]

          For L=0, D is 1.

          For L=1, the cartesian representation is trivial to obtain,
          noting that Y_11 ~ x, Y_10 ~ z, Y_1-1 ~ y.

          For higher L, we switch to complex spherical harmonics and
          apply a recursion using Clebsch-Gordan coefficients
          \f[
          D_{MN}^L = \frac{1}{C^L} \sum_{mn} D_{(M-m)(N-n)}^{L-1} D_{mn}^1
                     \langle LM|(L-1)(M-m) 1m\rangle 
                     \langle (L-1)(N-n) 1n|LN\rangle
          \f]
          where C is given by
          \f[
          C^L = \frac{(2(L-1)+1)(2\cdot1 +1)}{4\pi} Wigner3j(L-1,1,L,0,0,0)^2 
          \f]
          and results from the expansion
          \f[
          C^L Y_{LM} = \sum_m Y_{(L-1)(M-m)}Y_{1m} 
                     \langle(L-1)(M-m) 1m|LM\rangle 
          \f]

          The D-matrices for complex spherical harmonics are complex.
          These complex D-matrices are then reverted to those for the real
          spherical harmonics, which are real.

        */
      static SxArray<SxMatrix<Double> > 
      computeYlmRotMatrices (int lmax, const SxMatrix3<Double> &S,
            const SxClebschTable &complexClebschGordan);

      /// \brief Compute rotation group in real Ylm representation
      static SxYlmRotGroup
      computeYlmRotMatrices (const SxArray<SxMatrix3<Double> > &syms, int lmax);

      /** \brief returns \f$ \frac{n!}{m!} \f$
        */
      static double facFrac (int, int);
      static double faculty (int);

      static double jsb(int l, double x);
      // Only for l <=2
      static double jsbTaylor (int l, double x);

      /** \brief Get Fk(m)/m! for 
                 \f$ \frac{d^m}{dx^m} j_L(x) = \sum F^L_k(m) j_k(x)\f$
          @param L the L value of the spherical Bessel function to differentiate
          @param maxOrder the maximum order of the derivative
          @return FkL(:k,:m)

          The computation is done recursively via
          \f[
          \frac{d}{dx}j_k(x) = \frac{k}{2k+1}j_{k-1}(x) 
                             - \frac{k+1}{2k+1}j_{k+1}(x)
          \f]
      */
      static SxArray2<double> getJsbDerivCoeffs(int L, int maxOrder);
};
#endif // _SX_YLM_H_
