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

#ifndef _SX_GAUSS_QUAD_H_
#define _SX_GAUSS_QUAD_H_

#include <SxMath.h>
#include <SxVector.h>

/** \brief Gauss quadrature grids

  A Gauss quadrature formula replaces an integral by a weighted sum
  over special points. Both the integration points and weights are
  chosen such that polynomials up to a certain order are integrated
  exactly.
  Depending on an optional weight function, different types of integrations
  are available
  - Gauss-Legendre: weight 1, interval [-1,1]
  - Gauss-Laguerre: weight exp(-x), interval [0,\f$\infty\f$]

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_MATH SxGaussQuad
{
   public:
      /// Integration points
      SxVector<Double> nodes;
      
      /// Integration weights
      SxVector<Double> weights;

      /// Get size
      ssize_t getSize () const
      {
         SX_CHECK (nodes.getSize () == weights.getSize (),
                   nodes.getSize (), weights.getSize ());
         return nodes.getSize ();
      }

   protected:
      /** \brief Define recursion function pointer
        The recursion for monic orthogonal polynomials can be written as
        \f[
        p_{n+1}(x) + (b_n - x) p_n(x) + a_n p_{n-1}(x) = 0
        \f]
        This function sets a and b for a specific n >=0.
        Note that a0 should be 0.

        \note Monic orthogonal polynomials are the orthogonal
        polynomials normalized such that the highest power has
        a prefactor of one.

        Starting from a standard recursion formula
        \f[
        D_n P_{n+1}(x) = (C_n x - B_n) P_n(x) - A_n P_{n-1}(x)
        \f]
        the monic recursion coefficients become
        \f[
        a_n = \frac{A_n D_{n-1}}{C_n C_{n-1}}
        ; \quad 
        b_n = \frac{B_n}{C_n}
        \f]
        */
      typedef void (*Recursion)(int, double *a, double *b);

      /** \brief Golub-Welsch algorithm
        @param n              order of resulting quadrature
        @param recursion      recursion function for polynomial
        @param weightIntegral integral over the weight function

        The Golub-Welsch algorithm generates Gauss quadrature formulas
        from the recursion formula of the underlying polynomial.

        Basic idea of Golub-Welsch:
        We need the nodes of p_{n+1}(x) = 0.
        Regard this as a vector function \f$\mathbf p(x)\f$ and
        reinterpret recursion equations
        \f[
        p_{n+1} + b_n p_n + a_n p_{n-1} = x p_n 
        \f]
        as eigenvalue equation for tridiagonal matrix J.

        The eigenvalues can be stably computed from the symmetrized
        form (multiply J from left with diagonal matrix \f$1/\sqrt{a_n}\f$ and
        from right with diagonal matrix \f$\sqrt{a_{n-1}}\f$ )
        \f[
        \tilde J_{n,n-1}=J_{n-1,n}=\sqrt{a_n}
        \f]

        The integration weights are given by the first elements of the
        eigenvectors (don't know why...)
        \f[
        w_n = \left(\int dx w(x)\right) |\phi^{n}_1|^2 
        \f]

        */
      void GolubWelsch (int n, Recursion recursion, double weightIntegral);

      /// Recursion formula for monic Legendre polynomials
      static void recursionLegendre (int n, double *a, double *b);

      /// Recursion formula for monic Legendre polynomials
      static void recursionLaguerre (int n, double *a, double *b);

   public:
      /// Create a Gauss-Legendre quadrature
      void setupLegendre (int n)
      {
         GolubWelsch (n, &recursionLegendre, 2.); 
      }

      /// Create a Gauss-Laguerre quadrature
      void setupLaguerre (int n)
      {
         GolubWelsch (n, &recursionLaguerre, 1.); 
      }
};

#endif /* _SX_GAUSS_QUAD_H_ */
