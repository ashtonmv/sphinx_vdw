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

#include <SxGaussQuad.h>
#include <SxSymMatrix.h>

void SxGaussQuad::GolubWelsch (int n, Recursion recursion, 
                               double weightIntegral)
{
   SxSymMatrix<Double> symJ (n);
   symJ.set (0.);
   double an, bn;
   for (int i = 0; i < n; ++i)  {
      (*recursion)(i, &an, &bn);
      symJ (i, i) = bn;
      if (i > 0) symJ(i-1,i) = sqrt(an);
   }
   SxSymMatrix<Double>::Eigensystem eig = symJ.eigensystem ();
   nodes = eig.vals;
   weights = weightIntegral * eig.vecs.row (0).absSqr ();
}

void SxGaussQuad::recursionLegendre (int n, double *a, double *b)
{
   // The standard recursion is
   // (n+1) P_{n+1} = (2n + 1)x P_n - n P_{n-1}
   SX_CHECK (a);
   SX_CHECK (b);
   *a = double(n * n) / double((2*n + 1)*(2 * n - 1));
   *b = 0.;
}

void SxGaussQuad::recursionLaguerre (int n, double *a, double *b)
{
   // The standard recursion is
   // (n+1) L_{n+1} = (2n+1 - x)L_n - n L_{n-1} 
   SX_CHECK (a);
   SX_CHECK (b);
   *a = n * n;
   *b = 2 * n + 1;
}

