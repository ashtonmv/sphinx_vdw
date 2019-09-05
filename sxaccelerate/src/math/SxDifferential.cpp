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

#include <SxDifferential.h>
#include <SxMatrix.h>

SxDifferential::SxDifferential (int n)
{
   SX_CHECK (n > 0, n);
   setup (n);
}

void SxDifferential::setup (int n)
{
   // note: assymmetric formulas might be numerically worse than
   //       lower-N symmetric formula (except for very small n)
   //       possible solution:
   //       use assymmetric 3-order formula for i < 3     (or N-i < 3)
   //       use symmetric i-order formula for 3 <= i <= n
   //       use symmetric N-order formula for i >= n

   SX_CHECK (n > 0, n);
   int N = 2 * n + 1;
   preFac.resize (N);
   SxVector<Double> fac(N);
   fac(0) = 1.;
   for (int i = 1; i < N; ++i) fac(i) = fac(i-1) * i;
   // --- loop over start offset
   for (int o = 0; o < N; ++o)  {
      // f(x) = \sum_j x^j / j! d^j/d z^j f(z)|_{z=0}
      // \sum_n c_n f(n) = f'(0) => \sum_n c_n n^j = j! delta(j,1)
      // thus: D(n.j) = n^j/j!, c_n = D^{-1}(n,1) 
      SxMatrix<Double> D(N,N);
      for (int i = 0; i < N; ++i)
         for (int j = 0; j < N; ++j)
            D(j,i) = pow(double(i-o), double(j)) / fac(j); // (i-o)^j/j!
      preFac(o).copy (D.inverse ().colRef (1));
   }
}

