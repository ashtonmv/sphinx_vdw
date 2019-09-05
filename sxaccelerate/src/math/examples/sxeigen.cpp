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

#include <SxUtil.h>
#include <SxSymMatrix.h>
#include <SxMatrix.h>
#include <SxMathLib.h>
#include <SxVector.h>
#include <SxConstants.h>  /* for complex I */
#include <stdio.h>
#include <iostream>

#define TYPE_M1 Complex8
#define TYPE_M2 Complex16
#define TYPE_M3 Float
#define TYPE_M4 Double

/** 
  \example eigen.cpp
  \author  Sixten Boeck
  */
int main ()
{
   initSPHInXMath ();

   int n = 4;
   SxMatrix<TYPE_M1>    m1(n,n);
   SxSymMatrix<TYPE_M2> m2(n);
   SxMatrix<TYPE_M3>    m3(n,n);
   SxSymMatrix<TYPE_M4> m4(n);

   // --- set up matrix m1 and m3
   int r, c;
   for (r=0; r < n; r++)  {
      for (c=r; c < n; c++)  {
         m1(r,c) = TYPE_M2::Type(sqrt((double)(r*c+1)), 
                                 cbrt((double)(r*c+1)));
         m1(c,r) = m1(r,c).conj();
         m3(r,c) = ((double)(r+c+1));
         m3(c,r) = m3(r,c);
         if (r==c)  m1(r,c).im = 0.;
      }
   }
   // --- initialize matrix m2/4. m2/4 = m1/3, but trigonally packed
   //     (use only upper-right triangle)
   for (r=0; r < n; r++)
      for (c=r; c < n; c++)  {
         m2(r,c) = TYPE_M2::Type(sqrt((double)(r*c+1)), 
                                 cbrt((double)(r*c+1)));
         m4(r,c) = sqrt((double)(r*c+1));
         if (r==c)  m2(r,c).im = 0.;
      }
               
   // --- print matrices
   cout << "m1\n" << m1 << endl;
   cout << "m2\n" << m2 << endl;
   cout << "m3\n" << m3 << endl;
   cout << "m4\n" << m4 << endl;

   // --- TODO: support '^' for trigonal matrices!!!
   cout << "m1 ^ m1 =\n" << (m1 ^ m1) << endl; 
   cout << "m1 ^ m2 =\n" << (m1 ^ m2.expand()) << endl;           // TODO
   cout << "m3 ^ m1 =\n" << (m3 ^ m1) << endl;
   cout << "m2 ^ m4 =\n" << (m2.expand() ^ m4.expand()) << endl;  // TODO

   // --- compute eigensystems
   SxMatrix<TYPE_M1>::Eigensystem eig1;
   SxMatrix<TYPE_M3>::Eigensystem eig3;
   SxSymMatrix<TYPE_M2>::Eigensystem eig2;
   SxSymMatrix<TYPE_M4>::Eigensystem eig4;
   eig1 = m1.eigensystem ();
   eig2 = m2.eigensystem ();
   eig3 = m3.eigensystem ();
   eig4 = m4.eigensystem ();

   cout << "e1: "; eig1.print(true);
   cout << "e2: "; eig2.print(true);
   cout << "e3: "; eig3.print(true);
   cout << "e4: "; eig4.print(true);

   // --- check eigensystem
   printf ("Check eigensystem\n");
   for (r=0; r < n; r++)  {
      cout << "||m1^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m1^eig1.vecs.colRef(r)) 
            -  eig1.vals(r)*eig1.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m2^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m2.expand()^eig2.vecs.colRef(r)) 
            -  eig2.vals(r)*eig2.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m3^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m3^eig3.vecs.colRef(r)) 
            -  eig3.vals(r)*eig3.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m4^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m4.expand()^eig4.vecs.colRef(r)) 
            -  eig4.vals(r)*eig4.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   cout << endl;

   // --- check norm of eigenvectors
   printf ("Check norm\n");
   for (r=0; r < n; r++)  {
      cout << "||V("<<r<<")|| = "
        << ((eig1.vecs.colRef(r)) ^ eig1.vecs.colRef(r)).chop() << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||V("<<r<<")|| = "
        << ((eig3.vecs.colRef(r)) ^ eig3.vecs.colRef(r)).chop() << endl;
   }


   // --- Cholesky decomposition, M = L * L^  (M symmetric)
   m1 = SxMatrix<TYPE_M1> (3, 3,
           SxList<TYPE_M1::Type> () << 10 <<  6       << 5
                                    <<  6 << 22       << 1.+2.*I
                                    <<  5 <<  1.-2.*I << 3
        );
   cout << "m1 = \n" << m1;
   SxMatrix<TYPE_M1> L = m1.choleskyDecomposition();
   cout << "L = \n" << L;
   cout << "L ^ L^* - M = " 
        << ((L ^ L.adjoint()) - m1).absSqr().sum() << endl;


   return 0;
}

