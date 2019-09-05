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
#include <SxMatrix.h>
#include <SxSymMatrix.h>
#include <stdio.h>
#include <iostream>

#define TYPE_M1 Float
#define TYPE_M2 Double
#define TYPE_M3 Complex8

/** \example matrix.cpp 

  Using matrices in the SFHIngX Algebra Library

  \author Sixten Boeck
  */
int main ()
{
   // --- initialize vector/matrix class, setting up random number generator etc.
   initSPHInXMath ();

   int n = 5;
   SxMatrix<TYPE_M1> m1(n,n);
   SxMatrix<TYPE_M2> m2, i(3,3);
   SxMatrix<TYPE_M3> m3(3,3);

   // --- set up matrices: m1 in a for loop, m2 using lists
   int r, c;
   for (r=0; r < n; r++)
      for (c=0; c < n; c++)
         m1(r,c) = sqrt((double)(r*c));
   m2 = SxMatrix<TYPE_M2> ( 3, 3,
           SxList<TYPE_M2::Type> () << 1 << 2 << 3
                                    << 4 << 5 << 6
                                    << 7 << 8 << 9
        );

   // --- using identity matrix
   i=i.identity ();

   // --- if you have huge vectors ALWAYS use iterators for accesing them
   double x = 0.;
   SxMatrix<TYPE_M3>::Iterator it;
   for (it = m3.begin(); it != m3.end(); ++it, x+=1.5)  {
      *it = SxComplex8(sqrt(x), cbrt (x));
   }

   // --- print matrices
   printf ("m1\n"); m1.print (); cout << endl; // unformatted output
   printf ("m2\n"); m2.print (true);           // formatted output
   cout << "m3\n" << m3 << endl;               // streaming
   cout << endl;
   cout << "i \n" << i  << endl;
   cout << endl;

   // --- initialization of non-quadratic matrices
   m2 = SxMatrix<TYPE_M2> ( 2, 5,
           SxList<TYPE_M2::Type> () << 1 << 2 << 3 << 4 << 5 
                                    << 6 << 7 << 8 << 9 << 10
        );
   cout << " 2 x 5 matrix =\n" << m2 << endl;
   cout << endl;

   // --- in case of BLAS 1 matrices are nothing but vectors
   cout << "sum  (m1) = "  << m1.sum()    << endl; 
   cout << "prod (m2) = "  << m2.sum()    << endl;
   cout << "   |m3|^2 =\n" << m3.absSqr() << endl;
   cout << "  m3 * m3 =\n" << (m3 * m3)   << endl;  // don't forget '(' and ')'
   cout << "|m3*m3|^2 =\n" << (m3 * m3).absSqr() << endl;
   cout << endl;

   // --- usage of assignment operators
   cout << "Assignment operators:\n";
   cout << "m3\n"                << m3 << endl;
   m3 *= 10.;
   cout << "m3 = m3 * 10    =\n" << m3 << endl;
   m3 *= SxComplex16(2.,3.);
   cout << "m3 = m3 * (2,3) =\n" << m3 << endl;
   m3 /= SxComplex16(2.,3.);
   cout << "m3 = m3 / (2,3) =\n" << m3 << endl;
   cout << endl;

   // --- trigonometric functions
   cout << "Matrix functions\n";
   cout << "m2       =\n" << m2       << endl;
   cout << "sqrt(m2) =\n" << sqrt(m2) << endl;
   cout << "exp(m2)  =\n" << exp(m2)  << endl;
   cout << "erf(m2)  =\n" << derf(m2) << endl;
   cout << "sin(m2)  =\n" << sin(m2)  << endl;
   cout << "cos(m2)  =\n" << cos(m2)  << endl;
   cout << "sinh(m2) =\n" << sinh(m2) << endl;
   cout << "cosh(m2) =\n" << cosh(m2) << endl;
   cout << endl;

   // --- usage of the random number generator
   m3.randomize();
   cout << "Random numbers (normed columnwise)\n";
   cout << " m3  =\n" << m3 << endl;
   cout << "|m3| = "  << m3.absSqr().sum() << endl;
   cout << endl;

   // --- ectract rows and columns
   cout << "Rows, columns, diagonal\n";
   cout << "m1\n" << m1 << endl;
   cout << "row m1(1:c) = " << m1.row(1) << endl;
   cout << "col m1(r:2) = " << m1.colRef(2) << endl;
   cout << "dg  m1      = " << m1.diag() << endl;
   m1.row(2)(0) = 5.;
   cout << "m1 (element m1(0,2) is untouched)\n" << m1 << endl;
   m1.colRef(2)(0) = 5.;
   cout << "m1 (element m1(0,2)=5 now!!!)\n" << m1 << endl;
   cout << endl;

   // --- treat complex matrices
   cout << "Matrix operations\n";
   cout << "m3    =\n" << m3              << endl;
   cout << "m3^*  =\n" << m3.conj()       << endl;   
   cout << "m3^t  =\n" << m3.transpose()  << endl;
   cout << "m3^t* =\n" << m3.adjoint()    << endl;
   cout << endl;

   // --- BLAS level 3
   cout << "Matrix-Matrix operations\n";
   cout << "m1\n" << m1 << endl;
   cout << "m1^m1   =\n" << (m1^m1) << endl;  // don't forget '(' and ')'
   cout << "m3\n" << m3 << endl;
   cout << "m3^m3   =\n" << (m3^m3) << endl;  // don't forget '(' and ')'
   cout << "m3* ^m3 =\n" << (m3.adjoint()^m3) << endl; 

   m1 = SxMatrix<TYPE_M1> ( 3, 3, 
           SxList<TYPE_M1::Type> () <<  1 << 2 << 3 
                                    << -2 << 5 << 6
                                    << -3 << 8 << 9
        );
   cout << "m1      =\n" << m1         << endl;
   cout << "m1^-1   =\n" << m1.inverse() << endl;

   cout << "m3      =\n" << m3         << endl;
   cout << "m3^-1   =\n" << m3.inverse() << endl;
   cout << endl;

   // --- symmetric matrices
   cout << "Symmetric matrices:\n";
   SxSymMatrix<Double> symMat (    // input matrix
         SxList<double> () << 2. << 1. << 0. << 0.
                                 << 2. << 1. << 0.
                                       << 2. << 1.
                                             << 1.
   );
   SxSymMatrix<Double> symMatInv ( // expected result of symMat^-1
         SxList<double> () << 1. << -1. <<  1. << -1.
                                 <<  2. << -2. <<  2.
                                        <<  3. << -3.
                                               <<  4.
   );
   cout << "S:\n"; symMat.print(true);
   cout << "S^-1:\n";
   symMat.inverse().print(true);
   cout << "should be:\n";
   symMatInv.print(true);
   cout << "Error (should be 0): " 
        << (symMat.inverse() - symMatInv).absSqr().sum() << endl;

   return 0;
}

