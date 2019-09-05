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
#include <SxVector.h>
#include <SxMath.h>
#include <SxConstants.h>  /* for complex I */
#include <stdio.h>
#include <iostream>


#define TYPE_V1 Float
#define TYPE_V2 Double
#define TYPE_V3 Complex16

/** 
  \example vector.cpp
  \author  Sixten Boeck
  */
int main ()
{
   // --- initialize vector class, setting up random number generator etc.
   initSPHInXMath ();
   {

   int n = 3;
   SxVector<TYPE_V1> v1(n);
   SxVector<TYPE_V2> v2;
   SxVector<TYPE_V3> v3(3);

   // --- set up vectors: v1 in a for loop, v2 using lists
   for (int i=0; i < n; i++)
      v1(i) = sqrt((double)i);
   v2 = SxList<TYPE_V2::Type> () << 2 << 4 << 6;


   // --- if you have huge vectors ALWAYS use iterators for accessing them
   double x = 0.;
   SxVector<TYPE_V3>::Iterator it;
   for (it = v3.begin(); it != v3.end(); ++it, x+=1.5)  {
      *it = sqrt(x) + cbrt (x) * I;
   }

   // --- print vectors
   printf ("v1 = "); v1.print ();           // unformatted output
   printf ("v2 = "); v2.print (true);       // formatted output
   cout << "v3 = " << v3 << endl;           // streaming
   cout << endl;

	cout << "Select the first two elements: " << v1(SxIdx(0, 1)) << endl;
	cout << -v1 +v2 << endl;
	cout << v1 + v2 << endl;
	cout << sqr (v1(2)) << endl;
	cout << endl;

   // --- call some BLAS 1 routines and usage of temp. objects
   cout << "sum  (v1) = " << v1.sum ()     << endl; 
   cout << "prod (v2) = " << v2.product () << endl;
   cout << "   |v3|^2 = " << v3.absSqr();
   cout << "  v2 * v3 = " << (v2 * v3);  // don't forget '(' and ')'!!!
   cout << "|v2*v3|^2 = " << (v2 * v3).absSqr() << endl;

   // --- normalization
   cout << "Normalization\n";
   cout << "v1*v2            = " << (v1*v2) << endl;
   v1 = v1 * v2; v1.normalize ();
   cout << "v1*v2 / |v1*v2| = " << v1;
   cout << "           norm = " << v1.absSqr().sum() << endl;
   cout << endl;

   // --- usage of assignment operators
   cout << "Assignment operators:\n";
   cout << "v3              = " << v3 << endl;
   v3 *= 10.;
   cout << "v3 = v3 * 10    = " << v3 << endl;
   v3 *= SxComplex16(2.,3.);
   cout << "v3 = v3 * (2,3) = " << v3 << endl;
   v3 /= SxComplex16(2.,3.);
   cout << "v3 = v3 / (2,3) = " << v3 << endl;
   cout << endl;

   // --- trigonometric functions
   cout << "Vector functions\n";
   cout << "v2       = " << v2       << endl;
   cout << "sqrt(v2) = " << sqrt(v2) << endl;
	cout << "v2^(1/3.) = " << pow (v2,(1/3.)) << endl;
   cout << "exp(v2)  = " << exp(v2)   << endl;
   cout << "erf(v2)  = " << derf(v2)  << endl;
   cout << "sin(v2)  = " << sin(v2)   << endl;
   cout << "cos(v2)  = " << cos(v2)   << endl;
	cout << "atan(v2) = " << atan (v2) << endl;
   cout << "sinh(v2) = " << sinh (v2) << endl;
   cout << "cosh(v2) = " << cosh (v2) << endl;
   cout << endl;

   // --- usage of the random number generator
   v3.randomize();
   cout << "Random numbers\n";
   cout << " v3  = " << v3 << endl;
   cout << "|v3| = " << v3.absSqr().sum() << endl;
   cout << endl;

   // --- call some vector-vector operations
   cout << "Conjugate complex vector\n";
   cout << "v3^*    = " << v3.conj() << endl;
   cout << "Scalar product\n";
   cout << "<v3|v3> = " << (v3.conj() * v3).sum() << endl;
   cout << "        = " << (v3 ^ v3) << endl;  // conj is handled by '^'
   }

   return 0;
}

