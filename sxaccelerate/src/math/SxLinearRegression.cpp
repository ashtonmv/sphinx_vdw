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

#include <SxError.h>
#include <SxLinearRegression.h>

//------------------------------------------------------------------------------
// SxLinearRegression-class
//------------------------------------------------------------------------------
// Constructors
// -- Constructor
SxLinearRegression::SxLinearRegression ()
{
   // empty
}

SxLinearRegression::SxLinearRegression (const SxVector<Double> &x, const SxVector<Double> &y)
   : xVals (x)
{
   computeLeastSquares (x,y);
}

// -- Destructor
SxLinearRegression::~SxLinearRegression ()
{
   // empty
}

double SxLinearRegression::getVal (double x)
{
   return p2 * x + p1;
}

SxString SxLinearRegression::getPolynomial ()
{
   SxString funcStr ("g (x) = "), p1Str (p1), p2Str (p2), helpStr (" * x +");
   return funcStr + p2Str + helpStr + p1Str;
}

double SxLinearRegression::getCoeff ()
{
   return p2;
}

double SxLinearRegression::getTransl ()
{
   return p1;
}

// Methods
void SxLinearRegression::computeLeastSquares (const SxVector<Double> &x,const SxVector<Double> &y)
{
   ssize_t n (x.getSize ());
   double a12, a22, b1, b2;
   a12 = x.sum ();
   a22 = (x.sqr ()).sum ();

   b1 = y.sum ();
   //b2 = (x*y).sum ();
   b2 = dot (x,y);

   // Translation
   p2 = (b2 - 1./(double)n * a12 * b1) / (a22 - 1./(double)n * ::pow (a12, 2));
   // Regression coefficient
   p1 = 1. / (double)n * b1 - 1. / (double)n * a12 * p2; 
}

