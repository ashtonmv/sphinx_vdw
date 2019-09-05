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
#include <SxExponentialRegression.h>

//------------------------------------------------------------------------------
// SxExponentialRegression-class
//----------------------------------------------------------------------------
// Constructors
// -- Constructor
SxExponentialRegression::SxExponentialRegression ()
{
   
   // empty

}

SxExponentialRegression::SxExponentialRegression (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the exponential least square problem can be simplified to a
   // linear one (smallest y-value still bigger than 0)
   SX_CHECK (x.getSize () == y.getSize () && y.minval () > 0, x.getSize (), y.getSize (), y.minval ());
   
   // Converting the exponential least square problem in a linear one
   // f(x) = a * exp (b * x)
   // => log (f (x)) = log (a) + b * x
   // <=> p[0] + p[1] * x = y
   SxVector<Double> logarithmizedYVals;
   logarithmizedYVals = convertToLinearProblem (y); 
  
   // Solving the linear problem 
   computeLeastSquares (x, logarithmizedYVals);

   // Rechanging the problem to an exponential one again (changing the calculated coefficients)
   convertLinSolToExpProblem ();
   
}

// -- Destructor
SxExponentialRegression::~SxExponentialRegression ()
{

   // empty

}

SxVector<Double> SxExponentialRegression::convertToLinearProblem (const SxVector<Double> &y)
{

   SxVector<Double> ret (y.getSize ());

   const double factor = (1. / log(exp(1.)));

   for (int i = 0; i < ret.getSize (); ++i)  {

      ret(i) = factor * log (y(i));

   }
   
   return ret;

}

void SxExponentialRegression::convertLinSolToExpProblem ()
{

   p1 = exp (p1);

}

double SxExponentialRegression::getVal (double x)
{

   return exp (p2 * x) * p1;

}

SxString SxExponentialRegression::getFunc ()
{

   SxString funcStr ("g (x) = exp ("), p1Str (p1), p2Str (p2), helpStr (" * x) * ");
   return funcStr + p2Str + helpStr + p1Str;

}

// Methods
void SxExponentialRegression::computeLeastSquares (const SxVector<Double> &x,const SxVector<Double> &y)
{

   ssize_t n (x.getSize ());
   double a12, a22, b1, b2;

   a12 = x.sum ();
   a22 = (x.sqr ()).sum ();

   b1 = y.sum ();
   //b2 = (x*y).sum ();
   b2 = dot(x,y);

   // Translation
   p2 = (b2 - 1./(double)n * a12 * b1) / (a22 - 1./(double)n * ::pow (a12, 2));

   // Regression coefficient
   p1 = 1. / (double)n * b1 - 1. / (double)n * a12 * p2; 

}

bool SxExponentialRegression::testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the preconditions for solving the exponential least square
   // problem are fullfilled (smallest y-value still bigger than 0)
   if (x.getSize () != y.getSize () || y.minval () <= 0)  {
   
      return false;
   
   } else  {
   
      return true;
   
   }
   
}

SxString SxExponentialRegression::getErrStr ()
{

   return SxString ("The number of x-values must equal the number of y-values. All y-values must be positive.");
   
}

