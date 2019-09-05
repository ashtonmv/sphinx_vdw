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
#include <SxLogarithmicRegression.h>

//------------------------------------------------------------------------------
// SxLogarithmicRegression-class
//----------------------------------------------------------------------------
// Constructors
// -- Constructor
SxLogarithmicRegression::SxLogarithmicRegression ()
{
   
   // empty

}

SxLogarithmicRegression::SxLogarithmicRegression (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the logarithmic least square problem can be simplified to a
   // linear one (smallest x-value still bigger than 0)
   SX_CHECK (x.getSize () == y.getSize () && x.minval () > 0, x.getSize (), y.getSize (), x.minval ());
   
   // Converting the logarithmic least square problem in a linear one
   // f(x) = a + b * log (x)
   // => f (x) = a + b * y, y := log (x)
   SxVector<Double> logarithmizedXVals;

   logarithmizedXVals = convertToLinearProblem (x); 
  
   // Solving the linear problem 
   SxLinearRegression linReg (logarithmizedXVals, y);
   coeff = linReg.getCoeff ();
   transl = linReg.getTransl ();
  
}

// -- Destructor
SxLogarithmicRegression::~SxLogarithmicRegression ()
{

   // empty

}

// Methods
SxVector<Double> SxLogarithmicRegression::convertToLinearProblem (const SxVector<Double> &x)
{

   SxVector<Double> ret (x.getSize ());

   for (int i = 0; i < ret.getSize (); ++i)  {

      ret(i) = log (x(i));

   }
   
   return ret;

}

double SxLogarithmicRegression::getVal (double x)
{

   return transl + coeff * log (x);

}

SxString SxLogarithmicRegression::getFunc ()
{

   SxString funcStr ("log (x) * "), coeffStr (coeff), translStr (transl), helpStr (" + ");
   return funcStr + coeffStr + helpStr + translStr;

}

bool SxLogarithmicRegression::testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the logarithmic least square problem can be simplified to a
   // linear one (smallest x-value still bigger than 0)
   if (x.getSize () != y.getSize () || x.minval () <= 0)  {

      return false;

   } else   {
  
      return true;

   }
   
}

SxString SxLogarithmicRegression::getErrStr ()
{

   return SxString ("The number of x-values must equal the number of y-values. All x-values must be positive.");
   
}

