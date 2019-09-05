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
#include <SxPowerRegression.h>

//------------------------------------------------------------------------------
// SxPowerRegression-class
//----------------------------------------------------------------------------
// Constructors
// -- Constructor
SxPowerRegression::SxPowerRegression ()
{
   
   // empty

}

SxPowerRegression::SxPowerRegression (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the power least square problem can be simplified to a
   // linear one (all values positive)
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.minval () > 0  && y.minval () > 0, x.minval (), y.minval ());
  
   // Converting the power least square problem in a linear one
   // f (x) = a * x ^ b
   // <=> log (f (x)) = log (a) + b * log (x)
   SxVector<Double> logarithmizedXVals (x.getSize ()), logarithmizedYVals (y.getSize ());
   for (int i = 0; i < logarithmizedXVals.getSize (); ++i)  {
   
      logarithmizedXVals(i) = log10 (x(i));
      logarithmizedYVals(i) = log10 (y(i));
   
   }
  
   // Solving the linear problem 
   SxLinearRegression linReg (logarithmizedXVals, logarithmizedYVals);
   coeff = linReg.getCoeff ();
   transl = linReg.getTransl ();

   // Rechanging the linear problem to a power problem
   convertLinSolToPowProblem ();
  
}

// -- Destructor
SxPowerRegression::~SxPowerRegression ()
{

   // empty

}

// Methods
void SxPowerRegression::convertLinSolToPowProblem ()
{

   transl = ::pow (10., transl);

}

double SxPowerRegression::getVal (double x)
{

   return transl * ::pow (x, coeff);

}

SxString SxPowerRegression::getFunc ()
{

   SxString funcStr (" * x ^ "), coeffStr (coeff), translStr (transl);
   return translStr + funcStr + coeffStr;

}

bool SxPowerRegression::testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the power least square problem can be simplified to a
   // linear one (all values positive)
   if (x.getSize () != y.getSize () || (x.minval () <= 0  || y.minval () <= 0))  {

      return false;

   } else   {
  
      return true;

   }
   
}

SxString SxPowerRegression::getErrStr ()
{

   return SxString ("The number of x-values must equal the number of y-values. All x and y-values must be positive.");
   
}

