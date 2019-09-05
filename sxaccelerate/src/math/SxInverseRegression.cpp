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
#include <SxInverseRegression.h>

//------------------------------------------------------------------------------
// SxInverseRegression-class
//----------------------------------------------------------------------------
// Constructors
// -- Constructor
SxInverseRegression::SxInverseRegression ()
{
   
   // empty

}

SxInverseRegression::SxInverseRegression (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the inverse least square problem can be simplified to a
   // linear one (all y-values not 0)
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   for (int i = 0; i < y.getSize (); ++i)  {
      SX_CHECK (y(i) != 0, y(i));
   }
   
   // Converting the inverse least square problem in a linear one
   // f(x) = 1. / (a + b * x)
   // => 1. / f (x) = a + b * x
   SxVector<Double> inverseYVals;

   inverseYVals = convertToLinearProblem (y); 
  
   // Solving the linear problem 
   SxLinearRegression linReg (x, inverseYVals);
   coeff = linReg.getCoeff ();
   transl = linReg.getTransl ();
  
}

// -- Destructor
SxInverseRegression::~SxInverseRegression ()
{

   // empty

}

// Methods
SxVector<Double> SxInverseRegression::convertToLinearProblem (const SxVector<Double> &y)
{

   SxVector<Double> ret (y.getSize ());

   for (int i = 0; i < ret.getSize (); ++i)  {

      ret(i) = 1. / (y(i));

   }
   
   return ret;

}

double SxInverseRegression::getVal (double x)
{

   return 1. / (transl + coeff * x);

}

SxString SxInverseRegression::getFunc ()
{

   SxString funcStr ("1. / (x * "), coeffStr (coeff), translStr (transl), helpStr (" + "), helpStr2 (")");
   return funcStr + coeffStr + helpStr + translStr + helpStr2;

}

bool SxInverseRegression::testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y)
{

   // Checking if the inverse least square problem can be simplified to a
   // linear one (all y-values not 0)
   if (x.getSize () != y.getSize ())  {

      return false;

   } else   {
   
      for (int i = 0; i < y.getSize (); ++i)  {

         if (y(i) == 0) return false;

      }
   
      return true;

   }
   
}

SxString SxInverseRegression::getErrStr ()
{

   return SxString ("The number of x-values must equal the number of y-values. All y-values must not be zero.");
   
}

