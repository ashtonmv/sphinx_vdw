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
#include <SxNaturalCubicSpline.h>

//------------------------------------------------------------------------------
// SxNaturalCubicSpline-class
//------------------------------------------------------------------------------
// Constructors
// -- Constructor
SxNaturalCubicSpline::SxNaturalCubicSpline ()
{
   // empty
}

SxNaturalCubicSpline::SxNaturalCubicSpline (const SxVector<Double> &x, const SxVector<Double> &y)
   : xVals (x)
{
   compute (x,y);
}

// -- Destructor
SxNaturalCubicSpline::~SxNaturalCubicSpline ()
{
   // empty
}

void SxNaturalCubicSpline::print (const SxVector<Double> &diag1, const SxVector<Double> &diag2, const SxVector<Double> &diag3, const SxVector<Double> &b)
{
  ssize_t n = diag2.getSize ();
  
   // TESTOUTPUT
  for (int i = 0; i < n; ++i) {
     for (int iCol = 0; iCol < i - 1; ++iCol)  {
        cout << setw (9) << 0.;
     }
     // lower
     if (i-1 >= 0)  {
        cout << setw (9)<< diag3(i-1);
     }
     // mid
     cout << setw (9) << diag2(i);
     // upper
     if (i < n - 1) {
        cout << setw (9) << diag1(i);
     }
     for (int iCol = i + 2; iCol < n; ++iCol)  {
        cout << setw (9) << 0.;
     }
     cout << setw (9) << b(i) << endl;		
   }
}

//extern "C" {
//#include <f2c.h>
//#include <clapack.h>
//}
SxVector<Double> SxNaturalCubicSpline::symmetricTridiagonalGauss (
      const SxVector<Double> & upper,
      const SxVector<Double> &mid, 
      const SxVector<Double> &lower,
      const SxVector<Double> &rhs)
{
   /* --- LAPACK version
   integer N = mid.getSize (), nrhs=1, info = 0;
   //dgtsv_(&N, &nrhs, lower.elements, mid.elements, upper.elements, 
   //       rhs.elements, &N, &info);

   char fact='N', trans='N';
   SxVector<Double> dlf(N-1), df(N), duf(N-1), du2(N-2), res(N), work(3*N);
   integer *ipiv = new integer[N],
           *iwork = new integer[N];
   double rcond, ferr, berr;
   dgtsvx_ (&fact, &trans, &N, &nrhs,
            lower.elements, mid.elements, upper.elements,
            dlf.elements, df.elements, duf.elements, du2.elements, ipiv,
            rhs.elements, &N, res.elements, &N,
            &rcond, &ferr, &berr, work.elements, iwork, &info);
   delete iwork;
   delete ipiv; 

   if (info != 0) {
      cout << "dgtsv failed: info=" << info << endl;
      SX_EXIT;
   }
   return res;
   // return rhs;
   */
   SxVector<Double> diag1 (upper), diag2 (mid), diag3 (lower), b (rhs);
   ssize_t n = diag2.getSize ();
  
  // Overview:
  // ---------  
  // diag2(0) diag1(0) 0        0
  // diag3(0) diag2(1) diag1(0) 0
  // 0        ...
  // ...
  // 0        0        0        diag3(n-2) diag2(n-1) diag1(n-1)
  // 0        0        0        0          diag3(n-1) diag2 (n)
  
//  cout << "Before Gaussian Elimination: " << endl; 
//  print (diag1, diag2, diag3, b);
  	
  // Gaussian Elimination
  for (ssize_t i = 0; i < n - 1; ++i)  {
     diag2(i + 1) -= diag1 (i) * diag3 (i) / diag2(i);
     b(i + 1) -= b(i) * diag3(i) / diag2(i);
     diag3(i) = 0;
  }

//  cout << "Gaussian Elimination: " << endl;
//  print (diag1, diag2, diag3, b);
  

  // Backward substitution
  b(n-1) /= diag2(n-1);
  diag2(n-1) = 1;
  for (ssize_t i = n - 2; i >= 0; --i)  {
     b(i) -= diag1(i) * b(i + 1);
     diag1(i) = 0;
     b(i) /= diag2(i);
     diag2(i) = 1;
  } 
//  cout << "Backward substitution: " << endl;
//  print (diag1, diag2, diag3, b);

  return b;
}

double SxNaturalCubicSpline::getValXY (double x) const
{
   // Introducing epsilon which is supposed to guarantee that even if
   // floating point errors occured one would receive the corresponding
   // y-value for the last x-value
   SX_CHECK ((xVals(0) <= x && xVals(xVals.getSize () - 1) + 1e-8>= x), xVals(0), x, xVals(xVals.getSize () - 1));
   
   for (int i = 0; i< xVals.getSize (); ++i){
      if (i == xVals.getSize () - 1 || (xVals(i)<= x && x < xVals(i+1)))  {
         SxVector<Double> power (4);
         double dx = x - xVals(i);
         power(0) = 1;
         power(1) = dx;
         power(2) = dx * dx;
         power(3) = dx * dx * dx;
         //cout << "Hello!!!" << endl;         
         return dot (polyCoeff(i), power);
      }     
   }

   SX_EXIT;
   
}


// Methods
void SxNaturalCubicSpline::compute (const SxVector<Double> &x,const SxVector<Double> &y)
{
   // Testing if the sizes are the same
   SX_CHECK (x.getSize () == y.getSize ());
   
   ssize_t n (x.getSize () - 2);
  
   SxVector<Double> h (n + 1), gamma (n);
   for (int i = 0; i < n + 1; ++i)  {
      h(i) = x(i + 1) - x(i);
      if (i > 0 && i < n + 1) {
         gamma(i-1) = 6. * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i-1));
      }
   }
   
   SxVector<Double> diag1 (n-1), diag2 (n), diag3 (n-1), b (n);
   
   // Init vectors
  for (int i = 0; i < n - 1; ++i)  {
     diag1(i) = h(i + 1);
     diag2(i) = 2. * (h(i) + h (i + 1));
     diag3(i) = h(i + 1);
     b(i) = gamma(i);
  }
  b(n - 1) = gamma(n-1);
  diag2(n - 1) = 2. * (h(n - 1) + h(n));

  SxVector<Double> beta (n+2);
  beta(0) = 0.;
  beta(n + 1) = 0.;
  b = symmetricTridiagonalGauss (diag1, diag2, diag3, b);
  for (int i = 0; i < b.getSize (); ++i)  {
     beta(i+1) = b(i);
  }

/*  cout << "Beta: " << endl;
  for (int i = 0; i < beta.getSize (); ++i){
     cout << beta(i) << " | ";
  }
  cout << endl;
*/
  SxVector<Double> alpha (beta.getSize () - 1);

//  cout << "Alpha:" << endl;
  for (int i = 0; i < alpha.getSize (); ++i)  {
     alpha (i) = (y(i + 1) - y(i)) / h(i) - 1. / 3. * beta(i) * h(i) - 1. / 6. * beta(i+1) * h(i);
//     cout << alpha(i) << " | ";
  }
//  cout << endl;

  polyCoeff.resize (y.getSize ()); 
  //polyCoeff.resize (alpha.getSize ());
  //cout << "polyCoeff"<< polyCoeff.getSize () << endl;
  //cout << "y" << y.getSize () << endl;
  for (int i = 0; i < polyCoeff.getSize (); ++i)  {
     if (i < alpha.getSize ())  {
        polyCoeff(i).resize (4);
        polyCoeff(i)(0) = y(i);
        polyCoeff(i)(1) = alpha(i);
        polyCoeff(i)(2) = beta(i)/2.;
        polyCoeff(i)(3) = (beta(i + 1) - beta(i)) / (6. * h(i));   
     } else  {
        // Necessary to maintain all y-values
        polyCoeff(i).resize (1);
        polyCoeff(i)(0) = y(i);
     }  
  }
}
