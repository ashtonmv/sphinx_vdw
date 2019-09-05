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


//------------------------------------------------------------------------------
// SxCubicSpline-class
//------------------------------------------------------------------------------
// Constructors
// -- Constructor
template <class V>
SxCubicSpline<V>::SxCubicSpline ()
{
   // empty
}

template <class V>
SxCubicSpline<V>::SxCubicSpline (const V &x, const V &spline)
{
   SX_CHECK(x.getSize () > 0); 
   SX_CHECK(spline.getSize () > 0);
   SX_CHECK(spline.getSize () % 4 == 0);
   VALIDATE_VECTOR(x);
   VALIDATE_VECTOR(spline);
   xVals = x.getCopy();
   setSpline (spline);
}

// calculate spline constructor
template <class V>
SxCubicSpline<V>::SxCubicSpline (const V &x, const V &y, const enum SplineType splineTypeIn,
                                 const s slope, const s slope2)
{
   SX_CHECK(x.getSize () == y.getSize ());
   VALIDATE_VECTOR(x);
   VALIDATE_VECTOR(y);
   xVals = x.getCopy();
   yVals = y.getCopy();
   splineType = splineTypeIn;
   fitType = None;
   setH();
   leftHermite = slope;
   rightHermite = slope2;
   if (splineType == NaturalHermite) rightHermite = slope;
   computeSpline ();
}

template <class V>
SxCubicSpline<V>::SxCubicSpline (const V &xData, const V &yData, const V &basis, 
                                 const enum SplineType splineTypeIn,
                                 const enum FitType fitTypeIn,
                                 const s slope, const s slope2)
{
   SxVector<Int> pointsPerSpline = getPointsPerSpline (basis, xData);
   xVals = cleanZeroDepths(basis,pointsPerSpline);
   yVals.resize(xVals.getSize());
   yVals.set(0.0);
   xFit = xData.getCopy();
   yFit = yData.getCopy();
   SX_CHECK(xFit.getSize () == yFit.getSize ());
   SX_CHECK(xVals.getSize() > 0);
   VALIDATE_VECTOR(xFit);
   VALIDATE_VECTOR(yFit);
   VALIDATE_VECTOR(xVals);
   splineType = splineTypeIn;
   fitType = fitTypeIn;
   setH();
   leftHermite = slope;
   rightHermite = slope2;
   if (splineType == NaturalHermite) rightHermite = slope;
   computeLeastSquareFit ();
   if (xVals.getSize() != basis.getSize()) {
      SxCubicSpline<V> spline(xVals, yVals, splineType, leftHermite, rightHermite);
      xVals = basis.getCopy ();
      yVals = spline.getY(xVals);
   }
   computeSpline ();
}

template <class V>
SxCubicSpline<V>::SxCubicSpline (const V &x, const V &dx, const V &y, 
                                 enum SplineType splineTypeIn,
                                 const s slope, const s slope2)
{
   xVals = x.getCopy();
   yVals.resize(xVals.getSize());
   yVals.set(0.0);
   hVals = dx.getCopy();
   yVals = y.getCopy();
   SX_CHECK(xVals.getSize () == yVals.getSize ());
   SX_CHECK(xVals.getSize () == hVals.getSize ());
   VALIDATE_VECTOR(xVals);
   VALIDATE_VECTOR(hVals);
   VALIDATE_VECTOR(yVals);
   splineType = splineTypeIn;
   fitType = None;
   leftHermite = slope;
   rightHermite = slope2;
   if (splineType == NaturalHermite) rightHermite = slope;
   computeSpline();
}

// -- Destructor
template <class V>
SxCubicSpline<V>::~SxCubicSpline ()
{
   // empty
}

template <class V>
V SxCubicSpline<V>::getYFit () const
{
   SX_CHECK(yVals.getSize () > 0);
   VALIDATE_VECTOR(yVals);
   return yVals.getCopy ();
}

template <class V>
inline typename SxCubicSpline<V>::s SxCubicSpline<V>::getY (
      const typename SxCubicSpline::s x) const
{
  SX_CHECK(polyCoeff.getSize () > 0);
  SX_CHECK(polyCoeff.getSize () % 4 == 0);
  int iSpline = getSplineIdx(x, xVals);
  s diff = x - xVals(iSpline);
  SX_CHECK(polyCoeff.getSize () > 4*iSpline + 3, iSpline, polyCoeff.getSize ());
  s result = polyCoeff(4*iSpline + 0) + diff
              * ( polyCoeff(4*iSpline + 1) + diff 
              * ( polyCoeff(4*iSpline + 2) + diff 
              *   polyCoeff(4*iSpline + 3)));     
  SX_CHECK_NUM(result);    
  return result;
}

template <class V>
V SxCubicSpline<V>::getY (const V &x) const
{
   int dim = (int)x.getSize ();
   SX_CHECK(dim > 0);
   VALIDATE_VECTOR(x);
   V result (dim);
   for(int i = 0; i < dim; i++)   {
      result(i) = getY(x(i));
   }
   return result;
}

template <class V>
inline typename SxCubicSpline<V>::s SxCubicSpline<V>::getdYdX (
      const typename SxCubicSpline::s x) const
{
  SX_CHECK(polyCoeff.getSize () > 0);
  SX_CHECK(polyCoeff.getSize () % 4 == 0);
  int iSpline = getSplineIdx(x, xVals);
  s diff = x - xVals(iSpline);
  SX_CHECK(polyCoeff.getSize () > 4*iSpline + 3, iSpline, polyCoeff.getSize ());
  s result = polyCoeff(4*iSpline + 1) + diff
              * ( 2.0 * polyCoeff(4*iSpline + 2) + diff 
              *   3.0 * polyCoeff(4*iSpline + 3));     
  SX_CHECK_NUM(result);    
  return result;
}

template <class V>
V SxCubicSpline<V>::getdYdX (const V &x) const
{
   int dim = (int)x.getSize ();
   SX_CHECK(dim > 0);
   VALIDATE_VECTOR(x);
   V result (dim);
   for(int i = 0; i < dim; i++)   {
      result(i) = getdYdX(x(i));
   }
   return result;
}

// cubic spline interpolation
template <class V>
void SxCubicSpline<V>::computeSpline ()
{
   SX_CLOCK(Timer::computeSpline);
   VALIDATE_VECTOR (xVals);
   VALIDATE_VECTOR (yVals);
   int dim = (int)yVals.getSize ();
   SX_CHECK (dim == xVals.getSize (), dim, xVals.getSize ());
   setT(splineType);
   setBeta(splineType);
 
   //V polyCoeff2 = (TMat.inverse () ^ betaVals);

   // Use trigonal Matrix decomposition
   V polyCoeff2 = symmetricTridiagonalGauss (TMat, betaVals);

   polyCoeff.resize (4*dim);

   // Polynom coefficients 0 and 2 can be directly set
   for (int i = 0; i < dim; i++)   {
      polyCoeff(4*i+0) = yVals(i);
      polyCoeff(4*i+2) = polyCoeff2(i);
   }

   // Polynom coefficients 1 and 3 need higher 0 and 2
   for (int i = 0; i < dim - 1; i++)   {
      polyCoeff(4*i+3) 
         = (polyCoeff(4*(i+1)+2) - polyCoeff(4*i+2)) / (3.0 * hVals(i));
      polyCoeff(4*i+1) 
         = (polyCoeff(4*(i+1)+0) - polyCoeff(4*i+0)) / hVals(i) 
         - hVals(i) * (2.*polyCoeff(4*i+2) + polyCoeff(4*(i+1)+2)) / 3.0;
   }
   // outermost coefficients consistently set via interpolation
   // S(x) = a + b(x-x_n) + c(x-x_n)^2 + d(x-x_n^3)
   // d_n = d_{n+1}
   polyCoeff(4*(dim-1)+3) = polyCoeff(4*(dim-2)+3);
   //b_{n+1} = b_n + 2c_n*(x_{n+1}-x_n) + 3d_n(x_{n+1}-x_n)^2
   polyCoeff(4*(dim-1)+1) 
      = polyCoeff(4*(dim-2)+1) 
      + 2.*polyCoeff(4*(dim-2)+2)*hVals(dim-2)
      + 3.*polyCoeff(4*(dim-2)+2)*hVals(dim-2)*hVals(dim-2);
}

template <class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::calcDRDataDRGrid ()
{

   size_t nSpline = xVals.getSize ();
   size_t nData = xFit.getSize();
   SX_CHECK(nSpline > 0);
   SX_CHECK(nData > 0);

   setH ();
   setT(splineType);

   // calculate coefficients derivatives
   // y = S_n(x) = a + b*(x-x_n) + c*(x-x_n)^2 + d*(x-x_n)^3
   SxArray<M> dPolyCoeff (4);
   // da/dy 
   dPolyCoeff(0) = dady();
   // dc/dy
   dPolyCoeff(2) = dcdy(splineType);
   // db/dy
   dPolyCoeff(1) = dbdy(dPolyCoeff(0), dPolyCoeff(2));
   // dd/dy
   dPolyCoeff(3) = dddy(dPolyCoeff(2));
   
   M result (nData,nSpline);
   result.set(0.0);
   for(size_t iData = 0; iData < nData; iData++)   {
      int iSpline = getSplineIdx(xFit(iData), xVals);
      s diff = xFit(iData)-xVals(iSpline);
      for(size_t jSpline = 0; jSpline < nSpline; jSpline++)   {
      result(iData,jSpline)
         += dPolyCoeff(0)(iSpline,jSpline) 
         + diff * (dPolyCoeff(1)(iSpline,jSpline) 
         + diff * (dPolyCoeff(2)(iSpline,jSpline) 
         + diff * dPolyCoeff(3)(iSpline,jSpline)));
      }
   }

   return result;
}

// least square fit
template <class V>
void SxCubicSpline<V>::computeLeastSquareFit ()
{
   SX_CLOCK (Timer::leastSquareFit);

   enum SplineType internSType = splineType;
   double originalLeftHermite = leftHermite;

   // Mirror Datapoints if needed
   if (fitType == MirrorPlane || fitType == MirrorPoint) {
      mirror();
      if (splineType == Natural || splineType == HermiteNatural)
         internSType = Natural;
      else {
         internSType = Hermite;
         leftHermite = -rightHermite;
         if (fitType == MirrorPoint) leftHermite = rightHermite;
      }
   }

   ssize_t nSpline = xVals.getSize ();
   ssize_t nData = xFit.getSize();
   SX_CHECK(nData == yFit.getSize ());

   //setup T
   setT(internSType);

   // calculate coefficients derivatives
   // y = S_n(x) = a + b*(x-x_n) + c*(x-x_n)^2 + d*(x-x_n)^3
   SxArray<M> dPolyCoeff (4);
   // da/dy 
   dPolyCoeff(0) = dady();
   // dc/dy
   dPolyCoeff(2) = dcdy(internSType);
   // db/dy
   dPolyCoeff(1) = dbdy(dPolyCoeff(0), dPolyCoeff(2));
   // dd/dy
   dPolyCoeff(3) = dddy(dPolyCoeff(2));

   // Setup delta_i^q = (x_i - x_n)^q, q = 0..6
   // and setup resultVector for system of linear equation
   // fq = f(x_nu) * (x_nu-x_n)^k, k= 0..3
   SxArray<int> pointsPerSpline(nSpline);
   pointsPerSpline.set(0);
   SxArray<V> q (7);
   SxArray<V> fq (4);
   for (int i = 0; i < 7; i++)   {
      q(i).resize(nSpline);
      q(i).set(0.0);
      if (i < 4)   {
         fq(i).resize(nSpline);
         fq(i).set(0.0);
      }
   }
   for (ssize_t iVal = 0; iVal < nData; iVal++) {
      int iSpline = getSplineIdx(xFit(iVal), xVals);
      pointsPerSpline(iSpline)++;
      s diff = xFit(iVal)-xVals(iSpline);
      for (int i = 0; i < 7; i++)   {
         s powDiff = 1.0;
         for (int j = 1; j <= i; j++)   {
            powDiff *= diff;
         }
         q(i)(iSpline) += powDiff;
         if (i < 4) fq(i)(iSpline) += yFit(iVal) * powDiff;
      }
   }
   SxArray<SxList<int> > splineExtend = getSplineExtend (pointsPerSpline);

   for (ssize_t iVal = 0; iVal < nData; iVal++) {
      int iSpline = getSplineIdx(xFit(iVal), xVals);
      SxList<int>::Iterator it;
      for(it = splineExtend(iSpline).begin(); it != splineExtend(iSpline).end(); it++)  {
         s diff = xFit(iVal)-xVals(*it);
         for (int i = 0; i < 7; i++)   {
            s powDiff = 1.0;
            for (int j = 1; j <= i; j++)   {
               powDiff *= diff;
            }
            q(i)(*it) += powDiff;
            if (i < 4) fq(i)(*it) += yFit(iVal) * powDiff;
         }
      }
   }

   // Set up R_jk = da_k/df_m * delta_i^(j+k)
   SxArray<SxArray<M> > R (4);
   for (int iCoeff = 0; iCoeff < 4; iCoeff++)   {
      R(iCoeff).resize(4);
      for (int jCoeff = 0; jCoeff < 4; jCoeff++)   {
         R(iCoeff)(jCoeff).reformat(nSpline, nSpline);
         for (int iCol = 0; iCol < nSpline; iCol++)  {
            R(iCoeff)(jCoeff).colRef(iCol) <<= 
               dPolyCoeff(jCoeff).colRef(iCol) * q(iCoeff + jCoeff);
         }
      }
   }

   //Setup Mat = da_j/df_p * R_j,k
   M Mat (nSpline, nSpline);
   Mat.set(0.0);
   for(int iCoeff = 0; iCoeff < 4; iCoeff++)   {
      M dPCDagger = dPolyCoeff(iCoeff).adjoint();
      M RSum = R(iCoeff)(0);
      for(int jCoeff = 1; jCoeff < 4; jCoeff++)   {
         RSum += R(iCoeff)(jCoeff);
      }
      Mat += (dPCDagger ^ RSum);
   }
    
   // beta = da_k/df_m * fq
   V beta (nSpline);
   beta.set(0.0);
   for(int iCoeff = 0; iCoeff < 4; iCoeff++)   {
      VALIDATE_VECTOR(fq(iCoeff));
      VALIDATE_VECTOR(dPolyCoeff(iCoeff));
      beta += dPolyCoeff(iCoeff).adjoint() ^ fq(iCoeff);
   }

   // Solve linear equation
   yVals = Mat.solve(beta);
   VALIDATE_VECTOR(yVals);

   if (fitType == MirrorPlane || fitType == MirrorPoint) deMirror();
   leftHermite = originalLeftHermite;
}

// Setup differences in xVals
template <class V>
void SxCubicSpline<V>::setH ()
{
   SX_CHECK(xVals.getSize() > 0);
   VALIDATE_VECTOR(xVals);
   int nSpline = (int)xVals.getSize ();
   hVals.resize(nSpline);
   // last hVals is not needed since we have only nSpline - 1 intervals 
   hVals(nSpline-1) = 0.0;
   for (int i = 1; i < nSpline; i++)   {
      hVals(i-1) = xVals(i) - xVals(i-1);
      if (hVals(i-1) < 1e-12)   {
         cout << SX_SEPARATOR;
         cout << "Splinegrid has identical points!\nS/PHI/nX exits here!" 
              << endl;
         sxprintf ("i is %i: x(i) = %.6f, x(i-1) = %.6f\n",
               i,xVals(i),xVals(i-1));
         cout << SX_SEPARATOR;
         SX_EXIT;
      }
   }
}

// Setup coeff Matrix for Splinecoeff calculation
template <class V>
void SxCubicSpline<V>::setT (const enum SplineType sType)
{
   int nSpline = (int)xVals.getSize();
   // layout:
   // column 0: subdiagonal i+1,i
   // column 1: diagonal i,i
   // column 2: superdiagonal i,i+1
   TMat.reformat (nSpline, 3);
   for (int iRow = 0; iRow < nSpline; iRow++)   {
      // --- conditions for left bondary
      if (iRow == 0 )   {
         if (sType == Hermite || sType == HermiteNatural)   {
            // Hermite edge (f'(0)): 
            // 2h_i/3 c_i + h_i/3 c_{i+1} = ...
            TMat(iRow, 1) = 2. * hVals(iRow) / 3.;
            TMat(iRow, 2) = hVals(iRow) / 3.;
         } else if (sType == Natural || sType == NaturalHermite) {
            // cub spline: c_i = ...
            TMat(iRow,1) = 1.;
            TMat(iRow,2) = 0.;
         } else  {
            cout << "Unknown splinetype" << endl;
            SX_EXIT;
         }
      }
      // --- conditions for right bondary
      else if (iRow == (nSpline - 1)) {
         if (sType == Hermite || sType == NaturalHermite) {
            // Hermite
            TMat(iRow,1) = 2 * hVals(iRow-1) / 3.;
            TMat(iRow-1,0) = hVals(iRow-1) / 3.;
         } else if (sType == Natural || sType == HermiteNatural)  {
            // Natural
            TMat(iRow,1) = 1.;
            TMat(iRow-1,0) = 0.;
         } else  {
            cout << "Unknown splinetype" << endl;
            SX_EXIT;
         }
      }
      else {
         // Cubic Spline Textbook
         // h_{j-1}c_{j-1} + 2(h_{j-1} + h_j)c_j + h_jc_{j+1} = ...
         TMat(iRow-1,0) = hVals(iRow-1);
         TMat(iRow  ,1) = 2. * (hVals(iRow-1) + hVals(iRow));
         TMat(iRow  ,2) = hVals(iRow);
      }
   }
}

// Setup right hand side for spline coeff calculation
template <class V>
void SxCubicSpline<V>::setBeta (const enum SplineType sType)
{
   int nSpline = (int)xVals.getSize (); 
   betaVals.resize (nSpline);
   // left bondary
   // Hermite : beta(i) = [f(x_{i+1}) - f(x_i)] / h_i - f'(x_i)
   if (sType == Hermite || sType == HermiteNatural) 
      betaVals(0) = (yVals(1) - yVals(0)) / hVals(0) - leftHermite;
   // natural: beta(i) = 0.0;
   if (sType == Natural || sType == NaturalHermite) betaVals(0) = 0.;
   // right bondary
   if (sType == Natural || sType == HermiteNatural) 
      betaVals(nSpline-1) = 0.;
   if (sType == Hermite || sType == NaturalHermite) 
      betaVals(nSpline-1) 
         = rightHermite 
         - (yVals(nSpline-1) - yVals(nSpline-2)) / hVals(nSpline-2);
   // center
   for (int iRow = 1; iRow < nSpline - 1; iRow++)   {
      betaVals(iRow) 
         = 3.0 / hVals(iRow)   * (yVals(iRow+1) - yVals(iRow))
         - 3.0 / hVals(iRow-1) * (yVals(iRow) - yVals(iRow-1));
   }
}

template <class V>
int SxCubicSpline<V>::getSplineIdx (
      const typename SxCubicSpline<V>::S::Type x,
      const V& basis) const
{
   int nSpline = (int)basis.getSize ();
   // Check lower bound
   if (fabs(x-basis(0)) < 1e-10 || (x > basis(0) && x < basis(1)) ) 
      return 0;
   else if (x < basis(0))   {
      /*
      cout << SX_SEPARATOR;
      cout << "Warning!" << endl;
      cout << "Here is SxCubicSpline!" << endl;
      cout << "Value lies beyound left grid point!" << endl;
      cout << "Val: " << x << ", Netpoint: " << basis(0) << endl;
      cout << SX_SEPARATOR;
      */
      return 0;
   }
   // Check upper bound
   if (fabs(x-basis(nSpline-1)) < 1e-10 || (x < basis(nSpline-1) && x > basis(nSpline-2))) 
      return nSpline - 2;
   else if (x > basis(nSpline - 1))   {
      /*
      cout << SX_SEPARATOR;
      cout << "Warning !" << endl;
      cout << "Here is SxCubicSpline!" << endl;
      cout << "Value lies beyound right grid point!" << endl;
      cout << "Val: " << x << ", Netpoint: " << basis(nSpline - 1) << endl;
      cout << SX_SEPARATOR;
      */
      return nSpline - 2;
   }
   // intervall is half of nSpline when nSpline is even or [nSpline/2] when nSpline is odd
   int interval = nSpline / 2;
   int result = interval;
   bool found = false;
   int counter = 0;
   // find Idx via interval intersection
   while (!found) {
      if (counter > nSpline)    {
         cout << "Here is SxCubicSpline!" << endl;
         cout << "I get no spline Idx for x, is: " << x << endl;
         cout << "current idx : " << result << ", but x is not between " << basis(result)
              << " and " << basis(result+1) << endl;
         basis.print();
         SX_EXIT;
      }
      if (result + 1 >= nSpline) {
         cout << "Here is SxCubicSpline!" << endl;
         cout << "Somehow ran to upper bound!" << endl;
         cout << "I get no spline Idx for x, is: " << x << endl;
         cout << "current idx : " << result << ", but x is larger than " 
              << basis(result) << endl;
         basis.print();
         SX_EXIT;
      }
      if (x >= basis(result) && x <= basis(result+1)) found = true;
      else   {
         interval = interval / 2;
         if (interval == 0) interval = 1;
         if (x < basis(result)) result -= interval;
         else result += interval;
         counter++;
      }
   }
   return result;
}

template<class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::dady ()  
{
   ssize_t nSpline = xVals.getSize ();
   M result (nSpline,nSpline);
   result.set(0.0);
   // da_i / dy_j = d_{i,j}
   for (ssize_t i = 0; i < nSpline; i++)   {
      result(i,i) = 1.0;
   }
   return result;
}

template<class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::dbdy (
             const typename SxCubicSpline<V>::M &dA, 
             const typename SxCubicSpline<V>::M &dC)
{
   ssize_t nSpline = xVals.getSize ();
   SX_CHECK(nSpline == dA.colRef(0).getSize(),
            nSpline, dA.colRef(0).getSize());
   SX_CHECK(nSpline == dC.colRef(0).getSize(),
            nSpline, dC.colRef(0).getSize());
   SX_CHECK(nSpline == dA.row(0).getSize(),
            nSpline, dA.row(0).getSize());
   SX_CHECK(nSpline == dC.row(0).getSize(),
            nSpline, dC.row(0).getSize());
   M result (nSpline, nSpline);
   result.set(0.0);
   //db_i / dy_j = (da_{i+1} / dy_j + da_i / dy_j) / h_i - h_i / 3 * (2* dc_i / dy_j + dc_{i+1} / dy_j)
   for (ssize_t i = 0; i < nSpline - 1; i++)   {
      for (ssize_t j = 0; j < nSpline; j++)   {
         result(i,j) = (dA(i+1,j) - dA(i,j)) / hVals(i) 
                     - (2.0 * dC(i,j) + dC(i+1,j)) * hVals(i) / 3.0;
      }
   }
   return result;
}

template<class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::dcdy (const enum SplineType sType) 
{
   SX_CHECK(TMat.getSize() > 0);
   M dBeta = dbetady (sType);
   M result = symmetricTridiagonalGauss (TMat, dBeta);
   return result;
}

template<class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::dbetady (const enum SplineType sType) 
{
   ssize_t nSpline = xVals.getSize();
   M result (nSpline, nSpline);
   result.set(0.0);
   for (ssize_t i = 0; i < nSpline; i++)   {
      for (ssize_t j = 0; j < nSpline; j++)   {
         if (i == 0) {
            // Hermite edge
            if (sType == Hermite || sType == HermiteNatural)   {
               if (i+1 == j) result(i,j) += 1./hVals(i);
               if (i == j) result(i,j) -= 1./hVals(i);
            }
            // natural edge
            if (sType == Natural || sType == NaturalHermite) result(i,j) = 0.0;
         }         
         if (i == (nSpline-1)) {
            // natural edge
            if (sType == Natural || sType == HermiteNatural) result(i,j) = 0.0;
            // Hermite edge
            if (sType == Hermite || sType == NaturalHermite)   {
               if (i-1 == j) result(i,j) += 1./hVals(i-1);
               if (i == j) result(i,j) -= 1./hVals(i-1);
            }
         }
         if (i > 0 && i < (nSpline-1)) {
            if ( (i+1) == j) result(i,j) += 3.0 / hVals(i);
            if ( i == j) result(i,j) -= 3.0 / hVals(i) + 3.0 / hVals(i-1); 
            if ( (i-1) == j) result(i,j) += 3.0 / hVals(i-1);
         }
      }
   }
   return result;
}

template<class V>
typename SxCubicSpline<V>::M SxCubicSpline<V>::dddy (
      const typename SxCubicSpline<V>::M &dC) 
{
   ssize_t nSpline = xVals.getSize();
   SX_CHECK(nSpline == dC.colRef(0).getSize(),
            nSpline, dC.colRef(0).getSize());
   SX_CHECK(nSpline == dC.row(0).getSize(),
            nSpline, dC.row(0).getSize());
   M result (nSpline, nSpline);
   result.set(0.0);
   // dd_i/dy_j = (dc_{i+1}/dy_j - dc_i/dy_j) / (3h_i)
   for (ssize_t i = 0; i < nSpline - 1; i++)  {
      for (ssize_t j = 0; j < nSpline; j++)  {
         result(i,j) = (dC(i+1,j) - dC(i,j)) / (3.0 * hVals(i));
      }
   }
   return result;
}

template <class V>
void SxCubicSpline<V>::setSpline (const V &spline)
{
   SX_CHECK(spline.getSize() > 0);
   SX_CHECK(spline.getSize() % 4 == 0);
   VALIDATE_VECTOR(spline);
   polyCoeff = spline.getCopy ();
}

template <class V>
V SxCubicSpline<V>::getSpline () const
{
   SX_CHECK(polyCoeff.getSize () > 0);
   SX_CHECK(polyCoeff.getSize () % 4 == 0);
   VALIDATE_VECTOR(polyCoeff);
   return polyCoeff.getCopy ();
}

template <class V>
void SxCubicSpline<V>::mirror ()
{
   ssize_t nSpline = xVals.getSize();
   ssize_t nData   = xFit.getSize();
  
   if (nData > 0)  { 
      // Create symmetric dataset
      if (fabs(xFit(0)) < 1e-10)   {
         V xMirror(2*nData-1);
         V yMirror(2*nData-1);
         xMirror(nData-1) = xFit(0);
         yMirror(nData-1) = yFit(0);
         for (ssize_t i = 1; i < nData; i++)   {
            xMirror(i-1) = -xFit(nData-i);
            if (fitType == MirrorPlane) yMirror(i-1) = yFit(nData-i);
            else if (fitType == MirrorPoint) yMirror(i-1) = -yFit(nData-i);
            else {
               cout << SX_SEPARATOR;
               cout << "Unknown fit type in SxCubicSpline!\nS/PHI/nX exits here!" << endl;
               cout << SX_SEPARATOR;
               SX_EXIT;
            }
            xMirror(nData+i-1) = xFit(i);
            yMirror(nData+i-1) = yFit(i);
         }
         xFit = xMirror.getCopy ();
         yFit = yMirror.getCopy ();
      } else {
         V xMirror(2*nData);
         V yMirror(2*nData);
         for (ssize_t i = 0; i < nData; i++)   {
            xMirror(i) = -xFit(nData-1-i);
            if (fitType == MirrorPlane) yMirror(i) = yFit(nData-1-i);
            else if (fitType == MirrorPoint) yMirror(i) = -yFit(nData-1-i);
            else {
               cout << SX_SEPARATOR;
               cout << "Unknown fit type in SxCubicSpline!\nS/PHI/nX exits here!" << endl;
               cout << SX_SEPARATOR;
               SX_EXIT;
            }
            xMirror(nData+i) = xFit(i);
            yMirror(nData+i) = yFit(i);
         }
         xFit = xMirror.getCopy ();
         yFit = yMirror.getCopy ();
      }
   }

   // Mirror Basis
   if (nSpline > 0)  {
      if (fabs(xVals(0)) < 1e-10) {
         V xMirror(2*nSpline-1);
         xMirror(nSpline-1) = xVals(0);
         // Create symmetric spline mesh
         for (ssize_t i = 1; i < nSpline; i++)  {
            xMirror(i-1) = -xVals(nSpline-i);
            xMirror(nSpline+i-1) = xVals(i);
         }
         xVals = xMirror.getCopy ();
         yVals.resize(xVals.getSize());
         setH();

      } else {
         V xMirror(2*nSpline);
         // Create symmetric spline mesh
         for (ssize_t i = 0; i < nSpline; i++)  {
            xMirror(i) = -xVals(nSpline-1-i);
            xMirror(nSpline+i) = xVals(i);
         }
         xVals = xMirror.getCopy ();
         yVals.resize(xVals.getSize());
         setH();
      }
   }
}

template <class V>
void SxCubicSpline<V>::deMirror ()
{
   // DeMirror Vals
   // Backup
   ssize_t nSpline = xVals.getSize();   
   V xValsBackup = xVals.getCopy();
   V yValsBackup = yVals.getCopy();
   if ((nSpline & 1) == 1)   {
      nSpline = (nSpline + 1) / 2;
      xVals.resize(nSpline);
      yVals.resize(nSpline);
      for (ssize_t i = 0; i < nSpline; i++)   {
         xVals(i) = xValsBackup(nSpline-1+i);
         yVals(i) = yValsBackup(nSpline-1+i);
      }
   }
   else   {
      nSpline = nSpline / 2;
      xVals.resize(nSpline);
      yVals.resize(nSpline);
      for (ssize_t i = 0; i < nSpline; i++)   {
         xVals(i) = xValsBackup(nSpline+i);
         yVals(i) = yValsBackup(nSpline+i);
      }
   }

   //DeMirror Datapoints
   ssize_t nData = xFit.getSize();
   V xFitBackup = xFit.getCopy();
   V yFitBackup = yFit.getCopy();
   nData = nData / 2;
   xFit.resize(nData);
   yFit.resize(nData);
   for (ssize_t i = 0; i < nData; i++)   {
      xFit(i) = xFitBackup(nData+i);
      yFit(i) = yFitBackup(nData+i);
   }
}

template <class V>
SxArray<SxList<int> > SxCubicSpline<V>::getSplineExtend (SxArray<int> &pointsPerSpline)
{
   ssize_t nSpline = xVals.getSize();
   ssize_t nData = xFit.getSize();
   ssize_t minPoints = 1;
   if (nData < minPoints)  { 
      cout << "Too few points in spline fit." << endl;
      SX_EXIT;
   }

   SxArray<SxList<int> > result(nSpline);

   for (ssize_t iSpline = 0; iSpline < nSpline; iSpline++)  {
      result(iSpline).resize(0);
      ssize_t nPoints = pointsPerSpline(iSpline);
      for (ssize_t d = 1; nPoints < minPoints; ++d)  {
         if (iSpline - d >= 0)  {
            ssize_t jSpline = iSpline - d;
            result(iSpline) << int(jSpline);
            nPoints += pointsPerSpline(jSpline);
         }
         if (iSpline + d < nSpline)  {
            ssize_t jSpline = iSpline + d;
            result(iSpline) << int(jSpline);
            nPoints += pointsPerSpline(jSpline);
         }
      }
   }

   return result;    

}
template <class V>
SxVector<Int> SxCubicSpline<V>::getPointsPerSpline (const V& basis, const V& xData)
{

   ssize_t nSpline = basis.getSize();
   ssize_t nData = xData.getSize ();
   SxVector<Int> result (nSpline);
   result.set(0);
   for (ssize_t iVal = 0; iVal < nData; iVal++) {
      int iSpline = getSplineIdx(xData(iVal),basis);
      result(iSpline)++;
   }

   return result;
}

template <class V>
V SxCubicSpline<V>::cleanZeroDepths(const V& xData, const SxVector<Int> &pointsPerSpline)
{
   SX_CHECK(xData.getSize() == pointsPerSpline.getSize (), xData.getSize(), pointsPerSpline.getSize ());
   // keep outer right Point
   ssize_t fullSplines = 1;
   ssize_t nSplines = pointsPerSpline.getSize();

   for (ssize_t iSpline = 0; iSpline < nSplines - 1; iSpline++)  
      if (pointsPerSpline(iSpline) != 0) fullSplines++;

   if (fullSplines == nSplines) return xData;
   else  {
      V result(fullSplines);

      ssize_t iCounter = 0;
      for (ssize_t iSpline = 0; iSpline < nSplines - 1; iSpline++)  {
         if  (pointsPerSpline(iSpline) != 0)  {
            result(iCounter) = xData(iSpline);
            iCounter++;
         }
      }
      // keep Outer right Point
      result(iCounter) = xData(nSplines - 1);
      return result;
   }
}

template <class V>
V SxCubicSpline<V>::symmetricTridiagonalGauss (
      const M &Mat,
      const V &rhs)
{
   SX_CHECK (Mat.nCols () == 3, Mat.nCols ());
   SX_CHECK (Mat.nRows () == rhs.nRows (), Mat.nRows (), rhs.nRows ());
   int dim = (int)Mat.nRows ();
   V diag2, b;
   diag2.copy (Mat.colRef (1));
   b = rhs.transpose (); // improve data locality...
   b.reshape (rhs.nCols (), dim);
   V diag1 (dim-1), diag3 (dim - 1);

   for (int i = 0; i < dim - 1; i++)  {
      diag1 (i) = Mat(i,0);
      diag3 (i) = Mat(i,2);
   }
  
  // Overview:
  // ---------  
  // diag2(0) diag1(0) 0        0
  // diag3(0) diag2(1) diag1(0) 0
  // 0        ...
  // ...
  // 0        0        0        diag3(n-2) diag2(n-1) diag1(n-1)
  // 0        0        0        0          diag3(n-1) diag2 (n)
  
  //Mat.print(true); 
  //cout << "Before Gaussian Elimination: " << endl; 
  //print (diag1, diag2, diag3, b);
  	
  // Gaussian Elimination
  for (int i = 0; i < dim - 1; ++i)  {
     double x = diag3 (i) / diag2(i);
     diag2(i + 1) -= diag1 (i) * x;
     for (int j = 0; j < b.nRows (); j++) // b is rhs transposed
        b(j,i + 1) -= b(j,i) * x;
     diag3(i) = 0;
  }

  //cout << "Gaussian Elimination: " << endl;
  //print (diag1, diag2, diag3, b);
  

  // Backward substitution
  {
     double x = 1. / diag2(dim-1);
     for (int j = 0; j < b.nRows (); j++)
        b(j,dim - 1) *= x;
  }
  //diag2(dim-1) = 1;
  for (int i = dim - 2; i >= 0; --i)  {
     double x = 1. / diag2(i);
     for (int j = 0; j < b.nRows (); j++)  {
        b(j,i) -= diag1(i) * b(j,i + 1);
        b(j,i) *= x;
     }
     //diag1(i) = 0;
     //diag2(i) = 1;
  } 
  //cout << "Backward substitution: " << endl;
  //print (diag1, diag2, diag3, b);

  return b.transpose ();
}

template <class V>
void SxCubicSpline<V>::print (const V &diag1, const V &diag2, const V &diag3, const V &b)
{
  int dim = diag2.getSize ();
  
   // TESTOUTPUT
  for (int i = 0; i < dim; ++i) {
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

