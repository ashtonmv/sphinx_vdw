// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#include <SxRadRBasis.h>
#include <SxRadBasis.h>
#include <SxRadGBasis.h>
#include <SxGBasis.h>


SxRadRBasis::SxRadRBasis ()
{
   //empty
}

SxRadRBasis::~SxRadRBasis ()
{
   //empty/
}

SxRadRBasis::SxRadRBasis (double rMin, double rMax, int nPoints, enum gridType mode)
{
   set (rMin, rMax, nPoints, mode);
}

SxRadRBasis::SxRadRBasis (const SxDiracVec<Double> &basis)
{
   set (basis);
}

SxRadRBasis::SxRadRBasis (const SxDiracVec<TBasisType> &drIn, double factor, double r0)
{
   set (drIn, factor, r0);
}
   
void SxRadRBasis::set (double rMin, double rMax, int nPoints, enum gridType mode)
{
   r.resize(nPoints);
   dr.resize(nPoints);
   if (mode == Linear)   {
      // linear grid
      double deltaR = (rMax - rMin) / double (nPoints - 1);
      dr.set(deltaR);
      for(int i = 0; i < nPoints; i++)   {
         r(i) = rMin + double(i) * dr(i);
      }
   } else if (mode == Quadratic)   {
      // quadratic grid
      double deltaR = (rMax - rMin) / sqrt(1.0*(nPoints - 1.));
      for(int i = 0; i < nPoints; i++)   {
         r(i) = rMin + sqrt(1.0*i) * deltaR;
         if(i > 0) dr(i-1) = r(i) - r(i-1);
      }
   } else if (mode == Mixed)   {
      // quadratic grid + linear grid    1:2 
      double lenght = (rMax - rMin);
      double quadLenght = lenght / 3.0;
      double deltaRQuad = quadLenght / sqrt(nPoints/3.0 - 1.);
      double deltaRLin = lenght / (nPoints - 1);
      for(int i = 0; i < nPoints; i++)   {
         if (i < nPoints/3.) r(i) = rMin + sqrt(1.0*i) * deltaRQuad;
         else r(i) = rMin + i * deltaRLin;
         if(i > 0) dr(i-1) = r(i) - r(i-1);
      }
   } else if (mode == Hyperbel)   {
      // hyperbel grid
      double diff = rMax - rMin;
      double b = (nPoints - 1.) / (2.*diff);
      double a = 1. / (2.*diff);
      for(int i = 0; i < nPoints; i++)   {
         r(i) = rMin + 1.*i / (b+a*i);
         if(i > 0) dr(i-1) = r(i) - r(i-1);
      }
   } else if (mode == Logarithmic)  {
      SX_CHECK (fabs(rMin) > 1e-12, rMin); 
      double logDr = (log(rMax) - log(rMin)) / (nPoints - 1.0);
      for(int i = 0; i < nPoints; i++)   {
         r(i) = rMin * exp(i * logDr);
         if(i > 0) dr(i-1) = r(i) - r(i-1);
      }
   } else {
      cout << "UNKNOWN MODE IN SXRADGBASIS!" << endl;
      SX_EXIT;
   }

   // last dr is equal zero conventionally
   dr(nPoints-1) = 0.0;

   VALIDATE_VECTOR(r);
   VALIDATE_VECTOR(dr);
}

void SxRadRBasis::set (const SxDiracVec<TBasisType> &basis)
{
   r = basis.getCopy();
   int nPoints = (int)r.getSize();

   dr.resize(nPoints);
   dr(nPoints-1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      dr(i-1) = r(i) - r(i-1);
   }

   VALIDATE_VECTOR(dr);
}

void SxRadRBasis::set (const SxDiracVec<TBasisType> &drIn, double factor, double r0)
{
   int nPoints = (int)drIn.getSize();
   r.resize(nPoints);
   r(0) = r0;
   dr.resize(nPoints);
   dr(nPoints - 1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      if (drIn(i-1) > 1e-8)
         dr(i-1) = factor / drIn (i-1) /(1.0 * nPoints);
      else {
       cout << "i = " << i << " is " << dr(i-1) << endl;  
         SX_EXIT;
      }
      r(i) = r(i-1) + dr(i-1);
   }
   dr.print();
   VALIDATE_VECTOR(r);
}

const SxDiracVec<SxRadRBasis::TBasisType> & SxRadRBasis::getRadRFunc () const 
{ 
   return r;
}

const SxDiracVec<SxRadRBasis::TBasisType> & SxRadRBasis::getRadDR () const 
{ 
   return dr;
}

SxDiracVec<SxRadRBasis::TBasisType> SxRadRBasis::identity (const SxRadRBasis *basisPtr,
      const TPsi &vec) const
{
   // --- verify that this is really the identity projector, i.e.
   //     that the provided basis is the same as the vector's basis
   SX_CHECK (vec.getBasisPtr() == this);
   SX_CHECK (basisPtr == this);
   return vec;
}

SxDiracVec<TRadBasisType> SxRadRBasis::toRadBasis (const SxRadBasis *radBasisPtr,
      const TPsi &vec) const
{
   SX_CLOCK(Timer::radR2Rad);

   // Get aux data
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   const SxDiracVec<Double> &radFunc = radBasisPtr->radFunc(is);
   ssize_t newDim = radFunc.getSize ();
   ssize_t dim = r.getSize ();
   SxDiracVec<TBasisType> x = 1.0 * r,
                          y = 1.0 * vec;
   SxCubicSpline<TPsi> spline;

   // Check for boundary
   if (r(dim-1) < radFunc(newDim-1))  {
      double xLast = r(dim-1);
      double yLast = vec(dim-1);
      double xPreLast = r(dim-2);
      double yPreLast = vec(dim-2);
      double slope = (yLast - yPreLast)/(xLast - xPreLast);

      double xLastNew = radFunc(newDim-1);
      double slopeRight = 0.0;

     if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  { 
        // force to zero
        x.resize(dim+1,true);
        y.resize(dim+1,true);
        x(dim) = xLastNew;
        y(dim) = 0.0;
     } else {
        // exponential decay
        double alpha = slope / yLast;
        int nPoints = 100;
        x.resize(dim+nPoints,true);
        y.resize(dim+nPoints,true);
        double delta = (xLastNew - xLast) / nPoints;
        for (int i = 0; i < nPoints; i++)  {
           x(dim + i) = xLast + (i+1) * delta;
           y(dim + i) = yLast * exp(alpha * (i+1) * delta);
        }
        slopeRight = y(dim+nPoints-1) * slope;
     }
     if (l < 2) {
        spline = SxCubicSpline<SxDiracVec<TBasisType> > (
              x,y, 
              SxCubicSpline<SxDiracVec<TBasisType> >::NaturalHermite,
              slopeRight);
     } else  {
        spline = SxCubicSpline<SxDiracVec<TBasisType> > (
              x,y, 
              SxCubicSpline<SxDiracVec<TBasisType> >::Hermite,
              0.0,
              slopeRight);
     }
   } else  {
      if (l < 2) 
         spline = SxCubicSpline<SxDiracVec<Double> > (
               x, 
               y, 
               SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
      else
         spline = SxCubicSpline<SxDiracVec<Double> > (
               x, 
               y, 
               SxCubicSpline<SxDiracVec<Double> >::Hermite);
   }

   TPsi result = spline.getY(radFunc);

   // set aux data
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;
   result.setBasis (radBasisPtr);

   return result;
}

SxDiracVec<TRadBasisType> SxRadRBasis::toRadGBasis (const SxRadGBasis *radGBasisPtr,
      const TPsi &vec) const
{ 
   SX_CLOCK (Timer::radR2RadG);
   // Get aux data
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;
   
   const SxDiracVec<Double> &radGFunc = radGBasisPtr->getRadGFunc();
   int ng = (int)radGFunc.getSize ();
   SxDiracVec<TRadGBasisType> result(ng);


   SxCubicSpline<SxDiracVec<Double> > spline;
   if (vec.getSize() == 4 * getNElements ())  {
      spline = SxCubicSpline<SxDiracVec<Double> > (r, vec);
   } else if (vec.getSize() == getNElements ()) {
      if (l < 2)  {
         spline = SxCubicSpline<SxDiracVec<Double> > (r, vec,
               SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
      } else {
         spline = SxCubicSpline<SxDiracVec<Double> > (r, vec,
               SxCubicSpline<SxDiracVec<Double> >::Hermite);
      }
   } else  {
      cout << "Incompatible sizes between basis and vector!" << endl;
      SX_EXIT;
   }


   double drDense = 0.001 * TWO_PI / radGFunc(ng-1);
   double rMax = r(r.getSize()-1);
   int nRadPoints = max(int(rMax/drDense) + 1, (int)r.getSize ());
   nRadPoints = min(nRadPoints, 10 * (int)r.getSize ());
   SxRadRBasis denseRadR (r(0), rMax, nRadPoints);

   SxDiracVec<Double> denseVec = spline.getY(denseRadR.r);


   for (int ig = 0; ig < ng; ig++)  {
      SxDiracVec<Double> jl = jsb (l, radGFunc(ig) * denseRadR.r);
      // Psi(G) = int(jl(r*G) * Psi(r) * r^2dr)
      result(ig) = denseRadR.integrate(jl * denseVec, false);
   }
   result *= sqrt(2./PI);

   // set aux data
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;
   result.setBasis (radGBasisPtr);

   return result;
}

double SxRadRBasis::integrate (const TPsi &integrand, bool useSpline) const
{
   SX_CLOCK(Timer::integrateR);
   SX_CHECK(integrand.getSize() > 0);
   VALIDATE_VECTOR(integrand);

   // general Simpson integration sheme 3rd Order
   int dim = (int)r.getSize ();
   double result = 0.0;
   TPsi f = integrand * r * r;
   if (!useSpline)  {
      int i = 0;
      // general Simpson integration sheme 3rd Order
      while (i < dim - 1)  {
         if (i < dim - 3)  {
            result += simpsonTerm(r(i),r(i+1),r(i+2),r(i+3),f(i),f(i+1),f(i+2),f(i+3));
            i += 3;
         } else if (i == dim - 3)  {
            result += simpsonTerm(r(i),r(i+1),r(i+2),f(i),f(i+1),f(i+2));
            i += 2;
         } else if (i == dim - 2)  {
            result += simpsonTerm(r(i),r(i+1),f(i),f(i+1));
            i += 1;
         }
      }
   } else {
      TPsi coeffs = toSpline(f);
      // Integration scheme for cubic Splines S = a x^3 + b x^2 + c x + d
      // Setup Coefficients (integrated analytically S dx) a-> 1/4 a, b -> 1/3 b, c -> 1/2 c, d-> d
      for (int i = 0; i < dim - 1; i++)  {
         double h = r(i+1)-r(i); 
         result += coeffs(4*i+0) * h
                 + coeffs(4*i+1) * h*h / 2.0
                 + coeffs(4*i+2) * h*h*h / 3.0
                 + coeffs(4*i+3) * h*h*h*h / 4.0;
      }
   }
   
   return result;
}

double SxRadRBasis::simpsonTerm (double x1, double x2, double x3, double x4, 
      double f1, double f2, double f3, double f4) const
{
   double T1,T2,T3,T4, result;

   T1 = f1*(6.*x2*x3 + (2.*(x1-x2-x3)+x4)*x4-x1*(4.*x2+4.*x3-3.*x1)) / ((x2-x1)*(x1-x3));
   T2 = f2*((x1-x4)*(x1-x4)*(x1-2.*x3+x4))/((x2-x1)*(x2-x3)*(x2-x4));
   T3 = f3*((x1-x4)*(x1-x4)*(x1-2.*x2+x4))/((x1-x3)*(x2-x3)*(x3-x4));
   T4 = f4*(6.*x2*x3 + (2.*(x4-x2-x3)+x1)*x1-x4*(4.*x2+4.*x3-3.*x4)) / ((x4-x2)*(x3-x4));;

   result = (x1-x4)/12. * (T1+T2+T3+T4);

   return result;

}

double SxRadRBasis::simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const
{
   double T1,T2,T3, result;

   T1 = f1*(x2-x3)*(2.*x1-3.*x2+x3);
   T2 = f2*(x1-x3)*(x1-x3);
   T3 = f3*(x2-x1)*(x1-3.*x2+2.*x3);

   result = (x3-x1)/(6.*(x1-x2)*(x2-x3)) * (T1+T2+T3);

   return result;
}

double SxRadRBasis::simpsonTerm (double x1, double x2, double f1, double f2) const
{
   double result = 0.5 * (f1+f2) * (x2-x1);

   return result;
}


SxDiracVec<Double> SxRadRBasis::jsb (int l, const SxDiracVec<Double> &z) 
{

   SX_CHECK (z.getSize () > 0, z.getSize());

   int zSize = (int)z.getSize ();
   SxDiracVec<TReal8> vec(zSize);
   // --- use generic jsb
#pragma omp parallel for
   for (int i = 0; i < zSize; ++i)
      vec(i) = SxYlm::jsb(l, z(i));
   VALIDATE_VECTOR (vec);
   return vec;
}

SxDiracVec<Double> SxRadRBasis::toSpline (const SxDiracVec<Double> &vec) const
{
   SX_CLOCK(Timer::VecToSplineR);
   // get auxdata
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   SxCubicSpline<SxDiracVec<Double> > spline;
   if (l < 2) 
      spline = SxCubicSpline<SxDiracVec<Double> > (r, vec,
            SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
   else 
      spline = SxCubicSpline<SxDiracVec<Double> > (r, vec,
            SxCubicSpline<SxDiracVec<Double> >::Hermite);

   SxDiracVec<Double> result = spline.getSpline ();

   // set Basis and auxData
   result.setBasis(this);
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;

   return result;
}

SxDiracVec<Double> SxRadRBasis::toVec (const SxDiracVec<Double> &vec) const
{
   SX_CLOCK(Timer::SplineToVecR);
   // get auxdata
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   SxCubicSpline<SxDiracVec<Double> > spline (r, vec);

   SxDiracVec<Double> result = spline.getY (r);

   // set Basis and auxData
   result.setBasis(this);
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;

   return result;
}
