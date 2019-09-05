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

#include <SxRadGBasis.h>
#include <SxRadBasis.h>
#include <SxRadRBasis.h>
#include <SxGBasis.h>


SxRadGBasis::SxRadGBasis ()
{
   //empty
}

SxRadGBasis::~SxRadGBasis ()
{
   //empty
}

SxRadGBasis::SxRadGBasis (double gMin, double gMax, int nPoints, enum gridType mode)
{
   set (gMin, gMax, nPoints, mode);
}

SxRadGBasis::SxRadGBasis (const SxDiracVec<Double> &basis)
{
   set (basis);
}
SxRadGBasis::SxRadGBasis (const SxDiracVec<Double> &dgIn, double factor, double g0)
{
   set (dgIn, factor, g0);
}

void SxRadGBasis::set (double gMin, double gMax, int nPoints, enum gridType mode)
{
   g.resize(nPoints);
   dg.resize(nPoints);
   if (mode == Linear)   {
      // linear grid
      double deltaG = (gMax - gMin) / double (nPoints - 1);
      dg.set(deltaG);
      for(int i = 0; i < nPoints; i++)   {
         g(i) = gMin + double(i) * dg(i);
      }
   } else if (mode == Quadratic)   {
      // quadratic grid
      double deltaG = (gMax - gMin) / sqrt(1.0*(nPoints - 1.));
      for(int i = 0; i < nPoints; i++)   {
         g(i) = gMin + sqrt(1.0*i) * deltaG;
         if(i > 0) dg(i-1) = g(i) - g(i-1);
      }
   } else if (mode == Mixed)   {
      // quadratic grid + linear grid    1:2 
      double lenght = (gMax - gMin);
      double quadLenght = lenght / 3.0;
      double deltaGQuad = quadLenght / sqrt(nPoints/3.0 - 1.);
      double deltaGLin = lenght / (nPoints - 1);
      for(int i = 0; i < nPoints; i++)   {
         if (i < nPoints/3.) g(i) = gMin + sqrt(1.0*i) * deltaGQuad;
         else g(i) = gMin + i * deltaGLin;
         if(i > 0) dg(i-1) = g(i) - g(i-1);
      }
   } else if (mode == Hyperbel)   {
      // hyperbel grid
      double diff = gMax - gMin;
      double b = (nPoints - 1.) / (2.*diff);
      double a = 1. / (2.*diff);
      for(int i = 0; i < nPoints; i++)   {
         g(i) = gMin + 1.*i / (b+a*i);
         if(i > 0) dg(i-1) = g(i) - g(i-1);
      }
   } else if (mode == Logarithmic)  {
      SX_CHECK (fabs(gMin) > 1e-12, gMin); 
      double logDg = (log(gMax) - log(gMin)) / (nPoints - 1.0);
      for(int i = 0; i < nPoints; i++)   {
         g(i) = gMin * exp(i * logDg);
         if(i > 0) dg(i-1) = g(i) - g(i-1);
      }
   } else {
      cout << "UNKNOWN MODE IN SXRADGBASIS!" << endl;
      SX_EXIT;
   }

   // last dg is equal zero conventionally
   dg(nPoints-1) = 0.0;

   VALIDATE_VECTOR(g);
   VALIDATE_VECTOR(dg);
}

void SxRadGBasis::set (const SxDiracVec<Double> &basis)
{
   g = basis.getCopy();
   int nPoints = (int)g.getSize();

   dg.resize(nPoints);
   dg(nPoints-1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      dg(i-1) = g(i) - g(i-1);
   }

   VALIDATE_VECTOR(dg);
}

void SxRadGBasis::set (const SxDiracVec<Double> &dgIn, double factor, double g0)
{
   int nPoints = (int)dgIn.getSize();
   g.resize(nPoints);
   g(0) = g0;
   dg.resize(nPoints);
   dg(nPoints - 1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      if (dgIn(i-1) > 1e-8)
         dg(i-1) = factor / dgIn (i-1) / (1.0 * nPoints);
      else SX_EXIT;
      g(i) = g(i-1) + dg(i-1);
   }

   VALIDATE_VECTOR(g);
}

const SxDiracVec<Double> &SxRadGBasis::getRadGFunc () const 
{ 
   return g;
}

const SxDiracVec<Double> &SxRadGBasis::getRadDG () const 
{ 
   return dg;
}

SxDiracVec<SxRadGBasis::TBasisType> SxRadGBasis::identity (const SxRadGBasis *basisPtr,
                                                           const TPsi &vec) const
{
   // --- verify that this is really the identity projector, i.e.
   //     that the provided basis is the same as the vector's basis
   SX_CHECK (vec.getBasisPtr() == this);
   SX_CHECK (basisPtr == this);
   return vec;
}

SxDiracVec<TRadBasisType> SxRadGBasis::toRadBasis (const SxRadBasis *radBasisPtr,
                                                   const TPsi &vec) const
{ 
   SX_CLOCK (Timer::radG2Rad);
   SX_CHECK(radBasisPtr);
   
   int iSpecies = vec.handle->auxData.is;
   int iAtom    = vec.handle->auxData.ia;
   int n        = vec.handle->auxData.n;
   int l        = vec.handle->auxData.l;
   int m        = vec.handle->auxData.m;
   const SxDiracVec<Double> &rAbs = radBasisPtr->radFunc(iSpecies);
   int dimRad = (int)rAbs.getSize();
   SxDiracVec<TRadBasisType> result (dimRad);
   for (int ir = 0; ir < dimRad; ir++)   {
      SxDiracVec<Double> jl = jsb(l,g * rAbs(ir));
      // Psi(r) = int(jl(r*G) * Psi(G) * G^2dG)
      result(ir) = integrate(jl * vec, true);
   }
   result *= sqrt(2./PI);
   result.setBasis(*radBasisPtr);
   result.handle->auxData.is = iSpecies;
   result.handle->auxData.ia = iAtom;
   result.handle->auxData.n = n;
   result.handle->auxData.l = l;
   result.handle->auxData.m = m;

   return result;
}

SxDiracVec<TRadBasisType> SxRadGBasis::toRadRBasis (const SxRadRBasis *radRBasisPtr,
                                                   const TPsi &vec) const
{ 
   SX_CLOCK (Timer::radG2RadR);
   
   int iSpecies = vec.handle->auxData.is;
   int iAtom    = vec.handle->auxData.ia;
   int n        = vec.handle->auxData.n;
   int l        = vec.handle->auxData.l;
   int m        = vec.handle->auxData.m;
   const SxDiracVec<Double> &r = radRBasisPtr->getRadRFunc();
   int dimRad = (int)r.getSize();
   SxDiracVec<TRadBasisType> result (dimRad);

   SxCubicSpline<SxDiracVec<Double> > spline;
   if (vec.getSize() == 4 * getNElements ())  {
      spline = SxCubicSpline<SxDiracVec<Double> > (g, vec);
   } else if (vec.getSize() == getNElements ()) {
      if (l < 2)  {
         spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
               SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
      } else  {
         spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
               SxCubicSpline<SxDiracVec<Double> >::Hermite);
      }
   } else  {
      cout << "Incompatible sizes between basis and vector!" << endl;
      SX_EXIT;
   }


   SxRadGBasis denseRadG (g(0),g(g.getSize()-1), 10 * (int)g.getSize ());

   SxDiracVec<Double> denseVec = spline.getY(denseRadG.g);
   
   for (int ir = 0; ir < dimRad; ir++)   {
      SxDiracVec<Double> jl = jsb(l,denseRadG.g * r(ir));
      // Psi(r) = int(jl(r*G) * Psi(G) * G^2dG)
      result(ir) = denseRadG.integrate(jl * denseVec, false);
   }
   result *= sqrt(2./PI);

   //set aux data
   result.setBasis(*radRBasisPtr);
   result.handle->auxData.is = iSpecies;
   result.handle->auxData.ia = iAtom;
   result.handle->auxData.n = n;
   result.handle->auxData.l = l;
   result.handle->auxData.m = m;

   return result;
}


SxDiracVec<TGBasisType> SxRadGBasis::toGBasis (const SxGBasis *gBasisPtr,
                                               const TPsi &vec) const
{
   SX_CHECK (vec.getSize() > 0);
   SX_CLOCK (Timer::radG2Gk);
   SX_CHECK (vec.getBasisPtr () == this);
   SX_CHECK (gBasisPtr->structPtr);
   SX_CHECK (gBasisPtr->structPtr->cell.volume > 0.,
             gBasisPtr->structPtr->cell.volume);

   // Get auxData
   int iSpecies = vec.handle->auxData.is;
   int iAtom    = vec.handle->auxData.ia;
   int n        = vec.handle->auxData.n;
   int l        = vec.handle->auxData.l;
   int m        = vec.handle->auxData.m;


   // Interpolate to the gBasis points
   int lastRad = (int)g.getSize()-1;
   SxDiracVec<Double> gVec = sqrt(gBasisPtr->g2);
   int lastG = (int)gVec.getSize() - 1;
   if (g(0) > gVec(0)) {
      cout << "Cannot interpolate to outer left region." << endl;
      SX_EXIT;
   }
   if (g(lastRad) + 1e-10 < gVec(lastG) ) {
      cout << "Cannot interpolate to outer right region." << endl;
      sxprintf("max. g of radial G basis: %.15f\n", g(lastRad));
      sxprintf("max. g of Gk basis:       %.15f\n", gVec(lastG));
      SX_EXIT;
   }

   SxCubicSpline<SxDiracVec<Double> > spline;
   if (vec.getSize() == 4 * getNElements ())  {
      spline = SxCubicSpline<SxDiracVec<Double> > (g, vec);
   } else if (vec.getSize() == getNElements ()) {
      if (l < 2)  {
         spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
               SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
      } else  {
         spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
               SxCubicSpline<SxDiracVec<Double> >::Hermite);
      }
   } else  {
      cout << "Incompatible sizes between basis and vector!" << endl;
      SX_EXIT;
   }

   SX_START_TIMER(Timer::splineGetY);
   SxDiracVec<TGBasisType> result = spline.getY(gVec);
   SX_STOP_TIMER(Timer::splineGetY);

   // Multiply spherical harmonics if (l,m) is physical
   if (m != NONE_M)   { 
      result *= SxYlm::getYlmNormFactor(l,m) * gBasisPtr->getYlm(l,m);
   }
   
   // Shift to atom iAtom if specified
   if (iAtom != -1)  {
      SX_CLOCK (Timer::phaseFactors);
      result *= gBasisPtr->getPhaseFactors(iSpecies,iAtom);
   }
   // G+k normalization factor
   result *= sqrt(TWO_PI*TWO_PI*TWO_PI/gBasisPtr->structPtr->cell.volume);
    
   // set aux data
   result.setBasis (gBasisPtr);
   result.handle->auxData.is = iSpecies;
   result.handle->auxData.ia = iAtom;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;
   
   return result;
}

double SxRadGBasis::integrate (const TPsi &integrand, bool useSpline) const
{
   SX_CLOCK(Timer::integrateG);
   SX_CHECK(integrand.getSize() > 0);
   VALIDATE_VECTOR(integrand);
   // rectangular integration
   //return (integrand * g * g * dg).sum ();

   int dim = (int)g.getSize ();
   double result = 0.0;
   TPsi f = integrand * g * g;
   if (!useSpline)  {
      int i = 0;
      // general Simpson integration sheme 3rd Order
      while (i < dim - 1)  {
         if (i < dim - 3)  {
            result += simpsonTerm(g(i),g(i+1),g(i+2),g(i+3),f(i),f(i+1),f(i+2),f(i+3));
            i += 3;
         } else if (i == dim - 3)  {
            result += simpsonTerm(g(i),g(i+1),g(i+2),f(i),f(i+1),f(i+2));
            i += 2;
         } else if (i == dim - 2)  {
            result += simpsonTerm(g(i),g(i+1),f(i),f(i+1));
            i += 1;
         }
      }
   } else {
      TPsi coeffs = toSpline(f);
      // Integration scheme for cubic Splines S = a x^3 + b x^2 + c x + d
      // Setup Coefficients (integrated analytically S dx) a-> 1/4 a, b -> 1/3 b, c -> 1/2 c, d-> d
      for (int i = 0; i < dim - 1; i++)  {
         double h = g(i+1)-g(i); 
         result += coeffs(4*i+0) * h
                 + coeffs(4*i+1) * h*h / 2.0
                 + coeffs(4*i+2) * h*h*h / 3.0
                 + coeffs(4*i+3) * h*h*h*h / 4.0;
      }
   }
   return result;
}

double SxRadGBasis::simpsonTerm (double x1, double x2, double x3, double x4, 
      double f1, double f2, double f3, double f4) const
{
   double T1,T2,T3,T4, result;

   T1 = f1*(6.*x2*x3 + (2*(x1-x2-x3)+x4)*x4-x1*(4.*x2+4.*x3-3.*x1)) / ((x2-x1)*(x1-x3));
   T2 = f2*((x1-x4)*(x1-x4)*(x1-2.*x3+x4))/((x2-x1)*(x2-x3)*(x2-x4));
   T3 = f3*((x1-x4)*(x1-x4)*(x1-2.*x2+x4))/((x1-x3)*(x2-x3)*(x3-x4));
   T4 = f4*(6.*x2*x3 + (2*(x4-x2-x3)+x1)*x1-x4*(4.*x2+4.*x3-3.*x4)) / ((x4-x2)*(x3-x4));;

   result = (x1-x4)/12. * (T1+T2+T3+T4);

   return result;

}

double SxRadGBasis::simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const
{
   double T1,T2,T3, result;

   T1 = f1*(x2-x3)*(2.*x1-3.*x2+x3);
   T2 = f2*(x1-x3)*(x1-x3);
   T3 = f3*(x2-x1)*(x1-3.*x2+2.*x3);

   result = (x3-x1)/(6.*(x1-x2)*(x2-x3)) * (T1+T2+T3);

   return result;
}

double SxRadGBasis::simpsonTerm (double x1, double x2, double f1, double f2) const
{
   double result = 0.5 * (f1+f2) * (x2-x1);

   return result;
}


SxDiracVec<Double> SxRadGBasis::jsb (int l, const SxDiracVec<Double> &z)
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

SxDiracVec<Double> SxRadGBasis::toSpline (const SxDiracVec<Double> &vec) const
{
   SX_CLOCK(Timer::VecToSplineG);
   // get auxdata
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   SxCubicSpline<SxDiracVec<Double> > spline;
   if (l < 2) 
      spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
            SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
   else 
      spline = SxCubicSpline<SxDiracVec<Double> > (g, vec,
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

SxDiracVec<Double> SxRadGBasis::toVec (const SxDiracVec<Double> &vec) const
{
   SX_CLOCK(Timer::SplineToVecG);
   // get auxdata
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   SxCubicSpline<SxDiracVec<Double> > spline (g, vec);

   SxDiracVec<Double> result = spline.getY (g);

   // set Basis and auxData
   result.setBasis(this);
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;

   return result;
}
