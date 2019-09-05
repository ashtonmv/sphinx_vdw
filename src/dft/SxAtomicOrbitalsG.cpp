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

#include <SxAtomicOrbitalsG.h>
SxAtomicOrbitalsG::SxAtomicOrbitalsG ()
{
   // empty
}

SxAtomicOrbitalsG::SxAtomicOrbitalsG (
      const SxArray<SxArray<SxDiracVec<Double> > > &in,
      SxConstPtr<SxRadGBasis> radGBasisPtrIn,
      bool splineRepIn)
{
   // --- check dimensions
#  ifndef NDEBUG
   for (int is = 0; is < in.getSize(); is++)  {
      int nOrbTypes = (int)in(is).getSize();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         SX_CHECK(in(is)(iot).handle->auxData.is >= 0, in(is)(iot).handle->auxData.is);
         SX_CHECK(in(is)(iot).handle->auxData.n >= 0, in(is)(iot).handle->auxData.n);
         SX_CHECK(in(is)(iot).handle->auxData.l >= 0, in(is)(iot).handle->auxData.l);
         if (splineRepIn) {
            SX_CHECK (in(is)(iot).getSize() == 4 * radGBasisPtrIn->getNElements(),
                      in(is)(iot).getSize(), 4 * radGBasisPtrIn->getNElements());
         }
         else  {   
            SX_CHECK (in(is)(iot).getSize() == radGBasisPtrIn->getNElements(),
                      in(is)(iot).getSize(), radGBasisPtrIn->getNElements());
         }
      }
   }
#  endif /* NDEBUG */

   radGBasisPtr = radGBasisPtrIn;
   muSet = in;
   splineRep = splineRepIn;

   createFuncLMap ();
}

SxAtomicOrbitalsG::SxAtomicOrbitalsG (const SxAtomicOrbitalsG &in)
{
   (*this) = in;
}

SxAtomicOrbitalsG::~SxAtomicOrbitalsG ()
{
   // empty
}

void SxAtomicOrbitalsG::set(double val)
{

   for (int iSpecies = 0; iSpecies < muSet.getSize(); iSpecies++)   {
      for (int iot = 0; iot < muSet(iSpecies).getSize(); iot++)   {
         muSet(iSpecies)(iot).set(val);
      }
   }
}

void SxAtomicOrbitalsG::operator= (const SxAtomicOrbitalsG &in)
{
   radGBasisPtr = in.radGBasisPtr;
   muSet = in.muSet;
   splineRep = in.splineRep;
   funcLMap = in.funcLMap;
}

void SxAtomicOrbitalsG::operator+= (const SxAtomicOrbitalsG &in)
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(muSet.getSize() == in.muSet.getSize(),
            muSet.getSize(),
            in.muSet.getSize());
   for (int is = 0; is < in.muSet.getSize (); is++)   {
      SX_CHECK(muSet(is).getSize() == in.muSet(is).getSize(),
               muSet(is).getSize(),
               in.muSet(is).getSize());
      for (int iot = 0; iot < in.muSet(is).getSize (); iot++)   {
         SX_CHECK(muSet(is)(iot).getSize() == in.muSet(is)(iot).getSize(),
                  muSet(is)(iot).getSize(),
                  in.muSet(is)(iot).getSize());
         muSet(is)(iot) += in.muSet(is)(iot);
      }
   }
}

void SxAtomicOrbitalsG::operator-= (const SxAtomicOrbitalsG &in)
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(muSet.getSize() == in.muSet.getSize(),
            muSet.getSize(),
            in.muSet.getSize());
   for (int is = 0; is < in.muSet.getSize (); is++)   {
      SX_CHECK(muSet(is).getSize() == in.muSet(is).getSize(),
               muSet(is).getSize(),
               in.muSet(is).getSize());
      for (int iot = 0; iot < in.muSet(is).getSize (); iot++)   {
         SX_CHECK(muSet(is)(iot).getSize() == in.muSet(is)(iot).getSize(),
                  muSet(is)(iot).getSize(),
                  in.muSet(is)(iot).getSize());
         muSet(is)(iot) -= in.muSet(is)(iot);
      }
   }
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator+ (const SxAtomicOrbitalsG &in) const
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(radGBasisPtr->getNElements() == in.radGBasisPtr->getNElements());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes,
               in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) + in.muSet(is)(iot);
         resMuSet(is)(iot).handle->auxData.is = in.muSet(is)(iot).handle->auxData.is;
         resMuSet(is)(iot).handle->auxData.n  = in.muSet(is)(iot).handle->auxData.n;
         resMuSet(is)(iot).handle->auxData.l  = in.muSet(is)(iot).handle->auxData.l;
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, splineRep);

   return result;
}

void SxAtomicOrbitalsG::sumMPI() const
{
   ssize_t nSpecies = muSet.getSize();
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      ssize_t nOrbTypes = getNOrbTypes (is);
      resMuSet(is).resize(nOrbTypes);
      for (ssize_t iot = 0; iot < nOrbTypes; iot++) 
         SxLoopMPI::sum (muSet(is)(iot));
   }
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator- (const SxAtomicOrbitalsG &in) const
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(radGBasisPtr->getNElements() == in.radGBasisPtr->getNElements());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes,
               in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) - in.muSet(is)(iot);
         resMuSet(is)(iot).handle->auxData.is = in.muSet(is)(iot).handle->auxData.is;
         resMuSet(is)(iot).handle->auxData.n  = in.muSet(is)(iot).handle->auxData.n;
         resMuSet(is)(iot).handle->auxData.l  = in.muSet(is)(iot).handle->auxData.l;
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator* (double skalar) const
{
   SX_CHECK(radGBasisPtr->getNElements() > 0);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) * skalar;
         resMuSet(is)(iot).handle->auxData.is = muSet(is)(iot).handle->auxData.is;
         resMuSet(is)(iot).handle->auxData.n  = muSet(is)(iot).handle->auxData.n;
         resMuSet(is)(iot).handle->auxData.l  = muSet(is)(iot).handle->auxData.l;
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator/ (double skalar) const
{
   SX_CHECK(radGBasisPtr->getNElements() > 0);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) / skalar;
         resMuSet(is)(iot).handle->auxData.is = muSet(is)(iot).handle->auxData.is;
         resMuSet(is)(iot).handle->auxData.n  = muSet(is)(iot).handle->auxData.n;
         resMuSet(is)(iot).handle->auxData.l  = muSet(is)(iot).handle->auxData.l;
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator* (const SxArray<SxArray<double> > &skalar) const
{
   SX_CHECK(radGBasisPtr->getNElements() > 0);
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == skalar.getSize(), nSpecies, skalar.getSize());
   SxArray<SxArray<SxDiracVec<Double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      SX_CHECK(nOrbTypes == skalar(is).getSize(), nOrbTypes, skalar(is).getSize());
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = skalar(is)(iot) * muSet(is)(iot);
         resMuSet(is)(iot).handle->auxData.is = muSet(is)(iot).handle->auxData.is;
         resMuSet(is)(iot).handle->auxData.n  = muSet(is)(iot).handle->auxData.n;
         resMuSet(is)(iot).handle->auxData.l  = muSet(is)(iot).handle->auxData.l;
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsG SxAtomicOrbitalsG::operator* (const SxAtomicOrbitalsG &in) const
{
  SX_CHECK(radGBasisPtr.getPtr () != NULL);
   SX_CHECK(radGBasisPtr == in.radGBasisPtr);
   SX_CHECK(in.isSpline () == isSpline());
   bool wasSpline = false;
   if (in.isSpline ()) {
      in.toVec();
      toVec();
      wasSpline = true;
   }
   int nSpecies = getNSpecies ();
   SX_CHECK (nSpecies == in.getNSpecies ());
   SxArray<SxArray<SxRadGBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      SX_CHECK (nOrbTypes == in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = in.muSet(is)(iot) * muSet(is)(iot);
      }
   }
   SxAtomicOrbitalsG result = SxAtomicOrbitalsG(resMuSet, radGBasisPtr, 0);

   if (wasSpline)  {
      in.toSpline ();
      toSpline ();
      result.toSpline ();
   }

   return result;
} 

SxDiracVec<Double> & SxAtomicOrbitalsG::operator() (int is, int iot)
{
   return muSet(is)(iot);
}

const SxDiracVec<Double> & SxAtomicOrbitalsG::operator() (int is, int iot) const
{
   return muSet(is)(iot);
}

SxDiracVec<Double> SxAtomicOrbitalsG::operator() (int is, int ia, int iot, int l, int m) const
{
   SxDiracVec<Double> result = muSet(is)(iot).getCopy ();
   // Set auxData
   result.setBasis(*radGBasisPtr);
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = iot;
   result.handle->auxData.l  = l;
   result.handle->auxData.m = m;

   return result;
}

int SxAtomicOrbitalsG::getNSpecies () const
{
   return (int)muSet.getSize ();
}

int SxAtomicOrbitalsG::getNOrbTypes () const
{
   int result = 0;
   for(int is = 0; is < getNSpecies (); is++)   
      result += getNOrbTypes(is);
      
   return result;
}

int SxAtomicOrbitalsG::getNOrbTypes (int iSpecies) const
{
   return (int)muSet(iSpecies).getSize ();
}

void SxAtomicOrbitalsG::print (const SxString fileIn) const
{
   int nSpecies = (int)muSet.getSize ();
   for(int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize();
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         SxString file = fileIn + SxString(is) + SxString(iot) + ".dat";
         SxDiracVec<Double> vec;
         if (splineRep) vec = toVec(is,iot);
         else           vec = muSet(is)(iot);
         SxBinIO out;
         out.open(file, SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(toVector(radGBasisPtr->getRadGFunc()),toVector(vec));
         out.close();
      }
   }
}

double SxAtomicOrbitalsG::getNormSqr(int is, int iot) const
{
   SxDiracVec<Double> vec; 
   if (splineRep) vec = toVec(is,iot);
   else           vec = muSet(is)(iot);
   double result = radGBasisPtr->integrate(vec * vec);
   return result;
}

double SxAtomicOrbitalsG::getNormSqrSum() const
{
   double result = 0.0;
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += getNormSqr(is,iot);
      }
   }

   return result;
}
double SxAtomicOrbitalsG::dot (const SxAtomicOrbitalsG &orbitalsIn) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == orbitalsIn.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == orbitalsIn.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += dot(orbitalsIn,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitalsG::dot (const SxAtomicOrbitalsG &orbitalsIn, int is, int iot, int js, int jot) const
{
   SX_CHECK(radGBasisPtr == orbitalsIn.radGBasisPtr);
   SxDiracVec<Double> vec1; 
   if (splineRep) vec1 = toVec(is,iot);
   else           vec1 = muSet(is)(iot);
   SxDiracVec<Double> vec2; 
   if (orbitalsIn.splineRep) vec2 = orbitalsIn.toVec(js,jot);
   else           vec2 = orbitalsIn.muSet(js)(jot);
   double result = radGBasisPtr->integrate(vec1 * vec2);
   return result;
}

double SxAtomicOrbitalsG::sum (const SxAtomicOrbitalsG &orbitalsIn) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == orbitalsIn.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == orbitalsIn.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += sum(orbitalsIn,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitalsG::sum (const SxAtomicOrbitalsG &orbitalsIn, int is, int iot, int js, int jot) const
{
   SX_CHECK(radGBasisPtr == orbitalsIn.radGBasisPtr);
   SxDiracVec<Double> vec1; 
   if (splineRep) vec1 = toVec(is,iot);
   else           vec1 = muSet(is)(iot);
   SxDiracVec<Double> vec2; 
   if (orbitalsIn.splineRep) vec2 = orbitalsIn.toVec(js,jot);
   else           vec2 = orbitalsIn.muSet(js)(jot);
   double result = (vec1 * vec2).sum ();
   return result;
}

SxArray<SxArray<SxMatrix<Double> > > SxAtomicOrbitalsG::getOverlap () const
{
   int nSpecies = getNSpecies ();
   SxArray<SxArray<SxMatrix<Double> > > result(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      result(is).resize(lMax + 1);
      for (int l = 0; l <= lMax; l++)  {
         result(is)(l) = getOverlap (is, l);
      }
   }
   return result;
}

SxMatrix<Double>  SxAtomicOrbitalsG::getOverlap (int is, int l) const
{
   int nFL = getFuncPerL (is,l);
   SxMatrix<Double> result(nFL,nFL);
   for (int ifl = 0; ifl < nFL; ifl++)  {
      int iot = funcLMap(is)(l)(ifl);
      SxDiracVec<Double> vec1; 
      if (splineRep) vec1 = toVec(is,iot);
      else           vec1 = muSet(is)(iot);
      for (int jfl = ifl; jfl < nFL; jfl++)  {
         int jot = funcLMap(is)(l)(jfl);
         SxDiracVec<Double> vec2; 
         if (splineRep) vec2 = toVec(is,jot);
         else           vec2 = muSet(is)(jot);
         result(ifl,jfl) = result(jfl,ifl) = 
            radGBasisPtr->integrate(vec1 * vec2);
      }
   }
   return result;
}

void SxAtomicOrbitalsG::setBasis (SxConstPtr<SxRadGBasis> radGBasisPtrIn)
{
   radGBasisPtr = radGBasisPtrIn;
   for (int is = 0; is < muSet.getSize (); ++is)  {
      for (int iot = 0; iot < muSet(is).getSize (); ++iot)  {
         muSet(is)(iot).setBasis (radGBasisPtr.getConstPtr ());
      }
   }
}

void SxAtomicOrbitalsG::toSpline () const
{
   SX_CHECK (!splineRep);
   if (!splineRep)  {
      SxDiracVec<Double> radGFunc = radGBasisPtr->getRadGFunc();
      int nSpecies = getNSpecies ();
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbtypes = getNOrbTypes (is);
         for (int iot = 0; iot < nOrbtypes; iot++)   {
            SX_CHECK(is == muSet(is)(iot).handle->auxData.is,
                     is, muSet(is)(iot).handle->auxData.is);
            int ia = muSet(is)(iot).handle->auxData.ia;
            int n = muSet(is)(iot).handle->auxData.n;
            int l = muSet(is)(iot).handle->auxData.l;
            int m = muSet(is)(iot).handle->auxData.m;
            SxCubicSpline<SxDiracVec<Double> > spline;
            if (l < 2) 
               spline = SxCubicSpline<SxDiracVec<Double> > ( 
                     radGFunc, 
                     muSet(is)(iot), 
                     SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
            else
               spline = SxCubicSpline<SxDiracVec<Double> > ( 
                     radGFunc, 
                     muSet(is)(iot), 
                     SxCubicSpline<SxDiracVec<Double> >::Hermite);
            muSet(is)(iot) = spline.getSpline ();
            muSet(is)(iot).setBasis(&*radGBasisPtr);
            muSet(is)(iot).handle->auxData.is = is;
            muSet(is)(iot).handle->auxData.ia = ia;
            muSet(is)(iot).handle->auxData.n = n;
            muSet(is)(iot).handle->auxData.l = l;
            muSet(is)(iot).handle->auxData.m = m;
         }
      }
      splineRep = true;
   }
}

void SxAtomicOrbitalsG::toVec () const
{
   SX_CHECK (splineRep);
   if(splineRep)  {
      SxDiracVec<Double> radGFunc = radGBasisPtr->getRadGFunc();
      int nSpecies = getNSpecies ();
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbtypes = getNOrbTypes (is);
         for (int iot = 0; iot < nOrbtypes; iot++)   {
            SX_CHECK(is == muSet(is)(iot).handle->auxData.is,
                     is, muSet(is)(iot).handle->auxData.is);
            int ia = muSet(is)(iot).handle->auxData.ia;
            int n = muSet(is)(iot).handle->auxData.n;
            int l = muSet(is)(iot).handle->auxData.l;
            int m = muSet(is)(iot).handle->auxData.m;
            SxCubicSpline<SxDiracVec<Double> > spline (radGFunc, muSet(is)(iot));
            muSet(is)(iot) = spline.getY (radGFunc);
            muSet(is)(iot).setBasis(&*radGBasisPtr);
            muSet(is)(iot).handle->auxData.is = is;
            muSet(is)(iot).handle->auxData.ia = ia;
            muSet(is)(iot).handle->auxData.n = n;
            muSet(is)(iot).handle->auxData.l = l;
            muSet(is)(iot).handle->auxData.m = m;
         }
      }
      splineRep = false;
   } 
}

SxDiracVec<Double> SxAtomicOrbitalsG::toVec (int is, int iot) const
{
   SX_CHECK (splineRep);
   SxDiracVec<Double> radGFunc = radGBasisPtr->getRadGFunc();
   SxCubicSpline<SxDiracVec<Double> > spline (radGFunc, muSet(is)(iot));
   SxDiracVec<Double> result = spline.getY (radGFunc);
   return result;
}

void SxAtomicOrbitalsG::normalize ()
{
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         double norm = getNormSqr(is,iot);
         if (norm > 1e-12) muSet(is)(iot) /= sqrt(norm);
         else {
            cout << "WARNING: Try to normalize zero norm orbital." << endl;
            cout << "Norm is " << norm << endl;
            cout << "Orbital is set to zero." << endl;
            muSet(is)(iot).set(0.0);
         }
      }
   }
}

void SxAtomicOrbitalsG::createFuncLMap ()
{
   int nSpecies = (int)muSet.getSize ();
   funcLMap.resize(nSpecies);
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int lMax = getLMax (iSpecies);
      funcLMap(iSpecies).resize(lMax + 1);
      SxArray<int> funcPerL(lMax + 1);
      funcPerL.set(0);
      int nOrbTypes = (int)muSet(iSpecies).getSize ();
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         int l = muSet(iSpecies)(iot).handle->auxData.l;
         funcPerL(l)++;
      }
      SxArray<int> currentIFL(lMax + 1);
      currentIFL.set(0);
      for (int l = 0; l <= lMax; l++)  {
         funcLMap(iSpecies)(l).resize(funcPerL(l));
      }
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         int l = muSet(iSpecies)(iot).handle->auxData.l;
         funcLMap(iSpecies)(l)(currentIFL(l)) = iot;
         currentIFL(l)++;
      }
   }
}

SxRadGBasis::TPsi &SxAtomicOrbitalsG::getFuncL (int is, int l, int ifl)
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

const SxRadGBasis::TPsi &SxAtomicOrbitalsG::getFuncL (int is, int l, int ifl) const
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

SxArray<SxArray<SxMatrix<Double> > > SxAtomicOrbitalsG::orthogonalize ()
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitalsG org = 1.0 * *this;
   int nSpecies = (int)muSet.getSize ();
   SxArray<SxArray<SxMatrix<Double> > > result(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      result(is).resize(lMax + 1);
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL (is,l);
         if (nFL == 0) continue;
         SxMatrix<Double> S = getOverlap(is,l);
         result(is)(l) = S.choleskyDecomposition ().adjoint (). inverse ();
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            muSet(is)(iot).set(0.0);
            for (int jfl = 0; jfl < nFL; jfl++)  {
               muSet(is)(iot) += result(is)(l)(jfl,ifl) * org.getFuncL(is,l,jfl);
            }
         }
      }
   }

   // return CholeskyDecompositions
   return result;
}

int SxAtomicOrbitalsG::getLMax () const
{
   int nSpecies = (int)muSet.getSize ();
   int result = 0;
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(is)(iot).handle->auxData.l)
            result = muSet(is)(iot).handle->auxData.l;
      }
   }
   return result;
}

int SxAtomicOrbitalsG::getLMax (const int iSpecies) const
{
   int nOrbTypes = (int)muSet(iSpecies).getSize ();
   int result = 0;
   for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(iSpecies)(iot).handle->auxData.l)
            result = muSet(iSpecies)(iot).handle->auxData.l;
   }
   return result;
}

int SxAtomicOrbitalsG::getFuncPerL (int is, int l) const
{
   return (int)funcLMap(is)(l).getSize();
}

void SxAtomicOrbitalsG::orthogonalizeOn(SxAtomicOrbitalsG &basis)
{
   //Gram-Schmidt ortho
   bool switchThis = false;
   bool switchBasis = false;
   int nSpecies = (int)muSet.getSize ();
   SX_CHECK(nSpecies == basis.muSet.getSize());
   
   if (isSpline()) { 
      toVec ();
      switchThis = true;   
   }

   if (basis.isSpline())  {
      basis.toVec ();
      switchBasis = true;
   }

   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      SX_CHECK(lMax <= basis.getLMax(is));
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL(is,l);
         int nFLbasis = basis.getFuncPerL(is,l);
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            SxDiracVec<Double> vec = muSet(is)(iot).getCopy();
            for (int jfl = 0; jfl < nFLbasis; jfl++)  {
               int jot = basis.funcLMap(is)(l)(jfl);
               muSet(is)(iot) -= radGBasisPtr->integrate(vec * basis.muSet(is)(jot))
                               / basis.getNormSqr(is,jot) * basis.muSet(is)(jot);
            }
         }
      }
   }

   if (switchThis) toSpline ();
   if (switchBasis) basis.toSpline ();
}

void SxAtomicOrbitalsG::rotate(SxArray<SxArray<SxMatrix<Double> > > rotMat)
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitalsG org = 1.0 * *this;
   int nSpecies = (int)muSet.getSize ();
   SX_CHECK (nSpecies == rotMat.getSize(),nSpecies, rotMat.getSize());
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      SX_CHECK (lMax + 1 == rotMat(is).getSize(), lMax + 1, rotMat(is).getSize());
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL (is,l);
         SX_CHECK (nFL == rotMat(is)(l).row(0).getSize(), nFL, rotMat(is)(l).row(0).getSize());
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            muSet(is)(iot).set(0.0);
            for (int jfl = 0; jfl < nFL; jfl++)  {
               muSet(is)(iot) += rotMat(is)(l)(jfl,ifl) * org.getFuncL(is,l,jfl);
            }
         }
      }
   }
}

int SxAtomicOrbitalsG::getNOrbitals (
      SxAtomicStructure structure) const
{
   SX_CHECK(muSet.getSize() > 0);
   int nSpecies = (int)muSet.getSize ();
   int result = 0;
   for(int is = 0; is < nSpecies; is++)   {
      int nAtoms = structure.getNAtoms(is);
      int nOrbLocal = (int)muSet(is).getSize();
      for(int ia = 0; ia < nAtoms; ia++)   {
         for(int iot = 0; iot < nOrbLocal; iot++)   {
            int l = muSet(is)(iot).handle->auxData.l;
            result += 2*l+1;
         }
      }
   }

   return result;
}

SxArray<SxQuantumNumbers> SxAtomicOrbitalsG::getReducedOrbitalMap () const
{
   SX_CHECK(&muSet);
   int nSpecies = (int)muSet.getSize ();
   SxList<SxQuantumNumbers> list;
   for(int is = 0; is < nSpecies; is++)   {
      int nOrbLocal = (int)muSet(is).getSize();
      for(int iot = 0; iot < nOrbLocal; iot++)   {
         int l = muSet(is)(iot).handle->auxData.l;
         // SxQuantumNumberConstructor forbids uninitialized ia and m
         SxQuantumNumbers qNumbers (is,0,iot,l,0);
         // Therefore unset it here
         qNumbers.iAtom = -1; qNumbers.m = NONE_M; 
         list.append(qNumbers);
      }
   }
   SxArray<SxQuantumNumbers> result = SxArray<SxQuantumNumbers>(list);

   return result;
}


SxArray<SxQuantumNumbers> SxAtomicOrbitalsG::getOrbitalMap (
      SxAtomicStructure structure) const
{
   SX_CHECK(&muSet);
   int nSpecies = (int)muSet.getSize ();
   SxList<SxQuantumNumbers> list;
   for(int is = 0; is < nSpecies; is++)   {
      int nAtoms = structure.getNAtoms(is);
      int nOrbLocal = (int)muSet(is).getSize();
      for(int ia = 0; ia < nAtoms; ia++)   {
         for(int iot = 0; iot < nOrbLocal; iot++)   {
            int l = muSet(is)(iot).handle->auxData.l;
            for(int m = -l; m <= l; m++)   {
               list.append(SxQuantumNumbers(is,ia,iot,l,m));
            }
         }
      }
   }

   SxArray<SxQuantumNumbers> result = SxArray<SxQuantumNumbers>(list);

   return result;
}

int SxAtomicOrbitalsG::getIOT (int is, int n, int l) const
{
   int iot = 0;
   while ((muSet(is)(iot).handle->auxData.n != n) 
         || (muSet(is)(iot).handle->auxData.l != l))   {
      iot++;
      if (iot >= muSet(is).getSize()) SX_EXIT;
   }

   return iot;
}

int SxAtomicOrbitalsG::getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const
{
   for (int iOrbital = 0; iOrbital < map.getSize(); iOrbital++)  {
      if (map(iOrbital).iSpecies == is  &&
          map(iOrbital).iAtom    == ia  &&  
          map(iOrbital).n        == iot &&
          map(iOrbital).l        == l   &&
          map(iOrbital).m        == m) return iOrbital;
   }

   cout << "OrbitalIdx with is = " << is << ", ia = " << ia
        << ", n = " << iot << ", l = " << l << ", m = " << m
        << " not found!" << endl;
   SX_EXIT;
}
