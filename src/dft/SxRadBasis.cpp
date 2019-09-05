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

#include <SxRadBasis.h>
#include <SxGBasis.h>
#include <SxRadGBasis.h>
#include <SxRadRBasis.h>
#include <SxConstants.h>
#include <SxMathLib.h>
#include <stdio.h>
#include <math.h>
#include <SxYlm.h>
#include <SxPseudoPot.h>
#include <SxPAWPot.h>


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------


SxRadBasis::SxRadBasis ()
{
   registerMemoryObservers ();
}


SxRadBasis::SxRadBasis (const SxArray<SxDiracVec<TReal8> > &radFunc_,
                        const SxArray<Real8>               &logDr_)
{
   set(radFunc_, logDr_);
}

void SxRadBasis::set (const SxArray<SxDiracVec<TReal8> > &radFunc_,
                      const SxArray<Real8>               &logDr_)
{
   this->radFunc = radFunc_;
   this->logDr   = logDr_;

   registerMemoryObservers ();
}

SxRadBasis::SxRadBasis (double rMin, double rMax, int nPoints)
{
   set(rMin, rMax, nPoints);
}

SxRadBasis::SxRadBasis (const SxVector<Double> &rMin, 
                        const SxVector<Double> &rMax, 
                        const SxVector<Int> &nPoints)
{
   set(rMin, rMax, nPoints);
}

SxRadBasis::SxRadBasis (const SxSymbolTable *table)
{
   setup(table);
}

void SxRadBasis::set (double rMin, double rMax, int nPoints)
{
   SX_CHECK (nPoints > 1, nPoints);
   SX_CHECK (rMin > 0., rMin);
   SX_CHECK (rMax > rMin, rMax, rMin);

   radFunc.resize (1); // one species
   logDr.resize (1);
   SxDiracVec<TReal8> r(nPoints);
   double ldr = log(rMax/rMin) / double(nPoints - 1);
   for (int i = 0; i < nPoints; i++)
      r(i) = rMin * exp(ldr * i);
   SX_CHECK (fabs(r(nPoints-1)-rMax) < 1e-12 * rMax, r(nPoints-1), rMax);
   radFunc(0) = r;
   this->logDr(0) = ldr;
   registerMemoryObservers ();
}

void SxRadBasis::set (const SxVector<Double> &rMin, 
                      const SxVector<Double> &rMax, 
                      const SxVector<Int> &nPoints)
{
   SX_CHECK(rMin.getSize () == rMax.getSize ());
   SX_CHECK(rMin.getSize () == nPoints.getSize ());

   ssize_t nSpecies = rMin.getSize ();
   for (ssize_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      SX_CHECK (nPoints(iSpecies) > 1, nPoints(iSpecies));
      SX_CHECK (rMin(iSpecies) > 0., rMin(iSpecies));
      SX_CHECK (rMax(iSpecies) > rMin(iSpecies),
                rMax(iSpecies), rMin(iSpecies));

      radFunc.resize (nSpecies); 
      logDr.resize (nSpecies);
      int dim = nPoints(iSpecies);
      SxDiracVec<TReal8> r(dim);
      double ldr = log(rMax(iSpecies)/rMin(iSpecies)) / double(dim - 1);
      for (int i = 0; i < dim; i++)
         r(i) = rMin(iSpecies) * exp(ldr * i);
      SX_CHECK (fabs(r(dim-1)-rMax(iSpecies)) < 1e-12 * rMax(iSpecies),
                r(dim-1), rMax(iSpecies));
      radFunc(iSpecies) = r;
      this->logDr(iSpecies) = ldr;
   }
   registerMemoryObservers ();
}

SxRadBasis::SxRadBasis (const SxString &file)
{
   this->read(file);
}

SxRadBasis::SxRadBasis (const SxBinIO &io)
{
   this->read(io);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxRadBasis::~SxRadBasis ()
{
   deregisterAll ();
}

void SxRadBasis::read(const SxString &file)
{
try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}
void SxRadBasis::setup (const SxSymbolTable *table)
{
   SX_CHECK (table);
   const SxString specName = "species";
   const SxString basisName = "radBasis";
   const SxSymbolTable *species, *basis = NULL;
   int iSpecies=0;
   SxRadBasis result;

   try {
      species = table->getGroup (specName);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   try   {
      
      // get species
      for(species = table->getGroup (specName);
          species != NULL;
          species = species->nextSibling (specName), iSpecies++)
      {
         if (species->containsGroup (basisName))
            basis = species->getGroup (basisName);
         else  {
            cout << "No radial basis specified." << endl;
            cout << "SPHInX quits here!" << endl;
            SX_QUIT;
         }
         if(basis->containsGroup("fromFile"))   {
            const SxSymbolTable *mode = basis->getGroup ("fromFile");
            try {
               SxString fileName = mode->get("file")->toString();
               int is = mode->get("is")->toInt();
               SxRadBasis fileBasis;
               fileBasis.read(fileName);
               result.addMesh(fileBasis.radFunc(is),fileBasis.logDr(is));
            } catch (SxException e)   {
               e.print ();
               SX_EXIT;
            }
         } else if (basis->containsGroup("generate")) {
            const SxSymbolTable *mode = basis->getGroup ("generate");
            try {
               double rMin = mode->get("rMin")->toReal();
               double rMax = mode->get("rMax")->toReal();
               int nPoints = mode->get("nPoints")->toInt();
               SxRadBasis genBasis;
               genBasis.set(rMin, rMax, nPoints);
               result.addMesh(genBasis.radFunc(0),genBasis.logDr(0));
            } catch (SxException e)   {
               e.print ();
               SX_EXIT;
            }
         } else if (basis->containsGroup("fromPotential"))  {   
            const SxSymbolTable *mode = basis->getGroup ("fromPotential");
            SxString potFile = mode->get("file")->toString();
            SxParser potParser;
            SxConstPtr<SxSymbolTable> potTable = potParser.read (potFile);
            if (potTable->containsGroup("pseudoPot"))   {
               SxPseudoPot pot = SxPseudoPot(&*potTable);
               SxPtr<SxRadBasis> radBasisPotPtr 
                  = SxPtr<SxRadBasis>::create(pot.rad, pot.logDr);
               try {
                  int is = mode->get("is")->toInt();
                  result.addMesh(radBasisPotPtr->radFunc(is),
                        radBasisPotPtr->logDr(is));
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else if (potTable->containsGroup("pawPot"))  {
               SxPAWPot pot = SxPAWPot(&*potTable);
               SxPtr<SxRadBasis> radBasisPotPtr 
                  = SxPtr<SxRadBasis>::create(pot.rad,pot.logDr);
               try {
                  int is = mode->get("is")->toInt();
                  result.addMesh(radBasisPotPtr->radFunc(is),
                        radBasisPotPtr->logDr(is));
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else  {
               cout << "No known potential group found!" << endl;
               cout << "SPHInX quits here!" << endl;
               SX_QUIT;
            }
         } else {
            cout << "No known radial basis initialization specified!" << endl;
            cout << "SPHInX quits here!" << endl;
            SX_QUIT;
         }
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }

   this->set(result.radFunc,result.logDr);
}

void SxRadBasis::read (const SxBinIO &io)
{
   SxArray<SxDiracVec<Double> > radFunc_;
   SxArray<double> logDr_;

   try {
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      radFunc_.resize(nSpecies);
      logDr_.resize(nSpecies);
      for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
         SxString dimRadName = "dimRad-" + SxString(iSpecies);
         int dimRad = io.getDimension (dimRadName);
         radFunc_(iSpecies).resize(dimRad);
         SxDiracVec<Double> &Vec = radFunc_(iSpecies);
         SxString radialName = "radFunc-" + SxString(iSpecies); 
         io.read(radialName,&Vec,dimRad);
         SxString logDrName = "logDr-"+SxString(iSpecies);
         double value;
         io.read(logDrName, &value);
         logDr_(iSpecies) = value;
      }

   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }


   set(radFunc_, logDr_);
}

SxRadBasis SxRadBasis::readMesh(const SxString &file)
{
   SxRadBasis result;
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      result = readMesh (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return result;
}

SxRadBasis SxRadBasis::readMesh (const SxBinIO &io)
{
   SxArray<SxDiracVec<Double> > radFunc_;
   SxArray<double> logDr_;

   try {
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      radFunc_.resize(nSpecies);
      logDr_.resize(nSpecies);
      for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
         SxString dimRadName = "dimRad-" + SxString(iSpecies);
         int dimRad = io.getDimension (dimRadName);
         radFunc_(iSpecies).resize(dimRad);
         SxDiracVec<Double> &Vec = radFunc_(iSpecies);
         SxString radialName = "radFunc-" + SxString(iSpecies); 
         io.read(radialName,&Vec,dimRad);
         SxString logDrName = "logDr-"+SxString(iSpecies);
         double value;
         io.read(logDrName, &value);
         logDr_(iSpecies) = value;
      }

   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }


   SxRadBasis result(radFunc_, logDr_);
   return result;
}

int SxRadBasis::addMesh (const SxDiracVec<TReal8> &radFuncIn,
                               double             logDrIn)
{
   int is = int(radFunc.getSize ());
   int nSpecies = is + 1;
   radFunc.resize (nSpecies, true);
   logDr.resize (nSpecies, true);

   radFunc(is) = radFuncIn;
   logDr(is) = logDrIn;

   cashedRl.resize (0);
   gBases.resize (0);

   return is;
}

double SxRadBasis::getRMax () const
{
   SxVector<Double> rMax (getNSpecies ());

   for (int is = 0; is < getNSpecies (); is++) rMax(is) = getRMax(is);

   return rMax.maxval ();
}

SxDiracVec<SxRadBasis::TBasisType>
SxRadBasis::identity (const SxRadBasis *basis,
                       const SxDiracVec<TBasisType> &vec) const
{
   // --- verify that this is really the identity projector, i.e.
   //     that the provided basis is the same as the vector's basis
   SX_CHECK (vec.getBasisPtr() == this);
   SX_CHECK (basis == this);
   return vec;
}

SxDiracVec<SxRadBasis::TBasisType>
SxRadBasis::changeRadBasis (const SxRadBasis *basisPtr,
                            const SxDiracVec<TBasisType> &vec) const
{
   SX_CHECK(vec.getBasisPtr() == this);
   
   const SxRadBasis &basis = *basisPtr;
   
   int is = vec.handle->auxData.is;
   int l = vec.handle->auxData.l;

   // Check for identical basis
   if (basis.radFunc(is).getSize() == radFunc(is).getSize())
      if ((radFunc(is) - basis.radFunc(is)).norm() < 1e-10) return (1.0 * vec);

   SX_CHECK(is < basis.radFunc.getSize(), is, basis.radFunc.getSize());
   // Fourierinterpolation
   // result = (newRad|radG|vec);
   /*
   SxRadGBasis radGBasis(0.0, 30.0, 3000, SxRadGBasis::Linear);
   SxDiracVec<TBasisType> result = radGBasis.toRadBasis(&basis,toRadGBasis(&radGBasis,vec));
   // SplineInterpolation
   */
   ssize_t dim = radFunc(is).getSize ();
   ssize_t dimNew = basis.radFunc(is).getSize();
   SxDiracVec<TBasisType> x,y;
   SxCubicSpline<SxDiracVec<TBasisType> > spline;
   // Check for boundary
   if ((basis.radFunc(is)(dimNew-1) - radFunc(is)(dim-1)) > 1e-4)  {
      x = 1.0 * radFunc(is);
      y = 1.0 * vec;

      double xLast = radFunc(is)(dim-1);
      double yLast = vec(dim-1);
      double xPreLast = radFunc(is)(dim-2);
      double yPreLast = vec(dim-2);
      double slope = (yLast - yPreLast)/(xLast - xPreLast);

      double xLastNew = basis.radFunc(is)(dimNew-1);
      double slopeRight = 0.0;

     if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  { 
        // force to zero
        x.resize(dim+1,true);
        y.resize(dim+1,true);
        if (xLastNew - xLast > 5.0)
           x(dim) = xLastNew;
        else
           x(dim) = xLastNew + 5.0;
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
        spline 
           = SxCubicSpline<SxDiracVec<TBasisType> > (
         x,y, 
         SxCubicSpline<SxDiracVec<TBasisType> >::NaturalHermite,
         slopeRight);
     } else  {
        spline 
           = SxCubicSpline<SxDiracVec<TBasisType> > (
         x,y, 
         SxCubicSpline<SxDiracVec<TBasisType> >::Hermite,
         0.0,
         slopeRight);
     }
   } else  {
      x = radFunc(is);
      y = vec;
      if (l < 2) {
         spline = SxCubicSpline<SxDiracVec<TBasisType> >(
               x,y, 
               SxCubicSpline<SxDiracVec<TBasisType> >::NaturalHermite);
      } else  {
         spline  = SxCubicSpline<SxDiracVec<TBasisType> > (
               x,y, 
               SxCubicSpline<SxDiracVec<TBasisType> >::Hermite);
      }
   }
   
   SxDiracVec<TBasisType> result = spline.getY(basis.radFunc(is));

   cout << SX_SEPARATOR;
   cout << "WARNING: Change of RADBASIS by interpolation!" << endl;
   SxString file = "ORIG-" + SxString(is) + SxString(l) + ".dat";
   cout << "Check " << file << endl;
   SxBinIO out; 
   out.open(file, SxBinIO::ASCII_WRITE_ONLY);
   out.writeXYPlot(toVector(radFunc(is)),
                   toVector(vec));
   out.close();

   file = "INTERPOL-" + SxString(is) + SxString(l) + ".dat";
   cout << "Check " << file << endl;
   out.open(file, SxBinIO::ASCII_WRITE_ONLY);
   out.writeXYPlot(toVector(basis.radFunc(is)),
                   toVector(result));
   out.close();
   cout << SX_SEPARATOR;
   
   result.setBasis(basisPtr);
   result.handle->auxData.is = vec.handle->auxData.is;
   result.handle->auxData.ia = vec.handle->auxData.ia;
   result.handle->auxData.n  = vec.handle->auxData.n;
   result.handle->auxData.l  = vec.handle->auxData.l;
   result.handle->auxData.m  = vec.handle->auxData.m;

   return result;
}


SxDiracVec<TGBasisType>
SxRadBasis::toPWBasis (const SxGBasis *basis,
                       const SxDiracVec<TBasisType> &vec) const
{
   SX_CLOCK (Timer::RadTotal);
   SX_CHECK (basis->structPtr);
   SX_CHECK (basis->structPtr->cell.volume > 0.);
   double rNorm = FOUR_PI / sqrt (basis->structPtr->cell.volume);
   SxDiracVec<TGBasisType> res;

   int is = vec.handle->auxData.is;  // TODO: ugly
   int ia = vec.handle->auxData.ia;  // TODO: ugly
   int l  = vec.handle->auxData.l;   // TODO: ugly
   int m  = vec.handle->auxData.m;   // TODO: ugly

   bool useCache =  cashedVec.getSize () == vec.getSize () 
                 && cashedVec.handle->auxData.l  != l
                 && cashedVec.handle->auxData.is != is;
   if (useCache)  {
      // --- check cached vector elements
      for (int i = 0; i < vec.getSize (); ++i)  {
         if (vec(i) != cashedVec(i))  {
            useCache = false;
            break;
         }
      }
   }

   if (!useCache) {
      // --- clear cache
      for (int jk = 0; jk < cashedRl.getSize (); ++jk)
         cashedRl(jk) = SxDiracVec<Double> ();
      cashedVec = vec.getCopy ();
      cashedVec.handle->auxData.l = l;
      cashedVec.handle->auxData.is = is;
   }

   int ik = int(gBases.findPos (basis));
   if (ik == -1)  {
      // new G basis
      gBases.append (basis);
      int nk = int(gBases.getSize ());
      ik = nk - 1;
      cashedRl.resize (nk, true);

      UPDATE_MEMORY (cashedRl);
      UPDATE_MEMORY (cashedVec);
   }

   if (cashedRl(ik).getSize () > 0) {
      // use cashed result
      // cout << "Use cash: ik =" << ik << "; l=" << l << "; m=" << m << endl;
   } else {
      // cout << "Set cash: ik =" << ik << "; l=" << l << "; m=" << m << endl;
      cashedRl(ik) = toPWBasis (radFunc(is), vec, *basis, l, logDr(is));
   }

   // multiply with Ylm
   res = cashedRl(ik) * basis->getYlm (l,m);
   // normalization factors
   res *= SxYlm::getYlmNormFactor(l,m) * rNorm;

   // --- do we calculate <G|r><r|Psi> or <G|r><r|T*Psi>?
   if (ia >= 0)  {
      SX_CLOCK (Timer::Phase);
      res *= basis->getPhaseFactors(is,ia);
   }

   res.setBasis (basis);

   res.handle->auxData.is = is; // TODO: ugly
   res.handle->auxData.ia = ia; // TODO: ugly
   res.handle->auxData.l  = l;  // TODO: ugly
   res.handle->auxData.m  = m;  // TODO: ugly

   return res;

}
SxDiracVec<TRadGBasisType> SxRadBasis::toRadRBasis (const SxRadRBasis *radRBasisPtr,
                                                    const SxDiracVec<TBasisType> &vec) const
{
   SX_CLOCK(Timer::rad2radR);

   // Get aux data
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;

   const SxDiracVec<TRadRBasisType> &radRFunc = radRBasisPtr->getRadRFunc();

   ssize_t dim = radFunc(is).getSize ();
   ssize_t newDim = radRFunc.getSize ();
   SxDiracVec<TBasisType> x = 1.0 * radFunc(is),
                          y = 1.0 * vec;
   SxCubicSpline<TPsi> spline;

   // Check for boundary
   if ((radRFunc(newDim-1) - radFunc(is)(dim-1)) > 1e-4)  {
      double xLast = radFunc(is)(dim-1);
      double yLast = vec(dim-1);
      double xPreLast = radFunc(is)(dim-2);
      double yPreLast = vec(dim-2);
      double slope = (yLast - yPreLast)/(xLast - xPreLast);

      double xLastNew = radRFunc(newDim-1);
      double slopeRight = 0.0;

     if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  { 
        // force to zero
        x.resize(dim+1,true);
        y.resize(dim+1,true);
        if (xLastNew - xLast > 5.0)
           x(dim) = xLastNew;
        else
           x(dim) = xLastNew + 5.0;
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

   TPsi result = spline.getY(radRFunc);

   // set aux data
   result.handle->auxData.is = is;
   result.handle->auxData.ia = ia;
   result.handle->auxData.n  = n;
   result.handle->auxData.l  = l;
   result.handle->auxData.m  = m;
   result.setBasis (radRBasisPtr);

   return result;

}
SxDiracVec<TRadGBasisType> SxRadBasis::toRadGBasis (
      const SxRadGBasis *radGBasisPtr,
      const SxDiracVec<TBasisType> &vec) const
{
   SX_CLOCK (Timer::rad2radG);

   // Get aux data
   int is = vec.handle->auxData.is;
   int ia = vec.handle->auxData.ia;
   int n  = vec.handle->auxData.n;
   int l  = vec.handle->auxData.l;
   int m  = vec.handle->auxData.m;
   
   const SxDiracVec<Double> &radGFunc = radGBasisPtr->getRadGFunc();
   int ng = (int)radGFunc.getSize ();
   SxDiracVec<TRadGBasisType> result(ng);
   SxDiracVec<TBasisType> vec_rad3 = vec * radFunc(is).cub();

   for (int ig = 0; ig < ng; ig++)  {
      SxDiracVec<Double> jl = jsb (l, radGFunc(ig) * radFunc(is));
      result(ig) = (vec_rad3 * jl).integrate(logDr(is));
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

inline double weightSimpson (int i, int n)
{
   SX_CHECK (i >= 0 && i < n, i, n);
   if (i == 0) return 1. / 3.;
   int d = n - i;
   if (d > 4) return ((i & 1) ? 4. : 2. ) / 3.;
   if (n & 1)  {
      if (d == 1) return 1./3.;
      if (d & 1) return 2./3.;
      return 4. / 3.;
   } // else
   if (d == 1) return 3./8.;
   if (d == 4) return 17./24.;
   return 9. / 8.;
}


//------------------------------------------------------------------------------
// This function projects the radial vector onto plane waves.
//------------------------------------------------------------------------------
SxDiracVec<TRadBasisType> SxRadBasis::toPWBasis
     (const SxDiracVec<TReal8> &rad,
      const SxDiracVec<TReal8> &psi,
      const SxGBasis &pwBasis,
      int l, Real8 logDrIn) const
{
   SX_CHECK (rad.getSize() == psi.getSize() && rad.getSize() > 0,
             rad.getSize(), psi.getSize());
   SX_CLOCK (Timer::JsbInt);

   ssize_t ig, ng = pwBasis.ng;
   int nr = (int)psi.getSize ();
   SxDiracVec<TReal8>           psi_rad3, jl;
   SxDiracVec<TRadBasisType>    g;
   Real8                        gLast = 0., gAbsLast;

   g.resize (ng);

   // ---     radFunc * radX
   psi_rad3 = psi     * rad.cub();

   gAbsLast = -1.;
#ifdef USE_OPENMP
#pragma omp parallel
#pragma omp for firstprivate(gLast,gAbsLast)
#endif
   for (ig = 0; ig < ng; ig++)  {
      Real8 gAbs  = sqrt (pwBasis.g2(ig));
      if ( ig == 0 || fabs(gAbsLast - gAbs) > 1e-10 )  {
         gLast = 0.;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:gLast)
#endif
         for (int ir = 0; ir < nr; ++ir)
            gLast += weightSimpson (ir, nr) * psi_rad3(ir)
                     * SxYlm::jsb(l, gAbs * rad(ir));
         gLast *= logDrIn;
         g(ig) = gLast;

         gAbsLast = gAbs;
      }  else
         g(ig) = gLast;  // integrals for equal |G|s are equal!
   }

   g.setBasis (&pwBasis);
   return g;
}

Real8 SxRadBasis::tr (const SxDiracVec<Double> &vec) const
{
   SX_CHECK (vec.handle);
   int is = vec.handle->auxData.is;
   SX_CHECK (is >= 0 && is < radFunc.getSize (), is, radFunc.getSize ());
   // TODO: check if this is time critical
   // explicit for loop might be faster
   return (vec * radFunc(is).cub ()).integrate (logDr(is));
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
//
// Mathematica:
//
//    jSphericalBessel[n_Integer, z_] :=
//       fn[n, z] Sin[z] + (-1)^(n + 1) fn[-n - 1, z] Cos[z];
//       fn[n_Integer, z_] :=
//          Which[
//             n == 0,  1/z,
//             n == 1,  1/z^2,
//             n > 1, (2 n - 1)/z fn[n - 1, z] - fn[n - 2, z],
//             n < 0, (2 n + 3)/z fn[n + 1, z] - fn[n + 2, z]
//          ];
//    jsb = Table[jSphericalBessel[l, z], {l, 0, 3}] // Simplify
//
SxDiracVec<TReal8> SxRadBasis::jsb (int l, const SxDiracVec<TReal8> &z)
{

   SX_CHECK (z.getSize () > 0, z.getSize());

   int zSize = (int)z.getSize ();
   SxDiracVec<TReal8> vec(zSize);
   // --- use generic jsb
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
   for (int i = 0; i < zSize; i++)
      vec(i) = SxYlm::jsb(l, z(i));
   VALIDATE_VECTOR (vec);
   return vec;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
// return LegendreP [l, x]; x = cos(theta)
Real8 SxRadBasis::pl (int l, Real8 x)
{
   SX_CHECK (l >= SxQuantumNumbers::s && l <= SxQuantumNumbers::f, l);

   switch ( l )  {
      case SxQuantumNumbers::s :
         return 1.;
         break;
      case SxQuantumNumbers::p :
         return x;
         break;
      case SxQuantumNumbers::d :
         return 0.5 * ( 3.*x*x - 1. );
         break;
      case SxQuantumNumbers::f :
         return 0.5 * ( 5.*x*x*x - 3.*x );
         break;
      default:
         SX_EXIT;
   }
   return 0.;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Real8 SxRadBasis::cosTheta (int ig, int jg, const SxGBasis &g)
{
   Real8 giTgj = sqrt ( g.g2(ig) * g.g2(jg)  );
   if (giTgj <  1.e-7) return 1.;
   else return ( ( g.getG(ig) ^ g.getG(jg) ) / giTgj ) ;
}


// see numrec eq. 4.1.14
Real8 SxRadBasis::integrate (const SxDiracVec<TReal8> &vec) const
{
   SX_EXIT;
   SX_CHECK (vec.getSize() >= 10, vec.getSize());

   double res = 0.;

   const double C0=-31./48.;
   const double C1=+11./48.;
   const double C2= -5./48.;
   const double C3= +1./48.;


   double r1  = radFunc(0)(0);
   double dex = logDr(0);

   const double xexp = exp(dex);
   int nr = (int)vec.getSize();


   // --- integrate from 0 to first grid point
   double r2=r1*xexp;
   double r3=r2*xexp;
   double s21=(vec(2)-vec(1))/(r2-r1);
   double s31=(vec(3)-vec(1))/(r3-r1);
   double g1=(s21-s31)/(r2-r3);
   double b=s21-g1*(r1+r2);
   double a=vec(1)-b*r1-g1*r1*r1;
   res = a*r1 + .5*b*r1*r1 + 1./3.*g1*r1*r1*r1;


   // --- map onto linear grid
   SxDiracVec<TReal8> f (nr);
   double ri=r1/xexp;
   int ir;
   for (ir=0; ir < nr; ir++)  {
      ri=ri*xexp;
      f(ir)=vec(ir)*dex*ri;
   }

   // --- summation
   res=res+C0*f(1)+C1*f(2)+C2*f(3)+C3*f(4);
   for (ir=0; ir < nr; ir++)  {
      res=res+f(ir);
   }
   if (nr <= 4)
      res=res+C0*f(nr-1)+C1*f(nr-2)+C2*f(nr-3)+C3*f(nr-4);
   else
      res=res-0.5*f(nr-1);

   return res;


}

SxDiracMat<Double> SxRadBasis::realYlm (int lmax, 
                                        const SxGBasis &G,
                                        const Coord &dG)
{
   SX_CHECK (lmax >= 0, lmax);
   int nl = lmax + 1;
   int nLm = nl*nl;
   int ng = G.ng;
   SX_CHECK (ng >= 0, ng);
   SxVector<Double> Ylm(nLm);

   SxDiracMat<Double> result (ng, nLm);

   // set up Ylm for all G vectors
   for (int ig = 0; ig < ng; ig++)  {
      Coord Gk = G.getG (ig) + dG;
      if (Gk.normSqr () < 1e-12)  {
         result (ig, 0) = 1.;
         for (int lm = 1; lm < nLm; lm++)
            result(ig,lm) = 0.;
      } else {
         SxYlm::getYlmArray(lmax, Gk(0), Gk(1), Gk(2), &Ylm);
         for (int lm = 0; lm < nLm; lm++)
            result(ig,lm) = Ylm(lm);
      }
   }
   result.setBasis (&G);

   return result;
}


void SxRadBasis::registerMemoryObservers ()
{
   TRACK_MEMORY (radFunc);
   TRACK_MEMORY (cashedRl);
   TRACK_MEMORY (cashedVec);
   TRACK_MEMORY (cashYlm);
   TRACK_MEMORY (YlmNormFactors);
}

// --- Stand alone interpolation function
#include <SxNaturalCubicSpline.h>
SxDiracVec<Double> interpolateRad(const SxDiracVec<Double> &psi,
                                  double r0, double logDr,
                                  const SxDiracVec<Double> &newRad)
{
   SX_CHECK (logDr > 0., logDr);
   SX_CHECK (r0 > 0., r0);
   // natural cubic spline interpolation
   SxNaturalCubicSpline cubSpline (psi, true);
   // 1st derivative at r0
   double dPsi = (psi(1) - psi(0)) / ((exp(logDr) - 1.) * r0);

   int n = (int)newRad.getSize ();
   SxDiracVec<Double> res(n);
   for (int i = 0; i < newRad.getSize (); ++i)  {
      double r = newRad(i);
      if (r <= r0)  {
         // linear extrapolation below r0
         res(i) = dPsi * (r - r0) + psi(0);
      } else {
         double x = log(r/r0) / logDr;
         SX_CHECK ((x - (double)psi.getSize ()) < 1e-10, x,psi.getSize ());
         res(i) = cubSpline.getVal (x);
      }
   }
   VALIDATE_VECTOR(res);
   res.handle->auxData = psi.handle->auxData;
   return res;
}

