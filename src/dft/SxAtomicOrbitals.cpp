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

#include <SxAtomicOrbitals.h>
#include <SxNaturalCubicSpline.h>
#include <SxCubicSpline.h>
#include <SxRegex.h>
#include <SxFileIO.h>

SxAtomicOrbitals::SxAtomicOrbitals ()
{
   // empty
}

SxAtomicOrbitals::SxAtomicOrbitals (
      const SxArray<SxArray<SxRadBasis::TPsi> > &in,
      const SxConstPtr<SxRadBasis> radBasisPtrIn)
   : SxPsiSet (PW)
{
   SX_CHECK (radBasisPtrIn.getPtr () != NULL);

   radBasisPtr = radBasisPtrIn;
   muSet     =  in;

   // --- make sure that basis is set
   const SxRadBasis *radPtr = radBasisPtr.getPtr ();
   int iSpecies, nSpecies = (int)in.getSize();
   SX_CHECK (nSpecies <= radBasisPtrIn->radFunc.getSize(),
             nSpecies,   radBasisPtrIn->radFunc.getSize());
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int nl = (int)in(iSpecies).getSize();
      for (int l = 0; l < nl; l++)  {
         // check dimensions
         SX_CHECK (in(iSpecies)(l).getSize() ==
                   radBasisPtrIn->radFunc(iSpecies).getSize(),
                   in(iSpecies)(l).getSize(),
                   radBasisPtrIn->radFunc(iSpecies).getSize());
         // check basis
         if (muSet(iSpecies)(l).getBasisPtr () != radPtr)  {
            muSet(iSpecies)(l) = muSet(iSpecies)(l).getCopy ();
            muSet(iSpecies)(l).setBasis (radPtr);
         }
      }
   }

   createFuncLMap ();

   registerMemoryObservers ();
}

SxAtomicOrbitals::SxAtomicOrbitals (
      const SxArray<SxDiracMat<Double> > &in,
      const SxConstPtr<SxRadBasis> radBasisPtrIn)
   : SxPsiSet (PW)
{
   // --- check dimensions
   int nSpecies = (int)in.getSize();
   muSet.resize(nSpecies);
   SX_CHECK (nSpecies == radBasisPtr->radFunc.getSize(),
             nSpecies, radBasisPtr->radFunc.getSize());
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int nl = (int)in(iSpecies).row(0).getSize();
      muSet(iSpecies).resize(nl);
      for (int l=0; l < nl; l++)  {
         /*SX_CHECK (in(iSpecies).colRef(l).getSize() ==
                     radBasisPtr->radFunc(iSpecies).getSize(),
                     in(iSpecies).colRef(l).getSize(),
                     radBasisPtr->.radFunc(iSpecies).getSize());*/
         muSet(iSpecies)(l).copy(in(iSpecies).colRef(l));
      }
   }
   radBasisPtr = radBasisPtrIn;

   createFuncLMap ();

   registerMemoryObservers ();
}


SxAtomicOrbitals::SxAtomicOrbitals (const SxAtomicOrbitals &in)
   : SxPsiSet(PW)
{
   (*this) = in;
}

SxAtomicOrbitals::SxAtomicOrbitals (SxBinIO &io)
{

   radBasisPtr = SxPtr<SxRadBasis>::create(io);
   try  {
      read(io);
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
   createFuncLMap ();
}

void SxAtomicOrbitals::setup (const SxSymbolTable *table)
{
   
   SX_CHECK (table);
   radBasisPtr = SxPtr<SxRadBasis>::create(table);
   const SxString specName = "species";
   const SxString orbitalName = "orbital";
   const SxSymbolTable *species, *orbital;
   int iSpecies=-1;
   SxDiracVec<Double> newOrbital;

   try {
      species = table->getGroup (specName);
      int nSpecies = species->getNItems (specName);
      if (nSpecies != radBasisPtr->getNSpecies ())
      {
         cout << endl;
         cout << "ERROR: Mismatch in number of species" << endl;
         cout << "   Radial basis has " << radBasisPtr->getNSpecies () 
              << " species." << endl;
         cout << "   Orbital setup group '" << table->getName () 
              << "' has " << nSpecies << " species." << endl;
         SX_QUIT;
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   try   {
      
      // get species
      for(species = table->getGroup (specName);
          species != NULL;
          species = species->nextSibling (specName))
      {
         int iOrbital = -1;
         iSpecies++;
         for(orbital = species->getGroup (orbitalName);
             orbital != NULL;
             orbital = orbital->nextSibling (orbitalName))
         {
            iOrbital++;
            if (orbital->containsGroup("fromPotential"))   {
               const SxSymbolTable *mode = orbital->getGroup ("fromPotential");
               try {
                  SxParser parser;
                  SxString potFile = mode->get("file")->toString();
                  SxConstPtr<SxSymbolTable> potTable = parser.read (potFile);
                  SxAtomicOrbitals fromPot;
                  if (potTable->containsGroup("pseudoPot"))   {
                     SxPseudoPot pot = SxPseudoPot(&*potTable);
                     SxPtr<SxRadBasis> radBasisPotPtr 
                        = SxPtr<SxRadBasis>::create(pot.rad, pot.logDr);
                     fromPot  
                        = SxAtomicOrbitals(pot.getPseudoPsi (), radBasisPotPtr);
                  } else if (potTable->containsGroup("pawPot"))  {
                     SxPAWPot pot = SxPAWPot(&*potTable);
                     SxPtr<SxRadBasis> radBasisPotPtr 
                        = SxPtr<SxRadBasis>::create(pot.rad,pot.logDr);
                     fromPot 
                        = SxAtomicOrbitals(pot.getPhiPS (),radBasisPotPtr);
                  } else  {
                     cout << "No known potential group found!" << endl;
                     SX_QUIT;
                  }
                  int l = mode->get("l")->toInt();
                  int iot = mode->get("iot")->toInt();
                  int is = mode->get("is")->toInt();
                  newOrbital = fromPot(is,iot).getCopy ();
                  // Check Basis
                  if (fromPot.getBasis().getPtr() != radBasisPtr.getPtr())  {
                     newOrbital 
                        = fromPot.getBasis()
                        ->changeRadBasis(radBasisPtr.getPtr (),newOrbital);
                  }
                  if (l != newOrbital.handle->auxData.l) {
                     cout << "Inconsistent iot/l combination in fromPotential!"
                          << endl;
                     cout << "l is " << l 
                          <<", should be " << newOrbital.handle->auxData.l 
                          << endl;
                     SX_QUIT;
                  }
                  newOrbital.handle->auxData.is = iSpecies;
                  if (mode->contains("scale"))  {
                     double scale = mode->get("scale")->toReal();
                     SX_CHECK(scale > 1.0);
                     SxDiracVec<Double> scaledRad 
                        = scale * radBasisPtr->radFunc(iSpecies);
                     scaledRad(0) = 0.0;
                     SxCubicSpline<SxDiracVec<Double> > spline;
                     if (l < 2)
                        spline 
                           = SxCubicSpline<SxDiracVec<Double> > (scaledRad, 
                           newOrbital,
                           SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
                     else
                        spline 
                           = SxCubicSpline<SxDiracVec<Double> > (scaledRad, 
                           newOrbital,
                           SxCubicSpline<SxDiracVec<Double> >::Hermite);
                     newOrbital = spline.getY (radBasisPtr->radFunc(iSpecies));
                     if (l > 0) newOrbital(0) = 0;
                     newOrbital.handle->auxData.is = iSpecies;
                     newOrbital.handle->auxData.ia = -1;
                     newOrbital.handle->auxData.n = -1;
                     newOrbital.handle->auxData.l = l;
                     newOrbital.handle->auxData.m = NONE_M;
                  }
                  cout << "Add Orbital from Potential: is iot l " << is 
                       << " " << iot
                       << " " << l;
                  addOrbital(newOrbital);
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else if(orbital->containsGroup("fromFile"))   {
               const SxSymbolTable *mode = orbital->getGroup ("fromFile");
               try {
                  SxString fileName = mode->get("file")->toString();
                  int l = mode->get("l")->toInt();
                  int iot = mode->get("iot")->toInt();
                  int is = mode->get("is")->toInt();
                  // read Orbitalfile & check Basis
                  SxAtomicOrbitals fileOrbitals;
                  if (mode->contains("Siesta"))  {
                     fileOrbitals.readSiesta(fileName,radBasisPtr->radFunc(is));
                     newOrbital = fileOrbitals(is,iot).getCopy ();
                     newOrbital.setBasis(&*radBasisPtr);
                  } else  {
                     SxRadBasis vecBasis;
                     vecBasis.read(fileName);
                     fileOrbitals.read(fileName);
                     int nSpeciesFile = fileOrbitals.getNSpecies ();
                     if (is >= nSpeciesFile) {
                        cout << "Did not find species "
                             << is << " in quamol file " << fileName 
                             << " with " << nSpeciesFile << " species." << endl;
                        SX_QUIT;
                     }
                     int nOrbTypesFile = fileOrbitals.getNOrbTypes(is);
                     if (iot >= nOrbTypesFile) {
                        cout << "Did not find orbitaltype iot="
                             << iot << " in quamol file " << fileName 
                             << " with " << nOrbTypesFile << " orbital types." 
                             << endl;
                        SX_QUIT;
                     }
                     if (mode->contains("rCut"))  {
                        double rCut = mode->get("rCut")->toReal();
                        double rMin = vecBasis.getRMin(is);
                        double rMax = vecBasis.getRMax(is);
                        int nPoints = int ((double)vecBasis.getNPoints (is)
                                           * rCut / rMax);
                        SxRadRBasis cutRad (rMin,rCut,nPoints);
                        newOrbital = fileOrbitals(is,iot).getCopy ();
                        newOrbital = vecBasis.toRadRBasis (&cutRad, newOrbital);
                        newOrbital = cutRad.toRadBasis (&vecBasis, newOrbital);
                        newOrbital.setBasis(vecBasis);
                     } else {
                        newOrbital = fileOrbitals(is,iot).getCopy ();
                        newOrbital.setBasis(vecBasis);
                     }

                     if (fabs (vecBasis.radFunc(is)(0) 
                              - radBasisPtr->radFunc(iSpecies)(0)) > 1e-6 || 
                         fabs (vecBasis.logDr(is) 
                               - radBasisPtr->logDr(iSpecies)) > 1e-6)
                        newOrbital = vecBasis.changeRadBasis(
                              radBasisPtr.getPtr (),
                              newOrbital);
                  }
                  if (l != newOrbital.handle->auxData.l) {
                     cout << "Inconsistent iot/l combination in fromFile!" 
                          << endl;
                     cout << "l is " << l 
                          <<", should be " << newOrbital.handle->auxData.l 
                          << endl;
                     SX_QUIT;
                  }
                  newOrbital.handle->auxData.is = iSpecies;
                  if (mode->contains("scale"))  {
                     double scale = mode->get("scale")->toReal();
                     SX_CHECK(scale > 1.0);
                     SxDiracVec<Double> scaledRad 
                        = scale * radBasisPtr->radFunc(iSpecies);
                     scaledRad(0) = 0.0;
                     SxCubicSpline<SxDiracVec<Double> > spline;
                     if (l < 2)
                        spline
                           = SxCubicSpline<SxDiracVec<Double> > (scaledRad, 
                           newOrbital,
                           SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
                     else
                        spline
                           = SxCubicSpline<SxDiracVec<Double> > (scaledRad, 
                           newOrbital,
                           SxCubicSpline<SxDiracVec<Double> >::Hermite);
                     newOrbital = spline.getY (radBasisPtr->radFunc(iSpecies));
                     if (l > 0) newOrbital(0) = 0;
                     newOrbital.handle->auxData.is = iSpecies;
                     newOrbital.handle->auxData.ia = -1;
                     newOrbital.handle->auxData.n = -1;
                     newOrbital.handle->auxData.l = l;
                     newOrbital.handle->auxData.m = NONE_M;
                  }
                  cout << "Add Orbital from File: is iot l " << is 
                       << " " << iot
                       << " " << l;
                  addOrbital(newOrbital);
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else if(orbital->containsGroup("fromWaves"))  {
               const SxSymbolTable *mode = orbital->getGroup ("fromWaves");
               try {
                  SxString fileName = mode->get("file")->toString();
                  int l = mode->get("l")->toInt();
                  int iState = mode->get("iState")->toInt();
                  newOrbital = compressWave(fileName,iState,iSpecies,l);
                  newOrbital.handle->auxData.is = iSpecies;
                  newOrbital.handle->auxData.ia = -1;
                  newOrbital.handle->auxData.n = -1;
                  newOrbital.handle->auxData.l = l;
                  newOrbital.handle->auxData.m = NONE_M;
                  cout << "Add Orbital from wavestate: iState l "  
                  << iState << " " << l;
                  addOrbital(newOrbital);
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else if(orbital->containsGroup("generateGaussian"))  {
               const SxSymbolTable *mode 
                  = orbital->getGroup ("generateGaussian");
               try {
                  double beta = mode->get("beta")->toReal();
                  int l = mode->get("l")->toInt();
                  if (l > 0) beta = beta * sqrt(2./(1.* l));
                  const SxDiracVec<Double> &r = radBasisPtr->radFunc(iSpecies);
                  SxDiracVec<Double> rPowL = 1.0 * r;
                  if (l == 0) rPowL.set(1.0);
                  for (int il = 2; il <= l; il++) rPowL *= r;
                  newOrbital = rPowL * exp(-beta * r.sqr ());
                  newOrbital.handle->auxData.is = iSpecies;
                  newOrbital.handle->auxData.ia = -1;
                  newOrbital.handle->auxData.n = -1;
                  newOrbital.handle->auxData.l = l;
                  newOrbital.handle->auxData.m = NONE_M;
                  newOrbital.setBasis(radBasisPtr.getConstPtr());
                  if (l == 0)  {
                     cout << "Generate Gaussian Orbital: is l beta ";
                     cout << iSpecies; 
                     cout << " " << l;
                     cout << " " << beta;
                  } else  { 
                     cout << "Generate Gaussian Orbital: is l rMax "; 
                     cout << iSpecies; 
                     cout << " " << l;
                     cout << " " << beta / sqrt(2./(1.* l));
                  }
                  addOrbital(newOrbital);
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else if(orbital->containsGroup("gaussianBasis"))   {
               const SxSymbolTable *mode = orbital->getGroup ("gaussianBasis");
               try {
                  SxVector<Double> exps (mode->get("exponents")->toList());
                  int nGauss = (int)exps.getSize ();
                  SxVector<Double> coeffs (mode->get("coefficients")->toList());
                  int l = mode->get("l")->toInt();
                  if (nGauss != coeffs.getSize())  {
                     cout << "Wrong Number of coefficients: nGaussians is " 
                          << nGauss << ", but have " << coeffs.getSize () 
                          << "coefficients!" << endl;
                     SX_QUIT;
                  }
                  const SxDiracVec<Double> &r = radBasisPtr->radFunc(iSpecies);
                  int dim = (int)r.getSize();
                  newOrbital.resize (dim);
                  newOrbital.set(0.0);
                  for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
                     newOrbital 
                        += coeffs(iGauss) * exp(-r.sqr () * exps(iGauss));
                  }
                  newOrbital *= pow(r, double(l));
                  newOrbital.handle->auxData.is = iSpecies;
                  newOrbital.handle->auxData.ia = -1;
                  newOrbital.handle->auxData.n = -1;
                  newOrbital.handle->auxData.l = l;
                  newOrbital.handle->auxData.m = NONE_M;
                  newOrbital.setBasis(radBasisPtr.getConstPtr());
                  cout << "Generate orbital via gaussian basis: is l nGauss "; 
                  cout << iSpecies; 
                  cout << " " << l;
                  cout << " " << nGauss;
                  addOrbital(newOrbital);
               } catch (SxException e)   {
                  e.print ();
                  SX_EXIT;
               }
            } else {
               cout << "No known Initialization shema found!" << endl;
               cout << "SPHInX quits here!" << endl;
               SX_QUIT;
            }
         }
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }

   createFuncLMap ();
}

void SxAtomicOrbitals::addOrbital (const SxDiracVec<Double> &orbitalIN)
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   
   int nSpeciesOld = getNSpecies ();
   int nSpeciesNew = nSpeciesOld;
   int SpeciesIN = orbitalIN.handle->auxData.is;
   if (nSpeciesOld <= SpeciesIN) nSpeciesNew = nSpeciesOld+1;
   int lIN =  orbitalIN.handle->auxData.l;
   SxArray<SxArray<SxDiracVec<Double> > > muSetNew(nSpeciesNew);
   bool insert = false;
   for (int iSpecies = 0; iSpecies < nSpeciesOld; iSpecies++) {
      int nOrbTypes = getNOrbTypes(iSpecies);
      if (iSpecies == SpeciesIN) muSetNew(iSpecies).resize(nOrbTypes + 1);
      else muSetNew(iSpecies).resize(nOrbTypes);
      int iot = 0;
      for (int iotOrig = 0; iotOrig < nOrbTypes; iotOrig++)   {
         int l = muSet(iSpecies)(iotOrig).handle->auxData.l;
         // insert new Orbital if l < lMax
         if (!insert && 
             iSpecies == SpeciesIN && 
             l > lIN ) {
            muSetNew(iSpecies)(iot) = 1.0 * orbitalIN;
            muSetNew(iSpecies)(iot).handle->auxData.is = iSpecies;
            muSetNew(iSpecies)(iot).handle->auxData.ia = -1;
            muSetNew(iSpecies)(iot).handle->auxData.n = iot; 
            muSetNew(iSpecies)(iot).handle->auxData.l = lIN;
            muSetNew(iSpecies)(iot).handle->auxData.m = NONE_M;
            muSetNew(iSpecies)(iot).setBasis(radBasisPtr.getConstPtr ());
            cout << " as " << muSetNew(iSpecies)(iot).handle->auxData.is
              << " " << muSetNew(iSpecies)(iot).handle->auxData.n
              << " " << muSetNew(iSpecies)(iot).handle->auxData.l << endl;

            iot++;
            insert = true;
         }
         muSetNew(iSpecies)(iot) = 1.0 * muSet(iSpecies)(iotOrig);
         muSetNew(iSpecies)(iot).handle->auxData.n = iot;
         iot++;
      }
      // or append it if l = lMax
      if (iSpecies == SpeciesIN && !insert)   {
         muSetNew(iSpecies)(iot) = 1.0 * orbitalIN;
         muSetNew(iSpecies)(iot).handle->auxData.is = iSpecies;
         muSetNew(iSpecies)(iot).handle->auxData.ia = -1;
         muSetNew(iSpecies)(iot).handle->auxData.n = iot; 
         muSetNew(iSpecies)(iot).handle->auxData.l = lIN;
         muSetNew(iSpecies)(iot).handle->auxData.m = NONE_M;
         muSetNew(iSpecies)(iot).setBasis(radBasisPtr.getConstPtr());
         cout << " as " << muSetNew(iSpecies)(iot).handle->auxData.is
              << " " << muSetNew(iSpecies)(iot).handle->auxData.n
              << " " << muSetNew(iSpecies)(iot).handle->auxData.l << endl;
         insert = true;
         iot++;
      }
   }

   if (!insert)   {
      muSetNew(nSpeciesNew-1).resize(1);
      muSetNew(nSpeciesNew-1)(0) = 1.0 * orbitalIN;
      muSetNew(nSpeciesNew-1)(0).handle->auxData.is = nSpeciesNew-1;
      muSetNew(nSpeciesNew-1)(0).handle->auxData.ia = -1;
      muSetNew(nSpeciesNew-1)(0).handle->auxData.n = 0; 
      muSetNew(nSpeciesNew-1)(0).handle->auxData.l = lIN;
      muSetNew(nSpeciesNew-1)(0).handle->auxData.m = NONE_M;
      muSetNew(nSpeciesNew-1)(0).setBasis(radBasisPtr.getConstPtr());
      cout << " as " << muSetNew(nSpeciesNew-1)(0).handle->auxData.is
        << " " << muSetNew(nSpeciesNew-1)(0).handle->auxData.n
        << " " << muSetNew(nSpeciesNew-1)(0).handle->auxData.l << endl;
      insert = true;
   }
   *this = SxAtomicOrbitals(muSetNew,radBasisPtr);
}

SxAtomicOrbitals::~SxAtomicOrbitals ()
{
   // empty
}

void SxAtomicOrbitals::operator= (const SxAtomicOrbitals &in)
{
   if (&in == this) {
      cout << "a = a Error" << endl;
      SX_EXIT;
   }

   radBasisPtr = in.radBasisPtr;
   muSet = in.muSet;
   funcLMap = in.funcLMap; 
}

void SxAtomicOrbitals::operator+= (const SxAtomicOrbitals &in)
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   SX_CHECK(radBasisPtr == in.getRadBasisPtr());
   SX_CHECK(muSet.getSize() == in.muSet.getSize(),
            muSet.getSize(),
            in.muSet.getSize());
   for (int is = 0; is < in.muSet.getSize (); is++)   {
      SX_CHECK(muSet(is).getSize() == in.muSet(is).getSize(),
               muSet(is).getSize(),
               in.muSet(is).getSize());
      for (int l = 0; l < in.muSet(is).getSize (); l++)   {
         muSet(is)(l) += in.muSet(is)(l);
      }
   }
}

SxAtomicOrbitals SxAtomicOrbitals::operator+ (const SxAtomicOrbitals &in) const
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   SX_CHECK(radBasisPtr == in.getRadBasisPtr());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes,
               in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = in.muSet(is)(iot) + muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}
SxAtomicOrbitals SxAtomicOrbitals::operator- (const SxAtomicOrbitals &in) const
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   SX_CHECK(radBasisPtr == in.getRadBasisPtr());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes,
               in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) - in.muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}

SxAtomicOrbitals SxAtomicOrbitals::operator* (double skalar) const
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = skalar * muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}

SxAtomicOrbitals SxAtomicOrbitals::operator* (const SxArray<SxArray<double> > &skalar) const
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = skalar(is)(iot) * muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}

SxAtomicOrbitals SxAtomicOrbitals::operator* (const SxAtomicOrbitals &in) const
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   SX_CHECK(radBasisPtr == in.radBasisPtr);
   int nSpecies = getNSpecies ();
   SX_CHECK (nSpecies == in.getNSpecies ());
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      SX_CHECK (nOrbTypes == in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = in.muSet(is)(iot) * muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}

SxRadBasis::TPsi & SxAtomicOrbitals::operator() (int is, int iot)
{
   return muSet(is)(iot);
}

const SxRadBasis::TPsi & SxAtomicOrbitals::operator() (int is, int iot) const
{
   return muSet(is)(iot);
}

int SxAtomicOrbitals::getIOT (int is, int n, int l) const
{
   int iot = 0;
   while ((muSet(is)(iot).handle->auxData.n != n) 
         || (muSet(is)(iot).handle->auxData.l != l))   {
      iot++;
      if (iot >= muSet(is).getSize()) SX_EXIT;
   }

   return iot;
}


SxRadBasis::TPsi SxAtomicOrbitals::operator() (int is,
                                               int ia, 
                                               int n, 
                                               int l,
                                               int m) const
{
   int iot = getIOT (is,n,l);

   SxRadBasis::TPsi vec ( muSet(is)(iot) );
   vec.setBasis (radBasisPtr.getConstPtr ());

   vec.handle->auxData.is = is; vec.handle->auxData.ia = ia;
   vec.handle->auxData.n  = n;  vec.handle->auxData.l  = l; 
   vec.handle->auxData.m  = m;

   return vec;
}

int SxAtomicOrbitals::getNSpecies () const
{
   return (int)muSet.getSize ();
}

int SxAtomicOrbitals::getNOrbTypes (int iSpecies) const
{
   return (int)muSet(iSpecies).getSize ();
}

int SxAtomicOrbitals::getNOrbTypes () const
{
   int result = 0;
   for (int is = 0; is < getNSpecies (); is++)
      result += getNOrbTypes(is);

   return result;
}

void SxAtomicOrbitals::set(double val)
{

   for (int iSpecies = 0; iSpecies < muSet.getSize(); iSpecies++)   {
      for (int iot = 0; iot < muSet(iSpecies).getSize(); iot++)   {
         muSet(iSpecies)(iot).set(val);
      }
   }
}

SxQuantumNumbers SxAtomicOrbitals::getQuantumNumbers (int iSpecies, int idx) const
{
   SX_CHECK(iSpecies < muSet.getSize(), iSpecies, muSet.getSize());
   SX_CHECK(idx < muSet(iSpecies).getSize(), idx, muSet(iSpecies).getSize());
   int n = muSet(iSpecies)(idx).handle->auxData.n;
   int l = muSet(iSpecies)(idx).handle->auxData.l;

   SxQuantumNumbers result (iSpecies,n,l,0);

   return result;
}

SxArray<SxQuantumNumbers> SxAtomicOrbitals::getReducedOrbitalMap () const
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


SxArray<SxQuantumNumbers> SxAtomicOrbitals::getOrbitalMap (
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

int SxAtomicOrbitals::getNOrbitals (
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

int SxAtomicOrbitals::getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const
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


void SxAtomicOrbitals::createFuncLMap ()
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

SxRadBasis::TPsi SxAtomicOrbitals::getFuncL (int is, int l, int ifl)
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

const SxRadBasis::TPsi SxAtomicOrbitals::getFuncL (int is, int l, int ifl) const
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}


int SxAtomicOrbitals::getLMax () const
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

int SxAtomicOrbitals::getLMax (const int iSpecies) const
{
   int nOrbTypes = (int)muSet(iSpecies).getSize ();
   int result = 0;
   for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(iSpecies)(iot).handle->auxData.l)
            result = muSet(iSpecies)(iot).handle->auxData.l;
   }
   return result;
}

SxArray<SxArray<int> > SxAtomicOrbitals::getFuncPerL () const
{
   int nSpecies = (int)muSet.getSize ();
   SxArray<SxArray<int> > result(nSpecies);
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int lMax = getLMax (iSpecies);
      result(iSpecies).resize(lMax + 1);
      for(int l = 0; l <= lMax; l++)   {
         result(iSpecies)(l) = getFuncPerL(iSpecies,l);
      }
   }
   return result;   
}

int SxAtomicOrbitals::getFuncPerL (int is, int l) const
{
   return (int)funcLMap(is)(l).getSize();
}

void SxAtomicOrbitals::registerMemoryObservers ()
{
   TRACK_MEMORY (muSet);
}

void SxAtomicOrbitals::write (SxBinIO &io) const
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   const SxRadBasis &radBasis = *radBasisPtr;
   int nSpecies = getNSpecies();

   try {
      //create dimensions
      io.addDimension ("nSpecies", nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         io.addDimension ("dimOrbTypes-"+SxString(is), getNOrbTypes(is));
         io.addDimension ("dimRad-"+SxString(is),
                          (int)radBasis.radFunc(is).getSize());
      }
      //write data
      SxArray<SxVector<Int> > lNumbers (nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         SxString radialName= "radFunc-" + SxString(is);
         SxString dimRadName = "dimRad-" + SxString(is);
         io.write(radialName, radBasis.radFunc(is), dimRadName);
         SxString logDrName= "logDr-" + SxString(is);
         io.write(logDrName, radBasis.logDr(is));
         lNumbers(is).resize(getNOrbTypes(is));
         for (int iot=0; iot < getNOrbTypes(is); iot++)   {
            radialName= "radial-" + SxString(is) + "-" + SxString(iot);
            io.write(radialName, muSet(is)(iot), dimRadName);
            lNumbers(is)(iot) = muSet(is)(iot).handle->auxData.l;
         }
         SxString lName= "lNumbers-" + SxString(is);
         io.write(lName, lNumbers(is),"dimOrbTypes-"+SxString(is));
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
}

void SxAtomicOrbitals::write (SxString filename) const
{
   SxBinIO io(filename, SxBinIO::BINARY_WRITE_ONLY);
   write(io);
   io.setMode (SxBinIO::WRITE_DATA);
   write(io);
   io.close();
}

void SxAtomicOrbitals::read (SxBinIO &io)
{

   try {
      SxPtr<SxRadBasis> basisPtr = SxPtr<SxRadBasis>::create(io);
      radBasisPtr = basisPtr;
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      muSet.resize(nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbTypes = io.getDimension ("dimOrbTypes-"+SxString(is));
         muSet(is).resize(nOrbTypes);
         int dimRad = io.getDimension ("dimRad-"+SxString(is));
         SxDiracVec<Int> lVec (nOrbTypes);
         SxString lName= "lNumbers-" + SxString(is);
         io.read(lName,&lVec,nOrbTypes);  
         // TODO further check of SxRadBasis
         for (int iot=0; iot < nOrbTypes; iot++)   {
            muSet(is)(iot).resize(dimRad);
            SxString radialName= "radial-" + SxString(is) + "-" + SxString(iot);
            SxDiracVec<Double> &readVec = muSet(is)(iot);
            io.read(radialName,&readVec,dimRad);
            muSet(is)(iot).setBasis (radBasisPtr.getConstPtr());
            muSet(is)(iot).handle->auxData.is = is;
            muSet(is)(iot).handle->auxData.n = iot;
            muSet(is)(iot).handle->auxData.l = lVec(iot);
            muSet(is)(iot).handle->auxData.m = 0;
         }
         // TODO Inconsistent check --- is in file and is in RadBasis might differ
         /*
         if (dimRad != rBasis.radFunc(is).getSize())   {
            cout << "WARNING: Incompatible grid size in readAtomicOrbitals !" << endl;
            cout << "Grid size in File is " << dimRad << endl;
            cout << "Grid size in rad Basis is " 
                 << rBasis.radFunc(is).getSize() << endl;
            cout << "Interpolate now!" << endl;
            SxRadBasis fileBasis;
            fileBasis.read(io);
            for (int iot=0; iot < nOrbTypes; iot++)   { 
               SxNaturalCubicSpline interpol (toVector(fileBasis.radFunc(is)), 
                     toVector(muSet(is)(iot)));
               muSet(is)(iot).resize(rBasis.radFunc(is).getSize ());
               for (int i = 0; i < rBasis.radFunc(is).getSize (); i++)   {
                  int last = fileBasis.radFunc(is).getSize() - 1;
                  if (fileBasis.radFunc(is)(last) > rBasis.radFunc(is)(i))
                     muSet(is)(iot)(i) = interpol.getVal(rBasis.radFunc(is)(i));
                  else muSet(is)(iot)(i) = 0.0;
               }
            }
         }
         */
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }

   createFuncLMap ();
}

void SxAtomicOrbitals::read (const SxString &file)
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
void SxAtomicOrbitals::writeSiesta (const SxArray<SxString> &file, 
                                    const SxArray<SxDiracVec<Double> > &basis)
{
   int nSpecies = getNSpecies ();
   SX_CHECK (basis.getSize () == nSpecies);
   SX_CHECK (file.getSize () == nSpecies);

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      SxBinIO io;
      io.open(file(iSpecies), SxBinIO::ASCII_WRITE_ONLY);
      int oldL = -1;
      int z = 0;
      for (int iot = 0; iot < getNOrbTypes (iSpecies); iot++)  {
         // orbital header
         int l = muSet(iSpecies)(iot).handle->auxData.l;
         if (l == oldL)  z++;
         else {
            z = 0;
            oldL = l;
         }
         fprintf(io.fp,"%i ? %i ? 0.000000 #orbital l, n, z, is_polarized, population\n",l,z);
         int npts = (int)basis(iSpecies).getSize ();
         double delta = basis(iSpecies)(1) - basis(iSpecies)(0);
         double cutoff = basis(iSpecies)(npts-1);
         fprintf(io.fp,"%i %.12e %.12f # npts, delta, cutoff\n",npts, delta, cutoff);
         SxCubicSpline<SxDiracVec<Double> > spline (
               radBasisPtr->radFunc(iSpecies), 
               muSet(iSpecies)(iot), 
               SxCubicSpline<SxDiracVec<Double> >::Natural);
         SxDiracVec<Double> orbital = spline.getY(basis(iSpecies));
         orbital(npts-1) = 0.0;
         if (l != 0) orbital(0) = 0.0;
         // data
         io.writeXYPlot (basis(iSpecies), orbital);
      }
      io.close ();
   }
}


void SxAtomicOrbitals::readSiesta (const SxString &file, 
                                   const SxDiracVec<Double> &radFunc)
{
   SxString ionFile;
   SxList<SxString> cLine;
   SxString data;
   muSet.resize(1);

   try {
      ionFile = SxFileIO::readBinary (file,-1);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxList<SxString> ionFileTok = ionFile.tokenize ('\n');
   int line = 0;
   while (ionFileTok(line).contains("# Lmax for basis, no. of nl orbitals") == 0) 
      if (line >= ionFileTok.getSize()) {
         cout << "Expression not found" << endl;
         SX_EXIT;
      }
   data = ionFileTok(line).left("#").stripWhiteSpace ();
   //int lMax = (data.left(" ").stripWhiteSpace ()).toInt ();
   int nOrbTypes = (data.right(" ").stripWhiteSpace ()).toInt ();
   muSet(0).resize(nOrbTypes);
   for (int iot = 0; iot < nOrbTypes; iot++)  {
      while (ionFileTok(line).contains("#orbital") == 0) {
         line++;
         if (line >= ionFileTok.getSize()) {
            cout << "Expression not found" << endl;
            SX_EXIT;
         }
      }
      cLine = ionFileTok(line).left('#').stripWhiteSpace ().tokenize(' ');
      int l = cLine.first ().toInt ();
      line++;
      cLine = ionFileTok(line).left('#').stripWhiteSpace ().tokenize(' ');
      int nPoints = cLine.first ().toInt ();
      //int cutoff = cLine.last ().toInt ();
      SxDiracVec<Double> x (nPoints);
      SxDiracVec<Double> y (nPoints);
      line++;
      for (int iPoint = 0; iPoint < nPoints; iPoint++,line++)  {
         x(iPoint) = ionFileTok(line).stripWhiteSpace ().left (' ')
            .stripWhiteSpace ().toDouble();
         y(iPoint) = ionFileTok(line).stripWhiteSpace ().right(' ')
            .stripWhiteSpace ().toDouble();
      }
      SxCubicSpline<SxDiracVec<Double> > spline(x,y,
            SxCubicSpline<SxDiracVec<Double> >::Natural);
      muSet(0)(iot) = spline.getY(radFunc);
      // Set auxdata
      SxRadBasis::TPsi &mu = muSet(0)(iot);
      mu.handle->auxData.is = -1;
      mu.handle->auxData.ia = -1;
      mu.handle->auxData.n  = iot;
      mu.handle->auxData.l  = l;
      mu.handle->auxData.m  = NONE_M;
   }

   createFuncLMap ();

   registerMemoryObservers ();
}

void SxAtomicOrbitals::print (SxString file) const
{
   int nSpecies = getNSpecies();
   const SxRadBasis &rad = *radBasisPtr;
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         try {
         SxString outFile = file + SxString(is) + SxString(iot) + ".dat";
         SxBinIO out; 
         out.open(outFile, SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(toVector(rad.radFunc(is)),toVector(muSet(is)(iot)));
         out.close();
         } 
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }
   }
}

void SxAtomicOrbitals::print () const
{
   SxString file = "AtomicOrbitals";
   print (file);
}

void SxAtomicOrbitals::setBasis (SxConstPtr<SxRadBasis> radBasisPtrIn)
{
   radBasisPtr = radBasisPtrIn;
   for (int is = 0; is < muSet.getSize (); ++is)  {
      for (int iot = 0; iot < muSet(is).getSize (); ++iot)  {
         muSet(is)(iot).setBasis (radBasisPtr.getConstPtr ());
      }
   }
}

void SxAtomicOrbitals::readOrbital (FILE *fp, int is, int n, int l, int iot,
                                    bool ignoreComments)
{
   SX_CHECK (fp);
   SX_CHECK (is >= 0 && (   (is <= getNSpecies () && iot == -1)
             || (is < getNSpecies ())),
             is, getNSpecies (), iot);
   SX_CHECK (l >= 0, l);
   SX_CHECK ((iot >= 0 && iot < muSet(is).getSize ()) || (iot == -1),
             iot, muSet(is).getSize ());
   
   // --- read the file
   int N = 0;
   char c;
   SxStack<double> rFile, psiFile;
   while (!feof(fp))  {
      if (ignoreComments)  {
         while (fscanf (fp, " #%c", &c) == 1)  {
            // comment line
            while (c != '\n')  {
               if (feof(fp)) break;
               c = (char)fgetc (fp);
            }
         }
      }
      double r, f;
      if (fscanf(fp, "%lf %lf", &r, &f) == 2)  {
         rFile << r;
         psiFile << f;
         N++;
      } else {
         break;
      }
   }
   if (N == 0)  {
      cout << "Could not read orbital from text file!" << endl;
      SX_EXIT;
   }

   // --- transform data to vectors
   SxDiracVec<Double> r(rFile), psi(psiFile);

   // --- change radial grid?
   if (radBasisPtr)  {
      if (N != radBasisPtr->radFunc(is).getSize ()
          || (radBasisPtr->radFunc(is) - r).normSqr () > 1e-10)
      {
         // --- spline interpolation
         SxNaturalCubicSpline spline;
         if (l == 0)
            spline = SxNaturalCubicSpline (toVector(r), toVector(r * r * psi));
         else
            spline = SxNaturalCubicSpline (toVector(r), toVector(psi));
         // linear extrapolation below r0
         double rmax = r(r.getSize () - 1), rmin = r(0);
         double dpsi0 = (psi(1) - psi(0)) / (r(1) - r(0)),
                psi0 = psi(0);

         // --- get psi on rBasis radial grid
         psi = SxDiracVec<Double>((*radBasisPtr)(is));
         for (int ir = 0; ir < radBasisPtr->radFunc(is).getSize (); ++ir)  {
            double rr = radBasisPtr->radFunc(is)(ir);
            if (rr > rmax)  {
               psi(ir) = 0.;
            } else if (rr < rmin)  {
               psi(ir) = psi0 + (rr - rmin) * dpsi0;
            } else {
               psi(ir) = spline.getVal (rr);
               if (l == 0) psi(ir) /= rr * rr;
            }
         }
         VALIDATE_VECTOR (psi);
      }
   }

   // --- metadata
   if (iot != -1 && muSet(is)(iot).getSize () > 0)
      psi.handle->auxData = muSet(is)(iot).handle->auxData;
   psi.handle->auxData.is = is;
   psi.handle->auxData.n = n;
   psi.handle->auxData.l = l;

   // --- store away data in right place
   if (iot == -1)
      addOrbital (psi);
   else
      muSet(is)(iot) = psi;
}

void SxAtomicOrbitals::readOrbitals (FILE *fp, int is)
{
   char buffer[10240];
   int n, l;
   SxRegex lMatch("l *= *([0-9]+)");
   SxRegex nMatch("n *= *([0-9]+)");
   for ( ; !feof(fp); )  {
      buffer[0] = 0;
      fgets (buffer, 10240, fp);
      SxString line(buffer);
      line = line.stripWhiteSpace ();
      if (line.getSize () <= 1) continue;
      SxList<SxString> match;
      n = muSet.getSize () >= is ? 0 : (int)muSet(is).getSize ();
      try {
         match = lMatch.match (line);
         if (match.getSize () == 2)  {
            cout << line;
            cout << match(0) << endl;
            l = match(1).toInt ();
         } else {
            cout << match << endl;
            cout << "Missing l in '" << line << "'." << endl;
            SX_EXIT;
         }
         match = nMatch.match (line);
         if (match.getSize () == 2) n = match(1).toInt ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
      readOrbital (fp, is, n, l, -1, false);
   }
}

double SxAtomicOrbitals::getNormSqr(int is, int iot) const
{
   const SxDiracVec<Double> &vec = muSet(is)(iot);
   double logDr = radBasisPtr->logDr(is);
   const SxDiracVec<Double> &r = radBasisPtr->radFunc(is);
   double result = (vec * vec * r * r * r).integrate(logDr);
   return result;
}

double SxAtomicOrbitals::getNormSqrSum() const
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

double SxAtomicOrbitals::dot(const SxAtomicOrbitals &in) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == in.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == in.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += dot(in,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitals::dot(const SxAtomicOrbitals &in, int is, int iot, int js, int jot) const
{
   const SxDiracVec<Double> &vec1 = muSet(is)(iot);
   const SxDiracVec<Double> &vec2 = in.muSet(js)(jot);

   double logDr = radBasisPtr->logDr(is);
   const SxDiracVec<Double> &r = radBasisPtr->radFunc(is);
   SX_CHECK (logDr == in.radBasisPtr->logDr(is));
   SX_CHECK ((r(0) - in.radBasisPtr->radFunc(is)(0)) < 1e-6);
   double result = (vec1 * vec2 * r * r * r).integrate(logDr);
   return result;
}

void SxAtomicOrbitals::normalize ()
{
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         double norm = getNormSqr(is,iot);
         if (norm > 1e-6) muSet(is)(iot) /= sqrt(norm);
         else {
            cout << "WARNING: Try to normalize zero norm orbital." << endl;
            cout << "Norm is " << norm << endl;
            muSet(is)(iot).set(0.0);
         }
      }
   }
}

SxArray<SxArray<SxMatrix<Double> > > SxAtomicOrbitals::orthogonalize ()
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitals org = 1.0 * *this;
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
   return result;
}

SxArray<SxArray<SxMatrix<Double> > > SxAtomicOrbitals::getOverlap () const
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

SxMatrix<Double> SxAtomicOrbitals::getOverlap (int is, int l) const
{
   const SxDiracVec<Double> &r = radBasisPtr->radFunc(is);
   double logDr = radBasisPtr->logDr(is);
   int nFL = getFuncPerL (is,l);
   SxMatrix<Double> result(nFL,nFL);
   for (int ifl = 0; ifl < nFL; ifl++)  {
      for (int jfl = ifl; jfl < nFL; jfl++)  {
         result(ifl,jfl) = result(jfl,ifl) = 
            (getFuncL(is,l,ifl) * getFuncL(is,l,jfl)*r*r*r).integrate(logDr);
      }
   }
   return result;
}

void SxAtomicOrbitals::refine (SxArray<int> &factor)
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   const SxRadBasis &radBasis = *radBasisPtr;

   int nSpecies = getNSpecies ();
   SX_CHECK(factor.getSize() == nSpecies, factor.getSize(), nSpecies);
   SxArray<SxDiracVec<Double> > rad(nSpecies);
   SxArray<double> logDr(nSpecies);
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      SX_CHECK(factor(iSpecies) > 0, factor(iSpecies));
      if (factor(iSpecies) == 1) {
         cout << "Mesh for species " << (iSpecies + 1) 
              << " not refined." << endl;
         rad(iSpecies) = radBasis.radFunc(iSpecies).getCopy();
         logDr(iSpecies) = radBasis.logDr(iSpecies);
      }
      else  {
         // --- interpolate phi on finer mesh
         cout << "Refining mesh for species " << (iSpecies + 1) 
              << " by factor " << factor(iSpecies) << "." << endl;
         double logDrOrig = radBasis.logDr(iSpecies); 
         double r0 = radBasis.radFunc(iSpecies)(0);
         double r0New = r0;

         // --- define finer mesh
         int nrOrig = (int)radBasis.radFunc(iSpecies).getSize ();
         logDr(iSpecies) = logDrOrig / factor(iSpecies);
         int nExtra = (r0New < 0.) ? 0 : (int)(log(r0/r0New) / logDr(iSpecies));
         int nFine = (nrOrig-1)*factor(iSpecies) + 1 + nExtra;
         rad(iSpecies).resize(nFine);

         for (int i = 0; i < nFine; ++i)  {
            rad(iSpecies)(i) = r0 * exp(i * logDr(iSpecies));
         }

         // interpolate
         for (int iot = 0; iot < getNOrbTypes (iSpecies); iot++)  {
            // update rad,psi,logDr in psPot
            muSet(iSpecies)(iot) 
               = interpolateRad (muSet(iSpecies)(iot),
                     radBasis.radFunc(iSpecies)(0),
                     radBasis.logDr(iSpecies),
                     rad(iSpecies));
         }
      }
   }
   radBasisPtr = SxConstPtr<SxRadBasis>::create(rad,logDr);
   setBasis (radBasisPtr);
}

SxDiracVec<Double> SxAtomicOrbitals::compressWave(SxString &file, 
            int iState, 
            int iSpecies, 
            int l)
{
   SxPW waves;
   SxAtomicStructure structure;
   try {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      waves = SxPW (file, SxPW::InMemory);
      structure.read (io);
      io.close ();
   }
   catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   SxStack<double> xData;
   SxStack<double> yData;
   SxDiracVec<Double> yRad;
   SxDiracVec<Double> result;

   SxGkBasis &GkBasis = *waves.getGkBasisPtr ();
   GkBasis.changeTau(structure);
   double naNorm = 1.0/(structure.getNAtoms(iSpecies));
   double mNorm = 1.0/(2.0*l+1.0);
   double spinNorm = 1.0/(1.0 * waves.getNSpin ());

   for (int ik = 0; ik < waves.getNk(); ik++)  {
      yRad.resize(GkBasis(ik).g2.getSize());
      yRad.set(0.0);
      for (int iSpin = 0; iSpin < waves.getNSpin (); iSpin++)  {
         for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
            for (int m = -l; m <= l; m++)  {
               SxDiracVec<Double> wave = sqrt(waves(iState,iSpin,ik).absSqr ());
               SxDiracVec<Complex16> radial = wave
                  * GkBasis(ik).getPhaseFactors(iSpecies,iAtom).conj () //shift
                  * GkBasis(ik).getYlm(l,m) // project l
                  * SxYlm::getYlmNormFactor(l,m)
                  * sqrt(2.0 * structure.cell.volume/PI) 
                  * spinNorm
                  * naNorm 
                  * mNorm;
               yRad += 0.5 * (radial + radial.conj ());
            }
         }
      }
      for (int i = 0; i < GkBasis(ik).g2.getSize(); i++)  {
         xData.push(sqrt(GkBasis(ik).g2(i)));
         yData.push(yRad(i));
      }
   }

   SxDiracVec<Double> xToFit(xData);
   SxDiracVec<Double> yToFit(yData);

   SxCubicSpline<SxDiracVec<Double> > fit;
   SxRadGBasis radGBasis 
      = SxRadGBasis(0.0, GkBasis.getMaxGk(), 100, SxRadGBasis::Linear);
   if (l == 0) // ?(May have CUSP, no Mirror)?
      fit  = SxCubicSpline<SxDiracVec<Double> > (
            xToFit,
            yToFit,
            radGBasis.getRadGFunc(),
            SxCubicSpline<SxDiracVec<Double> >::Natural,
            SxCubicSpline<SxDiracVec<Double> >::MirrorPlane);
   else if (l & 1) // (l odd)
      fit  = SxCubicSpline<SxDiracVec<Double> > (
            xToFit,
            yToFit,
            radGBasis.getRadGFunc(),
            SxCubicSpline<SxDiracVec<Double> >::Natural,
            SxCubicSpline<SxDiracVec<Double> >::MirrorPoint);
   else // (l even and not zero)
      fit  = SxCubicSpline<SxDiracVec<Double> > (
            xToFit,
            yToFit,
            radGBasis.getRadGFunc(),
            SxCubicSpline<SxDiracVec<Double> >::Natural,
            SxCubicSpline<SxDiracVec<Double> >::MirrorPlane);
   SxDiracVec<Double> resultG = fit.getYFit ();
   resultG.setBasis(&radGBasis);
   resultG.handle->auxData.is = iSpecies;
   resultG.handle->auxData.ia = -1;
   resultG.handle->auxData.n = -1;
   resultG.handle->auxData.l = l;
   resultG.handle->auxData.m = NONE_M;

   //const SxRadBasis &radBasis = *radBasisPtr;
   //result = (radBasis | resultG);
   result = radGBasis.toRadBasis(&*radBasisPtr,resultG);

   result.setBasis(&*radBasisPtr);
   result.handle->auxData.is = iSpecies;
   result.handle->auxData.ia = -1;
   result.handle->auxData.n = -1;
   result.handle->auxData.l = l;
   result.handle->auxData.m = NONE_M;

   // ensure localization
   SxDiracVec<Double> rCut (radBasisPtr->radFunc(iSpecies).getSize ());
   rCut.set(10.0);
   result *= 1.0 / (1.0 + exp((radBasisPtr->radFunc(iSpecies) - rCut) / 1.0));

   result.normalize ();

   return result;
}

#ifdef USE_HDF5
void SxAtomicOrbitals::writeHDF5(const SxString &name)
{
   SxHDF5 file (name, SxHDF5::BINARY_CREATE);
   SxArray<SxArray<int> > funcPerL = getFuncPerL ();
   for (int is = 0; is < funcPerL.getSize (); is++)  {
      SxString isName = "Species-" + SxString(is);
      file.createGroup(isName);
      file.enterGroup(isName);
      SxString basisName = "radBasis";
      SxVector<Double> radBasis = toVector(radBasisPtr->radFunc(is));
      file.writeVector(basisName, radBasis);
      for (int l = 0; l < funcPerL(is).getSize(); l++)  {
         SxString lName = "Angularmomentum-" + SxString(l);
         file.createGroup(lName);
         file.enterGroup(lName);
         for (int il = 0; il < funcPerL(is)(l); il++)  {
            SxString radName = "Radial-" + SxString(il);
            SxVector<Double> radial = toVector(getFuncL(is,l,il));
            file.writeVector(radName, radial);
         }
         file.leaveGroup();
      }
      file.leaveGroup();
   }
}
#endif
