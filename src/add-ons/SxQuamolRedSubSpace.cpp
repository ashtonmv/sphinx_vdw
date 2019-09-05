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


#include <SxQuamolRedSubSpace.h>

#ifndef SX_STANDALONE
// Standart Constructor
SxQuamolRedSubSpace::SxQuamolRedSubSpace()
{
   GkBasisPtr = NULL;
}

void SxQuamolRedSubSpace::set (
      SxString fileName,
      SxConstPtr<SxRadBasis> radBasisPtrIn,
      SxPW &waveIn,
      SxGkBasis *GkPtrIn,
      SxConstPtr<SxRadGBasis> radGBasisPtrIn,
      SxPtr<SxOverlapBase> SPtrIn,
      SxAtomicStructure structureIn,
      double rCut)
{
   SX_CLOCK(Timer::setup);
   // set radBasis and Gk Basis
   radBasisPtr = radBasisPtrIn;
   GkBasisPtr = GkPtrIn;
   radGBasisPtr = radGBasisPtrIn;

   // set waves 
   waves = waveIn;

   // define Overlap Operator
   SPtr = SPtrIn;

   // define structure
   structure = structureIn;

   // first guess for new Basis
   improvedBasis.setBasis(radBasisPtrIn);
   improvedBasis.read(fileName);

   int lMax = improvedBasis.getLMax ();
   int nZeros = 10;
   
   // read in large Basisset
   basis = getJSB(rCut, lMax, nZeros);

   basisGk = radialsToOrbitals(basis);
   
   // get Overlap Matrix
   Sb = getOverlapB ();

   // get projections
   P = getProjections();

   // first guess for new Basis
   improvedBasis.setBasis(radBasisPtrIn);
   improvedBasis.read(fileName);


   // get coefficients
   beta = getCoefficients ();
}

// set Function
void SxQuamolRedSubSpace::set (
      SxString basisFile,
      SxAtomicOrbitals initGuess,
      SxConstPtr<SxRadBasis> radBasisPtrIn,
      SxPW &waveIn,
      SxGkBasis *GkPtrIn,
      SxConstPtr<SxRadGBasis> radGBasisPtrIn,
      SxPtr<SxOverlapBase> SPtrIn,
      SxAtomicStructure structureIn)
{
   SX_CLOCK(Timer::setup);
   // set radBasis and Gk Basis
   radBasisPtr = radBasisPtrIn;
   GkBasisPtr = GkPtrIn;
   radGBasisPtr = radGBasisPtrIn;

   // set waves 
   waves = waveIn;

   // define Overlap Operator
   SPtr = SPtrIn;

   // define structure
   structure = structureIn;
   
   // read in large Basisset
   basis.setBasis(radBasisPtrIn);
   basis.read(basisFile);

   basisGk = radialsToOrbitals(basis);
   
   // get Overlap Matrix
   Sb = getOverlapB ();

   // get projections
   P = getProjections();

   // first guess for new Basis
   improvedBasis.setBasis(radBasisPtrIn);
   improvedBasis = initGuess;

   checkInitialGuess();

   // get coefficients
   beta = getCoefficients ();
}

// Standart Destructor
SxQuamolRedSubSpace::~SxQuamolRedSubSpace()
{
   //empty
}

void SxQuamolRedSubSpace::info (SxString fileName, SxConstPtr<SxRadBasis> radBasisPtrIn)
{
   SxAtomicOrbitals quamols;
   quamols.setBasis(radBasisPtrIn);
   quamols.read(fileName);

   cout << SX_SEPARATOR;
   SxArray<SxArray<int> > funcPerL = quamols.getFuncPerL ();
   ssize_t nSpecies = funcPerL.getSize ();
   for(ssize_t is = 0; is < nSpecies; is++)   {
      ssize_t nL = funcPerL(is).getSize ();
      for(ssize_t iL = 0; iL < nL; iL++)   {
         sxprintf ("Species %d has %d functions in l-Channel %d.\n",
                 int(is), funcPerL(is)(iL), int(iL));
      }
   }
   cout << SX_SEPARATOR;
}

void SxQuamolRedSubSpace::checkInitialGuess()
{
   SxArray<SxArray<int> > funcPerLB = basis.getFuncPerL ();
   SxArray<SxArray<int> > funcPerLGuess = improvedBasis.getFuncPerL (); 
   ssize_t nSpecies = basis.getNSpecies ();
   if (nSpecies != funcPerLGuess.getSize() )   {
      cout << "Inconsistent number of species for improved Basis!" << endl;
      SX_QUIT;
   }
   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for(ssize_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      ssize_t nChannel = funcPerLB(iSpecies).getSize();
      if (nChannel < funcPerLGuess(iSpecies).getSize() )   {
         cout << "Inconsistent number of l-Channels for improved Basis!" << endl;
         SX_QUIT;
      }
   }
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::getOverlapB ()
{
   SX_CLOCK(Timer::overlapB);
   int nk = waves.getNk ();
   SxArray<SxDiracMat<Complex16> > result(nk);
   for(int ik = 0; ik < nk; ik++)   {
      SxDiracMat<Complex16> Snu = SPtr->apply(basisGk(ik));
      result(ik) = (basisGk(ik).adjoint () ^ Snu);
   }

   return result;
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::getCoefficients ()
{
   int nk = waves.getNk ();
   SxArray<SxDiracMat<Complex16> > result(nk);
   SxOrbitals orbitalsC = radialsToOrbitals(improvedBasis);

   for(int ik = 0; ik < nk; ik++)   {
      SxDiracMat<Complex16> Snu = SPtr->apply(orbitalsC(ik));
      SxDiracMat<Complex16> Pbc = (basisGk(ik).adjoint () ^ Snu);
      result(ik) = (Sb(ik).inverse () ^ Pbc);
   }

   return result;
}
SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::getOverlapC ()
{
   int nk = waves.getNk ();
   SxArray<SxDiracMat<Complex16> > result(nk);

   for(int ik = 0; ik < nk; ik++)   {
      result(ik) = (beta(ik).adjoint () ^ Sb(ik) ^ beta(ik));
   }

   return result;
}

SxArray<SxArray<SxDiracMat<Complex16> > > SxQuamolRedSubSpace::getProjections ()
{
   SX_CLOCK(Timer::projections);
   int nk = waves.getNk ();
   int nSpin = waves.getNSpin ();
   SxArray<SxArray<SxDiracMat<Complex16> > > result(nk);

   for(int ik = 0; ik < nk; ik++)   {
      result(ik).resize(nSpin);
         SxDiracMat<Complex16> Snu = SPtr->apply(basisGk(ik));
      for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
         result(ik)(iSpin) = (waves(iSpin,ik).adjoint () ^ Snu);
      }
   }
   return result;
}

// This routine calculates the Normmatrix

double SxQuamolRedSubSpace::calcNormB () 
{
   SX_CLOCK (Timer::normCalc);
   int nk = waves.getNk ();
   int nSpin = waves.getNSpin ();

   SxVector<Double> kWeight = GkBasisPtr->weights;

   double result = 0.0;
   for (int ik = 0; ik < nk; ik++)   {
      SxDiracMat<Complex16> invSb = Sb(ik).inverse();
      for (int iSpin = 0; iSpin < nSpin; iSpin++)   {
         SxComplex<double> normk = (P(ik)(iSpin) ^ invSb 
                                    ^ P(ik)(iSpin).adjoint ())
                                    .diag ().sum ();
         result += kWeight(ik) / double (nSpin) * normk.re;
      }
   }
   return result;
}

double SxQuamolRedSubSpace::calcNormC () 
{
   SX_CLOCK (Timer::normCalc);
   int nk = waves.getNk ();
   int nSpin = waves.getNSpin ();

   SxVector<Double> kWeight = GkBasisPtr->weights;

   double result = 0.0;
   SxArray<SxDiracMat<Complex16> > Sc = getOverlapC ();
   for (int ik = 0; ik < nk; ik++)   {
      SxDiracMat<Complex16> invSc = Sc(ik).inverse();
      SxDiracMat<Complex16> betaDagger = beta(ik).adjoint();
      for (int iSpin = 0; iSpin < nSpin; iSpin++)   {
         SxComplex<double> normk = (P(ik)(iSpin) ^ beta(ik) ^ invSc 
               ^ betaDagger ^ P(ik)(iSpin).adjoint ()).diag ().sum ();
         result += kWeight(ik) / double (nSpin) * normk.re;
      }
   }
   return result;
}

// Calculate Gradient
SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::calcGradient ()
{
   SX_CLOCK (Timer::gradCalc);
   SX_CHECK(GkBasisPtr);

   
   int nk = waves.getNk();
   int nSpin = waves.getNSpin();

   SxArray<SxDiracMat<Complex16> > result (nk);
   SxArray<SxDiracMat<Complex16> > Sc = getOverlapC ();

   for (int ik = 0; ik < nk; ik++)   {
      result(ik) = beta(ik).getCopy ();
      result(ik).set(0.0);
      SxDiracMat<Complex16> invSc = Sc(ik).inverse();
      SxDiracMat<Complex16> betaDagger = beta(ik).adjoint ();
      for (int iSpin = 0; iSpin < nSpin; iSpin++)   {
         SxDiracMat<Complex16> PDagger = P(ik)(iSpin).adjoint ();
         SxDiracMat<Complex16> gradk = (PDagger ^ P(ik)(iSpin) ^ beta(ik) ^ invSc)
                                   - (Sb(ik) ^ beta(ik) ^ invSc ^ betaDagger 
                                   ^ PDagger ^ P(ik)(iSpin) ^ beta(ik) ^ invSc);
         result(ik) += 1.0 / double(nSpin) * gradk;
      }
      result(ik) = keepSymmetry (result(ik));
   }

   SxDiracMat<Complex16> redResult = reduce (result);
   result = expand (redResult);
   
   return result;
}

SxDiracMat<Complex16> SxQuamolRedSubSpace::keepSymmetry (
      const SxDiracMat<Complex16> &MatIN)
{
   SxDiracMat<Complex16> result = MatIN.getCopy ();
   SxArray<SxQuantumNumbers> mapB = basis.getOrbitalMap(structure);
   SxArray<SxQuantumNumbers> mapC = improvedBasis.getOrbitalMap(structure);
   ssize_t nOrbitalsB = mapB.getSize ();
   ssize_t nOrbitalsC = mapC.getSize ();

   SX_CHECK (MatIN.colRef(0).getSize () == nOrbitalsB,
             MatIN.colRef(0).getSize (),
             nOrbitalsB);
   SX_CHECK (MatIN.row(0).getSize () == nOrbitalsC,
             MatIN.row(0).getSize (),
             nOrbitalsC);
   for (ssize_t iRow = 0; iRow < nOrbitalsB; iRow++)   {
      for (ssize_t iCol = 0; iCol < nOrbitalsC; iCol++)   {
         // No mixing between different species
         if (mapB(iRow).iSpecies != mapC(iCol).iSpecies)
            result(iRow,iCol) = 0.0;
         // No mixing between different l channels
         if (mapB(iRow).l != mapC(iCol).l)
            result(iRow,iCol) = 0.0;
      }
   }

   return result;
}


// Extract Radials to Orbitals in Gk(ik)
SxOrbitals SxQuamolRedSubSpace::radialsToOrbitals (const SxAtomicOrbitals &radialsIn)
{
   SX_CLOCK(Timer::rad2Orb);
   SX_CHECK(GkBasisPtr);
   int nk = waves.getNk ();
   SxOrbitals result(nk);
   SxGkBasis &GkBasis = *GkBasisPtr;

   int nSpecies = radialsIn.getNSpecies ();
   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = radialsIn.getNOrbTypes(is);
      muSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++) {
         muSet(is)(iot) = (*radGBasisPtr) | radialsIn(is,iot);
      }
   }

   // 0 means no spline representation
   SxAtomicOrbitalsG radialsG = SxAtomicOrbitalsG (muSet,radGBasisPtr,0);
   radialsG.toSpline ();
   
   for(int ik = 0; ik < nk; ik++)   {
      SxArray<SxQuantumNumbers> map = radialsIn.getOrbitalMap(structure);

      ssize_t gDim = GkBasis(ik).g2.getSize ();
      ssize_t nOrbitals = map.getSize ();
      result(ik).reformat(gDim,nOrbitals);

      for(int iOrbital = 0; iOrbital < nOrbitals; iOrbital++)   {
         int is = map(iOrbital).iSpecies;
         int ia = map(iOrbital).iAtom;
         int n = map(iOrbital).n;
         int l = map(iOrbital).l;
         int m = map(iOrbital).m;
         result(ik).colRef(iOrbital) <<= (GkBasis(ik)|radialsG(is,ia,n,l,m));
         result(ik).setBasis(&GkBasis(ik));
         result(ik).handle->auxData.ik = ik;
      }
   }
   return result;
}


int SxQuamolRedSubSpace::getOrbital(SxArray<SxQuantumNumbers> &map, int is, int n, int l) const
{
   int nOrbital = int(map.getSize ());
   for(int iOrbital = 0; iOrbital < nOrbital; iOrbital++)   {
      int js = map(iOrbital).iSpecies;
      int jn = map(iOrbital).n;
      int jl = map(iOrbital).l;
      if (js == is && jn == n && jl == l)   {
         return iOrbital;
      }
   }
   cout << SX_SEPARATOR;
   cout << "Orbital with " << endl 
        << "is = " << is << endl
        << "n  = " << n << endl
        << "l  = " << l << endl
        << "not found! Exit here!" << endl;
  SX_EXIT; 
}

SxDiracMat<Complex16> SxQuamolRedSubSpace::reduce (SxArray<SxDiracMat<Complex16> > &matIn) const
{
   SxArray<SxQuantumNumbers> mapB = basis.getOrbitalMap (structure);
   SxArray<SxQuantumNumbers> mapC = improvedBasis.getOrbitalMap (structure);
   SxArray<SxQuantumNumbers> redMapB = basis.getReducedOrbitalMap ();
   SxArray<SxQuantumNumbers> redMapC = improvedBasis.getReducedOrbitalMap ();

   ssize_t dimB = redMapB.getSize ();
   ssize_t dimC = redMapC.getSize ();

   SxDiracMat<Complex16> result (dimB, dimC);
   result.set(0.0);

   int nk = waves.getNk ();
   SxVector<Double> kWeight = GkBasisPtr->weights; 
   for(int ik = 0; ik < nk; ik++)   {
      SX_CHECK (mapB.getSize  () ==  matIn(ik).colRef(0).getSize (),
                mapB.getSize (), matIn(ik).colRef(0).getSize ());
      SX_CHECK (mapC.getSize () ==  matIn(ik).row(0).getSize (),
                mapC.getSize (), matIn(ik).row(0).getSize ());
      for(int iRow = 0; iRow < matIn(ik).colRef(0).getSize (); iRow++)   {
         for(int iCol = 0; iCol < matIn(ik).row(0).getSize (); iCol++)   {
            int isB = mapB(iRow).iSpecies;  int isC = mapC(iCol).iSpecies;
            int lB  = mapB(iRow).l;         int lC  = mapC(iCol).l;
            int mB  = mapB(iRow).m;         int mC  = mapC(iCol).m;
            int nB  = mapB(iRow).n;         int nC  = mapC(iCol).n;
            if (isB == isC && lB == lC && mB == mC)   {
               // TODO: Sum over all mm' or just over same mm' combination ?
               int jCol = getOrbital (redMapC, isC, nC, lC);
               int jRow = getOrbital (redMapB, isB, nB, lB);
               double aNorm = 1.0/structure.getNAtoms(isB);
               double mNorm = 1.0/(2.0 * lB + 1.0);
               result(jRow,jCol) += kWeight(ik) * aNorm * mNorm 
                  * matIn(ik)(iRow,iCol);
            }
         }
      }
   }
   return result;
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::expand (SxDiracMat<Complex16> &matIn) const
{
   SxArray<SxQuantumNumbers> mapB = basis.getOrbitalMap (structure);
   SxArray<SxQuantumNumbers> mapC = improvedBasis.getOrbitalMap (structure);
   SxArray<SxQuantumNumbers> redMapB = basis.getReducedOrbitalMap ();
   SxArray<SxQuantumNumbers> redMapC = improvedBasis.getReducedOrbitalMap ();

   ssize_t dimB = mapB.getSize ();
   ssize_t dimC = mapC.getSize ();

   int nk = waves.getNk ();

   SxArray<SxDiracMat<Complex16> > result (nk);
   result(0).reformat(dimB, dimC);
   result(0).set(0.0);

   SX_CHECK (redMapB.getSize () ==  matIn.colRef(0).getSize (),
             redMapB.getSize (), matIn.colRef(0).getSize ());
   SX_CHECK (redMapC.getSize () ==  matIn.row(0).getSize (),
             redMapC.getSize (), matIn.row(0).getSize ());
   for(int iRow = 0; iRow < mapB.getSize (); iRow++)   {
      for(int iCol = 0; iCol < mapC.getSize (); iCol++)   {
         int isB = mapB(iRow).iSpecies;  int isC = mapC(iCol).iSpecies;
         int iaB = mapB(iRow).iAtom;     int iaC = mapC(iCol).iAtom;
         int lB  = mapB(iRow).l;         int lC  = mapC(iCol).l;
         int nB  = mapB(iRow).n;         int nC  = mapC(iCol).n;
         int mB  = mapB(iRow).m;         int mC  = mapC(iCol).m;
         if (isB == isC && iaB == iaC && lB == lC && mB == mC)   {
            for (int jRow = 0; jRow < redMapB.getSize (); jRow++)   {
               for(int jCol = 0; jCol < redMapC.getSize (); jCol++)   {
                  int isBred = redMapB(jRow).iSpecies;  int isCred = redMapC(jCol).iSpecies;
                  int lBred  = redMapB(jRow).l;         int lCred  = redMapC(jCol).l;
                  int nBred  = redMapB(jRow).n;         int nCred  = redMapC(jCol).n;
                  if (isB == isBred && isC == isCred)   {
                     if (nB ==nBred && nC == nCred)   {
                        if (lC == lCred && lB == lBred)   {
                           result(0)(iRow,iCol) += matIn(jRow,jCol);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   for(int ik = 1; ik < nk; ik++)   {
      result(ik) = result(0).getCopy ();
   }

   return result;
}


void SxQuamolRedSubSpace::redBetaToRadial ()
{
   SxDiracMat<Complex16> redBeta = reduce (beta);

   SxArray<SxQuantumNumbers> redMapB = basis.getReducedOrbitalMap ();
   SxArray<SxQuantumNumbers> redMapC = improvedBasis.getReducedOrbitalMap ();

   ssize_t dimB = redMapB.getSize ();
   ssize_t dimC = redMapC.getSize ();

   // set improvedBasis to Zero
   improvedBasis = 0.0 * improvedBasis;

   for(int iRow = 0; iRow < dimB; iRow++)   {
      for(int iCol = 0; iCol < dimC; iCol++)   {
         int isB = redMapB(iRow).iSpecies, iotB = redMapB(iRow).n;
         int isC = redMapC(iCol).iSpecies, iotC = redMapC(iCol).n;
         if (isC == isB)  {
            improvedBasis(isC,iotC) += redBeta(iRow,iCol).re * basis(isB,iotB);
         }
      }
   }
}


double SxQuamolRedSubSpace::dotproductRad (const SxDiracVec<Double> &vec1, const SxDiracVec<Double> &vec2)
{
   // Check ingredients 
   SX_CHECK (vec1.handle->auxData.is == vec2.handle->auxData.is,
             vec1.handle->auxData.is,
             vec2.handle->auxData.is);
   SX_CHECK (vec1.getSize() ==  vec2.getSize(),
             vec1.getSize(),
             vec2.getSize());
   int is = vec1.handle->auxData.is;
   const SxRadBasis &radBasis = *radBasisPtr;
   int l1 = vec1.handle->auxData.l;
   int l2 = vec2.handle->auxData.l;

   double result = 0.0;

   if(l1 == l2)   {
      result = (vec1 * vec2 
             * radBasis.radFunc(is).cub()).integrate(radBasis.logDr(is));
   }
   return result;
}


SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::getCopy (
      const SxArray<SxDiracMat<Complex16> > &dataIN)
{
   ssize_t nk = dataIN.getSize ();
   SxArray<SxDiracMat<Complex16> > result(nk);

   for(ssize_t ik = 0; ik < nk; ik++)   {
      result(ik) = dataIN(ik).getCopy ();
   }

   return result;
}

void SxQuamolRedSubSpace::compute ()
{
   // Minimisation sheme CGLMin
   cout << "Start CGLMin" << endl;
   cout << SX_SEPARATOR;

   int nk = waves.getNk ();
   double norm = calcNormC ();

   SxArray<SxDiracMat<Complex16> > grad = calcGradient ();
   SxArray<SxDiracMat<Complex16> > dir = getCopy(grad);
   double oldNormAbs, newNormAbs, delta, sw = 0.5;
   int step = 0;
   restartCG = false;
   SxComplex<double> control;
   control.re = sw;
   control.im = norm;
   do   {
      step++;
      if(step == 1 || restartCG)   {
         grad = calcGradient ();
         dir = getCopy(grad);
         restartCG = false;
      } else {
         // Save old Gradient for CG
         SxArray<SxDiracMat<Complex16> > oldGrad = getCopy(grad);
         
         // Calculate new Gradient
         grad = calcGradient ();

         SxArray<SxDiracMat<Complex16> > oldDir = getCopy(dir);
         // Calculate new direction
         for(int ik = 0; ik < nk; ik++)   {
            double oldGradProd = oldGrad(ik).absSqr ().sum ();
                                                   
            double gamma = 0;
            if (oldGradProd > 1e-14)  {
               gamma =  grad(ik).absSqr ().sum () / oldGradProd;
               }
            dir(ik) = grad(ik) + gamma * oldDir(ik);
         } 
      }
      
      // Optimize
      oldNormAbs = norm / waves.getNStates();
      if (sw > 5 ) sw = 0.5;
      control = lineMin (dir, sw);
      sw = control.re;
      SxArray<SxDiracMat<Complex16> > toBeta = skalarMultBeta(sw, dir);
      setBeta (addBeta(beta, toBeta));
      norm = calcNormC (); 
      newNormAbs = norm / waves.getNStates();
      delta = newNormAbs - oldNormAbs;
 
      double spillage = 1 - newNormAbs; 
         cout << "Step: " << step
         << ", Norm: " << newNormAbs 
         << ", Delta: " << delta
         << ", Spillage: " << spillage
         << endl;
     
   }while((step < maxSteps) && (fabs(delta) > error));
   
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::getBeta ()
{
   ssize_t nk =  beta.getSize ();

   SxArray<SxDiracMat<Complex16> > result (nk);
   for(ssize_t ik = 0; ik < nk; ik++)   {
      result(ik) = beta(ik).getCopy();
   }

   return result;
}

void SxQuamolRedSubSpace::setBeta (
      const SxArray<SxDiracMat<Complex16> > &betaIN)
{
   ssize_t nk = betaIN.getSize ();

   for(ssize_t ik = 0; ik < nk; ik++)   {   
      beta(ik) = betaIN(ik).getCopy ();
   }
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::addBeta (
      const SxArray<SxDiracMat<Complex16> > &beta1,
      const SxArray<SxDiracMat<Complex16> > &beta2)
{ 
   ssize_t nk = beta1.getSize ();
   SX_CHECK (nk == beta2.getSize (), nk, beta2.getSize ());
   SxArray<SxDiracMat<Complex16> > result(nk);
   for(ssize_t ik = 0; ik < nk; ik++)   {
      SX_CHECK (beta1(ik).getSize () == beta2(ik).getSize (),
                beta1(ik).getSize (), beta2(ik).getSize ());
      result(ik) = beta1(ik) + beta2(ik);
   }
   return result;
}

SxArray<SxDiracMat<Complex16> > SxQuamolRedSubSpace::skalarMultBeta (
      const double skalar,
      const SxArray<SxDiracMat<Complex16> > &betaIN)
{
   ssize_t nk = betaIN.getSize ();
   SxArray<SxDiracMat<Complex16> > result(nk);
   for(ssize_t ik = 0; ik < nk; ik++)   {
      result(ik) = skalar * betaIN(ik);
   }
   return result;
}

SxComplex<double> SxQuamolRedSubSpace::lineMin (
      const SxArray<SxDiracMat<Complex16> > &dir, 
      double sw)
{
   // backup actual beta
   SxArray<SxDiracMat<Complex16> > backup = getBeta ();
   double normBackup = calcNormC ();

   //define points
   double x0 = 0.0;
   double N0 = normBackup;

   double x1 = sw;
   SxArray<SxDiracMat<Complex16> > toBeta = skalarMultBeta (x1, dir);
   setBeta (addBeta (backup, toBeta));
   double N1 = calcNormC ();

   double x2, N2;
   if (N1 > N0)   { 
      x2 = 2 * x1 - x0;
      toBeta = skalarMultBeta (x2, dir);
      setBeta (addBeta (backup, toBeta));
      N2 = calcNormC ();
   }
   else    {
      x2 = x1;
      N2 = N1;
      x1 = 0.5 * (x0 + x1);
      toBeta = skalarMultBeta (x1, dir);
      setBeta (addBeta (backup, toBeta));
      N1 = calcNormC ();
   }

   // parabel fit
   SxComplex<double> pFit = parabelFit(x0,x1,x2,N0,N1,N2);
   toBeta = skalarMultBeta (pFit.re, dir);
   setBeta(addBeta(backup, toBeta));
   double optNorm = calcNormC ();
   

   // Has lineminimazation a problem ?
   int lineMinSteps = 0;
    
   while ( (( optNorm < normBackup ) || 
           (fabs(optNorm - pFit.im) > 0.1 * fabs(optNorm - normBackup))) &&
           (lineMinSteps < 10) )   
   { 
      cout << "LMSteps: " << lineMinSteps << " " << pFit.re << " " << pFit.im << " " << optNorm << endl; 
      // Choose new fitting point
      // New point is old point : PROBLEM 
      if ( (fabs(pFit.re - x0) < 1e-10) || 
           (fabs(pFit.re - x1) < 1e-10) || 
           (fabs(pFit.re - x2) < 1e-10) )   {
         cout << "Identical points!" << endl;
         lineMinSteps = 11;
      }
      // New point is larger than x2 : Delete x0
      else if (pFit.re > x2)   {
         x0 = x1; N0 = N1;
         x1 = x2; N1 = N2;
         x2 = pFit.re; N2 = optNorm;
      }
      // New Point between x1 and x2 : Delete x0
      else if (pFit.re > x1)   {
         x0 = x1; N0 = N1;
         x1 = pFit.re; N1 = optNorm;
      }
      // New Point between X0 and x1: Delete x2
      else if (pFit.re > x0)   {
           x2 = x1; N2 = N1;
           x1 = pFit.re; N1 = optNorm;
      }
      //New point smaller x0: Delete x2
      else if (pFit.re < x0)   {
           x2 = x1; N2 = N1;
           x1 = x0; N1 = N0;
           x0 = pFit.re; N0 = optNorm;
      }
      //cout << SX_SEPARATOR;
      pFit = parabelFit(x0,x1,x2,N0,N1,N2);
      toBeta = skalarMultBeta (pFit.re, dir);
      setBeta(addBeta(backup, toBeta));
      optNorm = calcNormC ();

      lineMinSteps++;
      //restartCG = true;
   }

   // YES, we definitly have a problem
   if (lineMinSteps >= 10)   {
      cout << "No chance to find maximum... Scan for Maximum und restart CG" << endl;
      SxDiracVec<Double> swArray(101), N(101);
      double width = 10.0 / 100.0;
      for(int i = 0; i <= 100; i++)  {
         swArray(i) = width * i - 5.0;
         toBeta = skalarMultBeta (swArray(i), dir);
         setBeta(addBeta(backup, toBeta));
         N(i) = calcNormC ();
      }
      int i;
      N.maxval(&i);
      pFit.re = swArray(i);
      pFit.im = -N(i);
      toBeta = skalarMultBeta (pFit.re, dir);
      setBeta(addBeta(backup, toBeta));
      optNorm = calcNormC ();
      cout << pFit.re << " " << pFit.im  << " " << optNorm << endl;
      SxString file = "LineMin";
      file += ".dat";
      SxBinIO out;
      out.open(file, SxBinIO::ASCII_WRITE_ONLY);
      out.writeXYPlot(swArray,N);
      out.close();
      cout << SX_SEPARATOR;
      SX_EXIT;
   }

   setBeta(backup);

   return pFit;
}

SxComplex<double> SxQuamolRedSubSpace::parabelFit (double x0, double x1, double x2, double N0, double N1, double N2)
{
   //real part gives x Value / imag part gives y value
   SxComplex<double> result;

   double z1 = x1 - x0;
   double z2 = x2 - x0;
   
   double a = (z1 * (N0-N2) + z2 * (N1-N0))
         / (z1*z2*(z1-z2));
   double b = -(z1 * z1 * (N0-N2) + z2 * z2 * (N1-N0))
         / (z1*z2*(z1-z2));

   double opt = -0.5 * b / a;
   // Parabel has maximum? Go to it
   if (a < 0)   {
      result.re = opt + x0;
      result.im = a * opt * opt + b * opt + N0;
   } else if (a > 0) {
   // Parabel has minimum ? Go to maximal possible quantity
      double N = 0.1 * (waves.getNStates () - N0) + N0;
      opt = -b/(2.0*a) + sqrt((N-N0)/a + b*b/(4.0*a*a));
      result.re = opt + x0;
      result.im = a * opt * opt + b * opt + N0;
   } else {
      if (N0 < N1)  {
         if (N2 < N0)  {
            result.re = x2;
            result.im = N2;
         } else {
            result.re = x0;
            result.im = N0;
         }
      } else {
         if (N2 < N1)  {
            result.re = x2;
            result.im = N2;
         } else {
            result.re = x1;
            result.im = N1;
         }
      }
   }
   
   return result;
}
SxAtomicOrbitals SxQuamolRedSubSpace::getJSB (double rCut, int lMax, int nZeros)
{
   SX_CHECK(rCut > 0.);

   SxArray<SxVector<Double> > zeros = getJSBZeros(lMax, nZeros);

   /* Print roots for debug
   for(int l = 0; l <=lMax; l++)   {
      zeros(l).print();
   }
   */

   int nSpecies = int(radBasisPtr->radFunc.getSize());
   int nL = int(zeros.getSize());
   SxArray<SxArray<SxDiracVec<Double> > > muSet (nSpecies);
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int iot = 0;
      for (int l = 0; l < nL; l++)   {
         muSet(iSpecies).resize(nL * nZeros);
         for (int iZeros = 0; iZeros < nZeros; iZeros++)   {
            double alpha = zeros(l)(iZeros) / rCut;
            SxDiracVec<Double> r = radBasisPtr->radFunc(iSpecies);
            muSet(iSpecies)(iot) = jsb(l,alpha * r);
            int dim = int(muSet(iSpecies)(iot).getSize());
            for (int i = 0; i < dim; i++) 
               if (r(i) > rCut) 
                  muSet(iSpecies)(iot)(i) = 0.0; 
            muSet(iSpecies)(iot).setBasis(*radBasisPtr);
            muSet(iSpecies)(iot).handle->auxData.is = iSpecies;
            muSet(iSpecies)(iot).handle->auxData.ia = -1;
            muSet(iSpecies)(iot).handle->auxData.n = iot;
            muSet(iSpecies)(iot).handle->auxData.l = l;
            muSet(iSpecies)(iot).handle->auxData.m = NONE_M;
            iot++;
         }
      }
   }

   SxAtomicOrbitals result (muSet, radBasisPtr);

   return result;
}

SxDiracVec<Double> SxQuamolRedSubSpace::jsb (int l, const SxDiracVec<Double> &z) const
{

   SX_CHECK (z.getSize () > 0, z.getSize());

   ssize_t zSize = z.getSize ();
   SxDiracVec<TReal8> vec(zSize);
   // --- use generic jsb
   for (ssize_t i = 0; i < zSize; ++i)
      vec(i) = SxYlm::jsb(l, z(i));
   VALIDATE_VECTOR (vec);
   return vec;
}
SxArray<SxVector<Double> > SxQuamolRedSubSpace::getJSBZeros (int lMax, int nZeros)
{
   //TODO This is simply a test enviroment
   //     The routine should calculate an abitrary number JSB Zeros 
   //     for abitrary l and abitrary accuracy
   SxArray<SxList<double> > work (lMax+1);
   SxList<double>::Iterator it;
   double factor = 1.0;
   work(0).resize(nZeros);
   for (it = work(0).begin(); it != work(0).end(); it++)   {
      *it = factor * PI; // Roots from sin(x)/x are multiples of Pi
      factor += 1.0;
   }
   for (int l = 1; l <= lMax; l++)   {
      double guess = PI + l; // First root is near Pi + l
      work(l).resize(nZeros);
      for (it = work(l).begin(); it != work(l).end(); it++)   {
         *it = findRootJSB (l, guess);
         guess = *it + 3.0; // Next root is nearly 3.0 away from old one
      }
   }

   SxArray<SxVector<Double> > result(lMax+1);
   for(int l = 0; l <= lMax; l++)   {
      result(l) = SxVector<Double> (work(l));
   }

   return result;
}

double SxQuamolRedSubSpace::findRootJSB (int l, double guess)
{
   //findRoot algorithm by Mueller (parabolic)
   double p2 = guess;
   double p1 = p2 - 0.01;
   double p0 = p1 - 0.01;
   double fp2 = SxYlm::jsb(l, p2);
   double fp1 = SxYlm::jsb(l, p1);
   double fp0 = SxYlm::jsb(l, p0);
   
   int counter = 0;
   while ((fabs(fp2) > 1e-12) && counter < 1000)   {
      SX_CHECK (fabs(p0-p1) > 1e-10);
      SX_CHECK (fabs(p0-p2) > 1e-10);
      SX_CHECK (fabs(p1-p2) > 1e-10);
      double a = ((p1-p2)*(fp0-fp2)-(p0-p2)*(fp1-fp2)) 
               / ((p0-p2)*(p1-p2)*(p0-p1));
      double b = ((p0-p2)*(p0-p2)*(fp1-fp2)-(p1-p2)*(p1-p2)*(fp0-fp2))
               / ((p0-p2)*(p1-p2)*(p0-p1));
      double sigb = b > 0 ? 1.0 : -1.0;
      SX_CHECK(b*b-4*a*fp2 > 0);
      double p3 = p2 - 2*fp2 / (b+sigb*sqrt(b*b-4*a*fp2));
      p0 = p1; p1 = p2; p2 = p3;
      fp0 = fp1; fp1 = fp2;
      fp2 = SxYlm::jsb(l, p2);
      counter++;
   }
   SX_CHECK (counter < 1000);
   return p2;
}

#else // SX_STANDALONE

int main (int argc, char** argv)
{
   initSPHInXMath ();

   cout.precision(10);

   // Command line parsing
   SxCLI cli (argc,argv);

   // Define Author
   cli.authors = "B. Lange";

   // What does the program do ?
   cli.preUsageMessage = "Optimized Orbitalbasis reduction tool.";

   SxString inputFile = cli.option ("-i|--input", "file", "SPHInX input file")
                        .toString ("input.sx");
   SxString basisFile = cli.option ("-q|--quamol", "file", "Quamol input file")
                        .toString ("quamol.sxb");
   SxString guessFile = cli.option ("-g|--guess", "file", "Quamol input file")
                        .toString ("guess.sx");

   SxString wavesFile = cli.option ("-w|--waves", "file", "SPHInX waves file")
                        .toString ("waves.sxb");
   // Define several options
   double error = cli.option ("-e|--error",
                                  "Norm difference",
                                  "Convergence criterium for Norm").toDouble (1e-6,1e-10,1e-1);
   int maxSteps = cli.option ("--maxSteps",
                              "int","maximal optimization steps")
                              .toInt (100,0,1000);
   int nStates = cli.option ("--nStates",
                             "nStates",
                             "Number of Wave states for calculation")
                             .toInt ();
   bool noOpt = cli.option ("--noOpt",
                            "no optimization",
                            "no optimization is performed ").toBool();
   bool info = cli.option ("--info",
                            "info of Basis",
                            "gives information of Basis").toBool();
   bool jsbBasis = cli.option ("--jsbBasis",
                            "use jsb Basis",
                            "use spherical Besselfunctions as Basisfunctions").toBool();
   double rCut = cli.option ("--rCut",
                             "Localization cutoff",
                             "Localization cutoff").toDouble (15,5,70);
   int nPoints = cli.option ("--nPoints","Points for interpolation",
                             "Points for interpolation")
                            .toInt ();
   cli.finalize ();

    
   // --- read input file
   SxParser parser;
   SxConstPtr<SxSymbolTable> table = parser.read (inputFile);
   SxConstPtr<SxSymbolTable> guessTable = parser.read (guessFile);
   SxAtomicStructure structure (&*table);
   SxQuamolRedSubSpace redQuamol;

   if (info) {
      SxConstPtr<SxRadBasis> radBasisPtr = SxConstPtr<SxRadBasis>::create(basisFile);
      redQuamol.info(basisFile, radBasisPtr);
      return 0;
   }


   // --- Read in the PWaves and set GkBasis
   RelVec meshDim;
   SxPW waves;
   SxFermi fermi;
   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      waves.read (io);
      waves.getGkBasis().changeTau(structure);
      fermi.read (io);
      io.read ("meshDim", &meshDim);
      io.close ();

   }
   catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   SxGkBasis &GkBasis = waves.getGkBasis();

   // Generate Basis depending on the underlying |g+k| values
   SxPtr<SxRadGBasis> radGBasisPtr = 
      SxPtr<SxRadGBasis>::create(0.0, 
            GkBasis.getMaxGk (), 
            nPoints, 
            SxRadGBasis::Hyperbel);

   // Resize waves to needed value
   int nk = waves.getNk();
   SxArray<int> NperK (nk);
   for (int ik = 0; ik < nk; ik++)   {
      if (nStates > waves.getNStates(ik))   {
         sxprintf ("DFT run has %d states at maximal.\n", waves.getNStates(ik));
         sxprintf ("You choose %d states! SxQuamolRedSubSpace ends here!\n",
                 nStates);
         SX_QUIT;
      }
      NperK(ik) = nStates;
   }
   waves.setNStates(NperK);

   // Normconserving or PAW Potential ?
   if (table->containsGroup("pseudoPot"))   {
      SxConstPtr<SxRadBasis>radBasisFilePtr 
         = SxConstPtr<SxRadBasis>::create(basisFile);
      SxPtr<SxPseudoPot> psPotPtr = SxPtr<SxPseudoPot>::create (&*table);
      SxConstPtr<SxRadBasis> radBasisPotPtr 
         = SxConstPtr<SxRadBasis>::create 
         (psPotPtr->rad, psPotPtr->logDr);
      SxPtr<SxOverlapBase> SPtr = SxPtr<SxPWOverlap>::create ();
      SxAtomicOrbitals initGuess;
      try   {
         SxSymbolTable *guessGroup = guessTable->getGroup("Quamol");
         initGuess.setup(guessGroup);
      } catch (SxException e)   {
         e.print ();
         SX_EXIT;
      }
      if(jsbBasis) redQuamol.set(basisFile,
                                 radBasisFilePtr,
                                 waves,
                                 &GkBasis,
                                 radGBasisPtr,
                                 SPtr,
                                 structure,
                                 rCut);
      else         redQuamol.set(basisFile,
                                 initGuess,
                                 radBasisFilePtr,
                                 waves,
                                 &GkBasis,
                                 radGBasisPtr,
                                 SPtr,
                                 structure);
   } else if (table->containsGroup("pawPot"))   {
      SxPtr<SxRadBasis> radBasisFilePtr = SxPtr<SxRadBasis>::create(basisFile);
      SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);
      SxConstPtr<SxRadBasis> radBasisPotPtr = SxConstPtr<SxRadBasis>::create
         (pawPotPtr->rad, pawPotPtr->logDr);
         SxPtr<SxPartialWaveBasis> pBasis 
            = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
         pBasis->createProjBasis (GkBasis);
      SxPtr<SxOverlapBase> SPtr
        = SxPtr<SxPAWOverlap>::create (pBasis, pawPotPtr);
     
     SxAtomicOrbitals initGuess;
     try   {
         SxSymbolTable *guessGroup = guessTable->getGroup("Quamol");
         initGuess.setup(guessGroup);
      } catch (SxException e)   {
         e.print ();
         SX_EXIT;
      }
     if(jsbBasis) redQuamol.set(basisFile,
                                 radBasisFilePtr,
                                 waves,
                                 &GkBasis,
                                 radGBasisPtr,
                                 SPtr,
                                 structure,
                                 rCut);
     else         redQuamol.set(basisFile,
                                initGuess,
                                radBasisFilePtr,
                                waves,
                                &GkBasis,
                                radGBasisPtr,
                                SPtr,
                                structure);
   } else   {
      cout << "No known Potential Group found !" << endl;
      SX_QUIT;
   }
   
   redQuamol.error = error;
   redQuamol.maxSteps = maxSteps;

   cout << "States: " << nStates << endl;
   
   //calculate Norm
   double norm = redQuamol.calcNormC ();
   cout << endl;
   cout << "Fullnorm: " << norm / redQuamol.waves.getNStates() << endl;
   cout << "Spillage: " << 1 - norm / redQuamol.waves.getNStates() << endl;
   cout << SX_SEPARATOR;
   
   if (!noOpt) { redQuamol.compute ();}

   redQuamol.redBetaToRadial ();

   int ns = redQuamol.improvedBasis.getNSpecies ();
   for(int is = 0; is < ns; is++)   {
      int nOrbTypes = redQuamol.improvedBasis.getNOrbTypes (is);
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         SxString file = "Quamol";
         file += is;
         file += iot;
         file += ".dat";
         SxBinIO out; 
         out.open(file, SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(toVector(redQuamol.radBasisPtr->radFunc(is)),
               toVector(redQuamol.improvedBasis(is,iot)));
         out.close();
      }
   }
   
   SxString filename = "improvedQuamol.sxb";
   redQuamol.improvedBasis.write(filename);

   norm = redQuamol.calcNormC () / redQuamol.waves.getNStates();

   cout << "Final Norm: " << norm << endl;
   cout << "Final Spillage: " << 1 - norm << endl;
   double normB = redQuamol.calcNormB () / redQuamol.waves.getNStates();
   cout << "Spillage Basis Change: " << 1.0 - norm / normB << endl;  


   SxTimer::getGlobalTimer().print ();
   cout << "SxQuamolRedSubSpace completes successfully" << endl;
   return 0;
}

#endif /* SX_STANDALONE */
