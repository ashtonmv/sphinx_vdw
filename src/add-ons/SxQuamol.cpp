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


#include <SxQuamol.h>
#include <SxConstants.h>

#ifndef SX_STANDALONE
// Standart Constructor
SxQuamol::SxQuamol()
{
   // empty
}

// set Function
void SxQuamol::set (
      SxAtomicOrbitals &radialsIn,
      SxConstPtr<SxRadGBasis> radGBasisPtr,
      SxConstPtr<SxPW> wavePtrIn,
      SxConstPtr<SxFermi> fermiPtrIn,
      SxConstPtr<SxOverlapBase> SPtrIn)
{
   SX_CLOCK(Timer::setup);
   wavesPtr     = wavePtrIn;
   fermiPtr     = fermiPtrIn;
   SPtr         = SPtrIn;
   printLine = false;

   // Transform radials to radGBasis
   int nSpecies = radialsIn.getNSpecies ();
   SxArray<SxArray<SxDiracVec<Double> > > muSet (nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = radialsIn.getNOrbTypes(is);
      muSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++) {
         SxDiracVec<Double> &mu = muSet(is)(iot); 
         SxDiracVec<Double> &muIn = radialsIn(is,iot); 
         mu = ((*radGBasisPtr) | muIn);
      }
   }
   // 0 means no spline Representation
   radials = SxAtomicOrbitalsG (muSet,radGBasisPtr,0);
   radials = setOrthogonal(radials);

   functionalValue = 0.0;
   grad = 0.0 * radials;
   spillage = 0.0;
   spillageGrad = 0.0 * radials;
   eKin = 0.0;
   eKinGrad = 0.0 * radials;


   int nk = wavesPtr->getNk();
   dRgdRGk.resize(nk);
   for (int ik = 0; ik < nk; ik++) setDRDR (ik);

   ssize_t nGPoints = radGBasisPtr->getRadGFunc().getSize();
   double gMax = radGBasisPtr->getRadGFunc()(nGPoints-1);
   ssize_t nPoints = nGPoints * ssize_t (rMax / gMax);
   SxConstPtr<SxRadRBasis> radRBasisPtr 
      = SxConstPtr<SxRadRBasis>::create(0.0, rMax, (int)nPoints);

   computeLocRGGprime (radRBasisPtr);

}

// Standart Destructor
SxQuamol::~SxQuamol()
{
   //empty
}
void SxQuamol::setFixedList(const SxSymbolTable *table)
{
   SX_CHECK (table);
   const SxString specName    = "species";
   const SxString orbitalName = "orbital";
   const SxSymbolTable *species, *orbital;
   int nSpecies = radials.getNSpecies ();
   fixedOrbitals.resize(nSpecies);
   cout << SX_SEPARATOR;
   cout << "Variational Settings" << endl;
   cout << SX_SEPARATOR;
   try   {
      // get species
      int iSpecies = 0;
      for(species = table->getGroup (specName);
          species != NULL;
          species = species->nextSibling (specName))
      {
         int nOrbTypes = radials.getNOrbTypes(iSpecies);
         fixedOrbitals(iSpecies).resize(nOrbTypes);
         int iOrbital = 0;
         for(orbital = species->getGroup (orbitalName);
             orbital != NULL;
             orbital = orbital->nextSibling (orbitalName))
         {
            cout << "Species " << iSpecies 
                 << " Orbital " << iOrbital;
            if (orbital->hasAttribute ("fixed"))  {
               fixedOrbitals(iSpecies)(iOrbital) = true;
               cout << " is fixed." << endl;
            } else  {
               fixedOrbitals(iSpecies)(iOrbital) = false;
               cout << " is variational." << endl;
            }
            iOrbital++;
         }
         iSpecies++;
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
}

void SxQuamol::computeLocRGGprime (SxConstPtr<SxRadRBasis> radRBasisPtr)
{
   const SxRadRBasis &radRBasis = *radRBasisPtr;
   const SxRadGBasis &radGBasis = *radials.getRadGBasisPtr ();
   const SxDiracVec<Double> &r = radRBasis.getRadRFunc ();
   const SxDiracVec<Double> &dr = radRBasis.getRadDR ();
   const SxDiracVec<Double> &g = radGBasis.getRadGFunc ();
   SxDiracVec<Double> locR = r - rStart;
    for (int ir = 0; ir < locR.getSize (); ir++) 
      if (locR(ir) < 1e-12) locR(ir) = 0;
   SxDiracVec<Double> locR3 = locR * locR * locR;
   int lMax = radials.getLMax ();
   locRGGPrime.resize(lMax+1);
   for (int l = 0; l <= lMax; l++ ) {
      locRGGPrime(l).reformat(g.getSize(), g.getSize());
      SxDiracMat<Double> GToR (r.getSize (),g.getSize ());
      SxDiracMat<Double> locR3GToR (r.getSize (),g.getSize ());
      for (ssize_t ig = 0; ig < g.getSize (); ig++) {
         SxDiracVec<Double> delta (g.getSize());
         delta.set(0.0);
         delta(ig) = 1.0;
         delta.setBasis(radials.getRadGBasisPtr ().getPtr ());
         delta.handle->auxData.l = l;
         SxDiracVec<Double> deltaR = (radRBasis | delta);
         GToR.colRef(ig) << deltaR;
         locR3GToR.colRef(ig) << locR3 * r * r * dr * deltaR;
      }
      
      locRGGPrime(l) =  GToR.adjoint () ^ locR3GToR;
      //symmetrize
      locRGGPrime(l) = 0.5 * (locRGGPrime(l) + locRGGPrime(l).transpose ());
   }
}

void SxQuamol::setDRDR(int ik)
{
   SX_CLOCK(Timer::DRDataDRGrid);
   SX_CHECK(dRgdRGk.getSize() > ik);
   
   const SxGkBasis &gkBasis = *wavesPtr->getGkBasisPtr ();
   SxDiracVec<Double> radGFunc 
      = radials.getRadGBasisPtr ()->getRadGFunc ();
   SxCubicSpline<SxDiracVec<Double> > spline;
   spline.setXFit(sqrt(gkBasis(ik).g2));
   spline.setXVals(radGFunc);
   // for s and p orbitals (l = 0 or l = 1)
   spline.setSplineType(SxCubicSpline<SxDiracVec<Double> >::NaturalHermite);
   dRgdRGk(ik).resize(2);
   dRgdRGk(ik)(0) = spline.calcDRDataDRGrid ();
   // for d and beyond orbitals (l >= 2)
   spline.setSplineType(SxCubicSpline<SxDiracVec<Double> >::Hermite);
   dRgdRGk(ik)(1) = spline.calcDRDataDRGrid ();
}


void SxQuamol::computeFunctional (
      const SxAtomicOrbitalsG &functions, 
      bool calcGrad,
      bool calcKinLoc)
{
   SX_CLOCK(Timer::computeFunctional)

   functionalValue = 0.0;
   if (calcGrad) grad = 0.0 * grad;

   // Spillage
   computeSpillage(functions, calcGrad);
   
   functionalValue += sigma * spillage;
   if (calcGrad) grad += sigma * spillageGrad;

   // kinetic energy
   if(zeta > 1e-12 || calcKinLoc)  {
      computeKineticEnergy(functions, calcGrad);
      functionalValue += zeta * eKin;
      if (calcGrad) grad += zeta * eKinGrad;
   }

   // localization
   if(kappa > 1e-12 || calcKinLoc)  {
      computeLocalization(functions, calcGrad);
      functionalValue += kappa * locVal;
      if (calcGrad) grad += kappa * locGrad;
   }

   if (calcGrad)  {
      for(int is = 0; is < grad.getNSpecies (); is++)  {
         for(int iot = 0; iot < grad.getNOrbTypes(is); iot++)  {
            if (fixedOrbitals(is)(iot)) grad(is,iot).set(0.0);
         }
      }
   }
}

void SxQuamol::computeSpillage (
      const SxAtomicOrbitalsG &functions, 
      bool calcGrad)
{
   SX_CLOCK(Timer::computeSpillage);

   if (calcGrad) spillageGrad = 0.0 * spillageGrad;
   spillage = 0.0;
   double spaceNorm = 0.0;

   const SxPW &waves = *wavesPtr;
   const SxFermi &fermi = *fermiPtr;
   int nk            = waves.getNk ();
   int nSpin   = waves.getNSpin ();
   int nStates = waves.getNStates ();
   const SxGkBasis &gkBasis = *waves.getGkBasisPtr ();
   const SxAtomicStructure &structure = *gkBasis(0).structPtr;
   SxArray<SxQuantumNumbers> map = radials.getOrbitalMap (structure);
   ssize_t nOrbitals = map.getSize ();
   SxDiracVec<Double> radGFunc 
      = functions.getRadGBasisPtr ()->getRadGFunc ();

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      SxOrbitals muSet = expandRadialsG (functions,ik);
      SxOrbitals SMuSet = SPtr->apply(muSet);
      SxDMatC16 Smunu = muSet.adjoint () ^ SMuSet;
      SxDMatC16 invSmunu = Smunu.inverse ();
      SxOrbitals gradPW(muSet.nRows(),muSet.nCols());
      gradPW.set(0.0);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         SxOrbitals psiSet = waves(iSpin, ik);
         SxOrbitals SPsiSet = SPtr->apply(psiSet);
         SxDMatC16  C = invSmunu ^ (muSet.adjoint () ^ SPsiSet);
         SxOrbitals psiMuSet = muSet ^ C; 
         psiMuSet.setBasis(gkBasis(ik));
         psiMuSet.handle->auxData.ik = ik;
         psiMuSet.handle->auxData.iSpin = iSpin;
         SxOrbitals SPsiMuSet = SMuSet ^ C;
         SPsiMuSet.setBasis(gkBasis(ik));
         SPsiMuSet.handle->auxData.ik = ik;
         SPsiMuSet.handle->auxData.iSpin = iSpin;
         for (int iState = 0; iState < nStates; iState++)  {
            if (fabs(fermi.focc(iState,iSpin,ik)) > 1e-16)  { 
               SxWaveState psi = psiSet.colRef(iState);
               SxWaveState SPsi = SPsiSet.colRef(iState);
               SxWaveState SPsiMu = SPsiMuSet.colRef(iState);
               spaceNorm += gkBasis.weights(ik) * fermi.focc(iState,iSpin,ik) 
                  * dot(psi,SPsi).re;
               spillage += gkBasis.weights(ik) * fermi.focc(iState,iSpin,ik) 
                  * ( dot(psi,SPsi) - dot(psi,SPsiMu) ).re;

               if (calcGrad) {
                  SxWaveState dSPsiConj = (SPsi - SPsiMu).conj ();
                  for (ssize_t iOrbital = 0; iOrbital < nOrbitals; iOrbital++) {
                     SxComplex16 wgt = gkBasis.weights(ik)
                                     * fermi.focc(iState,iSpin,ik)
                                     * C(iOrbital,iState);
                     gradPW.colRef(iOrbital).plus_assign_ax (-wgt, dSPsiConj);
                  }
               }
            }
         }
      }
      if (calcGrad) {
         for (ssize_t iOrbital = 0; iOrbital < nOrbitals; iOrbital++) {
            int is = map(iOrbital).iSpecies;
            int ia = map(iOrbital).iAtom;
            int n = map(iOrbital).n;
            int l = map(iOrbital).l;
            int m = map(iOrbital).m;
            // Cubic Spline Edge conditions
            int mode = 0;
            if (l >= 2) mode = 1; 
            double vol = structure.cell.volume;
            SxDiracVec<Complex16> dMudR = sqrt(TWO_PI*TWO_PI*TWO_PI/vol)
               * SxYlm::getYlmNormFactor(l,m)
               * gkBasis(ik).getYlm(l,m)
               * gkBasis(ik).getPhaseFactors(is,ia);
            if (dRgdRGk(ik).getSize() == 0) setDRDR(ik);
            int iot = radials.getIOT(is,n,l);
            spillageGrad(is,iot) 
               += 2.0 * ((gradPW.colRef(iOrbital) * dMudR) 
                     ^ dRgdRGk(ik)(mode)).real ();
         }
      }
   }

   spaceNorm = SxLoopMPI::sum (spaceNorm);
   spillage = SxLoopMPI::sum (spillage);

   if (calcGrad) spillageGrad.sumMPI ();

   spillage /= spaceNorm;
   spillageGrad = spillageGrad / spaceNorm;
}

void SxQuamol::computeKineticEnergy (const SxAtomicOrbitalsG &functions, bool calcGrad)
{
   SX_CLOCK(Timer::computeKineticEnergy);
   const SxRadGBasis &radGBasis = *functions.getRadGBasisPtr ();
   if (calcGrad) eKinGrad = 0.0 * eKinGrad;
   eKin = 0.0;

   int nSpecies = functions.getNSpecies ();
   int nFunctions = functions.getNOrbTypes ();
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      int nOrbTypes = functions.getNOrbTypes(iSpecies);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         const SxDiracVec<Double> &g = radGBasis.getRadGFunc ();
         const SxDiracVec<Double> &dg = radGBasis.getRadDG ();
         const SxDiracVec<Double> &mu = functions(iSpecies,iot);
         SxDiracVec<Double> g2mu = g * g * mu;
         double norm = (mu * g2mu * dg).sum ();
         double eKinMu = (mu * g2mu * g * g * dg).sum () / norm;
         eKin += eKinMu / (1.0 * nFunctions);
         if (calcGrad) eKinGrad(iSpecies,iot) 
            = (2.0 * g2mu * g * g * dg / norm
            - 2.0 * eKinMu * g2mu * dg / norm) / (1.0 * nFunctions); 
      }
   }
}

void SxQuamol::computeLocalization (
      const SxAtomicOrbitalsG &functions, 
      bool calcGrad)
{
   SX_CLOCK(Timer::computeLocalization);
   locVal = 0.0;
   if (calcGrad) locGrad = 0.0 * functions;

   const SxRadGBasis &radGBasis = *functions.getRadGBasisPtr ();
   const SxDiracVec<Double> &g = radGBasis.getRadGFunc ();
   const SxDiracVec<Double> &dg = radGBasis.getRadDG ();
   SxDiracVec<Double> g2dg = g * g * dg;

   int nSpecies = radials.getNSpecies ();
   SxArray<SxArray<SxDiracVec<Double> > > muSetG (nSpecies);
   int nFunctions = functions.getNOrbTypes ();

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      int nOrbTypes = radials.getNOrbTypes(iSpecies);
      muSetG(iSpecies).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         const SxDiracVec<Double> &muG = functions(iSpecies,iot);
         int l = muG.handle->auxData.l;
         double normG = (muG * muG * g2dg).sum ();
         SxDiracVec<Double> locMuG = (locRGGPrime(l) ^ muG);
         double preFactor = 1.0 / normG / (1.0 * nFunctions);
         double locValMu = (muG * locMuG).sum ();
         locVal += preFactor * locValMu;

         if (calcGrad)  {
            muSetG(iSpecies)(iot) = 2.0 * preFactor 
               * (locMuG - locValMu/normG * muG * g2dg);
            muSetG(iSpecies)(iot).setBasis(&radGBasis);
            muSetG(iSpecies)(iot).handle->auxData.is = muG.handle->auxData.is;
            muSetG(iSpecies)(iot).handle->auxData.ia = muG.handle->auxData.ia;
            muSetG(iSpecies)(iot).handle->auxData.n  = muG.handle->auxData.n;
            muSetG(iSpecies)(iot).handle->auxData.l  = muG.handle->auxData.l;
            muSetG(iSpecies)(iot).handle->auxData.m  = muG.handle->auxData.m;

         } // calcGrad
      } // iot
   } // iSpecies

   if (calcGrad)  {
      locGrad = SxAtomicOrbitalsG (muSetG,functions.getRadGBasisPtr (),0); 
   } // calcGrad
}
/* old version: numerical gradient differs in smal region and end region
void SxQuamol::computeLocalization (
      const SxAtomicOrbitalsG &functions, 
      bool calcGrad)
{
   SX_CLOCK(Timer::computeLocalization);
   locVal = 0.0;
   if (calcGrad) locGrad = 0.0 * functions;

   SxConstPtr<SxRadGBasis> radGBasisPtr = functions.getRadGBasisPtr ();

   ssize_t nGPoints = radGBasisPtr->getRadGFunc().getSize();
   double gMax =  radGBasisPtr->getRadGFunc()(nGPoints-1);
   ssize_t nPoints = nGPoints * ssize_t(rMax / gMax);
   SxConstPtr<SxRadRBasis> radRBasisPtr 
      = SxConstPtr<SxRadRBasis>::create(0.0, rMax, nPoints);

   SxAtomicOrbitalsR radR = getOrbitals (functions, radRBasisPtr);

   int nSpecies = radials.getNSpecies ();
   
   SxAtomicOrbitalsR locGradR = 0.0 * radR;

   SxArray<SxArray<SxDiracVec<Double> > > muSetG (nSpecies);
   const SxDiracVec<Double> &r = radRBasisPtr->getRadRFunc ();
   SxDiracVec<Double> locR = r - rStart;
   for (int ir = 0; ir < locR.getSize (); ir++) 
      if (locR(ir) < 1e-12) locR(ir) = 0;
   const SxDiracVec<Double> &dr = radRBasisPtr->getRadDR ();
   SxDiracVec<Double> r2dr = r * r * dr;
   SxDiracVec<Double> locR3r2dr = locR * locR * locR * r2dr;

   int nFunctions = functions.getNOrbTypes ();

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      int nOrbTypes = radials.getNOrbTypes(iSpecies);
      muSetG(iSpecies).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         const SxDiracVec<Double> &muG = functions(iSpecies,iot);
         int l = muG.handle->auxData.l;
         const SxDiracVec<Double> &muR = radR(iSpecies,iot);
         double norm = (muR * muR * r2dr).sum ();
         double locValMu = (muR * muR * locR3r2dr).sum () / norm;
         locVal += locValMu / (1.0 * nFunctions);

         if (calcGrad)  {
            locGradR(iSpecies,iot) 
               = (1.0 / norm * muR * locR3r2dr // loc gradient
               - 1.0 * locValMu / norm * muR * r2dr) // keep Norm 
               / (1.0 * nFunctions); // average
             SxDiracMat<Double> trafoMat (muR.getSize (), muG.getSize ());
            for (int ig = 0; ig < muG.getSize (); ig++)  {
               SxDiracVec<Double> jl = SxRadGBasis::jsb(l,g(ig) * r);
               trafoMat.colRef(ig) 
                  << sqrt(2.0/PI) * jl * g(ig) * g(ig) * dg(ig);
            }
            muSetG(iSpecies)(iot) = 2.0 
               * (trafoMat.transpose () ^ locGradR(iSpecies,iot));
            muSetG(iSpecies)(iot).setBasis(radGBasisPtr.getPtr ());
            muSetG(iSpecies)(iot).handle->auxData.is = muR.handle->auxData.is;
            muSetG(iSpecies)(iot).handle->auxData.ia = muR.handle->auxData.ia;
            muSetG(iSpecies)(iot).handle->auxData.n  = muR.handle->auxData.n;
            muSetG(iSpecies)(iot).handle->auxData.l  = muR.handle->auxData.l;
            muSetG(iSpecies)(iot).handle->auxData.m  = muR.handle->auxData.m;

         } // calcGrad
      } // iot
   } // iSpecies

   if (calcGrad)  {
      locGrad = SxAtomicOrbitalsG (muSetG,radGBasisPtr,0); 
   } // calcGrad
}
*/

SxAtomicOrbitalsG SxQuamol::calcNumericGrad (const SxAtomicOrbitalsG &functions,
      double h)
{
   SX_CLOCK(Timer::calcNumericGrad);
   SxAtomicOrbitalsG trial = 1.0 * functions;
   SxAtomicOrbitalsG result = 0.0 * functions;

   int nSpecies = radials.getNSpecies();

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      int nOrbTypes = radials.getNOrbTypes(iSpecies);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         size_t dim = functions(iSpecies,iot).getSize ();
         for (size_t i = 0; i < dim; i++)  {
            trial(iSpecies,iot)(i) += h;
            computeFunctional (trial, false, false);
            double PH = functionalValue;
            trial(iSpecies,iot)(i) -= 2.0*h;
            computeFunctional (trial, false, false);
            double MH = functionalValue;
            trial(iSpecies,iot)(i) += h;
            result(iSpecies,iot)(i) = (PH-MH) / (2.0*h);
         }
      }
   }

   return result;
}

void SxQuamol::compute ()
{
   if (checkGrad)  {
      cout << "Numeric grad calculation" << endl;
      SxAtomicOrbitalsG numericGrad = calcNumericGrad (radials);
      cout << "Analytic grad calculation" << endl;
      computeFunctional (radials, true, false);
      SX_MPI_MASTER_ONLY
      
      {
         numericGrad.print("NumericGradient");
         grad.print("AnalyticGradient");
         (grad-numericGrad).print("GradientDiff");
      }
      
      return;
   }

   // Minimisation sheme CGLMin
   cout << "Start CGLMin" << endl;
   cout << SX_SEPARATOR;

   double sw = 0.1, oldsw = sw;
   SxVector3<Double> guess = SxVector3<Double> (0.0, sw, 2.0*sw);
   SxVector3<Double> yGuess = SxVector3<Double> (-1.0,-1.0,-1.0);
   int step = 0;
   double oldFuncValue = 0.0;
   double delta = functionalValue - oldFuncValue;
   int problem = 0;
   printLine = false;
   bool doNotStop = false;
   bool converged = false;
   double lastError = functionalValue;

   computeFunctional (radials, true, true);
   // Searching for Minimum
   grad = -1.0 * grad;
   //grad = setOrthogonal(grad, radials);
   SxAtomicOrbitalsG dir = 1.0 * grad;
   SxAtomicOrbitalsG oldGrad = 1.0 * grad;

   while (((step <= maxSteps) && !converged) || doNotStop) {
      SX_CLOCK(Timer::stepClock);
      step++;
      if (step == printStep) printLine = true;
      doNotStop = false;

      oldFuncValue = functionalValue;
      
      // Conjugate Gradient;
      if (step > 1)  {
         // Calculate new direction
         // Fletcher-Reeves
         //double numerator = grad.sum(grad);
         //double denominator = oldGrad.sum(oldGrad);
         // Polak-Ribiere
         double numerator = (grad-oldGrad).sum(grad);
         double denominator = oldGrad.sum(oldGrad);
         // Hestenes-Stiefel
         //double numerator = (grad-oldGrad).sum(grad);
         //double denominator = (grad-oldGrad).sum(dir);
         double gamma = 0;
         if (fabs(denominator) > 1e-12) gamma = numerator / denominator;

         //if (gamma > 0.9) problem++; // gradient is not decreasing
         if (gamma < 0.0) problem++; // gradient has significant component 
                                     // in old direction
         // restart CG every 25 steps
         if ((problem >= 9) || (step%25 == 0))  {
            //restart CG
            gamma = 0.;
            problem = 0;
            sw = 0.1;
            cout << "RESTART CG!" << endl;
            doNotStop = true;
            if((adaptiveZeta) &&(step%25 == 0)) {
               zeta = spillage / eKin;
               cout << "Adaptive Zeta is now " << zeta << endl;
            }
            if((adaptiveKappa)&& (step%25 == 0))   {
               kappa = spillage / locVal;
               cout << "Adaptive Kappa is now " << kappa << endl;
            }
         }
         //SD
         //dir = 1.0 * grad;
         //CG
         cout << "gamma = " << gamma << endl;
         dir = grad + gamma * dir;
         //dir = setOrthogonal(dir, radials);
      }
      double dirLength = dir.sum(dir);
      SxAtomicOrbitalsG dirNorm = 1.0 * dir;
      if (dirLength > 1e-12) dirNorm = dirNorm / sqrt(dirLength);
      oldGrad = 1.0 * grad;

      // --- Lineminimization
      double gradLength = grad.sum(dirNorm); 
      if (gradLength < 0) dirNorm = -1.0 * dirNorm;
      //if (sw > 1.0) guess = SxVector3<Double> (0.0,0.1,0.2);
      //else guess = SxVector3<Double> (0.0,sw,2.0*sw);
      
      guess = SxVector3<Double> (0.0,sw,2.0*sw);
      yGuess = SxVector3<Double> (-1.0,-1.0,-1.0);

      oldsw = sw;
      sw = lineMin (dirNorm, guess, yGuess);
      int iterSteps = 0;

      updateGuess(guess,yGuess,sw);
      computeFunctional(cos(sw) * radials + sin(sw) * dirNorm, true, false);
      grad = -1.0 * grad;
      cout << "Fallbackstep: " << iterSteps 
              << ", sw = " << sw
              << ", Functional = " << functionalValue 
              << ", gradiant along dir = " 
              << fabs(grad.sum(dirNorm)) / gradLength * 100
              << endl;
      // Stabilize Linimin for big steps
      bool funcConverged = false;
      while (fabs(sw) > 1e-12 
              && !funcConverged
              && fabs(grad.sum(dirNorm)/gradLength * 100) > relLineMin 
              && iterSteps <= 20)  {
         iterSteps++;
         sw = lineMin (dirNorm, guess, yGuess);
         updateGuess(guess,yGuess,sw);
         double funcOld = functionalValue;
         computeFunctional(cos(sw) * radials + sin(sw) * dirNorm, true, false);
         grad = -1.0 * grad;
         cout << "Fallbackstep: " << iterSteps 
              << ", sw = " << sw
              << ", Functional = " << functionalValue 
              << ", gradiant along dir = " 
              << fabs(grad.sum(dirNorm)) / gradLength * 100
              << endl;
         if (fabs(functionalValue - funcOld) < 0.01 * dF) {
            cout << "Functional converged: is "
                 << "| " << functionalValue << " - "
                 << funcOld << " | = "
                 << fabs(functionalValue - funcOld) 
                 << endl;
            funcConverged = true;
         }
      }

      if (iterSteps > 20 || fabs(sw) < 1e-12) problem = 10;

      // --- Optimize
      radials = cos(sw) * radials + sin(sw) * dirNorm;

      computeFunctional(radials, true, false);
      // Search for minimum
      grad = -1.0 * grad;
      cout << "Gradient lenght in search direction reduced to " 
           << fabs(grad.sum(dirNorm))/gradLength * 100 << " %" << endl; 
      delta = functionalValue - oldFuncValue;

      cout << "Step: " << step
           << ", stepwidth: " << sw
           << ", fallback: " << iterSteps
           << ", Functional: " << functionalValue
           << ", Delta: " << delta
           << ", Spillage: " << spillage;
      cout << ", eKin: " << zeta * eKin;
      cout << ", loc : " << kappa * locVal;
      cout << endl;
      cout << "<g|g> = " << grad.sum(grad) << endl;
      
      sw = 0.5*(sw+oldsw);
   
      if (((lastError < dF) && (fabs(delta) < dF)) || grad.sum(grad) < dRes) 
         converged = true;
   
      lastError = fabs(delta);
   }
}

void SxQuamol::updateGuess (SxVector3<Double> &x, SxVector3<Double> &y, double val)
{
   SxVector3<Double> xCopy = 1.0 * x;
   SxVector3<Double> yCopy = 1.0 * y;

   if (   fabs(xCopy(0) - val) < 1e-12
       || fabs(xCopy(1) - val) < 1e-12
       || fabs(xCopy(2) - val) < 1e-12)
      return;
   
   if (val < xCopy(0))  {
      x(0) = val;
      x(1) = xCopy(0);
      x(2) = xCopy(1);
      y(0) = -1.0;
      y(1) = yCopy(0);
      y(2) = yCopy(1);
      return;
   }

   if (val > xCopy(2))  {
      x(0) = xCopy(1);
      x(1) = xCopy(2);
      x(2) = val;
      y(0) = yCopy(1);
      y(1) = yCopy(2);
      y(2) = -1.0;
      return;
   }

   if (val < x(1))  {
      x(0) = xCopy(0);
      x(1) = val;
      x(2) = xCopy(1);
      y(0) = yCopy(0);
      y(1) = -1.0;
      y(2) = yCopy(1);
      return;
   }

   if (val > x(1))  {
      x(0) = xCopy(1);
      x(1) = val;
      x(2) = xCopy(2);
      y(0) = yCopy(1);
      y(1) = -1.0;
      y(2) = yCopy(2);
      return;
   }

   SX_EXIT;
}

double SxQuamol::lineMin (const SxAtomicOrbitalsG &dir,
                          const SxVector3<Double> &x,
                          SxVector3<Double> &y)
{
   SX_CHECK(fabs(x(1) - x(0)) > 1e-6);
   SX_CHECK(fabs(x(2) - x(0)) > 1e-6);
   SX_CHECK(fabs(x(2) - x(1)) > 1e-6);

   SX_CLOCK (Timer::lineMin);
   // Setup Guess functions
   for (int i = 0; i < 3; i++)  {
      if (y(i) < 0)  {
         SxAtomicOrbitalsG trial = cos(x(i)) * radials + sin(x(i)) * dir;
         computeFunctional(trial, false, false);
         y(i) = functionalValue;
      }
   }
   
   // parabel Fit
   SxComplex<double> pFit = parabelFit(x(0), x(1), x(2), y(0), y(1), y(2));
   // for debug
   if (printLine)  {
      cout << SX_SEPARATOR;
      cout << "Result of Linemin: " << endl;
      cout << pFit << endl;
      cout << SX_SEPARATOR;

      int dim = 51;
      SxVector<Double> nX(dim);
      SxVector<Double> nY(dim);
      SxVector<Double> tanY(dim);
      double delta = 0.5 * PI / dim;
      // I dont unterstand why, but I do not need g^2 in "integration"
      double proj = grad.sum(dir);
      for (int i = 0; i < dim; i++)  {
         nX(i) = i * delta;
         SxAtomicOrbitalsG trial = cos(nX(i)) * radials + sin(nX(i)) * dir;
         computeFunctional(trial, false, false);
         nY(i) = functionalValue;
         tanY(i) = y(0) - proj * sin(nX(i));
      }

      SX_MPI_MASTER_ONLY
      {
         SxString file = "NumericLine.dat";
         SxBinIO out;
         out.open(file,SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(nX,nY);
         out.close();
         file = "Tangent.dat";
         out.open(file,SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(nX,tanY);
         out.close();
      }
      SX_QUIT;
   }

   return pFit.re;
}

SxComplex<double> SxQuamol::parabelFit (double x0, double x1, double x2,
                                        double N0, double N1, double N2)
{
   // Points are sorted
   // x0 < x1 < x2
   SxComplex<double> result;

   double z1 = x1 - x0;
   double z2 = x2 - x0;

   SX_CHECK(fabs(z1) > 1e-6);
   SX_CHECK(fabs(z2) > 1e-6);
   SX_CHECK(fabs(z1-z2) > 1e-6);

   double a = (z1 * (N0-N2) + z2 * (N1-N0))
         / (z1*z2*(z1-z2));
   double b = -(z1 * z1 * (N0-N2) + z2 * z2 * (N1-N0))
         / (z1*z2*(z1-z2));

   double opt = 0.0;

   if (fabs(a) > 1e-6 * fabs(b) && a > 0)  {
      opt = -0.5 * b / a;
      // if Minimum left of initial point, decrase interval
      if (opt + x0 < 0) opt = 0.5 * z1;
   } else if (fabs(a) <= 1e-10)  {
      cout << "Parabolic fit: linear function ";
      // if increasing, shorten interval
      // if decreasing, enlarge interval
      if (N0 > N2) opt = 1.1 * z2;
      if (N0 < N2) opt = 0.5 * z1;
   } else {
      cout << "Parabolic fit: Negative a = " << a << endl;
      // Parabel has maximum ? 
      // if increasing side, shorten interval
      // if decreasing side, enlarge interval
      if ((N0 > N1) && (N1 > N2)) opt = 1.1 * z2;
      else if ((N0 < N1) && (N1 < N2)) opt = 0.5 * z1;
      else  {
         cout << "Unexpected shape in LineMin!" << endl;
         printLine = true;
      }
   }
   result.re = opt + x0;
   result.im = a * opt * opt + b * opt + N0;
   
   if (print || printLine)  {
      cout << SX_SEPARATOR;
      cout << "ParabelFit: " << endl;
      cout << "x0: " << x0 << ", N0: " << N0 << endl;
      cout << "x1: " << x1 << ", N1: " << N1 << endl;
      cout << "x2: " << x2 << ", N2: " << N2 << endl;
      cout << "a = " << a << endl;
      cout << "b = " << b << endl;
      cout << "sw = " << result.re << endl;
      cout << SX_SEPARATOR;
   }

   // for debug
   if (printLine) {
      int dim = 51;
      SxVector<Double> x (dim);
      SxVector<Double> y (dim);
      double delta = 0.5 * PI / dim;
      for (int i = 0; i < dim; i++)  {
         x(i) = i * delta + x0;
         y(i) = a * (x(i)-x0)*(x(i)-x0) + b * (x(i)-x0) + N0;
      }

      SX_MPI_MASTER_ONLY
      {
         SxString file = "FitLine.dat";
         SxBinIO out;
         out.open(file,SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(x,y);
         out.close();
      }
   }

   return result;
}

SxOrbitals SxQuamol::expandRadialsG (const SxAtomicOrbitalsG &functions, int ik)
{
   SX_CLOCK(Timer::radG2mu);

   const SxGkBasis &gkBasis = *wavesPtr->getGkBasisPtr ();
   const SxRadGBasis &radGBasis = *functions.getRadGBasisPtr ();
   SxArray<SxQuantumNumbers> map 
      = functions.getOrbitalMap(*gkBasis(ik).structPtr);

   //Expand the radialfunctions
   int nOrbitals = radials.getNOrbitals(*gkBasis(ik).structPtr);
   SxDiracVec<Double> gVec = sqrt(gkBasis(ik).g2);
   size_t dim = gkBasis(ik).g2.getSize ();
   SxOrbitals result(dim, nOrbitals);
   SxArray<SxArray<SxDiracVec<Complex16> > > funcGk (functions.getNSpecies ());
   for (int is = 0; is < functions.getNSpecies (); is++)  {
      int nOrbTypes = functions.getNOrbTypes(is);
      funcGk(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         SxDiracVec<Double> splineVec = radGBasis.toSpline(functions(is,iot));
         funcGk(is)(iot) = (gkBasis(ik) | splineVec);
      }
   }

   for (int io = 0; io < nOrbitals; io++)   {
      int is  = map(io).iSpecies;
      int ia  = map(io).iAtom;
      int iot = map(io).n;
      int l   = map(io).l;
      int m   = map(io).m;
      result.colRef(io) <<= funcGk(is)(iot);
      result.colRef(io) *= SxYlm::getYlmNormFactor(l,m) * gkBasis(ik).getYlm(l,m);
      result.colRef(io) *= gkBasis(ik).getPhaseFactors(is,ia);
   }

   result.setBasis(gkBasis(ik));
   result.handle->auxData.ik = ik;

   return result;
}

void SxQuamol::checkTrafo (const SxAtomicOrbitalsG &functionsIn)
{
   SX_CLOCK(Timer::checkTrafo);
   SxArray<SxArray<SxDiracVec<Double> > > muSet;
   cout << SX_SEPARATOR;
   cout << "Transformationtest" << endl;
   cout << SX_SEPARATOR;
   SxAtomicOrbitalsG functions = 1.0 * functionsIn;
   int nk = wavesPtr->getNk ();
   int nSpecies = functions.getNSpecies ();
   SxConstPtr<SxRadGBasis> radGBasisPtr = functions.getRadGBasisPtr ();
   const SxRadGBasis &radGBasis = *radGBasisPtr;
   const SxGkBasis &GkBasis = *wavesPtr->getGkBasisPtr ();
   const SxAtomicStructure &structure = *GkBasis(0).structPtr;
   SxArray<SxQuantumNumbers> map = radials.getOrbitalMap (structure);

   functions.normalize ();

   SX_MPI_MASTER_ONLY {functions.print("InitRadG");}
   
   SxArray<SxArray<SxDiracVec<Double> > > xData (nSpecies), yData(nSpecies);
   SxVector<Int> kDims(nk);
   for (int ik = 0; ik < nk; ik++)   {
      kDims(ik) = (int)GkBasis(ik).g2.getSize ();
   }

   muSet.resize(nSpecies);
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++ )  {
      int nOrbTypes = functions.getNOrbTypes(iSpecies);
      muSet(iSpecies).resize(nOrbTypes);
      xData(iSpecies).resize(nOrbTypes);
      yData(iSpecies).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         int l = functions(iSpecies,iot).handle->auxData.l;
         int nAtoms = structure.getNAtoms(iSpecies);
         xData(iSpecies)(iot).resize(kDims.sum() * nAtoms * (2*l+1));
         yData(iSpecies)(iot).resize(kDims.sum() * nAtoms * (2*l+1));
         xData(iSpecies)(iot).set(0.0);
         yData(iSpecies)(iot).set(0.0);

         SxDiracVec<Double> &function = functions(iSpecies,iot);
         double norm = radGBasis.integrate(function * function);
         cout << "iSpecies = " << iSpecies
              << ", iot = " << iot
              << ", norm = " << norm << endl;
      }
   }

   cout << SX_SEPARATOR;
   cout << "G+k space" << endl;
   cout << SX_SEPARATOR;

   for (int ik = 0; ik < nk; ik++)   {
      SxOrbitals orbitals = expandRadialsG (functions,ik);
      for (int iOrbital = 0; iOrbital < map.getSize(); iOrbital++)   {
         int is = map(iOrbital).iSpecies;
         int ia = map(iOrbital).iAtom;
         int n = map(iOrbital).n;
         int l = map(iOrbital).l;
         int m = map(iOrbital).m;
         int iot = radials.getIOT(is,n,l);
         int nAtoms = structure.getNAtoms(is);
         SxDiracVec<Double> data = (orbitals.colRef(iOrbital)
               * GkBasis(ik).getPhaseFactors(is,ia).conj() //shift
               * GkBasis(ik).getYlm(l,m) // project l
               * SxYlm::getYlmNormFactor(l,m)
               * sqrt(2.0 * structure.cell.volume/PI)).real ();
         int dim = kDims(ik);
         int kOffset = 0;
         if(ik > 0) kOffset = kDims(SxIdx(0,ik-1)).sum() * nAtoms * (2*l+1);
         int atomOffset = ia * (2*l+1) * dim;
         int mOffset = (m+l) * dim;
         int low = kOffset + atomOffset + mOffset;
         int high = low + dim - 1;
         SxIdx idx (low,high);
         xData(is)(iot)(idx) <<= sqrt(GkBasis(ik).g2);
         yData(is)(iot)(idx) <<= data;
      } // iOrbital
   } // ik

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++ )  {
      int nOrbTypes = functions.getNOrbTypes(iSpecies);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         SxDiracVec<Double> xVec = xData(iSpecies)(iot);
         SxDiracVec<Double> yVec = yData(iSpecies)(iot);
         int l = functions(iSpecies,iot).handle->auxData.l;
         yVec.handle->auxData.l = l;
         SxCubicSpline<SxDiracVec<Double> > fit;
         if (l == 0) // ?(May have CUSP, no Mirror)?
            fit  = SxCubicSpline<SxDiracVec<Double> > (
                  xVec,
                  yVec,
                  radGBasis.getRadGFunc(),
                  SxCubicSpline<SxDiracVec<Double> >::Hermite,
                  SxCubicSpline<SxDiracVec<Double> >::MirrorPlane);
         else if (l & 1) // (l odd)
            fit  = SxCubicSpline<SxDiracVec<Double> > (
                  xVec,
                  yVec,
                  radGBasis.getRadGFunc(),
                  SxCubicSpline<SxDiracVec<Double> >::Hermite,
                  SxCubicSpline<SxDiracVec<Double> >::MirrorPoint);
         else // (l even and not zero)
            fit  = SxCubicSpline<SxDiracVec<Double> > (
                  xVec,
                  yVec,
                  radGBasis.getRadGFunc(),
                  SxCubicSpline<SxDiracVec<Double> >::Hermite,
                  SxCubicSpline<SxDiracVec<Double> >::MirrorPlane);
         muSet(iSpecies)(iot) = fit.getYFit ();
         SxDiracVec<Double> &function = muSet(iSpecies)(iot);
         function.setBasis(&radGBasis);
         function.handle->auxData.is = iSpecies;
         function.handle->auxData.ia = -1;
         function.handle->auxData.n = iot;
         function.handle->auxData.l = l;
         function.handle->auxData.m = NONE_M;
         double norm = radGBasis.integrate(function * function);
         cout << "iSpecies = " << iSpecies
              << ", iot = " << iot
              << ", norm = " << norm << endl;
      } // iot
   } // iSpecies

   SxAtomicOrbitalsG GkTrafo(muSet,radGBasisPtr,0);
   
   SX_MPI_MASTER_ONLY {GkTrafo.print("GkTrafo");}

   cout << SX_SEPARATOR;
   cout << "RadialR space" << endl;
   cout << SX_SEPARATOR;

   int nGPoints = (int)radGBasisPtr->getRadGFunc().getSize();
   double gMax =  radGBasisPtr->getRadGFunc()(nGPoints-1);
   int nPoints = int (nGPoints * rMax / gMax);

   SxConstPtr<SxRadRBasis> radRBasisPtr 
      = SxConstPtr<SxRadRBasis>::create(0.0, rMax, nPoints);
   SxAtomicOrbitalsR radR = getOrbitals (functions,radRBasisPtr);

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++ )  {
      int nOrbTypes = functions.getNOrbTypes(iSpecies);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         SxDiracVec<Double> &function = muSet(iSpecies)(iot);
         function = (radGBasis | radR(iSpecies,iot));
         double norm = radGBasis.integrate(function * function);
         cout << "iSpecies = " << iSpecies
              << ", iot = " << iot
              << ", norm = " << norm << endl;
      } //iot
   } //iSpecies

   SxAtomicOrbitalsG RadRTrafo(muSet,radGBasisPtr,0);
   
   SX_MPI_MASTER_ONLY{RadRTrafo.print("RadRTrafo");}
   
   cout << SX_SEPARATOR;
}

SxAtomicOrbitals SxQuamol::getOrbitals (
      const SxAtomicOrbitalsG &functions, 
      SxConstPtr<SxRadBasis> radBasisPtr)
{
   SX_CLOCK(Timer::getROrbitals);
   const SxRadBasis &radBasis = *radBasisPtr;
   int nSpecies = functions.getNSpecies ();
   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = functions.getNOrbTypes(is);
      muSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++) {
        muSet(is)(iot) = (radBasis | functions(is,iot));
      }
   }
   
   SxAtomicOrbitals result (muSet,radBasisPtr);
   
   return result;
}

SxAtomicOrbitalsR SxQuamol::getOrbitals (
      const SxAtomicOrbitalsG &functions, 
      SxConstPtr<SxRadRBasis> radRBasisPtr)
{
   SX_CLOCK(Timer::getROrbitals);
   const SxRadRBasis &radRBasis = *radRBasisPtr;
   int nSpecies = functions.getNSpecies ();
   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = functions.getNOrbTypes(is);
      muSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++) {
        muSet(is)(iot) = (radRBasis | functions(is,iot));
        /*
        int l = functions(is,iot).handle->auxData.l;
        muSet(is)(iot) = (radGToRadR(l) ^ functions(is,iot));
        muSet(is)(iot).setBasis(radRBasisPtr.getPtr ());
        muSet(is)(iot).handle->auxData.is = functions(is,iot).handle->auxData.is;
        muSet(is)(iot).handle->auxData.ia = functions(is,iot).handle->auxData.ia;
        muSet(is)(iot).handle->auxData.n  = functions(is,iot).handle->auxData.n;
        muSet(is)(iot).handle->auxData.l  = functions(is,iot).handle->auxData.l;
        muSet(is)(iot).handle->auxData.m  = functions(is,iot).handle->auxData.m;
        */
      }
   }
   
   SxAtomicOrbitalsR result (muSet,radRBasisPtr,0);
   
   return result;
}

SxAtomicOrbitalsG SxQuamol::setOrthogonal (
      const SxAtomicOrbitalsG &functions, 
      const SxAtomicOrbitalsG &referenz)
{
   SxAtomicOrbitalsG result = 1.0 * functions;
   SxDiracVec<Double> dg = referenz.getRadGBasisPtr()->getRadDG();
   SxDiracVec<Double> g = referenz.getRadGBasisPtr()->getRadGFunc();
   for (int iSpecies = 0; iSpecies < result.getNSpecies (); iSpecies++)  {
      for (int iot = 0; iot < result.getNOrbTypes(iSpecies); iot++)  {
         SxDiracVec<Double> &res = result(iSpecies,iot);
         const SxDiracVec<Double> &func = functions(iSpecies,iot);
         for (int jot = 0; jot < referenz.getNOrbTypes(iSpecies); jot++)  {
            const SxDiracVec<Double> &ref = referenz(iSpecies,jot);
            if (res.handle->auxData.l == ref.handle->auxData.l)  {
               double normSqr = (ref * ref * dg).sum ();
               double proj = (ref * func * dg).sum ();
               res -= proj / normSqr * ref;
            }
         }
      }
   }

  return result; 
}

SxAtomicOrbitalsG SxQuamol::setOrthogonal (
      const SxAtomicOrbitalsG &functions)
{
   SxAtomicOrbitalsG result = 1.0 * functions;
   SxDiracVec<Double> dg = functions.getRadGBasisPtr()->getRadDG();
   SxDiracVec<Double> g = functions.getRadGBasisPtr()->getRadGFunc();
   for (int iSpecies = 0; iSpecies < result.getNSpecies (); iSpecies++)  {
      for (int iot = 0; iot < result.getNOrbTypes(iSpecies); iot++)  {
         const SxDiracVec<Double> &func = functions(iSpecies,iot);
         SxDiracVec<Double> &work = result(iSpecies,iot);
         for (int jot = 0; jot < iot; jot++)  {
            const SxDiracVec<Double> &res = result(iSpecies,jot);
            if (res.handle->auxData.l == work.handle->auxData.l)  {
               double normSqr = (res * res * g * g * dg).sum ();
               double proj = (res * func * g * g * dg).sum ();
               work -= proj / normSqr * res;
            }
         }
      }
   }

   result.normalize ();

   return result; 
}

void SxQuamol::completenessProfile (const SxAtomicOrbitalsG &functions)
{
   SxConstPtr<SxRadGBasis> radGBasisPtr = functions.getRadGBasisPtr ();
   int nGPoints = (int)radGBasisPtr->getRadGFunc().getSize();
   double gMax =  radGBasisPtr->getRadGFunc()(nGPoints-1);
   int nPoints = int (nGPoints * rMax / gMax);
   SxConstPtr<SxRadRBasis> radRBasisPtr 
      = SxConstPtr<SxRadRBasis>::create(0.0, rMax, nPoints);
   SxAtomicOrbitalsR radR = getOrbitals (functions,radRBasisPtr);
   // define exponents for test gaussians
   int xDim = 201;
   SxVector<Double> logX (xDim);
   double start = -10.0;
   double end = 10.0;
   double delta = (end - start) / (xDim - 1);
   for (int i = 0; i < xDim; i++)  {
      logX(i) = start + i*delta;
   }
   for (int is = 0; is < radials.getNSpecies(); is++)  {
      const SxDiracVec<Double> &r = radRBasisPtr->getRadRFunc ();
      int nOrbTypes = radR.getNOrbTypes(is);
      SxMatrix<Double> Smunu(nOrbTypes,nOrbTypes);
      Smunu.set(0.0);
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         int l = radR(is,iot).handle->auxData.l;
         for (int jot = 0; jot < nOrbTypes; jot++)  {
            if (radR(is,jot).handle->auxData.l == l)  {
               Smunu(iot,jot) += radR.dot (radR,is, iot, is, jot);
            }
         }
      }
      SxMatrix<Double> invSmunu = Smunu.inverse();

      for (int l = 0; l <= radR.getLMax(is); l++)  {
         SxVector<Double> result (xDim);
         result.set(0.0);
         for (int i = 0; i < xDim; i++)  {
            SxDiracVec<Double> probeGauss 
               = pow(r, double(l)) * exp(-r.sqr () * exp(logX(i)));
            probeGauss.setBasis(&*radRBasisPtr);
            double probeGaussNorm 
               = sqrt(radRBasisPtr->integrate(probeGauss * probeGauss));
            probeGauss /= probeGaussNorm;
            probeGauss.handle->auxData.is = is;
            probeGauss.handle->auxData.ia = -1;
            probeGauss.handle->auxData.n = -1;
            probeGauss.handle->auxData.l = l;
            probeGauss.handle->auxData.m = NONE_M;
            for (int iot = 0; iot < nOrbTypes; iot++)  {
               if (radR(is,iot).handle->auxData.l == l)  {
                  for (int jot = 0; jot < nOrbTypes; jot++)  {
                     if (radR(is,jot).handle->auxData.l == l)  {
                        result(i) += radRBasisPtr->integrate(probeGauss * radR(is,iot))
                                   * invSmunu(iot,jot)
                                   * radRBasisPtr->integrate(probeGauss * radR(is,jot));
                     }
                  }
               }
            }
         }

         SX_MPI_MASTER_ONLY
         {
            SxString file = "C-Profile-" + SxString(is) +"-" + SxString(l) +".dat";
            SxBinIO out;
            out.open(file,SxBinIO::ASCII_WRITE_ONLY);
            out.writeXYPlot(logX, result);
            out.close();
         }
      }
   }
}

void SxQuamol::getNormContribution (const SxAtomicOrbitalsG &functions,
      SxConstPtr<SxRadBasis> radBasisPtr)
{

   const SxPW &waves = *wavesPtr;
   const SxFermi &fermi = *fermiPtr;
   const SxGkBasis &gkBasis = *waves.getGkBasisPtr ();

   int nSpecies = functions.getNSpecies ();
   int nk = waves.getNk();
   int nSpin = waves.getNSpin();
   const SxVector<Double> &kWeight = gkBasis.weights;

   SxArray<SxArray<SxMatrix<Complex16> > > P (nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = radials.getLMax (is);
      P(is).resize(lMax+1);
      for(int l = 0; l <= lMax; l++)  {
         int nFL = functions.getFuncPerL (is,l);
         if (nFL == 0) continue;
          P(is)(l).reformat(nFL,nFL);
          P(is)(l).set(0.0);
      }
   }

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      SxOrbitals muSet = expandRadialsG (functions,ik);
      SxOrbitals SMuSet = SPtr->apply(muSet);
      SxDMatC16 Smunu = muSet.adjoint () ^ SMuSet;
      SxDMatC16 invSmunu = Smunu.inverse ();
      SxArray<SxQuantumNumbers> map 
         = functions.getOrbitalMap(*gkBasis(ik).structPtr);
      for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
         SxOrbitals psiSet = waves(iSpin, ik);
         int nStates = waves.getNStates ();
         SxIdx idx (0,nStates - 1);
         SxDiracVec<Double> focc = fermi.focc(iSpin,ik)(idx);
         // C(nOrbitals,nStates)
         SxDMatC16  C = invSmunu ^ (SMuSet.adjoint () ^ psiSet);
         for (int is = 0; is < nSpecies; is++)  {
            int nAtoms = gkBasis(ik).structPtr->getNAtoms(is);
            double atomFactor = 1.0 / nAtoms;
            int lMax = radials.getLMax (is);
            for(int l = 0; l <= lMax; l++)  {
               double lFactor = 1.0 / (2.*l+1.0);
               int nFL = functions.getFuncPerL (is,l);
               if (nFL == 0) continue;
               for (int ifl = 0; ifl < nFL; ifl++)  {
                  int iot = functions.funcLMap(is)(l)(ifl);
                  for (int jfl = ifl; jfl < nFL; jfl++)  {
                     int jot = functions.funcLMap(is)(l)(jfl);
                     for (int ia = 0; ia < nAtoms; ia++)  {
                        for (int m = -l; m <= l ; m++)  {
                           int iOrbital = 
                              functions.getOrbitalIdx(is,ia,iot,l,m,map);
                           int jOrbital = 
                              functions.getOrbitalIdx(is,ia,jot,l,m,map);

                           SxComplex16 cDot 
                              = ( focc 
                                * C.row(iOrbital).conj()
                                * C.row(jOrbital)).sum();
                           P(is)(l)(ifl,jfl) += 
                              atomFactor * lFactor * kWeight(ik) * cDot;
                           P(is)(l)(jfl,ifl) += 
                              atomFactor * lFactor * kWeight(ik) * cDot.conj ();
                        }
                     }
                  }
               }
            }
         }
      }
   }

   for (int is = 0; is < nSpecies; is++)  {
      int lMax = radials.getLMax (is);
      for(int l = 0; l <= lMax; l++)  {
         SxLoopMPI::sum (P(is)(l));
      }
   }

   
   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = functions.getLMax (is);
      for(int l = 0; l <= lMax; l++)  {
         int nFL = functions.getFuncPerL (is,l);
         if (nFL == 0) continue;
         SxMatrix<Complex16>::Eigensystem eig;
         eig = P(is)(l).eigensystem ();
         cout << SX_SEPARATOR;
         cout << "is = " << is
              << ", l = " << l
              << ", functions = " << nFL
              << endl;
         eig.vals.real ().print ();
         cout << "Total : " << eig.vals.sum().re << endl;

         muSet(is).resize(eig.vals.getSize());
         for (int iot = 0; iot < eig.vals.getSize(); iot++)  {
            int dim = int(radials.getFuncL(is,l,0).getSize());
            SxDiracVec<Double> newRadial (dim);
            newRadial.set(0.0);
            for (int ifl = 0; ifl < nFL; ifl++)  {
               newRadial += eig.vecs(ifl,iot) * functions.getFuncL(is,l,ifl);
            }
            newRadial.handle->auxData.is = is;
            newRadial.handle->auxData.ia = -1;
            newRadial.handle->auxData.n  = iot;
            newRadial.handle->auxData.l  = l;
            newRadial.handle->auxData.m  = NONE_M;
            newRadial.setBasis(functions.getRadGBasisPtr ().getPtr ());
            muSet(is)(iot) = newRadial.getCopy();
         }
      }
   }
   cout << SX_SEPARATOR;

   SxAtomicOrbitalsG resultG = 
      SxAtomicOrbitalsG(muSet,functions.getRadGBasisPtr (),0);
   SxAtomicOrbitals result = getOrbitals(resultG,radBasisPtr);
   result.normalize ();

   SX_MPI_MASTER_ONLY{result.write("normContrib.sxb");}
}

void SxQuamol::printSpillagePerState (const SxAtomicOrbitalsG &functions)
{
   SX_CHECK(wavesPtr->getGkBasisPtr ().getPtr () != NULL);

   const SxGkBasis &gkBasis = *wavesPtr->getGkBasisPtr ();
   const SxPW &waves = *wavesPtr;
   const SxFermi &fermi = *fermiPtr;


   int nk      = waves.getNk ();
   int nSpin   = waves.getNSpin ();
   int nStates = waves.getNStates ();

   double spaceNorm = 0.0;
   SxVector<Double> spillagePerState (nStates);
   spillagePerState.set(0.0);

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      SxOrbitals muSet = expandRadialsG (functions,ik);
      SxOrbitals SMuSet = SPtr->apply(muSet);
      SxDMatC16 Smunu = muSet.adjoint () ^ SMuSet;
      SxDMatC16 invSmunu = Smunu.inverse ();
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         SxDiracMat<Complex16> SPsiSet = SPtr->apply(waves(iSpin, ik));
         SxDMatC16  C = invSmunu ^ (muSet.adjoint () ^ SPsiSet);
         SxOrbitals psiMuSet = muSet ^ C;
         psiMuSet.setBasis(gkBasis(ik));
         psiMuSet.handle->auxData.ik = ik;
         psiMuSet.handle->auxData.iSpin = iSpin;
         SxOrbitals SPsiMuSet = SMuSet ^ C;
         SPsiMuSet.setBasis(gkBasis(ik));
         SPsiMuSet.handle->auxData.ik = ik;
         SPsiMuSet.handle->auxData.iSpin = iSpin;
         for (int iState = 0; iState < nStates; iState++)  {
            SxWaveState psi = waves(iState, iSpin, ik);
            SxWaveState SPsi = SPsiSet.colRef(iState);
            SxWaveState SPsiMu = SPsiMuSet.colRef(iState);
            SxComplex16 norm = dot(psi, SPsi);
            SxComplex16 projNorm = dot(psi, SPsiMu);
            spaceNorm += gkBasis.weights(ik) * fermi.focc(iState,iSpin,ik) 
               * norm.re;

            spillagePerState(iState) 
               += gkBasis.weights(ik) * fermi.focc(iState,iSpin,ik)
               * (norm - projNorm).re;
         }
      }
   }

   spaceNorm = SxLoopMPI::sum(spaceNorm);
   SxLoopMPI::sum (spillagePerState);
   spillagePerState /= spaceNorm;

   for (int iState = 0; iState < nStates; iState++)  {
      cout << "State: " << iState 
           << ", Spillage: " << spillagePerState(iState) << endl;
   }
}

void SxQuamol::calcResidues (const SxAtomicOrbitalsG &functions)
{
   SX_CHECK(wavesPtr->getGkBasisPtr ().getPtr() != NULL);
   int dim = 100;

   const SxPW &waves = *wavesPtr;
   const SxFermi &fermi = *fermiPtr;
   const SxGkBasis &gkBasis = *waves.getGkBasisPtr ();

   int nk = waves.getNk();
   int nStates = waves.getNStates();
   int nSpin = waves.getNSpin();
   SxVector<Double> kWeight = gkBasis.weights;
   
   // Assuming ik = 0 is reference structure
   SxRBasis rBasis (dim,dim,dim,gkBasis(0).structPtr->cell); 
   SxDiracVec<Complex16> Residuum(rBasis); //(nr)
   SxDiracVec<Complex16> signedResiduum(rBasis);
   Residuum.set(0.0);
   signedResiduum.set(0.0);

   for(int ik = 0; ik < nk; ik++)   {
      if (!SxLoopMPI::myWork(ik)) continue;
      SxOrbitals muSet = expandRadialsG(functions,ik);
      SxOrbitals SMuSet = SPtr->apply(muSet);
      SxDMatC16 Smunu = muSet.adjoint () ^ SMuSet;
      SxDMatC16 invSmunu = Smunu.inverse ();
      for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
         const SxDiracMat<Complex16> &Psi = waves(iSpin,ik);
         SxDMatC16  C = invSmunu ^ (SMuSet.adjoint () ^ Psi);
         SxDiracMat<Complex16> PsiMu = muSet ^ C;
         PsiMu.setBasis(gkBasis(ik));
         PsiMu.handle->auxData.ik = ik;
         PsiMu.handle->auxData.iSpin = iSpin;
         for(int iState = 0; iState < nStates; iState++)   {
            if (fabs((fermi.focc(iState, iSpin, ik))) > 1e-12)  { 
               SxDiracVec<Complex16> resDen 
                  = (PsiMu.colRef(iState) - Psi.colRef(iState)).absSqr ();
               Residuum += kWeight(ik) * fermi.focc(iState, iSpin, ik) 
                  * (rBasis|resDen); 
               SxDiracVec<Complex16> PWDen = Psi.colRef(iState).absSqr();
               SxDiracVec<Complex16> muDen = PsiMu.colRef(iState).absSqr();
               signedResiduum += kWeight(ik) * fermi.focc(iState, iSpin, ik)
                  * (rBasis|(muDen - PWDen));
            }      
         }//iState
      }//iSpin
   }//ik

   Residuum = rBasis.symmetrize(Residuum);
   SxString file = "Residues.sxb";
   rBasis.writeMesh3d (file, Residuum);
   signedResiduum = rBasis.symmetrize(signedResiduum);
   file = "signedResidues.sxb";
   rBasis.writeMesh3d (file, signedResiduum);

}

#else // SX_STANDALONE

int main (int argc, char** argv)
{
   initSPHInXMath ();
   
   SxLoopMPI::init (argc, argv);  // LoopMPI
   
   cout.precision(10);

   // Command line parsing
   SxCLI cli (argc,argv);

   // Define Author
   cli.authors = "B. Lange";

   cli.preUsageMessage = "Orbital Optimization tool.";

   SxString quamolFile = cli.option ("-i|--input", "file", "Quamol input file")
                                    .toString ("quamol.sx");
   SxString outputFile = cli.option ("-o|--out", "file", "Quamol SXB file")
                                    .toString ("quamol.sxb");

   // factor by which log-grid is interpolated
   SxArray<int> refine = cli.option ("--refine","int","radial grid refinement")
                           .toIntList ();
   int printStep = cli.option ("--printStep","int","line print at step")
                             .toInt (0,0);
   bool noOpt = cli.option ("--noOpt",
                            "no optimization is performed ").toBool();
   bool checkGrad = cli.option ("--checkGrad",
                            "just check gradient ").toBool();
   bool checkTrafo = cli.option ("--checkTrafo",
                            "just check transformation ").toBool();
   int nPoints = cli.option ("--nPoints","Points for interpolation",
                             "Points for interpolation")
                            .toInt (100,0);
   bool plotResiduum = cli.option ("--plotResiduum","plotResiduum").toBool();
   cli.finalize ();


   // --- Define necessary variables 
   SxPtr<SxPW> wavesPtr;
   SxPtr<SxFermi> fermiPtr;
   SxKPoints kPoints;
   SxPtr<SxAtomicStructure> structPtr = SxPtr<SxAtomicStructure>:: create ();
   int nSpin=-1, nStates=-1, nk=-1;
   double ekt = -1.0;


   //--- read Quamol steering file
   SxParser quamolParser;
   SxConstPtr<SxSymbolTable> quamolTable = quamolParser.read (quamolFile);
   SxSymbolTable *quamolGroup = quamolTable->getGroup("Quamol");

   //--- Define Subspace
   if (quamolGroup->containsGroup("subSpace"))  {
      SxSymbolTable *cmd = quamolGroup->getGroup("subSpace");
      if (cmd->contains("waveFile"))  {
         SxString wavesFile = cmd->get("waveFile")->toString();
         try {
            SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
            wavesPtr = SxPtr<SxPW>::create (wavesFile, SxPW::InMemory);
            structPtr->read (io);
            kPoints.read(io);
            io.close ();
         }
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
         nSpin   = wavesPtr->getNSpin ();
         nStates = wavesPtr->getNStates ();
         nk      = wavesPtr->getNk();
         ekt     = 0.0;
      } else { 
         cout << "No reference wave file specified. SxQuamol quits here!" << endl;
         SX_QUIT;
      }
      if (cmd->contains("nStates"))  { 
         //occupation by lowest states
         nStates = cmd->get("nStates")->toInt ();
         SxArray<int> NperK (nk);
         for (int ik = 0; ik < nk; ik++)   {
            if (nStates > wavesPtr->getNStates(ik))   {
               sxprintf ("DFT run has %d states at maximal.\n", 
                         wavesPtr->getNStates(ik));
               sxprintf ("You choose %d states! SxQuamol ends here!\n",nStates);
               SX_QUIT;
            }
            NperK(ik) = nStates;
         }
         wavesPtr->setNStates(NperK);
         Focc fakeFocc (nStates,nSpin, nk);
         fakeFocc.set(2.0);
         fermiPtr = SxPtr<SxFermi>::create(2.0*nStates, nStates, nSpin, kPoints);
         fermiPtr->setOccupencies(fakeFocc);
      } else if (cmd->contains("bands"))  { 
         //ocupation by bands
         SxVector<Int> occ = cmd->get("bands")->toIntList ();
         Focc fakeFocc (nStates, nSpin, nk);
         fakeFocc.set(0.0);
         for (int iState = 0; iState < occ.getSize(); iState++)  {
            for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
               for (int ik = 0; ik < nk; ik++)  {
                  if (occ(iState)-1 >= nStates) {
                     cout << "Requested band not in reference calculation!" 
                          << endl;
                     SX_QUIT;
                  }
                  fakeFocc(occ(iState)-1, iSpin, ik) = 2.0;
               }
            }
         }
         nStates = occ.maxval();
         double nElec = 2.0 * int(occ.getSize ());
         fermiPtr = SxPtr<SxFermi>::create(nElec, nStates, nSpin, kPoints);
         fermiPtr->setOccupencies(fakeFocc);
      } else if (cmd->containsGroup("window"))  {

         SxString wavesFile = cmd->get("waveFile")->toString();
         SxFermi fermiRef;
         try {
            SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
            fermiRef = SxFermi (io);
            fermiRef.kpPtr = &kPoints;
            io.close ();
         }
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }

         SxSymbolTable *mode = cmd->getGroup("window");
         double lowEnergy = mode->get("lowEnergy")->toReal () * EV2HA;
         double highEnergy = mode->get("highEnergy")->toReal () * EV2HA;
         ekt = mode->contains("ekt") ?
               mode->get("ekt")->toReal()*EV2HA : 0.0;

         Focc fakeFocc = fermiRef.getFoccByWindow(lowEnergy, highEnergy, ekt);

         // Calculate Number of electrons
         double nElec = 0.0;
         for (int iState = 0; iState < fermiRef.getNStates(); iState++)  {
            for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
               for (int ik = 0; ik < nk; ik++)  {
                  nElec += fakeFocc(iState,iSpin,ik);
               }
            }
         }
         nElec /= double(nSpin) * double(nk);


         fermiPtr = SxPtr<SxFermi>::create(nElec, nStates, nSpin, kPoints);
         fermiPtr->eps = fermiRef.eps;
         fermiPtr->setOccupencies(fakeFocc);
         nStates = fermiPtr->getHOMO () + 1;

      } else { 
         //Fermi or standart
         SxString wavesFile = cmd->get("waveFile")->toString();
         SxFermi fermiRef;
         double nElec = -1.0;
         try {
            SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
            fermiRef = SxFermi (io);
            fermiRef.kpPtr = &kPoints;
            io.close ();
         }
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
         if (cmd->containsGroup("fermi")) {
            SxSymbolTable *fermiGroup = cmd->getGroup("fermi");
            nElec = fermiGroup->contains("nElectrons") ?
               fermiGroup->get("nElectrons")->toReal() :
               2.0 * fermiRef.getNValenceBands ();
            ekt = fermiGroup->contains("ekt") ?
               fermiGroup->get("ekt")->toReal()*EV2HA :
               0.0;
         } else {
            double bandOcc = (nSpin == 1) ? 2.0 : 1.0;
            nElec = bandOcc * fermiRef.getNValenceBands ();
            ekt = 0.0;
         } 
         // check if enough nStates have been calculated
         if (2 * nStates < nElec) {
            cout << "Too few states in reference calculation!" << endl;
            cout << "need " << nElec << ", but have " << (2 * nStates) << endl;
            SX_QUIT;
         }
         // Distribute Electrons
         fermiPtr = SxPtr<SxFermi>::create(nElec, nStates, nSpin, kPoints);
         fermiPtr->eps = fermiRef.eps;
         fermiPtr->fermiDistribution(ekt);
         nStates = fermiPtr->getNValenceBands ();
      }
   } else {
      cout << "No sub space specified. SxQuamol quits here!" << endl;
      SX_QUIT;
   }

   SxPW &waves = *wavesPtr;
   SxFermi &fermi = *fermiPtr;
   fermi.printOccupation ();
   SxAtomicStructure &structure = *structPtr;

   SxGkBasis &GkBasis = *waves.getGkBasisPtr();
   GkBasis.changeTau(structure);
   fermiPtr->kpPtr = waves.getGkBasisPtr ().getPtr();

   // Generate Basis depending on the underlying |g+k| values
   SxPtr<SxRadGBasis> radGBasisPtr =
      SxPtr<SxRadGBasis>::create(0.0, GkBasis.getMaxGk (), 
            nPoints, SxRadGBasis::Linear);
   
   SxArray<int> NperK (nk);
   for (int ik = 0; ik < nk; ik++)  NperK(ik) = nStates;
   waves.setNStates(NperK);
   double nElectrons = fermi.nElectrons;
   cout << "nElectrons = " << nElectrons << endl;
   cout << "distributet in nStates = " << nStates << endl;
   cout << "Number of kPoints = " << nk << endl;

   //Prepare refining
   // --- Refine mesh
   if (refine.getSize() == 0)  {
      int refineVal = 1;
      refine.resize(structure.getNSpecies());
      refine.set(refineVal);
   }
   if (refine.getSize() == 1)  {
      int refineVal = refine(0);
      refine.resize(structure.getNSpecies());
      refine.set(refineVal);
   }
   if (refine.getSize() != structure.getNSpecies())  {
      cout << "Wrong number of refinements in refine Mesh!" << endl;
      SX_QUIT;
   }

   // Normconserving or PAW Potential ?
   SxPtr<SxOverlapBase> SPtr;
   SxArray<SxString> chemName;
   
   if (quamolGroup->contains("potFile"))  {
      SxParser parser;
      SxString potFile = quamolGroup->get("potFile")->toString();
      SxConstPtr<SxSymbolTable> table = parser.read (potFile);
      //--- Norm conserving pseudo potential
      if (table->containsGroup("pseudoPot"))   {
         SxPtr<SxPseudoPot> psPotPtr = SxPtr<SxPseudoPot>::create (&*table);

         // radial basis from potential
         SxPtr<SxRadBasis> radBasisPotPtr 
            = SxPtr<SxRadBasis>::create(psPotPtr->rad,psPotPtr->logDr);

         // potential overlap operator
         SPtr = SxPtr<SxPWOverlap>::create ();

         chemName = psPotPtr->chemName;

      //--- Projector augmented wave potential
      } else if (table->containsGroup("pawPot"))   {
         SxPtr<SxPAWPot> pawPotPtr = SxPtr<SxPAWPot>::create (&*table);

         // radial Basis from potential
         SxPtr<SxRadBasis> radBasisPotPtr 
            = SxPtr<SxRadBasis>::create(pawPotPtr->rad,pawPotPtr->logDr);
         
         //--- potential overlap Operator
         SxPtr<SxPartialWaveBasis> pBasis
            = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
         pBasis->createProjBasis (GkBasis);
         SPtr = SxPtr<SxPAWOverlap>::create (pBasis, pawPotPtr);

         chemName = pawPotPtr->chemName;

      //--- Other potentials not supported, yet   
      } else   {
         cout << "No known Potential Group found!" << endl;
         SX_QUIT;
      }
   } else {
      cout << SX_SEPARATOR;
      cout << "WARNING: No potential specified, " 
           << "SxQuamol assumes norm conserving pseudopotential here!" 
           << endl;
      cout << SX_SEPARATOR;
      // potential overlap operator
      SPtr = SxPtr<SxPWOverlap>::create ();
   }

   //--- Initialize basis and create quamol
   SxConstPtr<SxRadBasis> radBasisPtr;
   SxAtomicOrbitals initGuess;
   SxQuamol quamol;
   if (quamolGroup->containsGroup("basis"))  {
      SxSymbolTable *basisGroup = quamolGroup->getGroup("basis");
      radBasisPtr = SxConstPtr<SxRadBasis>::create(basisGroup);
      quamol.rMax = radBasisPtr->getRMax ();
      initGuess.setup(basisGroup);
      initGuess.refine(refine);

      quamol.set(initGuess, radGBasisPtr, wavesPtr, fermiPtr, SPtr);
   } else {
      cout << "No basis specified, SxQuamol quits here!" << endl;
      SX_QUIT;
   }

   // --- spillage
   if (quamolGroup->containsGroup("spillage"))  {
      SxSymbolTable *spillageGroup = quamolGroup->getGroup("spillage");
      quamol.sigma = spillageGroup->contains("sigma") ?
                      spillageGroup->get("sigma")->toReal () : 1.0;   
   } else {
      quamol.sigma = 1.0;
   }

   // --- Kinetic Energy
   if (quamolGroup->containsGroup("kineticEnergy"))  {
      SxSymbolTable *eKinGroup = quamolGroup->getGroup("kineticEnergy");
      quamol.zeta = eKinGroup->contains("zeta") ?
                      eKinGroup->get("zeta")->toReal () : 0.0; 
      quamol.adaptiveZeta = eKinGroup->contains("adaptive") ? true : false;  
   } else {
      quamol.zeta = 0.0;
   }

   // --- Localization
   if (quamolGroup->containsGroup("localization"))  {
      SxSymbolTable *locGroup = quamolGroup->getGroup("localization");
      quamol.kappa = locGroup->contains("kappa") ?
                      locGroup->get("kappa")->toReal () : 0.0;
      quamol.rStart = locGroup->contains("rStart") ?
                      locGroup->get("rStart")->toReal () : 0.0;
      quamol.adaptiveKappa = locGroup->contains("adaptive") ? true : false;
   } else {
      quamol.kappa = 0.0;
   }
   quamol.setFixedList(quamolGroup->getGroup("basis"));

   // --- Final Settings (main)
   if (quamolGroup->containsGroup("main"))   {
      SxSymbolTable *mainGroup = quamolGroup->getGroup("main");
      quamol.dF = mainGroup->contains("dF") ?
                      mainGroup->get("dF")->toReal () : 1e-6;
      quamol.dRes = mainGroup->contains("dRes") ?
                      mainGroup->get("dRes")->toReal () : 1e-6;
      quamol.relLineMin = mainGroup->contains("dRelLineMin") ?
                      mainGroup->get("dRelLineMin")->toReal () : 1;
      quamol.print = mainGroup->contains("print") ? true : false;
      quamol.maxSteps = mainGroup->contains("maxSteps") ?
                      mainGroup->get("maxSteps")->toInt () : 100;
      quamol.printStep = printStep;
      quamol.checkGrad = checkGrad;
   }

   if (checkTrafo) quamol.checkTrafo (quamol.radials);

   //calculate Spillage
   quamol.computeFunctional (quamol.radials, false, true);
   cout << SX_SEPARATOR;
   cout << "Functional: " << quamol.functionalValue << endl;
   cout << SX_SEPARATOR;
   cout << "Spillage: " << quamol.spillage << endl;
   cout << SX_SEPARATOR;
   cout << "Kinetic Energy: " << quamol.eKin << endl;
   cout << SX_SEPARATOR;
   cout << "Localization: " << quamol.locVal << endl;
   cout << SX_SEPARATOR;
      
   SxAtomicOrbitals radials = quamol.getOrbitals(quamol.radials,radBasisPtr);
   radials.normalize ();
   radials.print("InitialGuess");
   
   if (!noOpt) quamol.compute ();
   
   quamol.computeFunctional (quamol.radials, false, true);
   cout << SX_SEPARATOR;
   cout << "Functional: " << quamol.functionalValue << endl;
   cout << SX_SEPARATOR;
   cout << "Spillage: " << quamol.spillage << endl;
   cout << SX_SEPARATOR;
   cout << "Kinetic Energy: " << quamol.eKin << endl;
   cout << SX_SEPARATOR;
   cout << "Localization: " << quamol.locVal << endl;
   cout << SX_SEPARATOR;

   quamol.printSpillagePerState(quamol.radials);

   if (checkTrafo) quamol.checkTrafo (quamol.radials);
   if (plotResiduum) quamol.calcResidues (quamol.radials);

   radials = quamol.getOrbitals(quamol.radials,radBasisPtr);
   radials.normalize ();

   SX_MPI_MASTER_ONLY {
      radials.print("Quamol");
      quamol.radials.print("Quamol-G");
      radials.write(outputFile);
   }

   quamol.getNormContribution(quamol.radials,radBasisPtr);

#ifdef USE_HDF5
   SX_MPI_MASTER_ONLY {radials.writeHDF5("quamol.hdf5");}
#endif

   SxTimer::getGlobalTimer().print ();
   cout << "SxQuamol completes successfully" << endl;

   return 0;
}

#endif /* SX_STANDALONE */
