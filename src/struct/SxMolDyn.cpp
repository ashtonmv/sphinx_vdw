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

#include <SxMolDyn.h>
#include <SxPtr.h>
#include <SxHessianOps.h>
#include <SxTimer.h>
#include <SxTaylorExpPotential.h>
#include <SxTime.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------

// --- Embedded Class SxIntegrator 
SxMolDyn::SxIntegrator::SxIntegrator()
{
   // empty
}

   SxMolDyn::SxIntegrator::SxIntegrator 
(const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   //--- Parsing in input parameters
   
   if (cmd->containsGroup ("integrator")) {
      SxSymbolTable *initI = cmd->getGroup ("integrator");
      if (initI->contains ("scheme")) scheme 
         = initI->get("scheme") -> toString();
      if (initI->contains ("order")) order = initI->get("order")-> toInt();
      if (initI->contains ("dt")) dt = initI->get("dt")-> toReal();
      if (initI->contains ("timeSteps")) 
         timeSteps = initI -> get ("timeSteps") -> toInt ();
   } else {
      cout << "Please set up integrator !" << endl;
   }
   
   //--- Validation of input parameters
   validateInput ();
   
}

void SxMolDyn::SxIntegrator::validateInput ()
{
   if (    (scheme != SxString("PCGear"))
         && (scheme != SxString("PCABM"))
         && (scheme != SxString("BBK"))
         && (scheme != SxString("Verlet")) ) 
   {
      cout << "MD: Unknown Integrating scheme: " << scheme << endl;
      SX_EXIT;
   }
}

void SxMolDyn::SxIntegrator::printParameters ()
{
   cout << "**************************" << endl;
   cout << "Integrator --- Parameters" << endl;
   cout << "-------------------------" << endl;
   cout << "scheme         : " << scheme << endl;
   cout << "order          : " << order << endl;
   cout << "dt             : " << dt << endl;
   cout << "timeSteps      : " << timeSteps << endl;
   cout << "**************************" << endl;
   cout << endl;
}

void SxMolDyn::SxIntegrator::initialize (SxState *state, SxState *oldState) 
{
   int i, j;
   //TODO timeDamping
   double timeDamping = 1.;
   
   if (scheme == SxString("PCGear")) {
      
      gearPredictorMatrix = SxMatrix<Double> (order + 1, order + 1 );
      gearPredictorMatrix.set(0.);
      
      // setting up Gear - Predictor -  matrix
      for (i = 0; i <= order; i++) {
         for (j = 0; j <= order; j++) {
            gearPredictorMatrix(i, j) = 0.;
            if (i > j)  gearPredictorMatrix(i, j) 
               = 0.;
            if (i == j) gearPredictorMatrix(i, j) 
               = 1.;
            if (i < j)  gearPredictorMatrix(i, j) 
               = 1./(double) parent -> getFaculty (j - i )
                  * ::pow (dt/(timeDamping), (j - i));
         }
      }
      
      gearCoeff = SxVector<Double> (order + 1);
      // --- Setting up Gear Coefficients: ref (2) Table 4.1
      
      if (order == 3)  {
         gearCoeff(0) = 1./6.; gearCoeff(1) = 5./6.;
         gearCoeff(2) = 1.; gearCoeff(3) = 1./3.;
      }
      
      if (order == 4)  {
         gearCoeff(0) = 19./120.; gearCoeff(1) = 3./4.;
         gearCoeff(2) = 1.;gearCoeff(3) = 1./2.;
         gearCoeff(4) = 1./12.;
      }
      
      if (order == 5)  {
         // Haile reports: 	gearCoeff(0) = 3./16.;
         gearCoeff(0) = 3./20.; gearCoeff(1) = 251./360.;
         gearCoeff(2) = 1.; gearCoeff(3) = 11./18.;
         gearCoeff(4) = 1./6.; gearCoeff(5) = 1./60.;
      }
   }

   gearDerivatives = SxMatrix<Double> (order + 1, parent -> nDoF);
   gearDerivatives.set (0.);
   gearDerivatives = gearDerivatives.transpose ();
   gearDerivatives.colRef(0).set(parent -> state.X.coordRef());
   gearDerivatives.colRef(1).set(parent -> state.V.coordRef());
   gearDerivatives = gearDerivatives.transpose ();
   
   thermoDerivatives = SxVector<Double> (order + 1);
   thermoDerivatives.set (0);

   
   if (scheme == SxString("BBK")) {
     
      parent -> thermostat.apply (state);
      
      SxAtomicStructure Xold, X, V, F, dFdt, A, dAdt, Xnew, Vnew, XRand, XRandInv, deltaX; 
     
      //--- getting state variables 
      X.copy ((*state).X);
      V.copy ((*state).V);
      A.copy ((*state).A);
      XRand.copy ((*state).XRand);
      XRandInv.copy ((*state).XRandInv);
      Xold.copy ((*oldState).X);

      
      //--- getting random variables from thermostat
      double GLD = parent -> thermostat.GLD;
      double HLD = parent -> thermostat.HLD;
      double gammaLD = parent -> thermostat.gammaLD;
      double eMin = exp(-gammaLD*dt);
      
      
      Xnew =        X 
           +     (1./gammaLD * (1. - eMin) )*   V  
           +     (1./gammaLD/ gammaLD*(gammaLD*dt - (1. - eMin))) * A
           +       XRand;
           

      deltaX = Xnew - Xold;

      deltaX.set (parent -> randomFilter | deltaX.coordRef ());
      Xnew = Xold + deltaX;

      
      Vnew =  (HLD/dt)*(     Xnew - Xold
           +                   (1./ gammaLD/gammaLD*GLD) * A
           -                      XRand)   ;

      //cout << Vnew << endl;
      //cout << HLD << endl;
      
     
     (*oldState).set (*state);
     (*state).X.copy (Xnew);
     (*state).V.copy (Vnew);
      
      parent -> thermostat.apply (state);
      
   }
}



void SxMolDyn::SxIntegrator::initPCABM ()
{
   cout << "PCABM not implemented yet !" << endl;
   //TODO: make compatible to new structure class
   //	initTauVec = tau2vec (initTau);
   
   //	for (int i = 3; i >= 0; i--) {
   //		timeEvolStruc(i) = timeEvolStruc(i + 1);
   //		timeEvolVel(i) = timeEvolVel(i + 1);
   //		timeEvolForc(i) = timeEvolForc(i + 1);
   //	}
   SX_QUIT;
}

void SxMolDyn::SxIntegrator::initVerlet ()
{
   SX_QUIT;
}

void SxMolDyn::SxIntegrator::resetGearPredictorMatrix(double deltaTime) 
{
   // TODO: incorporate time-damping
   
   int i, j;
   double timeDamping = 1.; // just a dummy yet (see above)
   for (i = 0; i <= order; i++) {
      for (j = 0; j <= order; j++) {
         gearPredictorMatrix(i, j) = 0.;
         if (i > j)  gearPredictorMatrix(i, j)  = 0.;
         if (i == j) gearPredictorMatrix(i, j)  = 1.;
         if (i < j)  gearPredictorMatrix(i, j)  
            = 1./(double) parent -> getFaculty (j - i )
               * ::pow (deltaTime/(timeDamping), (j - i));
      }
   }
}
   
   void 
SxMolDyn::SxIntegrator::predictNewState (SxState *state, SxState *oldState)
{
   if (scheme == SxString("PCGear")) {
      //--- make prediction for positions, velocities 
      //    and higher order derivatives
      //    Ref: SFHIngX 1.0 manual Eq. 11.9-11.12 
      
      gearDerivatives = gearPredictorMatrix^gearDerivatives;
      thermoDerivatives = gearPredictorMatrix^thermoDerivatives;
      
      gearDerivatives = gearDerivatives.transpose ();
      (*state).A.set(gearDerivatives.colRef(2));
      (*state).V.set(gearDerivatives.colRef(1));
      (*state).X.set(gearDerivatives.colRef(0));
      gearDerivatives = gearDerivatives.transpose ();
      
      (*state).thermoA = thermoDerivatives(2);
      (*state).thermoV = thermoDerivatives(1);
      (*state).thermoX = thermoDerivatives(0);
   }

   
   if (scheme == SxString("BBK") ) {
      
      SxAtomicStructure Xold, X, V, F, dFdt, A, dAdt,
                        Xnew, Vnew, XRand, XRandInv, deltaX; 
      
      //--- getting state variables 
      X.copy ((*state).X);
      V.copy ((*state).V);
      A.copy ((*state).A);
      XRand.copy ((*state).XRand);
      XRandInv.copy ((*state).XRandInv);
      Xold.copy ((*oldState).X);

      
      //--- computing time derivative of accelerations 
      dAdt.copy ((*state).X);
      dAdt. set ((A.coordRef () - (*oldState).A.coordRef ())/dt);

      //--- getting random variables from thermostat
      double gammaLD = parent -> thermostat.gammaLD;
      double GLD = parent -> thermostat.GLD;
      double HLD = parent -> thermostat.HLD;

      //--- approximated coefficients, stable for gammaLD*dt > 0.05
      //double eMin = exp(-gammaLD*dt);
      //double eMin2 = (0.5*gammaLD*dt*(1. + eMin) - (1. - eMin)); 
      
     //--- power series expansion for coefficients, works also 
     //    for gammaLD*dt <0.05
      double x = gammaLD*dt;
      double eMin = -( x - 1./2*pow(x, 2.) + 1./6.*pow(x, 3.) 
                  +    1./24.*pow(x, 4.) - 1./120.*pow(x, 5.) - 1.);

      double eMin2 = 1./12.*pow(x, 3.) - 1./24.*pow(x, 4.) 
                   + 1./80.*pow(x, 5.) - 1./360.*pow(x, 6.);
      
      //cout << gammaLD << endl;
      //cout << X << endl;
      
      //--- computting new positions, velocities and accelerations
      Xnew =        (1. + eMin)             *X 
           -              eMin              *Xold 
           +      ( (1. - eMin)                        *dt/gammaLD) *A
           +      (eMin2*dt/gammaLD/gammaLD)   *dAdt
           +                                                XRand
           +         eMin                                 * XRandInv;


      deltaX = Xnew - Xold;

      //parent -> randomFilter.set (Xold, parent->performance.freezeMode, true);
      
      deltaX.set (parent -> randomFilter | deltaX.coordRef ());
      Xnew = Xold + deltaX;

      
      Vnew =  (HLD/dt)*(     Xnew - Xold
           +                   (1./ gammaLD/gammaLD*GLD) * A
           -        (1./gammaLD/gammaLD/gammaLD*GLD)     * dAdt
           +                     (XRandInv - XRand)   );

      //cout << Vnew << endl;
      //cout << HLD << endl;
      
     
     (*oldState).set (*state);
     (*state).X.copy (Xnew);
     (*state).V.copy (Vnew);
     
   }  

   
}

void 
SxMolDyn::SxIntegrator::correctNewState (SxState *state, SxState * /*oldState*/)
{
   SxAtomicStructure correction;
   double correctionThermo, predictedThermoA;
   SxAtomicStructure predictedA;
   SxVector<Double> correctedVec(parent -> nDoF);
   int i;
   
   //TODO: timeDamping
   double timeDamping = 1.;
   if (scheme == SxString("PCGear")) {
      
      // --- Gear - Corrector: ref(2) equ. 4.50 
      predictedA.copy   ((*state).A);
      (*state).A.set    ((*state).F.coordRef()/parent -> massVec);
      
      predictedThermoA = (*state).thermoA;
      (*state).thermoA = (*state).thermoF/(parent -> thermostat.mass);
      // ref(2) equ. 4.57
      correction =  (::pow(dt/timeDamping, 2.) / 2.)
         *((*state).A - predictedA);
      
      correctionThermo = ((*state).thermoA - predictedThermoA)
         * ::pow( dt/timeDamping, 2.) / 2.;
      
      
      // --- correction of predictions: ref (2) equ. 4.51 - 4.56
      
      //--- transposing gear derivatives is just a workaround (see above)
      gearDerivatives = gearDerivatives.transpose ();
      
      for (i = 0; i <= order; i++) {
         gearDerivatives.colRef(i).set 
            ( gearDerivatives.colRef(i) 
              + gearCoeff(i)*(double)(parent -> getFaculty(i))
              / ::pow (dt/(timeDamping), i) * correction.coordRef ());
         
         
         
         thermoDerivatives(i) += correctionThermo*gearCoeff(i)
            * (double) parent -> getFaculty(i)
            /  ::pow (dt/(timeDamping), i);
      }
      
      (*state).X.set(gearDerivatives.colRef(0));
      (*state).V.set(gearDerivatives.colRef(1));
      
      gearDerivatives = gearDerivatives.transpose ();
   }
   
   if (scheme == SxString("BBK")) {
      // --- only the accelerations need to be updated here
      (*state).A.set    ((*state).F.coordRef()/parent -> massVec);
   } 
}

SxMolDyn::SxIntegrator::~SxIntegrator()
{
   // empty:wq
   
}

// --- Embedded Class SxImpSampler

SxMolDyn::SxImpSampler::SxImpSampler()
{
   // empty
}

SxMolDyn::SxImpSampler::SxImpSampler (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   isOn = false;
   if (cmd->containsGroup ("impSampler")) {
      isOn = true;
      SxSymbolTable *iS = cmd -> getGroup("impSampler");
      nTrust = iS -> get("nTrust") -> toInt ();
   }
   
   
}

void SxMolDyn::SxImpSampler::initialize () 
{

   SxArray<const SxSymbolTable *> dummyT (1);
   dummyT(0) = NULL;
   SxAtomicStructure dummy;

   iSteps = 0;
   startState = SxState (parent);
   startOldState = SxState (parent);
   startState.set (parent -> state);
   startOldState.set (parent -> oldState);
   
   //--- getting energy of target potential
   ePotTS = parent -> getEPot ();

   //--- getting energy of reference potential
   dummy = parent -> tPotential -> getForces(startOldState.X, dummyT(0));
   ePotRS =  parent -> tPotential -> getEnergy ();

   //--- switching on reference potential
   lambda = (parent -> guess).lambda;
   (parent -> guess).lambda = 0.;
} 

void SxMolDyn::SxImpSampler::inc ()
{
   iSteps ++;
}

void SxMolDyn::SxImpSampler::apply ()
{
   SxArray<const SxSymbolTable *> dummyT (1);
   dummyT(0) = NULL;
   SxAtomicStructure dummy;
   
   if (iSteps >= nTrust) {
     //--- getting energy of target potential
      (parent -> guess).lambda = lambda; 
     dummy = parent -> getForces(parent -> state.X);
     ePotTF = parent -> getEPot ();
     (parent -> guess).lambda = 0.; 
     
     //--- getting energy of reference potential
     dummy = parent -> tPotential -> getForces(parent -> state.X, dummyT(0));
     ePotRF =  parent -> tPotential -> getEnergy ();

     double ddE = (ePotTF - ePotRF) -(ePotTS - ePotRS);
     double kBT = (parent -> thermostat.T)*KB*JL2HA;

     double boltzmann = exp(-ddE/kBT);
     
     
     cout << "TEST ImpSampler " << endl;
     cout << "Starting chain: target: " << ePotTS << "reference: " << ePotRS << endl;

     cout << "ending chain: target: " << ePotTF << "reference: " << ePotRF << endl;
     cout << "ddE : " << ddE << "kBT: " << kBT << "bm: " << boltzmann << endl; 


     if (boltzmann > 1.) boltzmann = 1.;

    // SxRandom random;
   //  double rand = random.get ();


     double random = rand()/double(RAND_MAX);

     if (random <= boltzmann) {
        cout << "VAL : " << (parent -> state.dUdL) << endl;
        ePotRS = ePotRF;
        ePotTS = ePotTF;
        startState.set (parent -> state);
        startOldState.set (parent -> oldState);
     } else {
        (parent -> state).set (startState);
        (parent -> oldState).set (startOldState);
     }
        
        

      

     iSteps = 0;
  }

}

SxMolDyn::SxImpSampler::~SxImpSampler()
{
   // empty
}

// --- Embedded Class SxTDIntegration

SxMolDyn::SxTDIntegration::SxTDIntegration()
{
   // empty
}

SxMolDyn::SxTDIntegration::SxTDIntegration (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   if (cmd->containsGroup ("tdIntegration")) {
   } 
   
   //--- Validation of input parameters
   validateInput ();
   
}


void SxMolDyn::SxTDIntegration::printParameters ()
{
   cout << "**************************" << endl;
   cout << "TDIntegration --- Parameters" << endl;
   cout << "-------------------------" << endl;
   cout << "**************************" << endl;
   cout << endl;
}

void SxMolDyn::SxTDIntegration::validateInput ()
{
   //--- empty
}


void SxMolDyn::SxTDIntegration::initialize () 
{
   //--- empty
}

void SxMolDyn::SxTDIntegration::switchOn () 
{
   switchedOff =  false;
   /*
   dUdLAvg = 0.;
   sumSteps = 0.;
   sampleLambda.resize (3);
   sampledUdL.resize (3);
   sampleLambda(0) = 1e-5;
   sampleLambda(1) = 0.5;
   sampleLambda(2) = (1 - 1e-5);
   nSamplePoints = 3;
   iSamplePoints = 0;
   */
}


void SxMolDyn::SxTDIntegration::goOn ()
{

   switchedOff = true;
   /*
  // --- sample on given mesh
  if (iSamplePoints < nSamplePoints) {
     sampledUdL(iSamplePoints) = dUdLAvg;
     dUdLAvg = 0.;
     sumSteps = 0.;
     iSamplePoints++;
  }
     
  
  else {
     // --- print out 

     sxprintf 
     
     // --- extend sample mesh
     iSamplePoints = 0;
     SxArray<double> newSampleLambda (nSamplePoints + 2);
     SxArray<double> newSampledUdL (nSamplePoints + 2);
     for (int i = 0; i < nSamplePoints; i++) {
        newSampleLambda(i) = sampleLambda(i);
        newSampledUdL (i) = sampledUdL(i);
     }

     // --- refine for sampling singularity at lambda = 1
     newSampleLambda (nSamplePoints) =
        (sampleLambda(nSamplePoints - 1) + sampleLambda(nSamplePoints - 2))/2.;
     newSampleLambda (nSamplePoints + 1) =
        (sampleLambda(nSamplePoints - 2) + sampleLambda(nSamplePoints - 3))/2.;

     nSamplePoints += 2;
     
     for (int i = 0; i < nSamplePoints; i++) {
        sampleLambda(i) = newSampleLambda(i);
        sampledUdL (i) = newSampledUdL(i);
     }

     bool unsorted = true;
     double cValue;
     
     while (unsorted) {
        unsorted = false;
        for (int i = 0; i < (nSamplePoints - 1); i++) {
           if (sampleLambda(i) > sampleLambda(i)) {
              cValue = sampleLambda(i);
              sampleLambda(i) = sampleLambda(i + 1);
              sampleLambda(i + 1) =  cValue;
              cValue = sampledUdL(i);
              sampledUdL(i) = sampledUdL(i + 1);
              sampledUdL(i + 1) =  cValue;
              unsorted = true;
           }
        }
     }   
     if (nSamplePoints > 11) switchedOff = true;
  }
  */
}

void SxMolDyn::SxTDIntegration::apply (SxGuess * /*guess */)
{
   //guess.lambda = sampleLambda (iSamplePoints);
}

void SxMolDyn::SxTDIntegration::takeMeasurement (const SxState &state)
{
   double dUdL = state.dUdL;
   sumSteps += 1.;
   dUdLAvg = (dUdLAvg * (sumSteps - 1.) + dUdL)/sumSteps;  
}

SxMolDyn::SxTDIntegration::~SxTDIntegration ()
{
   // empty
}

// --- Embedded Class SxThermostat

SxMolDyn::SxThermostat::SxThermostat()
{
   // empty
}

SxMolDyn::SxThermostat::SxThermostat (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   type = SxString("none");
   mass = 1.0;
   gammaLD = 1.0;
   nHarmEqu = -1;
   iHarmEqu = 0;
   gammaLDStore = 1.0;
   
   if (cmd->containsGroup ("thermostat")) {
      SxSymbolTable *initI = cmd->getGroup ("thermostat");
      
      if (initI->contains ("mass")) mass = initI->get("mass") -> toReal();
      if (initI->contains ("gammaLD")) {
         gammaLD = initI->get("gammaLD") -> toReal();
         gammaLDStore = gammaLD;
      }
      if (initI->contains ("scheme")) scheme = initI->get("scheme") -> toString();
      if (initI->contains ("temperature")) { 
         type = SxString("constantTemperature");
         T = initI->get("temperature")-> toReal();
         startTemp = T;
         //scheme = SxString ("NH");
      }
      if (initI->contains ("rescale")) 
         scheme = SxString("rescale");
      if (initI->containsGroup ("annealing")) {
         type = SxString("annealing");
         if (!(initI->contains ("rescale"))) {
            if (scheme == SxString("langevin")) {
               cout << "Annealing not tested for langevin dynamics yet. Please use NH-Thermostat " 
                    << "or contact SPHInX website." << endl;
               SX_QUIT;
            } else {
               scheme = SxString ("NH");
            }
         }
         tempSteps = initI -> getGroup ("annealing") 
            -> get ("tempSteps") -> toInt ();
         startTemp = initI -> getGroup ("annealing") 
            -> get ("startTemp") -> toReal ();
         endTemp = initI -> getGroup ("annealing") 
            -> get ("endTemp") -> toReal ();
         increment = initI -> getGroup ("annealing") 
            -> get ("increment") -> toString ();
      }
      
      if (initI->contains ("harmonicEquSteps")) {
         nHarmEqu = initI->get("harmonicEquSteps") -> toInt();
         if ((nHarmEqu > 0 ) && (nHarmEqu < 500))  {
            cout << "harmonicEquSteps should be 0 or > 500." << endl;
            SX_QUIT;
         }
      }
      
      
   } 
   
   //--- Validation of input parameters
   validateInput ();
   
}


void SxMolDyn::SxThermostat::printParameters ()
{
   cout << "**************************" << endl;
   cout << "Thermostat --- Parameters" << endl;
   cout << "-------------------------" << endl;
   cout << "modus          : " << type << endl;
   if (type == SxString("annealing")) {
      cout << "mass           : " << mass << endl;	
      cout << "startTemp      : " << startTemp << endl;
      cout << "endTemp        : " << endTemp << endl;
      cout << "tempSteps      : " << tempSteps << endl;
      cout << "increment      : " << increment << endl;
   } else {
      cout << "scheme         : " << scheme << endl;
      cout << "temperature    : " << T << endl;
   }
   cout << "**************************" << endl;
   cout << endl;
}

void SxMolDyn::SxThermostat::validateInput ()
{
   if (endTemp < 1e-3) endTemp = 1e-3;
   if (startTemp < 1e-3) startTemp = 1e-3;
}

double SxMolDyn::SxThermostat::getInitialEkin ()
{
   return (startTemp*KB*JL2HA/2.*L);
}

void SxMolDyn::SxThermostat::initialize () 
{
   L = parent->nDoF - 3;
   alpha = exp(log(endTemp/startTemp)/((double)tempSteps - 1.));
   dT = (endTemp-startTemp)/((double)tempSteps - 1.);

   // --- parameters for Langevin dynamics

   double dt = parent -> integrator.dt;

   if (scheme == SxString("langevin")) {

      sigmaSqrMin = SxVector<Double> (parent -> nDoF);
      sigmaSqrPlu = SxVector<Double> (parent -> nDoF);

      //--- approximated coefficients, stable for gammaLD*dt > 0.05
     /* 
      double eMin  = exp(-gammaLD*dt);
      double e2Min = exp(-2.*gammaLD*dt);
      double ePlu  = exp(gammaLD*dt);
      double e2Plu = exp(2.*gammaLD*dt);

      CLD = 2.*gammaLD*dt - 3. + 4.* eMin - e2Min;
      GLD = ePlu - 2.*gammaLD*dt  - eMin;
      ELD = 16.*(ePlu  + eMin)
          -  4.*(e2Plu + e2Min)
          - 24.
          - 4.*gammaLD*dt*(ePlu - eMin)
          + 2.*gammaLD*dt*(e2Plu - e2Min);
      HLD = gammaLD * dt /(ePlu - eMin);
     */
     //--- power series expansion for coefficients, works also 
     //    for gammaLD*dt <0.05

      
      double x = gammaLD*dt;
      CLD = 2./3.*pow(x,3.) -1./2.*pow(x, 4.) + 7./30.*pow(x, 5.)
          - 1./12.*pow(x, 6.) +31./1260*pow(x,7.)-1./160.*pow(x, 8.)
          +127./90720.*pow(x, 9.);

      GLD = 1./3.*pow(x, 3.) + 1./60.*pow(x, 5.);

      HLD = 1./2. - 1./12.*pow(x, 2.) + 7./720.*pow(x, 4.);
      double ELDbyCLD = 1./2.*pow(x, 3.) +3./8.* pow(x, 4.) + 29./160.*pow(x,5.)
               + 43./640.*pow(x, 6.) + 1831./89600.*pow(x, 7.)
               + 381./71680.*pow(x, 8.) + 235009./193536000.*pow(x,9.);

      ELD = ELDbyCLD*CLD;
  

   }

}

void SxMolDyn::SxThermostat::switchOn () 
{
   switchedOff =  false;
   iSteps = 0;
   if (type == SxString("annealing")) {
      T = startTemp;
   }
}


void SxMolDyn::SxThermostat::goOn ()
{
   if (type == SxString("annealing")) {
      if (increment == SxString("linear")) {
            T += dT;
         }
         if (increment == SxString("exponential")) {
            T *= alpha;
         }
         iSteps += 1;
         if (iSteps >= tempSteps) switchedOff = true;
   }
   if (type == SxString("constantTemperature")) {
      switchedOff = true;
   }
   if (type == SxString("none")) {
      switchedOff = true;
   }
}

void SxMolDyn::SxThermostat::apply (SxState *state)
{
   //if (parent -> guess.lambda > 1e-7) storeLambda = parent -> guess.lambda;
   if (nHarmEqu > 0) {
      if (iHarmEqu == 0) storeLambda =  parent -> guess.lambda;

      if (iHarmEqu < nHarmEqu) {
         (parent -> guess.lambda) = 0.; 
/*         
if (iHarmEqu < 15) {
            parent -> guess.lambda = 
               storeLambda * ((1. - double(iHarmEqu)/15.));
         }
*/         
         
         if (iHarmEqu > (nHarmEqu - 10)) {
        
            if (( parent -> getEKin ((*state).V)) > (0.5*L*KB*T*JL2HA)) {
               iHarmEqu -= 2;
               gammaLD *= 3.;
            } else {
               gammaLD = gammaLDStore;
            }

            cout << "TEST: " <<( parent -> getEKin ((*state).V)) 
                 << " " << (0.5*L*KB*T*JL2HA)  
                 << " " << (nHarmEqu - iHarmEqu) << endl;

        
            double fac = (1. - double((nHarmEqu -1) - iHarmEqu)/10) ;
            
            parent -> guess.lambda = 
               storeLambda * fac;
         }
         
         iHarmEqu ++;
      } else {
         parent -> guess.lambda = storeLambda;
         gammaLD = gammaLDStore;
      }
   }
   
   
   if (scheme == SxString("NH")) {
      // V velocities
      // Eq. (11.33)  xi is thermoV
      SxVector<Double> a = ((*state).V.coordRef () 
            * (*state).thermoV)*(parent -> massVec);
      // F forces; thermoF: forces on the additional degree of freedom
      a = -a + (*state).F.coordRef ();
      (*state).F.set (a);
      (*state).thermoF = 2.* (parent -> getEKin((*state).V)) 
                             - L*KB*T*JL2HA;
            (*state).T = T;
   }
   

   SxVector<Double> XRand (parent -> nDoF);
   SxVector<Double> XRandInv (parent -> nDoF);
   if (scheme == SxString("langevin")) {

      for (int i = 0; i < (parent -> nDoF); i++) {
         sigmaSqrMin (i) = ELD/CLD * T * KB*JL2HA / (gammaLD * gammaLD *(parent -> massVec (i)) );
         sigmaSqrPlu (i) = CLD * T * KB*JL2HA / (gammaLD * gammaLD *(parent -> massVec (i)) );
         
      }
      
      //SxRandom random;
      double rand; 
      //double sign, norm;
      SxVector<Double> gauss (1000);
      gauss.set (0.);
      //double dx = 1e-2;

      
      for (int i = 0; i < (parent -> nDoF); i++) {

// cout << "RAND " << rand << endl;
         /*
         norm = 1./sqrt(2.*PI*sigmaSqrMin(i));
         rand = (0.5 - random.get ());
         //cout << "RAND " << rand << endl;
         rand = rand /0.5 *norm;
         if (rand < 0) 
            sign = -1.;
         else sign = 1.;
         rand = rand * sign;
         */

         rand = parent -> gaussDev (sigmaSqrMin(i));
         
         XRandInv(i) = GLD/CLD*(*state).XRand.coordRef ()(i)
                     + rand;
                     //(gammaLD*parent -> massVec(i));
                           
         
         /*                  
         cout << "XR " << XRandInv(i) << endl;
         cout << "XR " << ((*state).XRand.coordRef ()(i) ) << endl;
         cout << "XR " << (-log(rand/norm) ) << endl;
         cout << "XR " << (sigmaSqrMin(i)) << endl;
*/
 /*                          
         norm = 1./sqrt(2.*PI*sigmaSqrPlu(i));
         rand = (0.5 - random.get ());
         rand = rand /0.5 *norm;
         if (rand < 0) 
            sign = -1.;
         else sign = 1.;
         rand = rand * sign;
   */
        rand = parent -> gaussDev (sigmaSqrPlu(i));
  
         
         XRand(i) = rand;
                  ///(gammaLD*parent -> massVec(i))));

         //int index = (int)(XRand(i)/dx )+ 500;
//         int index = (int)(500. + (parent -> gaussDev (0.3))/dx) ;
//         if ((index < 1000) && (index >= 0)) {
           // cout << "G " << (parent -> gaussDev ()) << endl;
//            gauss(index) += 1.;
//         }
         
      }
     
      
     // parent -> randomFilter.set (state -> X, parent->performance.freezeMode, true);
      /*
      XRand .set 
            (parent -> randomFilter | XRand);

      XRandInv .set 
            (parent -> randomFilter | XRandInv);
            */

      (*state).XRand.set(XRand);
      (*state).XRandInv.set(XRandInv);
      (*state).T = T;
   } 
}

SxMolDyn::SxThermostat::~SxThermostat()
{
   // empty
}

// --- Embedded Class SxPerformance

SxMolDyn::SxPerformance::SxPerformance()
{
   // empty
}

SxMolDyn::SxPerformance::SxPerformance (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   freezeMode = SxString("trans");
   extrapolateWaves = false;
   applyConstraints = false;
   restart = false;
   constraints = 0.*constraints;
   if (cmd->containsGroup ("performance")) {
      SxSymbolTable *initI = cmd->getGroup ("performance");
      
      if (initI->contains ("freezeMode")) 
         freezeMode = initI->get("freezeMode") -> toString();
      if (initI->contains ("extrapolateWaves"))
         extrapolateWaves = initI -> get ("extrapolateWaves") -> toBool();
      if (initI->contains ("restart")) 
         restart = initI -> get ("restart") -> toBool();
      if (initI->contains ("constraints")) {
         applyConstraints = true;
         constraints = SxMatrix3<Double>
            (initI -> get ("constraints")->toList()).transpose();
      }
   } 
   
   //--- Validation of input parameters
   validateInput ();
}

void SxMolDyn::SxPerformance::initialize () 
{
   parent -> forceFilter.set (parent -> initialTau, freezeMode, false);
   //--- 26.11.2007 L. Ismer (performance-related change)
   //parent -> velocityFilter.set (parent -> initialTau, SxString("xyz"), false);
   parent -> velocityFilter.set (parent -> initialTau, freezeMode, false);
   parent -> randomFilter.setMasses (parent -> massVec);
   parent -> randomFilter.set (parent -> initialTau, freezeMode, false);
   
   if (extrapolateWaves) {
      xOld.copy (parent -> initialTau);
      xAct.copy (parent -> initialTau);
      xNew.copy (parent -> initialTau);
      actWaves = SxPW (parent -> hamSolver -> getWaves ());
   }
}

void SxMolDyn::SxPerformance::update (const SxState &state)
{
   if (extrapolateWaves) {
      xOld.copy (xAct);
      xAct.copy (state.X);
      oldWaves = SxPW (actWaves);
      actWaves = SxPW (parent -> hamSolver -> getWaves ());
      
   }
}

double SxMolDyn::SxPerformance::getAlpha (const SxAtomicStructure &xNew_)
{
   SxAtomicStructure xPrime;
   double rAA, rNO, rAO, rAN, rOO;
   double alpha, Z, N, left, right, middle;
   
   rAA = (xAct.coordRef () ^ xAct.coordRef ()).chop (); 
   rNO = (xNew_.coordRef() ^ xOld.coordRef ()).chop (); 
   rAO = (xAct.coordRef () ^ xOld.coordRef ()).chop (); 
   rAN = (xAct.coordRef () ^ xNew_.coordRef()).chop ();
   rOO = (xOld.coordRef () ^ xOld.coordRef ()).chop ();
   
   Z = (-rAA - rNO + rAO + rAN);
   N = (rAA - 2.*rAO + rOO);
   
   if (N != 0.) 
      alpha = Z/N;
   else {
      alpha = 0.;
   }
   
   xPrime = xAct + alpha*(xAct - xOld);
   xPrime = xPrime - xNew_;
   middle = (xPrime.coordRef () ^ xPrime.coordRef ()).chop ();
   
   xPrime = xAct + alpha/10.*(xAct - xOld);
   xPrime = xPrime - xNew_;
   right = (xPrime.coordRef ()^ xPrime.coordRef ()).chop ();
   xPrime = xAct - alpha/10.*(xAct - xOld);
   xPrime = xPrime - xNew_;
   
   left = (xPrime.coordRef ()^xPrime.coordRef ()).chop ();
   
   if (!((left > middle) && (right > middle))) {
      alpha = 0.;
   }
   return alpha;
}

int SxMolDyn::SxPerformance::meditateRestart (int nIt, SxState *oldState, SxState *state)
{
        if (restart) {


      parent -> thermostat.nHarmEqu = 0;
      parent -> guess.lambda = parent -> thermostat.storeLambda; 
      int nSpecies = (*oldState).X.getNSpecies ();
      SxVector<Int> nAtoms (nSpecies);
      nAtoms.set ((*oldState).X.getNAtoms ());

      restart = false;
   //--- reads in whole file before pushing to trajectory (should be refined)
      ifstream inStream;
      string line;
      int fileLength;
      //char *buffer;
      SxAtomicStructure X;
      SxVector3<Double> coord;
      SxList<SxString> xStringTok, coords;
      SxList<SxString>::Iterator itX;
   
      inStream.open ("restart.dat");
      inStream.seekg (0,ios::end);
      fileLength = int(inStream.tellg ());
      if (fileLength <= 0)  {
         cout << "Fatal error for restart: restart.dat not found or empty"
              << endl;
         SX_QUIT;
      }
   
      SxArray<char> buffer(fileLength);
      inStream.seekg (0, ios::beg);
      inStream.read (buffer.elements, fileLength);
   
      SxString whole = SxString(buffer.elements);
   
      
      SxString struc = (whole.right ("structureOld#")).left ("\n#");
      SxList<SxString> strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X = SxAtomicStructure ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref (j) = coord;
         itX++;
      }
      (*oldState).X.copy(X);
      
      
      struc = (whole.right ("velocityOld#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X = SxAtomicStructure ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref (j) = coord;
         itX++;
      }
      (*oldState).V.copy(X);
      
      struc = (whole.right ("accelerationOld#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X = SxAtomicStructure ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref (j) = coord;
         itX++;
      }
      (*oldState).A.copy(X);
      
      
      struc = (whole.right ("structure#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref(j) = coord;
         itX++;
      }
      (*state).X.copy(X);
      
      
      struc = (whole.right ("velocity#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref(j) = coord;
         itX++;
      }
      (*state).V.copy(X);
      
      struc = (whole.right ("acceleration#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref(j) = coord;
         itX++;
      }
      (*state).A.copy(X);
      
      struc = (whole.right ("XRand#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref(j) = coord;
         itX++;
      }
      (*state).XRand.copy(X);
      
      struc = (whole.right ("XRandInv#")).left ("\n#");
      strucTok = struc.tokenize ('\n');
      itX = strucTok.begin ();
      X.copy ((*state).X);
      for (int j = 0; j < ((parent -> nDoF)/3); j++) {
         coords = (*itX).tokenize(' ');
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.ref(j) = coord;
         itX++;
      }
      (*state).XRandInv.copy(X);
      
       struc = (whole.right ("step#")).left ("\n#");

       int n = (int)(struc.toDouble ());

       cout << (*state).X.getNSpecies ();
       cout << (*state).X.getNAtoms ();
       
      
      return n; 
   } else {
      //cout << "NORESTART" << endl;
      return nIt;
   }
}


void SxMolDyn::SxPerformance::printParameters ()
{
   cout << "**************************" << endl;
   cout << "Performance --- Parameters" << endl;
   cout << "-------------------------" << endl;
   cout << "freezeMode        : " << freezeMode << endl;	
   cout << "extrapolateWaves : " << extrapolateWaves << endl;
   cout << "constraints      : " << constraints << endl;
   cout << "**************************" << endl;
   cout << endl;
}

void SxMolDyn::SxPerformance::validateInput ()
{
   if ((extrapolateWaves) && (parent -> hamSolver == NULL)) {
      cout << "extrapolating wavefunction does only work for "
           << " DFT potentials yet " << endl;
      SX_QUIT;
   }
   
}

void SxMolDyn::SxPerformance::extrapolateWavefunction 
(const SxAtomicStructure &newX)
{
   //empty
   double alpha;
   if (extrapolateWaves) {
      cout << "Extrapolating Wavefunction" << endl;
      alpha = getAlpha (newX);
      cout << "Alpha: " << alpha << endl;
      if (alpha > 1e-15) parent -> hamSolver -> extrapolateWaves 
         (alpha, oldWaves);
   }
}

SxMolDyn::SxPerformance::~SxPerformance()
{
   // empty
}

// --- Embedded Class SxPerturbation

SxMolDyn::SxPerturbation::SxPerturbation()
{
   // empty
}

SxMolDyn::SxPerturbation::SxPerturbation (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   applyPerturbation = false;
   
   if (cmd->containsGroup ("perturbation")) {
      SxSymbolTable *initI = cmd->getGroup ("perturbation");
      applyPerturbation = true;
      gas = initI->get("gas") -> toString ();
      dof = initI -> get ("dof") -> toInt ();
      interval = initI -> get ("interval") -> toInt ();
   }
   
   //--- Validation of input parameters
   validateInput ();
}

void SxMolDyn::SxPerturbation::printParameters ()
{
   if (applyPerturbation) { 
      cout << "**************************" << endl;
      cout << "Perturbation --- Parameters" << endl;
      cout << "-------------------------" << endl;
      cout << "gas              : " << gas << endl;
      cout << "dof              : " << dof << endl;
      cout << "interval         : " << interval << endl;
      cout << "**************************" << endl;
      cout << endl;
   }
}

void SxMolDyn::SxPerturbation::validateInput ()
{
   //empty
}

SxMolDyn::SxPerturbation::~SxPerturbation()
{
   // empty
}

// --- Embedded Class SxOutput

SxMolDyn::SxOutput::SxOutput()
{
   // empty
}

SxMolDyn::SxOutput::SxOutput (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   saveWaves = -1;
   saveHistory = true;
   saveRestartInfo = true;
   printSteps = -1;
   printStepsFreq = -1;
   inSectionsFreq = false;
   printStepsHessian = -1;
   inSectionsHessian = false;
   
   if (cmd->containsGroup ("output")) {
      SxSymbolTable *initI = cmd->getGroup ("output");
      if (initI -> contains ("saveWaves")) 
         saveWaves = initI->get("saveWaves") -> toInt ();
      if (initI -> contains ("saveHistory")) {
         saveHistory = initI->get("saveHistory") -> toBool ();
      }
      if (initI -> contains ("printSteps")) {
         printSteps = initI -> get("printSteps") -> toInt ();
      }
      if (initI -> contains ("saveRestartInfo")) 
         saveRestartInfo = initI->get("saveRestartInfo") -> toBool ();
      if (initI->containsGroup ("printGFS")) {
         SxSymbolTable *printGroup = initI->getGroup("printGFS");
         parent -> useTimeCorr = true;
         printStepsFreq = printGroup -> get("printSteps") -> toInt ();
         if (printGroup -> contains ("inSections")) {
            inSectionsFreq = printGroup -> get("inSections") -> toBool ();
         }
      }
      if (initI->containsGroup ("printHessian")) {
         SxSymbolTable *printGroup = initI->getGroup("printHessian");
         parent -> useTimeCorr = true;
         printStepsHessian = printGroup -> get("printSteps") -> toInt ();
         if (printGroup -> contains ("inSections")) {
            inSectionsHessian = printGroup -> get("inSections") -> toBool ();
         }
      }
      plotDih = false;
      if (initI->containsGroup ("plotDihedrals")) {
         SX_QUIT;
         SxSymbolTable *printGroup = initI->getGroup("plotDihedrals");
         C_one = printGroup -> get("C_one") -> toInt ();
         C_two = printGroup -> get("C_two") -> toInt ();
         N_one = printGroup -> get("N_one") -> toInt ();
         N_two = printGroup -> get("N_two") -> toInt ();
         C_alpha = printGroup -> get("C_alpha") -> toInt ();
         plotDih = true;
      }
   }
   
   //--- Validation of input parameters
   validateInput ();
   
   if ( ! (parent -> performance.restart) ) {
   	if ( !(fpMoldyn = fopen ("moldyn.dat", "w" )) ) {
         sxprintf ("Can't open file moldyn.dat");
         SX_EXIT;
      }
   
      fprintf (fpMoldyn, "# time(femtosec) ePot(H=hartree) eKin(H)     ");
      fprintf (fpMoldyn, "    eTot(H)         dUdL-harm(H)    dUdL-harm-E0(H)");
      fprintf (fpMoldyn, "  dudL-harm-E0(meV)  dUdL-harm-E0(meV/atom)\n");
      fprintf (fpMoldyn, "# -------------------------------------------");
      fprintf (fpMoldyn, "---------------------------------------------------");
      fprintf (fpMoldyn, "-----------------------------------------\n");
      fclose(fpMoldyn);
   
      if ( !(fpMoldynHist = fopen ("moldynHist.dat", "w")) ) {
         sxprintf ("Can't open file moldynHist.dat");
         SX_EXIT;
      }	
      fclose(fpMoldynHist);
   }
}

void SxMolDyn::SxOutput::printParameters ()
{
   cout << "**************************" << endl;
   cout << "Output --- Parameters" << endl;
   cout << "-------------------------" << endl;
   cout << "saveWaves        : " << saveWaves << endl;
   cout << "saveHistory      : " << saveHistory << endl;
   if (printStepsFreq != -1) {
      cout << "printGFS         : " << endl;
      cout << "printSteps       : " << printStepsFreq << endl;
      cout << "inSections       : " << inSectionsFreq << endl;
   }
   if (printStepsHessian != -1) {
      cout << "printHessian         : " << endl;
      cout << "printSteps       : " << printStepsHessian << endl;
      cout << "inSections       : " << inSectionsHessian << endl;
   }
   cout << "**************************" << endl;
   cout << endl;
}

void SxMolDyn::SxOutput::writeRestartInfo (SxState *oldState, SxState *state, int it)
{
   
   SxString help;
   int is, ia, nAtoms, counter;
   //int dummyInt;
   //double dummyDouble;

   SxString SOPTACC;
   
   SOPTACC = SxString("%20.12f");
   FILE *fpRestart;

   if ( !(fpRestart = fopen ("restart.dat", "w")) ) {
      sxprintf ("Can't open file restart.dat");
      SX_EXIT;
   } 
   
   fprintf (fpRestart, "#step#     ");
   help = SxString (it, SOPTACC);
   fprintf (fpRestart, "%s\n", help.ascii ());

   fprintf (fpRestart, "#structureOld#\n");
   counter = 0.;
   for (is = 0; is < oldState -> X.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = oldState -> X.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (oldState -> X.ref(counter)(0), SOPTACC);
         SxString numberTwo   (oldState -> X.ref(counter)(1), SOPTACC);
         SxString numberThree (oldState -> X.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#velocityOld#\n");
   counter = 0;
   for (is = 0; is < oldState -> V.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = oldState -> V.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (oldState -> V.ref(counter)(0), SOPTACC);
         SxString numberTwo   (oldState -> V.ref(counter)(1), SOPTACC);
         SxString numberThree (oldState -> V.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#accelerationOld#\n");
   counter = 0;
   for (is = 0; is < oldState -> A.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = oldState -> A.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (oldState -> A.ref(counter)(0), SOPTACC);
         SxString numberTwo   (oldState -> A.ref(counter)(1), SOPTACC);
         SxString numberThree (oldState -> A.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#structure#\n");
   counter = 0.;
   for (is = 0; is < state -> X.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> X.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> X.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> X.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> X.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#velocity#\n");
   counter = 0;
   for (is = 0; is < state -> V.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> V.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> V.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> V.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> V.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#acceleration#\n");
   counter = 0;
   for (is = 0; is < state -> A.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> A.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> A.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> A.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> A.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#XRand#\n");
   counter = 0;
   for (is = 0; is < state -> XRand.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> XRand.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> XRand.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> XRand.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> XRand.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpRestart, "#XRandInv#\n");
   counter = 0;
   for (is = 0; is < state -> XRandInv.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> XRandInv.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> XRandInv.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> XRandInv.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> XRandInv.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpRestart, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   fprintf (fpRestart, "#\n");
   
   fclose(fpRestart);
}

void SxMolDyn::SxOutput::validateInput ()
{
   //empty
}

void SxMolDyn::SxOutput::writeEnergiesVsTime (SxState *state)
{
   //int dummyInt = 0;
   //double dummyDouble = 0.;
   if ( !(fpMoldyn = fopen ("moldyn.dat", "a")) ) {
      sxprintf ("Can't open file moldyn.dat");
      SX_EXIT;
   }
   
   fprintf (fpMoldyn, "%12.2f   %14.10f  %14.10f  %14.10f  "
                      "%14.10f  %14.10f  %14.6f  %14.6f\n"
         , state -> t
         , state -> ePot
         , state -> eKin
         , state -> eTot
         , state -> dUdL
         , state -> dUdLRef
         , (state -> dUdLRef) * HA2MEV
         , 3 * (state -> dUdLRef) * HA2MEV /
               double((state -> X).coordRef().getSize()));
   fclose(fpMoldyn);
}

//--- output file-format taken from old sphinx
//    suitable for DFT- (and tight-binding) calculations
//    for empirical potential calculations eventually too slow

void SxMolDyn::SxOutput::writeState (SxState *state)
{
   SxString help;
   int is, ia, nAtoms, counter; // dummyInt
   double dummyDouble = 0.;

   SxString SOPTACC;
   
   SOPTACC = SxString("%20.12f");
   
   if ( !(fpMoldynHist = fopen ("moldynHist.dat", "a")) ) {
      sxprintf ("Can't open file moldynHist.dat");
      SX_EXIT;
   }
   
   fprintf (fpMoldynHist, "@\n");
   
   fprintf (fpMoldynHist, "#step#     ");
   help = SxString (state->iSteps, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#time#     ");
   help = SxString (state->t, SOPTACC);
   fprintf (fpMoldynHist, "%s\n",help.ascii ());
   
   fprintf (fpMoldynHist, "#ePot#     ");
   help = SxString (state->ePot, SOPTACC);
   fprintf (fpMoldynHist, "%s\n",help.ascii ());
   
   if (parent -> guess.isActive) {
      fprintf (fpMoldynHist, "#dUdL#     ");
      help = SxString (state->dUdL, SOPTACC);
      fprintf (fpMoldynHist, "%s\n",help.ascii ());
   }
   
   fprintf (fpMoldynHist, "#eKin#     ");
   help = SxString (state->eKin, SOPTACC);
   fprintf (fpMoldynHist, "%s\n",help.ascii ());
   
   fprintf (fpMoldynHist, "#eTot#     ");
   help = SxString (state->eTot, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#T#        ");
   help = SxString (state -> T, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#fricSum#  ");
   help = SxString (state -> thermoX, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#friction# ");
   help = SxString (state -> thermoV, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#fricVel#  ");
   help = SxString (state -> thermoA, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#fricEPot# ");
   help = SxString ( dummyDouble, SOPTACC) ;
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#fricEKin# ");
   help = SxString (dummyDouble, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());
   
   fprintf (fpMoldynHist, "#HamNH#    ");
   help = SxString (state -> hamNH, SOPTACC);
   fprintf (fpMoldynHist, "%s\n", help.ascii ());

   fprintf (fpMoldynHist, "#structure#\n");
   counter = 0.;
   for (is = 0; is < state -> X.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> X.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> X.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> X.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> X.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpMoldynHist, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpMoldynHist, "#velocity#\n");
   counter = 0;
   for (is = 0; is < state -> V.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> V.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> V.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> V.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> V.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpMoldynHist, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fprintf (fpMoldynHist, "#force#\n");
   counter = 0.;

   for (is = 0; is < state -> V.getNSpecies (); is++) {
      SxString elem = parent -> speciesData.chemName(is);
      if (elem.getSize() == 1) elem = elem + SxString (" ");
      for (ia = 0, nAtoms = state -> V.getNAtoms(is); ia < nAtoms; ia++)  {
         SxString numberOne   (state -> F.ref(counter)(0), SOPTACC);
         SxString numberTwo   (state -> F.ref(counter)(1), SOPTACC);
         SxString numberThree (state -> F.ref(counter)(2), SOPTACC);
         SxString tosave = elem + SxString(" ") + numberOne + SxString(" ")
            + numberTwo + SxString(" ")	+ numberThree;
         fprintf (fpMoldynHist, "%s\n", tosave.ascii ());
         counter++;
      }
   }
   
   fclose (fpMoldynHist);
}

void SxMolDyn::SxOutput::printGFS ()
{
   SxString fileName;
   SxArray<double> G;
   double dfreq;
   double freq = 0.;
   int startIndex, endIndex, i;
   
   if (((parent -> state.iSteps) % printStepsFreq) == 0) {
      
      if (inSectionsFreq) {
         fileName = SxString("freq-");
         startIndex = (parent -> state.iSteps) - printStepsFreq;
         fileName += SxString(startIndex);
         fileName += SxString("-to-");
         endIndex = parent -> state.iSteps - 1;
         fileName += SxString((parent -> state.iSteps));
         fileName += SxString(".dat");
         parent -> timeCorrVel.setObservationRange(startIndex, endIndex);
         G = parent -> timeCorrVel.getGeneralizedFrequencySpectrum ();
         dfreq = parent -> timeCorrVel. getDFreq  ();
      } else {
         startIndex = 0;
         endIndex = (parent -> state.iSteps) - 1;
         parent -> timeCorrVel.setObservationRange(startIndex, endIndex);
         G = parent -> timeCorrVel.getGeneralizedFrequencySpectrum ();
         dfreq = parent -> timeCorrVel. getDFreq  ();
         fileName = SxString ("freq.dat");
      }
      
      if ( !(fpFreq = fopen (fileName.ascii(), "w")) ) {
         sxprintf ("Can't open file freq.dat");
         SX_EXIT;
      }
      
      for (i = 0; i < G.getSize (); i++) {
         if (freq < 4000.) {
            fprintf (fpFreq, "%10.5f %e\n", freq, G(i));
            freq += dfreq;
         }
      }
      
      fclose(fpFreq);
   }
}

void SxMolDyn::SxOutput::printHessian ()
{
   SxString fileName;
   int startIndex, endIndex, i, j;
   int nDoF = parent -> nDoF;
   SxBinIO io;
   parent -> timeCorrVel.setMasses (parent -> massVec);
   parent -> timeCorrForce.setMasses (parent -> massVec);
   
   if (((parent -> state.iSteps) % printStepsHessian) == 0) {
      
      if (inSectionsFreq) {
         fileName = SxString("hessian-");
         startIndex = (parent -> state.iSteps) - printStepsHessian;
         fileName += SxString(startIndex);
         fileName += SxString("-to-");
         endIndex = parent -> state.iSteps - 1;
         fileName += SxString((parent -> state.iSteps));
         fileName += SxString(".sx");
         parent -> timeCorrVel.setObservationRange(startIndex, endIndex);
         parent -> timeCorrForce.setObservationRange(startIndex, endIndex);
      } else {
         startIndex = 0;
         endIndex = (parent -> state.iSteps) - 1;
         parent -> timeCorrVel.setObservationRange(startIndex, endIndex);
         parent -> timeCorrForce.setObservationRange(startIndex, endIndex);
         fileName = SxString ("hessian.sx");
      }
     
      SxMatrix<Double> hessian = parent -> timeCorrVel.getHessian 
         (parent -> timeCorrVel, parent -> timeCorrForce);
      
      io.open (fileName, SxBinIO::SX_WRITE_ONLY);
      io.write ("matrix", hessian, "dummy", "dummy");
      io.close();

      SxMatrix<Double>::Eigensystem eig;
      eig = hessian.eigensystem ();
      io.open (fileName, SxBinIO::SX_APPEND_ONLY);
      io.write ("eigenvalues", eig.vals, "dummy", 0);
      io.close();

      for (i = 0; i < nDoF; i++) {
         for (j = 0; j < nDoF; j++) {
            hessian(i, j) = hessian(i, j)
               /sqrt(parent -> massVec(i)* parent -> massVec(j));
         }
      }
      eig = hessian.eigensystem ();
      
      SxVector<Complex16> freqs (nDoF);
      for (i = 0; i < nDoF; i++) 
         freqs(i) = sqrt(eig.vals(i).re) * 5123.75;
      
      io.open (fileName, SxBinIO::SX_APPEND_ONLY);
      io.write ("frequencies", freqs, "dummy", 0);
      io.close();
   }
}

SxMolDyn::SxOutput::~SxOutput()
{
   // empty
}




// --- Embedded Class SxInit

SxMolDyn::SxInit::SxInit()
{
   // empty
}

SxMolDyn::SxInit::SxInit (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   fromTemperature = false;
   if (cmd->containsGroup ("initHistory")) {
      SxSymbolTable *initI = cmd->getGroup ("initHistory");
      
      if (initI->contains ("file")) {
         scheme = SxString("file");
         initHistoryFN = initI -> get("file") -> toString ();
      }
      
      if (initI -> containsGroup ("randomVel")) {
         if (initI -> getGroup("randomVel") -> contains("seed"))
            seed = initI -> getGroup("randomVel") -> get ("seed") -> toInt ();
         else
            seed = 1;
         if (initI -> getGroup("randomVel") -> contains("noRandom"))
            noRandom = initI -> getGroup("randomVel") -> get ("noRandom") -> toBool ();
         else
            noRandom = false;

         scheme = SxString("randomVel");
         if (initI -> getGroup ("randomVel") 
               -> containsGroup("initStructure")) 
            initStructureFN = initI -> getGroup("randomVel")
               -> getGroup("initStructure") 
               -> get("file") -> toString();
         else initStructureFN = SxString("none");
         if (initI -> getGroup ("randomVel") 
               -> containsGroup("initHessian")) 
            initHessianFN = initI -> getGroup("randomVel") 
               ->  getGroup("initHessian") -> 
               get("file") -> toString();
         else initHessianFN = SxString("none");
         gas = initI -> getGroup("randomVel") -> get ("gas") -> toString ();
         if (initI -> getGroup("randomVel") -> contains("initEkin"))
            initEkin = initI -> getGroup("randomVel") 
               -> get ("initEkin") -> toReal ();
         if (initI -> getGroup("randomVel") -> contains("fromTemperature")) 
            fromTemperature = true;
      }
      
      if (initI -> containsGroup ("devAtoms")) {
         scheme = SxString("devAtoms");
         initStructureFN = initI -> getGroup ("devAtoms") 
            -> getGroup("initStructure") -> 
            get("file") -> toString();
         if (initI -> getGroup ("devAtoms") -> containsGroup("initHessian")) 
            initHessianFN = initI -> getGroup ("devAtoms")
               ->  getGroup("initHessian") -> 
               get("file") -> toString();
         else initHessianFN = SxString("none");
         gas = initI -> getGroup ("devAtoms") -> get ("gas") -> toString ();
         deltaE = initI -> getGroup ("devAtoms") -> get ("deltaE") -> toReal ();
         dof = initI -> getGroup ("devAtoms") -> get ("dof") -> toInt ();
      }
      
   } else {
      cout << "Please set up initHistory !" << endl;
   }
   
}

void SxMolDyn::SxInit::validateInput ()
{
   //empty
}

void SxMolDyn::SxInit::printParameters ()
{
   cout << "**************************" << endl;
   cout << "Init --- Parameters" << endl;
   cout << "-------------------------" << endl;
   if (scheme == SxString ("file")) {
      cout << "init by          : file" << endl;
      cout << "filename         : " << initHistoryFN << endl;
   }
   if (scheme == SxString ("randomVel")) {
      cout << "init by           : random velocities" << endl;
      cout << "initial structure : " << initStructureFN << endl;
      cout << "initial hessian   : " << initHessianFN << endl;
      cout << "gas               : " << gas << endl;
      cout << "initEkin          : " << initEkin << endl;
   }
   if (scheme == SxString ("devAtoms")) {
      cout << "init by           : deviate Atoms" << endl;
      cout << "initial structure : " << initStructureFN << endl;
      cout << "initial hessian   : " << initHessianFN << endl;
      cout << "gas               : " << gas << endl;
      cout << "deltaE            : " << deltaE << endl;
      cout << "dof               : " << dof << endl;
   }
   cout << "**************************" << endl;
   cout << endl;
}

SxMolDyn::SxState SxMolDyn::SxInit::getInitialState ()
{
   double norm;
   int nAtoms;
   SxAtomicStructure initialV;
   SxState initialState (parent);
   SxRandom random;
   
   cout << endl << endl << endl;
   cout << "--------------------------------------" << endl;
   cout << "Entering Molecular Dynamics initialisation ...." << endl;
   
   if (scheme == SxString ("file")) { 
      cout << "Initialising by file ...." << endl;
      cout << "not implemented yet !!" << endl;
      SX_EXIT;
   } else {
      initialState.T = 0.;	
      initialState.E0Ref = SX_HUGE;
      if (scheme == SxString("randomVel")) {
         cout << "... initialising with random velocities ...." << endl;
         initialState.X.copy (parent -> initialTau);
         initialState.V.copy (parent -> initialTau);
         SxVector<TPrecTauR> helpVec = initialState.V.coordRef();

         //--- INITIALISATION OF PSEUDO RANDOM GENERATOR
		 if (!noRandom) {
			 double theTime = SxTime::getRealTime ();
			 // reinterprete theTime as integer value
          unsigned x = 0;
			 // make sure that all bits in theTime are used
			 for (unsigned i = 0; i < sizeof(double) / sizeof(unsigned); ++i)
				 x += reinterpret_cast<const unsigned *>(&theTime)[i];
			 srand(x * seed);
		 }
         
         if (fromTemperature) initEkin = parent->thermostat.getInitialEkin ();
         
         if (gas == SxString("ideal")) {
            cout << "... initialising in cartesian basis ...." << endl;
            
            for (int i = 0; i < parent -> nDoF; i++) {
               helpVec (i) = (0.5 - rand()/double(RAND_MAX))
               	/sqrt(parent -> massVec(i));
               //helpVec (i) = (0.5 - random.get ())
               //   /sqrt(parent -> massVec(i));
               
            }
            initialState.V.set(helpVec);
         }
         
         //--- initialisation in phonon basis, shall provide shorter 
         //    equilibration times
         //    is in testing phase
         //    TODO: more extended docu in case that it helps
         
         if (gas == SxString("phonon")) {
            cout << "Initialising in phonon basis ...." << endl;
            if (initHessianFN == "none") {
               cout << "Please provide an appropriate "
                  << "Filename for the Hessian Matrix" << endl;
            } else {
               SX_EXIT; // CF, 2014/08/01: seems to be a fragment
               /*
               SxBinIO io;
               SxMatrix<Double> initHessian 
                  (parent -> nDoF, parent -> nDoF);
               SxMatrix<Double>::Eigensystem eig;
               io.open (initHessianFN, SxBinIO::ASCII_READ_ONLY);
               io.read ("dummy", &initHessian, 
                     parent -> nDoF, parent -> nDoF, 0, 0);
               io.close();
               
               SxHessianOps hOps (initHessian, parent -> massVec);
               SxMatrix<Double> eigenVelocities = hOps.getEigenVelocities ();
               SxVector<Double> helpVec (parent -> nDoF);
               SxVector<Double> toAdd (parent -> nDoF);
               
               helpVec.set (0.);
               for (int i = 0; i < parent -> nDoF; i++) {
                  toAdd.set (eigenVelocities.colRef(i));
               toAdd *= (0.5 - rand()/double(RAND_MAX));
               //toAdd *= 	(0.5 - random.get ()/(double) 0x7fffffff);
               helpVec += toAdd;
               }
               */
            }
         }
         
         SxAtomicStructure P;
         P.copy (initialState.V);
         P.set  (P.coordRef () *parent -> massVec);
         P.set 
            ((parent -> velocityFilter | P.coordRef ()));
         
         SxVector3<Double> trans = P.sum ();
         nAtoms = initialState.X.getNAtoms ();
         trans = trans/nAtoms;
         for (int i = 0; i < nAtoms; i++) P.ref(i) = P(i) - trans;
         
         initialState.V.set (P.coordRef () / parent -> massVec);
         norm = parent -> getEKin (initialState.V);
         initialState.V = sqrt(initEkin/norm)*initialState.V;

         // in the first step we need a lower energy criterion
         // for an accurate reference energy
         parent->potential->dEnergyLow = true;
         initialState.F = parent -> getForces(initialState.X);
         parent->potential->dEnergyLow = false;
         
         initialState.A.set (initialState.F.coordRef () / parent -> massVec); 
         
      }
      
      if (scheme == ("devAtoms")) { 
         cout << "Initialising by displacing atoms ...." << endl;
         cout << "not implemented yet !!" << endl;
         SX_EXIT;
      }
   }
   
   initialState.XRand.copy (initialState.X);
   initialState.XRandInv.copy (initialState.X);
   initialState.XRand.coordRef ().set (0.);
   initialState.XRandInv.coordRef ().set (0.);
   return initialState;
}	

SxMolDyn::SxInit::~SxInit()
{
   // empty
}


// --- Embedded Class SxState

SxMolDyn::SxState::SxState()
{
   // empty
}

SxMolDyn::SxState::SxState (SxMolDyn *p)
{
   parent = p;
   // Just for allocation (a bit ugly)
   X.copy (parent -> initialTau);
   XRand.copy (parent -> initialTau);
   XRandInv.copy (parent -> initialTau);
   V.copy (parent -> initialTau);
   A.copy (parent -> initialTau);
   F.copy (parent -> initialTau);
   thermoX = thermoV = thermoA = thermoF = T = 0.;
   
   iSteps = -1;
   t = 0.;
}

SxMolDyn::SxState::~SxState()
{
   // empty
}

void SxMolDyn::SxState::set  (const SxState &in)
{
   ePot = in.ePot;
   eKin = in.eKin;
   eTot = in.eTot;
   thermoX = in.thermoX;
   thermoV = in.thermoV;
   thermoA = in.thermoA;
   thermoF = in.thermoF;
   T = in.T;
   X.copy (in.X);
   V.copy (in.V);
   A.copy (in.A);
   F.copy (in.F);

   XRand.copy (in.XRand);
   XRandInv.copy (in.XRandInv);
}

void SxMolDyn::SxState::updateVariables ()
{
   ePot = parent -> getEPot ();
   eKin = parent -> getEKin (V);
   eTot = ePot + eKin;
   t += parent -> integrator.dt;
   if (parent -> integrator.scheme == SxString("NH")) {
      hamNH = eTot + thermoV*thermoV*(parent -> thermostat.mass)
         + (parent -> thermostat.L)*KB*T*JL2HA;
   }
   if (parent -> integrator.scheme == SxString("none")) {
      hamNH = eTot;
   }
   
   iSteps += 1;
}

void SxMolDyn::SxState::plotDihedrals ()
{
   if (parent -> output.plotDih) {
   
   int C_one = parent -> output.C_one - 1;
   int C_two = parent -> output.C_two - 1;
   int N_one = parent -> output.N_one - 1;
   int N_two = parent -> output.N_two - 1;
   int C_alpha = parent -> output.C_alpha - 1;

   
   SxVector3<Double> b1 = X(N_one) - X(C_one);
   SxVector3<Double> b2 = X(C_alpha) - X(N_one);
   SxVector3<Double> b3 = X(C_two) - X(C_alpha);

   SxVector3<Double> b23;
   SxVector3<Double> b12;

   b23(0) = b2(1)*b3(2) - b2(2)*b3(1);
   b23(1) = b2(2)*b3(0) - b2(0)*b3(2);
   b23(2) = b2(0)*b3(1) - b2(1)*b3(0);

   b12(0) = b1(1)*b2(2) - b1(2)*b2(1);
   b12(1) = b1(2)*b2(0) - b1(0)*b2(2);
   b12(2) = b1(0)*b2(1) - b1(1)*b2(0);

   double b = b12(0)*b23(0) + b12(1)*b23(1) + b12(2)*b23(2);
   double a = sqrt(b2.absSqr ().sum ()) * (b1(0)*b23(0) + b1(1)*b23(1) + b1(2)*b23(2));
   
   double phi = atan2(a, b);
   phi = phi/2./PI*360.;
   
   //if (phi < -150) phi = phi + 360.; 

   cout << "PHI: " << phi << endl;
   
    b1 = X(C_alpha) - X(N_one);
    b2 = X(C_two) - X(C_alpha);
    b3 = X(N_two) - X(C_two);

  
   b23(0) = b2(1)*b3(2) - b2(2)*b3(1);
   b23(1) = b2(2)*b3(0) - b2(0)*b3(2);
   b23(2) = b2(0)*b3(1) - b2(1)*b3(0);

   b12(0) = b1(1)*b2(2) - b1(2)*b2(1);
   b12(1) = b1(2)*b2(0) - b1(0)*b2(2);
   b12(2) = b1(0)*b2(1) - b1(1)*b2(0);

    b = b12(0)*b23(0) + b12(1)*b23(1) + b12(2)*b23(2);
    a = sqrt(b2.absSqr ().sum ()) * (b1(0)*b23(0) + b1(1)*b23(1) + b1(2)*b23(2));
   
   double psi = atan2(a, b);
   psi = psi/2./PI*360.;

   //if (psi < -150) psi = psi + 360.; 

   cout << "PSI: " << psi << endl;
   }
}

//--- Embedded Class SxGuess
SxMolDyn::SxGuess::SxGuess ()
{
   // empty
}

SxMolDyn::SxGuess::SxGuess (const SxSymbolTable *cmd, SxMolDyn *p)
{
   parent = p;
   potential = NULL;
   
   
   if (cmd -> containsGroup ("guessedPotential") ) {
      isActive = true;
      SxSymbolTable *potGroup = cmd -> getGroup ("guessedPotential"); 
      if (potGroup->contains ("lambda")) 
         lambda = potGroup->get("lambda")->toReal ();
      else 
         lambda = 1.;
     
      //--- access on SxThermostat is a bit ugly;
      (parent -> thermostat).storeLambda = lambda;
         
      if (potGroup->containsGroup ("taylorExpPotential")) {
         cout << "before create" << endl;
         //SxAutoPointer<SxTaylorExpPotential> tPotential;
         parent->tPotential = SxPtr<SxTaylorExpPotential>::create (
               potGroup->getGroup("taylorExpPotential")
               );
         //--only for test purposes
         parent->tPotential2 = SxPtr<SxTaylorExpPotential>::create (
               potGroup->getGroup("taylorExpPotential")
               );

         

         SxArtifactFilter hessianFilter;
         
         hessianFilter.
            set (parent -> tPotential -> equTau, parent->performance.freezeMode, false);
         
         (parent -> tPotential -> hessian) = (hessianFilter | parent -> tPotential -> hessian);
      
         
         SxMatrix<Double>::Eigensystem eig;
         eig =  (parent -> tPotential -> hessian).eigensystem ();

         //--- only for test purposes
         (parent -> tPotential2 -> hessian).copy (parent -> tPotential -> hessian);
         
         (parent -> tPotential2 -> hessian) *= 0.8;



         //hessianFilter.
         //   set (tPotential -> equTau, parent -> performance.freezeMode, false);
         
        // (tPotential -> hessian) = (hessianFilter | tPotential -> hessian);

         //potential = &*tPotential;
      } else {
         cout << "No valid guess-potential found in group guess" << endl;
         SX_EXIT;
      }
   } else {
      isActive = false;
   }
}


SxMolDyn::SxGuess::~SxGuess ()
{
    // if (potential != NULL) delete potential;
}
// --- Outer Class SxMolDyn

SxMolDyn::SxMolDyn (const SxAtomicStructure &tauIn,  SxPotential *pot) 
   : initialTau (tauIn),
potential (pot)   
{
   int i, ia, is, counter;
   
   hamSolver = dynamic_cast<SxHamSolver *> (pot);
   speciesData = potential->getSpeciesData ();
   
   useTimeCorr = false;
   
   nDoF = 3*initialTau.getNAtoms ();
   // SIXTEN: masses need to implemented (in SxAtomicStructure ??) 
   massVec.resize (nDoF);
   
   // a bit ugly ...
   counter = 0.;
   for (is = 0; is < initialTau.getNSpecies (); is++) {
      for (ia = 0; ia < initialTau.getNAtoms (is); ia++) {
         for (i = 0; i < 3; i++) {
            massVec(counter) = speciesData.ionicMass(is);
            counter++;
         }
      }
   }
}

SxMolDyn::SxMolDyn ()
{
   // empty
}

SxMolDyn::~SxMolDyn ()
{
   // empty
}

//----------------------------------------------------------------------------
//    Interface to input file
//----------------------------------------------------------------------------
bool SxMolDyn::isRegistered (const SxSymbolTable *cmd)
{
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   return ( str == "molDyn" );
}


double SxMolDyn::getEPot ()
{
   double lambda = 0.;
   double ePot = 0.;
   if (guess.isActive) { 
      SxArray<const SxSymbolTable *> dummy (1);
      dummy(0) = NULL;
      lambda = guess.lambda;
      ePot =  (1. - lambda)*(tPotential -> getEnergy ());
   }
   if (fabs(lambda) > 1e-10) {
      ePot = ePot + lambda*(potential -> getPotentialEnergy ());
   }
   if (guess.isActive) {
      double EFull = potential -> getPotentialEnergy ();
      double EHarmRef = tPotential -> getEnergy ();
      state.dUdL = EFull - EHarmRef;
      if ( state.E0Ref == SX_HUGE ) state.E0Ref = EFull;
      state.dUdLRef = EFull - EHarmRef - state.E0Ref;
   } 
   return ePot;
}


double SxMolDyn::getEKin (SxAtomicStructure vel) 
{
   return (0.5*(vel.coordRef ().absSqr ()*massVec).sum ());
}


void SxMolDyn::print (const SxSymbolTable *cmd)
{
   execute (cmd, false);
}


int SxMolDyn::getFaculty (int number) 
{
   if (number == 0) return 1;
   else return number*getFaculty(number - 1);
}


double SxMolDyn::gaussDev (double sigmaSqr) 
{

 double fac,rsq,v1,v2, rand1, rand2;
 //SxRandom random;

 double sigma = sqrt (sigmaSqr);
 

 do {
    rand1 = rand()/double(RAND_MAX);
    rand2 = rand()/double(RAND_MAX);
   
   // v1 = 2.0*random.get() - 1.0;
   // v2 = 2.0*random.get() - 1.0;
    v1 = 2.0*rand1 - 1.0;
    v2 = 2.0*rand2 - 1.0;
    rsq = v1*v1 + v2*v2;
 } while (rsq >= 1.0 || rsq ==0.0);
 fac = sqrt(-2.0*log(rsq)/rsq);
 return (v2*fac*sigma);
}
 



void SxMolDyn::execute (const SxSymbolTable *cmd, bool calc)
{  
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   int i;
   
   if (str == "molDyn")  {
      cout << SX_SEPARATOR;
      cout << "| Molecular Dynamics\n";
      cout << SX_SEPARATOR;
      cout << "| Parameters: ...\n";
   }  else  {
      cout << "ERROR: Unknown command " << str << endl;
      SX_EXIT;
   }

         

   //--- Parsing in input-file
   integrator        = SxIntegrator   (cmd, this);
   init              = SxInit         (cmd, this);
   performance       = SxPerformance  (cmd, this);
   output            = SxOutput       (cmd, this);
   perturbation      = SxPerturbation (cmd, this);
   thermostat        = SxThermostat   (cmd, this);
   guess             = SxGuess        (cmd, this);
   // impSampler        = SxImpSampler   (cmd, this);
   
   elMinimCmds = potential->getMinimCmds (cmd);
   cout << "elMinimCmds " << elMinimCmds << endl;
   /*
   if (elMinimCmds.getSize() == 0)  {
      sxprintf ("No valid command found in group 'bornOppenheimer'\n");
      SX_EXIT;
   }
   */
   
   
   
   //--- Printing parameters for contol purposes
   init.printParameters ();
   integrator.printParameters ();
   thermostat.printParameters ();
   perturbation.printParameters ();
   performance.printParameters ();
   output.printParameters ();
   
   if (!calc)  return;
   
   //--- initialisation of system state
   
   performance.initialize ();
   thermostat.initialize ();
   state = SxState (this);
   state = init.getInitialState ();
   
   state.updateVariables ();
   state.t -= integrator.dt;
   //output.writeEnergiesVsTime (&state);
   oldState = SxState (this);
   oldState.set (state);

   if ((output.saveHistory) && (!(performance.restart))) 
      output.writeState (&state);
   integrator.initialize (&state, &oldState);
   performance.initialize ();
   // impSampler.initialize ();
   SxTimer timer;   
   for (thermostat.switchOn (); !(thermostat.switchedOff); thermostat.goOn ()) {
      for (TDInt.switchOn (); !(TDInt.switchedOff); TDInt.goOn ()) {
            
         timer.init(1);
         timer.start(0);
         for (i = 0; i < integrator.timeSteps; i++) {
            i = performance.meditateRestart(i, &oldState, &state);
            mdStep=i;
            
            integrator.predictNewState (&state, &oldState);

            state.F = getForces (state.X);
            thermostat.apply (&state);
            //--- workaround for equilibration (start)
            if (thermostat.iHarmEqu < thermostat.nHarmEqu) i = -1;
            //--- workaround for equlibration  (end);
            integrator.correctNewState (&state, &oldState);
      
            performance.update (state);
            state.updateVariables ();
            //TDInt.takeMeasurement (state);
            //TDInt.apply (&guess);
           // impSampler.inc ();
           // impSampler.apply ();
            if (useTimeCorr) {
               timeCorrVel.push (state.V, state.t);
               timeCorrForce.push (state.F, state.t);
            }

            state.plotDihedrals ();
            // only for test purposes, comment out for production
            //output.writeState (&state);
 
            if (i > -1 ) {
               if (((state.iSteps) % output.printSteps) == 0) {
                  //output.writeState (&state);
                  output.writeEnergiesVsTime (&state);
                  //output.printTDInt ();
                  if (output.saveHistory) output.writeState (&state);
                  if (output.saveRestartInfo) 
                     output.writeRestartInfo (&oldState, &state, i);
                  if (output.printStepsFreq != -1) output.printGFS ();
                  if (output.printStepsHessian != -1) output.printHessian ();
               }
            }
         }
         timer.stop(0);
         cout << "time: " << timer.getTime(0) << endl;
      }
   }
}


SxAtomicStructure SxMolDyn::getForces (const SxAtomicStructure &tau)  
{
   SxAtomicStructure forces;
   double lambda = 1.;
   if (guess.isActive) { 
      lambda = guess.lambda;
      if (lambda < 1.) {
         forces = tPotential -> getForces(tau);
         forces = (1. - lambda)*forces;
      }
   } 
   if (fabs(lambda) > 1e-10) {
      if ((hamSolver != NULL) && (performance.extrapolateWaves))
         performance.extrapolateWavefunction (tau);
      if (lambda < 1.) {
         if (elMinimCmds.getSize () == 0) {
            forces += lambda * (potential->getForces (tau));
         } else 
            forces += lambda * (potential->getForces (tau, elMinimCmds));
      } else {
         if (elMinimCmds.getSize () == 0) {
            forces.copy (potential->getForces (tau));
         } else 
            forces.copy (potential->getForces (tau, elMinimCmds));
      }
   } 
   forces.set ( (forceFilter | forces.coordRef ()));
   return forces;
}

