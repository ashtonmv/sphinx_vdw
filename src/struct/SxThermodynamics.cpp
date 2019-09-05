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

#include <SxThermodynamics.h>

SxThermodynamics::SxThermodynamics () 
{
   //empty
}

SxThermodynamics::SxThermodynamics 
(const SxVector<Complex16> &inputFreqs, int inDoF, 
     int inExtDoF, double inThreshold)
{
   //--- setting up frequency spectrum, threshold, norm and degrees 
   //    of freedom for a molecule or cluster
   setFrequencies (inputFreqs);
   nDoF = inDoF;
   setNExternalDoF(inExtDoF);
   setThreshold(inThreshold);
   setNorm (1.);

   if (nDoF != freqs.getSize ()) {
      cout << "SxThermodynamics:: error in frequency-spectrum. "
            <<   "Size does not correspond to nDoF" << endl;
         SX_EXIT;
      }
   // --- defaults for grueneisen-parameters;
   B0 = B0Prime = V0 = 0.;
   grueneisen =  false;
}

   
   SxThermodynamics::SxThermodynamics 
    (SxSecondaryStructure &peptideChain, const SxString &interpolation, 
     double inThreshold)
{
   //--- setting up frequency spectrum, threshold, norm and degrees 
   //    of freedom for a helix or linear chain
   
   double res = peptideChain.getPDRes ();
   nDoF = peptideChain.getNAtomsPP () * 3 *(int) (360.0/res + .5);
   setFrequencies 
      (peptideChain.getFrequencies (interpolation, res)) ;
   setNExternalDoF(4);
   setThreshold(inThreshold);
   setNorm 
      ((double) (int) (360.0/res + .5) /(double) peptideChain.getNPeptides ());
   
   if (nDoF != freqs.getSize ()) {
      cout << "SxThermodynamics:: error in frequency-spectrum. "
         <<   "Size does not correspond to nDoF" << endl;
              SX_EXIT;
   }

   // --- defaults for grueneisen-parameters;
   B0 = B0Prime = V0 = 0.;
   grueneisen =  false;
}

   SxThermodynamics::SxThermodynamics 
    (SxHessianOps &hessianOps, int inExtDoF, double inThreshold)
{
   //--- setting up frequency spectrum, threshold, norm and degrees 
   //    of freedom for a given input hessianOps object
   nDoF = hessianOps.getNDoF ();
   setFrequencies 
      (hessianOps.getEigenFrequencies ()) ;
   setNExternalDoF(inExtDoF);
   setThreshold(inThreshold);
   setNorm (1.);
   
   if (nDoF != freqs.getSize ()) {
      cout << "SxThermodynamics:: error in frequency-spectrum. "
         <<   "Size does not correspond to nDoF" << endl;
      SX_EXIT;
   }
   // --- defaults for grueneisen-parameters;
   B0 = B0Prime = V0 = 0.;
   grueneisen =  false;
}

SxThermodynamics::~SxThermodynamics ()
{
   //empty
}

void SxThermodynamics::setNorm (double normIn)
{
   norm = normIn;
}

void SxThermodynamics::setFrequencies (
      const SxVector<Complex16> &inputFrequencies)
{
   
   freqs.resize (inputFrequencies.getSize ());
   freqs.set (inputFrequencies);
   perm.resize (freqs.getSize ());
   
   // --- in order to provide correct calculation of thermodynamic properties
   //     the frequency-spectrum needs to be ordered with respect to 
   //     frequency value
   //     03/12/06 ... however, the according permutation needs to be stored
   //     in order to avoid inconsistencies with the grueneisen-parameters
   //     ... is stored in vector perm 
   bool changed = true;
   SxComplex16 help = 0.;
   int i, p;
   for (i = 0; i < nDoF; i++) perm (i) = i;
   while (changed) {
      changed = false;
      for (i = 0; i < (nDoF - 1); i++) {
         if (freqs(i).re > freqs(i + 1).re) {
            changed = true;
            help = freqs(i + 1);
            freqs(i + 1) = freqs (i);
            freqs(i) = help;
            p = perm(i+1);
            perm(i + 1) = perm (i);
            perm(i) = p;
         }
      }
   }
            
}

void SxThermodynamics::setGrueneisen ( const SxSymbolTable *inGroup )
{
   
   
   V0 = inGroup -> get("V0") -> toReal ();
   if (V0 < 1e-10) grueneisen = false;
   else {
      SxString fileName = inGroup -> get ("file") -> toString ();
      
      // --- loads in the grueneisen-coeff from a file
      FILE *fp = fopen (fileName.ascii (), "r");
      SX_CHECK (fp);
      char *bs = NULL;
      const int BUFLEN = 1024000;
      char buffer[BUFLEN];
      double f;
      SxList<double> ge;
      bool endOfFile = false;
	
      while (!endOfFile) {
         bs = fgets (buffer, 40, fp);
         SxString tk(bs);
         if (tk.getSize() > 1 ) {
            f = tk.toDouble();
            ge.append (f);
         } else 
            endOfFile = true;
      }
      fclose (fp);
   
      SxVector<Int> GECoeffUnsorted = SxVector<Double> (ge);
   
      if (nDoF != GECoeffUnsorted.getSize () ) {		
      cout << "Inconsistency for grueneisen parameter set" << endl;
      cout << GECoeffUnsorted.getSize () << " " << nDoF << endl;
      SX_EXIT;
      }

   

      //--- grueneisen-parameters need to be sorted consistent with 
      //    the frequency spectrum
      GECoeff = SxVector<Double> (nDoF);
      for (int i = 0; i < nDoF; i++) 
         GECoeff(perm(i)) = GECoeffUnsorted(i);
      // GECoeff(i) = GECoeffUnsorted(perm(i));

      B0 = inGroup -> get("B0") -> toReal ();
      B0Prime = inGroup -> get("B0Prime") -> toReal ();

      // --- Unit conversion : input of B0 is 
      //     expected to be in GPa 
      B0 = B0 * JL2KCALPM*1e9*(::pow (BOHRRADIUS, 3));
      //B0Prime = B0 * HA2KCALPM;
   }
}
   
void SxThermodynamics::setNExternalDoF (int inExtDoF)
{
   externalDoF = inExtDoF;
}

void SxThermodynamics::setThreshold (double inThreshold)
{
   threshold = inThreshold;
}

double SxThermodynamics::getSVib (double Vol, double Temp)
{
   int i = 0;
   double TS = 0;
   double fac1 = CM2KCALPM; 
   double fac3 = JL2HA*HA2KCALPM*KB;           
   double fac2 = CVEL*100.* HPLANCK/KB; 
   double w0, w;
   
   if (Temp == 0.) return 0.;
   else {
      for (i = externalDoF; i < nDoF; i++) { 
         if (freqs(i).re > threshold) {
            w0 = (double) freqs(i);
            if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            else w = w0;
            
            TS += 
               fac1*w*(1./(exp(fac2*w/Temp) - 1.))
               - Temp*fac3*log(1. - exp(-fac2*w/Temp));
         }
      }
      return ((TS/Temp)/norm);
   }
}

double SxThermodynamics::getSVibCl (double Vol, double Temp)
{
   int i = 0;
   double TS = 0;
   double fac3 = JL2HA*HA2KCALPM*KB;           
   double fac2 = CVEL*100.* HPLANCK/KB; 
   double w0, w;
   
   if (Temp == 0.) return 0.;
   else {
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) {
            w0 = (double) freqs(i);
            if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            else w = w0;
            
            TS += fac3*(1. -  log(fac2*w/Temp));
         }
      }
   }
   return (TS/norm);
}

double SxThermodynamics::getCVVib (double Vol, double Temp)
{
   int i = 0;
   double fac2 = CVEL*100.* HPLANCK/KB;
   double C, exponential; C = exponential = 0.;
   double w0, w;
   
   if (Temp == 0.) return 0.;
   else {
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) {
            w0 = (double) freqs(i);
            if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            else w = w0;
            
            exponential =  exp (-fac2*(double) w /Temp);
            C += JL2HA*HA2KCALPM*KB*::pow((fac2*(double) w /Temp), 2.) 
               * exponential/::pow((1. - exponential), 2.);
         
         }
      }
      return (C/norm);
   }
}

double SxThermodynamics::getCVVibCl (double /*Vol*/, double /*Temp*/)
{
   int i = 0;
   double fac = JL2HA*HA2KCALPM*KB;          
   double C, exponential; C = exponential = 0.;
   
   for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) 
            C += 1.;
   }
      return (C*fac/norm);
}

double SxThermodynamics::getFVib (double Vol, double Temp)
{
   return (getUVib(Vol, Temp) - Temp*getSVib(Vol, Temp));
}

double SxThermodynamics::getFVibCl (double Vol, double Temp)
{
   return (getUVibCl(Vol, Temp) - Temp*getSVibCl(Vol, Temp));
}

double SxThermodynamics::getUVib (double Vol, double Temp)
{
   int i;
   double fac2 = CVEL*100.* HPLANCK/KB;
   double U = 0.;
   double w, w0;
   
   if (Temp == 0.) return getZPV(Vol);
   else {
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) {
            w0 = (double) freqs(i);
            if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            else w = w0;
            U += (double)w/2.
               * (1. + 2./( exp (fac2*(double) w /Temp) - 1.));
         }
      }
      return (U*CM2KCALPM/norm);
   }
}

double SxThermodynamics::getUVibCl (double /*Vol*/, double Temp)
{
   int i;
   double fac = JL2HA*HA2KCALPM*KB;           
   double U = 0.;
   
   if (Temp == 0.) return 0.;
   else {
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) 
            U += Temp;
      }
      return (U*fac/norm);
   }
}

double SxThermodynamics::getZPV (double Vol)
{
   int i;
   double zerocorr = 0.;
   double w, w0;
   
   for (i = externalDoF; i < nDoF; i++) {
      if (freqs(i).re > threshold) {
         w0 = (double) freqs(i);
         if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
         else w = w0;
         zerocorr += w/2.;
      }
   }
   return (zerocorr*CM2KCALPM/norm);
}   

double SxThermodynamics::getRescaledTemp(double Vol, double Temp) 
{
   double Trescaled = 0;
   int i;
   double fac1 = CM2KCALPM; 
   double fac3 = JL2KCALPM*KB;           
   double fac2 = CVEL*100.* HPLANCK/KB; 
   double w, w0;
   
   if (Temp == 0.) {
      Trescaled = 
         getZPV(Vol)/(fac3*(double)(nDoF - externalDoF))*norm;
      return Trescaled;
   } else {
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) { 
            w0 = (double) freqs(i);
            if (grueneisen) w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            else w = w0;
            Trescaled += (w*
                  (0.5 + 1./(exp(fac2*w/Temp) - 1.)));
         }
      }
      return (Trescaled*fac1/fac3/(double)(nDoF - externalDoF));
   }
}

double SxThermodynamics::getEMurn (double Vol)
{
   double D = B0*Vol/(B0Prime*(B0Prime - 1.));
   double C = B0Prime*(1. - V0/Vol) + ::pow (V0/Vol, B0Prime) - 1.;
   return (D*C);
}

double SxThermodynamics::getdEdVMurn (double Vol)
{
   double returnValue;
   double D = B0/(B0Prime*(B0Prime - 1.));
   double C = B0Prime*(1. - V0/Vol) + ::pow (V0/Vol, B0Prime) - 1.;
   returnValue = D*C + D*Vol 
                     * (B0Prime*V0/Vol/Vol - B0Prime*::pow (V0/Vol, B0Prime - 1.)
                                       * V0/Vol/Vol);

   
   //--- comparison to numerical derivative
  // double toCompare = (getEMurn(Vol + 1e-5) - getEMurn(Vol - 1e-5))/2e-5;
  // cout << "dEdV: " << returnValue << " " << toCompare;
   return returnValue;
  // SX_EXIT;
}   

double SxThermodynamics::getdFvibdV (double Vol, double Temp)
{
   int i;
   double fac2 = CVEL*100.* HPLANCK/KB;
   double U = 0.;
   double w, w0;
   
      for (i = externalDoF; i < nDoF; i++) {
         if (freqs(i).re > threshold) {
            w0 = (double) freqs(i);
            w = w0*(1. - GECoeff(i)*(Vol - V0)/V0);
            //w = w0;
            if (Temp > 0.)
               U += (double)w/2.
                  * (1. + 2./( exp (fac2*(double) w /Temp) - 1.))
                  * GECoeff(i);
            else U+= (double)w/2.*GECoeff(i);
         }
      }
      return (-U*CM2KCALPM/norm/Vol);
}
  
   
double SxThermodynamics::getdFdV (double Vol, double Temp)
{
   return (getdEdVMurn (Vol) + getdFvibdV (Vol, Temp));
}

double SxThermodynamics::solveEquOfState (double p, double Temp)
{
   double V1, V2/*, V3*/;
   double dF1, dF2 /*, dF3, dFGuess */;
   double Vopt, dFopt, VoptOld/*, Vopt1, Vopt2*/;
   double l1, l2/*, l3*/;
   //double A;
   /* ############## Parabolic solver not working yet 
   // ---- Initialising parameters 
   Vopt = 1.;
   VoptOld = 0.5;
   V1 = 2*V0;
   V2 = 2*V0 + V0/1e4;
   V3 = 2*V0 - V0/1e4;
   dF1 = getdFdV (V1, Temp);
   dF2 = getdFdV (V2, Temp);
   dF3 = getdFdV (V3, Temp);

   // ---- loop to solve the equation of state p = -dF/dV
   while (fabs(Vopt - VoptOld) > 1e-7) {
      VoptOld = Vopt;
      //--- parabolic fit:
      //    determination of coefficients for parabola
      //    dF(V) = l2*V*V + l1*V + l3;
      A = (V2*V2 - V1*V1)/(V3*V3 - V1*V1);
      l1 = (A*(dF3 - dF1) - (dF2 - dF1)) / (A*(V3 - V1) - (V2 - V1));
      l2 = ((dF2 - dF1) - l1*(V2 - V1)) / (V2*V2 - V1*V1);
      l3 = dF1 - V1*l1 - V1*V1*l2;
      // check coefficients 
      //cout << dF1 << " " << (l2*V1*V1 + l1*V1 + l3) << endl;
      //cout << dF2 << " " << (l2*V2*V2 + l1*V2 + l3) << endl;
      //cout << dF3 << " " << (l2*V3*V3 + l1*V3 + l3) << endl;
   
      //--- solving equ. of state for parabola
      
      if ( ((l1*l1/4./l2/l2) > (l3 + p)/l2) ) {
         Vopt1 = - (l1/2./l2) - sqrt ((l1*l1/4./l2/l2) - (l3 + p)/l2);
         Vopt2 = - (l1/2./l2) + sqrt ((l1*l1/4./l2/l2) - (l3 + p)/l2);
      
         if (getEMurn (Vopt1) <= getEMurn(Vopt2)) Vopt = Vopt1;
         else Vopt = Vopt2;
      } else Vopt = -l1/l2;
      dFGuess = l2*Vopt*Vopt + l1*Vopt + l3;  
      
      //--- getting new derivative at new guess Vopt
      dFopt = getdFdV (Vopt, Temp);

      cout << "Guess: Vopt = " << Vopt << "; dFGuess = "
           << dFGuess << "; dF = " << dFopt << endl;

      //--- removing worst guess
      if ( (fabs (dF1 + p) >= fabs (dF2 + p)) &&
           (fabs (dF1 + p) >= fabs (dF3 + p)) ) {
         V1 = Vopt; dF1 = dFopt;
      } else {
         if ( (fabs (dF2 + p) >= fabs (dF1 + p)) &&
               (fabs (dF2 + p) >= fabs (dF3 + p)) ) {
            V2 = Vopt; dF2 = dFopt;
         } else 
            V3 = Vopt; dF3 = dFopt;
      }
   }
   
   */
   
   // ---- Initialising parameters for linear fit solver

   if (!grueneisen) return V0;

   Vopt = 1.;
   VoptOld = 0.5;
   V1 = V0;
   V2 = V0 + V0/1e4;
   dF1 = getdFdV (V1, Temp);
   dF2 = getdFdV (V2, Temp);
   int nSteps = 0;

   // ---- loop to solve the equation of state p = -dF/dV
   while (fabs(Vopt - VoptOld) > 1e-7) {
      nSteps ++;
      VoptOld = Vopt;
      //--- linear fit:
      //    determination of coefficients  
      //    dF(V)/dV = l2*V + l1;
      l2 = (dF2 - dF1) / (V2 - V1);
      l1 = dF1 - l2*V1;
      // solving equ. of state for linear-fit 
      Vopt = - (p + l1)/l2;
      
      //dFGuess = l2*Vopt + l1;  
      
      //--- getting new derivative at new guess Vopt
      dFopt = getdFdV (Vopt, Temp);
    //  cout << "Guess: Vopt = " << Vopt << "; dFGuess = "
    //       << dFGuess << "; dF = " << dFopt << endl;
    //--- removing worst guess
      if ( (fabs (dF1 + p) >= fabs (dF2 + p)) ) {
         V1 = Vopt; dF1 = dFopt;
      } else {
         V2 = Vopt; dF2 = dFopt;
      } 
   }
   /*
   cout << "Equ. of state: "
        << "T = " << Temp << "; V = " << Vopt << "; p = " << p
        << "; dF/dV = "<< dFopt << "; nSteps = " << nSteps << endl;
   */
//sxprintf ("Equ. of state: T = %5.1f ; V = %5.4f ; p = %5.1f ; dF/dV = %f ;  nSteps = %d \n", Temp, Vopt, p, dFopt, nSteps);
  
   sxprintf ("Equ. of state: T = %5.1f ; V = %7.6f ; p = %5.1f ; dF/dV = %f \n", Temp, Vopt, p, dFopt);
   fflush (stdout);
   return Vopt;
}
   
   void SxThermodynamics::print 
(const SxString &fileName, double startT, double endT, int samplePoints)
{
   int i;
   double CV, CVCl, SVib, SVibCl, HVib, HVibCl, GVib, GVibCl, 
   Trescaled, deltaT;

   double V;
   
   FILE *fp = NULL;
   SxString text;
/*
   V = 3.65;
   cout << V0 << endl;
   cout << (getFVib(V , 300.)) << endl;
   SX_EXIT;
  */       
   
   if ( !(fp = fopen (fileName.ascii (), "w")) ) {
      sxprintf ("Can't open file %s",fileName.ascii ()); 
      SX_EXIT;
   }     
   
   fprintf (fp,"------- Thermodynamic Analysis of System -----\n");
   fprintf (fp,"Model: System of quantum harmonic oscillators ");
   fprintf (fp,"(classical values are in parantheses)\n");
   fprintf (fp,"units: Energy - kcal/mol; Temperature - K\n");
   
   double zeroCorr = getZPV (V0);	
   fprintf (fp,"\nZero Point Vibrations in kcal/mol: %f\n\n" , zeroCorr); 
   
   SxArray<double> Temp (1);
   SxArray<double> P (1);

   P(0) = 0.;
   
   // --- an array of temperatures in generated correspondin to the 
   //     input parameters
   
   if ((fabs(startT-endT) <= 1e-10) || (samplePoints <= 1)) 
      Temp(0) = startT;
   else {
      Temp.resize(samplePoints);
      deltaT = (endT - startT)/(int)(samplePoints - 1);
      for (i = 0; i < samplePoints; i++) 
         Temp(i) =  startT + (double)i*deltaT;
   }
   
   // --- the thermodynamic properties are calculated for the given 
   //     temperature array and printed to an ascii-file
   for (i = 0; i < Temp.getSize (); i++) {
      V = solveEquOfState (P(0), Temp(i));
      fprintf (fp, 
            "-----------------------------------------------------------\n" );
      
      Trescaled = getRescaledTemp(V, Temp(i));
      fprintf(fp, 
            "Quantum mech. Temp.  %4.1f, rescaled Temp. %4.1f  \n\n"
            , Temp(i), Trescaled); 
      
      CV = getCVVib (V, Temp(i));
      CVCl = getCVVibCl (V, Temp(i));
      fprintf (fp, "CV :       %10.5f ( %10.5f )\n\n"
            , CV, CVCl);  
      
      HVib =  getUVib     (V, Temp(i));
      HVibCl =  getUVibCl (V, Temp(i));
      fprintf (fp, "    Hvib:   %10.3f ( %10.3f )\n",HVib, HVibCl);
      
      SVib = getSVib (V, Temp(i));
      SVibCl = getSVibCl (V, Temp(i));
      fprintf (fp, " -T*Svib:   %10.3f ( %10.3f )\n", 
            -Temp(i)*SVib, -Temp(i)*SVibCl);

      GVib =  HVib - Temp(i)*SVib;
      GVibCl =  HVibCl - Temp(i)*SVibCl;
      fprintf (fp, "    Gvib:   %10.3f ( %10.3f )\n", GVib, GVibCl);
   }
   fclose (fp);
}

