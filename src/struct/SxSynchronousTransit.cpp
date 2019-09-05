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

#include <SxSynchronousTransit.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------
SxSynchronousTransit::SxSynchronousTransit 
(const SxAtomicStructure &tauIn,  SxPotential *pot) 
   : inputTau (tauIn),
     potential (pot)
{
	int i, ia, is, counter;

   //hamSolver = dynamic_cast<SxHamSolver *> (pot);
   speciesData = potential->getSpeciesData ();
   
   nDoF = 3*inputTau.getNAtoms ();
	massVec.resize (nDoF);
	
	counter = 0.;
	for (is = 0; is < inputTau.getNSpecies (); is++) {
		for (ia = 0; ia < inputTau.getNAtoms (is); ia++) {
			for (i = 0; i < 3; i++) {
				massVec(counter) = speciesData.ionicMass(is);
				counter++;
			}
		}
	}
}

SxSynchronousTransit::SxSynchronousTransit ()
{
   // empty
}

SxSynchronousTransit::~SxSynchronousTransit ()
{
   // empty
}

//----------------------------------------------------------------------------
//    Interface to input file
//----------------------------------------------------------------------------
bool SxSynchronousTransit::isRegistered (const SxSymbolTable *cmd)
{
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   return ( str == "SynchronousTransit" );
}

double SxSynchronousTransit::getEPot ()
{
      return potential -> getPotentialEnergy ();
}

void SxSynchronousTransit::print (const SxSymbolTable *cmd)
{
   execute (cmd, false);
}

void SxSynchronousTransit::setElMinimCmds 
(const SxArray<const SxSymbolTable *> &in)
{
   elMinimCmds.resize (in.getSize ());
   for (int i = 0; i < elMinimCmds.getSize (); i++) 
      elMinimCmds(i) = in(i);
}

SxArray<double> SxSynchronousTransit::getEPots 
(const SxArray<SxAtomicStructure> &traj)
{
   //--dummy for forces here, but forces could also be of interest
   SxAtomicStructure dummy;
   SxArray<double> returnValue (traj.getSize ());
   for (int i = 0; i < traj.getSize (); i++) {
      dummy = getForces (traj(i));
      returnValue(i) = getEPot ();
   }
   return returnValue;

}
SxHessianOps SxSynchronousTransit::getRefinedDynamical 
         (const SxHessianOps &in, double refineUpTo, double dX)
{
   FILE *fp = NULL;
   int i;
   SxHessianOps returnValue (in);
   SxList<int> indices;
   SxList<SxComplex16> corrections;
   
   // lower barrier for frequencies to refine
   double refineHigherThan = 1.0;
   SxVector<Complex16> freqs = returnValue.getEigenFrequencies ();
   SxVector<Double> reducedMasses = returnValue.getReducedMasses ();
   SxMatrix<Double> eigenVelocities = returnValue.getEigenVelocities ();

   //--- control file
   if ( !(fp = fopen ("refine.dat", "w")) ) {
      sxprintf ("Can't open file refine.dat"); 
      SX_EXIT; 
   }     
   SxString line
      ("Orig.  Frequency    Correction           error         (unit 1/cm)\n");
   line += SxString
      ("------------------------------------------------------------------\n");
   fprintf (fp, "%s", line.ascii());
   
   cout << "Refining Frequencies ..." << endl;
   cout << "refine up to " << refineUpTo << " 1/cm " << endl;
   cout << "dX = " << dX << " Bohr " << endl;
   
   //--- indices of frequencies to correct are determined
   for (i = 0; i < nDoF; i++) {
      if ( (freqs(i).abs () > refineHigherThan) 
           && (freqs(i).abs () < refineUpTo) ) {
         indices.append(i);
      }
   }
   cout << "Frequencies to correct :" << endl;
   for (i = 0; i < indices.getSize (); i++) 
      cout << indices(i) << " " << freqs(indices(i)) << endl;
      
   double dS;
   int index, j;
   SxVector<Double> mode (nDoF);
   SxArray<SxAtomicStructure> traj (5);
   SxAtomicStructure displacedStructure (inputTau);
   SxArray<double> pes;
   double y1, x1, y2, x2, dx, dy, k1, k2, dxs, errorbar;
   SxComplex16 freq1, freq2, origFreq, correction;
          
   
   for (i = 0; i < indices.getSize (); i++) {
      index = indices(i);
      mode.set (sqrt(reducedMasses(index))
            *(eigenVelocities.colRef (index)));
      origFreq = freqs(index).abs ();

      //--- calculations of suitable displacements 
      if (origFreq.abs () < 10.0) origFreq = 10.; 
      dS = dX / (origFreq / 5123.75 * sqrt(reducedMasses(index)));
      
      for (j = -2; j <= 2; j++) {
         displacedStructure.copy (inputTau);
         displacedStructure.set (inputTau.coordRef () 
               + mode*(double)j/2.*dS);
         traj(j + 2 ).copy (displacedStructure);
      }
      pes = getEPots (traj);
      
      //--- curvature as calculated from 3-point harmonic fit 
      //    with middle point and outer two points 
      //    of the 5 point trajectory
      y1 = pes(4) - pes(2);
      y2 = pes(0) - pes(2);
      x1 = dS; x2 = -dS;
      dy = y1 - y2; dx = x1 - x2; dxs= x1*x1 - x2*x2;
      
      k1 = 0.5*( (y1 - x1*dy/dx) / (x1*x1 - dxs/dx) ); 
      
      //--- curvature as calculated from 3-point harmonic fit 
      //    with middle point and inner two points of the 
      //    5 point trajectory
      y1 = pes(3) - pes(2);
      y2 = pes(1) - pes(2);
      x1 = dS/2.; x2 = -dS/2.;
      dy = y1 - y2; dx = x1 - x2; dxs= x1*x1 - x2*x2;
      
      k2 = 0.5*( (y1 - x1*dy/dx) / (x1*x1 - dxs/dx) ); 
      
      //--- the frequencies are calculated
      
      freq2 = 2*sqrt(k2/reducedMasses(index))*5123.75;
      freq1 = 2*sqrt(k1/reducedMasses(index))*5123.75;
     
      //--- for the corrected value an average of the two 
      //    frequencies is taken
      correction = (freq2 + freq1)/2.;
      errorbar = ((freq2 - freq1).abs ()/2.);
      origFreq = freqs(index);
      
      corrections.append (correction);
      // --- control output
      for (j = 0; j < 5; j++) 
         sxprintf("PES%d %f\n",i , pes(j));
      
      if (origFreq.re > origFreq.im) 
        line = SxString(origFreq.re, "%14.4f") 
             + SxString ("  ") ;
      else 
        line = SxString(origFreq.im, "%14.4f") 
             + SxString ("*I");
      if (correction.re > correction.im) 
        line += SxString(correction.re, "%14.4f") 
              + SxString ("  ");
      else 
        line += SxString(correction.im, "%14.4f") 
              + SxString ("*I");
        line += SxString(errorbar, "%14.4f");
      fprintf (fp, "%s\n", line.ascii ());
      fflush(fp);

   }
   fclose (fp);

   //--- corrected values are incorporated to the dynamical matrix
   SxList<SxList<SxComplex16> > toReplace;
   SxList<SxComplex16> replaceEntry;
   for (i = 0; i < corrections.getSize (); i++) {
      replaceEntry.resize (0);
      replaceEntry.append ((double)(indices(i) + 1));
      replaceEntry.append (corrections(i));
      toReplace.append (replaceEntry);
   }
   returnValue.setFrequencies (toReplace);
   return returnValue;
}



void SxSynchronousTransit::execute (const SxSymbolTable *cmd, bool /*calc*/)
{
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   cout << "Synchronous Transit ...";
   elMinimCmds = potential->getMinimCmds (cmd);
         
   if (elMinimCmds.getSize() == 0)  {
      sxprintf ("No valid command found in group 'bornOppenheimer'\n");
      SX_EXIT;
   }

}


SxAtomicStructure SxSynchronousTransit::getForces (const SxAtomicStructure &tau)  
{
   return potential->getForces (tau, elMinimCmds);
}

