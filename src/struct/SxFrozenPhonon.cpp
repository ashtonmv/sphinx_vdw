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

#include <SxFrozenPhonon.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------
SxFrozenPhonon::SxFrozenPhonon 
(const SxAtomicStructure &tauIn,  SxPotential *pot) 
   : inputTau (tauIn),
     potential (pot)
{
   if (inputTau.cell.symGroupPtr->getNSymmorphic () > 1)  {
      cout << "SxFrozenPhonon does not support symmetries!" << endl;
      cout << "Specify symmetry {operator{S = [[1,0,0],[0,1,0],[0,0,1]];}} in the structure group!" << endl;
      SX_QUIT;
   }
   
   //hamSolver = dynamic_cast<SxHamSolver *> (pot);
   speciesData = potential->getSpeciesData ();
   
   nDoF = 3*inputTau.getNAtoms ();
	massVec.resize (nDoF);
	
	// a bit ugly ...
	int counter = 0.;
	for (int is = 0; is < inputTau.getNSpecies (); is++) {
		for (int ia = 0; ia < inputTau.getNAtoms (is); ia++) {
			for (int i = 0; i < 3; i++, counter++) {
				massVec(counter) = speciesData.ionicMass(is);
			}
		}
	}

   secondaryStructure = false;
}

SxFrozenPhonon::SxFrozenPhonon ()
{
   // empty
}

SxFrozenPhonon::~SxFrozenPhonon ()
{
   // empty
}

//----------------------------------------------------------------------------
//    Interface to input file
//----------------------------------------------------------------------------
bool SxFrozenPhonon::isRegistered (const SxSymbolTable *cmd)
{
   SX_CHECK (cmd);
   SxString str = cmd->getName ();
   return ( str == "frozenPhonon" );
}

double SxFrozenPhonon::getEPot ()
{
      return potential -> getPotentialEnergy ();
}

void SxFrozenPhonon::print (const SxSymbolTable *cmd)
{
   execute (cmd, false);
}

SxAtomicStructure SxFrozenPhonon::getResponse 
(const SxAtomicStructure &displacementIn, double dx, bool twoStep)
{
   SxAtomicStructure displacement;
   SxAtomicStructure force1, force2;
   SxAtomicStructure displaced;
   SxAtomicStructure returnValue;
   
   displacement.copy (displacementIn);
   double norm = sqrt(displacement.coordRef ().absSqr ().sum ());
   if (twoStep) {
      displaced = inputTau + dx*displacement;
      force1 = getForces (displaced);
      displaced = inputTau - dx*displacement;
      force2 = getForces (displaced);
      returnValue = (1./(2.*dx))*(force2 - force1);
   } else {
      displaced = inputTau + dx*displacement;
      force1 =  getForces (displaced) - inputForces;
      returnValue =  -(1./dx)*force1;
   }

      
   
   //--- response needs to be normed according to the 
   //    displacement length
   returnValue /= norm;
   return returnValue;
}


SxAtomicStructure SxFrozenPhonon::getDisplacement 
(int DoF, const SxString &basis)
{
   SxAtomicStructure returnValue;
   SxMatrix<Double> deviation (nDoF, nDoF);
   SxVector<Double> column(nDoF);
   SxVector<Double> orth(nDoF);
   returnValue.copy (inputTau);
   
   if (basis == "cartesian") {
      SxVector<Double> dispVec (nDoF);
      dispVec.set(0.);
      if (secondaryStructure) dispVec(peptideChain.getDoFInputOrder(DoF)) = 1.;
      else dispVec(DoF) = 1.;
      returnValue.set (dispVec);
   }
   
   if ((basis == "file") || (basis == "refine")) {
      if (secondaryStructure) {
         returnValue.set 
            (basisPeptideChain.getDisplacementVector (DoF));
      } else {
         returnValue.set (basisHOps.getDisplacementVector (DoF));
         
      }	
   }

   if (basis == "symmetry") {
      cout 
         << "SxFrozenPhonon:: Symmetry reduced basis not impl. yet" << endl;
   }
   return returnValue;
}

SxMatrix<Double> SxFrozenPhonon::getHessian (const SxMatrix<Double> &responses, const SxString &basis)
{
   SxMatrix<Double> returnValue (nDoF, nDoF);
   returnValue.set (0.);
   if (basis == "file") {
         returnValue.copy (responses);
   }
   if (basis == "cartesian") {
      if (secondaryStructure) {
         for (int i = 0; i < peptideChain.getNDoFPIU (); i++)  
            returnValue.colRef (peptideChain.getDoFInputOrder(i)).
               set (responses.colRef (i));
         returnValue =  returnValue.transpose ();
      } else 
         returnValue.copy (responses);
   }
   if (basis == "refine") {
      if (secondaryStructure) {
         peptideChain.setHessianFromRefinement 
            (responses.transpose (),basisHOps.getHessian (), 1000, false);
         returnValue.copy (peptideChain.getFullCartHessian ());
      } else {
         SxHessianOps hOps;
         hOps.setFromRefinement 
            (basisHOps, responses);
         returnValue.set (hOps.getHessian ());
         //returnValue.set (responses);
      }
   }
   if (basis == "symmetry") {
      cout << "SxFrozenPhonon:: Symmetry reduced basis not impl. yet" << endl;
   }
   return returnValue;
}
   

void SxFrozenPhonon::execute (const SxSymbolTable *cmd, bool /*calc*/)
{
   SX_CHECK (cmd);
   SxBinIO io;
   SxArtifactFilter f; 
   SxHessianOps hessianOps;
   SxString str = cmd->getName ();
   cout << "Frozen Phonon ..." << endl;
   
   elMinimCmds = potential->getMinimCmds (cmd);
         
   if (elMinimCmds.getSize() == 0)  {
      sxprintf ("No valid command found in group 'bornOppenheimer'\n");
      SX_EXIT;
   }

   SxString basis = SxString("cartesian");
   SxString fileName;
   if (cmd -> containsGroup ("basis")) {
      SxSymbolTable *basisGroup = cmd -> getGroup("basis");
      if (basisGroup -> contains("symmetry")) 
         basis = SxString("symmetry");
      if (basisGroup -> contains("cartesian")) 
         basis = SxString("cartesian");
      if (basisGroup -> contains("file")) {
         basis = SxString("file");
         fileName = basisGroup -> get ("file") -> toString ();
      }
      if (basisGroup -> containsGroup("secondaryStructure")) {
          SxSymbolTable *helixGroup = 
             basisGroup -> getGroup("secondaryStructure");
          SxString motif
             (helixGroup -> get("motif") -> toString ());
          bool toSymmetrize =  true;
          if (motif.contains("Sheet")) {
             toSymmetrize = false;
             cout << "TODO: FIX STRUCTURE-SYMMETRIZATION FOR BETA-SHEETS" << endl;
          }
          peptideChain = SxSecondaryStructure (inputTau, speciesData,  
                helixGroup -> get("periodicity") -> toString (),
                motif,
                helixGroup -> get("nTurns") -> toInt (), toSymmetrize);
          basisPeptideChain = SxSecondaryStructure (inputTau, speciesData,  
                helixGroup -> get("periodicity") -> toString (),
                motif,
                helixGroup -> get("nTurns") -> toInt (), toSymmetrize);
          secondaryStructure = true;
      }
   }

   bool twoStep = true; SxString freezeMode = SxString("notr"); 
   double deviation = 0.01; double refineDeviation = 0.;
   double refineFreqsDeviation = 0.; double refineUpTo = 0.;
   
   if (cmd -> containsGroup ("performance")) {
      SxSymbolTable *performanceGroup = cmd -> getGroup("performance");
      if (performanceGroup -> contains("freezeMode")) { 
         freezeMode = performanceGroup -> get("freezeMode") -> toString ();
      }
         
      if (performanceGroup -> contains("twoStep")) 
         twoStep = performanceGroup -> get("twoStep") -> toBool ();
      if (performanceGroup -> contains("deviation")) 
         deviation = performanceGroup -> get("deviation") -> toReal ();
      if (performanceGroup -> contains("refineDeviation")) 
         refineDeviation 
            = performanceGroup -> get("refineDeviation") -> toReal ();
      if (performanceGroup -> contains("refineFreqsDeviation")) {
         refineFreqsDeviation 
            = performanceGroup -> get("refineFreqsDeviation") -> toReal ();
         refineUpTo 
            = performanceGroup -> get("refineFreqsUpTo") -> toReal ();
      }
      cout << "in " << basis << " basis ...";
   }
   
   f.set (inputTau, freezeMode, false); 

   int startDof = 0; int endDof = nDoF - 1;
   if (secondaryStructure)   
      endDof = peptideChain.getNDoFPIU () - 1;

   if (cmd -> containsGroup ("dofRange") ) {
      SxSymbolTable *dofRangeGroup = cmd -> getGroup("dofRange");
      startDof = (dofRangeGroup -> get("startDof") -> toInt ()) - 1;
      endDof = (dofRangeGroup -> get("endDof") -> toInt ()) - 1;
   }
   
   SxAtomicStructure displacement;
   SxAtomicStructure response;
   SxMatrix<Double> responseMatrix(nDoF, nDoF);
   responseMatrix.set (0.);
   
   if (!twoStep) inputForces = getForces (inputTau);
   if (basis != SxString ("file")) {
      for (int i = startDof; i <= endDof; i++) {
         displacement = getDisplacement (i, basis);
         response = getResponse (displacement, deviation, twoStep);
         responseMatrix.colRef(i).set(response.coordRef ());
      }
   } else {
      responseMatrix = SxMatrix<Double> (nDoF, nDoF);
      io.open (fileName, SxBinIO::ASCII_READ_ONLY);
      io.read ("dummy", &responseMatrix, nDoF, nDoF, 0, 0);
      io.close();
   }   
  
   //--- in case that the user set's up the refineDeviation 
   //    parameter a refinement calculation on basis of the
   //    hessian as determined above
   if (refineDeviation > 1e-15) {
      basisHOps.set (getHessian(responseMatrix, basis), massVec);
      //basisHOps.set(responseMatrix.identity (), massVec);
      if (secondaryStructure) 
         basisPeptideChain.setHessian (basisHOps.getHessian (), 1000, false);
      basis = SxString("refine");
      for (int i = startDof; i <= endDof; i++) {
         displacement = getDisplacement (i, basis);
         response = getResponse (displacement, refineDeviation, twoStep);
         responseMatrix.colRef(i).set(response.coordRef ());
      }
   }

   // Hesse matrix has to be symmetric
   responseMatrix = 0.5 * (responseMatrix + responseMatrix.transpose ());

   SxMatrix<Double> hessian = getHessian (responseMatrix, basis);
   hessianOps = SxHessianOps (hessian, massVec);

   // --- a second refinement strategy can be applied
   //     (see SxSynchronousTransit)
   if (fabs(refineFreqsDeviation) > 1e-15)  {
      SxSynchronousTransit sTransit (inputTau, potential);
      sTransit.setElMinimCmds (elMinimCmds);
      hessianOps = sTransit.getRefinedDynamical 
         (hessianOps, refineUpTo, refineFreqsDeviation);
      hessian = hessianOps.getHessian ();
   }
   
   //--- ouptut of raw hessian for control purposes
   hessianOps.write ("hessian_raw.sx");
   
   if (secondaryStructure) {
         peptideChain.setHessian (hessian, 100, false);
         hessian = peptideChain.getFullCartHessian ();
         cout << "filtering artifacts ..."<< endl;
         f.setMasses (massVec);
         f.set (peptideChain.getCoords (), SxString("z"), false); 
         hessian = (f | hessian);
         peptideChain.setHessian (hessian,  100, false);
         hessian = peptideChain.getFullCartHessian ();
   }

   // ---applying artifact filter
   if (!secondaryStructure) hessian = (f | hessian);
   hessianOps = SxHessianOps (hessian, massVec);
   
      
   // ---standard output
   hessianOps.write ("hessian_end.sx");

   // ---optional output
   if (cmd -> containsGroup("output") ) {
      SxSymbolTable *outputGroup = cmd -> getGroup("output");
      if (outputGroup -> containsGroup("printDispersion")) {
         SxSymbolTable *dispGroup = outputGroup -> getGroup ("printDispersion");
         SxString dispersionFilename = 
            dispGroup -> get("file") -> toString();
         SxString interpolation = 
            dispGroup -> get("interpolation") -> toString ();
         
         peptideChain.printPhononDispersionRelation 
            (dispersionFilename, interpolation, 
             peptideChain.getPDRes() );
      }

      // --- output of thermodynamic data
      if (outputGroup -> containsGroup("printTD")) {

         //--- reading in output parameters
         SxSymbolTable *tdGroup = outputGroup -> getGroup ("printTD");
         SxString tdFilename = 
            tdGroup -> get("file") -> toString();
         SxSymbolTable *TGroup = tdGroup -> getGroup ("T");
         double startT = TGroup -> get("startT") -> toReal();
         double endT = TGroup -> get("endT") -> toReal();
         int tempSamplePoints = TGroup -> get("samplePoints") -> toInt();
         bool fromDispersion = tdGroup -> get("fromDispersion") -> toBool();

         //--- initialising thermodynamics
         //TODO: setup vibRots
         int vibRots = 3;
         SxThermodynamics td;
         if ((fromDispersion) && (secondaryStructure)) {
            SxString interpolation = SxString("spline");
            if (tdGroup -> contains("interpolation"))
               interpolation = tdGroup -> get("interpolation") -> toString ();
            td = SxThermodynamics (peptideChain, interpolation, 0);
         }
         else 
            td = SxThermodynamics (hessianOps, (6 - vibRots), 0.0001);
         
         //--- printing thermodynamic data to ascii-file
         td.print (tdFilename, startT, endT, tempSamplePoints);
      
      }
      
      //--- output of molden-file (normal-modes + frequencies)
      if (outputGroup -> containsGroup("printMolden")) {

         //--- reading in output parameters
         SxSymbolTable *moldenGroup = outputGroup -> getGroup ("printMolden");
         SxString moldenFilename = 
            moldenGroup -> get("file") -> toString();
         int replicZ = 1; 
         if (moldenGroup -> contains("replicZ") ) 
            replicZ = moldenGroup -> get("replicZ") -> toInt();
         
         hessianOps.printMolden (moldenFilename, inputTau, replicZ);
      }
      
      if (outputGroup -> containsGroup("printFC")) {

         //--- reading in output parameters
         SxSymbolTable *FCGroup = outputGroup -> getGroup ("printFC");
         SxString FCFilename = 
            FCGroup -> get("file") -> toString();
        if (!secondaryStructure) {
           cout 
              << "Force constant output only for secondary structure yet!"
              << endl;
           SX_QUIT;
        } else {
              peptideChain.printFC (FCGroup) ; 
        }
      }
   }
}


SxAtomicStructure SxFrozenPhonon::getForces (const SxAtomicStructure &tau)  
{
   return potential->getForces (tau, elMinimCmds);
}

