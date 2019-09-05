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

// sxvibrations by Lars Ismer !

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <SxError.h>
#include <SxPrecision.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSymbolTable.h>
#include <SxVector.h>
#include <SxAtomicStructure.h>
#include <SxBinIO.h>
#include <SxUtil.h>
#include <SxElemDB.h>
#include <SxSecondaryStructure.h>
#include <SxArtifactFilter.h>
#include <SxHessianOps.h>
#include <SxVibrations.h>
#include <SxThermodynamics.h>

#ifndef SX_STANDALONE

SxVibrations::SxVibrations ()
{
   //--- INPUT
   inputHessian           = "empty";
   inputStructure         = "empty";
   inputStructure2         = "empty";
   inputPeriodicity       = "empty";
   inputPeriodicity2      = "empty";
   inputEigenvalues       = "empty";
   inputHistory           = "empty";
   inputBasisHessian           = "empty";

   //---SYSTEM
   systemClass = SxString ("empty");
   sigma = 0;

   //---PERFORMANCE
   lowCoop =  false;
   overwriteFreezeMode = false;
   overwriteVibRots = false;
   applyPeriodicity = false;
   applyAveraging = false;
   setZero = 0.;
   supercellCutoff = 0;
   oldCalc = false;
   isDynamical = false;
   vdwCorrection = SxString("none");

   //---OUTPUT
   grueneisen = false;
   moldenFilename         = SxString("empty");
   anharmonicFilename     = SxString("empty");
   thermodynamicsFilename = SxString("empty");
   eigenvaluesFilename    = SxString("empty");
   forceConstantsFilename = SxString("empty");
   hessianFilename        = SxString("empty");
   dispersionFilename     = SxString("empty");
   //rotorType              = SxString("none");
   printMolden         = false;
   printAnharmonic     = false;
   printEigenvalues    = false;
   printForceConstants = false;
   printNuMu           = false;
   printThermodynamics = false;
   printHessian = false;
   printDispersion = false;
   printCouplings = false;
   fromDispersion =  false;
   fileBasis = false;
   interpolation = SxString ("fourier");
   tdInterpolation = SxString ("fourier");
   //incorporateRotor = false;

}
SxVibrations::~SxVibrations ()
{
   //empty
}

void SxVibrations::parseInputFile ()
{
   SxParser parser;
   table = parser.read (inputFile);

   try {

      //--- Group INPUT
      SxSymbolTable *input = table -> getGroup("input");

      //non - optional
      inputHessian = input -> get ("inputHessian") -> toString ();

      //optional

      if (input -> contains ("inputStructure"))
        inputStructure = input -> get ("inputStructure") -> toString();

      if (input -> contains ("inputStructure2"))
        inputStructure2 = input -> get ("inputStructure2") -> toString();

      if (input -> contains ("inputPeriodicity"))
        inputPeriodicity = input -> get ("inputPeriodicity") -> toString();

      if (input -> contains ("inputPeriodicity2"))
        inputPeriodicity2 = input -> get ("inputPeriodicity2") -> toString();

      if (input -> contains ("inputEigenvalues"))
        inputEigenvalues = input -> get ("inputEigenvalues") -> toString();

      if (input -> contains ("inputHistory"))
        inputHistory = input -> get ("inputHistory") -> toString();

      if (input -> contains ("inputBasisHessian")) {
        inputBasisHessian = input -> get ("inputBasisHessian") -> toString();
        fileBasis = true;
      }

      //---Group SYSTEM

      SxSymbolTable *system = table -> getGroup("system");

      //non - optional
      systemClass = system -> get ("class") -> toString ();
      if (systemClass == SxString ("helix")
            || systemClass == SxString("FES")) {
         nTurns = system -> get("nTurns") -> toInt ();
      }

      //optional
      if (system -> contains ("sigma"))
        sigma = system -> get ("sigma") -> toInt();


      //---Group PERFORMANCE

      if (table -> containsGroup ("performance")) {


         SxSymbolTable *performance = table -> getGroup("performance");

         //optional

         if (performance -> contains ("oldCalc"))
            oldCalc = performance -> get ("oldCalc") -> toBool();

         if (performance -> contains ("isDynamical"))
            isDynamical = performance -> get ("isDynamical") -> toBool();

         if (performance -> contains ("lowCoop"))
            lowCoop = performance -> get ("lowCoop") -> toBool();


         if (performance -> contains ("setZero"))
            setZero = performance -> get ("setZero") -> toReal();

         if (performance -> contains ("freezeMode")) {
            overwriteFreezeMode = true;
            freezeMode = performance -> get ("freezeMode") -> toString();
         }

         if (performance -> contains ("vibRots")) {
            overwriteVibRots = true;
            vibRots = performance -> get ("vibRots") -> toInt();
         }

         if (performance -> contains ("applyPeriodicity")) {
            applyPeriodicity
               = performance -> get ("applyPeriodicity") -> toBool();
         }

         if (performance -> contains ("applyAveraging")) {
            applyAveraging
               = performance -> get ("applyAveraging") -> toBool();
         }

         if (performance -> contains ("supercellCutoff")) {
            supercellCutoff
               = performance -> get ("supercellCutoff") -> toInt ();
         }

         if (performance -> contains ("vdwCorrection")) {
            vdwCorrection
               = performance -> get ("vdwCorrection") -> toString ();
         }


         if (performance -> contains ("replaceCurvature")) {
           replaceCurvature = true;
           SxList<double> list
              = performance -> get ("replaceCurvature") -> toList ();

           replaceList.resize (list.getSize () / 2);

           for (int i = 0; i < replaceList.getSize (); i++) {
              replaceList(i).resize (2);
              for (int j = 0; j < 2; j++) {
               replaceList(i)(0) = list(2*i);
               replaceList(i)(1) = list(2*i + 1);
              }
           }
         }

      }
      //---Group OUTPUT

      SxSymbolTable *output = table -> getGroup("output");

         //optional
      if (output -> containsGroup ("printMolden")) {
         SxSymbolTable *printGroup = output -> getGroup ("printMolden");
         printMolden = true;
         moldenFilename = printGroup -> get ("file") -> toString ();
         if (printGroup -> contains ("replicZ")) {
            replicZ = printGroup -> get ("replicZ") -> toInt ();
            if (replicZ > 1)
               zExt = printGroup -> get("zExt") -> toReal ();
            else zExt = 0.;

         }
      }

      if (output -> containsGroup ("printEV")) {
         SxSymbolTable *printGroup = output -> getGroup ("printEV");
         printEigenvalues = true;
         eigenvaluesFilename = printGroup -> get ("file") -> toString ();
      }

      if (output -> containsGroup ("printFC")) {
         FCGroup = output -> getGroup ("printFC");
         printForceConstants = true;
      }

      if (output -> containsGroup ("printNuMu")) {
         SxSymbolTable *printGroup = output -> getGroup ("printNuMu");
         printNuMu = true;
         nuMuFilename = printGroup -> get ("file") -> toString ();
      }

      if (output -> containsGroup ("printTD")) {
         SxSymbolTable *printGroup = output -> getGroup ("printTD");
         printThermodynamics = true;
         fromDispersion = printGroup -> get ("fromDispersion") -> toBool ();
         thermodynamicsFilename = printGroup -> get ("file") -> toString ();
         SxSymbolTable *TGroup = printGroup -> getGroup ("T");
         startT = TGroup -> get("startT") -> toReal();
         endT = TGroup -> get("endT") -> toReal();
         tempSamplePoints = TGroup -> get("samplePoints") -> toInt();
         if (printGroup -> contains ("p")) {
            p = printGroup -> get ("p") -> toReal ();
         if (printGroup -> contains ("interpolation"))
         tdInterpolation
            = printGroup -> get ("interpolation") -> toString ();
         }
         if (printGroup -> containsGroup ("grueneisen")) {
            grueneisenGroup = printGroup -> getGroup ("grueneisen");
            grueneisen = true;
         }
        /*
         SxSymbolTable *RGroup = printGroup -> getGroup ("rotor");
         rotorType = RGroup -> get("rotorType") -> toString();
         rotorBarrier = RGroup -> get("barrier") -> toReal();
         branchToReplace = RGroup -> get("branchToReplace") -> toInt ();
         */
      }


      if (output -> containsGroup ("printAH")) {
         SxSymbolTable *printGroup = output -> getGroup ("printAH");
         printAnharmonic = true;
         anharmonicFilename = printGroup -> get ("file") -> toString ();
         ahId = printGroup -> get ("ahId") -> toInt ();
      }

      if (output -> containsGroup ("printHessian")) {
         SxSymbolTable *printGroup = output -> getGroup ("printHessian");
         printHessian = true;
         hessianFilename = printGroup -> get ("file") -> toString ();
      }

      if (output -> containsGroup ("printDispersion")) {
         SxSymbolTable *printGroup = output -> getGroup ("printDispersion");
         printDispersion = true;
         dispersionFilename = printGroup -> get ("file") -> toString ();
         res = printGroup -> get ("nKPoints") -> toInt ();
         interpolation = printGroup -> get ("interpolation") -> toString ();
      }

      if (output -> containsGroup ("printCouplings")) {
         SxSymbolTable *printGroup = output -> getGroup ("printCouplings");
         printCouplings = true;
         couplingsFilename = printGroup -> get ("file") -> toString ();
      }
   } catch (SxException e) {
      e.print();
      SX_QUIT;
   }

}

void SxVibrations::printParameters ()
{
   cout << "---------------------------------------------" << endl;
   cout << "Input Parameters : " << endl << endl;
   cout << "Input-Structure  : " << inputStructure.ascii () << endl;
   cout << "Input-Structure 2  : " << inputStructure2.ascii () << endl;
   cout << "Input-Basis Hessian : " << inputBasisHessian.ascii () << endl;
   if (inputPeriodicity != SxString ("empty"))
      cout << "Input-Periodicty : " << inputPeriodicity.ascii () << endl;
   if (inputPeriodicity2 != SxString ("empty"))
      cout << "Input-Periodicty 2: " << inputPeriodicity2.ascii () << endl;
   if (inputEigenvalues != SxString ("empty"))
      cout << "Input-Eigenvalues : " << inputEigenvalues.ascii () << endl;
   if (inputHistory != SxString ("empty"))
      cout << "Input-History     : " << inputHistory.ascii () << endl;

   cout << "---------------------------------------------" << endl;
   cout << "System Parameters : " << endl << endl;
   cout << "System class:          " << systemClass.ascii () << endl;
   if ( sigma != 0)
      cout << "System symmetry-sigma: " << sigma << endl;

   cout << "---------------------------------------------" << endl;
   cout << "Performance Parameters : " << endl << endl;

   if  (setZero >= 0.)
      cout << "setZero: " << setZero << endl;


   if  (applyPeriodicity)
      cout << "Apply Periodicity" <<  endl;

   if  (lowCoop)
      cout << "Low Cooperativity Calc" <<  endl;

   if  (applyAveraging)
      cout << "Apply Averaging" <<  endl;

   if  (applyAveraging)
      cout << "Van der Waals Correction: " <<  vdwCorrection.ascii () << endl;

   cout << "---------------------------------------------" << endl;
   cout << "Output : " << endl << endl;
   if  (printMolden) {
      cout << "Save molden-format file to: " << moldenFilename.ascii () << endl;
      cout << "             z-replication: " << replicZ  << endl;
   }

   if  (printEigenvalues) {
      cout << "Save eigenvalues to       : " << eigenvaluesFilename.ascii ()
                                             << endl;
   }

   if  (printForceConstants) {
      cout << "Save force constants to     : "
           << forceConstantsFilename.ascii () << endl;
   }

   if  (printNuMu) {
      cout << "Save frequency vs reduced Mass infos to     : "
           << nuMuFilename.ascii () << endl;
   }

   if  (printThermodynamics) {
      cout << "Print harmonic thermodynamic analysis to : "
         << thermodynamicsFilename.ascii () << endl;
      if ( p >= 0.) {
      cout << "                          pressure (atm) :  "
           << p << endl;
      }
      cout << "                          temperature (K):  "
         << startT << " to "<< endT << endl;
      cout << "                          sample Points  :  "
         << tempSamplePoints << endl;
      /*
      cout << "Rotor: " << endl;
      cout << "   Type:            " << rotorType.ascii () << endl;
      cout << "   Barrier:         " << rotorBarrier << endl;
      cout << "   branchToReplace: " << branchToReplace << endl;
  */
  }

   if  (printAnharmonic) {
      cout << "Print anhharmonic analysis to : "
         << anharmonicFilename.ascii () << endl;
      cout << "  index of analysed eigenmode :  "
         << ahId << endl;
   }

   if  (printHessian) {
      cout << "Print Hessian to : "
         << hessianFilename.ascii () << endl;
   }

   if  (printDispersion) {
      cout << "Print phonon-dispersion to : "
         << dispersionFilename.ascii () << endl;
   }

   if  (printCouplings) {
      cout << "Print couplings to : "
         << couplingsFilename.ascii () << endl;
   }
}

void SxVibrations::flushTemplate (const SxString &filename)
{
   FILE *fp = NULL;

   if ( !(fp = fopen (filename.ascii (), "w")) ) {
      sxprintf ("Can't open file %s",filename.ascii ());
      SX_QUIT;
   }

   fprintf (fp,"format vibrations;\n\n");

   fprintf (fp,"input {\n");
   fprintf (fp,"   inputStructure    = \"tau_end.out\"; \n");
   fprintf (fp,"   inputHessian      = \"hessian_end.out\"; \n");
   fprintf (fp,"//   inputHistory      = \"relaxHist.out\"; \n");
   fprintf (fp,"//   inputStructure2    = \"tau_end.out\"; \n");
   fprintf (fp,"//   inputEigenvalues  = \"ewcm.out\";\n");
   fprintf (fp,"//   inputPeriodicity  = \"periodicity.dat\";\n");
   fprintf (fp,"//   inputPeriodicity2  = \"periodicity.dat\";\n");
   fprintf (fp,"//   inputBasisHessian = \"hessian_end.out\";\n");
   fprintf (fp,"}\n");

   fprintf (fp,"\nsystem {\n");
   fprintf (fp,"   class   = \"molecule\"; \n");
   fprintf (fp,"   sigma   = 0; \n");
   fprintf (fp,"//   nTurns  = 3;\n");
   fprintf (fp,"}\n");


   fprintf (fp,"\nperformance {\n");
   fprintf (fp,"   freezeMode = \"z\";\n");
   fprintf (fp,"   vibRots    = 2;\n");
   fprintf (fp,"   setZero    = 0; \n");
   fprintf (fp,"   applyPeriodicity = false;\n");
   fprintf (fp,"   applyAveraging   = false;\n");
   fprintf (fp,"   supercellCutoff = 5;\n");
   fprintf (fp,"   oldCalc = false;\n");
   fprintf (fp,"   isDynamical = false;\n");
   fprintf (fp,"   vdwCorrection = \"none\";\n");
   fprintf (fp,"   lowCoop = false;\n");
   fprintf (fp,"   replaceCurvature = [[1,0.1]];\n");
   fprintf (fp,"}\n\n");

   fprintf (fp,"output {\n");
   fprintf (fp,"  printMolden {\n");
   fprintf (fp,"     file = \"modes.molf\";\n");
   fprintf (fp,"     replicZ = 1;\n");
   fprintf (fp,"     zExt = 10.95;\n");
   fprintf (fp,"  }\n");
   fprintf (fp,"  //printEV {\n");
   fprintf (fp,"  //   file = \"ewcm.out\";\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printFC {\n");
   fprintf (fp,"  //   file = \"fc.out\";\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printNuMu {\n");
   fprintf (fp,"  //   file = \"numu.out\";\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printTD {\n");
   fprintf (fp,"  //   fromDispersion = false;\n");
   fprintf (fp,"  //   file = \"thermo.out\";\n");
   fprintf (fp,"  //   p = 1.0;\n");
   fprintf (fp,"  //   T {\n");
   fprintf (fp,"  //     startT = 0;\n");
   fprintf (fp,"  //     endT = 298.;\n");
   fprintf (fp,"  //     samplePoints = 298;\n");
   fprintf (fp,"  //   }\n");
   fprintf (fp,"  //rotor {");
   fprintf (fp,"  //   rotorType =\"none\"");
   fprintf (fp,"  //   barrier = 0.");
   fprintf (fp,"  //   branchToReplace = 0 ");
   fprintf (fp,"  //   }\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printAH {\n");
   fprintf (fp,"  //   file = \"anharmonic.out\";\n");
   fprintf (fp,"  //   ahId = 1;\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printHessian {\n");
   fprintf (fp,"  //   file = \"hessian_gen.out\";\n");
   fprintf (fp,"  //}\n");
   fprintf (fp,"  //printDispersion {\n");
   fprintf (fp,"  //   file = \"dispersion.out\";\n");
   fprintf (fp,"  //   nKPoints = 10;\n");
   fprintf (fp,"  //   supercellCutoff = 20;\n");
   fprintf (fp,"  //   interpolation = \"fourier\";\n");
   fprintf (fp," //}\n");
   fprintf (fp,"  //printCouplings {\n");
   fprintf (fp,"  //   file = \"couplings.out\";\n");
   fprintf (fp,"  //}\n");

   fprintf (fp,"}\n");
   fclose (fp);
}


void SxVibrations::validateInput ()
{
   if (     systemClass != SxString ("molecule")
         && systemClass != SxString ("bulk")
         && systemClass != SxString ("helix")
         && systemClass != SxString ("FES")
         && systemClass != SxString ("pSheet")
         && systemClass != SxString ("apSheet")
         && systemClass != SxString ("notr") ) {
      cout  << "Unknown system class: "
             << systemClass.ascii () << endl
            << "System class must be one of the following: "
            << "bulk, molecule, helix" << endl;
      exit (0);
   }

   if (       interpolation != SxString ("fourier")
         &&   interpolation != SxString ("linear")
         &&   interpolation != SxString ("spline") ) {
      cout << "Unknown interpolation type for phonon-dispersion relation:"
           << interpolation.ascii () << endl
           << "Known interpolation types are the following: "
           << "fourier, linear, spline" << endl;
   }


   if ((sigma <= 0) && (printThermodynamics)) {
      cout << "Sigma can't be <= 0,  must be an integer >= 1" << endl;
      cout << "exiting !" << endl;
      exit (0);
   }

   if ((printAnharmonic) && (inputHistory == SxString("empty"))) {
      cout << "Anharmonic analysis needs an history-file as input" << endl;
      cout << "Please provide an history-file from "
           << "a synchronous transit" << endl;
      exit (0);
   }

   if ((applyPeriodicity) && inputPeriodicity == SxString("empty")) {
      cout
         << "You must provide provide a file which contains the periodicity"
         << endl;
      exit (0);
   }

   if ((applyPeriodicity) && applyAveraging) {
      cout
         << "You have to choose between applyPeriodicity and applyAveraging"
         << endl;
      exit (0);
   }

   if (inputPeriodicity != SxString("empty")) {
      if (     (systemClass != "helix")
            && (systemClass != "FES")
            && (!(systemClass.contains("Sheet")))) {
         if (((supercellCutoff > peptideChain.getNPeptides ()/2)
                  || (supercellCutoff < 0))&&(!lowCoop)) {
            cout << "Supercell-Cutoff must be larger than 1 and it makes no "
               << "sense to set it larger then the number of irreducible "
               << "elements divided by 2" << endl;
            cout << "Number of Irreducible Elements divided by 2 "
               << peptideChain.getNPeptides ()/2
               << endl;
            exit (0);
         }
      }

      if (inputPeriodicity2 != SxString("empty")) {
      }

      if (lowCoop) {
         if (peptideChain.getNPeptides () != 5) {
            cout << "In the actual implementation 5 peptides are"
                 << " expected in the periodicity file, i.e. you must"
                 << " setup 5 columns!" << endl;
            SX_QUIT;
         }
      }
   }

   if ((fromDispersion) && !(printDispersion)) {
      cout << " You must setup your phonon dispersion"
           << " (printDispersion) when choosing the "
           << " fromDispersion flag !" << endl;
      SX_QUIT;
   }
  /*
   if ((rotorType == SxString("none")) || (rotorType == SxString("cosinus"))
         || rotorType == SxString("blocked")) {
      if (rotorType == SxString("none")) {
         incorporateRotor = false;
         branchToReplace = -1;
      }
   }
   else {
      cout << "Unknown Rotor Type:" << rotorType.ascii () << endl;
      cout << "Rotor type can be: none, cosinus."<< endl;
      SX_QUIT;
   }
      */
   if ( (vdwCorrection != SxString("Function_I")) &&
        (vdwCorrection != SxString("Function_II"))&&
        (vdwCorrection != SxString("none"))) {

      cout << "Unknown vdw correction:" << vdwCorrection.ascii () << endl;
      SX_QUIT;
   }

   if (lowCoop &&
         (  (inputPeriodicity == SxString("empty"))
          ||(inputPeriodicity2 == SxString("empty"))
          ||(inputStructure2 == SxString("empty"))   )) {
      cout << "You must provide two periodicity files and two structure files";
      cout << "for the low coop modus";
     SX_QUIT;
   }
}

void SxVibrations::initSystem ()
{
   cout << endl << endl << endl;
   cout << "*************************************" << endl;
   cout << "Initializing system ... " << endl;

   cout << "loading structure, " << endl;

   if (inputStructure == "empty") {
      cout << "loading structure from input.sx file !" << endl;
      // --- read from input file
      SxString inFile = SxString ("input.sx");

      SxParser parser;
      SxParser::Table table2 = parser.read (inFile, "std/structure.std");

      try  {
         tau = SxAtomicStructure(table2 -> getGroup("structure"));
         speciesData = SxSpeciesData(&*table2);

      } catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }
   } else  {
      tau = loadStructureFHI98 ();
   }

   nDoF = tau.nTlAtoms * 3;
   getMasses ();
   if ((inputHessian != SxString ("empty"))) {
      cout << "loading hessian, " << endl;
      hessian = SxMatrix<Double> (nDoF, nDoF);
      SxBinIO io (inputHessian, SxBinIO::ASCII_READ_ONLY);
      io.read ("dummy", &hessian, nDoF, nDoF, 0, 0);
   }

   if ((inputBasisHessian != SxString ("empty"))) {
      cout << "loading hessian, " << endl;
      basisHessian = SxMatrix<Double> (nDoF, nDoF);
      SxBinIO io (inputBasisHessian, SxBinIO::ASCII_READ_ONLY);
      io.read ("dummy", &basisHessian, nDoF, nDoF, 0, 0);
   }

   if (inputPeriodicity != SxString ("empty")) {
      if (inputStructure != SxString ("empty")) {
          bool toSymmetrize = true;
          if (systemClass.contains ("Sheet")) {
             //toSymmetrize = false;
             toSymmetrize = true;
             cout << "TODO: FIX STRUC-SYM for BETA-SHEETS" << endl;
          }
         peptideChain =
            SxSecondaryStructure
            (inputStructure, inputPeriodicity,
             systemClass, nTurns, toSymmetrize);
         tau = peptideChain.getCoords ();
      } else {
          bool toSymmetrize = true;
          if (systemClass.contains ("Sheet")) {
             //toSymmetrize = false;
             toSymmetrize = true;
             cout << "TODO: FIX STRUC-SYM for BETA-SHEETS" << endl;
          }

          peptideChain = SxSecondaryStructure (tau, speciesData,
                inputPeriodicity,
                systemClass,
                nTurns, toSymmetrize);
         tau = peptideChain.getCoords ();
      }
      if (peptideChain.getZPitch() <= 1e-15) {
         cout << " You must setup a valid supercell geometry"
              << " in the structure file !"
            << " z-Extension cannot be zero" << endl;
         SX_QUIT;
      }
   }

   cout << "setting up system properties, " << endl;
   if (systemClass == SxString("molecule")) {
      vibRots = 0;
      if (!overwriteFreezeMode) freezeMode = SxString("xyz");
   }

   if (systemClass == SxString("bulk")) {
      vibRots = 3;
      if (!overwriteFreezeMode) freezeMode = SxString("trans");
   }

   if ( (systemClass == SxString ("helix") ) ||
        (systemClass == SxString ("FES") )) {
      vibRots = 2;
      if (!overwriteFreezeMode) freezeMode = SxString("z");
   }

   if (systemClass == SxString("notr")) {
      vibRots = 6;
      if (!overwriteFreezeMode) freezeMode = SxString ("notr");
   }

   cout << "... ready!" << endl;

}

void SxVibrations::manipulateInput ()
{
   if (fileBasis) {
   //--- in case of refinement calculation, the input hessian must be
   //     transformed from normal to cartesian coordinates:
   //     treatment for linear chains:
      if ((applyPeriodicity)) {
         peptideChain.setHessianFromRefinement
            (hessian, basisHessian, supercellCutoff, applyAveraging);
         hessian = peptideChain.getFullCartHessian ();
      }
   } else {

      if ((applyPeriodicity) || (applyAveraging)) {
         peptideChain.setHessian (hessian, supercellCutoff, applyAveraging);
         hessian = peptideChain.getFullCartHessian ();
      }
   }

   cout << "...Symmetrizing Hessian," << endl;
   symmetrizeHessian (&hessian);


   cout << "setting entrys with absolute values < " << setZero
        << " to zero ... " << endl;
   setElementsToZero (&hessian, setZero);

   cout << "...ready!" << endl;

   if (replaceCurvature) {
     cout << "replacing curvatures ..."<< endl;
     hessianOps = SxHessianOps (hessian, masses.coordRef ());
     hessianOps.setCurvaturesDynBasis (replaceList);
     hessian = hessianOps.getHessian ();
   }

   cout << "filtering artifacts ..."<< endl;
   SxArtifactFilter f;
   f.set (tau, freezeMode, false);
   f.setMasses (masses.coordRef ());
   f.optimizeDMatrix (hessian);
   hessian = (f | hessian);
   if (applyPeriodicity || applyAveraging)
      peptideChain.setHessian (hessian, supercellCutoff, applyAveraging);
}

void SxVibrations::printOutput ()
{
     hessianOps = SxHessianOps (hessian, masses.coordRef ());

   if (printMolden)
      //saveEigenvecMolden
      //   (moldenFilename,tau, hessianOps, replicZ);
      hessianOps.printMolden (moldenFilename, tau, replicZ);

   if (printDispersion)
         printPhononDispersion ();


   if (printThermodynamics) {

      SxThermodynamics td;
      if (fromDispersion) {
         td = SxThermodynamics (peptideChain, tdInterpolation, 0);
      } else {
         td = SxThermodynamics (hessianOps, (6 - vibRots), 0.);
      }
      if (grueneisen) {
         cout << "Setting up Grueneisen Parameters ..." << endl;
         td.setGrueneisen (grueneisenGroup);
      }
         td.print (thermodynamicsFilename, startT, endT, tempSamplePoints);
   }

   //--- Comment: output of featured not yet transfered
   //    from rel-10 is commented out here

   if (printForceConstants)  {
      peptideChain.printFC (FCGroup) ;
   }

/*
   if (printEigenvalues)
         io.saveEwcm (eigenvaluesFilename, eig.vals);

   if (printNuMu)
         printFreqVsMass ();

   if (printAnharmonic)
      printAnharmonicCorrections ();

   if (printHessian)
         io.saveHessian(hessianFilename, hessian);

   if (printCouplings)
         printPhononCouplings ();

   fflush (stdout);
*/
   cout << " ready! " << endl;
}



void SxVibrations::printAnharmonicCorrections ()
{
        // -- empty: to be transfered from rel-10
}

void SxVibrations::printPhononDispersion ()
{
   peptideChain.printPhononDispersionRelation
      (dispersionFilename, interpolation, peptideChain.getPDRes());
}


void SxVibrations::printPhononCouplings ()
{
   // -- empty: to be transfered from rel-10
   SX_EXIT;
}

void   SxVibrations::printFreqVsMass ()
{
   // -- empty: to be transfered from rel-10
   SX_EXIT;
}

SxAtomicStructure SxVibrations::loadStructureFHI98 ()
{
   SxList<SxList<SxVector3<Double> > > tauList;
   SxMatrix3<Double> cell;
   int i, j;

   SxBinIO io;
   io.open (inputStructure, SxBinIO::ASCII_READ_ONLY);
   tauList = io.loadStructureFHI98 ();
   io.close ();
   io.open (inputStructure, SxBinIO::ASCII_READ_ONLY);
   speciesNameList = io.loadSpeciesFHI98 ();
   io.close ();
   io.open (inputStructure, SxBinIO::ASCII_READ_ONLY);
   cell = io.loadCellFHI98 ();
   io.close ();

   SxAtomicStructure tauFHI (cell);
   tauFHI.startCreation ();

   //Should be iterators (but not a bottleneck here)
   for (j = 0; j < tauList.getSize (); j++) {
      tauFHI.newSpecies ();
      for (i = 0; i < tauList(j).getSize (); i++) {
         tauFHI.addAtom (tauList(j)(i));
      }
   }
   tauFHI.endCreation ();
   return tauFHI;
}

void SxVibrations::getMasses ()
{
   int nAtoms = 0;
   int nSpecies = 0;
   int is, ia, i;
   SxVector3<Double> massVec;

   cout << "getting element table ..." <<  endl; fflush(stdout);
   cout << "ready" << endl; fflush (stdout);

   if (inputStructure != SxString("empty")) {
      SxList<SxList<SxVector3<Double> > > tauList;
      SxBinIO io;
      io.open (inputStructure, SxBinIO::ASCII_READ_ONLY);
      tauList = io.loadStructureFHI98 ();
      nSpecies = int(tauList.getSize ());

      SxList<SxString> speciesNames;
      io.open (inputStructure, SxBinIO::ASCII_READ_ONLY);
      speciesNames = io.loadSpeciesFHI98 ();
      io.close ();

      cout << "Getting atomic weights " << endl;
      masses = SxAtomicStructure ();
      masses.startCreation ();

      for (is = 0; is < nSpecies; is++) {
         nAtoms = int(tauList(is).getSize());
         for (ia = 0; ia < nAtoms; ia++) {
            for (i = 0; i < 3; i++) {
               massVec(i) = elemDB.getAtomicWeight(speciesNames(is));
            }
            masses.addAtom (massVec);
         }
      }
      masses.endCreation ();

   } else {
      masses = SxAtomicStructure ();
      masses.startCreation ();

      for (is = 0; is < speciesData.getNSpecies (); is++) {
         nAtoms = tau.getNAtoms (is);
         for (ia = 0; ia < nAtoms; ia++) {
            double mass = speciesData.ionicMass(is);
            const SxString &element = speciesData.chemName(is);
            if (speciesData.ionicMass(is) < 0.  && element.getSize () > 0)  {
               // fallback for uninitialized masses
               mass = elemDB.getAtomicWeight (element);
               cout << "Warning: ionic mass for species " << element
                    << " is taken from database as " << mass << endl;
            }
            for (i = 0; i < 3; i++) {
               massVec(i) = mass;
            }
            masses.addAtom (massVec);
         }
      }
      masses.endCreation ();
   }
}

void SxVibrations::setElementsToZero
(SxMatrix<Double> *matrix, double threshold)
{
   int i, j;

   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         if (fabs((*matrix)(i, j)) < threshold) (*matrix)(i, j) = 0.;
      }
   }

}

void SxVibrations::symmetrizeHessian (SxMatrix<Double> *matrix)
{
   int i, j;
   double helpi;

   for (i = 0; i < nDoF; i++) {
      for (j =0; j < nDoF; j++) {
         helpi = ((*matrix)(i, j) + (*matrix)(j, i))/2.;
         (*matrix)(i, j) = helpi;
         (*matrix)(j, i) = helpi;
      }
   }
}


#else /* SX_STANDALONE */

#include <SxCLI.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();

   SxVibrations sxv;

   SxCLI cli (argc, argv);
   sxv.inputFile = cli.option ("-i", "file", "sxvibrations input file")
                  .toString ("vibrations.sx");
   int tGroup = cli.newGroup ("template mode");
   cli.excludeGroup (cli.generalGroup);
   SxString outputFile = cli.option ("-o", "file", "sxvibrations template file name")
                         .toString ("vibrations.sx");
   cli.option ("--template", "generate sxvibrations input file template").toBool ();
   cli.finalize ();
   
   if (cli.groupAvailable (tGroup))  {
      sxv.flushTemplate (outputFile);
      exit (0);
   }

   sxv.parseInputFile ();
   sxv.validateInput ();
   sxv.printParameters ();
   sxv.initSystem ();
   sxv.manipulateInput ();
   sxv.printOutput ();
   return 0;
}
#endif /* SX_STANDALONE */
