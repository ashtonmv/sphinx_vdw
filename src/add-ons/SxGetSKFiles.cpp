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
#include <stdio.h>
#include <SxString.h>

#ifndef SX_STANDALONE
#ifdef WITH_TB

#include <SxGetSKFiles.h>

#include <SxError.h>
#include <iostream>

#include <SxTBSKData.h>

SxGetSKFiles::SxGetSKFiles (const SxSymbolTable *tableIn, 
                            double maxDIn, double stepIn, double sigmaIn)
   : structure (&*tableIn),
     SKData (structure, &*tableIn, sigmaIn) 
{

   lOrbs         = SKData.getlOrbs ();
   elementName   = SKData.getElementName ();
   chemName      = SKData.getChemName ();
   valenceCharge = SKData.getValenceCharge ();
   nSpecies      = SKData.nSpecies;
   maxDist       = maxDIn;
   strStep       = stepIn; 
   nStrSteps     = int (ceil (maxDist/strStep));
   sigma         = sigmaIn;

}

SxGetSKFiles::~SxGetSKFiles ()
{
   // empty
}

void SxGetSKFiles::setInitialStr ()
{
   // --- check number of atoms    
   if (structure.nTlAtoms != 2)  {
      cout << " use only two atoms in your input\n";
      SX_QUIT;
   }

   // --- check the cell dimensions    
   if ( structure.cell(0)(0) <  2.*maxDist || 
        structure.cell(0)(1) != 0  || 
        structure.cell(0)(2) != 0  || 

        structure.cell(1)(0) != 0  || 
       (structure.cell(1)(1) <  12. && structure.cell(1)(1) < 3.*sigma)  || 
        structure.cell(1)(2) != 0  || 

        structure.cell(2)(0) != 0  || 
        structure.cell(2)(1) != 0  || 
       (structure.cell(2)(2) <  12. && structure.cell(2)(2) < 3.*sigma) )  { 
      cout << " \n invalid cell for this application.\n";
      cout << " use a cell whose axes are along xyz axes respectively"
           << " and AT LEAST:\n"
           << " a1 = [2d  , 0.0        , 0.0        ]  \n"
           << " a2 = [0.0 , min{12, 3s}, 0.0        ]  \n" 
           << " a3 = [0.0 , 0.0        , min{12, 3s}]  (all in Bohr)\n"
           << " where: d = the maximum distance, s = sigma \n";
      SX_QUIT;
   }
         
   // --- initialize atomic positions
   structure.ref (0) = SxVector3<Double> (0, 0, 0);
   structure.ref (1) = SxVector3<Double> (0, 0, 0); 
    
}

void SxGetSKFiles::calculate (const SxSymbolTable *tableIn)
{
   // increase 2nd Atom x-component
   structure.ref (1)(0) += strStep;

   SKData.updateBasis (structure);

   SKData.tightBindingInit (&*tableIn);

   H = SKData.getHam ();
   S = SKData.getOvl ();

}

void SxGetSKFiles::setTo (const double num)
{
   // increase 2nd Atom x-component
   structure.ref (1)(0) += strStep;    
   
   SKData.setTo (num);

   H = SKData.getHam ();
   S = SKData.getOvl ();

}


void SxGetSKFiles::calcEigVals (const SxSymbolTable *tableIn)
{
   // put the 2nd atom as far as possible from the 1st atom
   // on the x-axis. Note: Calculating atomic eigenvaluse this way 
   // is dangerous if maxDist is small!
   structure.ref (1)(0) = maxDist/2.0;

   // update structure and basis in SKData
   SKData.updateBasis (structure);

   SKData.calcEigVals (&*tableIn);
   eigVals = SKData.getEigVals ();

   // return the structure to as before
   structure.ref (1)(0) = 0.;

}

void SxGetSKFiles::printHeaders (FILE *fp) 
{

   fprintf (fp," delta   nPoints   lMax1   lMax2\n");

   if (nSpecies == 2)  {

      fprintf (fp,"%6.3f    %d        %d      %d\n", strStep, nStrSteps, 
                    lOrbs(0).getSize()-1, lOrbs(1).getSize()-1);  

   }  else if (nSpecies == 1)  {

      fprintf (fp,"%6.3f     %d       %d      %d\n", strStep, nStrSteps, 
                    lOrbs(0).getSize()-1, lOrbs(0).getSize()-1);  
      fprintf (fp,"---------------+----------------------------+-----------\n");
      fprintf (fp," valenceCharge | eigenvalues                | Habbard_U \n");
      fprintf (fp,"    %d         ", int (valenceCharge(0)));

      for (int iOrb = 0; iOrb < lOrbs(0).getSize(); iOrb++)
         fprintf (fp, " %10.7f", eigVals(0)(iOrb)); 

      fprintf (fp, "        %5.2f   %5.2f ", 0., 0.); // TODO
      fprintf (fp, " // wrong hubbard_U \n");  

   }  else  {
      SX_EXIT;
   }
   
   fprintf (fp,"-------+------------------------------------------------\n");
   fprintf (fp," dist  | Hamiltonian and Overlap elements \n");
   fprintf (fp,"-------+------------------------------------------------\n");

}

void SxGetSKFiles::updateFile (FILE *fp) 
{

   fprintf (fp,"%7.4f",structure(1)(0));

   // the order is as needed by SxTBAtomAtom
   fprintf (fp,"%20.12le %20.12le %20.12le %20.12le %20.12le",
                H(0)(0,0), H(0)(0,1),H(0)(1,0), H(1)(1,1), H(0)(1,1));

   fprintf (fp,"%20.12le %20.12le %20.12le %20.12le %20.12le\n",
                S(0)(0,0), S(0)(0,1),S(0)(1,0), S(1)(1,1), S(0)(1,1));

}

#endif /* WITH_TB */

#else /* SX_STANDALONE */

#ifdef WITH_TB

#include <SxCLI.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   initSPHInXIO ();

   SxCLI cli (argc, argv);
   cli.preUsageMessage = " add-on to generate slater-koster tight-binding"
                         " data files";

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");
   double maxDist  = cli.option ("-d|--distance","number",
                     "maximum distance to reach between the\n"
                     "two atoms (Bohr)")
                     .toDouble (10.,0.4,20.);
   double strStep  = cli.option ("-g|--grid", "number", 
                     "grid spacing between structural steps (Bohr)")
                     .toDouble(0.02,0.001,1.0);
   double sigma    = cli.option ("-s|--sigma","number",
                     "the standard deviation of a Gaussian\n"
                     "distribution that is used to re-confine the\n"
                     "pseudo wavefunctions more. Thus, increasing\n"
                     "sigma reduces the effect of this Gaussian\n")
                     .toDouble(100,0.1,100);
   cli.finalize ();

   // --- print out input parameters
   cout << SX_SEPARATOR;
   cout << "| SxGetSKFiles \n";
   cout << SX_SEPARATOR;
   cout << "| options :\n"; 
   cout << "|    input file    : " << inFile << endl; 
   cout << "|    max. distance : " << maxDist << endl; 
   cout << "|    grid spacing  : " << strStep << endl; 
   cout << "|    sigma         : " << sigma   << endl; 
   cout << SX_SEPARATOR;

   // start timers
   initSPHInX();

   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   // call the constructor
   SxGetSKFiles getFiles (&*table, maxDist, strStep, sigma);

   // check the structure here, and initialize
   getFiles.setInitialStr ();
   
   const int nSpecies  = getFiles.nSpecies;
   const int nStrSteps = int (ceil (maxDist/strStep));
   int iStep;

   SxString fileName;
   FILE    *file;

   // TODO: do the calculations and save all elements and print once
   
   // --- open file
   if (nSpecies == 1)
      fileName = getFiles.chemName(0) + getFiles.chemName(0) + ".dat";
   else if (nSpecies == 2)
      fileName = getFiles.chemName(0) + getFiles.chemName(1) + ".dat";
   else
      SX_EXIT;
   file = fopen(fileName.ascii(),"w"); 

   // --- claculate eigenvalues and print headers
   getFiles.calcEigVals (&*table);
   getFiles.printHeaders (file);

   // --- start the structural loop
   // --- un-physical distances < 0.4 Bohr, set elements to a dummy number
   int nISteps = int (ceil(0.4/strStep)) - 1; 
   for (iStep = 0; iStep < nISteps; iStep++) {
      getFiles.setTo (0.1);   // dummy number
      getFiles.updateFile (file);
   }
   // --- now start the claculation loop for physical distances 
   for (iStep = 0; iStep < (nStrSteps - nISteps); iStep++) {
      getFiles.calculate (&*table);
      getFiles.updateFile (file);
   }
   // --- close file
   fclose (file); 

   cout << SX_SEPARATOR;
   cout << "| SxGetSKFiles exited normally.\n";
   cout << SX_SEPARATOR;

   return 0;
}
  

#else /* WITH_TB */

int main ()
{
   sxprintf ("ERROR: SPHInX compiled without TB support\n");
   sxprintf ("       sxgetskfiles not available.\n");
   return -1;
}
#endif /* WITH_TB */

#endif /* SX_STANDALONE */
