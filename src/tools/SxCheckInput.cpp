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

#include <SxCLI.h>
#include <SxUtil.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxNeighbors.h>
#include <SxConstants.h>
#include <SxElemDB.h>

//#define a0 0.5291772083 
#define a0 1.

// Functionheader
double getMaximalCovRadius (const SxArray<SxString> chemNames); 

//global variable
SxElemDB elemDB;

int main (int argc, char **argv)
{
   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on checks the input.sx file for inconsitencies! "
      "Up to now there is a dublet checker and a interatomic distance "
      "checker available.\n";
      
   cli.authors = "B. Lange";

   cli.finalize ();


   // --- we need the atomic structure out of input.sx
   SxParser parser;
   SxParser::Table table = parser.read ("input.sx");
   SxAtomicStructure structure (&*table);
   SxCell cell (&*table);
   SxArray<SxString> chemNames = SxSpeciesData(&*table).chemName;
   
   // --- get Number of total atoms
   int numberOfAtoms = structure.getNAtoms ();

   // --- map Atoms back into cell
   SxAtomicStructure mappedStructure = structure.getMapped ();

   // --- Check for dubletes
   cout << "Checking for dublets...";
   bool ok = true;
   for (int i = 0; i<numberOfAtoms-1; i++)   {

      for (int j = i+1; j<numberOfAtoms; j++)   {

         if (cell.getMapped (structure.getAtom(i) - structure.getAtom(j),
                             SxCell::Origin).norm () < 1e-4)
         {
            if (ok) cout << endl;

            cout << "WARNING !!!" << endl;
            cout << "Atomic Distance " << i << " - " << j << " is very small!" 
                 << endl;
            cout << "Atom seems to be a dublet!" << endl;
            ok = false;
         
         }
      }
   }
   if (ok) cout << "ok" << endl;

   // --- Check interatomic distances
   SxVector3<Int> gMesh 
      = SxGrid::suggestMesh (structure.cell, numberOfAtoms / 10 + 1);
   SxGrid grid (structure, gMesh);
   SxNeighbors neighbors;

   
   // --- Calculate greatest covradius in cell
   double maximalCovRadius = getMaximalCovRadius(chemNames);
   
   // --- Check cellsize
   cout << "Checking cell size." << endl;
   double norm0 = structure.cell(0).norm ();
   double norm1 = structure.cell(1).norm ();
   double norm2 = structure.cell(2).norm ();
   if (maximalCovRadius/a0 > norm0) {
      cout << "WARNING !!!" << endl;
      cout << "Cell Vector 1 is small!"<< endl;
      cout << "Norm is " << norm0 <<", should be at least " 
         << maximalCovRadius/a0 << "!\n"; 
      cout << "Check cell parameters" << endl;
   }
   if (maximalCovRadius/a0 > norm1) {
      cout << "WARNING !!!" << endl;
      cout << "Cell Vector 2 is small!"<< endl; 
      cout << "Norm is " << norm1 <<", should be at least " 
         << maximalCovRadius/a0 << "!\n"; 
      cout << "Check cell parameters" << endl;
   }
   if (maximalCovRadius/a0 > norm2) {
      cout << "WARNING !!!" << endl;
      cout << "Cell Vector 3 is small!"<< endl; 
      cout << "Norm is " << norm2 <<", should be at least " 
         << maximalCovRadius/a0 << "!\n"; 
      cout << "Check cell parameters" << endl;
   }



   // --- Create a cov. Sphere around each atom and check if 
   //     other atoms are inside this sphere
   ok = true;
   cout << "Checking for close atoms...";

   for (int i = 0; i < numberOfAtoms; i++)  {
      
      int mySpecies = structure.getISpecies(i);
      int myIdx = elemDB.getIdx(chemNames(mySpecies));
      double myCovRad = elemDB.getCovalentRadius(myIdx);  
      
      // a0 is due to convert Angstroem in Bohr
      double critDist = (myCovRad + maximalCovRadius) / a0;
      neighbors.compute (grid, structure, structure(i), critDist, 
                         SxNeighbors::StoreAll);
      
      for (int j = 0; j < neighbors.getSize (); j++)  {
      
         int yourSpecies = structure.getISpecies(neighbors.idx(j));
         int yourIdx = elemDB.getIdx(chemNames(yourSpecies));
         if (neighbors.idx(j) < i) continue; // has been checked already
         double yourCovRad = elemDB.getCovalentRadius(yourIdx);
      
         double errorDist = 0.7 * ((myCovRad + yourCovRad)/ a0);
         Coord whichCell = structure.cell.carToRel (neighbors.absPositions(j)
                                             - structure(neighbors.idx(j)));
         if (sqrt(neighbors.distSqr(j)) < errorDist) {
            if (ok) cout << endl;
            cout << "WARNING !!!" << endl;
            cout << "Atomic Distance " << (i+1) << " - " 
                 << (neighbors.idx(j)+1);
            if (whichCell.normSqr () > 1e-12)
               cout << "@" << neighbors.absPositions (j) << " in cell ";
            cout << " is small!" << endl;
            cout << "Distance is: " << sqrt(neighbors.distSqr(j)) 
                 << " bohr." << endl;
            cout << "Covalent radii are: "
                 << chemNames(mySpecies) << "=" << myCovRad << " / " 
                 << chemNames(yourSpecies) << "=" << yourCovRad
                 << " bohr." << endl;
            cout << "Warning distance is 70% of the covalent radii sum = "
                 << errorDist << " bohr." << endl;
            cout << "Check atom positions !" << endl;
            ok=false;
         }
      }
   }
   if (ok) cout << "ok" << endl;

   return 0;
}

double getMaximalCovRadius (const SxArray<SxString> chemNames)  {

   double result = 0.0;

   for (int i = 0; i < chemNames.getSize(); i++)  {
      int idx = elemDB.getIdx(chemNames(i));
      double covRad = elemDB.getCovalentRadius(idx);
      if(result < covRad)  result = covRad;
   }

   return result;
}
