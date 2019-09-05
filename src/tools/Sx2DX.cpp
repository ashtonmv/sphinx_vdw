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

#include <SxUtil.h>
#include <SxCLI.h>
#include <SxAtomicStructure.h>
#include <SxElemDB.h>
#include <SxSpeciesData.h>
#define a0 0.5291772083

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   
   // --- input and output files
   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();

   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   SxString outFile = cli.option("-o","file","DX output file to be written")
                      .toString("input.dx");
   SxVector3<Int> mult
              (cli.option ("-m","multiply","multiply charge density")
             .toIntList3 ("x,"));

   double delta = cli.option("-d","adjustment for bonds","additive constant for bonds in bohr")
                      .toDouble(0.5);

   cli.finalize ();
   
   // --- read from input file
   initSPHInXMath ();

   SxAtomicStructure structure;
   SxArray<SxString> chemName;
   SxElemDB elemDB;

   if (sxbFile)  {

      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         SxString chemNameList;
         io.read("chemNames", &chemNameList);
         chemName = chemNameList.tokenize (',');
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }

   } else  {
   
      try {
         SxParser parser;
         SxParser::Table table = parser.read (inFile, "std/structure.std");
   
         structure = SxAtomicStructure(&*table);
         chemName = SxSpeciesData::getElements (&*table);
      } catch (SxException e) {
         e.print ();
         SX_QUIT;
      }
   }

cout << inFile << " -> " << outFile << endl; 

   // --- multiply structure
   structure = structure.repeat(mult);

   // --- write down structure
   ofstream file;
   file.open(outFile.ascii ());
   int nTlAtoms = structure.nTlAtoms;
   file << "object 1 class array type float rank 1 shape 3 items " << nTlAtoms 
        << " data follows" << endl;
   for (int i = 0; i < nTlAtoms; i++)   {
      file << fixed << structure(i)(0) * a0 
           << " " << structure(i)(1) * a0
           << " " << structure(i)(2) * a0
           << endl;
   }
   file << endl;
   file << "object 2 class array type float rank 0 items " << nTlAtoms
        << " data follows" << endl; 
   int nSpecies = structure.getNSpecies ();
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int nAtom = structure.getNAtoms(iSpecies);
      int idx = elemDB.getIdx(chemName(iSpecies));
      double atRadius = elemDB.getAtomicRadius(idx); 
      for (int iAtom = 0; iAtom < nAtom; iAtom++)   {
         file << fixed << atRadius << endl;
      }
   }
   file << endl;

   file << "object 3 class array type float rank 1 shape 3 item " << nTlAtoms
        << " data follows" << endl;
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int nAtom = structure.getNAtoms(iSpecies);
      int idx = elemDB.getIdx(chemName(iSpecies));
      double colR = elemDB.getRed(idx); 
      double colG = elemDB.getGreen(idx); 
      double colB = elemDB.getBlue(idx); 
      for (int iAtom = 0; iAtom < nAtom; iAtom++)   {
         file << fixed << colR << " " << colG << " " << colB << endl;
      }
   }
   file << endl;

   // Calculate Bonds
   int nBonds = 0;
   SxList<SxArray<int> > bondList;
   for(int i = 0; i < nTlAtoms-1; i++)   {
      for (int j = i+1; j < nTlAtoms; j++)   {
         double dist = 
            (structure.coords.colRef(i) - structure.coords.colRef(j)).norm ();
         int iSpecies = structure.getISpecies(i);
         int jSpecies = structure.getISpecies(j);
         int idx = elemDB.getIdx(chemName(iSpecies));
         int jdx = elemDB.getIdx(chemName(jSpecies));
         double covDist = 
            elemDB.getCovalentRadius(idx) + elemDB.getCovalentRadius(jdx);
         if (dist < covDist + delta)   {
            SxArray<int> fromTo (2);
            fromTo(0) = i;
            fromTo(1) = j;
            bondList << fromTo;
            nBonds++;
         }
      }
   }

   if (nBonds > 0)   {
      file << "object 4 class array type int rank 1 shape 2 items " << nBonds 
         << " data follows" << endl;
      SxList<SxArray<int> >::Iterator it;
      for(it = bondList.begin(); it != bondList.end(); it++)   {
         file << (*it)(0) << " " << (*it)(1) << endl;
      }
      file << " attribute \"element type\" string \"lines\"" << endl;
      file << " attribute \"ref\" string \"positions\"" << endl;
      file << endl;
   }

   file << "object \"charge\" class field" << endl;
   file << " component \"positions\" value 1" << endl;
   file << " component \"data\" value 2" << endl;
   file << " component \"colors\" value 3" << endl;
   if (nBonds > 0)   {
      file << " component \"connections\" value 4" << endl;
   }
   file << "end" << endl;
   file.close();

   return 0;
}

