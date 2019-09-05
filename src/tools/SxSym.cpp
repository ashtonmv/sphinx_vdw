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
#include <SxSymFinder.h>
#include <SxString.h>
#include <SxParser.h>
#include <SxUniqueList.h>

int main (int argc, char** argv)
{
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   SxString inputFile
      = cli.option ("-i|--input", "file", "SPHInX input file")
        .toString ("input.sx");
   int minSymmetries 
      = cli.option ("-n", "number", "suppress output of translations that "
                    "produce less symmetries than this number")
        .toInt (0);
   SxString outputFile
      = cli.option ("-o|--output", "file", "save a high-symmetry structure")
        .toString ("");

   double epsStruct
      = cli.option ("--eps_struct", "number", "epsilon for structure "
                    "differences. Higher this, if affine symmetries are not "
                    "found.")
        .toDouble (1e-4);

   SxLinEquation::epsZero
      = cli.option ("--eps_lin", "number", "epsilon for linear dependency. "
                    "Higher this, if affine symmetries are found, but can't "
                    "be made symmorphic.")
        .toDouble (SxLinEquation::epsZero, 0);

   cli.authors = "Christoph Freysoldt";
   cli.finalize ();

   SxParser parser;
   SxParser::Table table = parser.read (inputFile, "std/structure.std");
   SxAtomicStructure structure (&*table);

   structure.epsEqual = epsStruct;
   //structure.print ();
   
   SxSymFinder symFinder;

   symFinder.compute (structure);

   // --- print result
   SxList<SxRotation>::Iterator rotIt;
   typedef SxUniqueList<SxRotation,SxRotation> RotationList;
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator it;
   SxMap<SxLinEquation, RotationList, SxNull>::Iterator itEnd;

   it = symFinder.equations.begin ();
   itEnd = symFinder.equations.end ();
   
   for (; it != itEnd; it++)
   {
      SX_CHECK (it.getKey().isSoluble ());
      
      // print only shifts with at least minSymmetries symmetries
      if (it.getValue().getSize () < minSymmetries) continue;
      
      cout << "The following " << it.getValue().getSize ();
      cout << " symmetries are available for" << endl;
      cout << it.getKey().getName ("translation") << endl;
      for (rotIt = it.getValue().begin ();
           rotIt != it.getValue().end (); ++rotIt)
      {
         cout << SxRotation::getName (*rotIt) << endl;
         // cout << *rotIt << endl;
      }
      cout << "---" << endl;
   }

   if (outputFile.getSize () != 0) {
      FILE *fp = fopen (outputFile.ascii (), "w");
      if (fp == NULL)  {
         cout << "Can't open '" << outputFile << "'" << endl;
         SX_EXIT;
      }
      // move atoms to high-symmetry position
      structure += symFinder.getHighSymShift ();
      structure %= structure.cell;
      
      structure.fprint (fp);
      fclose (fp);
   }


}

