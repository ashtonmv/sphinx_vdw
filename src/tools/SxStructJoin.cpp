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
#include <SxConfig.h>
#include <SxBinIO.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>

template <class T>
SxList<T> &operator<< (SxList<T> &list, const SxArray<T> &array)
{
   for (int i = 0; i < array.getSize (); ++i)
      list.append(array(i));
   return list;
}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on allows to join the atoms of two structures. The final "
      "structure will be the first structure plus all the atoms of the second "
      "structure, with their cartesian coordinates. The origin of the second "
      "structure can be optionally shifted with respect to the first one.";
   cli.authors = "C. Freysoldt";

   // ---  definition
   SxVector3<Double> trans(0.,0.,0.);
   
   SxCLI::CliArg *opt = &cli.option ("--by|--vector","vector", 
                                     "translation vector");
   if (opt->exists ())  {
      trans = SxVector3<Double>(opt->toList3 ());
   }
   opt->required (false);
   
   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   SxString substrateFile 
      = cli.option ("-i|--substrate", "input file", "S/PHI/nX structure file "
                    "for first structure (determines cell)")
        .toString ("input.sx");
   
   SxString adsorbateFile 
      = cli.option ("-j|--join|--adsorbate", "input file", 
                    "S/PHI/nX structure file for second structure")
        .toString ();
   
   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   bool outSxb = cli.option ("--outsxb", "output file is binary "
                              "rather than S/PHI/nX input file").toBool ();
   if (outSxb && outFile.getSize () == 0)
      outFile = "structure.sxb";

   cli.finalize ();

   // --- read input
   SxAtomicStructure substrate, adsorbate;
   SxArray<SxString> chemName, adsorbChemName;
   try {
      // --- get substrate structure
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (substrateFile, "std/structure.std");
      substrate = SxAtomicStructure (&*tree);
      chemName = SxSpeciesData::getElements (&*tree);
   } catch (SxException e) {
      e.print ();
      SX_QUIT;
   }
   try {
      // --- get adsorbate structure
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (adsorbateFile, "std/structure.std");
      adsorbate = SxAtomicStructure (&*tree);
      adsorbChemName = SxSpeciesData::getElements (&*tree);
      // shift origin of adsorbate structure
      adsorbate += trans;
   } catch (SxException e) {
      e.print ();
      SX_QUIT;
   }

   // --- join structures
   int is, ia, joinedSpecies;
   int nSpecies = substrate.getNSpecies ();
   SxAtomicStructure structure(substrate, SxAtomicStructure::Copy);
   SxList<SxString> chemNames;
   chemNames << chemName;
   structure.startCreation ();
   for (is = 0; is < adsorbate.getNSpecies (); ++is)  {
      cout << "Joining " << adsorbChemName(is) << "..." << endl;
      joinedSpecies = int(chemNames.findPos (adsorbChemName(is)));
      if (joinedSpecies < 0)  {
         joinedSpecies = nSpecies++;
         structure.newSpecies ();
         chemNames << adsorbChemName(is);
      }
      for (ia = 0; ia < adsorbate.getNAtoms(is); ++ia)
         structure.addAtom (joinedSpecies,adsorbate(is,ia));
   }
   structure.endCreation ();
   structure.atomInfo->meta.attach (SxAtomicStructure::Elements,
                                    SxArray<SxString> (chemNames));

   if (wrap) structure %= structure.cell;

   // --- output
   if (outSxb)  {
      structure.write (outFile);
   } else {
      FILE *file;
      if (outFile.getSize () > 0)  {
         file = fopen(outFile.ascii(),"w");
         if (file == NULL)  {
            cout << "Can't open '" << outFile << "'." << endl;
            SX_EXIT;
         }
      } else {
         file = stdout;
      }

      structure.fprint (file);
         
      if (file != stdout) fclose (file);

   }
}
