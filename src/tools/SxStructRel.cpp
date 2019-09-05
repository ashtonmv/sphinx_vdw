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

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on prints a structure in relative coordinates.";
   cli.authors = "C. Freysoldt";

   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();
   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";
   
   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure;
   SxArray<SxString> speciesInclude; 
   if (sxbFile)  {
      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (inFile, "std/structure.std");
      structure = SxAtomicStructure (&*tree);
      structure.readElements (&*tree);
   }

   
   // --- Relative coordinates
   if (! structure.isPeriodic ())  {
      cout << "Structure from file '" << inFile << "' is not periodic." << endl;
      cout << "Nonperiodic structures do not have relative coordinates.\n";
      SX_QUIT;
   }
   
   // --- output
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

   structure.fprintRel (file);
      
   if (file != stdout) fclose (file);

}

