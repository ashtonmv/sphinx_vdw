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
#include <SxStickyFilter.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on allows to shift all atoms within the structure.";
   cli.authors = "C. Freysoldt";

   // ---  definition
   SxVector3<Double> trans(cli.option ("--by|--vector","vector",
                                       "translation vector").toList3 ());
   
   bool relative = cli.option ("-r|--relative", "assume relative coordinates")
                   .toBool ();

   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   int binGroup = cli.newGroup ("binary files");
   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();
   cli.newGroup ("sx files");
   cli.excludeGroup (binGroup);
   bool keepMovable
      = cli.option("-m|--keep-movable","retain mobility of atoms as "
                   "it is defined in the input structure")
        .toBool ();
   cli.setGroup(SxCLI::generalGroup);
   
   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString (sxbFile ? "waves.sxb" : "input.sx");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";
   
   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   bool outSxb = cli.option ("--outsxb", "output file is binary "
                              "rather than S/PHI/nX input file").toBool ();
   cli.last ().defaultValue = "default: same format as input file";
   if (sxbFile) outSxb = true;
   if (outSxb && outFile.getSize () == 0)
      outFile = "structure.sxb";

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure;
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
      if (keepMovable)
         structure.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
               SxStickyFilter(tree->getGroup("structure")).getStickyArray ());
   }

   
   // --- Relative coordinates
   if (! structure.isPeriodic () && relative)  {
      cout << "Nonperiodic structures do not have relative coordinates.\n";
      cout << "--relative flag can't be used for the structure from." << endl;
      cout << "file '" << inFile << "'." << endl;
      SX_QUIT;
   }
   
   if (relative) structure.cell.changeToCar (&trans);
   
   // shift structure
   structure += trans;

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

