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
#include <SxAims.h>
#include <SxCLI.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

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
   SxString outFile = cli.option("-o","file","FHIaims structure file to "
                                 "be written")
                      .toString("geometry.in");

   cli.finalize ();
   
   // --- read from input file
   initSPHInXMath ();

   SxAtomicStructure structure;
   SxArray<SxString> chemName;

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
   
   SxAims aims (outFile);
   aims.write(structure, chemName);
   
   return 0;
}

