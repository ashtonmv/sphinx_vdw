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
#include <SxPoscar.h>
#include <SxCLI.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   
   // --- input and output files
   SxString inFile 
      = cli.option("-i","file","SPHInX structure file to read from (bohr)")
        .toString("struc.sx");
   SxString outFile 
      = cli.option("-o","file","POSCAR output file to be written (angstrom)")
        .toString("POSCAR");
   cli.finalize ();
   
   // --- read from input file
   initSPHInXMath ();
   SxParser parser;
   SxParser::Table table = parser.read (inFile.ascii (), "std/structure.std");
   SxAtomicStructure structure = SxAtomicStructure(&*table);
   SxArray<SxString> chemName = SxSpeciesData::getElements (&*table);

   SxPoscar (outFile).write (structure, chemName);
   return 0;
}

