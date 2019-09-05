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
      "This add-on computes the center of mass for some atoms.";
   cli.authors = "C. Freysoldt";

   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString ("input.sx");
   
   SxArray<int> atomList = cli.option ("-n|--atoms", "index list", 
                   "atoms to rotate").toIdxList();
   cli.finalize ();
   if (atomList.getSize () < 2) {
      cout << "Nothing to do..." << endl;
      return 0;
   }

   // --- read input
   SxAtomicStructure structure;
   SxParser parser;
   SxConstPtr<SxSymbolTable> tree;
   tree = parser.read (inFile, "std/structure.std");
   structure = SxAtomicStructure (&*tree);
   structure.readElements (&*tree);

   Coord centerOfMass;
   SxString label;
   SX_LOOP(iList)  {
      int iTlAtom = atomList(iList);
      if (iTlAtom < 0 || iTlAtom >= structure.getNAtoms ())  {
         cout << "Illegal atom number " << (iTlAtom+1) << endl;
         SX_QUIT;
      }
      if (iList == 0)  {
         if (structure.hasLabels ())
            label = structure.getLabels ()(iTlAtom);
         centerOfMass = structure.constRef(iTlAtom);
      } else  {
         if (structure.hasLabels ())  {
            const SxString &newLabel = structure.getLabels ()(iTlAtom);
            if (newLabel != label) label = "";
         }
         Coord relPos = structure.constRef(iTlAtom) - centerOfMass;
         // find nearest image
         structure.cell.map (&relPos, SxCell::WignerSeitz);
         cout << relPos << endl;
         // com := 1/(N+1) ( N * com + (com + relPos))
         centerOfMass += 1./int(iList + 1) * relPos;
      }
   }
   cout << "Center of mass:" << centerOfMass << endl;
   sxprintf ("      atom { coords = [% 12.8f, % 12.8f, % 12.8f];",
             centerOfMass(0), centerOfMass(1), centerOfMass(2));
   if (label.getSize () > 0) cout << " label = \"" << label << "\";";
   sxprintf (" }\n");
   cout << "atoms:" << endl;
   SX_LOOP(iList)  {
      int iTlAtom = atomList(iList);
      Coord relPos = structure.constRef(iTlAtom) - centerOfMass;
      // find nearest image
      structure.cell.map (&relPos, SxCell::WignerSeitz);
      sxprintf ("      // d=%.6f\n", relPos.norm ()); 
      sxprintf ("      atom { coords = [% 12.8f, % 12.8f, % 12.8f]\n"
                "                    + [% 12.8f, % 12.8f, % 12.8f]; }\n",
                centerOfMass(0), centerOfMass(1), centerOfMass(2),
                relPos(0), relPos(1), relPos(2));
   }
}

