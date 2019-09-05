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
#include <SxGrid.h>
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
   SxMatrix3<Int> repMatrix;

   int simpleRep = cli.newGroup ("simple repetition");

   SxVector3<Int> diag (cli.option ("-r|--repeat","3 numbers",
                                     "repetition along current cell vectors")
                            .toIntList3 ("x,"));

   int trueRep = cli.newGroup ("3dimensional repetition");
   cli.excludeGroup (simpleRep);
   SxVector3<Int> col1 (cli.option ("--a1","vector",
                                     "new a1 in relative coordinates")
                         .toIntList3 ()),
                  col2 (cli.option ("--a2","vector",
                                     "new a2 in relative coordinates")
                         .toIntList3 ()),
                  col3 (cli.option ("--a3","vector",
                                     "new a3 in relative coordinates")
                         .toIntList3 ());
   if (!cli.error)  {
      if (cli.groupAvailable (simpleRep))  {
         if (diag.product () <= 0)  {
            cout << "Illegal repetition factors " << diag << endl;
            cli.setError ();
         }
      } else if (cli.groupAvailable (trueRep))  {
         repMatrix = SxMatrix3<Int> (col1, col2, col3).transpose ();
         cout << repMatrix << endl;
         if (repMatrix.determinant () == 0)  {
            cout << "Illegal repetition matrix with zero determinant." << endl;
            cli.setError ();
         }
      } else {
         cout << "no repetition given!" << endl;
         cli.setError ();
      }
   }
   cli.setGroup (cli.generalGroup);
   
   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

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

   bool outSxb = cli.option ("--outsxb", "output file is binary "
                              "rather than S/PHI/nX input file").toBool ();
   bool labels = cli.option ("-l|--labels", "transfer labels").toBool ();
   bool movable = cli.option ("-m|--movable", "transfer movable").toBool ();
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
      if (movable)  {
         SxPtr<SxStickyFilter> stickyFilter;
         try {
            stickyFilter
               = SxPtr<SxStickyFilter>::create(tree->getGroup ("structure"));
         } catch (const SxException &e)  {
            e.print ();
            SX_EXIT;
         }
         // --- validate the sticky filter
         SX_CHECK (structure.cell.symGroupPtr);
         const SxSymGroup &S = *structure.cell.symGroupPtr;
         SxVector<Int> equivalentIdx (structure.getNAtoms(), 0);
         int           nSym = S.getSize ();
         cout << "| Validating sticky filter ....\n\n";
         for (int iSym = 0; iSym < nSym; iSym++) {
            structure.isEqual (S(iSym) ^ structure, &equivalentIdx);
            stickyFilter->validate (S(iSym).rot, equivalentIdx);
         }
         structure.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
                                          stickyFilter->getStickyArray ());
      }
   }

   SxAtomicStructure newStructure 
      = cli.groupAvailable (simpleRep) ? structure.repeat (diag)
                                       : structure.repeat (repMatrix);

   if (wrap) newStructure %= newStructure.cell;

   cout << "nSym = " << newStructure.cell.symGroupPtr->getSize () << endl;

   if (labels || movable)  {
      SxGrid grid(structure,10);
      SxConstPtr<SxAtomInfo> info = structure.match (grid, newStructure);
      newStructure.replaceInfo (info);
   }
   if (labels && structure.hasLabels ())  {
      const SxArray<SxString> &oldLabels = structure.getLabels ();
      SxArray<SxString> newLabels(newStructure.getNAtoms ());
      SX_CHECK (newStructure.atomInfo->parentMap.getSize () > 0);
      SX_LOOP(iNew)
         newLabels(iNew) = oldLabels(newStructure.atomInfo->parentMap(iNew));
      newStructure.atomInfo->meta.update (SxAtomicStructure::Labels, newLabels);
   }
   if (movable)  {
      SxConstPtr<SxAtomInfo> &info = newStructure.atomInfo;
      SX_CHECK (info->parentMap.getSize () > 0);
      SxArray<SxVector3<Int> > &oldMovable
         = structure.atomInfo->meta.get (SxAtomicStructure::StickyFilter);
      SxArray<SxVector3<Int> > newMovable(newStructure.getNAtoms ());
      SX_LOOP(iNew)
         newMovable(iNew) = oldMovable(info->parentMap(iNew));
      newStructure.atomInfo->meta.update (SxAtomicStructure::StickyFilter,
                                          newMovable);
   }

   // --- output
   if (outSxb)  {
      newStructure.write (outFile);
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

      newStructure.fprint (file);
         
      if (file != stdout) fclose (file);

   }
}

