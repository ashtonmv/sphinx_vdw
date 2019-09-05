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
#include <SxXYZ.h>
#include <SxCLI.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxSimpleParser.h>

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
   SxString outFile = cli.option("-o","file","XYZ output file to be written")
                      .toString("input.xyz");
   bool withAMat = cli.option ("-a", "prepend lattice vectors (a1, a2, a3)")
                   .toBool();
   bool trajectory = cli.option ("--trajectory", "file is a trajectory")
                     .toBool ();
   int interpolate =  cli.option ("--interpolate", "number of points", "trajectory interpolation with additional points").toInt (0,1);

   int simpleRep = cli.newGroup ("simple repetition");
   SxVector3<Int> repeat(0,0,0);
   cli.option ("-r|--repeat","3 numbers",
               "repetition along current cell vectors");
   if (cli.last ().exists ())
      repeat = SxVector3<Int> (cli.last ().toIntList3 ("x,"));

   int trueRep = cli.newGroup ("3dimensional repetition");
   SxMatrix3<Int> repMatrix;
   repMatrix.set (0);
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
         if (repeat.product () <= 0)  {
            cout << "Illegal repetition factors " << repeat << endl;
            cli.setError ();
         }
      } else if (cli.groupAvailable (trueRep))  {
         repMatrix = SxMatrix3<Int> (col1, col2, col3).transpose ();
         cout << repMatrix << endl;
         if (repMatrix.determinant () == 0)  {
            cout << "Illegal repetition matrix with zero determinant." << endl;
            cli.setError ();
         }
      }
   }
   //cli.setGroup (cli.generalGroup);
   bool doWrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   cli.finalize ();

   // --- read from input file
   initSPHInXMath ();

   SxAtomicStructure structure;
   SxArray<SxString> chemName;


   if (trajectory && !sxbFile)  {
      FILE *fp =sxfopen (outFile, "w");
      try {
         SxParser parser;
         SxParser::Table table = parser.read (inFile, "std/forces.std");

         chemName = SxSpeciesData::getElements (&*table);
         SX_LOOP(is) chemName(is) = chemName(is).toUpper ();

         int it = 0;
         SxAtomicStructure wrap, oldStructure;
         SYMBOLPARSE(&*table)  {
            FOREACH_SYMBOLGROUP("structure")  {
               structure = SxAtomicStructure(SYMBOLGROUP_TABLE);
               if (repeat.product () > 0)
                  structure = structure.repeat (repeat);
               else if (repMatrix.determinant () != 0)
                  structure = structure.repeat (repMatrix);
               if (doWrap)  {
                  if (it == 0)
                     wrap = (structure % structure.cell) - structure;
                  structure += wrap;
               }
               SX_CHECK(structure.getNSpecies () == chemName.getSize (),
                        structure.getNSpecies (), chemName.getSize ());
               if (interpolate > 0 && oldStructure.getNAtoms () > 0)  {
                  for (int ip = 1; ip <= interpolate; ++ip)  {
                     double x = double(ip)/(interpolate+1);
                     fprintf (fp, "%d\nStep %g\n", structure.getNAtoms (),
                              it + x);
                     SX_LOOP2(is,ia)  {
                        const Coord &newPos = structure.getAtom (is,ia);
                        Coord pos = oldStructure.getAtom (is,ia);
                        pos += x * (newPos - pos);
                        fprintf (fp, "%s %.8f %.8f %.8f\n", chemName(is).ascii (),
                                 pos(0) / A2B, pos(1) / A2B, pos(2) / A2B);
                     }
                  }
               }
               oldStructure = structure;
               fprintf (fp, "%d\nStep %d\n", structure.getNAtoms (), ++it);
               SX_LOOP2(is,ia)  {
                  const Coord &pos = structure.getAtom (is,ia);
                  fprintf (fp, "%s %.8f %.8f %.8f\n", chemName(is).ascii (),
                           pos(0) / A2B, pos(1) / A2B, pos(2) / A2B);
               }
            }
         }
      } catch (SxException e) {
         e.print ();
         SX_QUIT;
      }
      fclose (fp);
      return 0;
   }

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
   if (repeat.product () > 0)
      structure = structure.repeat (repeat);
   else if (repMatrix.determinant () != 0)
      structure = structure.repeat (repMatrix);
   if (doWrap) structure %= structure.cell;

   SxXYZ (outFile).write (structure, chemName, withAMat);
   return 0;
}

