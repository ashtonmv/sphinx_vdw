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
#include <SxSymFinder.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on allows to rotate a structure (cell, atoms).";
   cli.authors = "C. Freysoldt";

   // ---  definition
   SymMat rot;

   int explicitGroup = cli.newGroup ("Explicit rotation matrix");
   SxList<double> mElements = cli.option ("--matrix","matrix elements",
     "matrix elements row-wise, separated by colons ':'. Example:\n"
     "--matrix 1:0:0:0:0.707107:-0.707107:0:0.707107:0.707107\n"
     "would describe the matrix\n"
     "1    0         0       \n"
     "0    0.707107 -0.707107\n"
     "0    0.707107  0.707107\n").toDoubleList ();
   if (cli.groupAvailable (explicitGroup) && !cli.error)  {
      if (mElements.getSize () != 9)  {
         cout << "Number of matrix elements must be 9!" << endl;
         cli.setError ();
      } else {
         rot = SymMat (mElements);
         if (fabs(rot.determinant ()) < 1e-5)  {
            cout << "Rotation matrix has zero determinant!" << endl;
            cli.setError ();
         }
      }
   }
   int axisGroup = cli.newGroup ("Axis definition");
   cli.excludeGroup (explicitGroup);
   
   SxVector3<Double> axis (cli.option ("--axis", "vector", "rotation axis")
                            .toList3 () ); 
   
   SxString byString 
      = cli.option ("--by", "angle", "rotation angle (degrees) OR n/m")
        .toString ();

   int atomOnAxis = cli.option ("--center", "atom", "center of rotation axis")
                            .toInt (-1,-1); 

   double angle = 0.;
   if (cli.groupAvailable(axisGroup) && !cli.error)  {
      try {
         if (byString.contains ("/"))  {
            int n,m;
            n = byString.left("/").toInt ();
            m = byString.right("/").toInt ();
            angle = TWO_PI * double (n) / double (m);
         } else {
            angle = byString.toDouble () * PI / 180.;
         }
      } catch (SxException e)  {
         cout << "Error parsing angle option '" << byString << "':" << endl;
         e.print ();
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
   cli.last ().defaultValue = "default: same format as input file";
   if (sxbFile) outSxb = true;
   if (outSxb && outFile.getSize () == 0)
      outFile = "structure.sxb";
   
   SxArray<int> atomList = cli.option ("--atoms", "vector", 
                   "atoms to rotate").toIntList();
   cli.finalize ();

   // --- Set up rot
   if (cli.groupAvailable (axisGroup))  {
      if (axis.norm () < 1e-10)  {
         cout << "Invalid zero axis [0,0,0] (no direction)." << endl;
         SX_QUIT;
      }
      rot = SxRotation (axis, angle);
   } else if (!cli.groupAvailable (explicitGroup))  {
      cout << "No rotation given!" << endl;
      SX_QUIT;
   }

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
   }

   // shift structure if specified
   Coord shift;
   if (atomOnAxis >= structure.getNAtoms ()) {
      cout << "Rotation center is not in AtomGroup" << endl;
      SX_QUIT;
   }
   if (atomOnAxis > -1)  { 
      shift = structure(atomOnAxis);
      structure -= shift;
   }
   
   // rotate structure
   cout << "Rotation matrix:" << endl << rot << endl;
   if ((rot.transpose () - rot.inverse ()).absSqr ().sum () > 1e-5)  {
      cout << "WARNING: matrix is not unitary." << endl;
   }

   if (atomList.getSize () > 0)  {
      // --- subset rotation
      SX_LOOP(iAtom) {
         int ia = atomList(iAtom);
         structure.setAtom(ia, rot ^ structure(ia));
      }
   } else  {
      // --- full rotation
      structure ^= rot;
      structure.cell = rot ^ structure.cell;
   }
   
   if (atomOnAxis > -1)  { 
      structure += shift;
   }

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

