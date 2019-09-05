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
#include <SxBinIO.h>
#include <SxPtr.h>
#include <SxAtomicStructure.h>
#include <SxKPoints.h>
#include <SxParser.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   
   // --- parse command line
   SxCLI cli(argc, argv);
   cli.preUsageMessage = "Find the symmetrized k-point coordinates.";
   cli.authors = "C.Freysoldt";
   
   int wavesGroup = cli.newGroup ();
   SxString wavesFile 
      = cli.option ("-w|--waves", "file", "wave function sxb file")
        .toString("waves.sxb");
   
   cli.newGroup ();
   cli.excludeGroup (wavesGroup);
   
   SxString inputFile 
      = cli.option ("-i|--input", "file", "S/PHI/nX input file")
        .toString("input.sx");

   cli.setGroup (cli.generalGroup);
   float delta
      = cli.option ("-d|--delta", "float", "equality criterium")
        .toFloat (1e-6, 0.);
   bool xyOnly = cli.option ("--xy", "ignore z component").toBool ();

   bool allK = cli.option ("--all","list all equivalent k-points").toBool ();
   
   cli.finalize ();

   SxAtomicStructure structure;
   SxPtr<SxKPoints> kPoints;

   if (cli.groupAvailable (wavesGroup))  {
      try  {
         SxBinIO io(wavesFile, SxBinIO::BINARY_READ_ONLY);
         kPoints = SxPtr<SxKPoints>::create ();
         kPoints->read (io);
         structure.read (io);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else {
      SxParser parser;
      SxParser::Table table = parser.read (inputFile);
      structure = SxAtomicStructure (&*table);
      kPoints = SxPtr<SxKPoints>::create (structure.cell, &*table);
   }
   
   cout << "Enter k-coordinates:" << endl;
   SxVector3<TPrecG> kPoint;
   for (int dim = 0; dim < 3; dim++)
      cin >> kPoint(dim);


   // --- find kPoint
   SxArray<SymMat> symOps = structure.cell.symGroupPtr->getSymmorphic ();
   int nSym = int(symOps.getSize ());

   int ik, iSym;
   SxVector3<TPrecG> diffVec;
   bool found = false;
   SxCell recCell = structure.cell.getReciprocalCell ();

   for (ik = 0; ik < kPoints->nk && (!found || allK); ik++)  {
      for (iSym = 0; iSym < 2 * nSym; iSym++)  {
         if (iSym < nSym) 
            // test for S^kp = k
            diffVec = (symOps(iSym) ^ kPoint) - kPoints->getK(ik);
         else 
            // test for S^kp = -k 
            diffVec = (symOps(iSym-nSym) ^ kPoint) + kPoints->getK(ik);
         // map back to 1st cell
         recCell.map (&diffVec, SxCell::Origin);
         if (xyOnly) diffVec(2) = 0.;
         if (diffVec.normSqr () < delta)  {
            cout << "ik = " << (ik + 1) << ": " << kPoints->getK(ik) << endl;
            found = true;
            break;
         }
      }
   }
   if (!found)
      cout << endl << "Not found." << endl;
   
   return (found ? 0 : 2);
}

