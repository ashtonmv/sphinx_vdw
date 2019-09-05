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
#include <SxFourierInterpol.h>
#include <SxFileParser.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This interpolates a line (or any set of points) through a mesh.";

   int lGroup = cli.newGroup ("line");
   // --- line definition
   SxVector3<Double> start (cli.option ("--from","vector", "starting point")
                            .toList3 ());

   SxVector3<Double> end (cli.option ("--to","vector","end point").toList3 ());
   cout << "end = " << end << endl;

   int nPoints = cli.option ("-n|--nPoints", "number",
                             "number of points on line").toInt (100,2);
   int iGroup = cli.newGroup ("points from file");
   cli.excludeGroup (lGroup);
   SxString pointsFile = cli.option ("-p|--points", "file", "files with x y z")
                         .toString ("");
   cli.setGroup (cli.generalGroup);

   bool relative = cli.option ("--relative", "assume relative coordinates")
                   .toBool ();

   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   SxFFT::plannerCLI (cli);

   SxString meshFile =
      cli.option ("-i","mesh file","netcdf mesh file (one mesh only)")
      .toString ();

   cli.finalize ();

   // --- read input

   SxBinIO io;
   try {
      io.open (meshFile, SxBinIO::BINARY_READ_ONLY);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   SxFourierInterpol interpol(io);

   interpol.condense ();

   io.close ();

   SxAtomicStructure rVec;
   Coord dr;
   if (cli.groupAvailable (iGroup) && pointsFile.getSize () > 0)  {
      SxFileParser fp(pointsFile);
      while (!feof(fp.fp))  {
         SxVector3<Double> xyz;
         fp >> xyz(0) >> xyz(1) >> xyz(2);
         fp.skipWhite ();
         rVec.addAtom (xyz);
      }
      rVec.endCreation ();
      if (relative) rVec ^= interpol.cell;
   } else {
      cout << "start = " << start << endl;
      if (relative)  {
         interpol.cell.changeToCar(&start);
         interpol.cell.changeToCar(&end);
         cout << "start(abs) = " << start << endl;
         cout << "end(abs)   = " << end << endl;
      }

      rVec.resize (nPoints);

      dr = (end - start) / double (nPoints - 1);

      for (int i = 0; i < nPoints; i++)
         rVec.setAtom (i, start + i * dr);
   }

   SxVector<Double> result = interpol.interpolate (rVec);
//   result /= interpol.meshInG.getSize ();

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

   if (cli.groupAvailable (iGroup))  {
      SX_LOOP(i)  {
         const SxVector3<Double> &r = rVec.constRef (i);
         sxfprintf (file, "%.10f\t%.10f\t%.10f\t%.10f\n", r(0), r(1), r(2),
                    result(i));
      }
   } else {
      double length = sqrt (dr.absSqr ().sum ());
      for (int i = 0; i < nPoints; i++)
         sxfprintf (file, "%.10f %.10f\n", i * length, result(i));
   }

   if (file != stdout) fclose (file);

}

