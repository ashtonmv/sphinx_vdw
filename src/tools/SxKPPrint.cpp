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
#include <SxKPoints.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on prints out the k-points in sx-input style "
      "or the running path";
   int wGroup = cli.newGroup ("read waves file");
   SxString wavesFile = 
      cli.option ("-w|--waves","file","waves file")
      .toString ("waves.sxb");
   int iGroup = cli.newGroup ("read input file");
   cli.excludeGroup(wGroup);
   SxString inputFile = 
      cli.option ("-i|--input","file","input file")
      .toString ("input.sx");
   cli.setGroup(cli.generalGroup);
   SxString outputFile = 
      cli.option ("-o|--output","file","output file")
      .toString ("");
   cli.last ().defaultValue = "default: stdout";

   bool kpath = cli.option ("--path","print sum of k-point differences")
                .toBool ();
   cli.finalize ();

   // --- read input

   SxCell cell;
   SxKPoints kp;

   if (cli.groupAvailable(iGroup))  {
      SxParser parser;
      SxParser::Table table = parser.read(inputFile);
      SxAtomicStructure str(&*table);
      cell = str.cell;
      kp = SxKPoints (str.cell, &*table);
   } else {
      try {
         SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
         cell.read(io);
         kp.read (io);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } 

   // --- open output file
   FILE *fp;
   if (outputFile.getSize () > 0)  {
      fp = fopen(outputFile.ascii (), "w");
      if (!fp) {
         cout << "Error: '" << outputFile << "' can't be opened for writing.\n";
         SX_EXIT;
      }
   } else {
      fp = stdout;
   }
   
   // --- print out
   SxCell recCell = cell.getReciprocalCell ();
   double dk = 0.;
   if (kpath) sxfprintf(fp,"# dk\n");
   for (int ik = 0; ik < kp.nk; ++ik)  {
      if (kpath)  {
         if (ik > 0) dk += (kp.getK(ik) - kp.getK(ik-1)).norm ();
         sxfprintf (fp, "%8.6f\n", dk);
      } else {
         Coord relK = recCell.carToRel (kp.getK(ik));
         sxfprintf (fp, "      kPoint { coords = [%.8f,%.8f,%.8f]; relative;\n"
                        "               weight = %.8f; }\n",
                        relK(0), relK(1), relK(2), 
                        kp.weights(ik));
      }
   }

   // close output
   if (fp != stdout) fclose(fp);

}

