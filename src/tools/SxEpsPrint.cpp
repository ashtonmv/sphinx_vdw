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
#include <SxFermi.h>

int main (int argc, char **argv)
{
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage = "Print out eps.dat from waves file";
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString outFile
            = cli.option ("-o","file","output file").toString ("waves-eps.dat");

   bool printBandEnergy
      = cli.option ("--bandEnergy", " print band energy").toBool ();
   
   cli.finalize ();

   // --- read eps and weights from waves file
   SxFermi fermi;
   SxKPoints kp;
   try  {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      fermi.read (io);
      kp.read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   fermi.kpPtr = &kp;

   fermi.writeSpectrum (outFile, "");
   if (printBandEnergy)
      sxprintf ("band energy: %.12f", fermi.getEBand(SxFermi::UseFocc) * HA2EV);

   return 0;

}
