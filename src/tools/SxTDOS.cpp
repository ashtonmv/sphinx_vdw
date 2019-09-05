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
#include <SxSpectrum.h>

int main (int argc, char **argv)
{

   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates total density of states (TDOS)\n"
                 "  Gaussian broadening applied to each state").wrap ();
   
   int wGroup = cli.newGroup ("sxb");
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   int iGroup = cli.newGroup ("sx");
   cli.excludeGroup(wGroup);
   
   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("input.sx");
   
   SxString epsFile 
      = cli.option ("-e|--eps", "eps file", 
                    "use energies in <eps file> rather than those in "
                    "waves file").toString ();
   SxString weightsFile
      = cli.option ("-f|--weights", "weights file", 
                    "use weights in <weights file> (must have eps.dat format)"
                   ).toString ("");

   bool noWeights = cli.option ("--noweights", "do not use weights")
                    .toBool ();
 
   cli.setGroup (cli.generalGroup);

   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("tdos");

   // --- remove trailing ".dat"
   if (outFile.tail (4) == ".dat") {
      outFile = outFile.head (outFile.getSize () - 4);
   }

   double broadening = cli.option ("-b|--broad", "energy [eV]", 
                                   "Gaussian weight factor, eV")
                       .toDouble(0.1) / HA2EV;

   double shift = cli.option ("-s|--shift", "energy [eV]", 
          "shift borders of energy interval from (energyMin, energyMax) "
          "to (energyMin-shift, energyMax+shift), eV")
                  .toDouble(2. * broadening * HA2EV) / HA2EV;
   cli.last ().defaultValue = "default: 2 * broadening";

   SxList<double> range = cli.option("--range","range","Emin:Emax (eV)")
                          .toDoubleList ();
   if (range.getSize () != 2 && range.getSize () != 0 && !cli.error)  {
      cout << "Illegal range." << endl;
      cli.setError ();
   }
   
   bool useOcc = cli.option("--useocc",
                            "weight states with Fermi occupation numbers\n"
                            "otherwise all states have the same weight")
                 .toBool ();

   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);
   
   cli.finalize ();

   // --- read eps and weights from waves file
   int nk = -1;
   SxVector<Double> weights;
   SxFermi fermi;
   if (cli.groupAvailable(iGroup))  {
      if (!noWeights)  {
         SxParser parser;
         SxParser::Table table = parser.read(inputFile);
         SxAtomicStructure str(&*table);
         SxKPoints kp(str.cell, &*table);
         nk = kp.getNk ();
         weights = kp.weights;
         fermi = SxFermi (1. /* nElec */, 1 /* nStates */,1 /* nSpin */,kp);
      }
   } else {
      try  {
         SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
         
         nk = io.getDimension ("nk");
         
         weights.resize (nk);
         io.read ("kWeights", &weights, nk);
         
         fermi.read (io);
         
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

   if (epsFile.getSize () > 0)  {
      int kpInFile, statesInFile;
      SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
      if (kpInFile != weights.getSize () && !noWeights)  {
         cout << "'" << epsFile << "' has " << kpInFile;
         cout << " kpoints, but " << nk << " in " << inputFile << endl;
         SX_QUIT;
      }
      if (noWeights)
         fermi = SxFermi (statesInFile, 1, nk = kpInFile);
      else
         fermi.resize (statesInFile);
      if (weightsFile.getSize () > 0)  {
         // read focc's into fermi.eps
         fermi.readSpectrum (weightsFile);
         for (int iSpin = 0; iSpin < fermi.getNSpin (); ++iSpin)
            for (int ik = 0; ik < fermi.getNk (); ++ik)
               fermi.focc(iSpin,ik) = fermi.eps(iSpin,ik) * HA2EV;
         useOcc = true;
      }
      fermi.readSpectrum (epsFile);
   }

   // --- get minimum and maximum
   int nSpin = fermi.getNSpin ();
   double eMin = fermi.eps(0,0,0),
          eMax = fermi.eps(0,0,0);
   int i, iSpin, ik;
   if (range.getSize () == 0)  {
      double eMinK,eMaxK;
      for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (ik = 0; ik < nk; ++ik)  {
            eMinK = fermi.eps(iSpin,ik).minval ();
            if (eMinK < eMin) eMin = eMinK;
            eMaxK = fermi.eps(iSpin,ik).maxval ();
            if (eMaxK > eMax) eMax = eMaxK;
         }
      }
   } else {
      eMin = range(0) / HA2EV;
      eMax = range(1) / HA2EV;
      shift = 0.;
   }

   // --- calculate spectrum

   SxSpectrum tdos;
   tdos.nPerPeak = nPerPeak;
   double weight;

   for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // --- set up new spectrum
      tdos.set (eMin - shift, eMax + shift, broadening);
      
      // --- add peaks
      for (ik = 0; ik < nk; ++ik)  {
         if (noWeights)
            weight = 1.;
         else
            weight = 2./double(nSpin) * weights(ik);
         for (i = 0; i < fermi.getNStates(ik); ++i)  {
            if (useOcc) weight = fermi.focc(i,iSpin,ik) * weights(ik);
            tdos.addPeak (fermi.eps(i,iSpin,ik), weight);
         }
      }

      // compute spectrum
      tdos.compute ();
      // change from Hartree to eV 
      tdos.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
      tdos.energies *= HA2EV;
      
      // output
      try {
         SxString name = outFile;
         if (nSpin == 2)
            name += (SxList<SxString> () << ".up" << ".down")(iSpin);
         name += ".dat"; 
         SxBinIO io(name, SxBinIO::ASCII_WRITE_ONLY);
         tdos.fprint (io.fp);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

}
