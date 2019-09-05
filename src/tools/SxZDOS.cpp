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
#include <SxPW.h>
#include <SxProjector.h>
#include <SxTimer.h>

int main (int argc, char **argv)
{

   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates z-resolved density of states\n"
                 "  Gaussian broadening applied to each state").wrap ();
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString epsFile 
      = cli.option ("-e|--eps", "eps file", 
                    "use energies in <eps file> rather than those waves file")
        .toString ("");

   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("zdos");

   SxFFT::plannerCLI (cli);
   
   double broadening = cli.option ("-b|--broad", "energy [eV]", 
                                   "Gaussian weight factor, eV")
                       .toDouble(0.1) / HA2EV;

   int shiftGroup = cli.newGroup ("automatic energy range");
   
   double shift = cli.option ("-s|--shift", "energy [eV]", 
          "shift borders of energy interval from (energyMin, energyMax) "
          "to (energyMin-shift, energyMax+shift), eV")
                  .toDouble(2. * broadening * HA2EV) / HA2EV;
   cli.last ().defaultValue = "default: 2 * broadening";

   cli.newGroup ("explicit energy range");
   cli.excludeGroup (shiftGroup);
   
   SxList<double> range = cli.option ("-r|--range", "eMin:eMax", "energy range")
                          .toDoubleList ();

   cli.setGroup(cli.generalGroup);

   SxString potFile = cli.option("--sigfile","file",
                                 "apply model self energy (z-dependent)")
                      .toString ("");

   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);

   int dim = cli.option ("--a1|--a2|--a3","which axis to resolve").toChoice ();
   if (dim < 0) dim = 2;
   cli.last ().defaultValue = "default: a3";
   
   cli.finalize ();

   // --- read eps and weights from waves file
   int nk = -1;
   SxVector<Double> kWeights;
   SxFermi fermi;
   SxPtr<SxGkBasis> gkPtr;
   SxCell cell;
   try  {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      
      nk = io.getDimension ("nk");
      
      kWeights.resize (nk);
      io.read ("kWeights", &kWeights, nk);
      fermi.read (io);
      gkPtr = SxPtr<SxGkBasis>::create (io);
      fermi.kpPtr = &*gkPtr;

      cell.read (io);

      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (epsFile.getSize () > 0)  {
      int kpInFile, statesInFile;
      SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
      if (kpInFile != kWeights.getSize ())  {
         cout << "'" << epsFile << "' has " << kpInFile;
         cout << " kpoints, but " << nk << " in " << wavesFile << endl;
         SX_QUIT;
      }
      fermi.resize (statesInFile);
      fermi.readSpectrum (epsFile);
   }

   int nSpin = fermi.getNSpin ();
   double eMin = fermi.eps(0,0,0),
          eMax = fermi.eps(0,0,0),
          eMinK,eMaxK;

   int iSpin, ik;
   if (range.getSize () == 2)  {
      eMin = range(0) / HA2EV;
      eMax = range(1) / HA2EV;
   } else {
      for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (ik = 0; ik < nk; ++ik)  {
            eMinK = fermi.eps(iSpin,ik).minval ();
            if (eMinK < eMin) eMin = eMinK;
            eMaxK = fermi.eps(iSpin,ik).maxval ();
            if (eMaxK > eMax) eMax = eMaxK;
         }
      }
   }
   
   if (nSpin == 1) kWeights *= 2.; // to be consistent with TDOS


   // --- waves and R-meshes
   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   waves.setGkBasisPtr (gkPtr);
   SxMesh3D mesh = (*gkPtr)(0).fft3d(0).mesh;
   SxRBasis R (mesh, cell);

   int meshSize = mesh.product ();
   SxVector<Int> zOnly(meshSize);
   SxVector<Int>::Iterator zOnlyIt = zOnly.begin ();
   for (int iR = 0; iR < meshSize; ++iR)
      *zOnlyIt++ = mesh.getMeshVec(iR, SxMesh3D::Positive)(dim);

   SxMeshR rho;
   SxMeshR::Iterator rhoIt;
   int nPts = mesh(dim);
   SxVector<Double> weight(nPts);

   bool selfCorr = (potFile.getSize () > 0);
   SxVector<Double> sigma;
   if (selfCorr)  {
      FILE *fp = fopen(potFile.ascii (), "r");
      if (!fp)  {
         cout << "Invalid sigma file " << potFile << "." << endl;
         SX_QUIT;
      }
      sigma.resize (nPts);
      for (int i = 0; i <  nPts; ++i)  {
         fscanf(fp,"%lf",&sigma(i));
         cout << i << " " << sigma(i) << endl; 
      }
   }
      

   // --- calculate spectrum
   SxSpectrum zdos;
   zdos.nPerPeak = nPerPeak;
   int i;
   for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // --- set up new spectrum
      zdos.set (eMin - shift, eMax + shift, broadening, nPts);
      
      // --- add peaks
      for (ik = 0; ik < nk; ++ik)  {
         for (i = 0; i < fermi.getNStates(ik); ++i)  {
            if (fermi.eps(i,iSpin,ik) >= eMin && fermi.eps(i,iSpin,ik) <= eMax)
            {
               double eps = fermi.eps(i,iSpin,ik);
               rho = (R | waves(i,iSpin,ik)).absSqr ();
               rhoIt = rho.begin ();
               zOnlyIt = zOnly.begin ();
               weight.set (0.);
               for (int iR = 0; iR < meshSize; ++iR)
                  weight(*zOnlyIt++) += *rhoIt++;
               if (selfCorr)  {
                  eps += (weight * sigma).sum () * R.dOmega;
                  cout << eps << endl;
               }
               weight *= kWeights(ik) * R.dOmega;
               zdos.addPeak (eps, weight);
            }
         }
      }

      // compute spectrum
      zdos.compute ();
      // change from Hartree to eV 
      zdos.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
      zdos.energies *= HA2EV;
      
      // output
      try {
         SxString name = outFile;
         if (nSpin == 2)
            name += (SxList<SxString> () << ".up" << ".down")(iSpin);
         name += ".dat"; 
         SxBinIO io(name, SxBinIO::ASCII_WRITE_ONLY);
         zdos.fprint (io.fp);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

}
