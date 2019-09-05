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
#include <SxRho.h>
#include <SxPseudoPot.h>

int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates z-resolved density of states\n"
                 "  Gaussian broadening applied to each state").wrap ();
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString rhoFile
            = cli.option ("-r|--rho","file","input rho file")
             .toString   ("rho.sxb");
   
   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("mies");
   
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

   double rhoValue = cli.option("-d|--density","double","target density")
                     .toDouble(1e-3,0.);

   bool bottom = cli.option("--bottom","compute MIES for bottom surface")
                 .toBool ();

   SxList<double> range = cli.option ("-e|--range", "eMin:eMax", "energy range")
                          .toDoubleList ();

   cli.setGroup(cli.generalGroup);

   SxFFT::plannerCLI (cli);

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
   SxGkBasis gk;
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
   SxFFT3d &fft = gk(0).fft3d(0);
   SxRBasis R (fft.mesh, cell);

   int meshSize = fft.meshSize;

   int xySize = fft.meshSize / fft.mesh(dim);
   int zSize = fft.mesh(dim);

   SxVector<Double> weightXY(xySize);
   SxVector<Int> whichMeshIndex(xySize);
   try  {
      // load total rho
      SxBinIO io(rhoFile, SxBinIO::BINARY_READ_ONLY);
      SxRho rhoR(io);
      SxMeshR rhoTotal(meshSize);
      rhoTotal.set (0.);
      for (iSpin = 0; iSpin < nSpin; ++iSpin)
         rhoTotal += rhoR(iSpin);
      SxMeshR rhoZ(zSize);
      rhoZ.set (0.);
      // average rho over xy 
      for (int i = 0; i < meshSize; ++i)
         rhoZ(fft.mesh.getMeshVec(i,SxMesh3D::Positive)(dim)) += rhoTotal(i);
      double min = rhoZ(0);
      int start = 0;
      for (int i = 0; i < zSize; ++i)
         if (min > rhoZ(i)) min = rhoZ(start=i);
      // cout << "start=" << start << endl;

      SxVector3<Int> indexMap(0,1,dim);
      indexMap(dim) = 2;
      // set up map xy -> (z(DOS=rhoValue), DOS(x,y,z))
      ssize_t idxXY = 0, rhoIdx;
      SxVector3<Int> xyz;
      int &x = xyz(indexMap(0)), &y = xyz(indexMap(1)), &z = xyz(indexMap(2));
      for (x = 0; x < fft.mesh(indexMap(0)); x++)  {
         for (y = 0; y < fft.mesh(indexMap(1)); y++, idxXY++)  {
            for (z = start + (bottom ? 1 : -1); z != start; 
                 z+= (bottom ? 1 : -1))  {
               rhoIdx = fft.mesh.getMeshIdx(xyz, SxMesh3D::Unknown);
               if (rhoTotal(rhoIdx) > rhoValue)  {
                  weightXY(idxXY) = 1./rhoTotal(rhoIdx);
                  whichMeshIndex(idxXY) = (int)rhoIdx;
                  /*
                  cout << "x=" << x << "; y = " << y << "; z:=" << z
                       << ", rho=" << rhoTotal(rhoIdx) << endl;
                  */
                  break;
               }
            }
         }
      }
      SX_CHECK (idxXY == xySize);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   SxMeshR rho;

   // --- calculate spectrum
   SxSpectrum mies;
   mies.nPerPeak = nPerPeak;
   for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
      // --- set up new spectrum
      mies.set (eMin - shift, eMax + shift, broadening);
      
      // --- add peaks
      for (ik = 0; ik < nk; ++ik)  {
         for (int i = 0; i < fermi.getNStates(ik); ++i)  {
            if (fermi.eps(i,iSpin,ik) >= eMin && fermi.eps(i,iSpin,ik) <= eMax
                && fermi.focc(i,iSpin,ik) > 1e-5)
            {
               rho = (R | waves(i,iSpin,ik)).absSqr ();
               double weight = 0.;
               SxVector<Int>::Iterator idxIt = whichMeshIndex.begin ();
               SxVector<Double>::Iterator wIt = weightXY.begin ();
               for (int iXY = 0; iXY < xySize; ++iXY, ++idxIt, ++wIt)
                  // weight += rho(whichMeshIndex(iXY)) * weightXY(iXY);
                  weight += rho(*idxIt) * (*wIt);
               weight *= kWeights(ik) * R.dOmega * fermi.focc(i, iSpin, ik);
               mies.addPeak (fermi.eps(i,iSpin,ik), weight);
            }
         }
      }

      // compute spectrum
      mies.compute ();
      // change from Hartree to eV 
      mies.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
      mies.energies *= HA2EV;
      
      // output
      try {
         SxString name = outFile;
         if (nSpin == 2)
            name += (SxList<SxString> () << ".up" << ".down")(iSpin);
         name += ".dat"; 
         SxBinIO io(name, SxBinIO::ASCII_WRITE_ONLY);
         mies.fprint (io.fp);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }
   return 0;

}

