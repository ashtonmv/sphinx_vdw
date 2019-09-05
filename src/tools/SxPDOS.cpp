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
#include <SxPartialWaveBasis.h>
#include <SxAtomicStructure.h>
#include <SxPAWSet.h>
#include <SxProjector.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates projected density of states for PAW\n"
                 "  Gaussian broadening applied to each state"
                 " The first column of output is the total PDOS, and"
                 " proceeding data columns are PDOS of respective l.").wrap ();
   cli.authors = "Siyuan Zhang and Christoph Freysoldt";
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("input.sx");
   
   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("pdos");
   
   SxList<int> atoms = cli.option ("-a|--atoms", "atom list", 
                                   "atoms of interest for PDOS")
                      .toIdxList ();
   cout << "Entered atoms contain number: " << atoms << endl;
   cli.last ().required (false);
   cli.last ().defaultValue 
            = "default returns a total spectrum but no individuals";

   int wop = cli.option ("--wop", "1: write complete outputs; default: 0")
            .toInt(0);

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
   
   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);
   
   cli.finalize ();

   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read(inputFile);

   SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);
  
   // --- read waves file (not waves yet)
   int nk = -1;
   SxFermi fermi;
   SxVector<Double> weights;
   SxAtomicStructure structure;
   SxPtr<SxGkBasis> gkBasisPtr;
   SxPtr<SxPartialWaveBasis> pBasisPtr;
   SxBinIO io;
   try  {
      io.open (wavesFile, SxBinIO::BINARY_READ_ONLY);
      nk = io.getDimension ("nk");
      weights.resize (nk);
      io.read ("kWeights", &weights, nk);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   fermi.read (io);

   structure.read (io);
   gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
   gkBasisPtr->changeTau (structure);

   // create partial wave basis
   pBasisPtr = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
   pBasisPtr->createProjBasis (*gkBasisPtr);
   int nSpin = fermi.getNSpin ();
   // construct PAW waves container
   SxPAWSet waves (gkBasisPtr, pBasisPtr, fermi.getNStates (), nSpin);
   
   waves.read (io, SxPWSet::KeepGkBasis);
   io.close ();

   const SxPartialWaveBasis &pBasis = *pBasisPtr;
   const SxPAWPot &pawPot = *pawPotPtr;   
   SxDiracVec<Complex16> projpsi = pBasis | waves (0,0,0);

   // --- get minimum and maximum
   double eMin = fermi.eps(0,0,0),
          eMax = fermi.eps(0,0,0);
   int in, iSpin, ik;
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
   //int nPdosAtom = atoms.getSize ();
   int lMax = pawPot.getLMax ();

   SxSpectrum pdos;
   pdos.nPerPeak = nPerPeak;
   SxComplex16 kWeight, pWeight;
   SxVector<Double> specWeight (lMax + 2);

   for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
      pdos.set (eMin - shift, eMax + shift, broadening, lMax + 2);
      int ipBasis = 0;

      // --- check atom number
      int nAtom = structure.getNAtoms ();
      for (int iAtom = 0; iAtom < nAtom; ++iAtom)  {
         int iSpecies = structure.getISpecies (iAtom);
         int npt = pawPot.getNProjType (iSpecies);
         bool lastAtom = (iAtom == nAtom - 1);

         int aWrite = 0;
         double aWeight = 0.;
         if (atoms.getSize () == 0) aWeight = 1.;
         else if (atoms.contains (iAtom))  {
            aWrite = 1;
            aWeight = 1.;
         }

         SxSpectrum pdosA;
         pdosA.nPerPeak = nPerPeak;
         if (aWrite) 
            pdosA.set (eMin - shift, eMax + shift, broadening, lMax + 2);

         // ---integration
         SxDiracVec<Double> rad = pawPot.rad(iSpecies);
         double logdr = pawPot.logDr (iSpecies);
         double rPAW = pawPot.rc (iSpecies) * sqrt(-log(1e-4));
         cout << "Species " << iSpecies << " has rPAW=" << rPAW << endl;
         SxDiracVec<Double> rW (rad.getSize ());
         for (int ir = 0; ir < rad.getSize (); ++ir)  {
            if (rad(ir) < rPAW) rW(ir) = rad(ir) * rad(ir) * rad(ir);
            else rW(ir) = 0.;
         }

         for (int ipt = 0; ipt < npt; ++ipt)  {
            int li = pawPot.lPhi(iSpecies)(ipt);
            cout << "Atom " << iAtom << " type " << ipt << " has l=" 
                 <<  li << endl;
            int jpBasis = ipBasis;
            for (int jpt = ipt; jpt < npt; ++jpt)  {
               int lj = pawPot.lPhi(iSpecies)(jpt);
               if (lj == li)  {
                  double integral = (rW * pawPot.phiAE(iSpecies).colRef(jpt)
                  * pawPot.phiAE(iSpecies).colRef(ipt)).integrate (logdr);
                  cout << "ipt " << ipt << " jpt " << jpt 
                       << " integral=" << integral << endl;

                  // --- add peaks
                  for (ik = 0; ik < nk; ++ik)  {
                     kWeight = 2./double(nSpin) * weights(ik) * aWeight;
                     for (in = 0; in < fermi.getNStates(ik); ++in)  {

                        // --- calculate projected weight
                        projpsi = pBasis | waves (in, iSpin, ik);
                        pWeight = 0.;
                        specWeight.set (0.);
                        // --- sum <psi|p_j><p_i|psi> over m (im = m + l)
                        SxComplex16 sum = 0.;
                        for (int im = 0; im < 2 * li + 1; ++im)  {
                           SxComplex16 element = projpsi(jpBasis + im).conj () 
                                               * projpsi(ipBasis + im);
                           sum += element;
                        }
                        if (jpt == ipt) pWeight += integral * sum;
                        else pWeight += integral * (sum + sum.conj());
                        specWeight(0) = kWeight * pWeight;
                        specWeight(li + 1) = kWeight * pWeight;
                        pdos.addPeak (fermi.eps(in,iSpin,ik), specWeight);
                        if (aWrite)
                           pdosA.addPeak (fermi.eps(in,iSpin,ik), specWeight);
                     } //in
                  } //ik
                  jpBasis += 2 * lj + 1;
                  if (wop) cout << "jpBasis=" << jpBasis << endl;
               } //check jpt
            } //jpt
            ipBasis += 2 * li + 1;
            if (wop) cout << "ipBasis=" << ipBasis << endl;
         } //ipt
      
         // compute spectrum
         if (aWrite)  {
            pdosA.compute ();
            pdosA.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
            pdosA.energies *= HA2EV;
         }

         if (lastAtom)  {
            pdos.compute ();
            pdos.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
            pdos.energies *= HA2EV;
         }

         // output
         try {
            SxString ext = ".dat";
            if (nSpin == 2)  {
               if (iSpin == 0) ext = ".up";
               else ext = ".down";
            }

            if (aWrite)  {
               SxString name = outFile + "Atom" + SxString (iAtom) + ext;
               io.open (name, SxBinIO::ASCII_WRITE_ONLY);
               pdosA.fprint (io.fp);
               io.close ();
            }

            if (lastAtom)  {
               SxString name = outFile + "sum" + ext;
               io.open (name, SxBinIO::ASCII_WRITE_ONLY);
               pdos.fprint (io.fp);
               io.close ();
            }
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }

      } //atom
   } //spin
}
