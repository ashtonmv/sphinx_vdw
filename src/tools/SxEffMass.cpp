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
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxPseudoPot.h>

/// Wrapper class to access SxFermi internal with simplified interface
class SxFermiWrapper : public SxFermi
{
   public:
      double fermiFunction (double energy, double ekt)
      {
         beta = 1./ekt;
         return SxFermi::fermiFunction (energy, nElectrons, /*spin=*/0);
      }
      using SxFermi::operator=;
};

int main (int argc, char **argv)
{

   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "B.Lange";
   cli.preUsageMessage = 
      "This add-on calculates the effective Mass of a free electron\n";

   int wGroup = cli.newGroup ("sxb");
   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");
   int iGroup = cli.newGroup ("sx");
   cli.excludeGroup(wGroup);
   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("input.sx");
   
   SxString epsFile 
      = cli.option ("-e|--eps", "eps file", 
                    "use energies in <eps file> rather than those in "
                    "waves file").toString ();
   int bandGroup = cli.newGroup ("band");
   double eMin = cli.option("--VBM","Valence Band Maximum in eV").toDouble();
   double eMax = cli.option("--CBM","Conduction Band Minimum in eV").toDouble();
   cli.setGroup (cli.generalGroup);
   double ekt = cli.option("--ekt","electron temperature (eV)").toDouble(0.025);
   int size = cli.option("--size","number of data points").toInt(1000,10,100000);
   int test = cli.option("--test", "testmode #kPoints per dir").toInt(0,0,100000);
   int holeband = cli.option("--hband", "bandidx of holeband").toInt(-1);
   int elecband = cli.option("--eband", "bandidx of electronband").toInt(-1);
   
   cli.finalize ();
   
   ekt /= 27.2114;
   eMin /= 27.2114;
   eMax /= 27.2114;
   SxFermiWrapper fermi;
   SxKPoints kPoints;
   SxCell myCell;
   int nSpin = 1, nk = -1;
   double norm=1;
   if(!test)   {
      if (cli.groupAvailable(iGroup))  {
         SxParser parser;
         SxParser::Table table = parser.read(inputFile);
         myCell = SxCell(&*table);
         SxAtomicStructure str(&*table);
         SxSpeciesData mySpeciesData (&*table);
         double nElec = 0;
         for (int is = 0; is < str.getNSpecies() ; is++)  {
            double charge = mySpeciesData.valenceCharge(is);
            nElec += str.getNAtoms(is)*charge;
         }
         kPoints = SxKPoints(str.cell, &*table);
         nk = kPoints.getNk ();
         fermi = SxFermi (int(nElec), 1 /* nStates */, nSpin, kPoints);
         int kpInFile, statesInFile;
         SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
         fermi.resize (statesInFile);
         fermi.readSpectrum (epsFile);
         // Fermi Distribution of electrons
         // First use a higher temperature to put Fermienergy midgap
         fermi.fermiDistribution (1e-4);
         // Second use low temperure to avoid filling higher levels in semiconductors
         fermi.fermiFunction (fermi.eFermi, 1e-10);
         fermi.updateValCon();
      } else {
         try {
            SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
            fermi.read (io);
            kPoints.read(io);
            fermi.kpPtr = &kPoints;
            myCell.read(io);
            nSpin  = io.getDimension ("nSpin");
            nk = io.getDimension("nk");
            io.close ();
         }
         catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }
      cout << "#nElectrons = "<< fermi.nElectrons << endl;
   
      if(!(fermi.isSemiconductor ()))   {
         cout << "This seems not to be a semiconductor!" << endl;
         cout << "Tool works only for materials with a band gap." << endl;
         for(int ik = 0; ik < nk; ik++)   {
            cout << ik << "\t" << fermi.nValBands(ik)(0) << endl;
         }
         SX_QUIT;
      }
      cout << "Volume: " << myCell.volume << endl;
      norm = 1.0 / pow(myCell.volume,2.0/3.0);
      int NVB = fermi.getNValenceBands();
      cout << "#nValenceBands = "<< NVB << endl;

      if (!cli.groupAvailable(bandGroup)) {
         eMin = fermi.eps(NVB-2,0,0);
         eMax = fermi.eps(NVB+1,0,0);
         for(int ik = 0; ik < nk; ik++)   {
            for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
               int NValenceBands = fermi.getNValenceBands(iSpin, ik);
               if (eMin < fermi.eps(NValenceBands-1,iSpin,ik)) 
                  eMin = fermi.eps(NValenceBands-1,iSpin,ik); 
               if (eMax > fermi.eps(NValenceBands,iSpin,ik))
                  eMax = fermi.eps(NValenceBands,iSpin,ik);
            }
         }
      }

      fermi.fermiDistribution (ekt);
     
      int nState = fermi.getNStates ();
      /*
      cout << "Original" << endl; 
      for(int iState = 0; iState < nState; iState++)   {
         cout << iState << "\t" << fermi.eps(iState,0,0) << endl;
      }
      */

      if ((holeband >=0) || (elecband >=0))   {
         for(int iState = 0; iState < nState; iState++)   {
            if ((iState != holeband) && (iState != elecband))   {
               for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
                  for(int ik = 0; ik < nk; ik++)   {
                     if (fermi.eps(iState,iSpin,ik) < fermi.eFermi)   {
                        fermi.eps(iState,iSpin,ik) -= 0.1;
                     }
                     else   {
                        fermi.eps(iState,iSpin,ik) += 0.1;
                     }
                  }
               }
            }
         }
         /*
         cout << "Shifted" << endl;
         for(int iState = 0; iState < nState; iState++)   {
            cout << iState << "\t" << fermi.eps(iState,0,0) << endl;
         }
         */
      }
   } else   {
      nSpin = 1;
      kPoints.nk = nk = test;
      double volume = 6*PI*PI*(1.0 + 1.5/double(nk) + 0.5/double(nk*nk));
      cout << "Volume: " << volume << endl;
      norm = 1.0 / pow(volume,2.0/3.0);
      kPoints.kVec = SxArray<SxVector3<Double> >(nk);
      kPoints.weights.resize(nk); 
      fermi = SxFermi (2.0,2,1, kPoints);
      SxVector<Double> VB(nk);
      SxVector<Double> CB(nk);
      for(int ik = 0; ik < nk; ik++)   {
         kPoints.kVec(ik) = SxVector3<Double>(double(ik)/kPoints.nk,0,0);
         kPoints.weights(ik) = (ik+1.0) * (ik+1.0);
         VB(ik) = fermi.eps(0,0,ik) = -0.1*kPoints.kVec(ik).normSqr ();
         fermi.focc(0,0,ik) = 2.0;
         CB(ik) = fermi.eps(1,0,ik) = kPoints.kVec(ik).normSqr () + 0.1;
         fermi.focc(0,0,ik) = 0.0;
         fermi.nValBands(ik)(0)=0;
         fermi.nConBands(ik)(0)=1;
         
      }
      try{
         SxString file = "VB.dat";
         SxBinIO io (file, SxBinIO::ASCII_WRITE_ONLY);
         io.writeXYPlot(VB);
         io.close ();
      }
      catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }
      try{
         SxString file = "CB.dat";
         SxBinIO io (file, SxBinIO::ASCII_WRITE_ONLY);
         io.writeXYPlot(CB);
         io.close ();
      }
      catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }


      kPoints.weights /= kPoints.weights.sum ();
      fermi.kpPtr = &kPoints;
      fermi.nElectrons = 2.0;
      cout << "#nElectrons = "<< fermi.nElectrons << endl;
      
      if(!cli.groupAvailable(bandGroup)) {
         eMin = -0.05;
         eMax = 0.15;
         for(int ik = 0; ik < nk; ik++)   {
            for(int iSpin = 0; iSpin < nSpin; iSpin++)   {
               if (eMin < fermi.eps(0,iSpin,ik)) 
                  eMin = fermi.eps(0,iSpin,ik); 
               if (eMax > fermi.eps(1,iSpin,ik))
                  eMax = fermi.eps(1,iSpin,ik);
            }
         }
      }
     

      fermi.fermiDistribution (ekt);
   } 
   

   cout << "#eMin:   " << eMin*27.2114 << endl;
   cout << "#eMax:   " << eMax*27.2114 << endl;
   cout << "#eFermi: " << fermi.eFermi*27.2114 << endl;
   cout << "#ekt: " << ekt*27.2114 << endl;
   

   double stepwidth = 1.5*(eMax - eMin) / size;   
   SxVector<Double> energies (size+1);
   energies.set(0.0);
   SxVector<Double> ElecDiff (size+1);
   ElecDiff.set(0.0);
   // Conduction band
   double Deff = 0.0;
   for(int iStep=0; iStep<=size; iStep++)   {
      double energy = eMax - iStep * stepwidth;
      // Set CBM to zero
      energies(iStep) = energy - eMax;
      if(iStep == size)   {
         cout << "No sufficient exponential behavior found!"<< endl;
      }
     ElecDiff(iStep) += fermi.fermiFunction(energy, ekt);
     if ( ((ElecDiff(iStep)) < -1e-10) && (ElecDiff(iStep) > -1e-2))   {
        double DeffCheck = log(fabs(ElecDiff(iStep))) - energies(iStep)/ekt;
        if (fabs(Deff - DeffCheck) > 1e-4)   {
           Deff = DeffCheck;
        } else   {
           cout << "Exponential behavior at energy [eV]: " << energy*27.2114 << endl;
           break; // value found;
        }
      }
   }

   //  Get axisintersection
   double y0 = Deff - 1.5 * log(ekt);
   // Calculate Mass (scriptum)
   double m = exp(2.0/3.0 * y0 + 1.0/3.0 * log(2.0) + log(PI)) * norm;
   cout << "m_eff(Electrons): " << m << endl;

  try{
      SxString file = "nElecVsE-CB-";
      file += ekt*27.2114;
      file += "eV.dat";
      SxBinIO io (file, SxBinIO::ASCII_WRITE_ONLY);
      io.writeXYPlot(energies, -ElecDiff);
      io.close ();
   }
   catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   
   
   //Valence Band
   Deff = 0.0;
   energies.set(0.0);
   ElecDiff.set(0.0);
   for(int iStep=0; iStep<=size; iStep++)   {
      if(iStep == size)   {
         cout << "No sufficient exponential behavior found!"<< endl;
      }
      double energy = eMin + iStep * stepwidth;
      // set VBM to zero
      energies(iStep) = energy - eMin;
      ElecDiff(iStep) += fermi.fermiFunction(energy, ekt);
      if ( ((ElecDiff(iStep)) > 1e-10) && (ElecDiff(iStep) < 1e-2))   {
         double DeffCheck = log(fabs(ElecDiff(iStep))) + energies(iStep)/ekt;
         if (fabs(Deff - DeffCheck) > 1e-4)   {
           Deff = DeffCheck;
        } else   {
           cout << "Exponential behavior at energy [eV]: " << energy*27.2114 << endl;
           break; // value found;
        }
      }
   }
   //  Get axisintersection
   y0 = Deff - 1.5 * log(ekt);
   // Calculate Mass (scriptum)
   m = exp(2.0/3.0 * y0 + 1.0/3.0 * log(2.0) + log(PI)) * norm;
   cout << "m_eff(Holes): " << m << endl;

   try{
      SxString file = "nElecVsE-VB-";
      file += ekt*27.2114;
      file += "eV.dat";
      SxBinIO io (file, SxBinIO::ASCII_WRITE_ONLY);
      io.writeXYPlot(energies, ElecDiff);
      io.close ();
   }
   catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
    
}

