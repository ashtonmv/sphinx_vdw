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

#include <SxFermi.h>
#include <SxCLI.h>
#include <SxBinIO.h>

double getEffFermiDos(const SxFermi &fermi)
{
   SX_CHECK (fermi.kpPtr);
   double dos=0;
   SX_LOOP3(ik,iSpin(fermi.getNSpin ()),i(fermi.getNStates ()))
      dos += fermi.dFoccFermi(i,iSpin,ik) * fermi.kpPtr->weights(ik);
   return dos / fermi.beta;
}

void effdosplot (SxFermi &fermi)
{
   ofstream out("effDos.dat");
   for (double ekt = 0.005; ekt < 0.5; ekt += 0.005)  {
      fermi.fermiDistribution (ekt / HA2EV);
      out << ekt << ' ' << getEffFermiDos (fermi) / ekt
           << ' ' << fermi.eFermi * HA2EV << endl;
   }
}

void freeplot (SxFermi &fermi, double ekt0)
{
   fermi.fermiDistribution (ekt0);
   double Uref = 0.;
   SX_LOOP3(ik,iSpin,i)  {
      double kw = fermi.kpPtr->weights(ik);
      Uref += fermi.focc(i,iSpin,ik) * fermi.eps(i, iSpin, ik) * kw;
      double fn = fermi.focc(i,iSpin,ik) * fermi.getNSpin () / 2.;
   }
   double Sref = fermi.getEntropy ();
   //double Fref = Uref - ekt0 * Sref;

   ofstream out("free.dat");
   out << "# ekt U S F" << endl;
   for (double ekt = 0.01; ekt < 0.6; ekt += 0.01)  {
      fermi.fermiDistribution (ekt / HA2EV);
      double U = 0., S = 0.;
      SX_LOOP3(ik,iSpin,i)  {
         double kw = fermi.kpPtr->weights(ik);
         U += fermi.focc(i,iSpin,ik) * fermi.eps(i, iSpin, ik) * kw;
      }
      S = fermi.getEntropy ();
      double F = U - ekt/HA2EV * S;

      out << ekt << ' ' << (U-Uref) * HA2EV << ' ' << S << ' '
          << (F-Uref) * HA2EV << endl;
   }
}


int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates the Fermi energy for a given temperature.";
   SxString wavesFile = 
      cli.option ("-w|--waves","file","waves file")
      .toString ("waves.sxb");
   SxString epsFile =
      cli.option ("--eps","file","override energies from eps.dat file")
      .toString ("");

   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("");

   double nEl = cli.option ("-n|--nEl", "electrons", 
                            "number of Electrons")
                .toDouble (0.0,0);

   int orderMethfesselPaxton
      = cli.option ("-N", "int", "Methfessel-Paxton order").toInt (-1);
   cli.last ().defaultValue = "default: Fermi-Dirac smearing";
   int normal = cli.newGroup ("compute Fermi energy/occupations");
   double ekt = cli.option ("-e|--ekt", "energy", "temperature in eV")
                .toDouble () / HA2EV;
   bool print = cli.option ("-p|--print", "print occupation numbers").toBool ();

   cli.newGroup ("effektive DOS @ Fermi-Level kT-profile");
   cli.excludeGroup (normal);
   
   bool effDos = cli.option ("--dosprofile","compute eff. DOS (dN/dkT)")
                 .toBool ();
   cli.newGroup ("electronic free energy");
   bool free = cli.option ("--free", "compute free energy")
               .toBool ();

   cli.finalize ();

   // --- read input

   SxFermi fermi;
   SxCell cell;
   SxKPoints kp;

   try {
      if (epsFile.getSize () == 0 || inputFile.getSize () == 0)  {
         SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
         fermi.read (io);
         cell.read(io);
         kp.read (io);
      } else {
         SxParser parser;
         SxParser::Table table = parser.read(inputFile);
         SxAtomicStructure str(&*table);
         cell = str.cell;
         kp.read (cell.getReciprocalCell (), table->getGroup ("basis"),
                  "k", "kPoint", "kPoints");
         if (fabs(nEl) < 1e-16)  {
            cout << "You must specify the number of electrons" << endl;
            SX_QUIT;
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (epsFile.getSize () > 0)  {
      int kpInFile, statesInFile;
      SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
      if (kpInFile != kp.nk)  {
         cout << "'" << epsFile << "' has " << kpInFile;
         cout << " kpoints, but " << kp.nk << " in "
              << ((inputFile.getSize () > 0) ? inputFile : wavesFile) << endl;
         SX_QUIT;
      }
      fermi = SxFermi (statesInFile, 1, kpInFile);
      fermi.readSpectrum (epsFile, &cell, &kp);
   }
   fermi.orderMethfesselPaxton = orderMethfesselPaxton;

   // --- Recalculate Fermi distribution
   // get current number of electrons
   if (fabs(nEl) < 1e-6)  {
      nEl = 0.;
      for (int ik = 0; ik < fermi.getNk (); ik++)
         nEl += fermi.focc(0/* iSpin */,ik).sum () * kp.weights(ik);
   }
   cout << "nElectrons = " << nEl << endl;
   fermi.nElectrons = nEl;
   // set kpPtr
   fermi.kpPtr = &kp;
   if (effDos)  {
      effdosplot (fermi);
      return 0;
   }
   if (free)  {
      freeplot (fermi, ekt);
      return 0;
   }
   // recalculate
   fermi.fermiDistribution (ekt);
   if (print) fermi.printOccupation ();
   cout << "Fermi energy = " << fermi.eFermi * HA2EV << " eV." << endl;

   cout << "Material ";
   if (fermi.isSemiconductor ()) cout << "has a bandgap of " << fermi.getBandGap ()*HA2EV << " eV" << endl;
   else cout << "is a metal." << endl;
}

