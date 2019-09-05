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

#include <SxRBasis.h>
#include <SxPerturbK.h>
#include <SxFermi.h>
#include <SxPW.h>

#include <SxCLI.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();
   // --- init S/PHI/nX timers

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "P. Eggert, C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates an approximate photoemission spectrum.";
   SxString wavesFile = 
      cli.option ("-w|--waves","file","waves file with many unoccupied states")
      .toString ("waves.sxb");

   cli.newGroup ("nonlocal pseudopotential");
   bool nonlocalContributions =
      cli.option("--nl|--nonlocal","include non-local contributions").toBool ();
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file (for nl potentials)")
      .toString ("input.sx");

   cli.setGroup (cli.generalGroup);

   double eVac = cli.option ("--evac","energy","vacuum energy (eV)")
                .toDouble () / HA2EV;

   double hOmega = cli.option ("--omega","energy","excitation energy (eV)")
                .toDouble () / HA2EV;

   SxList<int> iks = cli.option ("--kpoints", "index list", 
                     "k-points as list, e.g. 1-4,7")
                     .toIdxList ();

   double gamma = cli.option ("--gamma","energy","Lorentz broadening (eV)")
                .toDouble () / HA2EV;
   
   double alpha = cli.option ("--alpha","energy","Gaussian broadening (eV)")
                .toDouble () / HA2EV;

   int nPoints = cli.option("--points", "number", "number of points to sample")
                 .toInt(100,2);

   SxList<int> iStates = cli.option ("--states|-n", "index list", 
                                     "states as list, e.g. 1-4,7")
                         .toIdxList ();
   cli.last ().defaultValue = "default: all occupied states";


   cli.version ("1.0");
   cli.finalize ();
   
   bool allOccupied = iStates.getSize () == 0;
   bool integrateK = iks.getSize () == 0;

   // --- read input

   SxPtr<SxGkBasis> gkBasisPtr;
   SxGBasis G;
   SxRBasis R;
   SxPW waves(wavesFile,SxPW::ReadOnDemand);
   SxFermi fermi;
   SxAtomicStructure structure;
   SxPerturbK kp;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io);
      waves.setGkBasisPtr (gkBasisPtr);
      fermi.read (io);
      fermi.kpPtr = &*gkBasisPtr;
      
      // Read cell including symmetries
      structure.read (io);
      
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxGkBasis &gkBasis = *gkBasisPtr;

   if (waves.getNSpin () == 2)  {
      cout << "No spin polarized version yet." << endl;
      SX_EXIT;
   }
   
   /*
   if (epsFile.getSize () > 0)  {
      int i, kpInFile, statesInFile, ik;
      int nSpin = 1, iSpin = 0;
      SxFermi::peekSpectrumFile (epsFile, &kpInFile, &statesInFile);
      SxFermi epsFromFile (0., statesInFile, nSpin, gkBasis);
      if (kpInFile != gkBasis.nk)  {
         cout << "'" << epsFile << "' has " << kpInFile;
         cout << " kpoints, but " << gkBasis.nk << " in " << wavesFile << endl;
         SX_QUIT;
      }
      epsFromFile.readSpectrum (epsFile, &structure.cell, &gkBasis);
      for (ik = 0; ik < gkBasis.nk; ik++)  {
         if (statesInFile < fermi.getNStates ())  {
            i = statesInFile - 1;
            double shift = epsFromFile.eps(i,iSpin,ik) 
                           - fermi.eps(i, iSpin, ik);
            cout << "ik = " << (ik + 1); 
            cout << ": shifting states " << (statesInFile + 1) << " to ";
            cout << fermi.getNStates () << " by " << (shift * HA2EV) << " eV\n";
            // shift higher DFT states like highest QP state known
            for (i = statesInFile; i < fermi.getNStates (); i++)
               fermi.eps(i, iSpin, ik) += shift;
         }
         // copy eps
         for (i = 0; i < statesInFile && i < fermi.getNStates (); i++)
            fermi.eps(i, iSpin, ik) = epsFromFile.eps (i,iSpin,ik);
      }
   }
   */

   // --- Recalculate Fermi distribution
//   if (cli.groupAvailable (fermiGroup) )  {
//      // get current number of electrons
//      double nEl = 0.;
//      for (int ik = 0; ik < fermi.getNk (); ik++)
//         nEl += fermi.focc(0/* iSpin */,ik).sum () * gkBasis.weights(ik);
//      cout << "nElectrons = " << nEl << endl;
//      fermi.nElectrons = nEl;
//      // set kpPtr
//      fermi.kpPtr = &gkBasis;
//      // recalculate
//      fermi.fermiDistribution (ekt);
//   }
   
   
   /*
   if (! fermi.isSemiconductor ()) {
      cout << endl << "This is not a semiconductor." << endl;
      cout << "The dielectric properties can be calculated only ";
      cout << "for semiconductors." << endl;
      SX_QUIT;
   }

   if (fermi.getNConductionBands (0,0) == 0)  {
      cout << endl << "ERROR: no unoccupied states." << endl;
      cout << "We need (many) unoccupied states, but there is none." << endl;
      SX_QUIT;
   }
   */

   // --- setup kp for pseudopotential
   if (nonlocalContributions)  {
      SxParser parser;
      SxParser::Table table = parser.read (inFile);
      SxPseudoPot psPot(&*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
   }

   // --- set up energy grid
   SxVector<Double> energy (nPoints), dEnergy;
   double eMin = -2. / HA2EV;          // little less than no kinetic energy
   double eMax = hOmega + 2. / HA2EV;  // all photon energy goes into kinetic energy
   int iE;
   for (iE = 0; iE < nPoints; iE++)
      energy(iE) = eMin + iE * (eMax - eMin) / double(nPoints - 1);

   int ik, ng, idir, i, iFinal;

   // --- compute
   PsiG wavesOcc, wavesVac, wavesRef, kpElements, kpAll;
   SxDiracVec<Double> epsVal, foccVal, deltaEps, lorentz;
   SxArray<SxVector<Double> > aKin(3); // idir, iE
   double ank, weight;
   int nVal, iSpin=0;
   int iFinalMin;
   if (integrateK)
      for (ik = 0; ik < gkBasis.getNk (); ++ik) iks.append (ik);

   for (int idxK = 0; idxK < iks.getSize (); idxK++)  {

      // --- collect initial states
      ik = iks(idxK);
      wavesRef = waves(0,ik);

      cout << "ik = " << (ik+1) << endl;
      ng = int(wavesRef.nRows ());
      if (allOccupied) {
         iStates.resize (0);
         for (i = 0; i < fermi.getNStates(ik); ++i)
            if (fermi.focc(i,iSpin,ik) > 1e-4) iStates.append(i);
      }
      nVal = int(iStates.getSize ());
      if (nVal == 0) continue;
      cout << nVal << " initial states." << endl;
      wavesOcc.reformat (ng, nVal);
      wavesOcc.setBasis (&gkBasis(ik));
      epsVal.resize (nVal);
      foccVal.resize (nVal);
      for (i = 0; i < nVal; i++)  {
         wavesOcc.colRef (i) <<= waves(iStates(i), 0, ik);
         epsVal(i) = fermi.eps(iStates(i),0,ik);
         foccVal(i) = allOccupied ? fermi.focc(iStates(i),0,ik) : 2.;
      }
      cout << "Collected initial states." << endl;
      
      // --- setup spectrum
      if (idxK == 0 || !integrateK)  {
         for (idir = 0; idir < 3; idir++)  {
            aKin(idir).resize (nPoints);
            aKin(idir) = 0.;
         }
      }
      
      // precalculate kp matrix elements
      cout << "Computing matrix elements..." << endl;
      for (iFinalMin = 0; iFinalMin < fermi.getNStates(ik); ++iFinalMin)
         if (fermi.eps(iFinalMin, iSpin, ik) >= eVac) break;
      
      {
         int nFinal =  int(wavesRef.nCols ()) - iFinalMin;
         PsiG wavesFinal = wavesRef(SxIdx(iFinalMin * ng, 
                                    int(wavesRef.getSize ())-1));
         wavesFinal.reshape (ng, nFinal);
         kpAll = kp.getMatrixElements (wavesFinal, wavesOcc);
         kpAll.reshape (nFinal, nVal * 3);
      }
      
      // --- loop over final states
      for (iFinal = iFinalMin; iFinal < fermi.getNStates (ik); iFinal++)  {
         if (fermi.eps(iFinal,0,ik) < eVac) continue;
         cout << "i=" << (iFinal+1) << "; eps= ";
         cout << fermi.eps(iFinal,0,ik) * HA2EV;
         // kpElements = kp.getMatrixElements (wavesOcc, waves(iFinal,0,ik));
         kpElements = kpAll.row(iFinal - iFinalMin);
         kpElements.reshape (nVal, 3);
         deltaEps = (fermi.eps(iFinal,0,ik) - hOmega) - epsVal;
         lorentz = (0.5 * gamma/PI)
                  / (deltaEps.sqr () + 0.25 * gamma * gamma);
         lorentz *= foccVal;
         dEnergy = fermi.eps(iFinal,0,ik) - eVac - energy;
         for (idir = 0; idir < 3; idir++)  {
            // Ank(iFinal, idir)
            ank = (lorentz * (kpElements.colRef(idir)).absSqr ()).sum ();
            if (! integrateK)
               cout << "; ank(" << idir << ")=" << ank; 

            weight = (2. - fermi.focc(iFinal,iSpin,ik))
                     / (sqrt(PI) * alpha); // Gaussian normalization
            if (integrateK) weight*= gkBasis.weights(ik);
            
//               aKin(idir) += (weight * ank * exp(-(dEnergy/alpha).sqr ()));
            aKin(idir).plus_assign_ax(weight * ank, exp(-(dEnergy/alpha).sqr ()));

         }
         cout << endl;
      }
      
      // output
      if ( (! integrateK) || (ik + 1 == gkBasis.getNk ()) )  {
         for (iE = 0; iE < nPoints; iE++)  {
            cout << (energy(iE) * HA2EV);
            for (idir = 0; idir < 3; idir++)
               cout << "\t" << (aKin(idir)(iE) / HA2EV);
            cout << "\n";
         }
      }
   }

}

