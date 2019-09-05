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
#include <SxPW.h>
#include <SxPerturbK.h>
#include <SxFermi.h>
#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxSpectrum.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "P. Eggert, C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates an approximate optical absorption spectrum.";
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

   SxList<int> iks = cli.option ("--kpoints", "index list", 
                     "k-points as list, e.g. 1-4,7")
                     .toIdxList ();

   int lorentzGroup = cli.newGroup ("Lorentz broadening");
   double gamma = cli.option ("--gamma","energy","Lorentz broadening (eV)")
                .toFloat () / HA2EV;
   int nPoints = cli.option("--points", "number", "number of points to sample")
                 .toInt(100,2);
   int gaussGroup = cli.newGroup ("Gaussian broadening");
   double broad = cli.option("--alpha", "energy", "Gaussian broadening (eV)")
                 .toFloat () / HA2EV;
   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);

   cli.setGroup (cli.generalGroup);
   bool useLorentz = cli.groupAvailable(lorentzGroup),
        useGauss   = cli.groupAvailable(gaussGroup);
   if (!(useLorentz || useGauss))  {
      cout << "Either Gaussian or Lorentzian broadening must be used." << endl;
      cli.setError ();
   }
   
   SxList<int> iStates = cli.option ("--initial|-I", "index list", 
                                     "initial states as list, e.g. 1-4,7")
                         .toIdxList ();
   cli.last ().defaultValue = "default: all occupied states";

   SxList<int> fStates = cli.option ("--final|-F", "index list", 
                                     "final states as list, e.g. 1-4,7")
                         .toIdxList ();
   cli.last ().defaultValue = "default: all unoccupied states";

   SxList<double> range = cli.option("--range","range","Emin:Emax (eV)")
                          .toDoubleList ();
   if (range.getSize () != 2 && !cli.error)  {
      cout << "Illegal range." << endl;
      cli.setError ();
   }

   bool debug = cli.option ("--details","request more output").toBool ();

   cli.version ("0.9");
   cli.finalize ();

   // --- read input


   SxPtr<SxGkBasis> gkBasisPtr;
   SxGBasis G;
   SxRBasis R;
   SxFermi fermi;
   SxAtomicStructure structure;
   SxPerturbK kp;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create(io);
      fermi.read (io);
      
      // Read cell including symmetries
      structure.read (io);
      
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   SxGkBasis &gkBasis = *gkBasisPtr;
   waves.setGkBasisPtr (gkBasisPtr);
   fermi.kpPtr = &*gkBasisPtr;

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

   if (nonlocalContributions)  {
      SxParser parser;
      SxParser::Table table = parser.read (inFile);
      SxPseudoPot psPot(&*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
   }

   // --- set up energies
   SxVector<Double> energy (nPoints), dEnergy;
   double eMin = range(0) / HA2EV;
   double eMax = range(1) / HA2EV;
   SxSpectrum gaussSpectrum;
   if (useGauss) gaussSpectrum.nPerPeak = nPerPeak;
   if (useLorentz)  {
      for (int iE = 0; iE < nPoints; iE++)
         energy(iE) = eMin + iE * (eMax - eMin) / double(nPoints - 1);
   }

   int ik, ng, idir, iInitial, iFinal;
   int nVal = fermi.getNValenceBands (0,0);
   int nInitial, nFinal;
   // --- default for initial states: all occupied
   if (iStates.getSize () == 0)  {
      for (int i = 0; i < nVal; i++)
         iStates.append(i);
   }
   nInitial = int(iStates.getSize ());
   
   // --- default for final states: all unoccupied
   if (fStates.getSize () == 0)  {
      for (int i = nVal; i < fermi.getNStates(); ++i)  {
         fStates.append (i);
      }
   }
   nFinal = int(fStates.getSize ());

   PsiG wavesRef, wavesInitial, wavesFinal, kpAllDir;
   SxArray<PsiG> kpElements(3);
   SxVector<Double> deltaEps, lorentz;
   SxVector<Double> gaussPeak(3);
   SxArray<SxVector<Double> > spectrum(3); // idir, iE
   double kpSquare, dEps;
   bool allK = iks.getSize () == 0;
   if (allK)
      for (ik = 0; ik < gkBasis.getNk (); ++ik) iks << ik;

   SxArray3<double> sym(3,3,3);
   if (allK)  {
      sym.set (0.);
      const SxSymGroup &syms = *structure.cell.symGroupPtr;
      int nSym = syms.getNSymmorphic ();
      for (int iSym = 0; iSym < nSym; ++iSym) {
         SymMat S = syms.getSymmorphic (iSym);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               for (int k = 0; k < 3; ++k)
                  sym(i,j,k) += S(i,j) * S(i,k) / double(nSym);
      }
      cout << "Symmetrizer=" << sym << endl;
   }

   for (int idxK = 0; idxK < iks.getSize (); idxK++)  {
      ik = iks(idxK);
      ng = gkBasis(ik).ng;

      // --- check epsSortIdx
      for (int i = 0; i < nInitial; ++i)
         if (iStates(i) != fermi.epsSortList (iStates(i),0,ik))
            cout << "Warning: " << (iStates(i)+1) << " is state "
                 << (fermi.epsSortList (iStates(i),0,ik))
                 << "in energetic order." << endl;
      for (int i = 0; i < nFinal; ++i)
         if (fStates(i) != fermi.epsSortList (fStates(i),0,ik))
            cout << "Warning: " << (fStates(i)+1) << " is state "
                 << (fermi.epsSortList (fStates(i),0,ik))
                 << "in energetic order." << endl;

      // --- set up initial states
      wavesInitial.reformat (ng, nInitial);
      wavesInitial.setBasis (&gkBasis(ik));
      for (int i = 0; i < nInitial; ++i)
         wavesInitial.colRef(i) <<= waves(iStates(i), 0, ik);
      // --- set up final states
      wavesFinal.reformat (ng, nFinal);
      wavesFinal.setBasis (&gkBasis(ik));
      for (int i = 0; i < nFinal; ++i)
         wavesFinal.colRef(i) <<= waves(fStates(i), 0, ik);
      
      // --- calculate k.p
      kpAllDir = kp.getMatrixElements (wavesInitial, wavesFinal);
      kpAllDir.reshape (nInitial * nFinal, 3);
      for (idir = 0; idir < 3; ++idir)  {
         kpElements(idir) = kpAllDir.colRef(idir);
         kpElements(idir).reshape (nInitial, nFinal);
      }

      // --- set up spectrum
      if (!allK || ik == 0)  {
         if (useLorentz)  {
            for (idir = 0; idir < 3; ++idir)  {
               spectrum(idir).resize (nPoints);
               spectrum(idir).set (0.);
            }
         } else {
            gaussSpectrum.set(eMin, eMax, broad, 3);
         }
      }
      
      // --- double sum over states
      for (int i = 0; i < nInitial; ++i)  {
         iInitial = iStates(i);
         for (int j = 0; j < nFinal; ++j)  {
            iFinal = fStates(j);
            dEps = fermi.eps(iFinal,0,ik) - fermi.eps(iInitial,0,ik);
            if (useLorentz)  {
               deltaEps = dEps - energy;
               lorentz = (0.5 * gamma/PI)
                        / (deltaEps.sqr () + 0.25 * gamma * gamma);
            }
            if (debug) 
               cout << "|<" << (iInitial+1) << "|p|" << (iFinal+1) << ">|^2 @"
                    << (dEps*HA2EV) << "=";
            
            for (idir = 0; idir < 3; idir++)  {
               if (allK)  {
                  kpSquare = 0.;
                  // symmetrization
                  for (int jdir = 0; jdir < 3; ++jdir)
                     for (int kdir = 0; kdir < 3; ++kdir)
                        kpSquare += sym(idir, jdir, kdir)
                                  * (kpElements(jdir)(i,j)
                                  *  kpElements(kdir)(i,j).conj ()).re;
                  kpSquare *= gkBasis.weights(ik);
               } else {
                  kpSquare = kpElements(idir)(i,j).absSqr ();
               }
               if (useLorentz)  {
                  // spectrum(idir) += kpSquare * lorentz;
                  spectrum(idir)
                     .plus_assign_ax(kpSquare / dEps, lorentz);
               } else {
                  gaussPeak(idir) = kpSquare / dEps;
               }
               if (debug) cout << "\t" << kpSquare;
            }
            if (useGauss) {
               gaussSpectrum.addPeak (dEps, gaussPeak);
            }
            if (debug) cout << endl;
         }
      }
      cout << "ik = " << (ik+1);
      cout << endl;
      
      // output
      if (!allK || ik+1 == gkBasis.getNk ())  {
         if (useLorentz)  {
            for (int iE = 0; iE < nPoints; iE++)  {
               cout << (energy(iE) * HA2EV);
               for (idir = 0; idir < 3; idir++)
                  cout << "\t" << (spectrum(idir)(iE) / HA2EV);
               cout << "\n";
            }
         } else {
            gaussSpectrum.compute ();
            gaussSpectrum.spectra  /= HA2EV;
            gaussSpectrum.energies *= HA2EV;
            gaussSpectrum.fprint(stdout);
         }
      }
   }

}

