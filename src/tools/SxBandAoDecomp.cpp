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
#include <SxProjector.h>
#include <SxAOBasis.h>
#include <SxPseudoPot.h>
#include <SxPAWPot.h>
#include <SxFermi.h>
#include <SxPWOverlap.h>
#include <SxPAWOverlap.h>


int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This is the Mulliken-like band analysis tool.";

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   SxList<int> states, kPoints;
   // --- selection by state number
   states = cli.option ("-n|--states","list",
                        "states (starting from 1) to be selected. "
                        "They are given as a ','-separated list of numbers or "
                        "ranges a-b, e.g. '-n 1,5-8,10' would select the states"
                        " 1, 5, 6, 7, 8, and 10.")
            .toIdxList ();
   cli.last ().defaultValue = "default: all states";
   kPoints = cli.option ("-k|--kpoints","list",
                         "k-points (starting from 1) to be "
                         "selected. The list is given like for the -n option")
             .toIdxList ();
   cli.last ().defaultValue = "default: all k-points";

   bool readStructFromInput = cli.option("--counterpoise",
         "Read structure from input file (allows for "
         "counterpoise-like corrections of the BSSE)").toBool ();

   double filter = cli.option("--filter","min occ",
                              "set minimum occupation to print")
                   .toDouble (0.01);

   SxString outFile = cli.option("-o|--output","file","(packed) output file")
                      .toString ("");

   bool recompPhaseFactors
      = cli.option ("--recompute-phases",
                    "recompute phase factors: slower, but less memory")
        .toBool ();
   bool averageK
      = cli.option ("--averageK",
                    "average over k-points")
        .toBool ();

   SxString basisFile = cli.option ("--basis", "file", "AOBasis file")
                        .toString ("");

   cli.version ("0.2");
   cli.finalize ();

   SxPtr<SxGkBasis> gkBasisPtr;
   SxFermi fermi;

   FILE *output = NULL;
   if (outFile.getSize () > 0)  {
      if ((output = fopen(outFile.ascii (), "w")) == NULL)  {
         cout << "Can't open output file '" << outFile << "'." << endl;
         SX_QUIT;
      }
   }
   SxAtomicStructure structure;
   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      if (!readStructFromInput) structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxPW waves(wavesFile, SxPW::ReadOneByOne);
   waves.setGkBasisPtr (gkBasisPtr);
   SxGkBasis &gkBasis = *gkBasisPtr;

   // read pseudopotential data
   cout << "Read PP data" << endl;
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   if(readStructFromInput) structure = SxAtomicStructure(&*table);

   // setup phase factors
   cout << "Set up phase factors" << endl;
   if (recompPhaseFactors)
      for (int ik = 0; ik < gkBasis.getNk (); ++ik)
         gkBasis(ik).memMode = SxGBasis::SaveMemory;
   gkBasis.changeTau(structure);

   cout << "Setting up ao basis" << endl;
   SxArray<SxDiracVec<TReal8> > pseudoRad;
   SxArray<Real8> logDr;
   SxArray<SxArray<SxDiracVec<TReal8> > > pseudoPsi;
   SxVector<Int> nMu;
   SxArray<SxVector<Int> > lMu(structure.getNSpecies ());

   SxConstPtr<SxRadBasis> radBasisPtr;
   SxPtr<SxSpeciesData> sData;
   SxPtr<SxOverlapBase> SPtr;
   SxPtr<SxAOBasis> aoBasisPtr;

   // Normconserving or PAW Potential ?
   if (table->containsGroup("pseudoPot"))   {
      SPtr = SxPtr<SxPWOverlap>::create ();

      if (basisFile.getSize () == 0)  {
         SxPtr<SxPseudoPot> potPtr = SxPtr<SxPseudoPot>::create (&*table);
         sData = potPtr;
         radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
         pseudoRad = potPtr->rad;
         logDr = potPtr->logDr;
         pseudoPsi = potPtr->getPseudoPsi ();
         nMu = potPtr->lMax + 1;
         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            lMu(is).resize (nMu(is));
            for (int l = 0; l < nMu(is); ++l)  {
               lMu(is)(l) = l;
            }
         }
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr, pseudoPsi,
                                         SPtr);
      } else {
         sData = SxPtr<SxSpeciesData>::create (&*table);
         SxParser aoParser;
         SxConstPtr<SxSymbolTable> aoTable = aoParser.read(basisFile);
         SxSymbolTable *aoGroup = aoTable->getGroup("AOBasis");
         SxAtomicOrbitals orbs;
         radBasisPtr = SxPtr<SxRadBasis>::create(aoGroup);
         orbs.setup(aoGroup);
         nMu.resize (structure.getNSpecies ());
         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            nMu(is) = orbs.getNOrbTypes (is);
            lMu(is).resize (nMu(is));
            for (int iMu = 0; iMu < nMu(is); ++iMu)  {
               lMu(is)(iMu) = orbs.getL(is, iMu);
            }
         }
         radBasisPtr = orbs.getRadBasisPtr ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, orbs);
      }
   } else if (table->containsGroup("pawPot"))   {
      SxPtr<SxPAWPot> potPtr =SxPtr<SxPAWPot>::create (&*table);
      sData = potPtr;
      SxPtr<SxPartialWaveBasis> pBasis 
         = SxPtr<SxPartialWaveBasis>::create (potPtr, structure);
      pBasis->createProjBasis (gkBasis);
      SPtr = SxPtr<SxPAWOverlap>::create (pBasis, potPtr);
      if (basisFile.getSize () == 0)  {
         radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
         pseudoRad = potPtr->rad;
         logDr = potPtr->logDr;
         pseudoPsi = potPtr->getPhiPS ();
         nMu.resize (structure.getNSpecies ());
         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            nMu(is) = potPtr->getNProjType (is);
            lMu(is).resize (nMu(is));
            for (int iMu = 0; iMu < nMu(is); ++iMu)  {
               lMu(is)(iMu) = potPtr->lPhi(is)(iMu);
            }
         }
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr, pseudoPsi,
                                         SPtr);
      } else  {
         SxParser aoParser;
         SxConstPtr<SxSymbolTable> aoTable = aoParser.read(basisFile);
         SxSymbolTable *aoGroup = aoTable->getGroup("AOBasis");
         SxAtomicOrbitals orbs;
         radBasisPtr = SxPtr<SxRadBasis>::create(aoGroup);
         orbs.setup(aoGroup);
         nMu.resize (structure.getNSpecies ());
         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            nMu(is) = orbs.getNOrbTypes (is);
            lMu(is).resize (nMu(is));
            for (int iMu = 0; iMu < nMu(is); ++iMu)  {
               lMu(is)(iMu) = orbs.getL(is, iMu);
            }
         }
         radBasisPtr = orbs.getRadBasisPtr ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, orbs, SPtr);
      }
   } else   {
      cout << "No known Potential Group found !" << endl;
      SX_QUIT;
   }

   SxAOBasis &aoBasis = *aoBasisPtr;

   // setup radial basis and aoBasis
   aoBasis.setOverlapCaching (SxAOBasis::CacheCurrentK);
   aoBasis.setInvOverlapCaching (SxAOBasis::CacheCurrentK);

   // --- setup map (is)(ia)(l) -> isal
   SxArray<SxArray<SxArray<int> > > mapSal;
   mapSal.resize (structure.getNSpecies ());
   int nSal = 0;
   int nOrb = aoBasis.getNOrb ();
   {
      int is, ia, nAtoms, nl, l;
      for (is = 0; is < structure.getNSpecies (); ++is)  {
         nAtoms = structure.getNAtoms(is);
         nl = nMu(is);
         mapSal(is).resize (nAtoms);
         for (ia = 0; ia < nAtoms; ++ia)  {
            mapSal(is)(ia).resize(nl);
            {
               SxArray<int> &a = mapSal(is)(ia);
               for (l = 0; l < nl; ++l) a(l) = nSal++;
            }
         }
      }
   }

   // --- setup map io -> isal
   SxArray<int> map(nOrb); // map orbital to m-summed orbitals
   for (int io = 0, isal = 0; io < nOrb; ++io)  {
     SxAOBasis::OrbitalIndex &sao = aoBasis.orbitalMap(io);
     map(io) = isal;
     //cout << io << "-> " << isal << endl;
     int l = aoBasis.refOrbMap(sao.is)(sao.io).l;
     int m = aoBasis.refOrbMap(sao.is)(sao.io).m;
     if (l == m) isal++;
   }

   int iSpin = 0, nk = gkBasis.getNk ();
   int nStates = waves.getNStates ();
   PsiG psi, psiAO, allPsi, allPsiAO;
   SxDiracVec<TReal8> normSal(nSal);

   int idxK, ik, idxState, iState;

   if (kPoints.getSize () == 0)
      for (ik = 0; ik < nk; ++ik) kPoints << ik;
   if (states.getSize () == 0)
      for (iState = 0; iState < nStates; ++iState) states << iState;

   if (!averageK)  {
      for (idxK = 0; idxK < kPoints.getSize (); ++idxK)  {
         ik = kPoints(idxK);
         cout << "ik = " << (ik+1) << endl;
         allPsi.reformat(gkBasis(ik).ng, states.getSize ());
         allPsi.setBasis(&gkBasis(ik));
         for (idxState = 0; idxState < states.getSize (); ++idxState)
            allPsi.colRef(idxState) <<= waves(states(idxState),iSpin,ik);
         allPsiAO = (aoBasis | allPsi);
         allPsiAO.reshape(nOrb, states.getSize ());

         for (idxState = 0; idxState < states.getSize (); ++idxState)  {
            iState = states(idxState);
            cout << "i = " << (iState + 1) << endl;
            if (output) fprintf(output,"%u",iState+1);
            // psi = waves(iState,iSpin,ik);
            // psiAO = (aoBasis | psi);
            psiAO = allPsiAO.colRef(idxState);
            normSal.set (0.);
            {
               PsiG normOrb;
               normOrb = 
                  psiAO.conj () * (aoBasis.getInverseOverlap (ik) ^ psiAO);
               cout << normOrb.real () << endl;
               for (int io = 0; io < nOrb; ++io)
                  normSal(map(io)) += normOrb(io).re; // imaginary parts sum to 0
            }
            // --- output
            {
               int is, ia, iMu, iSal = 0;
               SxString lNames="spdfghijklm";
               double shellNorm;
               for (is = 0; is < structure.getNSpecies (); ++is)  {
                  for (ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                     for (iMu = 0; iMu < nMu(is); ++iMu, iSal++)  {
                        shellNorm = normSal(iSal);
                        if (output) fprintf(output,"\t%8.6f",shellNorm);
                        if (fabs(shellNorm) >= filter)  {
                           cout << sData->chemName(is) << (ia + 1);
                           cout << ' ' << lNames(lMu(is)(iMu)) << "-shell: "
                              << shellNorm << endl;
                        }
                     }
                  }
               }
               cout << "Total AO norm: " << normSal.sum () << endl;
               if (output) fprintf(output,"\n");

            }

         }
      }
   } else  {
      for (idxState = 0; idxState < states.getSize (); ++idxState)  {
         iState = states(idxState);
         cout << "i = " << (iState + 1) << endl;
         if (output) fprintf(output,"%u",iState+1);
         normSal.set (0.);
         for (idxK = 0; idxK < kPoints.getSize (); ++idxK)  {
            ik = kPoints(idxK);
            psi = waves(states(idxState),iSpin,ik);
            psiAO = (aoBasis | psi);
            PsiG normOrb =
                  psiAO.conj () * (aoBasis.getInverseOverlap (ik) ^ psiAO);
               for (int io = 0; io < nOrb; ++io)
                  normSal(map(io)) += gkBasis.weights(ik) * normOrb(io).re;
         }
         // --- output
         {
            int is, ia, iMu, iSal = 0;
            SxString lNames="spdfghijklm";
            double shellNorm;
            for (is = 0; is < structure.getNSpecies (); ++is)  {
               for (ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                  for (iMu = 0; iMu < nMu(is); ++iMu, iSal++)  {
                     shellNorm = normSal(iSal);
                     if (output) fprintf(output,"\t%8.6f",shellNorm);
                     if (fabs(shellNorm) >= filter)  {
                        cout << sData->chemName(is) << (ia + 1);
                        cout << ' ' << lNames(lMu(is)(iMu)) << "-shell: "
                           << shellNorm << endl;
                     }
                  }
               }
            }
            cout << "Total AO norm: " << normSal.sum () << endl;
            if (output) fprintf(output,"\n");
         }
      }
   }

   if (output) fclose(output);

   return 0;

}

