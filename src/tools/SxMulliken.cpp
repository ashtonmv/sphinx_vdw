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
#include <SxKGrid.h>
#include <SxPseudoPot.h>
#include <SxPAWOverlap.h>
#include <SxPWOverlap.h>
#include <SxFermi.h>

SxDiracMat<TAOBasisType>
getDensityMatrix (const SxAOBasis          &ao,
                  const PsiG               &waves,
                  const Focc               &focc)
{
   int iSpin = waves.handle->auxData.iSpin;
   int ik = waves.handle->auxData.ik;

   int nStates = int(waves.nCols ());
   if (nStates == 0) nStates = 1;

   // --- get number of occupied states
   int nOccStates = 0;
   PrecFocc lastFocc = focc(0,iSpin,ik);
   bool sorted = true;
   for (int i = 0; i < nStates; ++i)  {
      if (focc(i,iSpin,ik) > lastFocc + 1e-6) sorted = false;
      if (fabs(lastFocc = focc(i,iSpin,ik)) > 1e-5) ++nOccStates;
   }
   SX_CHECK (nOccStates > 0);
   
   // --- get waves/focc for occupied states
   PsiG occWaves;
   SxDiracVec<TPrecFocc> myFocc;
   int ng = int(waves.nRows ());
   if (sorted)  {
      // standard case: packed storage
      occWaves = waves(SxIdx(0, ng * nOccStates - 1));
      occWaves.reshape (ng, nOccStates);
      occWaves.handle->auxData = waves.handle->auxData;
      myFocc = focc(iSpin,ik)(SxIdx(0,nOccStates-1));
   } else {
      if (2 * nOccStates < nStates)  {
         // pack occupied states
         occWaves.reformat (ng, nOccStates);
         occWaves.handle->auxData = waves.handle->auxData;
         myFocc.resize (nOccStates);
         int iState = 0;
         for (int i = 0; i < nStates; ++i)  {
            if (focc(i,iSpin,ik) > 1e-5)  {
               occWaves.colRef(iState) <<= waves.colRef(i);
               myFocc(iState++) = focc(i,iSpin,ik);
            }
         }
      } else {
         // do not care about the unnecessary states
         occWaves = waves;
         nOccStates = nStates;
         myFocc = focc(iSpin,ik);
      }
   }
   
   // get (S^-1) <mu|psi>
   SxAOBasis::TPsi aoPsi;
   //               S^{-1}           < mu | psi >
   aoPsi = ao.getInverseOverlap(ik) ^ (ao | occWaves);

   // release occWaves memory
   occWaves = PsiG ();

   // multiply aoPsi with sqrt(focc(iState))
   for (int i = 0; i < nOccStates; ++i)
      aoPsi.colRef (i) *= (myFocc(i) > 0. ? SxComplex16(1.) : I )
                        * sqrt(fabs(myFocc(i)));
   
   // return <mu|psi(i)>focc(i)<psi(i)|mu>
   return (aoPsi ^ aoPsi.adjoint ());
}

int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This is the Mulliken analysis tool.";

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   enum SpinMode { Alpha, Beta, Rho, Spin } mode = Rho;
   cli.option ("--alpha|--beta|--rho|--spin", "mode for spin");
   cli.last ().defaultValue = "rho";
   if (cli.last ().exists ())  {
      mode = SpinMode (cli.last ().toChoice ());
   };
   bool recompPhaseFactors
      = cli.option ("--recompute-phases",
                    "recompute phase factors: slower, but less memory")
        .toBool ();

   SxString basisFile = cli.option ("--basis", "file", "AOBasis file")
                        .toString ("");

   cli.version ("1.0");
   cli.finalize ();

   SxPtr<SxGkBasis> gkBasisPtr;
   SxFermi fermi;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      fermi.read (io);
      fermi.kpPtr = &*gkBasisPtr;
      structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   waves.setGkBasisPtr (gkBasisPtr);
   SxGkBasis &gkBasis = *gkBasisPtr;

   // read pseudopotential data
   cout << "Read input file" << endl;
   SxParser parser;
   SxParser::Table table = parser.read (inFile);

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
      SxPtr<SxPseudoPot> potPtr = SxPtr<SxPseudoPot>::create (&*table);
      radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
      sData = potPtr;
      pseudoRad = potPtr->rad;
      logDr = potPtr->logDr;
      pseudoPsi = potPtr->getPseudoPsi ();
      
      if (basisFile.getSize () == 0)  {
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
      radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
      pseudoRad = potPtr->rad;
      logDr = potPtr->logDr;
      pseudoPsi = potPtr->getPhiPS ();
      
      if (basisFile.getSize () == 0)  {
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
      } else {
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

   SxAOBasis::TPsi psiAO, psi;
   int iSpin = 0, nk = gkBasis.getNk ();
   int nOrb = aoBasis.getNOrb ();

   SxDiracMat<TAOBasisType> P, S;
   SxDiracMat<TReal8> M(nOrb, nOrb);
   M.set (0.);
   cout << "Computing full Mulliken matrix..." << endl;
   for (int ik = 0; ik < nk; ++ik)  {
      cout << "Computing S..." << endl;
      S = aoBasis.getOverlap(ik);
      cout << "Computing P..." << endl;
      for (iSpin = 0; iSpin < waves.getNSpin (); ++iSpin)  {
         if (mode == Alpha && iSpin == 1) continue;
         if (mode == Beta  && iSpin == 0) continue;
         P = getDensityMatrix (aoBasis, waves(iSpin, ik), fermi.focc);
         if (mode == Spin && iSpin == 1)  {
            M.plus_assign_ax(-gkBasis.weights(ik),S.conj () * P);
         } else {
            M.plus_assign_ax(gkBasis.weights(ik),S.conj () * P);
         }
      }
   }

   cout << "Condensing full Mulliken matrix to l-channels..." << endl;
   SxArray<SxArray<SxArray<int> > > mapSal;
   mapSal.resize (structure.getNSpecies ());
   int nSal = 0;
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

   SxArray<int> map(nOrb); // map orbital to m-summed orbitals
   for (int io = 0, isal = 0; io < nOrb; ++io)  {
     SxAOBasis::OrbitalIndex &sao = aoBasis.orbitalMap(io);
     map(io) = isal;
     //cout << io << "-> " << isal << endl;
     int l = aoBasis.refOrbMap(sao.is)(sao.io).l;
     int m = aoBasis.refOrbMap(sao.is)(sao.io).m;
     if (l == m) isal++;
   }

   SxDiracMat<Double> MsumM(nSal, nSal);
   MsumM.set (0.);
   
   {
      int iSal,jSal,iOrb,jOrb;
      SxAOBasis::OrbitalIndex orb;
      for (iOrb = 0; iOrb < nOrb; ++iOrb)  {
         iSal = map(iOrb);
         for (jOrb = 0; jOrb < nOrb; ++jOrb)  {
            jSal = map(jOrb);
            MsumM(iSal,jSal) += M(iOrb,jOrb);
         }
      }
   }

   // --- output
   {
      double totalN = 0.;
      SxString lNames="spdfghijklm";
      for (int is = 0, iSal = 0; is < structure.getNSpecies (); ++is)  {
         for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
            double charge = sData->valenceCharge(is), nMull;
            for (int iMu = 0; iMu < nMu(is); ++iMu, iSal++)  {
               charge -= (nMull = MsumM.colRef(iSal).sum ());
               cout << sData->chemName(is) << (ia + 1);
               cout << ' ' << lNames(lMu(is)(iMu)) << "-shell: "
                    << SxString(nMull, "%.3lf") << endl;;
               totalN += nMull;
            }
            cout << sData->chemName(is) << (ia + 1);
            if (mode == Rho)  {
               cout << SxString(charge, " Mulliken charge: %.3lf") << endl;
            } else {
               cout << " total occupation: " 
                    << SxString(sData->valenceCharge(is) - charge, "%.3lf")
                    << endl;
            }
         }
      }
      cout << "Total nEl " << SxString(totalN, "%.3lf") << endl;
   }

   return 0;

}
