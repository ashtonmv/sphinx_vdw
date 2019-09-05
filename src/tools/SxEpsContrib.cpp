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
#include <SxBinIO.h>
#include <SxProjector.h>
#include <SxPWHamiltonian.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates the various contributions to the eigenvalues. "
      "The Hamiltonian can be split into the following contributions:\n"
      "\n   H = T + V_el + V_xc + V_nl \n\n"
      "What is calculated by this add-on is the expectation value of each "
      "of these operators for each wavefunctions, which sum to the actual "
      " eigenvalue (which equals the expectation value of the total "
      "Hamiltonian). This separation may be physically misleading, but it "
      "looks impressive ;-).";
   SxString rhoFile = 
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxString wavesFile =
      cli.option ("-w|--waves","file","waves file")
      .toString ("waves.sxb");

   SxList<int> states =
      cli.option ("-n|--states","list","states to consider")
      .toIdxList ();
   cli.last ().defaultValue = "(default: all states)";

   cli.version ("0.1");
   cli.finalize ();

   // The Hamiltonian
   SxPWHamiltonian ham;
   
   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   
   // read xc information
   ham.read (&*table);
   // read pseudopotential
   ham.psPot = SxPseudoPot(&*table);

   // --- read waves file
   SxPW waves(wavesFile, SxPW::ReadOneByOne);
   SxBinIO io;
   try  {
      io.open (wavesFile, SxBinIO::BINARY_READ_ONLY);
      waves.setGkBasisPtr (SxPtr<SxGkBasis>::create (io));
      ham.structure.read(io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   ham.wavesPtr = &waves;
   SxGkBasis &gkBasis = waves.getGkBasis();
   
   // --- read density file
   SxVector3<Int> mesh;
   int &nSpin = ham.xcPtr->nSpin;

   // read mesh ...
   try  {
      io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
      io.read("dim", &mesh);
      nSpin = io.getDimension("nMeshes");
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   
   // --- setup G and R bases
   SxRBasis R(mesh, ham.structure.cell);
   SxGBasis G(mesh, ham.structure,
              SxGBasis::getGCut(SxGBasis::getECut(&*table)));
   G.registerRBasis (R);
   R.registerGBasis (G);
   ham.gBasisPtr = &G;
   ham.rho.rBasisPtr = &R;

   // ... and read rho
   try  {
      ham.rho.rhoR.resize (nSpin);
      ham.rho.readRho (io);
      io.close ();
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   // --- setup necessary parts of Hamiltonian
   ham.vEffR.resize (nSpin);
   int nR = int(ham.rho(0).getSize ());
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      ham.vEffR(iSpin).resize (nR);

   ham.computePhiLocPseudo(G);
   gkBasis.changeTau (ham.structure);
   ham.computePhiNonLocal (gkBasis);
   ham.computePhiGauss(G);
   ham.computeRhoGauss(G);

   ham.xcPtr->init (); 
   if (ham.psPot.nlcc) 
      ham.xcPtr->computePhiCore (ham.psPot, &G);
   if (ham.contrib & SxPWHamiltonian::CALC_X) ham.xcPtr->enableExchange ();
   if (ham.contrib & SxPWHamiltonian::CALC_C) ham.xcPtr->enableCorrelation ();

   SxPWHamiltonian::Contrib fullContrib = ham.contrib;

   ham.contrib = SxPWHamiltonian::CALC_EFF;
   ham.calcForces = false;
   ham.compute (SxFermi (), true);
   
   SxMeshR vElStat = ham.vHartreeR + ham.vLocR;

   // --- evaluate expectation values
   if (states.getSize () == 0)  {
      // --- all states
      for (int i = 0; i < waves.getNStates (); ++i) states << i;
   }
   int ik, iSpin, iBand;
   SxList<int>::Iterator iState;
   int nk = waves.getNk ();
   int nBands = int(states.getSize ());
   PsiG psi;
   SxFermi kin   (1.,nBands,nSpin,gkBasis),
           elstat(1.,nBands,nSpin,gkBasis),
           xc    (1.,nBands,nSpin,gkBasis), 
           nl    (1.,nBands,nSpin,gkBasis);
   for (ik = 0; ik < nk; ++ik)  {
      for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (iState = states.begin (), iBand = 0; 
              iState != states.end (); 
              ++iState,++iBand)  {
            if (*iState >= waves.getNStates(ik)) continue;
            psi = waves(*iState, iSpin, ik);
            // --- kinetic energy
            ham.contrib = SxPWHamiltonian::CALC_KIN;
            kin.eps(iBand, iSpin, ik) = (psi ^ (ham * psi)).chop ();
            // --- electrostatic (incl local potential) energy
            ham.contrib = static_cast<SxPWHamiltonian::Contrib>
                        (SxPWHamiltonian::CALC_HARTREE
                         + SxPWHamiltonian::CALC_LOC);
            ham.vEffR(iSpin) = vElStat;
            elstat.eps(iBand, iSpin, ik) = (psi ^ (ham * psi)).chop ();
            // --- xc energy
            ham.contrib = static_cast<SxPWHamiltonian::Contrib>
                        (SxPWHamiltonian::CALC_XC & fullContrib);
            ham.vEffR(iSpin) = ham.xcPtr->vXc(iSpin);
            xc.eps(iBand, iSpin, ik) = (psi ^ (ham * psi)).chop ();
            // --- nonlocal pseudopotential energy
            ham.contrib = SxPWHamiltonian::CALC_NL;
            nl.eps(iBand, iSpin, ik) = (psi ^ (ham * psi)).chop ();
            (cout << '.').flush ();
         }
         (cout << '\n').flush ();
      }
   }
   kin.writeSpectrum("eps_kin","dat");
   elstat.writeSpectrum("eps_elstat","dat");
   xc.writeSpectrum("eps_xc","dat");
   nl.writeSpectrum("eps_nl","dat");
}

