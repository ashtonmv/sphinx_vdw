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
      "This add-on calculates the electrostatic, exchange-correlation, and "
      "effective potentials.";
   SxString rhoFile = 
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxFFT::quickFFTPlanner ();
   SxFFT::plannerCLI (cli);
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxString structFile =
      cli.option ("-s|--structure","file","S/PHI/nX input file for structure")
      .toString ("");
   cli.last ().defaultValue = 
      "default: from input file = no structure optimization";

   bool elStat = 
      cli.option ("--elstat", "write electrostatic potential to vElStat.sxb")
      .toBool ();

   bool xc = 
      cli.option ("--xc", "write exchange-correlation potential to vXC.sxb")
      .toBool ();
   
   int xcDebugGroup = cli.newGroup ("Exchange or correlation only");
   bool onlyX = 
      cli.option ("--only-exchange", "write exchange potential to vXC.sxb")
      .toBool ();
   
   bool onlyC = 
      cli.option ("--only-correlation", 
                  "write correlation potential to vXC.sxb").toBool ();
   cli.newGroup("Effective potential");
   cli.excludeGroup (xcDebugGroup);

   bool eff = cli.option ("--eff", "write effective potential to vEff.sxb")
              .toBool ();
   cli.setGroup (cli.generalGroup);

   bool enforceDipole =
      cli.option ("--dipole","use dipole correction").toBool ();

   cli.version ("1.1");
   cli.finalize ();

   if (onlyX || onlyC) xc = true;
   if (onlyX && onlyC)  {
      cout << "Cannot write both only exchange and only correlation." << endl;
      SX_QUIT;
   }

   // The Hamiltonian
   SxPWHamiltonian ham;
   
   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   
   // read structure
   if (structFile.getSize () == 0)
      ham.structure = SxAtomicStructure(&*table);
   else
      ham.structure = SxAtomicStructure(&*
                      SxParser ().read (structFile,"std/structure.std") );
   // read pseudopotential
   ham.psPot = SxPseudoPot(&*table);

   // --- read density file
   SxVector3<Int> mesh;
   int &nSpin = ham.xcPtr->nSpin;
   SxBinIO io;

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
   // read xc information
   ham.read (&*table);
   // switch on dipole correction if wanted
   if (enforceDipole) ham.dipoleCorrection = true;
   
   // --- setup necessary parts of Hamiltonian
   ham.vEffR.resize (nSpin);
   int nR = int(ham.rho(0).getSize ());
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      ham.vEffR(iSpin).resize (nR);
   
   ham.computeESelf ();
   ham.computePhiLocPseudo(G);
   ham.computePhiGauss(G);
   ham.computeRhoGauss(G);

   if (xc || eff)  {
      ham.xcPtr->init (); 
      if (ham.psPot.nlcc) 
         ham.xcPtr->computePhiCore (ham.psPot, &G);
      if (!onlyC) ham.xcPtr->enableExchange ();
      if (!onlyX) ham.xcPtr->enableCorrelation ();
   } else {
      // avoid access to uninitialized values in compute
      ham.xcPtr->nlcc = false;
      ham.xcPtr->eXc = 0.;
   }
   
   // --- calculate potentials
   ham.contrib = ham.CALC_NONE;
   ham.calcForces = false;

   SxFermi fermi; // just a placeholder
   if (eff)  {
      // --- effective potential
      ham.contrib = ham.CALC_EFF;
      ham.compute (fermi, true);
      try  {
         // --- write vEff
         io.open("vEff.sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (ham.vEffR, R.cell, mesh);
         io.setMode (SxBinIO::WRITE_DATA);
         io.writeMesh (ham.vEffR, R.cell, mesh);
         io.close ();

         // --- write vXC
         if (xc)  {
            io.open ("vXC.sxb", SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (ham.xcPtr->vXc, R.cell, mesh);
            io.setMode (SxBinIO::WRITE_DATA);
            io.writeMesh (ham.xcPtr->vXc, R.cell, mesh);
            io.close ();
         }
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else if (xc)  {
      // --- xc only
      if (onlyC)
         ham.contrib = ham.CALC_C;
      else if (onlyX)
         ham.contrib = ham.CALC_X;
      else
         ham.contrib = ham.CALC_XC;
      ham.compute (fermi, true);
      // --- write out
      try  {
         io.open ("vXC.sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (ham.xcPtr->vXc, R.cell, mesh);
         io.setMode (SxBinIO::WRITE_DATA);
         io.writeMesh (ham.xcPtr->vXc, R.cell, mesh);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }

   if (elStat)  {
      // --- electrostatic only
      ham.contrib = (SxPWHamiltonian::Contrib)(ham.CALC_HARTREE + ham.CALC_LOC);
      ham.compute (fermi, true);
      // --- write out
      try  {
         io.open ("vElStat.sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (ham.vEffR, R.cell, mesh);
         io.setMode (SxBinIO::WRITE_DATA);
         io.writeMesh (ham.vEffR, R.cell, mesh);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }
   ham.printHartree = true;
   ham.printEnergies ();

   return 0;
}

