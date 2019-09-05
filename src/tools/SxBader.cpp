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
#include <SxConfig.h>
#include <SxBinIO.h>
#include <SxTypes.h>
#include <SxMesh3D.h>
#include <SxNeighbors.h>
#include <SxRBasis.h>
#include <SxRho.h>
#include <SxPAWRho.h>
#include <SxProjector.h>
#include <SxPAWHamiltonian.h>
#include <SxSpaceSeparator.h>

SxVector<Double> computeNHatR (SxPAWPot &pawPot, SxPAWRho &pawRho, SxAtomicStructure &structure)
{
   int nSpin = pawRho.getNSpin ();
   int nAtoms = structure.getNAtoms ();
   SxVector<Double> result(nAtoms);
   result.set(0.0);
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      int ia;
      int is = structure.getISpecies(iAtom,&ia); 
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (int ipt = 0; ipt < pawPot.getNProjType (is); ++ipt)  {
            int offsetI = pawPot.offset(is)(ipt);
            for (int jpt = 0; jpt < pawPot.getNProjType (is); ++jpt)  {
               if (pawPot.lPhi(is)(ipt) == pawPot.lPhi(is)(jpt))  {
                  int offsetJ = pawPot.offset(is)(jpt);
                  int nm = 2 * pawPot.lPhi(is)(ipt) + 1;
                  for (int m = 0; m < nm; ++m)  {
                     result(iAtom) += 
                        pawRho.Dij(iSpin,is,ia)(offsetI + m, offsetJ + m)
                        * pawPot.deltaS(is)(ipt,jpt);
                  }
               }
            }
         }
      }
   }

   return result;
}

SxDiracVec<Double> computeRhoCore (SxPAWPot &pawPot, SxPAWRho &pawRho, 
      SxAtomicStructure &structure)
{
   const SxRBasis &R = pawRho.pwRho(0).getBasis<SxRBasis> ();
   const SxGBasis &G = R.getGBasis ();
   
   int nSpecies = structure.getNSpecies ();

   SxDiracVec<Complex16> rhoCoreG(G);
   rhoCoreG.set(0.0);
   for (int is = 0; is < nSpecies; is++)  {
      SxDiracVec<Double> rhoPseudoCore = pawPot.rhoCorePS(is);
      rhoPseudoCore.handle->auxData.is = is;
      rhoPseudoCore.handle->auxData.l = 0;
      rhoPseudoCore.handle->auxData.m = 0;
      SxDiracVec<Complex16> rhoPseudoCoreG = (G | rhoPseudoCore);
      rhoCoreG += G.structureFactors(is) * rhoPseudoCoreG;
   }
   
   SxDiracVec<Double> result = R.symmetrize(R | rhoCoreG);

   RhoR density(1);
   density(0) = 1.0 * result;
   SxString file = "rhoCore.sxb";
   SxRho(density).writeRho(file);

   return result;
}

SxVector<Double> computeNPseudoCore (SxPAWPot &pawPot,SxAtomicStructure &structure)
{
   const SxRadBasis &rad = pawPot.getRadBasis ();
   int nAtoms = structure.getNAtoms ();
   int nSpecies = structure.getNSpecies ();

   SxVector<Double> result(nAtoms);

   for (int is = 0; is < nSpecies; is++)  {
      SxDiracVec<Double> r = rad.radFunc(is);
      double nPseudoCore = sqrt(FOUR_PI) 
         * (pawPot.rhoCorePS(is) * r.cub ())
         .integrate(rad.logDr(is));
      for (int ia = 0; ia < structure.getNAtoms(is); ia++)  {
         int iAtom = structure.getIAtom(is,ia);
         result(iAtom) = nPseudoCore;
      }
   }
   
   return result;
}

int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "B. Lange";
   cli.preUsageMessage = "This is the Bader analysis tool.";

   SxString rhoFile = cli.option ("-r|--rho","file","rho file")
                        .toString("rho.sxb");

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   bool writeSxb = cli.option ("--printSxb", "flag", "print filters")
                     .toBool ();


   cli.version ("1.0");
   cli.finalize ();

   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   SxAtomicStructure structure = SxAtomicStructure (&*table);

   SxVector3<Int> dims;
   SxBinIO io;
   SxString fileName;
   io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
   io.read("dim", &dims);
   sxprintf("Mesh: %i x %i x %i\n",dims(0), dims(1), dims(2));
   SxMesh3D mesh(dims);
   SxRBasis R(mesh, structure.cell);
   SxGBasis G(mesh, structure, 4. * SxGBasis::getECut(&*table));
   R.registerGBasis (G);
   G.registerRBasis (R);

   int nAtoms = structure.getNAtoms();
   int nSpecies = structure.getNSpecies ();
   SxVector<Double> nValence (nAtoms);
   SxVector<Double> nPseudoCore (nAtoms);
   SxVector<Double> nComp (nAtoms);
   SxVector<Double> compRadii (nSpecies);
   SxArray<SxCubicSpline<SxDiracVec<Double> > > coreAE(nSpecies);
   SxDiracVec<Double> rhoPseudoTot (R);
   SxDiracVec<Double> rhoPseudoCore (R);
   SxDiracVec<Double> rhoCompensation (R);
   SxDiracVec<Double> rhoPseudoValence (R);

   SxSpaceSeparator separator (structure);

   // Normconserving or PAW Potential ?
   if (table->containsGroup("pseudoPot"))   {
      SxPtr<SxPseudoPot> potPtr = SxPtr<SxPseudoPot>::create (&*table);
      SxRho rho (io, &R);
      int nSpin = rho.getNSpin ();
      rhoPseudoTot = (nSpin == 2) ? rho.rhoR(0) + rho.rhoR(1) : rho.rhoR(0);
      rhoPseudoTot.setBasis(R);

      for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
         int is = structure.getISpecies(iAtom);
         nValence(iAtom) = potPtr->valenceCharge(is);
      }
      nComp.set(0.0);
      rhoPseudoCore.set(0.0);
      nPseudoCore.set(0.0);
      cout << "Total charge is " << rho.getNorm () << endl;
   } else if (table->containsGroup("pawPot"))   {
      SxPtr<SxPAWPot> potPtr =SxPtr<SxPAWPot>::create (&*table);
      SxPAWRho rho (potPtr);
      // pw density need R Basis
      rho.pwRho = SxRho (io, &R);
      // setup paw density
      rho.readRho(rhoFile);
      int nSpin = rho.getNSpin ();
      rhoPseudoTot = (nSpin == 2) ? rho.pwRho.rhoR(0) + rho.pwRho.rhoR(1) 
                          : rho.pwRho.rhoR(0);
      rhoPseudoTot.setBasis(R);
      rhoPseudoCore = computeRhoCore (*potPtr, rho, structure);
      nPseudoCore = computeNPseudoCore(*potPtr, structure);
      nComp = computeNHatR (*potPtr, rho, structure);

      for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
         separator.PAWDist(iSpecies) = potPtr->rCore(iSpecies);
         compRadii(iSpecies) = potPtr->rc(iSpecies);
         for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
            int iTlAtom = structure.getIAtom(iSpecies,iAtom);
            nValence(iTlAtom) = potPtr->valenceCharge(iSpecies);
         }
      }
      
   } else   {
      cout << "No known Potential Group found !" << endl;
      SX_QUIT;
   }

   rhoPseudoValence = rhoPseudoTot - rhoPseudoCore;

   
   cout << SX_SEPARATOR;
   cout << "Tot pseudo electrons =  " << rhoPseudoTot.sum () * R.dOmega << endl;
   cout << "Core charges are ";
   nPseudoCore.print ();
   cout << "Comp charges are "; 
   nComp.print ();
   cout << "Rho valence electrons are " 
        << rhoPseudoTot.sum () * R.dOmega + nComp.sum () - nPseudoCore.sum ()
        << endl;
   cout << "Compensation radii are ";
   compRadii.print();
   cout << "PAW Dists are ";
   separator.PAWDist.print();
   
   cout << SX_SEPARATOR;

   SxArray<SxDiracVec<Double> > voronoi = separator.voronoi (R);

   cout << SX_SEPARATOR;
   cout << "Voronoi charges" << endl;
   cout << SX_SEPARATOR;
   double nSum = 0.0;
   double nPseudoElec = 0.0;
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      nSum += voronoi(iAtom).sum ();
      nPseudoElec += (voronoi(iAtom) * rhoPseudoTot).sum () * R.dOmega;
      double charge = (voronoi(iAtom) * rhoPseudoTot).sum () * R.dOmega
         - nPseudoCore(iAtom)
         + nComp(iAtom)
         - nValence(iAtom);
      sxprintf("Atom %i: charge = %2.8f e-, nPoints = %.2f\n",
            iAtom, charge, voronoi(iAtom).sum ());
      if (!writeSxb) continue;
      fileName = "voronoiFilter-" + SxString(iAtom) + ".sxb";
      io.open      (fileName, SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (voronoi(iAtom), structure.cell, mesh);
      io.setMode   (SxBinIO::WRITE_DATA);
      io.writeMesh (voronoi(iAtom), structure.cell, mesh);
      io.close();
   }
   sxprintf("nPoints(TOT) = %.2f\n",nSum);
   sxprintf("nPseudoElec(TOT) = %.2f\n",nPseudoElec);
   cout << SX_SEPARATOR;

   cout << SX_SEPARATOR;
   cout << "Bader charges" << endl;
   cout << SX_SEPARATOR;
   
   SxArray<SxDiracVec<Double> > baderTrinkle 
      = separator.baderTrinkle(rhoPseudoTot);
   
   if (writeSxb)  {
      for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
         SxDiracVec<Double> partialRho 
            = baderTrinkle(iAtom) * rhoPseudoTot;
         fileName = "baderTrinkleRho-" + SxString(iAtom) + ".sxb";
         io.open      (fileName, SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (partialRho, structure.cell, mesh);
         io.setMode   (SxBinIO::WRITE_DATA);
         io.writeMesh (partialRho, structure.cell, mesh);
         io.close();
      }
   }
   nSum = 0.0;
   nPseudoElec = 0.0;
   nValence.print ();
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      nSum += baderTrinkle(iAtom).sum ();
      nPseudoElec += (baderTrinkle(iAtom) * rhoPseudoTot).sum () * R.dOmega;
      double charge = (baderTrinkle(iAtom) * rhoPseudoTot).sum () * R.dOmega 
         - nPseudoCore(iAtom)
         + nComp(iAtom) 
         - nValence(iAtom);
         if (writeSxb) {
            fileName = "baderTrinkleFilter-" + SxString(iAtom) + ".sxb";
            io.open      (fileName, SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (baderTrinkle(iAtom), structure.cell, mesh);
            io.setMode   (SxBinIO::WRITE_DATA);
            io.writeMesh (baderTrinkle(iAtom), structure.cell, mesh);
            io.close();
         }
      sxprintf("Atom %i: charge = %2.8f e- nPoints = %.2f\n", iAtom, charge,
            baderTrinkle(iAtom).sum ());
   }
   sxprintf("nPoints(TOT) = %.2f\n",nSum);
   sxprintf("nPseudoElec(TOT) = %.2f\n",nPseudoElec);
   sxprintf("nPseudoElec(TOT) = %.2f\n",rhoPseudoTot.sum () * R.dOmega);

   SxTimer::getGlobalTimer().print ();
   cout << "SxBader completes successfully" << endl;

   return 0;

}
