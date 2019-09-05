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
// Conversion program for GW program. 
// Reads S/PHI/nX files and generates input for GW
// AUTHOR: C. Freysoldt
//----------------------------------------------------------------------------

#include <SxDirac.h>
#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxFFT3d.h>
#include <SxList.h>
#include <SxArray.h>
#include <SxAtomicStructure.h>
#include <SxXC.h>
#include <SxTimer.h>
#include <SxProjector.h>
#include <SxPseudoPot.h>
#include <SxPW.h>

void inconsistencyExit (const SxString& what,
                        const SxString& where1,
                        const SxString& where2)
{
   cout << "Inconsistency!" << endl;
   cout << what << " differs between '";
   cout << where1 << "' and '" << where2 << "'.";
   cout << endl;
   SX_EXIT;
}

// TODO: write for metals
int main (int argc, char ** argv)  {

   initSPHInXMath (); // init math libraries

   // --- command line interface
   SxCLI cli (argc, argv);
   cli.preUsageMessage = 
      "Conversion program for GW program.\n"
      "Reads S/PHI/nX files and generates input for GW\n";
   cli.authors = "C. Freysoldt";

   SxString ibzWavesFile
      = cli.option ("--ibz","file","waves file for IBZ integration")
        .toString ();
   
   SxString bandWavesFile 
      = cli.option ("--band","file", "waves file for band structure")
        .toString ("");
   bool bandstructure = (bandWavesFile.getSize () != 0);
   
   double maxg2 = 0.;
   bool autoG2;
   autoG2 = ! cli.option("--ecut","determine energy cutoff").exists ();
   if (!autoG2) maxg2 = cli.last ().toDouble ();
   cli.last ().defaultValue = "default: enclosing sphere";
   
   int xcGroup = cli.newGroup ("read xc potential");
   SxString vxcFile 
      = cli.option ("--vxc","file","xc potential file").toString ();

   int rhoGroup = cli.newGroup ("calculate xc potential");
   cli.excludeGroup (xcGroup);
   SxString rhoFile
      = cli.option ("--rho","file","rho file for calculation of xc potential")
        .toString ();

   SxXC::XCFunctional xcFunctional = 
         (SxList<SxXC::XCFunctional> ()
          << SxXC::LDA 
          << SxXC::PBE 
          << SxXC::PBE_LDA)
         (cli.option("--lda|--pbe|--pbe_lda","xc functional to use")
             .required ().toChoice ());

   cli.setGroup (SxCLI::generalGroup);
   SxString outputFile
      = cli.option ("-o|--output","file","output file name")
        .toString ();
   bool progress = cli.option ("--progress","show conversion progress")
                   .toBool ();

   bool metal = cli.option ("-m|--metal","put occupation numbers into file")
                   .toBool ();

   cli.version ("$Id: Sx2gwst.cpp 4803 2007-10-11 20:31:43Z freysoldt.svn $");
   cli.finalize ();
   
   bool calculateVxc = cli.groupAvailable (rhoGroup);
   if (! (calculateVxc || cli.groupAvailable (xcGroup)) )  {
      cout << "Rho or xc potential must be provided." << endl;
      SX_QUIT;
   }

   // --- variables we'll need 
   SxString title = "GW data file", title2 = "S/PHI/nX generated";

   int ng;
   int nk = -1; // Monkhorst-Pack
   int nStates = -1; // Monkhorst-Pack
   int nk2 = -1; // band structure
   int nStates2 = -1; // band structure
   int nOccStates = -1; // nr of occupied states TODO: metals
   double nElec = 0.;

   SxMatrix<Double> kVecList, kVecList2;

   SxFermi fermi, fermi2;

   SxDiracVec<Complex16> vXcInG;

   SxVector3<Int> meshDim;
   SxVector<Int> n123IBZ, n123Band;
   SxVector<Int> nGk, nGk2;
   SxBinIO io;
   double eCut;

   SxAtomicStructure structure;

   // --- read input
   try {
      // --- ibz waves file
      io.open (ibzWavesFile, SxBinIO::BINARY_READ_ONLY);
      io.read ("meshDim", &meshDim);
      int n123Size = io.getDimension ("nAllGk");
      n123IBZ.resize (n123Size);
      io.read ("fftIdx", &n123IBZ, n123Size);
      io.read ("eCut", &eCut);
      int nSpin = io.getDimension ("nSpin");
      if (nSpin != 1)  {
         cout << "Can't do spin GW calculations!"  << endl;
         SX_QUIT;
      }
      nk = io.getDimension ("nk");
      nGk.resize (nk);
      io.read ("nGk", &nGk, nk);
      SxList<int> nPerK;
      io.read ("nPerK", &nPerK, nk);
      nStates = nPerK(0);
      for (int ik = 0; ik < nk; ik++)
         if (nPerK(ik) != nStates)  {
            cout << ibzWavesFile << ": ";
            cout << "All k-points must have same number of states!" << endl;
            SX_QUIT;
         }
      kVecList.reformat (nk, 3);
      io.read ("kVec", &kVecList, nk, 3);
      fermi.read (io);
      if (!metal)  {
         // --- check occupation numbers
         int tempNOcc;
         for (int ik = 0; ik < nk; ik++)  {
            tempNOcc = 0;
            for (int i = 0; i < nStates; i++)
               if (fabs (2. - fermi.focc(i, 0 /*iSpin*/, ik)) < 1e-3)
                  tempNOcc++;
            if (ik == 0) nOccStates = tempNOcc;
            if (nOccStates != tempNOcc)  {
               cout << "Number of occupied states differs at ik = " << ik << endl;
               cout << "Occupied states must have (2.0 - focc) < 1e-3." << endl;
               cout << "Lower ekt in case you have an insulator." << endl;
               SX_QUIT;
            }
         }
      } else {
         nElec = 0.;
         SxDiracVec<Double> weights(nk);
         io.read("kWeights", &weights,nk);
         for (int ik = 0; ik < nk; ++ik)
            for (int iSpin = 0; iSpin < nSpin; ++iSpin)
               nElec += weights(ik) * fermi.focc(iSpin,ik).sum ();
      }
      structure.read(io);
      io.close ();
      
      // --- band waves file ----------------------------------
      if (bandstructure)  {
         io.open (bandWavesFile, SxBinIO::BINARY_READ_ONLY);
         nSpin = io.getDimension ("nSpin");
         if (nSpin != 1)  {
            cout << "Can't do spin GW calculations!"  << endl;
            SX_EXIT;
         }
         SxVector3<Int> tempMeshDim;
         io.read ("meshDim", &tempMeshDim);
         for (int dim = 0; dim < 3; dim++)  {
            if (meshDim(dim) != tempMeshDim(dim))
               inconsistencyExit ("meshDim", ibzWavesFile, bandWavesFile);
         }
         double tempECut;
         io.read ("eCut", &tempECut);
         if (fabs (eCut - tempECut) > 1e-10)
            inconsistencyExit ("eCut", ibzWavesFile, bandWavesFile);
         SxMatrix3<Double> tempCell;
         io.read ("cell", &tempCell);
         int dim2;
         for (int dim = 0; dim < 3; dim++)
            for (dim2 = 0; dim < 3; dim++)
               if (structure.cell(dim, dim2) != tempCell(dim, dim2))
                  inconsistencyExit ("cell", ibzWavesFile, bandWavesFile);
         n123Size = io.getDimension ("nAllGk");
         n123Band.resize (n123Size);
         io.read ("fftIdx", &n123Band, n123Size);
         nk2 = io.getDimension ("nk");
         nGk2.resize (nk2);
         io.read ("nGk", &nGk2, nk2);
         io.read ("nPerK", &nPerK, nk2);
         nStates2 = nPerK(0);
         for (int ik = 0; ik < nk2; ik++)
            if (nPerK(ik) != nStates2)  {
               cout << bandWavesFile << ": ";
               cout << "All k-points must have same number of states!" << endl;
               SX_EXIT;
            }
         kVecList2.reformat (nk2, 3);
         io.read ("kVec", &kVecList2, nk2, 3);
         fermi2.read (io);
         io.close ();
      } // if bandstructure
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // set up rel->abs in 
   SxCell recCell (TWO_PI * structure.cell.inv.transpose ());

   // --- set up FFT
   SxFFT3d fft (SxFFT::Reverse, meshDim, structure.cell.volume);

   // --- get enclosing shell
   
   if (autoG2) {
      cout << "Look for enclosing sphere...";
      SxVector<Int>::Iterator n123It;
      SxArray<bool> isUsed (fft.meshSize);
      isUsed.set (false);
      
      // set true, if meshPosition is used
      for (n123It = n123IBZ.begin (); n123It != n123IBZ.end (); n123It++)
         isUsed(*n123It) = true;
      if (bandstructure)
         for (n123It = n123Band.begin (); n123It != n123Band.end (); n123It++)
            isUsed(*n123It) = true;

      // run through mesh and get largest G-Vector
      int x, y, z;
      SxVector3<Int> meshBound = meshDim / 2;
      double g2;
      for (x = meshBound(0) - meshDim(0); x < meshBound(0); x++)  {
         for (y = meshBound(1) - meshDim (1); y < meshBound(1); y++)  {
            for (z = meshBound(2) - meshDim (2); z < meshBound (2); z++)  {
               if (isUsed(fft.mesh.getMeshIdx(x,y,z,SxMesh3D::Origin)))  {
                  g2 = recCell.relToCar(SxVector3<Int> (x,y,z)).normSqr ();
                  if (g2 > maxg2) maxg2 = g2;
               }
            }
         }
      }

      maxg2 *= 1.0001; // gCut must be larger than largest G+k
   }
   cout << endl << "gCut = " << maxg2 << endl;

   // set up G-basis
   SxGBasis globalGBasis;
   globalGBasis.set (meshDim, structure.cell, maxg2);
   ng = globalGBasis.ng;

   // the G-vectors are expected in relative (integer) coordinates
   // --- construct relative coordinates from n123
   SxMatrix<Int> gVecRel (ng, 3);
   SxVector3<Int> G;
   int ig = 0;
   for (SxVector<TPrecFFTIdx>::Iterator n123It = globalGBasis.n123(0).begin ();
        n123It != globalGBasis.n123(0).end ();
        ++n123It, ig++)  {
      G = fft.mesh.getMeshVec (*n123It, SxMesh3D::Origin);
      gVecRel(ig,0) = G(0);
      gVecRel(ig,1) = G(1);
      gVecRel(ig,2) = G(2);
#ifndef NDEBUG
      double diff = (SxVector<Double>(recCell.relToCar(G))
              - toVector (globalGBasis.gVec.row(ig)) ).absSqr (). sum ();
      SX_CHECK (diff < 1e-10, ig, diff);
#endif
   }
   SX_CHECK (ig == ng, ig, ng);

   if (calculateVxc)  {
      // --- read rho  
      RhoR rhoR(1 /* nSpin */);
      if (progress)
         cout << "Reading rho file '" << rhoFile << "'." << endl;
      try {
         io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
         SxMatrix3<Double> tempCell;
         io.read ("cell",&tempCell);
         if ((tempCell - structure.cell).absSqr ().sum () > 1e-8)
            inconsistencyExit ("cell", rhoFile, ibzWavesFile);
         SxVector3<Int> tempMeshDim;
         io.read ("dim", &tempMeshDim);
         if ((meshDim - tempMeshDim).normSqr () != 0)
            inconsistencyExit ("dim/meshDim", rhoFile, ibzWavesFile);
         if (io.getDimension ("nMeshes") != 1)  {
            cout << rhoFile << " contains spin-polarized rho." << endl;
            SX_EXIT;
         }
         // read rho (readMesh gives SxArray<MeshR>)
         rhoR(0) = io.readMesh ()(0);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }

      // --- calculate vXC and Fourier transform it
      if (progress)
         cout << "Calculating vXC..." << endl;
      SxRBasis rBasis (meshDim, structure.cell);
      rBasis.registerGBasis (globalGBasis);
      globalGBasis.registerRBasis (rBasis);
      rhoR(0 /* iSpin */).setBasis (&rBasis);
      SxXC xc(/*nSpin*/1);
      xc.xcFunctional = xcFunctional;
      xc.computeXC ();
      xc.updateXC (rhoR);
      vXcInG = globalGBasis | xc.vXc(0);
   } else {
      // --- read vXC and Fourier transform it
      if (progress)
         cout << "Reading xc potential from '" << vxcFile << "'." << endl;
      try  {
         io.open (vxcFile, SxBinIO::BINARY_READ_ONLY);
         SxMatrix3<Double> tempCell;
         io.read ("cell",&tempCell);
         if ((tempCell - structure.cell).absSqr ().sum () > 1e-8L)
            inconsistencyExit ("cell", vxcFile, ibzWavesFile);
         SxVector3<Int> tempMeshDim;
         io.read ("dim", &tempMeshDim);
         if ((meshDim - tempMeshDim).absSqr ().sum () != 0)
            inconsistencyExit ("dim/meshDim", vxcFile, ibzWavesFile);
         if (io.getDimension ("nMeshes") != 1)  {
            cout << vxcFile << " contains spin-polarized xc potential." << endl;
            SX_EXIT;
         }
         SxMeshR vXC = io.readMesh ()(0 /* iSpin */);
         SxRBasis rBasis (meshDim, structure.cell);
         vXC.setBasis (&rBasis);
         vXcInG = globalGBasis | vXC;
         io.close ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   }
   // gw code has different prefactor convention
   vXcInG /= (fft.scaleRev * fft.meshSize);

   // --- put band energies in matrix (spin dimension removed)
   if (progress)
      cout << "Preparing band energies..." << endl;
   SxMatrix<Double> epsIBZ (nk, nStates);
   SxMatrix<Double> foccIBZ (nk, nStates);
   int iState;
   for (int ik = 0; ik < nk; ik ++)  {
      for (iState = 0; iState < nStates; iState++)  {
         epsIBZ (ik, iState) = fermi.eps  (iState, 0/*iSpin*/, ik);
         foccIBZ(ik, iState) = fermi.focc (iState, 0/*iSpin*/, ik);
      }
   }


   SxMatrix<Double> epsBS;
   if (bandstructure)  {
      epsBS.reformat (nk2, nStates2);
      for (int ik = 0; ik < nk2; ik ++)
         for (iState = 0; iState < nStates2; iState++)
            epsBS(ik, iState) = fermi2.eps (iState, 0/*iSpin*/, ik);
   }
   
   // --- the kVectors are expected in relative coordinates
   if (progress)
      cout << "Transforming k-vectors to relative coordinates ..." << endl;
   SxVector3<Double> kVecCart, kVecRel;
   int dim;
   for (int ik = 0; ik < nk; ik ++)  {
      for (dim = 0; dim < 3; dim++)
         kVecCart(dim) = kVecList(ik, dim);
      kVecRel = recCell.carToRel(kVecCart);
      for (dim = 0; dim < 3; dim++)
         kVecList(ik, dim) = kVecRel(dim);
   }
   if (bandstructure)  {
      for (int ik = 0; ik < nk2; ik ++)  {
         for (dim = 0; dim < 3; dim++)
            kVecCart(dim) = kVecList2(ik, dim);
         kVecRel = recCell.carToRel(kVecCart);
         for (dim = 0; dim < 3; dim++)
            kVecList2(ik, dim) = kVecRel(dim);
      }
   }
   
   // --- write output ------------------------------
   if (progress)
      cout << "Writing output '" << outputFile << "' ..." << endl;
   try  {
      SxBinIO output (outputFile, SxBinIO::BINARY_WRITE_LARGE);

      // --- WRITE_HEADER
      title.resize (80, true);
      title2.resize (80, true);

      output.write ("title", title);
      output.write ("title2", title2);

      output.addDimension ("xyz", 3);
      output.write ("bMat", SxMatrix3<Double> (), "xyz");

      int nSym = structure.cell.symGroupPtr->getNSymmorphic ();
      output.addDimension ("nSym", nSym);
      output.write ("symOps", structure.cell.symGroupPtr->getSymmorphic (),
                    "nSym", "xyz");
      
      output.addDimension ("ng", ng);
      output.write ("gVec", gVecRel, "ng", "xyz");
      
      output.addDimension ("nk-IBZ", nk);
      output.write ("kVec-IBZ", SxMatrix<Double> (), "nk-IBZ", "xyz");

      output.addDimension ("nStates-IBZ", nStates);
      if (metal)
         output.addDimension ("nElectrons", (int)lround(nElec));
      else
         output.addDimension ("nOccStates-IBZ", nOccStates);
      output.write ("eps-IBZ", epsIBZ, "nk-IBZ", "nStates-IBZ");
      output.write ("focc-IBZ", foccIBZ, "nk-IBZ", "nStates-IBZ");

      output.addDimension ("nCoeff-IBZ", nk * nStates * ng);
      output.write ("psi-IBZ", PsiG (), "nCoeff-IBZ");

      output.write ("vXC", vXcInG, "ng");

      if (bandstructure)  {
         output.addDimension ("nk-BS", nk2);
         output.write ("kVec-BS", SxMatrix<Double> (), "nk-BS", "xyz");

         output.addDimension ("nStates-BS", nStates2);
         output.write ("eps-BS", epsBS, "nk-BS", "nStates-BS");

         output.addDimension ("nCoeff-BS", nk2 * nStates2 * ng);
         output.write ("psi-BS", PsiG (), "nCoeff-BS");
      }
      
      // --- WRITE_DATA ---------------

      output.setMode (SxBinIO::WRITE_DATA);
      
      output.write ("title", title);
      output.write ("title2", title2);
      output.addDimension ("xyz", 3);

      // bMat is written as xyz:dim (dim counts basis vectors)
      output.write ("bMat", recCell, "xyz");

      output.addDimension ("nSym", nSym);
      {  
         SxArray<SymMat> symOp (nSym);
         SymMat E (1.,0.,0.,
                   0.,1.,0.,
                   0.,0.,1.);
         SxCell &cell = structure.cell;
         bool unityNotFirst = (cell.symGroupPtr->getSymmorphic(0) - E)
                              .absSqr ().sum () > 1e-10;
         for (int iSym = 0; iSym < nSym; iSym++)  {
            symOp(iSym) = cell.carToRel(cell.symGroupPtr->getSymmorphic(iSym));
            if (unityNotFirst && (symOp(iSym) - E).absSqr ().sum () < 1e-10)  {
               symOp(iSym) = symOp(0);
               symOp(0)    = E;
               cout << "Warning: exchanging symmetries 1 and " << (iSym + 1);
               cout << " in order to put identity in first place." << endl;
            }
         }
         output.write ("symOps", symOp, "nSym", "xyz");
      }

      output.addDimension ("ng", ng);

      output.write ("gVec", gVecRel, "ng", "xyz");

      output.addDimension ("nk-IBZ", nk);
      
      SX_CHECK (kVecList.nRows () == nk, kVecList.nRows (), nk);
      SX_CHECK (kVecList.nCols () == 3, kVecList.nCols ());
      output.write ("kVec-IBZ", kVecList, "nk-IBZ", "xyz");

      output.addDimension ("nStates-IBZ", nStates);
      if (metal)
         output.addDimension ("nElectrons", (int)round(nOccStates));
      else
         output.addDimension ("nOccStates-IBZ", nOccStates);

      output.write ("eps-IBZ", epsIBZ, "nk-IBZ", "nStates-IBZ");
      output.write ("focc-IBZ", foccIBZ, "nk-IBZ", "nStates-IBZ");

      output.addDimension ("nCoeff-IBZ", nk * nStates * ng);

      SxPtr<SxPW> waves;
      // --- read waves (ibz)
      if (progress)
         cout << "Reading IBZ waves '" << ibzWavesFile << "' ..." << endl;
      waves = SxPtr<SxPW>::create(ibzWavesFile,SxPW::ReadOnDemand);
      /*
      try {
         io.open (ibzWavesFile, SxBinIO::BINARY_READ_ONLY);
         waves.read (io);
         io.close ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
      */

      // --- sort and write waves (ibz)
      
      SxVector<Int> sortIdx (ng);
      SxVector<Int>::Iterator sortIdxIt;
      SxVector<TPrecFFTIdx>::Iterator globalN123It;
      SxVector<Int> globalMesh (fft.meshSize);

      PsiG psiGlobal (ng), psiI;
      PsiG::Iterator psiIt;
      
      int psiOffset = 0; // for writing
      
      // iterator runs over all k-points and all G-vectors
      SxVector<Int>::Iterator n123It = n123IBZ.begin (); 

      if (progress)
         cout << "Writing ibz waves ..." << endl;
      for (int ik = 0; ik < nk; ik++)  {
         if (progress)
            cout << "ik= " << (ik+1);
         // --- get sortIdx
         globalMesh.set (-1);
         for (ig = 0; ig < nGk(ik); ig++)
            globalMesh (*n123It++) = ig;
         sortIdxIt = sortIdx.begin ();
         globalN123It = globalGBasis.n123(0).begin ();
         for (ig = 0; ig < ng; ig++)
            *sortIdxIt++ = globalMesh(*globalN123It++);
         
         for (iState = 0; iState < nStates; iState++)  {
            if (progress && iState % (nStates / 20 + 1) == 0)
               (cout << '.').flush (); // one dot for about 5%
            psiI = (*waves) (iState, 0 /*iSpin*/, ik);
            // --- sort ... 
            psiIt = psiGlobal.begin ();
            sortIdxIt = sortIdx.begin ();
            for (ig = 0; ig < ng; ig++, sortIdxIt++)
               if (*sortIdxIt != -1)
                  *psiIt++ = psiI (*sortIdxIt);
               else
                  *psiIt++ = 0.;
            
            // --- ... and write waves
            output.write ("psi-IBZ", psiGlobal, "nCoeff-IBZ", psiOffset);
            psiOffset += ng;
         }
         if (progress) cout << endl;
      }
      SX_CHECK (psiOffset == nk * nStates * ng, psiOffset, nk * nStates * ng);

      output.write ("vXC", vXcInG, "ng");

      if (bandstructure)  {
         output.addDimension ("nk-BS", nk2);
         output.write ("kVec-BS", kVecList2, "nk-BS", "xyz");
         output.addDimension ("nStates-BS", nStates2);
         output.write ("eps-BS", epsBS, "nk-BS", "nStates-BS");
      
         output.addDimension ("nCoeff-BS", nk2 * nStates2 * ng);

         SX_CHECK (kVecList2.nRows () == nk2, kVecList2.nRows (), nk2);
         SX_CHECK (kVecList2.nCols () == 3, kVecList2.nCols ());

         // --- read waves (band)
         if (progress)
            cout << "Reading BS waves '" << bandWavesFile << "' ..." << endl;
         waves = SxPtr<SxPW>::create(bandWavesFile,SxPW::ReadOnDemand);
         /*
         try {
            io.open (bandWavesFile, SxBinIO::BINARY_READ_ONLY);
            waves.read (io);
            io.close ();
         } catch (SxException e) {
            e.print ();
            SX_EXIT;
         }
         */

         // --- sort & write waves (band)
         n123It = n123Band.begin ();
         psiOffset = 0;

         if (progress)
            cout << "Writing BS waves ..." << endl;
         for (int ik = 0; ik < nk2; ik++)  {
            if (progress)
               cout << "ik= " << (ik+1);
            // --- get sortIdx
            globalMesh.set (-1);
            for (ig = 0; ig < nGk2(ik); ig++)
               globalMesh (*n123It++) = ig;
            sortIdxIt = sortIdx.begin ();
            globalN123It = globalGBasis.n123(0).begin ();
            
            for (ig = 0; ig < ng; ig++)
               *sortIdxIt++ = globalMesh(*globalN123It++);
            for (iState = 0; iState < nStates2; iState++)  {

               if (progress && iState % (nStates2 / 20 + 1) == 0)
                  (cout << '.').flush (); // one dot for about 5%
               
               psiI = (*waves) (iState, 0 /*iSpin*/, ik);
               // --- sort ... 
               psiIt = psiGlobal.begin ();
               sortIdxIt = sortIdx.begin ();
               for (ig = 0; ig < ng; ig++, sortIdxIt++)
                  if (*sortIdxIt != -1)
                     *psiIt++ = psiI (*sortIdxIt);
                  else
                     *psiIt++ = 0.;

               // --- ... and write waves
               output.write ("psi-BS", psiGlobal, "nCoeff-BS", psiOffset);
               psiOffset += ng;
            }
            if (progress) cout << endl;
         }
         SX_CHECK (psiOffset == nk2 * nStates2 * ng,
                   psiOffset, nk2 * nStates2 * ng);
      }

      output.close ();
   } catch (SxException e)  {
      cout << endl;
      e.print ();
      SX_EXIT;
   }

}
