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
#include <SxPAWRho.h>
#include <SxCLI.h>
#include <SxGrid.h>

int main(int argc, char **argv)
{
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   cli.authors = "C. Freysoldt";

   cli.preUsageMessage 
      = "This add-on creates a supercell mesh file.";

   SxString meshIn 
      = cli.option ("-i|--input","file","file").toString ("rho.sxb");

   SxString outFile
      = cli.option ("-o|--output","file","file").toString ();

   // ---  definition
   SxMatrix3<Int> repMatrix;

   int simpleRep = cli.newGroup ("simple repetition");

   SxVector3<Int> diag (cli.option ("-r|--repeat","3 numbers",
                                     "repetition along current cell vectors")
                            .toIntList3 ("x,"));
   SxMesh3D finalMesh(0,0,0);
   if (cli.option ("-m|--mesh","3 numbers", "final mesh").exists ())  {
      finalMesh = SxVector3<Int>(cli.last ().toIntList3 ("x,"));
   }

   /*
   int trueRep = cli.newGroup ("3dimensional repetition");
   cli.excludeGroup (simpleRep);
   SxVector3<Int> col1 (cli.option ("--a1","vector",
                                     "new a1 in relative coordinates")
                         .toIntList3 ()),
                  col2 (cli.option ("--a2","vector",
                                     "new a2 in relative coordinates")
                         .toIntList3 ()),
                  col3 (cli.option ("--a3","vector",
                                     "new a3 in relative coordinates")
                         .toIntList3 ());
   */
   if (!cli.error)  {
      if (cli.groupAvailable (simpleRep))  {
         if (diag.product () <= 0)  {
            cout << "Illegal repetition factors " << diag << endl;
            cli.setError ();
         }
      /*
      } else if (cli.groupAvailable (trueRep))  {
         repMatrix = SxMatrix3<Int> (col1, col2, col3).transpose ();
         cout << repMatrix << endl;
         if (repMatrix.determinant () == 0)  {
            cout << "Illegal repetition matrix with zero determinant." << endl;
            cli.setError ();
         }
      */
      } else {
         cout << "no repetition given!" << endl;
         cli.setError ();
      }
   }
   cli.newGroup ("PAW");
   bool paw = cli.option ("--paw", "PAW density").toBool ();
   SxString structFile = cli.option ("--structure", "sx file",
                          "repeated structure file for PAW. Needed for atom sorting")
                         .toString ();

   cli.finalize ();

   SxRBasis R;
   SxBinIO io;
   SxMesh3D mesh;
   try  {
      io.open (meshIn, SxBinIO::BINARY_READ_ONLY);
      io.read ("dim", &mesh);
      R.cell.read (io);
      R.fft3d.mesh = mesh;
      R.fft3d.meshSize = mesh.product ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxRho rho(io, &R);

   int nMeshes = int(rho.rhoR.getSize ());

   // --- expand mesh data
   SxMesh3D newMesh = mesh * diag;
   for (int iMesh = 0; iMesh < nMeshes; ++iMesh)  {
      SxMeshR newRho (newMesh.product ());
      for (int x = 0; x < newMesh(0); ++x)  {
         for (int y = 0; y < newMesh(1); ++y)  {
            for (int z = 0; z < newMesh(2); ++z)  {
               newRho (newMesh.getMeshIdx (x, y, z, SxMesh3D::Positive))
                  = rho(iMesh)(mesh.getMeshIdx (x, y, z, SxMesh3D::Unknown));
            }
         }
      }
      rho(iMesh)=newRho;
   }

   R.cell = SxCell (R.cell.basis(0) * diag(0),
                    R.cell.basis(1) * diag(1),
                    R.cell.basis(2) * diag(2));
   R.fft3d.mesh = newMesh;
   R.fft3d.meshSize = newMesh.product ();
   if (finalMesh.product () > 0 && finalMesh != newMesh)  {
      SxFFT3d repFFT(SxFFT::Reverse, newMesh, R.cell.volume),
              finalFFT(SxFFT::Forward, finalMesh, R.cell.volume);
      int nMesh1 = newMesh.product (),
          nMesh2 = finalMesh.product ();
      SX_LOOP(iMesh)  {
         cout << "Interpolating mesh " << (iMesh+1) << "..." << endl;
         SxDiracVec<Complex16> rhoR = rho(iMesh);
         SxDiracVec<Complex16> rhoG1(nMesh1), rhoG2(nMesh2);
         repFFT.fftReverse (nMesh1, rhoR.elements, rhoG1.elements);
         rhoG2.set (0.);
         for (int x = -finalMesh(0)/2; 2 * x < finalMesh(0); x++)  {
            if (abs(x) * 2 >= newMesh(0)) continue;
            for (int y = 0; y < finalMesh(1); y++)  {
               if (abs(y) * 2 >= newMesh(1)) continue;
               for (int z = 0; z < finalMesh(2); z++)  {
                  if (abs(z) * 2 >= newMesh(2)) continue;
                  rhoG2(finalMesh.getMeshIdx(x,y,z,SxMesh3D::Origin))
                     = rhoG1(newMesh.getMeshIdx(x,y,z,SxMesh3D::Origin));
               }
            }
         }
         rhoR.resize (nMesh2);
         finalFFT.fftForward (nMesh2, rhoG2.elements, rhoR.elements);
         rho(iMesh) = rhoR;
      }
      R.fft3d.mesh = finalMesh;
      R.fft3d.meshSize = nMesh2;
   }

   if (paw)  {
      SxPAWRho pawRho;
      SxAtomicStructure strNew, strOrig;
      try {
         SxParser parser;
         strNew = SxAtomicStructure(&*parser.read (structFile, "std/structure.std"));

         pawRho.Dij.readRho (io);
         strOrig.read (io);
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
      pawRho.pwRho = rho;
      if ((strNew.cell - R.cell).absSqr ().sum () > 1e-10)  {
         cout << "Mismatch between cell in " << structFile << " and repeated rho:"
              << endl;
         cout << structFile << ": " << strNew << endl;
         cout << "repeated rho: " << R.cell << endl;
         SX_QUIT;
      }
      SxGBasis G;
      G.fft3d.resize (1);
      G.fft3d(0).mesh.set (0);
      G.structPtr = &strNew;
      R.registerGBasis (G);
      SxGrid grid(strOrig, 10);
      SxRadMat newDij;
      newDij.resize (strNew.atomInfo, nMeshes);
      SX_LOOP2(is,ia)  {
         int origIdx = strOrig.find (strNew.getAtom(is,ia), grid);
         if (origIdx < 0) {
            cout << "Atom not found: " << strNew.getAtom(is,ia) << endl;
            SX_QUIT;
         }
         int iaOrig, isOrig = strOrig.getISpecies (origIdx, &iaOrig);
         SX_LOOP(iSpin)
            newDij(iSpin, is,ia) = pawRho.Dij(iSpin, isOrig,iaOrig);
      }
      pawRho.Dij = newDij;
      pawRho.writeRho (outFile);
   } else {
      rho.writeRho (outFile);
   }
   io.close ();

   return 0;
}

