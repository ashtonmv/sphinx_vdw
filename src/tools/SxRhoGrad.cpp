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
#include <SxFourierInterpol.h>
#include <SxRBasis.h>
#include <SxRho.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This compute a mesh gradient.";

   SxString outFile
      = cli.option ("-o","filename", "output file name")
        .toString ("gradient.sxb");
   
   SxFFT::plannerCLI (cli);

   SxString meshFile = 
      cli.option ("-i","mesh file","netcdf mesh file (one mesh only)")
      .toString ();

   cli.finalize ();

   // --- read input

   SxCell cell;
   SxVector3<Int> mesh;
   SxDiracVec<Double> realSpace;
   int meshSize;
   try {
      SxBinIO io;
      io.open (meshFile, SxBinIO::BINARY_READ_ONLY);
      if (! io.contains ("cell"))  {
         cout << endl << "'cell' is missing in netCDF file '";
         cout << io.filename << "'." << endl;
         cout << "This doesn't seem to be a netCDF mesh file." << endl;
         SX_QUIT;
      }
      cell.read (io);
      io.read ("dim", &mesh);

      meshSize = mesh.product ();
      realSpace.resize (meshSize);
      io.read ("mesh", &realSpace, meshSize);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   // --- Fourier transform to reciprocal space
   SxFFT3d fft(SxFFT::Both, mesh(0), mesh(1), mesh(2), cell.volume);
   SxDiracVec<SX_T_FFT_COMPLEX> in (realSpace), out(meshSize);
   fft.fftReverse (meshSize, in.elements, out.elements);
   SxDiracVec<Complex16> meshInG = out;

   SxDiracVec<TPrecCoeffG> gradMesh;
   gradMesh.reformat (meshSize, 3);

   // --- set up G-vectors and compute I G rho(G)
   SxCell recCell = cell.getReciprocalCell ();
   SX_LOOP(i)  {
      SxVector3<Double> g = recCell.relToCar (fft.mesh.getMeshVec ((int)i, SxMesh3D::Origin));
      SX_LOOP(iDir) gradMesh(i,iDir) = meshInG(i) * g(iDir) * I;
   }
   SxRBasis rBasis(mesh, cell);
   RhoR gradR(3);
   SX_LOOP(iDir)  {
      fft.fftForward (meshSize, gradMesh.colRef(iDir).elements, out.elements);
      gradR(iDir) = out;
      gradR(iDir).setBasis (&rBasis);
   }

   SxRho(gradR).writeRho (outFile);
}

