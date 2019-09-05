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

#include <SxRho.h>
#include <SxProjector.h>
#include <SxCLI.h>

int main (int argc, char **argv)  { 
   
   // --- parse command line arguments
   SxCLI cli(argc, argv);
   cli.preUsageMessage 
      = SxString("Create a saw-tooth-potential");
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   double eCut = cli.option ("--ecut", "number", "cutoff energy (Ry)")
                 .toDouble ();

   double height = cli.option ("--height", "number", "energy height (eV)")
                   .toDouble () / HA2EV;

   SxString outFile = cli.option ("-o|--output","file"," output file name")
                      .toString   ("sawpot.sxb");
   
   int dir = cli.option ("--dir","direction",
                         "potential sampling direction (1,2, or 3)")
             .toInt (3,1,3) - 1;

   double alpha = cli.option ("--alpha", "broadening", "apply Gaussian broadening").toDouble (0.);
   cli.last ().defaultValue = "default: no broadening";

   int asymGroup = cli.newGroup ("Asymmetric saw");
   double cutPos = cli.option ("--cut", "number", 
                               "cutting position (in bohr) along sampling axis")
                   .toDouble ();
   bool asym = cli.groupAvailable (asymGroup);

   cli.finalize ();

   SxCell cell;
   try {
      cell = SxCell(&*SxParser ().read (inFile));
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   SxMesh3D dim = SxGBasis::getMeshSize (eCut, cell);

   cout << dim << endl;

   SxVector3<Int> x;

   RhoR pot(1);
   pot(0).resize (dim.product ());

   double saw;
   for (x(0) = 0; x(0) < dim(0); ++x(0))  {
      for (x(1) = 0; x(1) < dim(1); ++x(1))  {
         for (x(2) = 0; x(2) < dim(2); ++x(2))  {
            if (asym)  {
               saw = x(dir) * height / dim(dir);
               if (x(dir) > cutPos/cell.basis(dir).norm () * dim(dir))
                  saw -= height;
            } else {
               saw = 2. * x(dir) * height / dim(dir);
               if (2 * x(dir) >= dim(dir)) saw = 2.*height-saw;
            }
            pot(0)(dim.getMeshIdx (x,SxMesh3D::Positive)) = saw;
         }
      }
   }
   if (alpha > 0.)  {
      SxFFT::quickFFTPlanner ();
      SxGBasis G;
      G.set (dim, cell, 2.*eCut);
      SxRBasis R(dim, cell);
      pot(0).setBasis (R);
      SxDiracVec<Double> gauss = exp(-0.5*alpha*alpha * G.g2);
      pot(0) = R | ( (G | pot(0)) * gauss);
   }
   
   // --- write down charge densities
   try {
      SxBinIO io;
      io.open      (outFile, SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (pot, cell, dim);
      io.setMode   (SxBinIO::WRITE_DATA);
      io.writeMesh (pot, cell, dim);
      io.close();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return 0;
}
