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
#include <SxCLI.h>
#define a0 0.5291772083

/**
  @ingroup Add-Ons
  @author  Alexey Dick, dick@fhi-berlin.mpg.de
**/

int main (int argc, char **argv)  { 
   
   // --- parse command line arguments
   SxCLI cli(argc, argv);
   cli.preUsageMessage 
      = SxString(" Convert a chargedensity sxb file in"
                 " the format needed for openDX.").wrap ();
   SxString inFile
      = cli.option ("-i","file","input charge density file.")
             .toString ();
   SxString outFile
            = cli.option ("-o","file","output charge density file")
             .toString   ("rho.dx");
   SxVector3<Int> mult
              (cli.option ("-m","multiply","multiply charge density")
             .toIntList3 ("x,"));
   cli.finalize ();
   cout << endl;
   
   // --- read in meshes
   SxMatrix3<Double> cell;
   SxVector3<Int>    dim;
   SxBinIO io;
   io.open (inFile, SxBinIO::BINARY_READ_ONLY);
   RhoR rhoIn = io.readMesh (&cell, &dim);
   

   cout << " "  << inFile << " -> "<< outFile << endl;
  
   // TODO different Spin Channels
   if(rhoIn.getSize() > 1) {
      cout << "Different Spin channels not yet implemented!" << endl;
      SX_EXIT;
   } 
   // --- write down charge densities
   ofstream file;
   file.open(outFile.ascii ());
   int meshSize = dim(0)*dim(1)*dim(2);
   int meshSizeMult = mult(0) * mult(1) * mult(2) * meshSize;
   SxMesh3D mesh3D (dim);
   SX_CHECK(meshSize == rhoIn(0).getSize());
   file << "object 1 class array items " << meshSizeMult << " data follows" << endl;
   int counter = 0;
   for (int x = 0; x < mult(0)*dim(0); x++)  {
      for (int y = 0; y < mult(1)*dim(1); y++)  {
         for (int z = 0; z < mult(2)*dim(2); z++)  {
            int xCell = x % dim(0);
            int yCell = y % dim(1);
            int zCell = z % dim(2);
            ssize_t i = mesh3D.getMeshIdx(xCell,yCell,zCell,SxMesh3D::Positive);
            file << scientific << rhoIn(0)(i) << " ";
            if ((counter+1) % 8 == 0) file << endl;
            counter++;
         }
      }
   }
   file << endl;
   file << "attribute \"dep\" string \"positions\"" << endl << endl;
   file << "object 2 class gridpositions counts " 
        << mult(0)*dim(0) << " " << mult(1)*dim(1) << " " << mult(2)*dim(2) << endl;
   file << "origin " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
   file << "delta  " 
        << cell(0,0)/double(dim(0)) * a0 << " " 
        << cell(1,0)/double(dim(1)) * a0 << " " 
        << cell(2,0)/double(dim(2)) * a0 << endl;
   file << "delta  " 
        << cell(0,1)/double(dim(0)) * a0 << " " 
        << cell(1,1)/double(dim(1)) * a0 << " " 
        << cell(2,1)/double(dim(2)) * a0 << endl;
   file << "delta  " 
        << cell(0,2)/double(dim(0)) * a0 << " " 
        << cell(1,2)/double(dim(1)) * a0 << " " 
        << cell(2,2)/double(dim(2)) * a0 << endl << endl;

   file << "object 3 class gridconnections counts " 
        << mult(0)*dim(0) << " " << mult(1)*dim(1) << " " << mult(2)*dim(2) << endl;
   file << " attribute \"element type\" string \"cubes\"" << endl;
   file << " attribute \"ref\" string \"positions\"" << endl << endl;

   file << "object \"electron density\" class field" << endl;
   file << " component \"data\" 1" << endl;
   file << " component \"positions\" 2" << endl;
   file << " component \"connections\" 3" << endl;
   file.close();
   
   return 0;
}
