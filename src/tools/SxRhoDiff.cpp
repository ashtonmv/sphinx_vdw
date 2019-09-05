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

/**
  @ingroup Add-Ons
  @author  Alexey Dick, dick@fhi-berlin.mpg.de
**/

int main (int argc, char **argv)  { 
   
   // --- parse command line arguments
   SxCLI cli(argc, argv);
   cli.preUsageMessage 
      = SxString(" Add-on creates charge density file that contains"
                 " difference of two input charge densities").wrap ();
   SxString inFile1
      = cli.option ("--i1","file","input charge density file.")
             .toString ();
    SxString inFile2
      = cli.option ("--i2","file","input charge density file.")
             .toString ();
   SxString outFile
            = cli.option ("-o","file","output charge density file")
             .toString   ("difrho.sxb");
   cli.finalize ();
   cout << endl;
   
   // --- read in meshes
   SxMatrix3<Double> cell1, cell2;
   SxVector3<Int>    dim1, dim2;
   SxBinIO io;
   io.open (inFile1, SxBinIO::BINARY_READ_ONLY);
   RhoR rhoIn1 = io.readMesh (&cell1, &dim1);
   io.open (inFile2, SxBinIO::BINARY_READ_ONLY);
   RhoR rhoIn2 = io.readMesh (&cell2, &dim2);
   

   // --- perform checks
   if ( !(cell1 == cell2) )  { 
      cout << " Real space cells do not match:" << endl;
      cout << " " << inFile1 << " -> " << cell1 << endl;
      cout << "  vs. " << endl;
      cout << " " << inFile2 << " -> " << cell2 << endl;
      return 1;
   }
   if ( !(dim1 == dim2) )  { 
      cout << " FFT meshes do not match:" << endl;
       cout << " " << inFile1 << " -> " << dim1 << endl;
       cout << "  vs. " << endl;
       cout << " " << inFile2 << " -> " << dim2 << endl;
       return 1;
   }
   if (rhoIn1.getSize() != rhoIn2.getSize())  {
      cout << " Number of meshes do not match:" << endl;
      cout << " " << inFile1 << " -> " << rhoIn1.getSize() << endl;
      cout << "  vs. " << endl;
      cout << " " << inFile2 << " -> " << rhoIn2.getSize() << endl;
      return 1;
   }
    
   cout << " "  << inFile1 << " - " << inFile2;
   cout << " -> "<< outFile << endl;
   
   // --- write down charge densities
   RhoR rhoDiff = rhoIn1 - rhoIn2;
   io.open      (outFile, SxBinIO::BINARY_WRITE_ONLY);
   io.writeMesh (rhoDiff, cell1, dim1);
   io.setMode   (SxBinIO::WRITE_DATA);
   io.writeMesh (rhoDiff, cell1, dim1);
   io.close();

   return 0;
}
