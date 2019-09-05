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
#include <SxLocpot.h>

enum FileType {
   None,
   sxb,
   VASP_LOCPOT
};

// Periodic average over w many points (w may be non-integer)
SxVector<Double> average (const SxVector<Double> &x, double w)
{
   int nAvg = int(floor(0.5 * w)), n = int(x.getSize ()), idx;
   double fw = w - (2*nAvg-1), sum;
   SxVector<Double> res(x.getSize ());
   for (int i = 0; i < n; ++i)  {
      // lower boundary
      idx = (i - nAvg) % n;
      if (idx < 0) idx += n;
      sum = 0.5 * x(idx);
      // upper boundary
      idx = (i + nAvg) % n;
      sum += 0.5 * x(idx);
      sum *= fw;

      // integer average
      for (int j = 1-nAvg; j < nAvg; ++j)  {
         idx = (i + j) % n;
         if (idx < 0) idx += n;
         sum += x(idx);
      }
      res(i) = sum / w;
   }
   return res;
}

/**
  @ingroup Add-Ons
  @author Alexey Dick, dick@fhi-berlin.mpg.de
**/

int main (int argc, char **argv)  
{
   // --- initialize S/PHI/nX utilities
   initSPHInXMath ();

   // --- parse input
   SxCLI cli(argc, argv);
   cli.preUsageMessage 
      = SxString(" Add-on averages mesh along one of the primitive cell axis")
        .wrap ();

   cli.authors = "Alexey Dick, Bjoern Lange";
  

   FileType fileType = sxb;
   if (cli.option ("--vaspLOCPOT", "potentials are VASP LOCPOT files").toBool ())  {
      fileType = VASP_LOCPOT;
   }

   SxString inFile
            = cli.option ("-i|--input","file","input mesh file")
             .toString ("rho.sxb");
   SxString weightFile
            = cli.option ("--weight","file","weight mesh file")
             .toString ("");
   double avgWidth = cli.option ("--average", "length", "local average (bohr)")
                     .toDouble (0., 0.);
   SxString outFile
            = cli.option ("-o|--output","file"," output file name")
             .toString   ("meshline.dat");
   if (cli.option ("-d|--dat",
                   "get output file name by substituting .sxb by .dat")
       .toBool ())  {
      if (!cli.error && inFile.contains (".sxb"))  {
         outFile = inFile.substitute (".sxb",".dat");
      } else {
         cout << "--dat option is not applicable:" << endl;
         cout << "Input file '" << inFile << "' has no .sxb" << endl;
         cli.setError ();
      }
   }
   SxString normalAxis
            = cli.option ("-a|--axis","string",
                          "one of primitive cell vectors along which averaging "
                          "is performed (a1, a2 or a3)")
             .toString   ("a1");
   int     iMesh
            = cli.option ("-m|--mesh","int","mesh number")
             .toInt   (1) - 1;
   

   // --- choose correct axis
   int nAxis;          // axis normal to plane
   int pAxis1, pAxis2; // axis parallel to plane
   nAxis = pAxis1 = pAxis2 = -1;

   if (normalAxis == "a1")  {
      nAxis  = 0;
      pAxis1 = 1;
      pAxis2 = 2;
   } else if (normalAxis == "a2")  {
      nAxis  = 1; 
      pAxis1 = 2;
      pAxis2 = 0;
   } else if (normalAxis == "a3")  {
      nAxis  = 2; 
      pAxis1 = 0;
      pAxis2 = 1;
   } else  {
      sxprintf ("Wrong normal axis\n");
      cli.setError ();
   }

   cli.finalize ();
   
   bool weighted = weightFile.getSize () > 0;
   // --- read in meshes
   SxMatrix3<Double> cell;
   SxMesh3D mesh;
   SxBinIO io;
   SxMeshR data, weight;
   RhoR meshesIn;

   if (fileType == sxb)  {
      try  {
         io.open (inFile, SxBinIO::BINARY_READ_ONLY);
         meshesIn = io.readMesh (&cell, &mesh);
         io.close ();
         if (weighted)  {
            io.open (weightFile, SxBinIO::BINARY_READ_ONLY);
            weight = io.readMesh ()(0);
            io.close ();
         }
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
   } else if (fileType == VASP_LOCPOT)  {
      meshesIn.resize(1);
      SxLocpot locpot(inFile);
      locpot.read ();
      cell = locpot.getCell ();
      mesh = locpot.getMesh ();
      meshesIn(0) = locpot.getPotential ();
   } 
      
   if (iMesh < meshesIn.getSize())  {  
      data = meshesIn(iMesh);
   } else {
      cout << "Invalid mesh number " << (iMesh+1) << endl;
      cout << '\'' << inFile << "' contains " << meshesIn.getSize ();
      cout << " meshes." << endl;
      SX_QUIT;
   }

   cout << " " << inFile << "(mesh:" << (iMesh+1) << ")";
   cout << " averaged along \"" << normalAxis << "\" axis"; 
   if (weighted)
      cout << " weighted with mesh 1 from '" << weightFile << "'";
   cout << " -> " << outFile << endl;
   

   // --- write plane to file (xfarbe format)
   double sum, wsum, w = 1.;
   SxVector3<Int> i; // fft index
   SxVector<Double> x(mesh(nAxis));
   SxVector<Double> y(mesh(nAxis));
   int idx;
   for (i(nAxis)=0; i(nAxis) < mesh(nAxis); i(nAxis)++)  {  
      sum = 0.;
      wsum = 0.;
      for (i(pAxis1)=0; i(pAxis1) < mesh(pAxis1); i(pAxis1)++)  {
         for (i(pAxis2)=0; i(pAxis2) < mesh(pAxis2); i(pAxis2)++)  {
            idx = (int)mesh.getMeshIdx(i, SxMesh3D::Positive);
            if (weighted) wsum += (w = weight(idx));
            sum += w * data(idx);

         }
      }
      if (weighted) sum /= wsum;
      x(i(nAxis)) = i(nAxis);
      y(i(nAxis)) = sum;
   }
   x /= double(mesh(nAxis)) / cell.col(nAxis).norm ();
   if (!weighted)
      y /= double(mesh(pAxis1) * mesh(pAxis2));

   if (avgWidth > 1e-16)  {
      double nAvg = avgWidth / cell.col(nAxis).norm () * mesh(nAxis);
      cout << "Averaging along line (" << nAvg << " points)" << endl;
      y = average (y, nAvg);
   }
   
   try  {
      io.open(outFile, SxBinIO::ASCII_WRITE_ONLY); 
      io.writeXYPlot (x, y); 
      io.close();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return 0;
}
