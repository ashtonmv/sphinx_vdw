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

int main(int argc, char **argv)
{
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   cli.authors = "C. Freysoldt";

   cli.preUsageMessage 
      = "This add-on splits mesh files into single-mesh files.";

   SxString meshIn 
      = cli.option ("-i|--input","file","multi-mesh file").toString ("rho.sxb");

   int explGroup = cli.newGroup ("explicit file names");

   SxList<SxString> meshOut 
      = cli.option ("-o|--output", "output files (multiple option)").toList ();
   cli.last ().defaultValue = "automatic file names";

   cli.newGroup ("automatic file names");
   cli.excludeGroup (explGroup);

   SxString stem
      = cli.option ("-s|--stem","name", 
                    "stem for output file names, the output files will be\n"
                    "<name>-<i>.sxb\n"
                    "for each mesh <i>.").toString ("rho");

   cli.version ("1.0");
   cli.finalize ();

   SxRBasis R;
   SxBinIO io;
   try  {
      io.open (meshIn, SxBinIO::BINARY_READ_ONLY);
      SxVector3<Int> mesh;
      io.read ("dim", &mesh);
      R.cell.read (io);
      R.fft3d.mesh = mesh;
      R.fft3d.meshSize = mesh.product ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxRho rho(io, &R);
   io.close ();

   int nMeshes = int(rho.rhoR.getSize ());

   SxString fileName; 
   for (int i = 0; i < nMeshes; i++)  {
      if (meshOut.getSize () == 0)  {
         fileName = stem + "-" + (i+1) + ".sxb";
      }  else  {
         if (i >= meshOut.getSize ())  {
            cout << "Can't write mesh " << (i+1) << endl;
            cout << "Too few filenames provided. Use --stem to autoname files";
            cout << endl;
            SX_QUIT;
         }
         fileName = meshOut(i);
      }
      (cout << "Writing file '" << fileName << "'...").flush ();
      SxRho(rho.rhoR(i)).writeRho (fileName);
      cout << endl;
   }
   return 0;
}

