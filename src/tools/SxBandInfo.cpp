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

#include <SxFermi.h>
#include <SxCLI.h>
#include <SxBinIO.h>

int main (int argc, char **argv)
{
   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on composes an xmgrace-compatible band structure + info file";
   SxString epsFile =
      cli.option ("--eps","file","read energies from this file")
      .toString ("eps.dat");
   SxString outFile = cli.option("-o|--output","file","output file")
      .toString ("bandinfo.dat");
   SxString infoFile = cli.option("-i|--info","file",
         "info input file: lines of <band index> <info>").toString ();
   bool sort = cli.option("-s|--sort","sort eigenenergies at each k-point")
               .toBool ();

   cli.finalize ();

   // --- read input


   int nk, nStates;
   SxFermi::peekSpectrumFile (epsFile, &nk, &nStates);
   SxKPoints dummyKP; // ugly: make better SxFermi constructor
   dummyKP.nk = nk;
   SxFermi fermi(1., nStates, 2, dummyKP);
   fermi.readSpectrum (epsFile);

   FILE *info = fopen(infoFile.ascii(), "r");
   if (info == NULL) {
      cout << "Can't open info file '" << infoFile << "'." << endl;
      SX_QUIT;
   }
   
   SxVector<Int> ik(nStates);
   ik.set (0);
   SxArray<SxList<double> > bandInfo;
   bandInfo.resize (nStates);

   int varRead, idx;
   double infoVal;
   int line = 0;
   while (!feof (info))  {
      line++;
      varRead = fscanf(info,"%d%lf", &idx, &infoVal);
      if (feof(info)) break;
      if (varRead != 2)  {
         cout << "line=" << line << endl;
         cout << "Warning: couldn't parse this part:" << endl;
         char stuff[81];
         fscanf(info,"%80s", stuff);
         cout << stuff;
         SX_QUIT;
      }
      idx--; // internal: start at 0 rather than 1
      if (idx >= nStates || idx < 0) {
         cout << "Ignoring unknown band " << idx << endl;
         continue;
      }
      if (ik(idx) < nk)  {
         bandInfo(idx).append(infoVal); 
         ik(idx)++;
      } else {
         cout << "Idx " << (idx) << " info exceeds number of k-points" << endl;
      }
   }
   fclose (info);

   FILE *output = fopen(outFile.ascii (), "w");
   if (output == NULL) {
      cout << "Can't open output file '" << outFile << "'." << endl;
      SX_QUIT;
   }
   for (int i = 0; i < nStates; i++)  {
      if (ik(i) == 0) continue;
      for (int k = 0; k < ik(i); ++k)  {
         // this will crash if number of k-points varies
         idx = sort ? fermi.epsSortList(i,0,k) : i;
         fprintf(output, "%d\t%f\t%f\n", k+1, fermi.eps(idx,0,k) * HA2EV,
                 bandInfo(idx)(k));
      }
      fprintf(output,"&\n");
   }
   fclose(output);
     
   return 0;

}

