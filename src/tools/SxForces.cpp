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
#include <SxForceSym.h>
#include <SxStickyFilter.h>
#include <SxHamSolver.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxLoopMPI::init (argc, argv);

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates the forces from the waves and density";
   SxString rhoFile = 
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxString outFile =
      cli.option ("-o|--output","file","S/PHI/nX format output file")
      .toString ("");

   SxString wavesFile =
      cli.option ("-w|--waves","file","waves file")
      .toString ("waves.sxb");

   bool filter = cli.option("--filter","apply sticky filter").toBool ();

   SxFFT::plannerCLI (cli);

   cli.newGroup ("pulay");
   bool usePulay = cli.option ("--pulay",
                    "use Pulay Forces")
                   .toBool ();
   SxString basisFile = cli.option ("-b|--basis","file","waves file")
      .toString ("quamol.sxb");
   cli.version ("0.1");
   cli.finalize ();

   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   // initialize MPI level "waves-k"
   SxKPoints (SxAtomicStructure (&*table).cell, &*table);

   SxHamSolver potential (wavesFile, rhoFile, table);

   if (usePulay) {
      // read Basis
      SxBinIO ioBasis;
      ioBasis.open(basisFile, SxBinIO::BINARY_READ_ONLY);
      SxAtomicOrbitals TBOrbitals(ioBasis);
      ioBasis.close ();
      // construct SxMuPW
      SxConstPtr<SxMuPW> muBasisPtr = SxConstPtr<SxMuPW>::create (TBOrbitals, potential.wavesPtr->getGkBasisPtr ());
      //set  SxMuPWPtr in HamSolver
      potential.muBasisPtr = muBasisPtr;
   }
   
   const SxAtomicStructure &structure = potential.structure;

   // --- setup filter
   SxOperator<SxAtomicStructure> F;
   F.append (SxForceSym(structure));
   if (filter)  {
      cout << "Adding sticky filter ..." << endl;
      F.append(SxStickyFilter (table->getGroup("structure")) );
   }

   // compute and filter forces
   SxAtomicStructure fSym = F | potential.getForces ();

   // --- output
   cout << SX_SEPARATOR;
   cout << "| STRUCTURE:" << endl;
   cout << SX_SEPARATOR;
   structure.print (*potential.potPtr);
   cout << SX_SEPARATOR;
   cout << "| FORCES:" << endl;
   cout << SX_SEPARATOR;
   fSym.print (*potential.potPtr);
   
   // --- output in sx format
   if (outFile.getSize () > 0)  {
      FILE *fp = fopen (outFile.ascii (), "w");
      if (!fp)  {
         cout << "Cannot open " << outFile << " for writing " << endl;
         SX_QUIT;
      }
      structure.fprint (fp, fSym);
      fclose (fp);
   }
   SxTimer::getGlobalTimer().print ();

}

