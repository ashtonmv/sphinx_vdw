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

#include<SxCLI.h>
#include<SxFermi.h>
#include<SxParser.h>
#include<SxHamiltonian.h>

int main (int argc, char **argv)
{
   SxCLI cli(argc, argv);

   SxString inputFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                        .toString ("input.sx");
   int maxPerLine = cli.option ("-n", "number", 
                                "max. number of states per line").toInt (10,1);
   bool annotate = cli.option ("--annotate", "annotate occupations").toBool ();
   cli.finalize ();


   SxParser parser;
   SxParser::Table tree = parser.read (inputFile);
   const SxSymbolTable* table = &*tree;

   SxAtomicStructure str(table);
   double nElectrons 
      = SxHamiltonian::getNElectrons (table, 
                                      SxSpeciesData (table).valenceCharge,
                                      str);
                          
   int nSpin = SxHamiltonian::getNSpin (table);
   int nStates = SxHamiltonian::getNStates (table, nElectrons);
   SxKPoints kp(str.cell,table);
    
   SxFermi fermi(nElectrons, nStates, nSpin, kp);
   const SxSymbolTable *gGroup = table->getGroup("initialGuess");
   if (gGroup->containsGroup("occupations"))
      fermi.readOccupations (gGroup);
   else
      fermi.fermiDistribution (SxHamiltonian::getEkt(table));

   sxprintf("Note: filter with | sed -e'1,/occupations/d'\n");
   sxprintf ("   occupations  {\n      values = [");
   const char *nl =         "\n                ";
   int nk = fermi.getNk ();
   int n = nStates * nSpin * nk;
   for (int ik = 0; ik < nk; ++ik)  {
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         if (annotate && (ik > 0 || iSpin > 0))
            sxprintf("// ik = %d, iSpin = %d%s", ik+1, iSpin+1,nl);
         for (int i = 0; i < nStates; ++i)  {
            sxprintf (--n ? "%.6g," : "%.6g];\n   }\n", fermi.focc(i,iSpin,ik));
            if (i % maxPerLine % 10 == 4) sxprintf(" ");
            if ((i+1) % maxPerLine == 0 && n)  {
               if (annotate) sxprintf(" // <-- i = %d", i+1);
               sxprintf ("%s",nl);
            }
         }
         sxprintf ("%s",nl);
      }
   }
}
