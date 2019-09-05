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
#include <SxAtomicOrbitals.h>
#include <SxPseudoPot.h>
#include <SxPAWPot.h>
#include <SxParser.h>

int main (int argc, char **argv)
{

   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "B.Lange";
   cli.preUsageMessage = 
      "This add-on calculates the OverlapMatrix for several Quamolfunctions\n";

   int radBasisPot = cli.newGroup ("fromPot");
   SxString inputFile 
      = cli.option("--input","file","specify sphinx input file")
        .toString ("input.sx");
   int radBasisFile = cli.newGroup ("fromSxbFile");
   cli.excludeGroup(radBasisPot);
   SxString quamolSxbFile 
      = cli.option("--quamolsxb","file","specify rad basis")
        .toString ("quamol.sxb");
   cli.setGroup (cli.generalGroup);
   SxString quamolFile 
      = cli.option("--quamol","file","specify functions")
        .toString ("quamol.sx");
   double limit = cli.option ("-l|--limit","limit value","limit value for basisfunctions").toDouble(1e-4);
   SxString outputFile = cli.option ("-o|--out", "file", "Quamol SXB file")
                                    .toString ("optQuamol.sxb");
   
   cli.finalize ();

   // Read in Quamols
   SxParser parser, quamolParser;
   SxConstPtr<SxSymbolTable> inputTable = parser.read (inputFile);
   SxConstPtr<SxSymbolTable> quamolTable = quamolParser.read (quamolFile);
   SxConstPtr<SxRadBasis> radBasisPtr;
   SxAtomicOrbitals functions;

   if (cli.groupAvailable(radBasisPot))  {
      if (inputTable->containsGroup("pseudoPot"))  {
         SxPtr<SxPseudoPot> psPotPtr = SxPtr<SxPseudoPot>::create (&*inputTable);
         radBasisPtr = SxConstPtr<SxRadBasis>::create(psPotPtr->rad, 
               psPotPtr->logDr);
      } else if (inputTable->containsGroup("pawPot"))  {
         SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*inputTable);
         pawPotPtr->extendRad(70.0);
         radBasisPtr = SxConstPtr<SxRadBasis>::create(pawPotPtr->rad, 
               pawPotPtr->logDr);
      } else   {
         cout << "No known Potential Group found !" << endl;
         SX_QUIT;
      }
   } else if (cli.groupAvailable(radBasisFile))  {
      try {
         SxBinIO io (quamolSxbFile, SxBinIO::BINARY_READ_ONLY);
         radBasisPtr = SxConstPtr<SxRadBasis>::create(io);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else  {
      cout << "No radial basis specified!" << endl;
      SX_QUIT;
   }

   try  {
      SxSymbolTable *quamolGroup= quamolTable->getGroup("Quamol");
      functions.setup(quamolGroup);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   functions.normalize ();

   // Calc Overlap matrices l dependent
   int nSpecies = functions.getNSpecies ();
   SxArray<SxList<SxDiracVec<Double> > > muSetList(nSpecies);

   for (int is = 0; is < nSpecies; is++)  {
      int lMax = functions.getLMax (is);
      for(int l = 0; l <= lMax; l++)  {
         SxMatrix<Double> S = functions.getOverlap(is,l);
         int nFL = int(S.nRows ());
         SxMatrix<Double>::Eigensystem eig;
         eig = S.eigensystem ();
         cout << SX_SEPARATOR;
         cout << "is = " << is 
              << ", l = " << l 
              << ", functions = " << nFL 
              << endl;
         eig.vals.print();
         int iot = 0;
         for (int i = 0; i < eig.vals.getSize(); i++)  {
            if (eig.vals(i).re > limit)  {
               int dim = int(functions.getFuncL(is,l,0).getSize());
               SxDiracVec<Double> newRadial (dim);
               newRadial.set(0.0);
               for (int ifl = 0; ifl < nFL; ifl++)  {
                 newRadial += eig.vecs(ifl,i) * functions.getFuncL(is,l,ifl);
               }
               newRadial.handle->auxData.is = is;
               newRadial.handle->auxData.ia = -1;
               newRadial.handle->auxData.n  = iot;
               newRadial.handle->auxData.l  = l;
               newRadial.handle->auxData.m  = NONE_M;
               newRadial.setBasis(*radBasisPtr);
               muSetList(is) << newRadial;
               iot++;
            }
         }
      }
   }
   cout << SX_SEPARATOR << endl;

   SxArray<SxArray<SxDiracVec<Double> > > muSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int iot = 0;
      int nOrbTypes = int(muSetList(is).getSize());
      muSet(is).resize(nOrbTypes);
      SxList<SxDiracVec<Double> >::Iterator it;
      for (it = muSetList(is).begin(); it != muSetList(is).end(); it++)  {
         muSet(is)(iot) = *it;
         iot++;
      }
   }

   // Write reduced functions
   SxAtomicOrbitals result (muSet,radBasisPtr);
   result.normalize ();
   for(int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = result.getNOrbTypes (is);
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         SxString file = "Quamol";
         file += is;
         file += iot;
         file += ".dat";
         SxBinIO out; 
         out.open(file, SxBinIO::ASCII_WRITE_ONLY);
         out.writeXYPlot(toVector(radBasisPtr->radFunc(is)),
                         toVector(result(is,iot)));
         out.close();
      }
   }

   result.write(outputFile);
}

