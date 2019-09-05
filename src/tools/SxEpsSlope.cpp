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
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxPseudoPot.h>
#include <SxPerturbK.h>


int main (int argc, char **argv)
{

   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This add-on calculates d(eps)/dk.\n";
   /*
      We do this to test against numerical derivatives:
      sxepsslope ... --path | sed -e'1,/Start/d' | xmgrace -nxy -
      Data->transformations->Differences->Central differences for first set 

      Should match if k-direction is along x, y, or z.
      */

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   int n = cli.option ("-n","index","state").toInt (false,1) - 1;
   bool shortTangent = cli.option ("--shortTangent", "request short tangents")
                       .toBool ();
   bool longTangent  =  cli.option ("--longTangent", "request long tangents")
                       .toBool ();
   bool printTangent = shortTangent || longTangent;

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   bool nonLocal 
      = !cli.option ("--kin-only",
                     "omit non-local pseudopotential contributions").toBool ();

   int blockSize
      = cli.option ("--blocksize","set blocksize. Smaller blocksizes reduce "
                    "the memory consumption, but may slow down the "
                    " <phi|psi> part").toInt (64);

   bool path
      = cli.option ("--path", "first col in output is path length").toBool ();

   int precision
      = cli.option ("--prec", "digits", "set number of digits for output")
        .toInt (6);

   cli.version ("1.2");
   cli.finalize ();

   cout.precision (precision);


   SxPW waves;
   SxFermi fermi;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      waves.read (io);
      fermi.read (io);
      structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxGkBasis &gkBasis = waves.getGkBasis();
   fermi.kpPtr = &gkBasis;
   
   int ik, idir, jk;
   // int nStates = waves.getNStates ();
   int nk = waves.getNk ();
   int iSpin = 0;

   SxDiracMat<TPrecCoeffG> npn(nk,3);

   SxPerturbK kp;
   if (nonLocal)  {
      SxParser parser;
      SxParser::Table table = parser.read(inFile);
      SxPseudoPot psPot(&*table);
      gkBasis.changeTau(structure);
      kp.set (psPot, gkBasis, structure);
      kp.blockSize = blockSize;
   }
   
   SxDiracMat<TPrecCoeffG> kpMatElem;
   for (ik = 0; ik < nk; ik++)  {
      (cout << "ik = " << ik << endl).flush ();
      kpMatElem = kp.getMatrixElements (waves(n,iSpin,ik), 
                                        waves(n,iSpin,ik));
//                                        waves(iSpin,ik));
      for (idir = 0; idir < 3; idir++)  {
          npn(ik,idir) = kpMatElem.colRef(idir).chop ();
//         npn(ik,idir) = kpMatElem.colRef(idir)(n);
      }
      
   } // ik

   if (nonLocal)
      kp.printTimer ();
   
   // --- calculate lines E_nj(k) = <nj|p|nj> * (k - k_j) + eps(n,jk)
   SxArray<SxDiracVec<Double> > kVec(nk);
   for (ik = 0; ik < nk; ik++)
      kVec(ik) = SxDiracVec<Double> (gkBasis.kVec(ik));

   SxDiracVec<TPrecCoeffG> dk;
   double tangentVal;
   cout << "Start" << endl;
   double pathLength = 0.;
   for (ik = 0; ik < nk; ik++)  {
      if (path)  {
         cout << pathLength;
         if (ik + 1 < nk)
            pathLength += sqrt((kVec(ik+1) - kVec(ik)).absSqr ().sum ());
      } else
         cout << ik;
      cout << "\t" << fermi.eps(n,iSpin,ik) * HA2EV;
      if (printTangent)  {
         for (jk = 0; jk < nk; jk++)  {
            if (shortTangent && abs(ik-jk) > 1)  {
               cout << "\t0.0";
            } else {
               dk = kVec(ik) - kVec(jk);
               tangentVal = (npn.row(jk) ^ dk).chop ().re 
                            + fermi.eps(n,iSpin,jk); 
               cout << "\t" << (tangentVal * HA2EV);
            }
         }
      } else {
         for (idir = 0; idir < 3; idir++)
            cout << "\t" << (npn(ik, idir).re * HA2EV);
      }
      cout << endl;
   }
   
}

