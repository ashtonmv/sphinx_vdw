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
#include <SxRedundantCoords.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);

   SxString inFile 
      = cli.option ("-i|--input", "input file", 
                    "take original input file")
        .toString ("input.sx");
   double maxDist = cli.option("-c", "cutoff", "max distance").toDouble (10.);
   double primarySetLimit = cli.option ("-d", "delta", "primary set limit")
                           .toDouble (0.05);
   double rmsThreshold    = cli.option ("-f", "rms factor", "join rms factor")
                            .toDouble (3.);
   double planeCutLimit   = cli.option ("-p", "polyhedron scale", "polyhedron scaling factor").toDouble (0.95);
   bool angles = cli.option ("-a", "with angles").toBool ();

   SxArray<int> bAtoms = cli.option ("-o", "atom list",
                                     "bond-orthogonal terms for these atoms")
                         .toIdxList ();
   cli.finalize ();

   // --- read input
   SxParser parser;
   SxConstPtr<SxSymbolTable> tree;
   tree = parser.read (inFile, "std/structure.std");
   SxAtomicStructure structure (&*tree);
   structure.readElements (&*tree);

   SxRedundantCoords ric;
   ric.verbose = true;
   ric.maxDist = maxDist;
   ric.primarySetLimit = primarySetLimit;
   ric.rmsThreshold = rmsThreshold;
   ric.planeCutLimit = planeCutLimit;
   ric.setup (structure);
   if (angles) ric.getAngles (structure);
   if (bAtoms.getSize () > 0) ric.getBornVonKarmanAngles (bAtoms);

   ric.param.resize (ric.getNParam ());
   SX_LOOP(i) ric.param(i) = 1./(1. + sqrt((double)(i+1)));

   int nDof = structure.getNAtoms () * 3;
   SxVector<Double> x(nDof), y, x2;
   x.set (0.);
   SxMatrix<Complex16> H(nDof, nDof);
   SX_LOOP(i)  {
      x(i) = -1.;
      H.colRef(i) <<= ric.applyH (structure, x);
      x(i) = 0.;
   }
   SxMatrix<Complex16>::Eigensystem eig = H.eigensystem ();
   cout << eig.vals << endl;
   SxVector<Double> complete(nDof);
   complete.set (0.);
   int nZero = 0;
   SX_LOOP(i)  {
      if (fabs (eig.vals(i).re < -1e-7))  {
         cout << "Warning: H has negative eigenvalue: " << eig.vals(i) << endl;
      }
      if (fabs (eig.vals(i).re) < 1e-7)  {
         nZero++;
         SxVector<Double> u = eig.vecs.colRef(i).getCopy ();
         for (int xyz = 0; xyz < 3; xyz++)  {
            double c = 0.;
            for (int j = xyz; j < nDof; j +=3) c += u(j);
            c /= (double)structure.getNAtoms ();
            for (int j = xyz; j < nDof; j +=3) u(j) -= c;
         }
         complete += u.absSqr ();
      }
   }
   //complete -= 3. / nDof;
   cout << complete << endl;
   cout << SX_SEPARATOR;
   cout << "Hessian will have " << nZero << " zero eigenvalues" << endl;
   SX_LOOP(iDof)  {
      if (complete(iDof) > 0.01)  {
         cout << "Incomplete configurational space: atom " << (iDof / 3 + 1)
              << " " << char((iDof % 3) + 'x') 
              << ":" << (1. - complete(iDof)) << endl;
      }
   }
   cout << endl;
   cout << "      ric {" << endl;
   cout << "         maxDist         = " << maxDist << endl;
   cout << "         typifyThreshold = " << primarySetLimit << endl;
   cout << "         rmsThreshold    = " << rmsThreshold << endl;
   cout << "         planeCutLimit   = " << planeCutLimit << endl;
   if (angles)
      cout << "         withAngles;" << endl;
   if (bAtoms.getSize () > 0)  {
      cout << "         bvkAtoms = [" << (bAtoms(0) + 1);
      for (int i = 1; i < bAtoms.getSize (); ++i)
         cout << "," << (bAtoms(i) + 1);
      cout << "];" << endl;
   }
   cout << "      }" << endl;

}

