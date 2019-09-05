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
#include <SxConfig.h>
#include <SxBinIO.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxGrid.h>

SxMatrix3<Double> symmetrizeCell (const SxCell &cell)
{
   SxArray<SymMat> latSym = cell.getLatticeSymmetries ();
   int nSym = (int)latSym.getSize ();
   cout << "Cell has " << nSym << " lattice symmetries with epsSym = "
        << cell.epsSym << endl;
   SxMatrix<Double> U(9,6);
   U.set (0.);
   // free parameters: diagonal elements
   for (int i = 0; i < 3; ++i) U(i+3*i,i) = 1.;
   // free parameters: off-diagonal elements
   U(1+2*3,3) = U(2+1*3,3) = sqrt(0.5);
   U(0+2*3,4) = U(2+0*3,4) = sqrt(0.5);
   U(0+1*3,5) = U(1+0*3,5) = sqrt(0.5);
   SxMatrix<Double> M(9,9);
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      SymMat &S = latSym(iSym);
      //cout << (S.transpose () ^ S) << endl;
      SxMatrix3<Double> N = cell.inverse () ^ S ^ cell;
      //SX_LOOP2(i,j) cout << "delta N =" << N(i,j) - round(N(i,j)) << endl;
      SX_LOOP2(i,j) N(i,j) = round(N(i,j));
      S = cell ^ N ^ cell.inverse ();
      //cout << (S.transpose () ^ S) << endl;
      for (int j1 = 0; j1 < 3; ++j1)  {
         for (int k1 = 0; k1 < 3; ++k1)  {
            for (int j2 = 0; j2 < 3; ++j2)  {
               for (int k2 = 0; k2 < 3; ++k2)  {
                  M(j1 + 3 * k1, j2 + 3 * k2) = -S(j2, j1) * S(k2, k1);
               }
            }
            // delta term
            M(j1 + 3 * k1, j1 + 3 * k1) += 1.;
         }
      }
      SxMatrix<Double>::Eigensystem eig = SxMatrix<Double>(U.adjoint () ^ M ^ U).eigensystem ();
      eig.vecs = U ^ eig.vecs;
      int nDof = 0;
      for (int i = 0; i < eig.vecs.nCols (); ++i)
         if (fabs (eig.vals(i).re) < 1e-9)
            U.colRef (nDof++) <<= eig.vecs.colRef (i);
      if (nDof == 0) {
         cout << "Symmetry catastrophy: no degrees of freedom left..." << endl;
         SX_EXIT;
      }
      U.reshape (U.getSize ());
      U.resize (9 * nDof, true);
      U.reshape (9, nDof);
      //cout << "---" << endl;
   }
   cout << "Computing symmetry constraints from " << nSym << " symmetries: "
        << U.nCols () << " free parameters" << endl;
   //cout << M.eigenvalues ().real () << endl;
   SxVector<Double> b(9);
   b.set (0.); b(0) = b(4) = b(8) = 1.;
   b = U ^ (U.transpose () ^ U).inverse () ^ U.transpose () ^ b;
   b.reshape (3,3);
   SxMatrix<Double>::Eigensystem eig = SxMatrix<Double> (b).eigensystem ();
   b.set (0.); SX_LOOP(i) b(i,i) = sqrt(eig.vals(i));
   b = eig.vecs ^ b ^ eig.vecs.transpose ();
   // conserve volume
   b /= cbrt(SxMatrix3<Double>(b).determinant ());
   cout << "Transformation matrix: " << b << endl;
   return SxMatrix3<Double> (b) ^ cell;
}


int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage = "This symmetrizes an almost symmetric structure";
   cli.authors = "C. Freysoldt";

   SxString inFile
      = cli.option ("-i|--input", "input file", "S/PHI/nX input file")
        .toString ("input.sx");

   bool sxbFile = cli.option ("-b|--sxb", "input file is binary waves file "
                              "rather than S/PHI/nX input file")
                  .toBool ();

   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   bool outSxb = cli.option ("--outsxb", "output file is binary "
                              "rather than S/PHI/nX input file").toBool ();
   cli.last ().defaultValue = "default: same format as input file";
   if (sxbFile) outSxb = true;
   if (outSxb && outFile.getSize () == 0)
      outFile = "structure.sxb";

   double epsCell = -1;
   cli.option ("--eps_cell", "number", "epsilon for lattice symmetry"
               "detection. Increase this, if you are Ali ;-).")
      .required (false);
   if (cli.last ().exists ())
      epsCell = cli.last ().toDouble ();

   double epsStruct
      = cli.option ("-e|--eps", "number", "epsilon for structure "
                    "differences. Increase this, if symmetries are not "
                    "found.")
        .toDouble (1e-4);

   bool nonSymmorphic
      = cli.option ("--nonsymmorphic", "include non-symmorphic symmetries")
        .toBool ();

   bool reduceToPrimitive
      = cli.option ("--primitive","find primitive cell").toBool ();
   bool cellSym = cli.option ("--symmetrize-cell",
                              "symmetrize the unit cell vectors").toBool ();
   bool reorient
      = cli.option ("--reorient",
                    "rotate into standard orientation (if possible)").toBool ();

   int printOptions = SxAtomicStructure::DefaultPrint;
   if (cli.option("--printsym", "print symmetries").toBool ())
      printOptions |= SxAtomicStructure::PrintSymmetries;

   cli.authors = "C. Freysoldt";
   cli.finalize ();

   // --- read input
   SxAtomicStructure structure;
   if (sxbFile)  {
      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (inFile, "std/structure.std");
      structure = SxAtomicStructure (&*tree);
      structure.readElements (&*tree);
   }

   // --- get symmetries with reduced equality criterium
   if (epsCell > 0)
      structure.cell.epsSym = epsCell;
   if (cellSym)  {
      SxCell &cell = structure.cell;
      // change coordinates to relative coordinates
      SX_LOOP(ia) structure.ref(ia) = cell.carToRel (structure.ref(ia));
      // update the cell
      cell = symmetrizeCell (cell);
      cell.epsSym = EPS_SYM_DEFAULT;
      cell.setup ();
      // change coordinates back to cartesian
      SX_LOOP(ia) structure.ref(ia) = cell.relToCar (structure.ref(ia));
   }
   if (reorient)  {
      SxCell &cell = structure.cell;
      SxCell::CellType(cell).print ();
      SymMat U = cell.getStandardRot ();
      cell = U ^ cell;
      cell.setup ();
      structure = U ^ structure;
   }
   structure.epsEqual = epsStruct;
   if (nonSymmorphic)
      structure.updateSymGroup ();
   else
      structure.updateSymmetries ();
   //structure.print ();
   structure.cell.symGroupPtr->print ();

   // --- symmetrize structure
   int nSym = structure.cell.symGroupPtr->getSize ();
   SxAtomicStructure dVec, delta = structure.getNewStr ();
   delta.set (0., 0., 0.);
   SxGrid grid(structure, 10);
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      // rotate
      dVec = (*structure.cell.symGroupPtr)(iSym) ^ structure;
      // reorder
      dVec.replaceInfo (structure.match(grid,dVec));

      // determine structural difference
      dVec -= structure;
      dVec.map (SxCell::Origin);

      // average: sum differences...
      delta += dVec;
   }
   // ... and divide by total number
   delta /= (double)nSym;

   cout << "Center of mass shifted by " << (-delta.sum ()) << endl;

   // shift
   structure += delta;

   if (reduceToPrimitive)  {
      // --- change to primitive cell (after symmetrization!)
      SxCell primCell = structure.getPrimitiveCell ();
      structure = structure.getPrimStr (primCell);
   }


   // --- output
   if (outSxb)  {
      structure.write (outFile);
   } else {
      FILE *file;
      if (outFile.getSize () > 0)  {
         file = fopen(outFile.ascii(),"w");
         if (file == NULL)  {
            cout << "Can't open '" << outFile << "'." << endl;
            SX_EXIT;
         }
      } else {
         file = stdout;
      }

      structure.fprint (file, printOptions);

      if (file != stdout) fclose (file);

   }
}

