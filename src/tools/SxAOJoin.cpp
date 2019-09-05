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

#include <SxParser.h>
#include <SxRadBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxCLI.h>
#include <SxProjector.h>
#include <SxRadialAtom.h>

int main (int argc, char** argv)
{
   initSPHInXMath ();

   cout.precision(10);

   // Command line parsing
   SxCLI cli (argc,argv);

   // Define Author
   cli.authors = "B. Lange";

   // What does the program do ?
   cli.preUsageMessage = "AO Basis joining tool.";

   SxString basisFile = cli.option ("-b|--basis", "file", "AOBasis input file")
                                    .toString ("basis.sx");
   SxString outputFile = cli.option ("-o|--out", "file", "AOBasis sxb file")
                                    .toString ("basis.sxb");
   enum OrthoMode { DiagLaplace = 0, DiagR2 = 1 };
   int orthoMode = cli.option ("--laplace|--r2",
                               "diagonalize l-channel Laplacian/r^2 matrix")
                   .toChoice ();
   cli.finalize();

   SxConstPtr<SxRadBasis> radBasisPtr;
   SxParser parser;
   SxConstPtr<SxSymbolTable> table = parser.read (basisFile);
   SxAtomicOrbitals orbitals;
   
   try   {
         SxSymbolTable *basisGroup = table->getGroup("AOBasis");
         radBasisPtr = SxConstPtr<SxRadBasis>::create(basisGroup);
         orbitals.setup(basisGroup);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   
   if (orthoMode > -1)  {
      for (int is = 0; is < orbitals.getNSpecies (); ++is)  {
         SxDiracVec<Double> r2 = radBasisPtr->radFunc(is).sqr ();
         for (int l = 0; l < orbitals.getLMax (is); ++l)  {
            int n = orbitals.getFuncPerL(is,l);
            if (n <= 1) continue;
            SxMatrix<Double> S(n,n), R(n,n);
            SxDiracMat<Double> phi(r2.getSize (), n);
            for (int i = 0; i < n; ++i)  {
               phi.colRef (i) <<= orbitals.getFuncL(is,l,i);
               for (int j = 0; j < n; ++j)  {
                  S(j,i) = tr(orbitals.getFuncL(is,l,i) * orbitals.getFuncL(is,l,j));
                  if (orthoMode == DiagR2)
                     R(j,i) = tr(orbitals.getFuncL(is,l,i) 
                                 * orbitals.getFuncL(is,l,j) * r2);
                  else {
                     SX_CHECK (orthoMode == DiagLaplace);
                     SxDiracVec<Double> L;
                     L = SxRadialAtom::laplace (orbitals.getFuncL(is,l,j));
                     R(j,i) = -tr(orbitals.getFuncL(is,l,i) * L);
                  }
               }
            }
            cout << "is=" << is << "; l = " << l << "; n=" << n << endl;
            cout << "S=" << S << endl;
            SxMatrix<Double> L = S.choleskyDecomposition ().inverse ();
            cout << "R=" << R << endl;
            R = L.adjoint () ^ R ^ L;
            cout << "R(ortho)=" << R << endl;
            SxMatrix<Double>::Eigensystem eig = R.eigensystem ();
            eig.vecs = L ^ eig.vecs;
            cout << "eigenvals R: " << eig.vals << endl;
            phi = phi ^ toVector(eig.vecs);
            for (int i = 0; i < n; ++i)  {
               orbitals.getFuncL(is,l,i) <<= phi.colRef (i);
            }
         }
         for (int ipt = 0; ipt < orbitals.getNOrbTypes(is); ++ipt)  {
            SxString filename = "orbital-" + SxString(is) + "-" + SxString(ipt) + ".dat";
            FILE *fp = fopen (filename.ascii (), "w");
            if (fp)  {
               for (int ir = 0; ir < radBasisPtr->radFunc(is).getSize (); ++ir)  {
                  fprintf (fp, "%14.10f %14.10f\n", radBasisPtr->radFunc(is)(ir),
                           orbitals(is, ipt)(ir));
               }
               fclose (fp);
            }
         }
      }
   }

   orbitals.write(outputFile);
}
