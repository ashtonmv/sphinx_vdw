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

#include <SxWannierHam.h>
#include <SxString.h>
#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxLatticeShells.h>


#ifndef SX_STANDALONE


SxWannierHam::SxWannierHam ()
{
   // empty
}

SxWannierHam::~SxWannierHam ()
{
   // empty
}


#else /* SX_STANDALONE */


int main (int argc, char **argv)
{

   // --- parse command line
   SxCLI cli (argc, argv);

   cli.preUsageMessage =
      "Computes the Hamilton operator within the Wannier function basis set.";

   SxString wavesFile = cli.option
      ("-w|--waves", "file", "wavesfile").toString ("waves.sxb");

   SxString trafoFile = cli.option
      ("-t|--trafo", "file",
       "sxb-file containing the unitary trafo of the Bloch functions")
      .toString ("ukmn.sxb");

   SxString bsFile = cli.option  // "bs" for bandstructure
      ("-i|--input", "file",
       "input file for bandstructure, optional [if omitted, only the matrix "
       "elements <0n|H|Rm> are calculated]")
      .toString ("");

   int nShells = cli.option
      ("-s|--shells", "number",
       "how many shells of lattice vectors (levels of neighbourhood: "
       "1st nearest, 2nd nearest, and so on ...) shall be considered")
      .toInt (-1);

   // --- for compatibility with old S/PHI/nX only
   SxString prepFile = cli.option
      ("-p|--prep", "file",
       "file containing initial preparation of WF's (\"prep.sxb\") required "
       "[in case of using old S/PHI/nX output only]")
      .toString ("");

   bool fromOldSPHInX = cli.option
      ("--old", "use output from old S/PHI/nX version").toBool ();

   cli.authors = "Matthias Wahn";
   cli.version (1.0);
   cli.finalize ();

   // --- validate data
   if (fromOldSPHInX && prepFile == "")  {
      cout << SX_SEPARATOR;
      sxprintf ("| ERROR: When using output from old S/PHI/nX, you have to "
                "provide an sxb-file\n"
                "|        containing the initial preparation of the "
                "WF's (\"prep.sxb\"). Use\n"
                "|        option -p | --prep.\n");
      cout << SX_SEPARATOR;
      SX_QUIT;
   }

   if (!fromOldSPHInX && prepFile != "")  {
      cout << SX_SEPARATOR;
      sxprintf ("| ERROR: When providing a preparation file (\"prep.sxb\"), "
                "you have to use\n"
                "|        option --old as well.\n");
      cout << SX_SEPARATOR;
      SX_QUIT;
   }

   if (nShells < 1)  {
      cout << SX_SEPARATOR;
      sxprintf ("| ERROR: Option -s | -shells must specify a non-negative "
                "integer number.\n"
                "         You forgot to set a value or the value you have "
                "chosen is not valid.\n"
                "         Sorry.\n");
      cout << SX_SEPARATOR;
      SX_QUIT;
   }

   // --- read input file(s)
   SxGkBasis  gk;
   SxFermi    fermi;
   SxCell     cell;

   try  {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gk.read (io, false);  // true --> initialization of SxGBasis::FFT3d
      fermi.read (io);
      cell.read (io);
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   SxDirTensor3<TPrecCoeffG> D;

   if (!fromOldSPHInX)  {
      D.readDirTensor3 (trafoFile);  // D_mu,v^k := sum_m P_mu,m^k * U_m,n^k
   }  else  {
      SxDirTensor3<TPrecCoeffG> P;   // sum_v <Bloch_mu,k|Gauss_v> *
                                     //                      * (S^-1/2)_v,mu
      SxDirTensor3<TPrecCoeffG> U;   // the proper unitary trafo U

      P.readTensor3 (prepFile);
      U.readUkmn (trafoFile);

      // --- Ueberschiebung without summing
      D = P * U;
   }

//D.writeDirTensor3 ("D.sxb");
//D.print ("D - Tensor", 3);

   // --- validate data from input file(s)
   if (gk.getNk() != D.n1)  {
      cout << SX_SEPARATOR;
      sxprintf ("| ERROR: No consistency between waves file and trafo file:\n"
                "|        number of k-points differs. Sorry.\n");
      cout << SX_SEPARATOR;
      SX_QUIT;
   }

   if (D.n2 != D.n3)  {
      cout << SX_SEPARATOR;
      sxprintf ("| ERROR: The input file(s) do(es) not specify a set of "
                "quadratic\n"
                "         matrices. Sorry.\n");
      cout << SX_SEPARATOR;
      SX_QUIT;
   }

   // --- print eigenenergies
   int ik, i;
   int nk = D.n1;
   int nBands = D.n2;

   cout << "| eigenenergies [Hartree]:\n|\n";

   for (ik = 0; ik < nk; ik++)  {
      sxprintf ("|   ik = %3d :  ", ik);

      for (i = 0; i < nBands; i++)  {  // TODO: should go from iBottom to iTop
         sxprintf ("%13.6f ", fermi.eps (i,0,ik));
      }
      sxprintf ("\n");
   }

   // --- get nearest neighbour shells of lattice vectors
   SxLatticeShells shells (cell, nShells);
   shells.find ();

   // --- get matrix elements < 0n | H | Rm >
   cout << SX_SEPARATOR;
   cout << endl << SX_SEPARATOR;
   cout << "| Wannier representation of H -- < 0n | H | Rm >" << endl;

   SxDirTensor3<TPrecCoeffG> H;           // Hamiltonian in WF basis
   SxDiracMat<TPrecCoeffG>   HR;          // Hamiltonian at lattice vector R
   double                    d3kVBZ;      // V.E. in k-space times vol. of BZ
   Coord                     k, R;        // k-point, lattice vector
   int                       iR, nR;      // neighboured lattice vectors
   int                       i1, i2, mu;  // band indeces
   PrecCoeffG                elem;        // elem = <w1|H|w2>
   double                    phase;       // - k R
   SxComplex16               phaseFac;    // exp{- k R }

   H.reformat (nk, nBands, nBands);
   HR.reformat (nBands, nBands);
   d3kVBZ = 1. / nk;
   nR = shells.getNAllVecs ();

   for (iR = 0; iR < nR; iR++)  {
      R = Coord (shells.getNeighbVec (iR));
      HR.set (0.);

      for (i1 = 0; i1 < nBands; i1++)  {
         for (i2 = 0; i2 < nBands; i2++)  {

            // --- sum and integral
            elem.set (0.,0.);
            for (mu = 0; mu < nBands; mu++)  {
               for (ik = 0; ik < nk; ik++)  {
                  k         = gk.kVec(ik);
                  phase     = - (k * R).sum ();
                  phaseFac  = SxComplex16 ( cos(phase), sin(phase) );
                  elem     += D(ik,mu,i1).conj()
                            * D(ik,mu,i2)
                            * fermi.eps(mu,0,ik)
                            * phaseFac
                            * d3kVBZ;
               }  // :ik
            }  // :mu

            HR(i1,i2) = elem;
         }  // :i2
      }  // :i1

      H.set (HR, iR);  // = < 0n | H | Rm >

      cout << SX_SEPARATOR;
      cout << "< 0n | H | Rm >  with  R = " << R << endl << endl;
      SxWannierHam::printMatrix (HR);
   }  // :iR

   cout << SX_SEPARATOR;
   cout << "| Calculation of matrix elements finished." << endl;
   cout.flush ();
   cout << SX_SEPARATOR;
            
   // --- get bandstructure from < 0n | H Rm >, if wanted
   if (bsFile != "")  {

      // --- read k-points for bandstructure from input file
      SxParser parser;
      const SxParser::Table table = parser.read (bsFile);
      const SxAtomicStructure str (&*table);
      SxKPoints kPoints (str.cell, &*table);
      int nkBS = kPoints.getNk ();

      cout << endl;
      cout << SX_SEPARATOR;
      cout << "| k-points used for the bandstructure "
           << "[in Cartesian coordinates]:" << endl;
      cout << SX_SEPARATOR;
      cout << "|  -ik-     -x-      -y-       -z-      -label-" << endl;

      SxString x, y, z, l;
      for (ik = 0; ik < nkBS; ik++)  {
         k = kPoints.getK (ik);
         x = SxString (k(0), "%8.6f");
         y = SxString (k(1), "%8.6f");
         z = SxString (k(2), "%8.6f");
         if (   kPoints.kLabels.getSize() == nkBS
             && kPoints.kLabels(ik) != "")
            l = kPoints.kLabels(ik);
         else
            l = "";
         sxprintf ("| %4d: %9s %9s %9s    %s\n", 
                   ik+1, x.ascii(), y.ascii(), z.ascii(), l.ascii());
      }
   
      // --- calculate bandstructure
      cout << SX_SEPARATOR;
      cout << "| Bandstructure obtained from WF's [in eV]:" << endl;
      cout << SX_SEPARATOR;  cout.flush ();

      SxDirTensor3<TPrecCoeffG>             E(nkBS, nBands, nBands);
      SxDirTensor3<Double>                  bandstr(nkBS, 1, nBands);
      SxDiracMat<TPrecCoeffG>::Eigensystem  eigSysE;
      SxDiracMat<TPrecCoeffG>               Z(nBands, nBands);

      for (ik = 0; ik < nkBS; ik++)  {
         k = kPoints.getK (ik);

         Z.set (0.);
         for (iR = 0; iR < nR; iR++)  {
            R         = Coord (shells.getNeighbVec (iR));
            phase     = (k * R).sum ();
            phaseFac  = SxComplex16 ( cos(phase), sin(phase) );
            Z        += H(iR) * phaseFac;
         }
         E.set (Z, ik);

         eigSysE = E(ik).eigensystem ();
         bandstr.set (SxDiracMat<Double> (eigSysE.vals) * HA2EV, ik);
      }

      // --- dump bandstructure
      SxWannierHam::printEps (bandstr);
      SxWannierHam::printEps (bandstr, "epsWannier.dat");

      cout << SX_SEPARATOR;
   }

   sxprintf ("| Program exited normally.\n");
   cout << SX_SEPARATOR;

   return 0;
}

void SxWannierHam::printMatrix (const SxDiracMat<Complex16> &M)
{
   int i, nLines = int(M.nRows ());
   int j, nCols  = int(M.nCols ());

   for (i = 0; i < nLines; i++)  {
      for (j = 0; j < nCols; j++)  {
         sxprintf ("(%9.6f, %9.6f)", M(i,j).re, M(i,j).im);
         if (j < nCols-1)  {
            sxprintf (", ");
         }  else  {
            sxprintf ("\n");
         }
      }
   }
}

void SxWannierHam::printEps (const SxDirTensor3<Double> &bs,
                             const SxString &file)
{
   SX_CHECK (bs.n2 == 1, bs.n2);  // bs must contain effective line vectors

   FILE *fp = NULL;

   fp = fopen (file.ascii(), "w");
   if (!fp)  {
      sxprintf ("Can't open file \"%s\". Sorry.\n", file.ascii());
      SX_QUIT;
   }

   // --- print header lines
   fprintf (fp, "# Eigenspectrum obtained by Wannier functions.\n");
   fprintf (fp, "# All energies in eV.\n");

   // --- print data
   int ik, im;
   for (ik = 0; ik < bs.n1; ik++)  {
      fprintf (fp, "%d\t", ik+1);               // write k-point, starting at 1

      for (im = 0; im < bs.n3; im++)  {
         fprintf (fp, "%12.6f\t", bs(ik,0,im)); // write eigenenergy
      }
      fprintf (fp, "\n");
   }
   fclose (fp);
}

void SxWannierHam::printEps (const SxDirTensor3<Double> &bs)
{
   SX_CHECK (bs.n2 == 1, bs.n2);  // bs must contain effective line vectors

   sxprintf ("|  -ik-          -energies-\n");

   int ik, im;
   for (ik = 0; ik < bs.n1; ik++)  {
      sxprintf ("| %4d:\t", ik+1);            // write ik, starting at 0

      for (im = 0; im < bs.n3; im++)  {
         sxprintf ("%12.6f\t", bs(ik,0,im));  // write eigenenergy
      }
      sxprintf ("\n");
   }
}


#endif /* SX_STANDALONE */
