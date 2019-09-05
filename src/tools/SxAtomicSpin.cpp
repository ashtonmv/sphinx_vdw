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
#include <SxBinIO.h>
#include <SxMesh3D.h>
#include <SxPAWRho.h>
#include <SxRBasis.h>
#include <SxRho.h>
#include <iostream>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   
   SxCLI cli(argc, argv);

   SxString rhoFile
            = cli.option("-r|--rho","file","input rho file")
            .toString   ("rho.sxb");

   SxString inputFile
            = cli.option("--input","file", "read potential from input file")
            .toString   ("input.sx");
   SxString outFile
            = cli.option("-o|--output", "file", "prefix of the output file;"
                         ".dat is always added to the end of the output file")
            .toString ("atomicSpin");

     cli.finalize ();


     SxParser parser;
     SxParser::Table table = parser.read(inputFile);

     SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);

     SxAtomicStructure structure = SxAtomicStructure (&*table);

     SxBinIO io;
     SxVector3<Int> dims;
     io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
     io.read("dim", &dims);

     SxCell cell;
     cell.read(io);
     sxprintf("Mesh: %i x %i x %i\n",dims(0), dims(1), dims(2));
     SxMesh3D mesh(dims);
     SxRBasis R(mesh, cell);
     ofstream out, out2;
     SxString ext = ".dat";
     SxString name = outFile + ext;
     out.open (name.ascii());

     // --- read rho file

     SxPAWRho pawRho(pawPotPtr);

     pawRho.pwRho = SxRho(io, &R);
     pawRho.readRho (rhoFile);
     
     SxPAWPot pawPot(&*table);


     // get atoms from structure and species
     int nAtom = structure.getNAtoms();

     for (int iAtom = 0; iAtom < nAtom; iAtom++) {
        int is = structure.getISpecies(iAtom);
        int ia = iAtom - structure.atomInfo->offset(is);
        int npt = pawPot.getNProjType (is);
        double spin = 0;
        // integration
        for (int ipt = 0; ipt < npt; ipt++) {
           int li = pawPot.lPhi(is)(ipt);
           int offsetI = pawPot.offset(is)(ipt) + li;
           for (int jpt = 0; jpt < npt; jpt++) {
              int lj = pawPot.lPhi(is)(jpt);
              int offsetJ = pawPot.offset(is)(jpt) + lj;
              if ( li == lj) {
                 double logdr = pawPot.logDr (is);

                 SxDiracVec<Double> rad = pawPot.rad(is);
                 double rPAW = pawPot.rc (is) * sqrt(-log(1e-4)); // calculation of radius of the PAW sphere
                 SxDiracVec<Double> rW (rad.getSize ());
                 for (int ir = 0; ir < rad.getSize (); ++ir)  {
                    if (rad(ir) < rPAW) rW(ir) = rad(ir) * rad(ir) * rad(ir);
                    else rW(ir) = 0.;
                 }

                 double integral = (rW * pawPot.phiAE(is).colRef(jpt) * pawPot.phiAE(is).colRef(ipt)).integrate (logdr);
                 for (int m=-li; m<=li; m++) {
                    spin += integral * (pawRho.Dij(0,is,ia)(offsetI+m,offsetJ+m)
                          -pawRho.Dij(1,is,ia)(offsetI+m,offsetJ+m));
                 }
              }
           }
        }
        cout << "The spin of Atom " << iAtom << " is " << spin << endl;
        // output
           out << "Atomnumber " << iAtom << "   atomspin " << spin << endl;
     }
     out.close();
}
