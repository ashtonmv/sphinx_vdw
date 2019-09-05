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
#include <SxFermi.h>
#include <SxHDF5.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "B. Lange";
   cli.preUsageMessage = "HDF5 to SXB Wave Convertor";

#ifdef USE_HDF5
   SxString wavesIn 
      = cli.option ("-i", "hdf5 waves file",
                     "wave function to read in")
        .toString ("waves.h5");
   SxString wavesOut
      = cli.option ("-o", "sx waves file",
                    "waves function to be written out")
        .toString ("waves.sxb");
   bool showReading
      = cli.option ("--show", "show Reading").toBool ();
   double nElectrons
      = cli.option ("--nElec", "Number of Electrons")
      .toDouble (-1.0);
   SxList<int> meshsize 
      = cli.option  ("--meshsize","XxYxZ meshsize")
      .toIntList ();
    bool ortho
      = cli.option ("--ortho", "show Reading").toBool ();

   cli.version ("0.1");
   cli.finalize ();


   SxHDF5 waves (wavesIn);

   // get eCut (Hartree -> Rydberg)
   double eCut 
      = 2.0 * waves.getData<double>("kinetic_energy_cutoff");
   // get gCut = 4.0 * eCut;
   double gCut = 4.0 * eCut;
   int nk = waves.getData<int>("number_of_kpoints");
   int nSpins = waves.getData<int>("number_of_spins");
   int nSpecies = waves.getData<int>("number_of_species");
   double TFile = waves.getData<double>("kinetic_energy");
   // double aLat = waves.getData<double>("lattice_constant");
   SxMatrix3<Double> cellMat 
      = waves.getMatrix3<Double>("primitive_vectors");
   SxVector3<Int> meshDims 
      = waves.getVector3<Int>("grid_size");
   SxMesh3D meshFile (meshDims);
   SxCell cell (cellMat);
   cell.setup ();
   
   waves.enterGroup("kpoints");
   SxList<double> weights; 
   SxVector<Double> weightsVec 
      = waves.getVector<Double>("weights",nk);
   if (fabs(weightsVec.sum ()-1.0) > 1e-10)  {
      cout << SX_SEPARATOR;
      cout << "WARNING: k-point weights do not sum to 1 but to "
           << weightsVec.sum () << "! Delta is " 
           << fabs(weightsVec.sum ()-1.0) << "!" << endl;
      cout << "Normalization ist forced!" << endl;
      cout << SX_SEPARATOR;
   }
   weightsVec /= weightsVec.sum ();
   SxList<Coord> kCoords; 
      for(int ik = 0; ik < nk; ik++)  {
         SxCell recCell = cell.getReciprocalCell ();
         SxString coordName = "reduced_coordinates";
         Coord relCoord = waves.getCoord(coordName,ik);
         Coord absCoord = (recCell.relToCar(relCoord));
         kCoords.append(absCoord);
         weights.append(weightsVec(ik));
      }
   waves.leaveGroup();

   waves.enterGroup("species");
   SxArray<SxString> chemNames(nSpecies);
   SxArray<SxArray<Coord> > coords (nSpecies);
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
     chemNames(iSpecies) = waves.getGroupNameByIdx (iSpecies);
     waves.enterGroupByIdx(iSpecies);
     int nAtoms = waves.getData<int>("number_of_atoms");
     for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
        SxString coordName = "reduced_atom_positions";
        Coord relCoord =  waves.getCoord(coordName,iAtom);
        Coord absCoord = (relCoord ^ cell);
        coords(iSpecies).append(absCoord);
     }
     waves.leaveGroup();
   }
   waves.leaveGroup();

   waves.enterGroup("waves");
   SxVector<Int> nStatesPerK (nk);
   SxArray<int> nCoeffs (nk);
   SxArray<SxArray<Coord> > redWaveCoords (nk);
   SxString nStatesName = "number_of_states";
   SxString nCoeffsName = "number_of_coefficients";
   SxString redCoordsPWName 
      = "reduced_coordinates_of_planewaves";
   for(int ik = 0; ik < nk; ik++)  {
     waves.enterGroupByIdx(ik);
     nStatesPerK(ik) = waves.getData<int>(nStatesName);
     nCoeffs(ik) = waves.getData<int>(nCoeffsName);
     SxMatrix<Int> waveCoords;
     waveCoords 
        = waves.getMatrix<Int>(redCoordsPWName,nCoeffs(ik),3);
     redWaveCoords(ik).resize(nCoeffs(ik));
     for (int iCoeff = 0; iCoeff < nCoeffs(ik); iCoeff++)
        redWaveCoords(ik)(iCoeff) 
           = Coord(waveCoords(iCoeff,0),
                   waveCoords(iCoeff,1),
                   waveCoords(iCoeff,2));
     waves.leaveGroup();
   }
   
   int nStates = nStatesPerK.minval ();

   Eps eigenvals (nStates,nSpins,nk);
   Focc foccs (nStates,nSpins,nk);
   Eps eKinFile (nStates,nSpins,nk);
   SxArray<double> eKinTotalFile (nk);

   for(int ik = 0; ik < nk; ik++)  {
      waves.enterGroupByIdx(ik);
      eKinTotalFile(ik) 
         = waves.getData<double>("total_kinetic_energy");
      for(int iSpin = 0; iSpin < nSpins; iSpin++)  {
         SxArray<hsize_t> lower (2), upper (2);
         lower(0) = 0; lower(1) = 0;
         upper(0) = 1; upper(1) = nStates;
         eigenvals(iSpin,ik) 
            = waves.getVBlock<Double>("eigenvalues",
                  lower, upper, nStates);
         foccs(iSpin,ik) 
            = waves.getVBlock<Double>("occupation",
                  lower, upper, nStates);
         eKinFile(iSpin,ik)
            = waves.getVBlock<Double>("kinetic_energy",
                  lower, upper, nStates);
      }
      waves.leaveGroup();
   }

   SxArray<SxArray<SxDiracMat<Complex16> > > coeffs (nk);
   SxString coeffName = "coefficients_of_wavefunctions";
   for(int ik = 0; ik < nk; ik++)  {
      waves.enterGroupByIdx(ik);
      coeffs(ik).resize(nSpins);
      for(int iSpin = 0; iSpin < nSpins; iSpin++)  {
         coeffs(ik)(iSpin).reformat(nCoeffs(ik), nStates);
         for(int iState = 0; iState < nStates; iState++)  {
            coeffs(ik)(iSpin).colRef(iState) 
               << waves.getWaveState(coeffName, iSpin,
                                     iState);
         }
      }
      waves.leaveGroup();
   }
   waves.leaveGroup();

   if (showReading)  {
      cout << "Result" << endl;
      cout << SX_SEPARATOR;
      cout << "Energy cutoff = " 
           << eCut << " Rydberg" << endl;
      cout << "Number of Spins = " << nSpins << endl;
      cout << "Number of k-points = " << nk << endl;
      cout << "Number of Species = " << nSpecies << endl;
      cout << "Number of States = " << nStates << endl;
      cout << "Cell = "; cellMat.print();
      
      cout << SX_SEPARATOR;
      
      cout << "k-points:" << endl;
      for (int ik = 0; ik < nk; ik++)  {
         cout << "\tNumber: " << ik+1 << endl;
         cout << "\tabs. Coords: ";
         kCoords(ik).print();
         cout << "\tweight: " << weights(ik) << endl;
      }

      cout << SX_SEPARATOR;
      
      cout << "Species:" << endl;
      for (int iSpecies = 0; 
           iSpecies < nSpecies; iSpecies++)  {
         cout << "\tSpecies name: " 
              << chemNames(iSpecies) << endl;
         size_t nAtoms = coords(iSpecies).getSize();
         cout << "\tNumber of Atoms = " << nAtoms << endl;
         cout << endl;
         for (size_t iAtom = 0; iAtom < nAtoms; iAtom++)  {
            cout << "\tabsolute atomic positions: "; 
            coords(iSpecies)(iAtom).print();
         }
         cout << endl;
         for (size_t iAtom = 0; iAtom < nAtoms; iAtom++)  {
            cout << "\trelative atomic positions: "; 
            (coords(iSpecies)(iAtom) 
             ^ cellMat.inverse()).print();
         }
      }

      cout << SX_SEPARATOR;

      cout << "Waves:" << endl;
      for(int ik = 0; ik < nk; ik++)  {
         cout << "\tk-point: "<< ik+1 << endl;
         cout << "\tNumber of states = " 
              << nStatesPerK(ik) << endl;
         cout << "\tNumber of coefficients = " 
              << nCoeffs(ik) << endl;
         cout << "\tEigenvalues :" << endl;
         for (int iSpin = 0; iSpin < nSpins; iSpin++)
            for (int iState = 0; iState < nStates; iState++)
               cout << eigenvals(iState,iSpin,ik) << endl;
         cout << "\tCoefficients :" << endl;
         for (int iSpin = 0; iSpin < nSpins; iSpin++)
            for (int iState = 0; iState < nStates; iState++)
               for (int iCoeff = 0; 
                    iCoeff < nCoeffs(ik); iCoeff++) {
                  cout << "\t" 
                     << iSpin << ","
                     << iState << ","
                     << iCoeff << ":"
                     << coeffs(ik)(iSpin)(iCoeff,iState) 
                     << " : " 
                     << "[" << redWaveCoords(ik)(iCoeff)(0)
                     << "," << redWaveCoords(ik)(iCoeff)(1) 
                     << "," << redWaveCoords(ik)(iCoeff)(2)
                     << "]" << endl;
               }
      }
      
      cout << SX_SEPARATOR;
   }

   SxAtomicStructure structure(cell);
   structure.startCreation ();
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      structure.newSpecies ();
      size_t nAtoms = coords(iSpecies).getSize ();
      for (size_t iAtom = 0; iAtom < nAtoms; iAtom++) 
            structure.addAtom (coords(iSpecies)(iAtom));
   }
   structure.endCreation ();
   cell.symGroupPtr 
      = SxPtr<SxSymGroup>::create(structure.getSymGroup ());

   SxVector3<Int> folding (1,1,1);
   SxKPoints kpoints(cell, kCoords, weights, folding);
   SxFermi fermi(nElectrons, nStates, nSpins, kpoints);
   for (int ik = 0; ik < nk; ik++)  {
      for (int iSpin = 0; iSpin < nSpins; iSpin++)  {
            fermi.eps(iSpin,ik) = eigenvals(iSpin,ik);
            fermi.focc(iSpin,ik) = foccs(iSpin,ik);
      }
   }

   cout << "MeshFile: "; meshFile.print();
   SxMesh3D mesh;
   if (meshsize.getSize() > 0)
      mesh = SxMesh3D (meshsize(0),meshsize(1),meshsize(2));
   else 
      mesh = SxGBasis::getCommensurableMesh(eCut, cell);

   cout << "predicted Mesh: "; mesh.print();
   SxGBasis gBasis (mesh, structure, gCut);
   SxPtr<SxGkBasis> gkBasisPtr 
      = SxPtr<SxGkBasis>::create(kpoints, gBasis, gCut);
   SxGkBasis &gkBasis = *gkBasisPtr;
   SxPW planeWaves(nStates,nSpins,gkBasisPtr);
   
   // Fill in waveCoeffs using the reduced (G+k) coordinates
   SxVector<Complex16> waveStore (mesh.getSize ());
   for (int ik = 0; ik < nk; ik++)   {
      size_t ngFile = redWaveCoords(ik).getSize();
      size_t ngWaves = gkBasis(ik).ng;
      for (int iSpin = 0; iSpin < nSpins; iSpin++)  { 
         for (int iState = 0; iState < nStates; iState++)  {
            planeWaves(iState,iSpin,ik).set(0.);
            waveStore.set(0.0);
            for (size_t ig = 0; ig < ngFile; ig++)  {
               int i1 = int(redWaveCoords(ik)(ig)(0));
               int i2 = int(redWaveCoords(ik)(ig)(1));
               int i3 = int(redWaveCoords(ik)(ig)(2));
               size_t iComp 
                  = mesh.getMeshIdx (i1,i2,i3,SxMesh3D::Unknown);
               waveStore(iComp) 
                  = coeffs(ik)(iSpin).colRef(iState)(ig);
            } //ig
            for (size_t ig = 0; ig < ngWaves; ig++)  {
               int iComp = gkBasis(ik).n123(0)(ig);
               planeWaves(iState,iSpin,ik)(ig) 
                  = waveStore(iComp);
            } //ig
         } //iState
      } //iSpin
   } //ik

   if (ortho) planeWaves.orthonormalize ();
   
   planeWaves.setGkBasisPtr(gkBasisPtr);
   planeWaves.writeWavesFile (wavesOut,fermi,structure,chemNames,true);

   // Test: Calculate <Psi|Psi>
   // Test: Calculate <G^2>

   SxComplex16 norm = 0.0;
   SxComplex16 T = 0.0;
   Eps TAll (nStates, nSpins, nk);
   Eps NAll (nStates, nSpins, nk);
   for (int ik = 0; ik < nk; ik++)  {
      for (int iSpin = 0; iSpin < nSpins; iSpin++)  {
         for (int iState = 0; iState < nStates; iState++)  {
            const SxGBasis &g = planeWaves.getGkBasis ()(ik);
            const PsiG &pw = planeWaves(iState,iSpin,ik);
            double normk = pw.absSqr ().sum ();
            double Tk = 0.5 * (pw.absSqr () * g.g2).sum ();
            double focc = fermi.focc(iState,iSpin,ik);
            double wk = gkBasis.weights(ik);
            TAll (iState,iSpin,ik) = Tk * focc;
            NAll (iState,iSpin,ik) = normk * focc;
            norm += wk * focc * normk;
            T += wk * focc * Tk;
         }
      }
   }


   if ((T-TFile).absSqr () > 1e-12)  {

      cout << "Norm is " << norm.re << endl;
      cout << "eKin is " << T.re << endl;
      cout << "eKinFile is " << TFile << endl;
      cout << "Delta is " << (T.re - TFile) << endl;

      for (int ik = 0; ik < nk; ik++)  {
         for (int iSpin = 0; iSpin < nSpins; iSpin++)  {
            for (int iState = 0; iState < nStates; iState++)  {
               double delta = TAll(iState,iSpin,ik)-eKinFile(iState,iSpin,ik);
               if ((delta*delta) > 1e-12)  {
                  cout << "\t" 
                     << ik << ","
                     << iSpin << ","
                     << iState << " : "
                     << NAll(iState,iSpin,ik) << " : "
                     << TAll(iState,iSpin,ik) 
                     << " should be " 
                     << eKinFile(iState,iSpin,ik) 
                     << ", Delta is " << delta
                     << endl;
               }
            }
            double delta = (TAll(iSpin,ik).sum () - eKinTotalFile(ik));
            if ((delta*delta) > 1e-12)  {
               cout << "Sum is : " 
                  << NAll(iSpin,ik).sum () << " : "
                  << TAll(iSpin,ik).sum () << " should be "
                  << eKinTotalFile(ik)
                  << ", Delta is " << delta
                  << endl;
            }
         }
      }
   }

#else
   
   cout << "Tool is only working with HDF5 support!" << endl;
   SX_QUIT;

#endif /* USE_HDF5 */   
}
