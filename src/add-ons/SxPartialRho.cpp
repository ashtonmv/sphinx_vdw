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
#include <SxPartialRho.h>

#ifndef SX_STANDALONE
#include <SxProjector.h>
#include <SxConfig.h>

SxPartialRho::SxPartialRho (const SxPtr<SxPW>    &wavesIn, 
                            const SxPtr<SxFermi> &fermiIn)
: SxRho ()
{
   waves  = wavesIn;
   fermi  = fermiIn;
   rBasisPtr = &waves->getGkBasis ()(0).getRBasis ();
   useFocc   = true;
   symmetrize = false;
   SX_CHECK (fermi->getNSpin () == waves->getNSpin (),
             fermi->getNSpin (), waves->getNSpin ());
}


void SxPartialRho::compute(double eMin, double eMax) 
{
   int iSpin, nSpin = fermi->getNSpin ();
   int i, ik, nStates;
   int nk = fermi->getNk ();
   SX_CHECK (nk == waves->getNk (), nk, waves->getNk ());
   SX_CHECK (rBasisPtr);
   const SxRBasis &R = *rBasisPtr;
 
   rhoR.resize(nSpin);
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      rhoR(iSpin).resize(R.getMeshSize ());
      rhoR(iSpin).set (0.);
      rhoR(iSpin).setBasis (rBasisPtr);
   }
   
   PrecEps eig;
   PrecFocc focc;
   double weight;
   SxDiracVec<Double> psiSquare;
   for (ik=0; ik < nk; ik++)  {
      nStates = minimum (waves->getNStates (ik), fermi->focc.getNStates(ik));
      
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (i=0; i < nStates; i++)  {
            eig = fermi->eps(i,iSpin,ik);
            focc = useFocc ? fermi->focc(i,iSpin,ik)
                           : PrecFocc(3-nSpin);
            // --- include only energies from given interval 
            if ((eig >= eMin) && (eig <= eMax))  {
               weight = waves->getGkBasis ().weights(ik) * focc;
               if (weight > 0.00001)  {
                  psiSquare = (R | (*waves)(i,iSpin,ik)).absSqr();
                  // rhoR(iSpin) += weight * psiSquare;
                  rhoR(iSpin).plus_assign_ax (weight, psiSquare);
               }
            }
         }
      }
   }
   // --- symmetrize charge density
   if (symmetrize)  {
      cout << "Symmetrizing ..." << endl;
      for (iSpin=0; iSpin < nSpin; iSpin++)
         rhoR(iSpin) = R.symmetrize (rhoR(iSpin));
   }
}

void SxPartialRho::compute(const SxList<int> &kPoints,
                           const SxList<int> &states)
{
   SX_CHECK( kPoints.getSize() > 0, kPoints.getSize() );
   SX_CHECK( states.getSize() > 0, states.getSize() );
   
   int             iSpin, nSpin = waves->getNSpin();
   int             iStates, iKpoints; // for the lists 
   int             i,ik; // actual i, k from the list

   const SxRBasis &R = *rBasisPtr;

   rhoR.resize(nSpin);
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      rhoR(iSpin).resize(R.getMeshSize ());
      rhoR(iSpin).set (0.);
      rhoR(iSpin).setBasis (rBasisPtr);
   }
   
   
   PrecFocc occupation(3-nSpin);
   double weight;
   SxDiracVec<Double> psiSquare;
   for (iKpoints = 0; iKpoints < kPoints.getSize(); iKpoints++)  {
      ik = kPoints(iKpoints);
      if (ik >= waves->getNk () || ik < 0)  {
         cout << "Warning: Ignoring illegal k-point " << (ik+1);
         (cout << ".\n").flush ();
         continue;
      }
      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (iStates = 0; iStates < states.getSize(); iStates++)  {
            i = states(iStates);
            if (i >= waves->getNStates (ik) || i < 0)  {
               cout << "Warning: Ignoring illegal state " << (i+1);
               (cout << ".\n").flush ();
               continue;
            }
            // user means sorted index, not internal one:
            i = fermi->epsSortList(states(iStates), iSpin, ik);
            // set occupation if focc is used
            if (useFocc) occupation = fermi->focc(i,iSpin,ik);
            weight = waves->getGkBasis ().weights(ik) * occupation;
            if (weight > 0.00001)  {
               psiSquare = (R | (*waves)(i,iSpin,ik)).absSqr();
               // rhoR(iSpin) += weight * psiSquare;
               rhoR(iSpin).plus_assign_ax (weight, psiSquare);
            }
         }
      }
   }
   // --- symmetrize charge density
   if (symmetrize)  {
      cout << "Symmetrizing ..." << endl;
      for (iSpin = 0; iSpin < nSpin; iSpin++)
         rhoR(iSpin) = R.symmetrize (rhoR(iSpin));
   }
}

#else // SX_STANDALONE
#include <SxCLI.h>

// --------------------------------------------------------------------------
int main (int argc, char **argv)  
{    
   initSPHInXMath ();

   // --- parse command line arguments
   SxCLI cli(argc, argv);

   cli.preUsageMessage = SxString(
   "This add-on calculates partial electron densities.\n"
   "The states to be used for the calculation can be selected in two ways, "
   "either by energy or by state and kpoint index. The standard occupations "
   "(Fermi) can be overridden by the --full option.").wrap ();
   cli.authors = "A.Dick, C.Freysoldt";

   SxString inWaves 
      = cli.option ("-w|--waves","file","waves file")
        .toString ("waves.sxb");
   SxString outRho
      = cli.option ("-o","file","partial rho output file")
        .toString ("prho.sxb");
   
   SxFFT::plannerCLI (cli);
   // --- selection by energy
   int byEnergy = cli.newGroup ("Selection by energy");
   double eMin = cli.option ("--min", "energy", 
                             "minimum energy in eV")
                 .toDouble () / HA2EV ;
   double eMax = cli.option ("--max", "energy", 
                             "selection by energy: maximum energy in eV")
                 .toDouble () / HA2EV;

   bool ikStyle = !cli.groupAvailable (byEnergy);

   cli.newGroup ("Selection by index");
   cli.excludeGroup (byEnergy);
   SxList<int> states, kPoints;
   // --- selection by state number
   states = cli.option ("-n|--states","list",
                        "states (starting from 1) to be selected. "
                        "They are given as a ','-separated list of numbers or "
                        "ranges a-b, e.g. '-n 1,5-8,10' would select the states"
                        " 1, 5, 6, 7, 8, and 10.")
            .toIdxList ();
   cli.last ().defaultValue = "default: all states";
   kPoints = cli.option ("-k|--kpoints","list",
                         "k-points (starting from 1) to be "
                         "selected. The list is given like for the -n option")
             .toIdxList ();
   cli.last ().defaultValue = "default: all k-points";
   
   cli.setGroup (cli.generalGroup);
   bool useFocc 
      = !(cli.option("--full","selection by index: set occupation "
                     "to full occupation for all states (useful for unoccupied "
                     "states)")
        .toBool ());
   int symmetrize
      = cli.option ("--symmetrize|--nosymmetrize", 
                    "force (or suppress) symmetrization")
        .required(false)
        .toChoice ();
   enum { AUTO=-1, SYM=0, NOSYM=1};
   cli.last ().defaultValue = "default: symmetrize if all k-points are used";
   SxVector3<Int> newMesh(0,0,0); 
   int meshGroup = cli.newGroup("set mesh size explicitly");
   SxCLI::CliArg *opt = &cli.option("--mesh","mesh dimensions",
                                    "use these mesh dimensions");
   cli.last ().defaultValue = "default: use standard mesh from waves file";
   if (opt->exists ()) newMesh = SxVector3<Int> (opt->toIntList3 ());

   cli.newGroup("set mesh size via cutoff");
   cli.excludeGroup(meshGroup);
   double meshCut = cli.option("--ecut","energy",
                               "use mesh appropriate for this energy")
                    .toDouble (-1.);
   cli.last ().defaultValue = "default: use standard mesh from waves file";

   cli.finalize ();

   // --- read bases
   SxRBasis rBasis;
   SxPtr<SxFermi> fermi = SxPtr<SxFermi>::create ();
   SxCell cell;
   SxPtr<SxGkBasis> gkPtr = gkPtr.create ();
   try  {
      SxBinIO io (inWaves, SxBinIO::BINARY_READ_ONLY);
      fermi->read (io);
      SxVector3<Int> mesh;
      io.read("meshDim", &mesh);
      cell.read(io);
      rBasis.set (mesh, cell);
      gkPtr->read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
      
   SxPtr<SxPW> waves;
   waves = SxPtr<SxPW>::create (inWaves, SxPW::ReadOneByOne);
   waves->setGkBasisPtr (gkPtr);
   SxGkBasis &gk = *gkPtr;
   for (int ik = 0; ik < gk.getNk (); ++ik) gk(ik).registerRBasis (rBasis);

   // --- replace mesh if necessary
   if (meshCut >= 0.) newMesh = SxGBasis::getCommensurableMesh (meshCut, cell);
   if (newMesh.normSqr () > 0)  {
      cout << "Using mesh: " << newMesh << endl;
      
      // whether all k-points are used
      bool allK = kPoints.getSize () == 0 
               || kPoints.getSize () == waves->getNk ();
      // whether we need symmetry 
      bool needSymmetry = (symmetrize == SYM) || (symmetrize == AUTO && allK);
      // check mesh for symmetry
      if (   needSymmetry 
          && ! SxGBasis::isCommensurableMesh(newMesh, cell,true))
      {
         cout << "The supplied mesh is incompatible with symmetry.\n";
         SX_QUIT;
      }
      // --- replace FFTs
      SxFFT3d newFFT(SxFFT3d::Forward, newMesh, cell.volume);
      for (int ik = 0; ik < gk.getNk (); ++ik)
         gk(ik).replaceMesh (newFFT);
      rBasis.set (newMesh, cell);
   }

   SxPartialRho partialRho (waves,fermi);
   partialRho.useFocc = useFocc;

   // --- compute partial charge density
   if (!ikStyle)  {
      cout << "Selection by energy: " << endl;
      cout << "eMin = " << (eMin * HA2EV) << endl;
      cout << "eMax = " << (eMax * HA2EV) << endl;
      // default: symmetrize
      partialRho.symmetrize = (symmetrize != NOSYM);
      partialRho.compute (eMin, eMax);
   } else {
      // fill empty k-point and state lists (no list => ALL k-points / states)
      if (kPoints.getSize() == 0)
         for (int i=0; i<waves->getNk (); i++) kPoints << i;
      if (states.getSize() == 0)
         for (int i=0; i<waves->getNStates (); i++) states << i;
      
      cout << "Selection by index: " << endl << "k-points:";
      for (int i = 0; i < kPoints.getSize (); i++)
        cout << " " << (kPoints(i) + 1);
      cout << endl << "states:";
      for (int i = 0; i < states.getSize (); i++)
        cout << " " << (states(i) + 1);
      cout << endl;
      if (symmetrize == AUTO && kPoints.getSize () == waves->getNk () )
         symmetrize = SYM;
      partialRho.symmetrize = (symmetrize == SYM);
      partialRho.compute (kPoints, states);
   }
   // --- write partial charge density to file
   partialRho.writeRho (outRho);

}

#endif // SX_STANDALONE
