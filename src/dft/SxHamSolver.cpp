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


#undef INDIVIDUAL_K_TIMING

/**
 * Perl script to create files containing the individual timings per iteration
 * in a post processing step.  Analysis yields considerable load imbalance in k!
 *
#!/usr/bin/perl -w
use strict;
my $fi; my $fo;
my $fn="sphinx.log";
my $line; my $counter=0;
open ($fi, $fn) or die $!;
open ($fo, ">/dev/null") or die $!;
while ($line=<$fi>)
{
   chomp ($line);
   if ($line eq "BEGIN_K_TIMING")
   {
      close ($fo);
      open ($fo, ">k_it_$counter.dat") or die $!;
      next;
   }
   elsif ($line eq "END_K_TIMING")
   {
      close ($fo);
      open ($fo, ">/dev/null") or die $!;
      $counter++;
   }
   print $fo "$line\n";
}
 * --- end of the Perl script.
 */


#include <SxHamSolver.h>
#include <SxPW.h>
#include <SxMuPW.h>
#include <SxProjector.h>
#include <SxSymMatrix.h>
#include <SxRhoMixer.h>
#include <SxMatrix.h>
#include <SxConstants.h>
#include <SxFSAction.h>
#include <SxSpeciesRef.h>
#include <SxPWHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxPAWBasis.h>
#include <SxPAWSet.h>
#include <SxOverlap.h>
#include <SxParser.h>
#include <SxLoopMPI.h>  // LoopMPI
#include <SxFileParser.h>
#include <SxSimpleParser.h>
#include <SxInitPAW.h>
#include <SxInitPW.h>
#include <SxInitKP.h>
#include <SxFileIO.h>



#ifdef INDIVIDUAL_K_TIMING
#include <SxTime.h>
#endif



// ref1: Ab-initio Molekular Dynamik fuer fluessige Metalle, 
//       G. Kresse, 1993, TU Wien, Dissertation
// ref2: Numerical Receipes
// ref3: G. Kresse, J. Furthmueller, Phys. Rev. B 54, 11169 (1996)


SxHamSolver::SxHamSolver()
{
   // empty
}


SxHamSolver::SxHamSolver (const SxAtomicStructure &str, 
                          const SxSymbolTable *table)
   : energy (-1.),
     mesh (SxGBasis::getMesh(table, str.cell)),
     structure(str,SxAtomicStructure::Copy),
     G  (mesh, str, SxGBasis::getGCut(table)),       // |G>
     R  (G.fft3d(0).mesh, str.cell)                  // |R>
{
   nConvAbs = nCosAcc = nConvIt = nConvOcc = nConvRel = 0;
   nConvCycle = -1;
   nStepsAbs = nStepsCos = nStepsIt = nStepsOcc = nStepsRel = 0;
   dEnergyLow = false;

   cout << SX_SEPARATOR;
   cout << "| Energy cut-off : " << SxGBasis::getECut(table) << " Ry\n";
   cout << "| FFT mesh       : " << mesh << endl;
   cout << "| Size of G basis: " << G.ng << endl;
   cout << SX_SEPARATOR;
   {
      double cutVol = G.ng * str.cell.getReciprocalCell ().volume;
      double eCutEff = pow(cutVol * 3. / FOUR_PI, 2./3.);
      cout << "effective G basis cutoff: " << eCutEff << endl;
   }

   // --- clean up files produced by SxHamSolver
   SxBinIO::deleteFile ("energy.dat");
   SxBinIO::deleteFile ("residue.dat");
   SxBinIO::deleteFile ("spins.dat");

   printHartree = true;

   const SxSymbolTable *top = table->topLevel();
   potPtr = SxHamInit::setupPotential (top);
   structure.atomInfo->meta.update (SxAtomicStructure::Elements,
                                    potPtr->chemName);

   int nEmptyStates, nEmptyStatesChi;
   nSpin            = SxHamiltonian::getNSpin (top);
   nEmptyStates     = SxHamiltonian::getNEmptyStates (top);
   nEmptyStatesChi  = SxHamiltonian::getNEmptyStatesChi (top);
   ektDefault       = SxHamiltonian::getEkt (top);
   nElectrons       = SxHamiltonian::getNElectrons (top, potPtr->valenceCharge,
                                                    str);
   nStates          = SxHamiltonian::getNStates (top, nElectrons);
   nStatesChi = (nEmptyStatesChi >= 0)
              ?  nStates - nEmptyStates + nEmptyStatesChi : -1;
   nValStates = -1;

   // --- set number of components in G/R basis
   {
      int nComp = SxHamInit::getNComp (top);
      if (nComp > 1)  {
         G.setNComp (nComp);
         R.setNComp (nComp);
      }
   }

   // --- cross linking of dual bases R and G
   R.registerGBasis (G);
   G.registerRBasis (R);

   str.print (*potPtr);
}

SxHamSolver::SxHamSolver(const SxString &wavesFile, 
                         const SxString &rhoFile, 
                         SxConstPtr<SxSymbolTable> tablePtr, 
                         const SxString &tmpDir, 
                         bool saveMemory)
{
   // Read information from wavesFile
   try  {
      SxBinIO io;
      io.open (wavesFile, SxBinIO::BINARY_READ_ONLY);
      SxPtr<SxGkBasis> gkPtr = SxPtr<SxGkBasis>::create(io, true, saveMemory);
      fermi.read(io);
      fermi.kpPtr = gkPtr.getPtr ();
      nStates = fermi.getNStates ();
      nSpin = fermi.getNSpin ();
      wavesPtr = SxPtr<SxPW>::create(nStates,nSpin,gkPtr,tmpDir);
      wavesPtr->read(io,SxPW::KeepAll);
      structure.read(io);
      io.close();
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   init(rhoFile, tablePtr);
   setupHam(rhoFile, tablePtr);
}

SxHamSolver::SxHamSolver(const SxString &rhoFile, SxConstPtr<SxSymbolTable> tablePtr)
{
   init(rhoFile, tablePtr);
}

void SxHamSolver::init (const SxString &rhoFile,
                        SxConstPtr<SxSymbolTable> tablePtr)
{
   if (structure.getNAtoms() == 0) structure = SxAtomicStructure(&*tablePtr);

   // --- get plane-wave cutoff for big G basis
   double gCut;
   if (wavesPtr)  {
      SX_CHECK(wavesPtr->getGkBasisPtr ());
      gCut = SxGBasis::getGCut(wavesPtr->getGkBasisPtr()->gkCut);
   } else  {
      gCut = SxGBasis::getGCut(SxGBasis::getECut(&*tablePtr));
   }

   // Read information from rhoFile
   SxRho rho;
   try  {
      SxBinIO io;
      io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
      io.read("dim", &mesh);
      // --- setup G and R bases
      R.set(mesh, structure.cell);
      G.set(mesh, structure, gCut);
      // --- read in rho
      // TODO Why do I need RBasis ?
      rho = SxRho(io, &R);
      io.close();
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   // --- cross linking of dual bases R and G
   R.registerGBasis (G);
   G.registerRBasis (R);
   

   nConvAbs = nCosAcc = nConvIt = nConvOcc = nConvRel = 0;
   nConvCycle = -1;
   nStepsAbs = nStepsCos = nStepsIt = nStepsOcc = nStepsRel = 0;
   dEnergyLow = false;  

   printHartree = true;

   // --- set up potPtr
   if (   tablePtr->containsGroup ("PWHamiltonian") 
        || tablePtr->containsGroup ("Hamiltonian"))
   {
      potPtr = SxPtr<SxPseudoPot>::create (&*tablePtr);
   } else if (tablePtr->containsGroup ("PAWHamiltonian"))  {
      potPtr = SxPtr<SxPAWPot>::create (&*tablePtr);
   } else  {
      cout << "Constructor only for PP or PAW Hamiltonian!" << endl;
      SX_QUIT;
   }
   if (structure.getNSpecies () != potPtr->getNSpecies ())  {
      cout << "Inconsistency detected in species: " << endl;
      cout << "potential in input file has "
           << potPtr->getNSpecies () << " species." << endl;
      cout << "structure has               "
           << structure.getNSpecies () << " species." << endl;
      SX_QUIT;
   }
   structure.atomInfo->meta.update (SxAtomicStructure::Elements,
                                    potPtr->chemName);

   // --- get various dimensions
   if (fermi.getNSpin() > 0)  {
      nSpin            = fermi.getNSpin();
      ektDefault       = SxHamiltonian::getEkt (&*tablePtr);
      nElectrons       = fermi.nElectrons;
      nStates          = fermi.getNStates();;
      nValStates       = fermi.getNValenceBands();
   } else  {
      nSpin      = SxHamiltonian::getNSpin (&*tablePtr);
      ektDefault = SxHamiltonian::getEkt (&*tablePtr);
      nElectrons = SxHamiltonian::getNElectrons (&*tablePtr,
                                                 potPtr->valenceCharge,
                                                 structure);
      nStates    = SxHamiltonian::getNStates (&*tablePtr, nElectrons);

      int nEmptyStates     = SxHamiltonian::getNEmptyStates (&*tablePtr);
      int nEmptyStatesChi  = SxHamiltonian::getNEmptyStatesChi (&*tablePtr);
      nStatesChi = (nEmptyStatesChi >= 0)
                 ?  nStates - nEmptyStates + nEmptyStatesChi : -1;
      nValStates = -1;
   }
}


void SxHamSolver::setupHam(const SxString &rhoFile,
                           SxConstPtr<SxSymbolTable> tablePtr,
                           SxPtr<SxGkBasis> gkPtr)
{
   if (!gkPtr)  {
      SX_CHECK(wavesPtr);
      SX_CHECK(wavesPtr->getGkBasisPtr());
      gkPtr = wavesPtr->getGkBasisPtr();
   }

   SxGkBasis &Gk = *gkPtr;
   Gk.changeTau(structure);

   if (tablePtr->containsGroup ("PAWHamiltonian"))  {
      hamPtr = SxPtr<SxPAWHamiltonian>::create (G, potPtr, Gk, rhoFile,
                                                &*tablePtr);
   } else if (SxPWHamiltonian::getHamiltonianGroup (&*tablePtr)) {
      hamPtr = SxPtr<SxPWHamiltonian>::create (G, Gk, potPtr, rhoFile,
                                               &*tablePtr, wavesPtr);
   } else {
      cout << "Invalid Hamiltonian, must be PAW or PW" << endl;
      SX_QUIT;
   }
   structure.atomInfo->meta.update (SxAtomicStructure::Elements,
                                    potPtr->chemName);

}

SxHamSolver::~SxHamSolver ()
{
   // empty
}




bool SxHamSolver::isRegistered (const SxSymbolTable *cmd) const
{
   if (!cmd) return false;
   SxString str = cmd->getName();
   return (   str == "initialGuess" 
           || str == "SD"    
           || str == "scfDiag" 
           || str == "subspaceDiag" 
           || str == "CCG"    
           || str == "directDiag"
           || str == "bandStructure"
           || str == "hSqrCG"
           || str == "diffEqtn"
           || str == "linGrad"
           || str == "lcao"
           || str == "evalForces");
}


void SxHamSolver::execute (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (cmd);

   SxString str = cmd->getName();
   SxString indent = " "; // = cmd->level * ' '; 
   SxString prefix = "|" + indent;
   bool     storeWaves = true, storeRho = true;

   rhoMixing = foccMixing = deltaT = dPsiConv = dEnergy = dEps = dRelEps
             = ekt = dRelR = 0.;
   maxSteps  = 1;

   try  {
      storeWaves =  (cmd->contains("noWavesStorage"))
                 ? !(cmd->get("noWavesStorage")->toAttribute())
                 :   true;
      storeRho   =  (cmd->contains("noRhoStorage"))
                 ? !(cmd->get("noRhoStorage")->toAttribute())
                 :   true;
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   if      (str == "initialGuess")  initialGuess (cmd, calc);
   else if (str == "SD")            steepestDescent (cmd, calc);
   else if (str == "scfDiag")       scfDiagonalization (cmd, calc);
   else if (str == "subspaceDiag")  subspaceDiagonalization (cmd, calc);
   else if (str == "CCG")           allStateCG (cmd, calc);
   else if (str == "bandStructure") bandStructureCG (cmd, calc);
   else if (str == "hSqrCG")        hSqrCG (cmd, calc);
   else if (str == "diffEqtn")      diffEqtn (cmd, calc);
   else if (str == "directDiag")    directDiagonalization (cmd, calc);
   else if (str == "linGrad")       linGradTest (cmd, calc);
   else if (str == "lcao")          tightBindingPAW (cmd, calc);
   else if (str == "evalForces") {
      SxAtomicStructure force;
      const SxSymbolTable *hamGroup =
        SxHamiltonian::getHamiltonianGroup(cmd->topLevel());
      if (hamGroup->containsGroup("vdwCorrection")) {
         SxPWHamiltonian &H = *dynamic_cast<SxPWHamiltonian *> (hamPtr.getPtr());
         SxVDW VDW(structure, hamGroup, H.rho.rhoR);
         force = getSymForces (structure, NULL, VDW);
      } else force = getSymForces (structure, NULL);
      SX_MPI_MASTER_ONLY {
         SxString fname;
         SYMBOLPARSE(cmd) {
            fname = SYMBOLGET("file") || "forces.sx";
         }
         deleteFileOnce(fname);

         FILE * forceFile = sxfopen(fname, "a");
         structure.fprint (forceFile, force);
      }
   } else {
      SX_EXIT;
   }

   // --- write data to disc
   if (calc)  {
      writeData (storeWaves, storeRho);
      if (str != "initialGuess") printTiming (/* restart= */true,true); 
      else printTiming (/* restart= */true,false);
   }
}

void SxHamSolver::writeData (bool storeWaves, bool storeRho) const
{
   cout << SX_SEPARATOR;
   cout << "| Saving data\n";
   cout << SX_SEPARATOR;
   if (   dynamic_cast<SxPWHamiltonian *> (hamPtr.getPtr ())
       || dynamic_cast<SxPAWHamiltonian *> (hamPtr.getPtr ()))
       
   {
      if (storeRho)  {
         cout << "|   Charge density ...   "; cout.flush ();
         SX_MPI_MASTER_ONLY {
            hamPtr->getRho ().writeRho ("rho.sxb");
         }
         cout << "done\n";
      }  else  {
         cout << "storage omitted\n";
      }  
   }
   hamPtr->writePotentials ();

   cout << "|   Wavefunctions ...    "; cout.flush ();
   if (storeWaves)  {
#ifndef USE_PARALLEL_NETCDF4
      SX_MPI_MASTER_ONLY {
#endif
         write ("waves.sxb");
#ifndef USE_PARALLEL_NETCDF4
      }
#endif
      cout << "done\n";
   }  else  {
      cout << "storage omitted\n";
   }

   cout << "|   Eigenspectrum ...    "; cout.flush ();
   SX_MPI_MASTER_ONLY {
      /* not written to NetCDF hence no attempt for parallelization was made */
      fermi.writeSpectrum ("eps", "dat");
   }
   cout << "done\n";

   const SxPWHamiltonian *H 
      = dynamic_cast<const SxPWHamiltonian *>(hamPtr.getPtr ());
   if (H)  {
      cout << "|   Electrostatic potential ...    "; cout.flush ();
      SX_MPI_MASTER_ONLY {
         /* TODO   parallelize using NetCDF4 */
         H->writeVElStat ("vElStat-eV.sxb");
      }
      cout << "done\n";
   }

   cout << SX_SEPARATOR;
   cout.flush ();
}

void SxHamSolver::setupAOBasis (const SxString &TBInFile,
                                const SxGkBasis &gk)
{
   SX_CHECK (!aoBasisPtr);
   if (dynamic_cast<SxPAWHamiltonian *>(hamPtr.getPtr ()))   {
      // get PAW Hamiltonian
      SxPAWHamiltonian &hamPAW = dynamic_cast<SxPAWHamiltonian &>(*hamPtr);
      SxPAWPot &potPAW = *(hamPAW.potPtr);
      SxConstPtr<SxOverlapBase> SPtr = hamPAW.getS ();
      SxOverlap S(SPtr);
      SxAtomicOrbitals TBOrbitals;
      SxConstPtr<SxRadBasis> radBasisPtr;

      if (TBInFile.getSize () > 0 && SxFileInfo(TBInFile).exists ())  {
         cout << "Taking ao basis from file " << TBInFile << endl;
         radBasisPtr = SxConstPtr<SxRadBasis>::create(TBInFile);
         TBOrbitals.setBasis (radBasisPtr);
         TBOrbitals.read (TBInFile);
         TBOrbitals.print();
      } else {
         cout << "Taking ao basis from potential" << endl;
         cout << "Take the following PhiPS:" << endl;
         cout << "is\tnProjLocal\tlPhi" << endl;
         // get RadBasis
         radBasisPtr = potPAW.getBasisPtr ();
         // --- collect LCAO orbitals with initial occupation > 0
         SxArray<SxArray<SxDiracVec<Double> > > usedPhiPS(structure.getNSpecies());
         for (int is = 0; is < structure.getNSpecies(); is++)   {
            int nProjLocal = (int)potPAW.foccInit(is).getSize();
            SxList<int> usedProjectors;
            for (int ipl = 0; ipl < nProjLocal; ipl++)   
               if (potPAW.foccInit(is)(ipl) > 1e-10) usedProjectors.append(ipl);
            int nUsedProjectors = int(usedProjectors.getSize());
            usedPhiPS(is).resize(nUsedProjectors);
            for (int iupl = 0; iupl < nUsedProjectors; iupl++)   {
               int ipl = usedProjectors(iupl);
               usedPhiPS(is)(iupl) = potPAW.phiPS (is).colRef (ipl);
               usedPhiPS(is)(iupl).setBasis (radBasisPtr.getConstPtr());
               usedPhiPS(is)(iupl).handle->auxData.is = is;
               usedPhiPS(is)(iupl).handle->auxData.n = ipl;
               usedPhiPS(is)(iupl).handle->auxData.l = potPAW.lPhi(is)(ipl);
               usedPhiPS(is)(iupl).handle->auxData.m = 0;
               cout << is << "\t" << ipl << "\t" << potPAW.lPhi(is)(ipl) << endl;
            }
         }
         TBOrbitals = SxAtomicOrbitals(usedPhiPS, radBasisPtr);
         TBOrbitals.print();
      }

      // TB orbitals in G+k basis
      aoBasisPtr = SxPtr<SxAOBasis>::create (gk, TBOrbitals, SPtr);
   } else {
      cout << "Support for auxiliary AO basis currently only for PAW" << endl;
      SX_EXIT;
   }
}
         
void SxHamSolver::moveWaves (const SxAtomicStructure &newStr)
{
   if ((newStr - structure).absSqr ().sum () < 1e-12) return;
   if (aoBasisPtr)  {
      const SxAOBasis &aoBasis = *aoBasisPtr;
      SxOverlap S = hamPtr->getS ();
      cout << "Translating waves with ao basis set..." << endl;
      // orthonormalize waves
      SX_NEW_LOOP (*wavesPtr);
      for (int ik = 0; ik < wavesPtr->getGkBasis ().getNk (); ik++)  {
         SX_MPI_LEVEL("waves-k");
         if (SxLoopMPI::myWork(ik)) {
            SxGBasis &gk = wavesPtr->getGkBasis ()(ik);
            const SxBasis &waveBasis = wavesPtr->getBasis (ik);

            for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
               gk.changeTau (structure);
               SxDiracMat<Complex16> aoS = aoBasis.calculateOverlap (ik);
               SxDiracMat<Complex16> invS = aoS.inverse ();

               PsiGI psiG = gk | (*wavesPtr)(iSpin, ik);
               PsiGI psiAO = invS ^ (aoBasis | psiG);
               psiAO.setBasis (aoBasis);

               psiG -= gk | psiAO;

               cout << "Avg. square residue at ik=" << (ik+1) 
                  << ": " << ((psiG.conj () * (S * psiG)).sum ().re / (double)psiG.nCols ())
                  << endl;

               //cout << "Cross-check: " << (aoBasis | psiG).normSqr () << endl;

               gk.changeTau (newStr);

               psiG += gk | psiAO;

               //cout << "Before ortho: " << (psiG.adjoint () ^ ( S * psiG)) << endl;

               // --- recompute projections
               (*wavesPtr)(iSpin, ik) <<= waveBasis | psiG;
               S.orthonormalize (&(*wavesPtr)(iSpin,ik), Loewdin);
            }
         }
      }
   }
   structure <<= newStr;
}

SxAtomicStructure SxHamSolver::getForces (const SxAtomicStructure &x, 
                                          const SxSymbolTable *cmd)
{
   bool Hirshfeld = false;
   if (cmd && cmd->contains("Hirshfeld"))  {
      Hirshfeld = cmd->get("Hirshfeld")->toAttribute ();
   }
   if (SxLoopMPI::nr () > 1)  {
      // make 100% sure that all tasks have the same structure
      SxLoopMPI::bcast (const_cast<SxAtomicStructure&>(x).coords, 0);
   }
   bool tauChanged = false;
   if (x.getSize () > 0)  {
      tauChanged = (x - structure).absSqr ().maxval () > 1e-20;
   }
   if (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()))   {
      SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr); 
      if (tauChanged)  {
         // --- get shape of ideal atomic densities in G space
         //     and overlap of all those densities
         const SxPseudoPot& pot = *SxPtr<SxPseudoPot>(potPtr);
         SxArray<PsiG> atomRho(potPtr->getNSpecies ());
         PsiG allAtom(G);
         allAtom.set (0.);
         SX_LOOP(is)  {
            atomRho(is) = pot.getAtomRhoG (G, (int)is);
            allAtom += atomRho(is) * G.structureFactors(is);
         }
         SxDiracVec<Double> allAtomR = R.symmetrize (R | allAtom);
         // switch to the deformation density: subtract ideal atoms
         SX_LOOP(iSpin)
            H.rho(iSpin).plus_assign_ax(-1./nSpin, allAtomR);

         if (Hirshfeld) H.rho.displaceHirshfeld (x, atomRho, &allAtomR);

         // --- now change the structure
         structure.copy (x);
         H.structure = structure;
         G.changeTau (structure);
         wavesPtr->getGkBasis().changeTau (structure);

         // switch back from deformation density to full density
         // by adding the atomic contributions
         allAtom.set (0.);
         SX_LOOP(is) allAtom += atomRho(is) * G.structureFactors(is);
         allAtomR = R | allAtom;
         SX_LOOP(iSpin)
            H.rho(iSpin).plus_assign_ax(+1./nSpin, allAtomR);

         H.update (fermi); // recompute Hamiltonian
      }
      H.calcForces = true;    // TODO: ugly!!!
      H.compute (fermi, true, true);
      H.calcForces = false;   // TODO: ugly!!!

      if (cmd)  {
         execute (cmd);
         H.calcForces = true; // TODO ugly
         H.update(fermi);
         H.calcForces = false; // TODO ugly
      } else {
         H.printEnergies ();
      }

      if (muBasisPtr.getPtr())  {
         //TODO
         //calculate Pulay forces
         cout << "Not yet implemented" << endl;
         SX_EXIT;
      }

      return H.fTotal;
   } else if (dynamic_cast<SxPAWHamiltonian *>(hamPtr.getPtr ()))   {

      // get PAW Hamiltonian
      SxPAWHamiltonian &hamPAW = dynamic_cast<SxPAWHamiltonian &>(*hamPtr);
      
      // get structure and update Gk and G Basis
      SX_CHECK(&R == hamPAW.rPtr);
      SX_CHECK(&G == &hamPAW.rPtr->getGBasis ());
      if (tauChanged)  {
         SxPAWRho atomRho;
         atomRho.atomicChargeDensity (structure, hamPAW.potPtr, R,
                                      hamPAW.pawRho.getNSpin ());
         hamPAW.pawRho -= atomRho;

         if (Hirshfeld)  {
            SxArray<PsiG> atomRhoG (structure.getNSpecies ());
            SX_LOOP(is)  {
               atomRhoG(is) = hamPAW.potPtr->getAtomRhoG(G, (int)is);
            }
            hamPAW.pawRho.pwRho.displaceHirshfeld (x, atomRhoG);
         }

         moveWaves (x);
         structure.copy (x);
         hamPAW.structure = structure;
         G.changeTau (structure);
         wavesPtr->getGkBasis().changeTau (structure);


         atomRho.atomicChargeDensity (structure, hamPAW.potPtr, R, hamPAW.pawRho.getNSpin ());
         hamPAW.pawRho += atomRho;

         // orthonormalize waves
         SxOverlap S = hamPAW.getS ();
         SX_NEW_LOOP (*wavesPtr);
         for (int ik = 0; ik < wavesPtr->getGkBasis ().getNk (); ik++)  {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik)) {
               const SxGBasis &gk = wavesPtr->getGkBasis ()(ik);
               const SxBasis &waveBasis = wavesPtr->getBasis (ik);
               for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
                  // --- recompute projections
                  (*wavesPtr)(iSpin, ik) <<= 
                     waveBasis | (gk | (*wavesPtr)(iSpin, ik));
                  S.orthonormalize (&(*wavesPtr)(iSpin,ik));
               }
            }
         }
      }

      // update RhoTerms and Struct
      if (aoBasisPtr)  {
         // assume that waves have been translated
         hamPAW.compute(*wavesPtr, fermi, true);
      } else {
         /*
         if (hamPAW.pawRho.psCore.getSize () > 0)  {
            SX_CHECK (hamPAW.pawRho.psCore.getSize () == pawRho.pwRho(0).getSize (), );
            // update pseudo-core position
            cout << "Shifting pseudo-cores" << endl;
            // new cores
            SxMeshR newCore = hamPAW.computeCorePS ();
            for (int iSpin = 0; iSpin < hamPAW.nSpin; ++iSpin)  {
               hamPAW.pawRho.pwRho(iSpin).plus_assign_ax (1./nSpin,
                                                          newCore - hamPAW.psCore);
            }
            hamPAW.pawRho.psCore = newCore;
         }
         */
         // --- keep density 
         hamPAW.computeRhoTerms ();
      }

      // perform electronic loop
      if (cmd) execute (cmd);

      // calculate Forces
      hamPAW.computeForces(*wavesPtr, fermi);

      // calculate Pulay forces
      if (muBasisPtr.getPtr())  {
         SxAtomicStructure pulay = structure.getNewStr ();
         pulay.set(Coord(0.0,0.0,0.0));
         cout << SX_SEPARATOR;
         cout << "Calculate Pulay Forces!" << endl;
         cout << SX_SEPARATOR;
         SxOverlap S = hamPAW.getS ();
         int nk = wavesPtr->getGkBasis ().getNk ();
         int nOrbs = muBasisPtr->getNStates ();
         int nWStates = wavesPtr->getNStates ();
         const SxGkBasis &gk = wavesPtr->getGkBasis ();
         const SxPWSet &waves = *wavesPtr;
         const SxMuPW &aoBasis = *muBasisPtr;
         for (int ik = 0; ik < nk; ik++)  {
            SxPAWBasis pawBasis (gk(ik).getThis (), hamPAW.pBasis);
            for (int iSpin = 0; iSpin < nSpin; iSpin++)  {

               // construct eigenvalues
               SxDiracVec<Double> epsilon (nWStates);
               for (int iState = 0; iState < nWStates; iState++)  {
                  PsiG psi = gk(ik) | waves (iState, iSpin, ik);
                  PsiG HPsi = hamPAW * psi; 
                  PsiG SPsi = S * psi; 
                  epsilon (iState) = dot(psi, HPsi) / dot(psi, SPsi);
               }

               // construct overlap
               SxDiracMat<Complex16> Snm (nOrbs, nOrbs);
               for (int iOrbital = 0; iOrbital < nOrbs; iOrbital++)  {
                  PsiG phiI = aoBasis(iOrbital,iSpin,ik);
                  PsiG SPhiI = S * phiI;
                  for (int jOrbital = iOrbital; jOrbital < nOrbs; jOrbital++)  {
                     PsiG phiJ = aoBasis(jOrbital,iSpin,ik);
                     Snm(jOrbital,iOrbital) = dot (phiJ, SPhiI);
                     Snm(iOrbital,jOrbital) = Snm(jOrbital,iOrbital).conj();
                  }
               }
                     
               // construct inverse of Overlap
               SxDiracMat<Complex16> invS = Snm.inverse ();

               // construct coefficients
               SxDiracMat<Complex16> Cnmu(nStates, nOrbs);
               PsiG SPsi = S * (gk(ik) | waves(iSpin,ik));
               for (int iOrbital = 0; iOrbital < nOrbs; iOrbital++)  {
                  Cnmu.colRef(iOrbital) <<= (SPsi.adjoint() ^ aoBasis(iOrbital,iSpin,ik));
               }
               Cnmu = (Cnmu ^ invS).conj ();

               // check PW representation
               /*
               cout << "Checking Plane-Wave representation:" << endl;
               for (int iState = 0; iState < nWStates; iState++)  {
                  PsiG diff = 1.0 *  (gk(ik) | waves (iState, iSpin, ik));
                  for (int iOrbital = 0; iOrbital < nOrbs; iOrbital++)  {
                     diff -= Cnmu(iState,iOrbital) * aoBasis(iOrbital,iSpin,ik);
                  }
                  cout << iState << ": diff is " << diff.norm () << endl;
               }
               */

               // calculate Pulayforces
               for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)  {
                  for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
                     for (int iState = 0; iState < nWStates; iState++)  {
                        PsiG psi = gk(ik) | waves(iState,iSpin,ik);
                        PsiG paw = pawBasis | psi;
                        PsiG HPsi = hamPAW * paw;
                        SPsi = S * paw;
                        for (int iOrbital = 0; iOrbital < nOrbs; iOrbital++)  {
                           SxQuantumNumbers thisOrbital = aoBasis.psiOrb(iOrbital);
                           if ( (thisOrbital.iSpecies == iSpecies) && (thisOrbital.iAtom == iAtom))  {
                              PsiG phiI = aoBasis(iOrbital,iSpin,ik);
                              for (int dir = 0; dir < 3; dir++)  {
                                 pulay.ref(iSpecies, iAtom)(dir) 
                                    += 2.0 * (gk.weights(ik) * fermi.focc(iState, iSpin, ik)
                                          * Cnmu(iState,iOrbital).conj () 
                                          * dot(I * gk(ik).gVec.colRef(dir) * phiI, HPsi - epsilon(iState) * SPsi)).re;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

         // Symmetrize forces
         SxForceSym forceSym;
         forceSym.setup(structure);
         pulay = forceSym * pulay;

         cout << SX_SEPARATOR;
         /*
         cout << "Current Structure :" << endl;
         for (int iSpecies = 0; iSpecies < x.getNSpecies (); iSpecies++)  {
            for (int iAtom = 0; iAtom < x.getNAtoms(iSpecies); iAtom++)  {
               cout << "Species " << iSpecies << ", Atom : " << iAtom << " : " << structure(iSpecies,iAtom) << endl;
            }
         }
         cout << "dist : " << (structure(0,1) - structure(0,0)).norm() << endl;
         cout << "Hellman-Feynman Forces :" << endl;
         for (int iSpecies = 0; iSpecies < x.getNSpecies (); iSpecies++)  {
            for (int iAtom = 0; iAtom < x.getNAtoms(iSpecies); iAtom++)  {
               cout << "Species " << iSpecies << ", Atom : " << iAtom << " : " << hamPAW.forces(iSpecies,iAtom) 
                  << "\t" << hamPAW.forces(iSpecies,iAtom).norm () << endl;
            }
         }
         */
         cout << "Pulay Forces :" << endl;
         for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)  {
            for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
               cout << "Species " << iSpecies << ", Atom : " << iAtom << " : " << pulay(iSpecies, iAtom)
                  << "\t" << pulay(iSpecies,iAtom).norm () << endl;
         }
         }
         cout << "Total Forces :" << endl;
         for (int iSpecies = 0; iSpecies < structure.getNSpecies (); iSpecies++)  {
            for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
               cout << "Species " << iSpecies << ", Atom : " << iAtom << " : " << hamPAW.forces(iSpecies,iAtom) + pulay(iSpecies, iAtom)
                  << "\t" << (hamPAW.forces(iSpecies,iAtom) + pulay(iSpecies,iAtom)).norm () << endl;
            }
         }
         cout << SX_SEPARATOR;

         return hamPAW.forces + pulay;
      }
      
      return hamPAW.forces;

   } else   {
      cout << "unknown Hamiltonian !" << endl;
      SX_EXIT;
   }
}


SxSpeciesData SxHamSolver::getSpeciesData () const
{
   return *potPtr;
}


PrecEnergy SxHamSolver::getEnergy () const
{
   return energy;
}

void SxHamSolver::initialGuess (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK(!hamPtr);
   SxSymbolTable *rhoGroup      = cmd->containsGroup ("rho") 
                                ? cmd->getGroup("rho")
                                : NULL;
   SxSymbolTable *wavesGroup    =  cmd->containsGroup ("waves")
                                ? cmd->getGroup("waves")
                                : NULL;
   SxSymbolTable *exchangeGroup =  cmd->containsGroup ("exchange")
                                ? cmd->getGroup("exchange")
                                : NULL;
   SxSymbolTable *xcPotGroup    = (cmd->containsGroup("xcPotential"))
                                ?  cmd->getGroup("xcPotential")
                                :  NULL;
   SxSymbolTable *occGroup      = cmd->containsGroup("occupations") 
                                ? cmd->getGroup("occupations")
                                : NULL;
   const SxSymbolTable *topLvl = cmd->topLevel ();


   // Setup GkBasis
   SxPtr<SxGkBasis> gkPtr = SxPtr<SxGkBasis>::create(G, topLvl);

   SxPtr<SxPAWHamiltonian> pawHam;
   SxPtr<SxHamInit> initializer;
   /* --- CF 2017-06-02
      It is tempting to put the initializer setup into a
      separate routine. However, for the time being, the
      different SxHamInit specializations have very different
      constructors, and I don't want to provide a wrapper
      interface with loads of mostly unused parameters
      */
   if (dynamic_cast<SxPAWPot*>(potPtr.getPtr ()))  {
      // --- PAW 
      SxPtr<SxInitPAW> initPAW;
      initPAW = initPAW.create (cmd, potPtr, gkPtr, this);
      initializer = initPAW;
      hamPtr = pawHam = initPAW->pawHam; 
   } else if (   topLvl->containsGroup ("PWHamiltonian") 
              || topLvl->containsGroup ("Hamiltonian")) 
   {
      SX_CHECK (dynamic_cast<SxPseudoPot*> (potPtr.getPtr ()));
      initializer = SxPtr<SxInitPW>::create (R, nSpin, nElectrons, potPtr);
   } else if (topLvl->containsGroup ("KP8x8Hamiltonian")) {
      initializer = SxPtr<SxInitKP>::create (SxInitKP::KdotP8x8, &R);
      gkPtr->setNComp(G.getNComp());
   } else if (topLvl->containsGroup ("kpHamiltonian")) {
      cout << "n-band k.p Hamiltonian used here...\n";
      initializer = SxPtr<SxInitKP>::create (SxInitKP::KdotPgeneral, &R);
      gkPtr->setNComp(G.getNComp());
   } else if (topLvl->containsGroup ("StrainField")) {
      initializer = SxPtr<SxInitKP>::create (SxInitKP::StrainField, &R);
      gkPtr->setNComp(G.getNComp());
   } else  {
      sxprintf ("No suitable Hamiltonian found in input file!\n");
      SX_QUIT;
   }
   // --- end of Hamiltonian initializer setup

   // setup spin constraints
   if (topLvl->containsGroup ("spinConstraint")) {
      spinConstraint=SxPtr<SxSpinConstraint>::create(potPtr);
      spinConstraint->read (topLvl, structure);
      pawHam->nuA.resize (structure.getNAtoms ());
      pawHam->nuA.set (0.);
   }

   if (topLvl->containsGroup ("aoBasis"))  {
      const SxSymbolTable *aoGroup = topLvl->getGroup ("aoBasis");
      if (aoGroup->contains ("file"))
         setupAOBasis (aoGroup->get ("file")->toString (), *gkPtr);
      else if (aoGroup->contains ("fromPotential") 
            && aoGroup->get ("fromPotential")->toAttribute ())
         setupAOBasis ("", *gkPtr);
      else
         setupAOBasis ("quamol.sxb", *gkPtr); // falls back to potential if no quamol.sxb
   }

   // --- initialize rho
   if (rhoGroup && rhoGroup->hasAttribute("random"))  {
      initializer->randomRho ();
   } else if (rhoGroup && rhoGroup->hasAttribute("fromWaves"))  {
      if (!wavesGroup || !wavesGroup->contains("file"))  {
         cout << "ERROR for fromWaves: no waves to construct rho from" << endl;
         cout << "The density can only be initialized from waves if they are\n";
         cout << "read from a file." << endl; 
         SX_QUIT;
      }
      // will be performed below
   } else if (!rhoGroup || rhoGroup->hasAttribute("atomicOrbitals"))  {
      if (calc) initializer->atomicRho (rhoGroup, structure);
   } else if (rhoGroup->contains("file"))  {
      SxString file = rhoGroup->get("file")->toString();
      initializer->readRho (file);
   } else {
      SX_EXIT;
   }
   if (rhoGroup && rhoGroup->containsGroup ("charged"))  {
      PsiG extraCharge(G);
      extraCharge.set (0.);
      SYMBOLPARSE(rhoGroup)  {
         FOREACH_SYMBOLGROUP("charged")  {
            double beta = SYMBOLGET("width") || 1;
            double charge = SYMBOLGET("charge");
            if (HAVE_SYMBOL("z"))  {
               // --- charge layer at specific z
               double z = SYMBOLGET("z");
               cout << "Adding Gaussian charge (Q=" << charge
                    << ") layer at z=" << z << endl;
               const SxMesh3D &meshG = G.fft3d(0).mesh;
               SX_LOOP(ig)  {
                  SxVector3<Int> idx3 = meshG.getMeshVec(G.n123(0)(ig),
                                                         SxMesh3D::Positive);
                  if (idx3(0) == 0 && idx3(1) == 0)  {
                     // --- gx = gy = 0
                     double gz = G.gVec((ssize_t)ig,2);
                     extraCharge(ig) += charge
                                      * exp(-0.5 * sqr(gz * beta))
                                      * exp(-I * gz * z);
                  }
               }
            } else {
               // --- Gaussian "point" charge at some position
               Coord pos(SYMBOLGET("coords")->toList ());
               if (SYMBOLGET("relative").toBool ())
                  structure.cell.changeToCar (&pos);
               cout << "Adding Gaussian charge (Q=" << charge
                    << ") at " << pos << endl;
               PsiG T = G.getPhaseFactors (pos);
               SxDiracVec<Double> gauss = exp((-0.5*beta*beta) * G.g2);
               extraCharge.plus_assign_ax (charge, gauss * T);
            }
         }
      }
      extraCharge /= sqrt(R.cell.volume);
      initializer->addCharge (R | extraCharge);
   }

   // --- initialize waves
   fermi = SxFermi (nElectrons, nStates, nSpin, *gkPtr);
   SYMBOLPARSE(SxHamiltonian::getHamiltonianGroup(topLvl))  {
      fermi.orderMethfesselPaxton = SYMBOLGET("MethfesselPaxton") || -1;
   }
   if (wavesGroup && wavesGroup->contains ("file"))  {
      SxString file = wavesGroup->get("file")->toString();
      wavesPtr = initializer->setupWaves (wavesGroup, gkPtr, nStates, nSpin);
      try {
         SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
         wavesPtr->read (io, SxPWSet::KeepAll);
         bool recomputeFermi = wavesGroup->contains ("changeCharge")
                             && wavesGroup->get("changeCharge")->toAttribute ();
         double nElec = fermi.nElectrons;
         if (recomputeFermi) fermi.nElectrons = -1.;
         fermi.read (io, true);
         fermi.kpPtr = &*gkPtr;
         if (recomputeFermi)  {
            cout << "Recomputing occupations..." << endl;
            fermi.nElectrons = nElec;
            fermi.fermiDistribution (ekt);
         }
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
   } else if (pawHam && (!wavesGroup || wavesGroup->containsGroup ("lcao"))) {
      if (calc)  {
         wavesPtr = initializer->setupWaves (wavesGroup, gkPtr, nStates, nSpin);
      }
      SxSymbolTable *lcao = wavesGroup ?  wavesGroup->getGroup ("lcao")
                                       : NULL ; 
      tightBindingPAW ( lcao, calc );
      fermi.fermiDistribution (ekt);
   } else if (!wavesGroup || wavesGroup->containsGroup("lcao")) {
      if (calc)  {
         wavesPtr = initializer->setupWaves (wavesGroup, gkPtr, nStates, nSpin);
      }
   } else if (wavesGroup && wavesGroup->hasAttribute("random")) {
      if (calc)  {
         gkPtr->setNComp(G.getNComp());
         wavesPtr = initializer->setupWaves (wavesGroup, gkPtr, nStates, nSpin);
         initializer->randomize (wavesPtr);
      }
   }
   if (spinConstraint && !dynamic_cast<SxPAWSet*>(wavesPtr.getPtr ()))  {
      cout << "Spin constraints need PAW basis" << endl;
      SX_QUIT;
   }

   // --- read occupations
   if (occGroup) fermi.readOccupations(occGroup);
   
   // --- construct rho from waves
   if (rhoGroup && rhoGroup->hasAttribute("fromWaves"))  {
      if (calc) initializer->rhoFromWaves (wavesPtr, fermi.focc);
   }

   // --- now both waves and rho are initialized
   //     thus, the Hamiltonian can be setup now
   if (calc)  {
      hamPtr = initializer->setupHam (wavesPtr, potPtr, structure);
      if (nStatesChi >= 0)
         SxPtr<SxPWHamiltonian>(hamPtr)->validateNStates (nStatesChi);
   }
   if (!pawHam) hamPtr->read (topLvl); // PAW reads topLvl in SxInitPAW
   if (fabs(SxHamiltonian::getNExcessElectrons(topLvl)) > 1e-10)  {
      // charged systems with dipole correction may have a finite net force
      if (pawHam)  {
         noNetForce = !bool(pawHam->dipoleCorr);
      } else {
         SX_CHECK (dynamic_cast<SxPWHamiltonian*>(hamPtr.getPtr ()));
         noNetForce = !(SxPtr<SxPWHamiltonian> (hamPtr)->dipoleCorrection);
      }
   }

   SxXC::XCFunctional xcBackup = SxXC::Unknown;
   if (xcPotGroup)  {
      if (!dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()))
      {
         cout << "unexpected group xcPotential" << endl;
         SX_QUIT;
      }
      SxPWHamiltonian &H = *dynamic_cast<SxPWHamiltonian *> (hamPtr.getPtr());
      SX_CHECK (H.xcPtr);
      SxXC &xc = *H.xcPtr;
      xcBackup = xc.xcFunctional;
      xc.xcFunctional = SxXC::READ_VXC;
      
      // --- read exchange potential
      if (xcPotGroup->contains("file"))  {
         xc.vXcFile = xcPotGroup->get("file")->toString();
      }
      if (xcPotGroup->containsGroup("compute"))  {
         SxSymbolTable *computeGroup = xcPotGroup->getGroup("compute"),
                       *item;
         SxList<SxSymbolTable *> cmdList;

         // get compute commands
         for (item = computeGroup->begin ();
              item != NULL;
              item = item->nextSibling ())
         {
            if (isRegistered (item)) cmdList << item;
         }
         
         // get xc functional
         xc.xcFunctional = (SxXC::XCFunctional)computeGroup->get("xc")->toInt();
         if (xc.xcFunctional >= SxXC::EXX)  {
            cout << "Illegal xc functional for potential initialization!\n";
            cout << "xc = " << (int) xc.xcFunctional << endl;
            cout << "Just try LDA." << endl;
            SX_QUIT;
         }

         // run DFT calculation
         if (cmdList.getSize () > 0)  {
            // compute something: setup waves if lcao
            if (wavesGroup->containsGroup("lcao"))
               tightBinding (wavesGroup->getGroup("lcao"), calc);
            // run compute commands
            for (int i = 0; i < cmdList.getSize (); ++i)
               execute(cmdList(i), calc);
            // stop further initializations
            wavesGroup = NULL;
            cout << SX_SEPARATOR;
            cout << "| WARNING: rho/waves are taken from potential calculation"
                    " above." << endl;
            cout << SX_SEPARATOR;
         }
      }
   }

   if (wavesGroup)  {
      const SxSymbolTable *lcaoGroup = wavesGroup->containsGroup("lcao")
                                     ? wavesGroup->getGroup("lcao")
                                     : NULL;
      if (!pawHam && lcaoGroup)  {
         SX_CHECK(fermi.getNk() == fermi.kpPtr->nk, fermi.getNk(), fermi.kpPtr->nk);
         tightBinding (lcaoGroup, calc);
         if (occGroup) {
            cout << SX_SEPARATOR;
            cout << "| Setting occupation numbers to initial values." << endl;
            fermi.readOccupations(occGroup);
            cout << SX_SEPARATOR;
         }

         // if rhoGroup->hasAttribute(readRho -> readRho
         if (   rhoGroup  && rhoGroup->contains("file")
             && lcaoGroup->get("rhoMixing")->toReal () > 1e-7)
         {
            SxString file = rhoGroup->get("file")->toString();

            // --- warn that lcao changed rho
            cout << SX_SEPARATOR;
            cout << "| WARNING: lcao rhoMixing ignored - rereading '";
            cout << file << "'." << endl;
            cout << SX_SEPARATOR;

            if (calc)  {
               hamPtr->getRho().readRho (file);
            }
         }
      }

      if (calc)  {
         if (pawHam)  {
            pawHam->computeRhoTerms ();
            if (SxXC::isHybrid (pawHam->xcPtr->xcFunctional))  {
               if (!exchangeGroup)  {
                  pawHam->exchangePtr->wavesPtr = wavesPtr;
                  pawHam->exchangePtr->setFocc (fermi.focc);
                  pawHam->exchangePtr->computeXij (pawHam->pawRho.Dij);
               }
            }
         } else {
            hamPtr->compute (fermi, true, true);
         }
      }
   } else if (!pawHam)  {
      // --- default: tb initialization
      tightBinding (NULL, calc);
      if (occGroup) {
         cout << SX_SEPARATOR;
         cout << "| Setting occupation numbers to initial values." << endl;
         fermi.readOccupations(occGroup);
         cout << SX_SEPARATOR;
      }
      if (calc)  {
         hamPtr->compute (fermi, true, true);
      }
   }

   // restore xc functional
   if (xcPotGroup)  {
      //H.xcFunctional = xcBackup;
      SxPWHamiltonian *pwHam = dynamic_cast<SxPWHamiltonian*>(hamPtr.getPtr ());
      SX_CHECK (pwHam);
      pwHam->xcPtr->xcFunctional = xcBackup;
   }
   // --- print results
   cout << SX_SEPARATOR;
   cout << "| Initialization done\n";
   cout << SX_SEPARATOR;
   hamPtr->printEnergies (); fermi.printOccupation (); 
   cout << SX_SEPARATOR;


}



void SxHamSolver::tightBinding (const SxSymbolTable *cmd, bool calc)
{
   SX_ALLOC_CACHE;
   SX_CLOCK (Timer::LcaoTotal);

   SX_CHECK (calc);
 
   SxString TBInFile;
   if (cmd)  {
      try {
         maxSteps   = cmd->contains ("maxSteps")
                    ? cmd->get("maxSteps")->toInt()
                    : 1;
         rhoMixing  = cmd->contains ("rhoMixing")
                    ? cmd->get("rhoMixing")->toReal()
                    : 0.;
         dEnergy = cmd->contains ("dEnergy")
                    ? cmd->get("dEnergy")->toReal()
                    : 1e-8;
         foccMixing = 0.05;
         if (cmd->contains("foccMixing"))
            foccMixing = cmd->get("foccMixing")->toReal();
         if (cmd->contains("ekt"))  ekt = cmd->get("ekt")->toReal()/HA2EV;
         else                       ekt = ektDefault;
         TBInFile = cmd->contains ("file")
                  ? cmd->get("file")->toString()
                  : SxString("");
      } catch (const SxException &e) {
         e.print ();
         SX_EXIT;
      }
   } else {
      // --- default run
      maxSteps   = 1;
      rhoMixing  = 0.;
      foccMixing = 1.;
      ekt = ektDefault;
      TBInFile = "";
   }

   if (!dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()))  {
      cout << "lcao initialization only for PW Hamiltonians!" << endl;
      SX_QUIT;
   }
   SX_CHECK (dynamic_cast<SxPW *> (wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &> (*wavesPtr);
   SxPWHamiltonian &hamPW = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxPseudoPot*>(potPtr.getPtr ()));
   SxPseudoPot &psPot     = dynamic_cast<SxPseudoPot&>(*potPtr);

   // --- basis sets and atomic orbitals
   SxConstPtr<SxRadBasis> rPtr;
   SxGkBasis &Gk = waves.getGkBasis();
   SxPtr<SxGkBasis> gkPtr = waves.getGkBasisPtr();
   SX_CHECK(fermi.getNk() == fermi.kpPtr->nk, fermi.getNk(), fermi.kpPtr->nk);
   SxAtomicOrbitals TBOrbitals;
   SxPtr<SxMuPW> muPtr;
   if(TBInFile != "")  {
      try {
         cout << "Taking TBOrbitals from " << TBInFile << endl; 
         SxBinIO io(TBInFile, SxBinIO::BINARY_READ_ONLY);
         TBOrbitals = SxAtomicOrbitals(io);
         io.close ();
      } catch (const SxException &e) {
         e.print();
         SX_EXIT;
      }
      muPtr = muPtr.create (TBOrbitals, gkPtr);
   } else {
      cout << "Taking TBOrbitals from Potential." << endl;
      rPtr = SxConstPtr<SxRadBasis>::create(psPot.rad, psPot.logDr);
      TBOrbitals = SxAtomicOrbitals(psPot.pseudoPsi, rPtr);
      muPtr = muPtr.create (TBOrbitals, psPot.pseudoFocc, gkPtr);
   }
   SxMuPW &mu = *muPtr; // sum_r<G+k|r><r|mu>

   // --- check number of required orbitals to initialize the Bloch states
   int nOrb      = mu.getNStates ();
   {
      SX_MPI_SOURCE(TopLevel, TaskGroupMaster);
      SX_MPI_TARGET(TopLevel, TaskGroupAll);
      nOrb = SxLoopMPI::bcast (nOrb, 0);
   }
   int nPWStates = waves.getNStates ();
   int nFromTB   = minimum (nOrb, nPWStates);

   // --- print information
   cout << SX_SEPARATOR;
   cout << "| Tight Binding Initialization\n";
   cout << SX_SEPARATOR;
   sxprintf ("| Number of atomic orbitals:   %d\n", nOrb);
   sxprintf ("| Number of Bloch states:      %d\n", nPWStates);
   if (nOrb < nPWStates)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Warning: %d Bloch states will be initialized with "
                        "random numbers!\n", nPWStates - nFromTB);
      sxprintf ("|          If you need a complete LCAO basis just\n");
      sxprintf ("|             - reduce number of empty states, or\n");
      sxprintf ("|             - use more localized orbitals\n");
      waves.randomize ();
   }
   cout << SX_SEPARATOR;

   // --- Mixer and TB Hamiltonian = PW Hamiltonian with |mu>
   SxRhoMixer mixer (cmd, nSpin ==2);
   mixer.print ();
   cout << SX_SEPARATOR;
   SxPWHamiltonian::Contrib origContrib, potContrib , nonPotContrib;
   origContrib   = hamPW.contrib;
   potContrib    = SxPWHamiltonian::Contrib(origContrib
                 & (SxPWHamiltonian::CALC_EFF | SxPWHamiltonian::CALC_EXT));
   nonPotContrib = SxPWHamiltonian::Contrib(origContrib
                 & (SxPWHamiltonian::CALC_KIN | SxPWHamiltonian::CALC_NL));
   SxPWHamiltonian::Contrib statContrib, nonstatContrib;
   if (hamPW.rsProjector)  {
      // all but nl
      nonstatContrib = SxPWHamiltonian::Contrib(origContrib
                     & (  SxPWHamiltonian::CALC_EFF 
                        | SxPWHamiltonian::CALC_EXT
                        | SxPWHamiltonian::CALC_KIN ));
      statContrib    = SxPWHamiltonian::CALC_NONE;
   } else {
      nonstatContrib = potContrib;
      statContrib    = nonPotContrib;
   }
   hamPW.contrib = potContrib;
   fermi.resize (nOrb);
   hamPW.set (mu, fermi);
   fermi.resize (nPWStates);
   hamPW.contrib = origContrib;

   double foccMix;
   PsiG muI, muJ, H_muI, Hstat_muI;
   int nk       = waves.getNk();
   SX_CHECK(nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SxMatrix<TPrecCoeffG> H(nOrb,nOrb), S, L;
   SxSymMatrix<TPrecCoeffG> Hsym(nOrb);
   SxSymMatrix<TPrecCoeffG>::Eigensystem eig;
   
   SxArray<SxMatrix<TPrecCoeffG> > Lk(nk), Hkstat(nk);

   bool computeStatic, computeL;
   double deltaE = 1.0, newEnergy = 1000;
   SX_MPI_LEVEL("waves-k");
   for (int it = 0; it < maxSteps; it++)  {
      SX_CLOCK (Timer::LcaoStep);
      mixer.addRhoIn (hamPW.rho);
      FILE *SEigFilePtr;
      SxString SEigFile = "OverlapValuesOverK.dat";
      if ( (SEigFilePtr = fopen(SEigFile.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << SEigFile.ascii() << endl;
         SX_EXIT;
      }
      for (int ik = 0; ik < nk; ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
            computeStatic = (it == 0 && iSpin == 0);
            computeL      = (it == 0 && iSpin == 0);

            if (computeStatic) 
               Hkstat(ik) = SxMatrix<TPrecCoeffG> (nOrb, nOrb);
            if (computeL) S.reformat (nOrb, nOrb);
            
            if (hamPW.rsProjector && computeStatic)  {
               // --- compute nl-part via fast EES-G

               computeStatic = false;

               const SxArray<PsiG> &muRef = mu.muPhi(ik);
               Hkstat(ik).set (0.);

               int nMu = int(muRef.getSize ());
               SxArray<SxVector<TPrecCoeffG> > pPhiMu(nMu);
               SxArray<SxConstPtr<SxAtomInfo> > info(nMu);
               Coord kVec = Gk.getK (ik);
               SxAtomicStructure atPhi, atMu;
               for (int iPhi = 0; iPhi < hamPW.phiNl(ik).getSize (); ++iPhi) {
                  atPhi = SxSpeciesRef(hamPW.phiNl(ik)(iPhi)
                                       .handle->auxData.is  ) | structure;

                  // --- set up <phi @ ka|mu @ ia>
                  for (int iMu = 0; iMu < nMu; ++iMu)  {
                     atMu = SxSpeciesRef(muRef(iMu).handle->auxData.is)
                            | structure;
                     
                     // --- set up relative distances
                     SxAtomicStructure distPhiMu;
                     for (int ia = 0; ia < atMu.getNAtoms (); ++ia)  {
                        distPhiMu.newSpecies ();
                        for (int ka = 0; ka < atPhi.getNAtoms (); ++ka)  {
                           distPhiMu.addAtom (atPhi(ka) - atMu(ia));
                        }
                     }
                     distPhiMu.endCreation ();
                     SX_CHECK (distPhiMu.getNSpecies() == atMu.getNAtoms(),
                               distPhiMu.getNSpecies(), atMu.getNAtoms ());
                     // store info
                     info(iMu) = distPhiMu.atomInfo;
                     
                     // get projection
                     pPhiMu(iMu) = hamPW.rsProjector->project 
                                   (hamPW.phiNl(ik)(iPhi),
                                    muRef(iMu), 
                                    distPhiMu);

                     // --- add mu phase exp(-ikr_mu) for consistency
                     //     with mu(i,iSpin,ik)
                     SxVector<TPrecCoeffG>::Iterator pIt = pPhiMu(iMu).begin ();
                     SxComplex16 phase;
                     for (int ia = 0; ia < atMu.getNAtoms (); ++ia)  {
                        phase = SxEESGProj::getPhase (-(kVec ^ atMu(ia)));
                        for (int ka = 0; ka < atPhi.getNAtoms (); ++ka)
                           *pIt++ *= phase;
                     }
                     // phi phase exp(ikr_phi) cancels out
                  }

                  // --- add phi contribution to Hamiltonian
                  for (int iOrb = 0; iOrb < nOrb; ++iOrb)  {
                     int iMu = mu.getRefIdx (iOrb);
                     SX_CHECK (iMu != -1);
                     int ia = mu.psiOrb(iOrb).iAtom;
                     for (int jOrb = 0; jOrb < nOrb; ++jOrb)  {
                        int jMu = mu.getRefIdx(jOrb);
                        SX_CHECK (jMu != -1);
                        int ja = mu.psiOrb(jOrb).iAtom;
                        for (int ka = 0; ka < atPhi.getNAtoms (); ++ka)  {
                           int pIdxI = info(iMu)->offset(ia) + ka;
                           int pIdxJ = info(jMu)->offset(ja) + ka;

                   // sum_ka <mu @ ia|phi @ ka> E^{-1}(phi) <phi @ ka| mu @ ja>
                           Hkstat(ik)(iOrb, jOrb) += pPhiMu(iMu)(pIdxI).conj ()
                                                   * pPhiMu(jMu)(pIdxJ) 
                                                   / hamPW.eKB(iPhi);
                        }
                     }
                  }
               }
            }
            
            for (int i = 0; i < nOrb; i++)  {
               muI = mu(i,iSpin,ik);                 //  |mu_i>

               // --- H|i>
               if (computeStatic)  {
                  hamPW.contrib = statContrib;
                  Hstat_muI = hamPW * muI;           // H_stat|mu_i>
               }
                  
               hamPW.contrib = nonstatContrib;
               H_muI = hamPW * muI;                  // H_nonstat|mu_i>

               // --- compute <j|H|i>
               
               if (hamPW.rsProjector)  {
                  // --- using fast EES-G techniques to project H|i> onto |j>
                  const SxArray<PsiG> &muRef = mu.muPhi(ik);
                  Coord kVec = Gk.getK (ik);
                  SxArray<SxVector<TPrecCoeffG> > pMuJ(muRef.getSize ());

                  // collect projections
                  for (int iMu = 0; iMu < muRef.getSize (); ++iMu)
                     pMuJ(iMu) = hamPW.rsProjector->project(muRef(iMu), H_muI,
                           SxSpeciesRef(muRef(iMu).handle->auxData.is) | structure);
                  // --- put them into the right place within H
                  for (int j = i; j < nOrb; ++j)  {
                     H(j, i) = Hkstat(ik)(j,i) 
                                + pMuJ(mu.getRefIdx(j))(mu.psiOrb(j).iAtom) 
                                // missing muJ phase
                                * SxEESGProj::getPhase(kVec ^ structure.getAtom(
                                         mu.psiOrb(j).iSpecies,
                                         mu.psiOrb(j).iAtom));
                     H(i, j) = H(j, i).conj ();
                  }
                     
               } else {
                  for (int j = i; j < nOrb; j++)  {
                     muJ = mu(j,iSpin,ik);               //   |mu_j>
                     // --- static matrix element
                     if (computeStatic)  {
                        Hkstat(ik)(j, i) = (muJ ^ Hstat_muI).chop ();
                        Hkstat(ik)(i, j) = Hkstat(ik)(j, i).conj ();
                     }
                     
                     // --- full matrix element
                     H(j, i) = Hkstat(ik)(j,i) + (muJ ^ H_muI).chop ();
                     H(i, j) = H(j, i).conj ();

                     // --- overlap matrix element
                     if (computeL && !hamPW.rsProjector)  {
                        S(j,i) = (muJ ^ muI).chop ();
                        S(i,j) = S(j,i).conj ();
                     }
                  }
               }
            }
            if (computeL && !hamPW.rsProjector)   {
               // print out TightBinding Hamiltonian and Overlap
               /*
               cout << SX_SEPARATOR;
               cout << "Tight Binding Hamiltonian: " << endl;
               H.print(true);
               cout << "Tight Binding Overlap: " << endl;
               S.print(true);
               cout << SX_SEPARATOR;
               */
               fprintf(SEigFilePtr, "%i\t",ik);
               SxMatrix<Complex16>::Eigensystem Seig;
               Seig = S.eigensystem ();
               for (int i = 0; i < nOrb; i++) {
                  fprintf(SEigFilePtr, "%f\t",Seig.vals(i).re);
               }
               fprintf(SEigFilePtr, "\n");
            }
            // restore original contrib
            hamPW.contrib = origContrib;

            if (computeL)  {
               if (hamPW.rsProjector)  {
                  // --- compute overlap matrix via EES-G
                  const SxArray<PsiG> &muRef = mu.muPhi(ik);
                  Coord kVec = Gk.getK (ik);
                  int nMu = int(muRef.getSize ());

                  // --- setup lookup table (iMu, ia)->(iOrb)
                  SxArray<SxArray<int> > muOrbs(nMu);
                  for (int iMu = 0; iMu < nMu; ++iMu) 
                     muOrbs(iMu).resize(structure.getNAtoms(muRef(iMu).handle->auxData.is));
                  for (int iOrb = 0; iOrb < nOrb; ++iOrb)
                     muOrbs(mu.getRefIdx(iOrb))(mu.psiOrb(iOrb).iAtom) = iOrb;
                  
                  SxAtomicStructure atMu, atNu;
                  SxVector<TPrecCoeffG> muNu;
                  SxComplex16 ovlp;
                  // --- loop over orbital types
                  for (int iNu = 0; iNu < nMu; ++iNu)  {
                     atNu = SxSpeciesRef(muRef(iNu).handle->auxData.is) 
                            | structure;

                     for (int iMu = iNu; iMu < nMu; ++iMu)  {
                        atMu = SxSpeciesRef(muRef(iMu).handle->auxData.is)
                               | structure;
                        
                        // --- set up relative distances
                        SxAtomicStructure distMuNu (  atMu.getNAtoms () 
                                                    * atNu.getNAtoms ());
                        int iMuNu = 0;
                        for (int ia = 0; ia < atMu.getNAtoms (); ++ia)
                           for (int ja = 0; ja < atNu.getNAtoms (); ++ja, ++iMuNu)
                              distMuNu.ref(iMuNu) = atMu(ia) - atNu(ja);
                        
                        // get projection
                        muNu = hamPW.rsProjector ->project (muRef(iMu),
                                                            muRef(iNu),
                                                            distMuNu);

                        // --- transfer muNu to resulting overlap matrix
                        iMuNu = 0;
                        for (int ia = 0; ia < atMu.getNAtoms (); ++ia)  {
                           int iOrb = muOrbs(iMu)(ia);
                           for (int ja = 0; ja < atNu.getNAtoms (); ++ja, ++iMuNu) {
                              int jOrb = muOrbs(iNu)(ja);
                              // --- add mu phases exp(-ikr_mu) for consistency
                              //     with mu(i,iSpin,ik)
                              ovlp = muNu(iMuNu)
                                   * SxEESGProj::getPhase(kVec^distMuNu(iMuNu));
                              S(iOrb, jOrb) = ovlp;
                              S(jOrb, iOrb) = ovlp.conj ();
                           }
                        }
                     }
                  }
               }
               // --- compute Cholesky decomposition
               SxMatrix<Complex16> SInv = S.inverse ();
               S = SxMatrix<Complex16> ();
               Lk(ik) = L = SInv.choleskyDecomposition ().adjoint ();
               
               // --- test L
               SxComplex16 llMat = ((L.adjoint()^L) - SInv ).normSqr ();
               if (fabs(llMat.re) > 1e-10 || fabs(llMat.im) > 1e-10)  {
                  sxprintf ("| Warning: L^t*L - S^-1 = (%e,%e) (should be 0)\n", 
                          llMat.re, llMat.im);
               }
            }

            // --- solve eigenproblem
            L = Lk(ik);
            H = (L ^ H ^ L.adjoint());  
            if (!H.isHermitian())  sxprintf ("| TB: H is not hermitian\n");
            for (int i = 0; i < nOrb; i++)
               for (int j = i; j < nOrb; j++)
                  Hsym(i,j) = H(i,j);
            eig = Hsym.eigensystem ();

            eig.vecs = L.adjoint() ^ eig.vecs;

            // --- change to tight-binding basis
            for (int i = 0; i < nFromTB; i++)  {
               fermi.eps(i,iSpin,ik) = eig.vals(i);
               waves(i,iSpin,ik).set (0.);
            }
            // --- higher States are large in energy
            for (int j = 0; j < nOrb; j++)  {
               if (!hamPW.psPot.realSpace (mu.psiOrb(j).iSpecies))  {
                  // --- reciprocal-space summation
                  muJ = mu(j,iSpin,ik);
                  for (int i = 0; i < nFromTB; i++)  {
                     //waves(i,iSpin,ik) += muJ * eig.vecs(j,i);
                     waves(i,iSpin,ik).plus_assign_ax(eig.vecs(j,i), muJ);
                  }
               }
            }
            // --- add real-space orbitals via EES-G
            if (hamPW.rsProjector)  {
               const SxArray<PsiG> &muRef = mu.muPhi(ik);
               int nMu = int(muRef.getSize ());
               SxArray<SxVector<TPrecCoeffG> > eigVecRS(nMu);
               
               // --- resize eigVecRS
               for (int iMu = 0; iMu < nMu; ++iMu)  {
                  int is = muRef(iMu).handle->auxData.is;
                  if (hamPW.psPot.realSpace(is))
                     eigVecRS(iMu).resize (structure.getNAtoms(is));
               }

               Coord tau;
               SxComplex16 phase;
               int refIdx;
               // --- do the wave setup
               for (int i = 0; i < nFromTB; ++i)  {
                  // --- collect relevant eigenvector elements
                  for (int j = 0; j < nOrb; ++j)  {
                     const SxQuantumNumbers &orb = mu.psiOrb(j);
                     if (hamPW.psPot.realSpace(orb.iSpecies))   {
                        tau = structure.constRef(orb.iSpecies, orb.iAtom);
                                                 
                        // missing global mu phase exp(-i k tau_mu)
                        phase = SxEESGProj::getPhase (-(Gk.getK(ik) ^ tau));
                        refIdx = mu.phiIdx(j);
                        eigVecRS(refIdx)(orb.iAtom) = eig.vecs(j,i) * phase;
                     }
                  }
                  // --- sum local orbitals using EES-G
                  for (int iMu = 0; iMu < nMu; ++iMu)  {
                     int is = muRef(iMu).handle->auxData.is;
                     if (hamPW.psPot.realSpace(is))  {
                        SxAtomicStructure atoms = SxSpeciesRef(is) | structure;
                        waves(i,iSpin,ik) += hamPW.rsProjector
                                             ->gradient (eigVecRS(iMu),
                                                         muRef(iMu), atoms);
                     }
                  }
               }
            }
         }
         if (it + 1 >= maxSteps)  {
            // --- free memory
            Lk(ik)     = SxMatrix<Complex16> ();
            Hkstat(ik) = SxMatrix<Complex16> ();
         }
      }
      fermi.eps.synMPI ();
      hamPW.contrib = origContrib;
      waves.orthonormalize ();

      // set first mixing to 1. if mixing is not zero
      foccMix = (it == 0 && foccMixing > 1e-7) ? 1. : foccMixing;
      SX_CHECK(fermi.getNk() == fermi.kpPtr->nk, fermi.getNk(), fermi.kpPtr->nk);
      SX_CHECK(fermi.getNk() == Gk.getNk(), fermi.getNk(), Gk.getNk());
      fermi.fermiDistribution (ekt, foccMix);
      if (it % 10 == 0) fermi.printOccupation (); 

      if (rhoMixing > 1e-10) {
         hamPW.computeRho (fermi.focc, waves);
         mixer.addRhoOut (hamPW.rho);
         hamPW.rho = mixer.getMixedRho ();
         hamPW.set (waves, fermi);
         if (it % 10 == 0) hamPW.printEnergies ();

         //cout << it+1 << ": eTot=" << hamPW.eTotal << endl;
      }

      // SCF LCAO 
      double oldEnergy = newEnergy;
      newEnergy = hamPW.getEnergy();
      deltaE = fabs(newEnergy - oldEnergy);
      if (it == 0 && maxSteps > 1)  {
         cout << SX_SEPARATOR;
         cout << "Selfconsistent LCAO run, maxSteps = " << maxSteps << endl;
         cout << SX_SEPARATOR;
         cout << it+1 << " Step: Etot = " << newEnergy << endl;
      }
      else cout << it+1 << " Step: Etot = " << newEnergy << " , Delta = " << deltaE << endl;
      if (deltaE < dEnergy) break;
   }

   cout << SX_SEPARATOR;
   cout << endl << endl << endl;
   hamPW.contrib = SxPWHamiltonian::CALC_NONE;
   hamPW.set (waves, fermi);
   hamPW.contrib = origContrib;
}

void SxHamSolver::tightBindingPAW (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SxGkBasis &Gk = waves.getGkBasis();
   SxPtr<SxGkBasis> gkPtr = waves.getGkBasisPtr ();

   SX_CLOCK (Timer::LcaoTotal);
   SX_CHECK (calc == true);

   SxString TBInFile;
   int maxStepsTB = 0;
   bool basisDecomp = false;
   bool pulayForces = false;
   bool initializeRho = false;
   if (cmd)  {
      try {
         TBInFile = cmd->contains ("file")
                  ? cmd->get("file")->toString()
                  : SxString("");
         maxStepsTB = cmd->contains ("maxSteps")
                  ? cmd->get("maxSteps")->toInt()
                  : 0;
         dEnergy = cmd->contains ("dEnergy")
                  ? cmd->get("dEnergy")->toReal()
                  : 1e-8;
         basisDecomp = cmd->hasAttribute ("basisDecomposition");
         pulayForces = cmd->hasAttribute ("pulayForces");
         if (cmd->contains ("atomicOrbitals", true))
            initializeRho = cmd->get ("atomicOrbitals", true)->toAttribute ();
      } catch (const SxException &e) {
         e.print ();
         SX_EXIT;
      }
   } else {
      // --- default run
      TBInFile = "";
   }

   // --- get PAW Hamiltonian, overlap operator, and potential
   SX_CHECK (hamPtr);
   SX_CHECK (dynamic_cast<SxPAWHamiltonian *>(hamPtr.getPtr ()));
   SxPAWHamiltonian &hamPAW = dynamic_cast<SxPAWHamiltonian &> (*hamPtr);
   SxOverlap S = hamPAW.getS ();
   SxPAWPot &potPAW = *(hamPAW.potPtr);

   if (initializeRho)  {
      cout << "Reinitializing the density from atoms" << endl;
      hamPAW.pawRho.atomicChargeDensity (structure, potPtr, R, nSpin);
   }
   // --- set up intermediate waves for exchange
   bool updateXWaves = hamPAW.exchangePtr && !hamPAW.exchangePtr->wavesPtr;
   if (updateXWaves)  {
      // get RadBasis
      const SxRadBasis &rad = *potPAW.getBasisPtr ();
      // --- collect LCAO orbitals with initial occupation > 0
      SxArray<SxArray<SxDiracVec<Double> > > usedPhiPS(structure.getNSpecies());
      for (int is = 0; is < structure.getNSpecies(); is++)   {
         int nProjLocal = (int)potPAW.foccInit(is).getSize();
         SxList<int> usedProjectors;
         for (int ipl = 0; ipl < nProjLocal; ipl++)   
            if (potPAW.foccInit(is)(ipl) > 1e-10) usedProjectors.append(ipl);
         int nUsedProjectors = int(usedProjectors.getSize());
         usedPhiPS(is).resize(nUsedProjectors);
         for (int iupl = 0; iupl < nUsedProjectors; iupl++)   {
            int ipl = usedProjectors(iupl);
            usedPhiPS(is)(iupl) = potPAW.phiPS (is).colRef (ipl);
            usedPhiPS(is)(iupl).setBasis (rad);
            usedPhiPS(is)(iupl).handle->auxData.is = is;
            usedPhiPS(is)(iupl).handle->auxData.n = ipl;
            usedPhiPS(is)(iupl).handle->auxData.l = potPAW.lPhi(is)(ipl);
            usedPhiPS(is)(iupl).handle->auxData.m = 0;
            cout << is << "\t" << ipl << "\t" << potPAW.lPhi(is)(ipl) << endl;
         }
      }
      SxAtomicOrbitals xOrbitals(usedPhiPS, potPAW.getBasisPtr ());
      SxMuPW xAoBasis (xOrbitals, gkPtr);
      int nExStates = xAoBasis.getNStates (), nk = Gk.getNk ();
      SxPtr<SxPW> xWaves = SxPtr<SxPW>::create (nExStates, nSpin, gkPtr);
      Focc focc(nExStates, nSpin, nk);
      for (int ik = 0; ik < nk; ++ik)  {
         SX_MPI_LEVEL("waves-k");
         if (! SxLoopMPI::myWork(ik)) continue;
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            for (int i = 0; i < nExStates; ++i)  {
               PsiG ao = xAoBasis(i,iSpin,ik);
               (*xWaves)(i,iSpin,ik) << ao;
               int is  = ao.handle->auxData.is,
                   ipl = ao.handle->auxData.n,
                   l   = ao.handle->auxData.l;
               focc(i, iSpin, ik) = potPAW.foccInit(is)(ipl) / (2 * l + 1);
            }
            S.orthonormalize (&(*xWaves)(iSpin, ik), Loewdin);
         }
      }

      hamPAW.exchangePtr->wavesPtr = xWaves;
      hamPAW.exchangePtr->setFocc (focc);
      hamPAW.exchangePtr->computeXij (hamPAW.computeDij (*xWaves, focc));
   }

   cout << SX_SEPARATOR;
   cout << "Setup initial Guess "<< endl;
  
   // --- define TB orbitals in radial basis
   SxAtomicOrbitals TBOrbitals;
   SxConstPtr<SxRadBasis> radBasisPtr;
   if(TBInFile.getSize () > 0)  {
      // --- read from file
      cout << "Taking TBOrbitals from File " << TBInFile << endl;
      radBasisPtr = SxConstPtr<SxRadBasis>::create(TBInFile);
      TBOrbitals.setBasis (radBasisPtr);
      TBOrbitals.read (TBInFile);
      SX_MPI_MASTER_ONLY {
         TBOrbitals.print();
      }
   } else   {
      cout << "Taking TBOrbitals from Potential" << endl;
      cout << "Take the following PhiPS:" << endl;
      cout << "is\tnProjLocal\tlPhi" << endl;

      // get RadBasis
      radBasisPtr = potPAW.getBasisPtr ();
      // --- collect LCAO orbitals with initial occupation > 0
      SxArray<SxArray<SxDiracVec<Double> > > usedPhiPS(structure.getNSpecies());
      for (int is = 0; is < structure.getNSpecies(); is++)   {
         int nProjLocal = (int)potPAW.foccInit(is).getSize();
         SxList<int> usedProjectors;
         for (int ipl = 0; ipl < nProjLocal; ipl++)   
            if (potPAW.foccInit(is)(ipl) > 1e-10) usedProjectors.append(ipl);
         int nUsedProjectors = int(usedProjectors.getSize());
         usedPhiPS(is).resize(nUsedProjectors);
         for (int iupl = 0; iupl < nUsedProjectors; iupl++)   {
            int ipl = usedProjectors(iupl);
            usedPhiPS(is)(iupl) = potPAW.phiPS (is).colRef (ipl);
            usedPhiPS(is)(iupl).setBasis (radBasisPtr.getConstPtr());
            usedPhiPS(is)(iupl).handle->auxData.is = is;
            usedPhiPS(is)(iupl).handle->auxData.n = ipl;
            usedPhiPS(is)(iupl).handle->auxData.l = potPAW.lPhi(is)(ipl);
            usedPhiPS(is)(iupl).handle->auxData.m = 0;
            cout << is << "\t" << ipl << "\t" << potPAW.lPhi(is)(ipl) << endl;
         }
      }
      TBOrbitals = SxAtomicOrbitals(usedPhiPS, radBasisPtr);
      SX_MPI_MASTER_ONLY {
         TBOrbitals.print();
      }
   }

   // TB orbitals in G+k basis
   SxAOBasis initBasis(*gkPtr, TBOrbitals);

   if (pulayForces) muBasisPtr = SxConstPtr<SxMuPW>::create (TBOrbitals, gkPtr);

   // Setup PAW Hamiltonian from PAW density
   hamPAW.computeRhoTerms ();

   // --- check Number of required orbitals to initialize the Bloch states
   int nOrb      = initBasis.getNOrb ();
   int nPWStates = wavesPtr->getNStates ();
   int nFromTB   = min(nStates,nOrb);

   // --- print information
   cout << SX_SEPARATOR;
   cout << "| Tight Binding Initialization\n";
   cout << SX_SEPARATOR;
   sxprintf ("| Number of atomic orbitals:   %d\n", nOrb);
   sxprintf ("| Number of Bloch states:      %d\n", nPWStates);
   if (nOrb < nPWStates)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Warning: %d Bloch states will be initialized with "
                        "random numbers!\n", nPWStates - nFromTB);
      sxprintf ("|          If you need a complete LCAO basis just\n");
      sxprintf ("|             - reduce number of empty states, or\n");
      sxprintf ("|             - use more localized orbitals\n");
   }
   cout << SX_SEPARATOR;

   SxRhoMixer mixer (cmd,nSpin == 2);
   double deltaE = 1.0;
   if (mixer.rhoMixing < 1e-10) maxStepsTB = 0;

   for(int iMixSteps = 0; iMixSteps < maxStepsTB || iMixSteps == 0; iMixSteps++)
   {
      SX_NEW_LOOP (waves);
      
      for (int ik = 0; ik < Gk.getNk (); ik++)
      {
         for(int iSpin = 0; iSpin < nSpin; iSpin++)
         {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik))  {
               SxPAWBasis pawBasis (Gk(ik).getThis (), hamPAW.pBasis);
               SxDiracMat<Complex16> Hnm(nOrb,nOrb), Snm(nOrb,nOrb);
               // --- construct Hamiltonian Hmn and Overlap Smn
               // H = H.adjoint, S = S.adjoint ()
               PsiG HPsiJs, SPsiJs;
               int blockSize = min(64, nOrb);
               for (int jState = 0, iB = 0; jState < nOrb; jState++, iB++)   {
                  if (iB == blockSize)  {
                     iB = 0;
                     blockSize = min(blockSize, nOrb - jState);
                  }
                  if (iB == 0)  {
                     PsiG psiJs = SxDiracMat<Complex16> (Gk(ik).ng, blockSize);
                     psiJs.handle->auxData.ik = ik;
                     psiJs.handle->auxData.iSpin = iSpin;
                     psiJs.setBasis (Gk(ik));
                     for (int k = 0; k < blockSize; ++k)
                        psiJs.colRef (k) <<= initBasis.getAOinG(ik, jState+k);
                     psiJs = pawBasis | psiJs;
                     HPsiJs = initBasis | (hamPAW * psiJs);
                     SPsiJs = initBasis | (S * psiJs);
                  }

                  Hnm.colRef(jState) <<= HPsiJs.colRef(iB);
                  Snm.colRef(jState) <<= SPsiJs.colRef(iB);
               }
               // transform generalized eigenvalue problem to standard
               // eigenvalue problem via cholesky decomposition
               SxDiracMat<Complex16> Lnm = Snm.choleskyDecomposition();
               SxDiracMat<Complex16> invLnm = Lnm.inverse();
               // transform Hnm
               Hnm = invLnm ^ Hnm ^ invLnm.adjoint();

               // solve standard eigenvalue problem
               SxDiracMat<Complex16>::Eigensystem eig;
               eig = Hnm.eigensystem ();

               // --- set up waves
               eig.vecs = invLnm.adjoint() ^ eig.vecs;
               // note: waveBasis might be G+k or PAW
               const SxBasis &waveBasis = waves.getBasis (ik);
               blockSize = 64;
               for (int iState = 0; iState < nFromTB; iState += blockSize)   {
                  blockSize = min(64, nFromTB - iState);
                  PsiG vecsBlock = eig.vecs(SxIdx(iState * nOrb, 
                                                  (iState + blockSize)*nOrb-1));
                  vecsBlock.reshape (nOrb, blockSize);
                  vecsBlock.setBasis (initBasis);
                  waves.getBlock (iState, blockSize, iSpin, ik)
                     <<= waveBasis | (Gk(ik) | vecsBlock);
               }
               for (int iState = 0; iState < nFromTB; ++iState)
                  fermi.eps(iState,iSpin, ik) = eig.vals(iState).re;

               // --- Setup higher waves
               if (nStates > nOrb)   {
                  int nHigh = nStates-nFromTB;
                  SxDiracMat<TPrecCoeffG> highGk(Gk(ik).ng, nHigh),
                        highStates, tbStates;
                  tbStates   = waves.getBlock (0, nFromTB, iSpin, ik);
                  highStates = waves.getBlock (nFromTB, nHigh, iSpin, ik);

                  // fill with random numbers
                  highGk.randomize ();
                  // --- project to wave basis
                  highGk.handle->auxData.ik = ik;
                  highGk.handle->auxData.iSpin = iSpin;
                  highGk.setBasis (&Gk(ik));
                  highStates <<= waveBasis | highGk;

                  // --- orthonormalization
                  S.setOrthogonal (&highStates, tbStates);
                  S.orthonormalize (&highStates);

                  highGk = Gk(ik) | highStates;

                  // --- diagonalize high state subspace
                  SxDiracMat<TPrecCoeffG>::Eigensystem eigHigh;
                  eigHigh = SxDiracMat<TPrecCoeffG> (highGk.adjoint ()
                        ^ (hamPAW * highStates)).eigensystem();
                  highStates.rotate (eigHigh.vecs);
                  for (int i = 0; i < nHigh; ++i)
                     fermi.eps (nFromTB + i, iSpin, ik) = eigHigh.vals (i).re;

               }
            }  // LoopMPI
         }
      }
      fermi.eps.synMPI ();

      if (updateXWaves && maxStepsTB > 0)  {
         hamPAW.exchangePtr->wavesPtr = wavesPtr;
         hamPAW.exchangePtr->pBasisX = hamPAW.pBasis;
         hamPAW.exchangePtr->setFocc (fermi.focc);
      }

      if (maxStepsTB == 0) break;
      double oldEnergy = hamPtr->getEnergy (*wavesPtr, fermi);

      if (deltaE < dEnergy) break;
      fermi.fermiDistribution(ektDefault);
      if (hamPAW.exchangePtr) hamPAW.exchangePtr->setFocc (fermi.focc);
      mixer.addRhoIn (hamPtr->getRho());
      hamPtr->computeRho (fermi.focc, *wavesPtr);
      mixer.addRhoOut (hamPtr->getRho ());
      hamPtr->getRho () = mixer.getMixedRho ();
      double newEnergy = hamPtr->getEnergy (*wavesPtr, fermi);
      deltaE = fabs (newEnergy - oldEnergy);
      cout << iMixSteps << " Mixing Step: " << newEnergy 
           << ", Delta : " << deltaE 
           << ", Residuum : " << mixer.getNormR () << endl;

      SX_MPI_MASTER_ONLY
      {
         SxFileIO::appendToFile
         (   SxString(iMixSteps)                            + "\t" // iteration
           + SxString(GETTIME(Timer::LcaoTotal), "%12.8f")  + "\t" // time
           + SxString(newEnergy, "%15.12f")                 + "\t" // E_tot  [H]
           + SxString(deltaE,  "%15.12f")                   + "\n" // dEnergy [H]
         , "energy.dat");
         SxFileIO::appendToFile
         (   SxString (iMixSteps)                           + "\t" // iteration
           + SxString(mixer.getNormR())                            // |R|
           + ((nSpin == 1) ? SxString("\n")
                           : ("\t" + SxString(mixer.getNormRPol()) + "\n")) // |Rpol|
         , "residue.dat");
      }
   }
   energy = hamPtr->getEnergy () - ekt * fermi.getEntropy ();

   // Print Basis composition
   if (basisDecomp)  {
      FILE *filePtr;
      FILE *HFilePtr;
      FILE *SFilePtr;
      FILE *SEigFilePtr;
      FILE *SVecFilePtr;
      SxString fileName = "BasisDecomposition.dat";
      SxString HFile = "HamiltonianOverK.dat";
      SxString SFile = "OverlapOverK.dat";
      SxString SEigFile = "OverlapValuesOverK.dat";
      SxString SVecFile = "OverlapEigVecOverK.dat";
      if ( (filePtr = fopen(fileName.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << fileName.ascii() << endl;
         SX_EXIT;
      }
      if ( (HFilePtr = fopen(HFile.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << HFile.ascii() << endl;
         SX_EXIT;
      }
      if ( (SFilePtr = fopen(SFile.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << SFile.ascii() << endl;
         SX_EXIT;
      }
      if ( (SEigFilePtr = fopen(SEigFile.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << SEigFile.ascii() << endl;
         SX_EXIT;
      }
      if ( (SVecFilePtr = fopen(SVecFile.ascii(), "w")) == NULL)  {
         cout << "Cannot open " << SVecFile.ascii() << endl;
         SX_EXIT;
      }
      for (int ik = 0; ik < Gk.getNk (); ik++)  {
         fprintf(HFilePtr, "%i\t",ik);
         fprintf(SFilePtr, "%i\t",ik);
         fprintf(SEigFilePtr, "%i\t",ik);
         SxPAWBasis pawBasis (Gk(ik).getThis (), hamPAW.pBasis);
         int iSpin = 0;
         SxMatrix<Complex16> Hnm(nOrb,nOrb), Snm(nOrb,nOrb);
         // --- construct Hamiltonian Hmn and Overlap Smn
         // H = H.adjoint, S = S.adjoint ()
         for (int jState = 0; jState < nOrb; jState++)   {
            PsiG psiJ = initBasis.getAOinG(ik, jState);
            psiJ.handle->auxData.iSpin = iSpin;
            PsiG pawJ = pawBasis | psiJ;
            PsiG HPsiJ = hamPAW * pawJ,
                 SPsiJ = S * pawJ;

            for (int iState = jState; iState < nOrb; iState++)   {
               PsiG psiI = initBasis.getAOinG(ik,iState);
               Hnm(iState,jState) = dot(psiI,HPsiJ);
               Hnm(jState,iState) = Hnm(iState,jState).conj();
               Snm(iState,jState) = dot (psiI, SPsiJ );
               Snm(jState,iState) = Snm(iState,jState).conj();
               fprintf(HFilePtr, "%f\t",Hnm(iState,jState).absSqr());
               fprintf(SFilePtr, "%f\t",Snm(iState,jState).absSqr());
            }
         }
         fprintf(HFilePtr, "\n");
         fprintf(SFilePtr, "\n");
         SxMatrix<Complex16>::Eigensystem Seig;
         Seig = Snm.eigensystem ();
         for (int i = 0; i < nOrb; i++) {
            fprintf(SEigFilePtr, "%f\t",Seig.vals(i).re);
            fprintf(SVecFilePtr, "%i\t%i\t",ik,i);
            for (int j = 0; j < nOrb; j++)  {
               fprintf(SVecFilePtr, "%f\t",Seig.vecs.colRef(i)(j).re);
            }
            fprintf(SVecFilePtr, "\n");
         }
         fprintf(SEigFilePtr, "\n");
         fprintf(SVecFilePtr, "\n");
         
         // transform generalized eigenvalue problem to standard
         // eigenvalue problem via cholesky decomposition
         SxMatrix<Complex16> Lnm = Snm.choleskyDecomposition();
         SxMatrix<Complex16> invLnm = Lnm.inverse(); 
         // transform Hnm
         Hnm = invLnm ^ Hnm ^ invLnm.adjoint();

         // solve standard eigenvalue problem
         SxMatrix<Complex16>::Eigensystem eig;
         eig = Hnm.eigensystem ();

         int nEigenStates = (int)eig.vecs.row(0).getSize();
         eig.vecs = invLnm.adjoint() ^ eig.vecs;

         fprintf(filePtr,"%i kPoint \n\n",ik);
         for(int iState = 0; iState < nEigenStates; iState++)  {
            double sum = eig.vecs.colRef(iState).absSqr().sum();
            for(int iOrb = 0; iOrb < nOrb; iOrb++)  {
               PsiG Psi = initBasis.getAOinG(ik,iOrb);
               int is = Psi.handle->auxData.is;
               int iot = Psi.handle->auxData.n;
               int l = Psi.handle->auxData.l;
               fprintf(filePtr, "State = %i, Basisorbital = %i, Species = %i, iot = %i, l = %i, Coefficient = %f\n",
                                 iState,iOrb,is,iot,l,eig.vecs(iOrb, iState).absSqr()/sum);
            }
            fprintf(filePtr, "\n");
         }
         fprintf(filePtr,"\n-------------------------------------\n");
      }
      fclose (filePtr);
   }
}

void SxHamSolver::setNStates (const SxArray<int> &nPerK)
{
//   SX_EXIT;  // FIXME: to be tested
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   dynamic_cast<SxPW*>(wavesPtr.getPtr ())->setNStates (nPerK);
   fermi.resize (nPerK);
}

void SxHamSolver::setNStates (int nStatesNew)
{
   SxArray<int> nPerK(wavesPtr->getNk());
   if (nStatesNew > -1)  {
      nPerK.set (nStatesNew);
      setNStates (nPerK);
   }  else  {
      for (int ik = 0; ik < wavesPtr->getNk(); ik++)  {
         nPerK(ik) = wavesPtr->getNStates (ik);
      }
      setNStates (nPerK);
   }
}


void SxHamSolver::steepestDescent (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (wavesPtr);
   SX_NO_MPI;
   SxPWSet &waves = *wavesPtr;
   SxGkBasis &Gk = waves.getGkBasis();
   SX_CLOCK (Timer::ElMinim);
   SX_CHECK (calc == true);

   int printSteps = 1;

   try  {
      ekt     = (cmd->contains("ekt")) 
              ? cmd->get("ekt")->toReal() / HA2EV 
              : ektDefault;
      keepRho = (cmd->contains("keepRhoFixed")) 
              ? cmd->get("keepRhoFixed")->toAttribute ()
              : false;
      keepOcc = (cmd->contains("keepOccFixed")) 
              ? cmd->get("keepOccFixed")->toAttribute ()
              : false;
      deltaT   = cmd->get("deltaT")->toReal() / 33.3333333; // eMass
      maxSteps = cmd->get("maxSteps")->toInt ();
      if (dEnergyLow) {
         dEnergy = (cmd->contains("dEnergyLow"))
            ? cmd->get("dEnergyLow")->toReal () 
            : 1e-8;
      } else dEnergy  = cmd->get("dEnergy")->toReal ();
      dPsiConv = cmd->get("dPsi")->toReal();
      printSteps = cmd->get("printSteps")->toInt();
      if (cmd->contains("spinMoment"))  {
         fermi.setSpinMoment (cmd->get("spinMoment")->toReal (), true);
         fermi.fermiDistribution (ekt);
      }
      if (cmd->contains("keepSpinFixed"))  {
         fermi.setSpinMoment (fermi.getSpinMoment (), 
                              cmd->get("keepSpinFixed")->toAttribute ());
      }
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   cout << SX_SEPARATOR;
   cout << "| Steepest Descent\n";
   cout << SX_SEPARATOR;
   cout << "|   deltaT:                 " << 33.3333333 * deltaT << endl;
   cout << "|   max. steps:             " << maxSteps << endl;
   cout << "|   dEnergy:                " << dEnergy  << " H\n";
   cout << "|   dPsi:                   " << dPsiConv << endl;
   cout << "|   details printed every:  " << printSteps  << " iteration(s)\n";
   if (keepOcc)  cout << "|   occupation:             " << "fixed\n";
   else  {
      cout << "|   ekt:                    " << ekt * HA2EV << " eV" << endl;
      if (ekt * HA2EV > 1e-5)  {
         cout << "|   smearing scheme:        ";
         fermi.printSmearing ();
      }
   }
   if (fermi.keepSpinMoment)
      cout << "|   spin moment:       " << fermi.spinMoment << " (fixed)\n";
   cout << "|   type of calculation:    ";
   if (keepRho)  cout << "band structure run\n";
   else          cout << "SCF minimization\n";
   cout << SX_SEPARATOR;

   if (!calc)  return;


   double foccMix = 0.05;
   //double foccMix = 1.;

   int ik, nk         = waves.getNk ();
   int iSpin;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   int i;
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   int nPlWaves;
   SxDiracVec<TPrecCoeffG> psiG, g, X, speed;
   double eps, eKinC, eBand, lastEnergy = 1e50;  // just a large value
   double UNIT = printHartree ? 1. : HA2EV;
   for (ik=0, nPlWaves=0; ik < nk; ik++)  nPlWaves += Gk(ik).ng;

   std::streamsize oldPrec = cout.precision ();
   cout.precision(10);
   for (int it=0; it < maxSteps; it++)  {
      eKinC = 0.;
      for (ik=0; ik < nk; ik++)  {
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            for (i=0; i < nStates; i++)  {

               psiG    = waves(i,iSpin,ik);
               speed.copy (psiG);

               // TODO: ugly
               psiG.handle->auxData.i     = i;
               psiG.handle->auxData.iSpin = iSpin;
               psiG.handle->auxData.ik    = ik;
               
               g    = *hamPtr | psiG;
               eps  = (psiG ^ g).chop();
               X    = deltaT * (g - psiG*eps);


               psiG -= X;
               fermi.eps(i,iSpin,ik) = eps;

               // --- compute change of the wavefunction coefficients
               speed -= psiG;
               eKinC += (speed ^ speed).chop().re;
            }
            SxPWOverlap ().orthonormalize (&waves(iSpin,ik), GramSchmidt);
         }
      }
      eKinC = 1000. * sqrt (eKinC / (nStates * nPlWaves * deltaT*deltaT) );

      if (!keepOcc)  fermi.fermiDistribution (ekt, foccMix);

      if (!keepRho)  {
         hamPtr->computeRho (fermi.focc, waves);
      }
      // TODO: should be compute(...)
      hamPtr->getEnergy (waves, fermi);

      // --- print results
      eBand  = fermi.getEBand (keepRho ? SxFermi::NoFocc 
                                       : SxFermi::UseFocc);
      energy = (keepRho) ? eBand : hamPtr->getEnergy ();
      sxprintf ("e%s(%d)=%15.12f, dPsi=%g\n", 
              (keepRho) ? "Band" : "Tot", 
               it+1, energy*UNIT, eKinC); fflush (stdout);
      double entropy = fermi.getEntropy ();
      printEnergyDatLine (it, energy, energy - ekt * entropy, eBand, entropy);

      // --- convergence?
      if (   fabs (energy - lastEnergy) < dEnergy   // energy convergence
          && eKinC < dPsiConv )                     // wavefunction conv.
      {
         sxprintf ("Convergence reached.\n");
         break;
      }
      if ( it+1 >= maxSteps )  {
         sxprintf ("WARNING: Maximum number of steps exceeded.\n");
         sxprintf ("         Convergence not yet reached.\n");
         break;
      }

      if ( !(it % printSteps) )  {
         hamPtr->printEnergies (); fermi.printOccupation ();
      }


      lastEnergy = energy;

   }
   cout.precision (oldPrec);
   // final printout
   hamPtr->printEnergies ();
   fermi.printOccupation (true);
   energy = hamPtr->getEnergy ();
}

SxVector<Double> SxHamSolver::lineFit4th (double f0, double fp0,
                                          double x1, double f1, double fp1,
                                          double x2, double f2)
{
   // Find 4th-order polynomial
   // f(0)   = f0   = a0
   // f'(0)  = fp0  = a1
   // f(x1)  = f1   = a0 + a1 x1 +   a2 x1^2 +   a3 x1^3 +   a4 x1^4
   // f'(x1) = fp1  =      a1    + 2 a2 x1   + 3 a3 x1^2 + 4 a4 x1^3
   // f(x2)  = f2   = a0 + a1 x2 +   a2 x2^2 +   a3 x2^3 +   a4 x2^4

   // vector contains a0..a4
   SxVector<Double> res(5);
   res(0) = f0;
   res(1) = fp0;
   SxMatrix<Double> m(3,3);
   SxVector<Double> rhs(3);
   m(0,0) = 1.;
   m(1,0) = 2.;
   m(2,0) = 1.;
   rhs(0) = (f1 - f0 - fp0 * x1) / (x1 * x1);
   m(0,1) = x1;
   m(1,1) = 3. * x1;
   m(2,1) = x2;
   rhs(1) = (fp1 - fp0) / x1;
   m(0,2) = x1 * x1;
   m(1,2) = 4. * x1 * x1;
   m(2,2) = x2 * x2;
   rhs(2) = (f2 - f0 - fp0 * x2) / (x2 * x2);
   res(SxIdx(2,4)) = m.inverse () ^ rhs;
   return res;
}

void SxHamSolver::allStateCG (const SxSymbolTable *cmd, bool calc)
{
   SxHamiltonian &H = *hamPtr;
   SxPWSet &waves = *wavesPtr;
   SxGkBasis &Gk = waves.getGkBasis();
   SxPtr<SxGkBasis> gkPtr = waves.getGkBasisPtr();

   SX_CLOCK (Timer::ElMinim);
   bool testLineMinim = false;
   bool keepOccFixed = false;
   int  printSteps=1, minSteps=2;
   double kappa = 0.1;
   bool fixKappa = false;
   bool finalDiag = false, initialDiag = true;
   double lMinAvg = 1.;
   double ltsFac = 0.25;
   //bool calcForces = false; // TODO: need interface in SxHamiltonian
   try  {
      maxSteps = cmd->contains ("maxSteps")
               ? cmd->get("maxSteps")->toInt()
               : 100;
      if (dEnergyLow) {
         dEnergy = (cmd->contains("dEnergyLow"))
            ? cmd->get("dEnergyLow")->toReal ()
            : 1e-8;
      }
      else {
         dEnergy  = cmd->contains ("dEnergy")
            ? cmd->get("dEnergy")->toReal()
            : 1e-8;
      }
      printSteps = cmd->contains ("printSteps")
                 ? cmd->get("printSteps")->toInt()
                 : 10;
      if (cmd->contains("keepOccFixed"))
         keepOccFixed = cmd->get("keepOccFixed")->toAttribute();
      ekt       = (cmd->contains("ekt")) 
                ? cmd->get("ekt")->toReal() / HA2EV 
                : ektDefault;
      if (cmd->contains("kappa")) kappa = cmd->get("kappa")->toReal ();
      fixKappa = kappa < 0.;
      kappa = fabs(kappa);
            
      H.setXCMeshDensity (1);
      if (cmd->contains("xcMeshDensity"))  {
         if (!dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ())) {
            cout << "Finer mesh for XC contributions yet only implemented for PW Basis! " << endl; SX_EXIT;
         } else {
            H.setXCMeshDensity (cmd->get("xcMeshDensity")->toInt ());
         }
      }
      if (cmd->contains ("initialDiag"))  {
         initialDiag = cmd->get("initialDiag")->toAttribute ();
      }
      if (cmd->contains ("finalDiag"))  {
         finalDiag = cmd->get("finalDiag")->toAttribute ();
      }
      if (cmd->contains ("testLineMinim"))  {
         testLineMinim = cmd->get ("testLineMinim")->toAttribute ();
      }
      if (cmd->contains ("testLineStep"))  {
         ltsFac = cmd->get ("testLineStep")->toReal ();
      }
      if (cmd->contains ("trialStep0"))  {
         lMinAvg = cmd->get ("trialStep0")->toReal () / 1.5;
      }
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   bool changedHam = hamPtr->rereadTable (cmd);

   bool isMetal = false;
   double s = 3. - nSpin;
   // --- check that all states are fully occupied
   int i, iSpin, ik, nk = waves.getNk();
   SX_CHECK(nStates == waves.getNStates(), nStates, waves.getNStates());
   SX_CHECK(nSpin == waves.getNSpin (), nSpin, waves.getNSpin ()) ;

   for (ik = 0; ik < nk && !isMetal; ik++)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork (ik))  {
         for (iSpin = 0; iSpin < nSpin && !isMetal; iSpin++)  {
            for (i=0; i < nStates && !isMetal; i++)  {
               if ( fabs(fermi.focc(i,iSpin,ik) - s) > 1e-5 )  {
                  isMetal = true;
               }
            }
         }
      }
   }
   isMetal = SxLoopMPI::lor(isMetal);

   cout << SX_SEPARATOR;
   cout << "| All-state conjugate gradient\n";
   cout << SX_SEPARATOR;
   cout << "|   max. steps:             " << maxSteps      << endl;
   cout << "|   min. steps:             " << minSteps      << endl;
   cout << "|   dEnergy:                " << dEnergy       << " H\n";
   cout << "|   ekt:                    " << (ekt * HA2EV) << " eV\n";
   if (HA2EV * ekt > 1e-5 && isMetal)  {
      cout << "|   smearing scheme:        ";
      fermi.printSmearing ();
   }
   cout << "|   details printed every:  " << printSteps    << " iteration(s)\n";
   if (keepOccFixed)
      cout << "|   occupation numbers:     " << "fixed\n";
   else
      cout << "|   type of calculation:    " << "SCF run\n";
   if (isMetal)  {
      cout << "|   kappa (eps precond.):   " << kappa;
      if (fixKappa) cout << " (fixed)" << endl;
      else          cout << " (auto-adjust)" << endl;
   }
   cout << SX_SEPARATOR;

   if (!calc)  return;

   // --- check that we have more than 1 Bloch state (current vector class
   //     cannot handle 1x1 matrices correctly)  // TODO: XPress
   if (waves.getNStates () == 1)  {
      sxprintf ("The current implementation of the CCG scheme cannot handle\n"
              "systems with only 1 state.\n"
              "Please provide an empty state.\n");
      SX_EXIT;
   }
   if (initialDiag)  {
      // perform an iterative diagonalization with current Hamiltonian
      dRelEps=0.1;
      blockStateByStateCG (5, 64, false, SCFrun, 0, false);
   }

   PsiG  K;
   // transfer StorageModel from waves
   bool keepWavesOnDisk = waves.wavesOnDisk ();
   SxPtr<SxPW> dPsiPtr;
   SxPtr<SxPWSet> psiTrialPtr,
                  XPtr = waves.getNew ();
   
   SxPWSet &X        = *XPtr;

   SxFermi fermiTrial (fermi);
   SxArray<SxArray<SxDiracMat<Complex16> > > dEta, xEta;
   double trace, traceOld = 1., traceR = 0., traceEps = 0., traceEpsOld = 0.;
   double sigma = 0., sigmaEps = 0., gamma = 0., curv = 0.;
   double lambdaMin = 0., lambdaTMin, lambdaTMax, lambdaT=1., lambdaTStep = 1.;
   lambdaTMin = lambdaTMax = lambdaT;

   double dMu, dNElec, fermidos;
   double eTot, eTrial, eFit;
   double freeEnergy = 0., fTrial = 0., fOld = 0.;

   SxOrthoMethod orthoMethod = isMetal ? GramSchmidt : Loewdin;
   SxOverlap S = H.getS ();
   SxString filename;
   
   if (isMetal && ekt < 1e-5) ekt = 1e-5;
   if (isMetal && !initialDiag)  {
      for (ik = 0; ik < nk; ++ik)  {
         SX_MPI_LEVEL("waves-k");
         if (SxLoopMPI::myWork (ik))  {
            for (iSpin = 0; iSpin < nSpin; ++iSpin)
               fermi.eps(iSpin,ik) = subspaceDiagonalization(&waves(iSpin,ik));
         }
      }
      // put SX_MPI_SOURCE into SxBundle3::synMPI
      fermi.eps.synMPI ();
   }
   if (!keepOccFixed) fermi.fermiDistribution (ekt);
   for (ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork (ik))  {
         for (iSpin = 0; iSpin < nSpin; ++iSpin)
            S.orthonormalize (&waves(iSpin,ik), orthoMethod);
      }
   }
   H.computeRho (fermi.focc, waves);

   // --- check PAW exchange projector basis
   {
      SxPAWHamiltonian* pawHam = dynamic_cast<SxPAWHamiltonian*>(&*hamPtr);
      if (pawHam
          && pawHam->exchangePtr
          && pawHam->exchangePtr->pBasisX.getPtr () != pawHam->pBasis.getPtr ())
      {
         cout << "Updating PAW exchange projector basis" << endl;
         pawHam->exchangePtr->pBasisX =  pawHam->pBasis;
      }
   }

   eTot = H.getEnergy (waves, fermi);
   freeEnergy = eTot - ekt * fermi.getEntropy ();
   cout << SX_SEPARATOR;
   cout << "Etot after TB Initialization" << endl;
   cout << "Etot = " << eTot << ", free energy = " << freeEnergy << endl;
   cout << SX_SEPARATOR;
   if (isMetal)  {
      dEta.resize (nk);
      xEta.resize (nk);
      for (ik = 0; ik < nk; ++ik)  {
         dEta(ik).resize (nSpin);
         xEta(ik).resize (nSpin);
      }
   }

   SxDiracMat<Complex16> ham;
   SxDiracMat<Complex16>::Eigensystem eig;
   double xEtaFactor = 1.;
   int slowPoints = 0, slowLimit = 10;
   int lastRestart = 0;

   bool checkGradient = false, quadMinOK = false;
   bool etaOnly = false;

   for (int it=0; it < maxSteps; it++)  {
      if (!isMetal)  {
         for (ik = 0; ik < nk; ik++)
            for (iSpin = 0; iSpin < nSpin; iSpin++)
               fermi.eps(iSpin,ik).set (0.);  // no eigenvalues available
         if (!keepOccFixed)  fermi.fermiDistribution (ekt);
      }

      trace = 0.;
      traceR = 0.;
      dNElec = 0.; fermidos = 0.;
      if (!dPsiPtr)
         dPsiPtr = SxPtr<SxPW>::create (nStates,nSpin,gkPtr,
                                        keepWavesOnDisk ? "." : "");
      for (ik = 0; ik < nk; ik++)  {
         SX_MPI_LEVEL("waves-k");
         if (!SxLoopMPI::myWork (ik)) continue;

         const SxGBasis &gk = waves.getGkBasis()(ik);
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            
            PsiG psiGI = gk | waves(iSpin, ik);
            PsiG dPsi  = (*dPsiPtr)(iSpin,ik);

            dPsi <<= H | waves(iSpin, ik);
            if (isMetal)  {
               dEta(ik)(iSpin) = (psiGI.adjoint () ^ dPsi);
               dPsi -= S | ((waves(iSpin, ik) ^ dEta(ik)(iSpin))
                           .setAux(waves(iSpin,ik).handle->auxData)); // TODO: ugly
               for (i = 0; i < nStates; ++i)  {
                  dEta(ik)(iSpin)(i,i) -= fermi.eps(i,iSpin,ik);
               }

               dEta(ik)(iSpin) *= -1.;
            }
            else {
               dPsi -= S | ((waves(iSpin,ik) ^ (psiGI.adjoint() ^ dPsi))
                           .setAux(waves(iSpin,ik).handle->auxData));
            }

            // --- compute trace and alpha
            for (i=0; i < nStates; i++)  {
               //trace += 2. * dPsiSet(i,iSpin,ik).normSqr () * Gk.weights(ik)
               //         * fermi.focc(i,iSpin,ik);
               //K = H.preconditioner (dPsiSet(i,iSpin,ik), SxPWHamiltonian::Arias);
               K = H.preconditioner (gk | waves(i,iSpin,ik));
               trace += 2.  * Gk.weights(ik) * fermi.focc(i,iSpin,ik)
                  * dot(dPsi.colRef(i), K * dPsi.colRef(i)).re;
               if (isMetal)  {
                  // --- contribution from subspace rotation
                  for (int j = 0; j < i; ++j)  {
                     double deltaEps = fermi.eps(i,iSpin,ik)
                        - fermi.eps(j,iSpin,ik);
                     double deltaFocc = fermi.focc(i,iSpin,ik)
                        - fermi.focc(j,iSpin,ik);
                     if (fabs(deltaEps) > 1e-10)
                        traceR -= 2. * Gk.weights(ik) * deltaFocc / deltaEps
                           * dEta(ik)(iSpin)(j,i).absSqr ();
                  }
                  // --- derivative of Fermi energy along dEta
                  double ff = fermi.dFoccFermi(i,iSpin,ik)
                            * Gk.weights(ik);
                  fermidos += ff;
                  dNElec += ff * dEta(ik)(iSpin)(i,i).re;
               }
            }
         }
      }

   
      {
         SX_MPI_SOURCE("waves-k", TaskGroupMaster);
         SX_MPI_TARGET(TopLevel, TaskGroupAll);
         // SX_MPI_SUM_VARS(trace,traceR,dNElec,fermidos);
         trace    = SxLoopMPI::sum (trace);
         traceR   = SxLoopMPI::sum (traceR);
         dNElec   = SxLoopMPI::sum (dNElec);
         fermidos = SxLoopMPI::sum (fermidos);
      }

      if (isMetal)  {
         dMu = (fabs(fermidos) > 1e-16) ? -dNElec/fermidos : 0.;
         double traceF = 0.;
         // --- contribution from occupation changes
         if (!keepOccFixed)  {
            for (ik = 0; ik < nk; ik++)  {
               SX_MPI_LEVEL("waves-k");
               if (!SxLoopMPI::myWork (ik)) continue;
               for (iSpin = 0; iSpin < nSpin; iSpin++)  {
                  for (i=0; i < nStates; i++)  {
                     traceF += Gk.weights(ik)
                        * fermi.dFoccFermi(i,iSpin,ik)
                        * (dEta(ik)(iSpin)(i,i).re + dMu)
                        * dEta(ik)(iSpin)(i,i).re;
                  }
               }
            }
            SX_MPI_SOURCE("waves-k", TaskGroupMaster);
            SX_MPI_TARGET(TopLevel, TaskGroupAll);
            traceF = SxLoopMPI::sum (traceF);
         } else {
            dMu = 0.;
         }
         cout << "trace= " << trace
              << "; traceR= " << traceR
              << "; traceF= " << traceF << endl;

         traceEpsOld = traceEps;
         traceEps = kappa * (traceR + traceF);
         if (traceEps > trace) slowPoints++;
         trace += traceEps;

      }
      // unload waves - currently not needed
      waves.memMinimize ();

      if (checkGradient)  {
         checkGradient = false;
         double sigOld = 0.;
         fermidos=dNElec=0.;
         double sigOldE = 0.;
         for (ik = 0; ik < nk; ik++)  {
            SX_MPI_LEVEL("waves-k");
            if (!SxLoopMPI::myWork (ik)) continue;
            const SxGBasis &gk = waves.getGkBasis()(ik);
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               PsiG dPsi = (*dPsiPtr)(iSpin,ik);
               for (i=0; i < nStates; i++)  {
                  sigOld -= ( 2. * Gk.weights(ik) * fermi.focc(i,iSpin,ik) )
                          * dot(gk | X(i,iSpin,ik), dPsi.colRef(i) ).re;
                  if (isMetal)  {
                     // --- derivative from subspace rotation
                     for (int j = 0; j < i; ++j)  {
                        double deltaEps = fermi.eps(i,iSpin,ik)
                                        - fermi.eps(j,iSpin,ik);
                        double deltaFocc = fermi.focc(i,iSpin,ik)
                                         - fermi.focc(j,iSpin,ik);
                        if (fabs(deltaEps) > 1e-10)  {
                           sigOldE += 2. * Gk.weights(ik) * deltaFocc / deltaEps
                                     * (dEta(ik)(iSpin)(j,i)
                                     *  xEta(ik)(iSpin)(i,j)).re
                                     * xEtaFactor;
                        }
                     }
                     // --- derivative of Fermi energy
                     double ff = fermi.dFoccFermi(i,iSpin,ik) * Gk.weights(ik);
                     fermidos += ff;
                     dNElec += ff * xEta(ik)(iSpin)(i,i).re;
                  }
               }
            }
         }

         {
            SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
            SX_MPI_TARGET (TopLevel, TaskGroupAll);
            sigOld   = SxLoopMPI::sum (sigOld);
            if (isMetal) {
               // SX_MPI_SUM_VARS(fermidos,dNelec);
               fermidos = SxLoopMPI::sum (fermidos);
               dNElec   = SxLoopMPI::sum (dNElec);
            }
         }
         dMu = (fabs(fermidos) > 1e-16) ? -dNElec/fermidos : 0.;
         //cout << "dMu sig = " << dMu << endl;

         if ((isMetal) && (!keepOccFixed))  {
            // --- derivative from occupation changes
            for (ik = 0; ik < nk; ik++)  {
               SX_MPI_LEVEL("waves-k");
               if (!SxLoopMPI::myWork (ik)) continue;
               for (iSpin = 0; iSpin < nSpin; iSpin++)  {
                  for (i=0; i < nStates; i++)  {
                     sigOldE -= Gk.weights(ik)
                              * fermi.dFoccFermi(i,iSpin,ik)
                              * (xEta(ik)(iSpin)(i,i).re + dMu)
                              * dEta(ik)(iSpin)(i,i).re
                              * xEtaFactor;
                  }
               }
            }
         } else {
            dMu = 0.;
         }
         if (isMetal) {
            SX_MPI_SOURCE("waves-k", TaskGroupMaster);
            SX_MPI_TARGET(TopLevel, TaskGroupAll); 
            sigOldE  = SxLoopMPI::sum (sigOldE);
         }

         if (!fixKappa && quadMinOK && xEtaFactor > 0.9) {
            // --- check gradient of F_min / d kappa
            double dsdk = sigmaEps / kappa;
            double dedk = sigOldE / kappa;
            double dFdkappa = -2. * sigma * dsdk * (freeEnergy - fOld)
                            + sqr (sigma) * (dedk + dsdk * lambdaMin);
            dFdkappa /= sqr (2. * curv * lambdaMin);
            //cout << "dF/dkappa=" << dFdkappa << endl;
            // dF / d (ln kappa) / (dF/d lambda)
            double dfnorm = kappa * dFdkappa / trace;
            //cout << "dF/d(ln kappa) / dF/d eta=" << dfnorm << endl;
            double kappaOld=kappa;
            if (fOld - freeEnergy > 1e-10
                && fabs(trace) < 0.9 * fabs(traceOld) /* gamma < 0.9 */ )
            {
               if (dfnorm < -5.)
                  kappa *= 5.;
               else if (dfnorm < -1.)
                  kappa *= 2.;
               else if (dfnorm < -0.5)
                  kappa *= 1.5;
               else if (dfnorm < -0.1)
                  kappa *= 1.1;
               if (dfnorm > 5.)
                  kappa /= 5.;
               else if (dfnorm > 1.)
                  kappa /= 2.;
               else if (dfnorm > 0.5)
                  kappa /= 1.5;
               else if (dfnorm > 0.1)
                  kappa /= 1.1;
            }
            cout << "kappa=" << kappa << endl;
            // --- update traces
            trace -= traceEps;
            traceEps *= kappa/kappaOld;
            trace += traceEps;
            traceOld -= (1. - kappa/kappaOld) * traceEpsOld;
         }

         sigOld += sigOldE;
         cout << "Check gradient: reduced to " << (sigOld / sigma) << endl;
         SxVector<Double> coeff = lineFit4th (freeEnergy, sigOld,
                                              -lambdaMin, fOld, sigma,
                                              lambdaT - lambdaMin, fTrial);
         if (testLineMinim)  {
            filename = "lineminim";
            filename += SxString(".it=") + (it-1) + SxString(".dat");
            FILE *fp = fopen (filename.ascii(), "a");
            fprintf (fp,"&\n");
            for (double l = lambdaTMin; l <= lambdaTMax; l+= lambdaTStep)  {
               double x = l - lambdaMin;
               fprintf (fp, "%16.12f\t%16.12f\t%16.12f\n", l,
                        coeff(0) + x * (coeff(1) + x * (coeff(2) + x * (
                        coeff(3) + x * coeff(4)))),
                        freeEnergy + x * sigOld);
            }
            fclose (fp);
         }
         if (fabs(sigOld) > 0.15 * fabs(sigma))  {
            // --- find 4-th order local minimum and jump to it

            // unload dPsi - no longer needed
            dPsiPtr->memMinimize ();

            // --- Newton zero search for gradient m
            double x, m, mp;
            x = (coeff(2) > 0.) ? 0. : -lambdaMin;
            for (int itNewton = 0; itNewton < 100; itNewton++)  {
               m = coeff(1) + x * (2. * coeff(2) + x * (3. * coeff(3)
                          + x * 4. * coeff(4)));
               mp = 2. * coeff(2) + x * (6. * coeff(3) + x * 12. * coeff(4));
               if (fabs(m) < 1e-9 * fabs(sigma)) break;
               if (mp <= 0.) break;
               x -= m / mp;
            }
            double freeGuess = coeff(0) + x * (coeff(1) + x * (coeff(2) +
                               x * (coeff(3) + x * coeff(4))));

            if (mp > 0.)  {
               cout << "4-th order step: " << x << endl;
               // --- jump to 4th-order minimum
               for (ik = 0; ik < nk; ik++)  {
                  SX_MPI_LEVEL("waves-k");
                  if (!SxLoopMPI::myWork (ik)) continue;
                  for (iSpin = 0; iSpin < nSpin; iSpin++)  {
                     waves(iSpin,ik) -= (x * X(iSpin,ik));
                     S.orthonormalize (&waves(iSpin,ik), orthoMethod);
                     if (isMetal)  {
                        ham = (-x * xEtaFactor) * xEta(ik)(iSpin);
                        for (i = 0; i < nStates; ++i)
                           ham(i,i) += fermi.eps(i,iSpin,ik);
                        eig = ham.eigensystem ();

                        fermi.eps(iSpin, ik) = eig.vals.real ();
                        waves(iSpin, ik).rotate (eig.vecs);
                        X(iSpin, ik).rotate (eig.vecs);
                        xEta(ik)(iSpin) = eig.vecs.adjoint () ^ xEta(ik)(iSpin)
                                          ^ eig.vecs;
                     }
                  }
               }
               fermi.eps.synMPI ();
               //syncWaves ();

               X.memMinimize (); // X currently not needed
            
               // --- SCF update
               //double fermiOld = fermi.eFermi;
               if (!keepOccFixed) fermi.fermiDistribution (ekt);
               //double deltaMu = fabs((fermi.eFermi - fermiOld) 
               //                      - dMu * x * xEtaFactor);
               //cout << "deltaMu=" << deltaMu << endl;
               H.computeRho (fermi.focc, waves);
               //if (calcForces) H.calcForces = true;
               eTot = H.getEnergy (waves, fermi);
               /*
               if (calcForces)  {
                  H.calcForces = false;
                  cout << "Forces="; 
                  for (int ia = 0; ia < H.fTotal.getNAtoms (); ++ia)  {
                     for (int d = 0; d < 3; ++d)
                        sxprintf ("\t%.16f", H.fTotal(ia)(d));
                  }
                  cout << endl;
               }
               */
               freeEnergy = eTot - ekt * fermi.getEntropy ();
               cout << "error in 4-th order refinement: "
                    << (freeEnergy - freeGuess) / fabs(fOld - freeEnergy)
                    << endl;
               //sxprintf ("F = %.12f\n", freeEnergy);
               sxprintf ("\n e(it=%d)=%20.18f H, sigma=%g\n",
                       it, eTot, sigOld);
               fflush (stdout);
               printEnergyDatLine (it, eTot, freeEnergy,
                                   fermi.getEBand (SxFermi::UseFocc),
                                   fermi.getEntropy ());

               // go to next iteration
               continue;
            } else {
               cout << "4-th order refinement skipped: "
                    << "(d^2 F)/(d lambda)^2 = " << mp << endl;
            }
         }
      }

      //double goodGamma = 0.5;
      double goodGamma = 0.4;
      if (it > 0)  gamma = trace / traceOld;
      //if (eTot) gamma = 0.;
      traceOld = trace;
      if (gamma > 1.)
         slowPoints += 3;
      else if (gamma > 1.7 * goodGamma)
         slowPoints += 2;
      else if (gamma > 1.4 * goodGamma)
         slowPoints++;
      else if (gamma < goodGamma)  {
         slowPoints = 0;
         if (slowLimit > 10) slowLimit--;
      }
      if (slowPoints > slowLimit && fOld > freeEnergy)  {
         cout << "Slow convergence detected. Trying steepest descent" << endl;
         gamma = 0.;
         slowPoints = 0;
         slowLimit += 5;
      }
      //if (it % 20 == 0) gamma = 0.;
      //gamma = 0.;

      sigma = 0.; sigmaEps = 0.;
      dNElec = 0.; fermidos = 0.;
      // --- conjugate direction (search vector X)
      for (ik = 0; ik < nk; ik++)  {
         SX_MPI_LEVEL("waves-k");
         if (!SxLoopMPI::myWork (ik)) continue;
         const SxGBasis &gk = waves.getGkBasis()(ik);
         const SxBasis &waveBasis = waves.getBasis (ik);
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            PsiG dPsi = (*dPsiPtr)(iSpin,ik);
            if (isMetal)  {
               if ( it == 0 )  {
                  xEta(ik)(iSpin) = kappa * dEta(ik)(iSpin);
               } else {
                  xEta(ik)(iSpin) = kappa * dEta(ik)(iSpin) 
                                  + gamma * xEta(ik)(iSpin);
               }
               //xEta(ik)(iSpin).set (0.);
            }
            {
               // --- get preconditioned search direction
               SxDiracMat<TPrecCoeffG> KdPsi;
               KdPsi.copy (dPsi);
               for (i=0; i < nStates; i++)
                  KdPsi.colRef(i) *= H.preconditioner (gk | waves(i,iSpin,ik));

               if ( it == 0 )  {
                // 1st step: optimized steepest-descent
                  X(iSpin,ik) <<= waveBasis | KdPsi;
               }  else  {
                  //X(iSpin, ik) <<= (waveBasis | KdPsi) + gamma * X(iSpin, ik);
                  X(iSpin, ik) *= gamma;
                  X(iSpin, ik) += waveBasis | KdPsi;
               }
            }
            // search direction must be orthogonal to current waves
            S.setOrthogonal (&X(iSpin,ik), waves(iSpin,ik));

            // --- compute gradients along search direction
            for (i=0; i < nStates; i++)  {
               // --- compute directional derivative of E at psi, see Ref1, fig4
               if (Gk.weights(ik) * fabs(fermi.focc(i,iSpin,ik)) > 1e-16)  {
                  sigma -= ( 2. * Gk.weights(ik) * fermi.focc(i,iSpin,ik) )
                         * dot(gk | X(i,iSpin,ik), dPsi.colRef(i) ).re;
               }
               if (isMetal)  {
                  // --- derivative from subspace rotation
                  for (int j = 0; j < i; ++j)  {
                     double deltaEps = fermi.eps(i,iSpin,ik)
                                     - fermi.eps(j,iSpin,ik);
                     double deltaFocc = fermi.focc(i,iSpin,ik)
                                      - fermi.focc(j,iSpin,ik);
                     if (fabs(deltaEps) > 1e-10)  {
                        sigmaEps += 2. * Gk.weights(ik) * deltaFocc / deltaEps
                                    * (dEta(ik)(iSpin)(j,i)
                                    *  xEta(ik)(iSpin)(i,j)).re;
                     } else {
                        cout << "@";
                     }
                  }
                  // --- derivative of Fermi energy
                  double ff = fermi.dFoccFermi(i,iSpin,ik) * Gk.weights(ik);
                  fermidos += ff;
                  dNElec += ff * xEta(ik)(iSpin)(i,i).re;
               }
            }
         }
      }

      {
         SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
         SX_MPI_TARGET (TopLevel, TaskGroupAll);
         // SX_MPI_SUM_VARS (sigma, sigmaEps,dNElec,fermidos);
         sigma    = SxLoopMPI::sum (sigma);
         sigmaEps = SxLoopMPI::sum (sigmaEps);
         dNElec   = SxLoopMPI::sum (dNElec);
         fermidos = SxLoopMPI::sum (fermidos);
      }
      dMu = (fabs(fermidos) > 1e-16) ? -dNElec/fermidos : 0.;

      if ((isMetal) && (!keepOccFixed))  {
         double sigmaF = 0.;
         // --- derivative from occupation changes
         for (ik = 0; ik < nk; ik++)  {
            SX_MPI_LEVEL("waves-k");
            if (!SxLoopMPI::myWork (ik)) continue;
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               for (i=0; i < nStates; i++)  {
                  sigmaF -= Gk.weights(ik)
                          * fermi.dFoccFermi(i,iSpin,ik)
                          * (xEta(ik)(iSpin)(i,i).re + dMu)
                          * dEta(ik)(iSpin)(i,i).re;
               }
            }
         }
         SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
         SX_MPI_TARGET (TopLevel, TaskGroupAll);
         sigmaF = SxLoopMPI::sum (sigmaF);
         cout << "sigma=" << sigma
              << "; sigmaR=" << sigmaEps
              << "; sigmaF=" << sigmaF << endl;
         sigmaEps += sigmaF;
           
      } else {
         dMu = 0.;
      }

      // unload dPsi - no longer needed
      dPsiPtr->memMinimize ();
      if (!keepWavesOnDisk) dPsiPtr = SxPtr<SxPW> (); // destroy
      // do not destroy waves kept on disk, because when recreating,
      // the complete scratch file would be rewritten with random numbers

      //lambdaT = 1.0;
      lambdaT = 1.5 * lMinAvg;
      xEtaFactor = 1.;
      double xNorm = 0., eNorm = 0.;
      if ((isMetal))  {
         // --- find xEta trial step within linear regime
         for (int iTry = 0; iTry < 20; ++iTry)  {
            bool smallChanges = true;
            for (ik = 0; ik < nk; ik++)  {
               SX_MPI_LEVEL("waves-k");
               if (!SxLoopMPI::myWork (ik)) continue;
               for (iSpin = 0; iSpin < nSpin; iSpin++)  {
                  xNorm += Gk.weights(ik) * xEta(ik)(iSpin).normSqr ();
                  eNorm += Gk.weights(ik) * fermi.eps(iSpin,ik).normSqr ();
                  ham = (-lambdaT * xEtaFactor) * xEta(ik)(iSpin);
                  for (i = 0; i < nStates; ++i)
                     ham(i,i) += fermi.eps(i,iSpin,ik);
                  eig = ham.eigensystem ();
                  fermiTrial.eps(iSpin, ik) = eig.vals.real ();
                  if (keepOccFixed)  {
                     for (i = 0; i < nStates; ++i)  {
                        double smallEps = 1e-2 * fabs(xEta(ik)(iSpin)(i,i).re);
                        double deltaEps = fermiTrial.eps(i,iSpin,ik)
                                        + lambdaT * xEtaFactor * xEta(ik)(iSpin)(i,i).re
                                        - fermi.eps(i,iSpin,ik);
                        if (smallEps > 1e-6 && fabs(deltaEps) > smallEps)
                           smallChanges = false;
                     }
                  }
               }
            }
            {
               SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
               SX_MPI_TARGET (TopLevel, TaskGroupAll);
               xNorm = SxLoopMPI::sum (xNorm);
               eNorm = SxLoopMPI::sum (eNorm);
               smallChanges = SxLoopMPI::sum (smallChanges ? 0 : 1) == 0;
               fermiTrial.eps.synMPI ();
            }
            if (xNorm > 10. * eNorm) break;
            if (keepOccFixed)  {
               if (smallChanges)  {
                  fermiTrial.focc = fermi.focc;
                  break;
               }
            } else {
               fermiTrial.fermiDistribution (ekt);

               // --- check if occupation number changes are small
               double deltaMu = (fermiTrial.eFermi - fermi.eFermi) 
                              / (xEtaFactor * lambdaT);
               //cout << "dMu real=" << deltaMu << endl;
               
               // less than 10% difference from predicted dMu is ok
               if (fabs(deltaMu - dMu) < 0.1 * fabs(dMu)) break;

               // small absolute change is ok
               double smallE = max(1e-3 * ekt, 1e-8);
               if (fabs(deltaMu) < smallE && fabs(dMu) < smallE) break;

               // negligible changes in occupations is ok
               double deltaOcc = 0.;
               for (ik = 0; ik < nk; ++ik)  {
                  for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
                     deltaOcc += (fermiTrial.focc(iSpin,ik) 
                                  - fermi.focc(iSpin,ik)).absSqr ().sum ()
                                 * Gk.weights(ik);
                  }
               }
               // if (sqrt(deltaOcc) < 1e-6) break;
               if (deltaOcc < 1e-12) break;
                     
            }
            // --- reduce xEta step
            xEtaFactor *= 0.5;
            if (iTry % 3 == 0) slowPoints++;
         }
         //cout << "xEtaFactor =" << xEtaFactor << endl;
         //cout << "lambdaT=" << lambdaT << endl;
         sigma += xEtaFactor * sigmaEps;
      } else {
         fermiTrial = fermi;
      }
      
      if (xNorm > 10. * eNorm && lastRestart < it-2)  
      {
         cout << "eps has gone crazy" << endl;
         lastRestart = it;
         slowPoints = 0;
         cout << "Restarting from subspace diagonalization" << endl;
         for (ik = 0; ik < nk; ik++)  {
            SX_MPI_LEVEL("waves-k");
            if (!SxLoopMPI::myWork (ik)) continue;
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               ham = (-1.) * dEta(ik)(iSpin);
               for (i = 0; i < nStates; ++i)
                  ham(i,i) += fermi.eps(i,iSpin,ik);
               eig = ham.eigensystem ();
               fermi.eps(iSpin,ik) = eig.vals;
               waves(iSpin,ik).rotate (eig.vecs);
               xEta(ik)(iSpin).set (0.);
               X(iSpin,ik).set (0.);
            }
         }
         fermi.eps.synMPI ();
         if (!keepOccFixed) fermi.fermiDistribution (ekt);
         H.computeRho (fermi.focc, waves);
         eTot = H.getEnergy (waves, fermi);
         fOld = freeEnergy;
         freeEnergy = eTot - ekt * fermi.getEntropy ();
         sxprintf ("F = %.12f\n", freeEnergy);
         sxprintf ("\n e(it=%d)=%20.18f H, gamma=%g,sigma=%g\n",
                 it, eTot, gamma, sigma);
         fflush (stdout);
         printEnergyDatLine (it, eTot, freeEnergy,
                             fermi.getEBand (SxFermi::UseFocc),
                             fermi.getEntropy ());
         if ( !(it % printSteps) )  {
            H.printEnergies ();
            if (isMetal) fermi.printOccupation ();
         }
         continue;
      }

      if (xEtaFactor < 5e-3 )  {
         cout << "eta downscaling has become critical" << endl;
         cout << "Switching off psi search" << endl;
         for (ik = 0; ik < nk; ++ik)  {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork (ik))  {
               for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
                  X(iSpin,ik).set (0.);
               }
            }
         }
         sigma = xEtaFactor * sigmaEps;
         traceOld = traceEps;
         etaOnly = true;
      } else {
         etaOnly = false;
      }

      psiTrialPtr = waves.getNew ();

      // --- get trial waves
      for (ik = 0; ik < nk; ik++)  {
         SxPWSet &psiTrial = *psiTrialPtr;
         SX_MPI_LEVEL("waves-k");
         if (!SxLoopMPI::myWork (ik)) continue;
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            psiTrial(iSpin,ik) <<= (waves(iSpin,ik) - lambdaT * X(iSpin,ik));
            S.orthonormalize (&psiTrial(iSpin,ik), orthoMethod);
            if ((isMetal) && (!keepOccFixed))  {
               // --- subspace rotation
               ham = (-lambdaT * xEtaFactor) * xEta(ik)(iSpin);
               for (i = 0; i < nStates; ++i)
                  ham(i,i) += fermi.eps(i,iSpin,ik);
               eig = ham.eigensystem ();
               psiTrial(iSpin, ik).rotate (eig.vecs);
            }
         }
      }

      // --- quadratic line minimization with trial step
      H.computeRho (fermiTrial.focc, *psiTrialPtr);
      eTrial = H.getEnergy (*psiTrialPtr, fermiTrial);
      fTrial = eTrial - ekt * fermiTrial.getEntropy ();
      sxprintf ("Ftrial = %.12f\n", fTrial);
      
      lambdaMin = lineMinimization (freeEnergy, lambdaT, fTrial, sigma, &curv);
      if (fabs(sigma) < 0.1 * fabs(eTot) && !etaOnly) {
         double g = fabs(lambdaMin) / lMinAvg;
         if (g > 2.) g = 1.2;
         if (g < 0.5) g = 0.8;
         lMinAvg *= pow(g, 0.7);
         //cout << "lMinAvg = " << lMinAvg << endl;
      }

      if (eTrial > eTot)
         sxprintf ("SxHamSolver: eTrial(%12.8f) > eTot(%12.8f)\n", 
                 eTrial, eTot);
      sxprintf ("\n eGuess(it=%d,lMin=%g)=%20.18f H\n",
            it, lambdaMin,
            (eTot + sigma*lambdaMin + curv*lambdaMin*lambdaMin));
      fflush (stdout);

      if (testLineMinim)  {
         SX_NO_MPI;
         SxPWSet &psiTrial = *psiTrialPtr;
         double lambdaTorig = lambdaT;
         filename = "lineminim";
         filename += SxString(".it=") + it + SxString(".dat");
         FILE *fpTest = fopen (filename.ascii(), "w");
         fprintf (fpTest, "# ---------------------------------------\n");
         fprintf (fpTest, "# lambdaT     eTot[H]     eFit[H]     eLin[H]\n");
         lambdaTMin = -0.5 * lMinAvg; lambdaTMax = 1.5001 * lMinAvg;
         //lambdaTStep = (lambdaTMax - lambdaTMin) / 17.;;
         lambdaTStep = ltsFac * lMinAvg;

         // --- loop for creating lambdaT-test plot
         for (lambdaT=lambdaTMin; lambdaT<=lambdaTMax; lambdaT+=lambdaTStep)  {
            for (ik = 0; ik < nk; ik++)
               for (iSpin = 0; iSpin < nSpin; iSpin++)  {
                  psiTrial(iSpin,ik) <<= (waves(iSpin,ik) - lambdaT * X(iSpin,ik));
                  S.orthonormalize (&psiTrial(iSpin,ik), orthoMethod);
                  if (isMetal)  {
                     ham = (-lambdaT * xEtaFactor) * xEta(ik)(iSpin);
                     for (i = 0; i < nStates; ++i)
                        ham(i,i) += fermi.eps(i,iSpin,ik);
                     eig = ham.eigensystem ();

                     fermiTrial.eps(iSpin, ik) = eig.vals.real ();
                     psiTrial(iSpin, ik).rotate (eig.vecs);
                  }
               }

            if (keepOccFixed) fermiTrial.focc=fermi.focc;
            else fermiTrial.fermiDistribution (ekt);
            H.computeRho (fermiTrial.focc, psiTrial);
            eTrial = H.getEnergy (psiTrial, fermiTrial);
            
            // --- create testbed output file
            //   = e0   + s l           + c^2 l^2                + ...
            eFit = freeEnergy + sigma*lambdaT + curv*lambdaT*lambdaT;
            if (isMetal)  {
               fprintf (fpTest, "%20.16f\t%20.16f\t%20.16f\t%20.16f\n", 
                     lambdaT, eTrial - ekt * fermiTrial.getEntropy (),
                     eFit, freeEnergy + sigma*lambdaT);
            } else {
               fprintf (fpTest, "%20.16f\t%20.16f\t%20.16f\t%20.16f\n", 
                     lambdaT, eTrial, eFit, eTot + sigma*lambdaT);
            }
            fflush(fpTest);
         }
         fclose (fpTest);
         lambdaT = lambdaTorig;
      }
      // unload psiTrial - no longer needed
      psiTrialPtr->memMinimize ();
      if (!keepWavesOnDisk) psiTrialPtr = SxPtr<SxPWSet> ();

      checkGradient = true;
      if (fabs(sigma) > 0.1 * fabs(eTot)) checkGradient=false;
      if (isMetal && (lambdaMin > lambdaT || lambdaMin < 0.)
          && !etaOnly)  {
         bool smallChanges = true;
         for (ik = 0; ik < nk; ik++)  {
            SX_MPI_LEVEL("waves-k");
            if (!SxLoopMPI::myWork (ik)) continue;
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               ham = (-lambdaMin * xEtaFactor) * xEta(ik)(iSpin);
               for (i = 0; i < nStates; ++i)
                  ham(i,i) += fermi.eps(i,iSpin,ik);
               eig = ham.eigensystem ();
               fermiTrial.eps(iSpin, ik) = eig.vals.real ();
               if (keepOccFixed)  {
                  for (i = 0; i < nStates; ++i)  {
                     double smallEps = 1e-2 * fabs(xEta(ik)(iSpin)(i,i).re);
                     double deltaEps = fermiTrial.eps(i,iSpin,ik)
                                     + lambdaMin * xEtaFactor * xEta(ik)(iSpin)(i,i).re
                                     - fermi.eps(i,iSpin,ik);
                     if (smallEps > 1e-6 && fabs(deltaEps) > smallEps)
                        smallChanges = false;
                  }
                  if (!smallChanges) break;
               }
            }
         }
         {
            SX_MPI_SOURCE(TopLevel, TaskGroupAll);
            SX_MPI_TARGET(TopLevel, TaskGroupAll);
            smallChanges = SxLoopMPI::sum (smallChanges ? 0 : 1) == 0;
            fermiTrial.eps.synMPI ();
         }
         if (!keepOccFixed)  {
            fermiTrial.fermiDistribution (ekt);
            double deltaMu = (fermiTrial.eFermi - fermi.eFermi) 
                           / (xEtaFactor * lambdaMin);
            double smallE = max(1e-3 * ekt, 1e-8);
            bool bigShift = fabs(dMu * xEtaFactor * lambdaMin) > 10. * ekt;
            bool tinyShift = fabs(dMu) < smallE;
            if (!smallChanges
                || (bigShift   && fabs(deltaMu - dMu) > 0.01 * fabs(dMu))
                || (!tinyShift && fabs(deltaMu - dMu) > 0.1 * fabs(dMu) )
                || (tinyShift  && fabs(deltaMu - dMu) > smallE          ))
            {
               // don't jump into dangerous regime, make trial step instead
               cout << "deltaMu (lMin) = " << (deltaMu - dMu) 
                    << "/" << dMu << endl;
               cout << "Changing step from " << lambdaMin;
               if (lambdaMin > 0. || eTrial < eTot)  {
                  lambdaMin = lambdaT;
               } else if (gamma > 0.1) {
                  lambdaMin *= 0.01; // just hope that this is ok...
                  traceOld *= 1000.;
               }
               cout << " to " << lambdaMin << endl;
               // 4-th order step might make things even worse
               checkGradient = false;
            }
         }
      }
      if (lambdaMin < 0.) slowPoints++;

      // --- compute new psi
      if ( fabs(lambdaMin) > 1000. )  {
         cout << "Gigantic lMin" << endl;
         cout << "reducing lMin = " << lambdaMin << " to " << lMinAvg
              << endl;
         lambdaMin = lMinAvg;
      }
      for (ik = 0; ik < nk; ik++)  {
         SX_MPI_LEVEL("waves-k");
         if (!SxLoopMPI::myWork (ik)) continue;
         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            waves(iSpin,ik) -= (lambdaMin * X(iSpin,ik));
            S.orthonormalize (&waves(iSpin,ik), orthoMethod);
            if (isMetal)  {
               ham = (-lambdaMin * xEtaFactor) * xEta(ik)(iSpin);
               for (i = 0; i < nStates; ++i)
                  ham(i,i) += fermi.eps(i,iSpin,ik);
               eig = ham.eigensystem ();

               fermi.eps(iSpin, ik) = eig.vals.real ();
               waves(iSpin, ik).rotate (eig.vecs);
               X(iSpin, ik).rotate (eig.vecs);
               xEta(ik)(iSpin) = eig.vecs.adjoint () ^ xEta(ik)(iSpin) ^ eig.vecs;
            }
         }
      }
      fermi.eps.synMPI ();

      X.memMinimize (); // X currently not needed
      
      // --- SCF update
      //double fermiOld = fermi.eFermi;
      if (!keepOccFixed) fermi.fermiDistribution (ekt);
      //double deltaMu = fabs((fermi.eFermi - fermiOld) 
      //                      - dMu * lambdaMin * xEtaFactor);
      H.computeRho (fermi.focc, waves);
      //if (calcForces) H.calcForces = true;
      eTot = H.getEnergy (waves, fermi);
      /*
      if (calcForces)  {
         H.calcForces = false;
         cout << "Forces="; 
         for (int ia = 0; ia < H.fTotal.getNAtoms (); ++ia)  {
            for (int d = 0; d < 3; ++d)
               sxprintf ("\t%.16f", H.fTotal(ia)(d));
         }
         cout << endl;
      }
      */
      fOld = freeEnergy;
      freeEnergy = eTot - ekt * fermi.getEntropy ();
      sxprintf ("F = %.12f\n", freeEnergy);
      sxprintf ("\n e(it=%d)=%20.18f H, gamma=%g,sigma=%g\n",
              it, eTot, gamma, sigma);
      fflush (stdout);
      printEnergyDatLine (it, eTot, freeEnergy,
                          fermi.getEBand (SxFermi::UseFocc),
                          fermi.getEntropy ());

//    fprintf (fp, "%d\t%12.8f\t%20.15f\n", 
//             it, 0.0, eTot); fflush  (fp);
//      if (lambdaMin < 0.)  { sxprintf ("negative lambdaMin\n"); return; }
      SX_MPI_MASTER_ONLY {
         if (nSpin == 2 && dynamic_cast<SxPAWHamiltonian*>(&*hamPtr))  {
            SxArray<double> spins = SxPtr<SxPAWHamiltonian> (hamPtr)->pawRho.getSpinMom (structure);
            FILE *spinFile = fopen ("spins.dat", "a");
            fprintf(spinFile, "%d ", it);
            SX_LOOP(ia) fprintf(spinFile, " %.6f", spins(ia));
            fprintf(spinFile, "\n");
            fclose (spinFile);
         }
      }
     
      double fGuess = fOld + lambdaMin * sigma + curv * sqr(lambdaMin);
      quadMinOK = fabs(fGuess - freeEnergy) < 3e-2 * fabs(fOld - fGuess);
      if (quadMinOK)  {
         cout << "Quadratic approximation error in line minimisation: " 
              << (fGuess - freeEnergy)/(fOld - fGuess) << endl;
         if (gamma < 0.3 || fixKappa || xEtaFactor < 0.9 || !isMetal)  {
            /// no need to adjust kappa
            checkGradient = false;
         }
      } 
      if (fabs(lambdaMin  - lambdaT) < 0.05 * fabs(lambdaT))
         checkGradient=false;
      if (keepOccFixed) checkGradient=false;
      if (checkGradient && !quadMinOK)  {
         cout << "will check gradient in it=" << (it+1) << ": "
              << 100. * (fGuess - freeEnergy)/(fOld - fGuess)
              << "%"" deviation " << endl;
      }
      if (fabs(lambdaMin) < 0.05 * lMinAvg) slowPoints +=5;

      // --- check convergence criteria
      if ((it + 1) >= minSteps)  {
         if (fabs(fOld - freeEnergy) < dEnergy || (fabs(lambdaMin) == 0.))  {
            sxprintf ("\nConvergence reached.\n");
            break;
         }
      }

      if ( !(it % printSteps) )  {
         H.printEnergies ();
         if (isMetal) fermi.printOccupation ();
      }

      if ( it+1 >= maxSteps )  {
         sxprintf ("WARNING: Maximum number of steps exceeded.\n");
         sxprintf ("         Convergence not yet reached.\n");
      }
   }

   if (finalDiag)  {
      for (ik = 0; ik < nk; ++ik)  {
         SX_MPI_LEVEL("waves-k");
         if (SxLoopMPI::myWork (ik))  {
            for (iSpin = 0; iSpin < nSpin; ++iSpin)  {
               fermi.eps(iSpin,ik) = subspaceDiagonalization(&waves(iSpin,ik));
            }
         }
      }
      fermi.eps.synMPI ();
   }
   
   H.printEnergies (); 
   if (isMetal || finalDiag) fermi.printOccupation (true);
   energy = H.getEnergy () - ekt * fermi.getEntropy ();

   if (changedHam) hamPtr->backToDefault (cmd);
}

void SxHamSolver::scfDiagonalization (const SxSymbolTable *cmd, bool calc)
{
   SX_CLOCK (Timer::ElMinim);

   double dRelRes, maxResidue = 1.;
   int printSteps;
   SX_CHECK (nSpin == wavesPtr->getNSpin(), nSpin, wavesPtr->getNSpin ());
   double dumpTime = 12 * 3600; // 12h
   bool storeRho, storeWaves, propWaves = false;
   bool calcForces = false;
   SxList<SxSymbolTable*> elMinim;
   bool printVarEnergy = false;
   bool useKSenergy = false;
   double epsilon2M = 1e-16;
   try  {
      if (dEnergyLow) {
         dEnergy = (cmd->contains("dEnergyLow"))
            ? cmd->get("dEnergyLow")->toReal ()
            : 1e-8;
      } else {
         dEnergy  = cmd->contains ("dEnergy")
            ? cmd->get("dEnergy")->toReal()
            : 1e-8;
      }
      if (cmd->contains ("maxResidue"))
         maxResidue = cmd->get ("maxResidue")->toReal ();

      ekt       = (cmd->contains("ekt")) 
                ? cmd->get("ekt")->toReal() / HA2EV 
                : ektDefault;
      keepOcc = cmd->contains("keepOccFixed")
              ? cmd->get("keepOccFixed")->toAttribute ()
              : false;
      keepRho   = (cmd->contains("keepRhoFixed")) 
                ? cmd->get("keepRhoFixed")->toAttribute ()
                : false;
      if (cmd->contains("propagateWaves"))
         propWaves = cmd->get("propagateWaves")->toAttribute ();
      if (cmd->contains("useKS"))
         useKSenergy = cmd->get("useKS")->toAttribute ();

      if (keepRho)  dEps = dEnergy;
      dRelRes       = (cmd->contains("dRelRes"))
                    ?  cmd->get("dRelRes")->toReal()
                    :  1e-2;    // from SFHIngX 1.0
      epsilon2M     = (cmd->contains("dSpinMoment"))
	            ? sqr (cmd->get("dSpinMoment")->toReal())
		    : 1e-16;
      maxSteps      = cmd->contains("maxSteps", true)
                    ? cmd->get("maxSteps")->toInt ()
                    : 100;
      printSteps    = cmd->contains("printSteps")
                    ? cmd->get("printSteps")->toInt()
                    : 15;
      if (cmd->contains("dumpTime"))
         dumpTime = cmd->get("dumpTime")->toReal ();
      storeWaves =  (cmd->contains("noWavesStorage"))
                 ? !(cmd->get("noWavesStorage")->toAttribute())
                 :   true;
      storeRho   =  (cmd->contains("noRhoStorage"))
                 ? !(cmd->get("noRhoStorage")->toAttribute())
                 :   true;
      if (cmd->contains ("calcForces"))
         calcForces = cmd->get("calcForces")->toAttribute ();

      if (cmd->contains ("printVarEnergy"))
         printVarEnergy = cmd->get ("printVarEnergy")->toAttribute ();

      if (cmd->contains("spinMoment"))  {
         fermi.setSpinMoment (cmd->get("spinMoment")->toReal (), true);
         fermi.fermiDistribution (ekt);
      }
      if (cmd->contains("keepSpinFixed"))  {
         fermi.setSpinMoment (fermi.getSpinMoment (), 
                              cmd->get("keepSpinFixed")->toAttribute ());
      }

      SxSymbolTable *group;
      for (group = cmd->begin (); group; group=group->nextSibling ())  {
         SxString name = group->getName ();
         if (   name == "CCG" || name == "blockCCG" || name == "rmmDiis"
             || name == "blockRmmDiis" || name == "blockRmmCG")
            elMinim.append (group);
      }
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   if (elMinim.getSize () == 0)  {
      cout << "scfDiag: blockCCG, CCG, rmmDiis, blockRmmDiis, "
           << "or blockRmmCG group must be used!" << endl;
      SX_EXIT;
   }
   
   if (!calc)  return;

   // --- setup mixer
   SxRhoMixer mixer (cmd, nSpin == 2);
   {
      // important: dipole correction may cause convergence problems
      // with Pulay mixer, so switch on mixer correction:
      // TODO: ugly
      //What Hamiltonian?
      SxPWHamiltonian *pwHamPtr
         = dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ());
      if (pwHamPtr)
         mixer.linearVacMixing |= pwHamPtr->dipoleCorrection;
      
      SxPAWHamiltonian *pawHamPtr
         = dynamic_cast<SxPAWHamiltonian *>(hamPtr.getPtr ());
      if (pawHamPtr)  
         if (pawHamPtr->dipoleCorr) mixer.linearVacMixing = true;
   }

   // --- print input parameter
   cout << SX_SEPARATOR;
   cout << "| SCF calculation\n";
   cout << SX_SEPARATOR;
   cout << "|   max. steps:                " << maxSteps << endl;
   cout << "|   details printed every:     " << printSteps << " iteration(s)\n";
   cout << "|   dEnergy:                   " << dEnergy  << " H\n";
   cout << "|   dEps/R(rho) factor:        " << dRelRes << "\n";
   cout << "|   total energy functional:   " 
        << (useKSenergy ? "Kohn-Sham" : "Harris-Foulkes") << endl; 
   if (elMinim.getSize () > 1) cout << SX_SEPARATOR;
   for (SxList<SxSymbolTable*>::Iterator elIt = elMinim.begin ();
        elIt != elMinim.end (); ++elIt)
   {
      SxString name = (*elIt)->getName (); 
      if (name == "blockCCG")  {
         blockStateByStateCG(*elIt, PrintOnly);
      } else if (name == "rmmDiis")  {
         rmmDiis (*elIt, !keepRho, true);
      } else if (name == "blockRmmDiis")  {
         blockRmmDiis (*elIt, !keepRho, true);
      } else if (name == "blockRmmCG")  {
         blockRmmCG (*elIt, !keepRho, true);
      } else if (name == "CCG")  {
         stateByStateCG(*elIt, true, /* print= */ true);
      } else {
         // never should arrive here
         SX_EXIT;
      }
      if ((*elIt)->contains ("dEnergy", true))  {
         double dE = (*elIt)->get ("dEnergy")->toReal ();
         cout << "|      until dE falls below:   " << dE << " H\n";
         cout << SX_SEPARATOR;
      }
   }
   if (!keepRho)  {
      cout << "|   mixer:\n";
      mixer.print ();
   }
   cout << "|   type of calculation:       ";
   if (keepRho)  {
      cout << "band structure run\n";
      if (propWaves)
         cout << "|   initial waves:             from previous k-point\n";
   } else {
      cout << "SCF minimization\n";
      if (keepOcc)
         cout << "|   occupation:                " << "fixed\n";
      else  {
         cout << "|   ekt:                       " << ekt * HA2EV << " eV"
              << endl;
         if (ekt * HA2EV > 1e-5)  {
            cout << "|   smearing scheme:           ";
            fermi.printSmearing ();
         }
      }
   }
   if (fermi.keepSpinMoment)
      cout << "|   spin moment:               " << fermi.spinMoment
           << " (fixed)\n";
   if (spinConstraint)  {
      cout << "|   atomic spin constraints:   enabled" << endl;
      cout << "|   dSpinMoment:               " << sqrt(epsilon2M) << endl;
   }
   cout << "|   data dumping:              ";
   if (dumpTime <= 0.)
      cout << "every iteration.";
   else if (dumpTime <= 10000.)
      cout << int(dumpTime / 60.) << " min";
   else
      cout << int(dumpTime / 3600.) << " hours";
   cout << endl;
   cout << SX_SEPARATOR;

   // --- self consistent loop
   SxList<SxSymbolTable*>::Iterator elIt = elMinim.begin ();
   bool changedHam = hamPtr->rereadTable (*elIt);
   PrecEnergy eBand, lastEnergy = 1e50;    // just a large value
   double dEpsFromdRes = dRelRes;
   double lastR, firstR = 1.;
   double avgRho = fermi.nElectrons / structure.cell.volume;
   int lastRefineForStep = -1;
   int nUnConv = wavesPtr->getNStates ();
   SxTimer lastDump(1);
   lastDump.start ();


   double eTot = -1.; // note: will contain last computed energy
   for (int it=0; it < maxSteps ; it++)  {

      // --- dump data if it becomes time 
      {
         lastDump.stop ();
         bool dumpNow = lastDump.getTime () > dumpTime;
         dumpNow = SxLoopMPI::lor (dumpNow);
         if (dumpNow)  {
            writeData (storeWaves, storeRho);
            lastDump.reset ();
         }
         lastDump.start ();
      }
      
      lastR = mixer.getNormR ();
      mixer.addRhoIn (hamPtr->getRho());
      // --- get dEps
      if (!keepRho)  {
         if (it > 0)  {
            double res = mixer.getNormR () / avgRho;
            dEps = dEpsFromdRes * pow(res, 1.5);
         } else {
            dEps = 1e-6;
         }
         //cout << "dRelRes = " << dEpsFromdRes << endl;
         //cout << "dEps = " << dEps << endl;
      }


      // ############## khr: most of the computing time goes here
      //
      // --- diagonalize H[rhoIn]|Psi> = e|Psi> and get rhoOut
      if ((*elIt)->getName () == "CCG")  {
         stateByStateCG (*elIt, !keepRho);
      } else if ((*elIt)->getName () == "blockCCG")  {
         if (keepRho && propWaves && it == 0)
            blockStateByStateCG (*elIt, PropagateRun);
         else
            nUnConv = blockStateByStateCG (*elIt, keepRho ? BandStructureRun 
                                                          : SCFrun);
      } else if ((*elIt)->getName () == "rmmDiis")  {
         rmmDiis (*elIt, !keepRho);
      } else if ((*elIt)->getName () == "blockRmmDiis")  {
         blockRmmDiis (*elIt, !keepRho);
      } else if ((*elIt)->getName () == "blockRmmCG")  {
         blockRmmCG (*elIt, !keepRho);
      } else  {
         SX_EXIT;
      }
      if (spinConstraint) 
      {
         SxPAWHamiltonian  &pawHam = *SxPtr<SxPAWHamiltonian>(hamPtr);
         spinConstraint->setNuA (0.);
         pawHam.nuA -= spinConstraint->computeNu(wavesPtr,fermi,ekt,epsilon2M); // wave functions get rotated
      }
      // ############## khr: most of the computing time goes here


      // ### some 5% are below
      if (!keepOcc)  fermi.fermiDistribution (ekt);
      double doubleCounting = hamPtr->getDoubleCounting ();
      if (keepRho)  {
         // TODO: should be compute(...)
         //hamPtr->getEnergy (*wavesPtr, fermi); // why? CF 2009-07-25
      }  else  {
         {
            SX_CLOCK (Timer::khrDebug1);
            hamPtr->computeRho (fermi.focc, *wavesPtr);
            if ( dynamic_cast<SxPAWHamiltonian*>(&*hamPtr))  {
               double spin = SxPtr<SxPAWHamiltonian>(hamPtr)->pawRho.getSpin ();
	       sxprintf("The spin is %.16f\n", spin);
            }
               
            // --- the variational energy
            if (printVarEnergy)  {
               hamPtr->getEnergy (*wavesPtr, fermi);
               sxprintf ("variational free energy: F(%d)=%.16f\n",
                        it+1, hamPtr->getEnergy () - ekt * fermi.getEntropy ());
            }
         }
         // --- mix densities rhoIn and rhoOut
         {
            SX_CLOCK (Timer::khrDebug2);
            mixer.addRhoOut (hamPtr->getRho ());
         }
         {
            SX_CLOCK (Timer::khrDebug3);
            hamPtr->getRho () = mixer.getMixedRho ();
         }
         if (calcForces)
            if (SxPtr<SxPWHamiltonian> H = hamPtr) H->calcForces = true;

         // TODO: should be compute(...)
         {
            SX_CLOCK (Timer::khrDebug4);
            hamPtr->getEnergy (*wavesPtr, fermi);
         }
         if (calcForces) {
            SxAtomicStructure F;
            if (SxPtr<SxPAWHamiltonian> pawHam = hamPtr)  {
               pawHam->computeForces (*wavesPtr, fermi);
               F = pawHam->forces;
            } else if (SxPtr<SxPWHamiltonian> H = hamPtr) {
               H->calcForces = false;
               F = H->fTotal;
            }
            if (F.getSize () > 0)  {
               int iDof = 0;
               const char *xyz = "xyz";
               SX_LOOP3(is,ia,iDir)
                  sxprintf("%d %3s%d %c %.16f\n", ++iDof,
                           potPtr->chemName(is).ascii (), int(ia+1),
                           xyz[iDir], F(is,ia)(iDir));
            }
         }
         // --- refine dRelRes
         if (it > 0)  {
            double newR = mixer.getNormR ();
            if (lastR > newR)  {
               // roughen calculation if last step was better than average
               if (log(firstR/lastR) < it * log(lastR/newR))
                  dEpsFromdRes *= 2.;
               // refine calculation if last step was much worse than average
               if (log(firstR/lastR) > 1.3 * it * log(lastR/newR))
                  dEpsFromdRes *= 0.5;
               if (dEpsFromdRes > dRelRes) dEpsFromdRes = dRelRes;
            } else {
               // refine calculation when residue grows, but not
               // twice in a row
               if (lastRefineForStep < it) {
                  dEpsFromdRes *= .1;
                  lastRefineForStep = it + 1;
               }
            }
         } else {
            firstR = mixer.getNormR ();
         }
      }

      // --- print results
      eBand  = fermi.getEBand (keepRho ? SxFermi::NoFocc 
                                       : SxFermi::UseFocc);
      eTot = hamPtr->getEnergy (); // Kohn-Sham functional
      double eHarris = (keepRho ? fermi.getEBand (SxFermi::UseFocc) : eBand)
                     + doubleCounting;
      if (!useKSenergy) {
         eTot = eHarris;
      }
      
      energy = (keepRho) ? eBand : eTot;
      sxprintf ("eTot(%d)=%15.12f, eBand=%15.12f,  |R|=%15.12f\n", 
               it+1, eTot, eBand, mixer.getNormR());
      if (nSpin == 2)  {
         sxprintf ("                                        "
                 "        |RPol|=%15.12f\n",
                  mixer.getNormRPol());
      }
      double entropy = fermi.getEntropy ();
      sxprintf ("F(%d)=%15.12f\n", it+1, eTot - ekt * entropy);
      fflush (stdout);

      SX_MPI_MASTER_ONLY
      {
         printEnergyDatLine (it, eTot, eTot - ekt * entropy,
                             eBand, entropy);
         SxFileIO::appendToFile
         (   SxString (it)                                + "\t" // iteration
           + SxString(mixer.getNormR())                          // |R|
           + ((nSpin == 1) ? SxString("\n")
                           : ("\t" + SxString(mixer.getNormRPol()) + "\n")) // |Rpol|
         , "residue.dat");
      }

      SX_MPI_MASTER_ONLY {
         if (nSpin == 2 && dynamic_cast<SxPAWHamiltonian*>(&*hamPtr))  {
            SxArray<double> spins = SxPtr<SxPAWHamiltonian> (hamPtr)->pawRho.getSpinMom (structure);
            FILE *spinFile = fopen ("spins.dat", "a");
            fprintf(spinFile, "%d ", it);
            SX_LOOP(ia) fprintf(spinFile, " %.6f", spins(ia));
            fprintf(spinFile, "\n");
            fclose (spinFile);
         }
      }
 



      // --- convergence?
      if ( (fabs (energy - lastEnergy) < dEnergy && mixer.getNormR () < maxResidue)
          || (keepRho && nUnConv <= 0))
      {
         sxprintf ("Convergence reached.\n");
         break;
      }
      if ( it+1 >= maxSteps )  {
         sxprintf ("WARNING: Maximum number of steps exceeded.\n");
         sxprintf ("         Convergence not yet reached.\n");
         break;
      }

      if ((*elIt)->contains ("dEnergy", true))  {
         double dE = (*elIt)->get("dEnergy")->toReal ();
         // proceed to next minimizer ?
         if (fabs(energy - lastEnergy) < dE) {
            cout << SX_SEPARATOR;
            cout << "| Energy limit of " << dE << " Hartree reached." << endl;
            if (*elIt != elMinim.last ()) {
               elIt++;
               cout << "| Switching to " 
                    << (*elIt)->getName () << "." << endl;
               // --- allow for Hamiltonian changes
               if (hamPtr->rereadTable (*elIt)) changedHam = true;
               cout << SX_SEPARATOR;
               printTiming (/* restart= */true);
            } else {
               cout << "| No other minimizer available..." << endl;
               cout << SX_SEPARATOR;
            }
         }
      }

      if ( !(it % printSteps) )  {
         hamPtr->printEnergies ();
         fermi.printOccupation ();
      }
      cout.flush ();

      lastEnergy = energy;

   }


   hamPtr->printEnergies ();
   fermi.printOccupation (true);
   energy = eTot - ekt * fermi.getEntropy ();

   if (changedHam) hamPtr->backToDefault (cmd);
}

void SxHamSolver::subspaceDiagonalization (const SxSymbolTable *cmd, bool calc)
{
   bool fermiDist = false;
   try  {
      ekt = cmd->contains("ekt") ? cmd->get("ekt")->toReal() / HA2EV
                                  : ektDefault; 
      if (cmd->contains ("fermiDistribution"))
         fermiDist = cmd->get("fermiDistribution")->toAttribute ();
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   cout << SX_SEPARATOR;
   cout << "| Subspace diagonalization" << endl;
   if (fermiDist)  {
      cout << "| with subsequent recalculation of occupations" << endl;
      if (fabs(ekt - ektDefault) > 1e-10)
         cout << "| with ekt= " << (ekt * HA2EV) << " eV" << endl;
      cout << "| smearing scheme: "; fermi.printSmearing ();
   }
   cout << SX_SEPARATOR;
   if (!calc) return;
   for (int ik = 0; ik < wavesPtr->getNk (); ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         fermi.eps(iSpin, ik) = subspaceDiagonalization (&(*wavesPtr)(iSpin, ik));
      }
   }
   fermi.eps.synMPI ();
   if (fermiDist)
      fermi.fermiDistribution (ekt);
}

SxDiracVec<TPrecEps> 
SxHamSolver::subspaceDiagonalization (PsiGI *wavesPtrIN,
                                      int subDiagSize,
                                      int subDiagOvlp)
{
   SX_CHECK (wavesPtrIN);
   PsiGI &psiAll = *wavesPtrIN;
   // --- subspace diagonalization ---
   
   int ng = (int)psiAll.nRows ();
   SX_CHECK(nStates == psiAll.nCols (), nStates, psiAll.nCols ());
   
   SxDiracVec<TPrecEps> epsAll(nStates);
   int from = 0;
   int to = subDiagSize > 1 ? subDiagSize 
                            : nStates;

   PsiGI dPsiAll;
   // calculate dPsi for all states
   {  SX_CLOCK (Timer::SubspaceMatrix);
      dPsiAll = *hamPtr | psiAll;
   }
   int ng2 = (int)dPsiAll.nRows ();
   
   int nBlock, nextFrom;
   PsiGI psiGI, dPsiGI;
   SxDiracMat<TPrecCoeffG> h;
   SxDiracSymMat<TPrecCoeffG> hSym;
   SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;
   double eps;
   for (int iSubDiag = 0; /* empty */; iSubDiag++)  {

      if (iSubDiag > nStates)  {
         cout << "number of diagonalization exceeds reasonable limit.";
         cout << " There's something wrong." << endl;
         cout << epsAll << endl;
         SX_EXIT;
      }
      if (to > nStates) to = nStates;

      nBlock = to - from;
      //cout << "Diagonalization " << (iSubDiag+1);
      //cout << " from " << from;
      //cout << " to " << (to-1) << endl;

      // --- set up Hamilton matrix for block
      psiGI = PsiGI ();
      psiGI = psiAll(SxIdx(from*ng, to*ng - 1));
      psiGI.reshape (ng, nBlock);
      { 
         SX_CLOCK (Timer::SubspaceMatrix);
         /*
         // TODO
         // h = H.getMatrix (psiGI);
         dPsiGI.reformat (ng, nBlock);
         for (i = from; i < to; i++)
           dPsiGI.colRef(i-from) <<= (H | psiAll(i,iSpin,ik));
         */
         dPsiGI = PsiGI ();
         dPsiGI = dPsiAll (SxIdx(from * ng2, to * ng2 - 1));
         dPsiGI.reshape (ng2, nBlock);
         h   = psiGI.overlap (dPsiGI, ng2);
         /*
         if (verbose)  {
            for (i = from; i < to; i++)  {
               for (j = from; j < to; j++)  {
                  cout << "<" << i << "|" << j << ">=";
                  cout << (h(i-from,j-from) * HA2EV) << endl;
               }
            }
         }
         */
      }

      // --- diagonalize matrix
      { SX_CLOCK (Timer::SubspaceDiag);
         
         // TODO: ugly!!! (with X-Press not necessary anymore)
         // free space, so we can reuse it if we need the same amount
         // this makes use of the memory allocation caching and
         // is a workaround for the missing SxSYMMAT::reformat ()
         hSym = SxDiracSymMat<TPrecCoeffG> ();
         hSym = SxDiracSymMat<TPrecCoeffG> (nBlock);
         { 
            int i,j;
            for (i=0; i < nBlock; ++i)
               for (j=i; j < nBlock; ++j)
                  hSym(i,j) = h(i,j);
         }

         eig = hSym.eigensystem();
      }

      // copy eigenvalues
      epsAll(SxIdx(from,to-1)) <<= eig.vals;
      // copy eigenfunctions
      psiGI <<= (psiGI ^ eig.vecs);

      if (to == nStates) break;

      nextFrom = to - subDiagOvlp;
      // include degenerate states
      eps = epsAll(nextFrom);
      while (nextFrom > 0 && fabs(eps - epsAll(nextFrom - 1)) < 1e-4)
      {
         nextFrom--;
      }

      // --- update dPsiAll for overlapping states
      {
         PsiGI dPsiOvlp, eigVecOvlp;
         int overlap = to - nextFrom;
         dPsiOvlp = dPsiAll(SxIdx(ng2*nextFrom, ng2 * to - 1));
         dPsiOvlp.reshape (ng2, overlap);
         eigVecOvlp = eig.vecs(SxIdx(nBlock * (nBlock - overlap),
                                     (int)eig.vecs.getSize () - 1));
         eigVecOvlp.reshape(nBlock, overlap);
         dPsiOvlp <<= (dPsiGI ^ eigVecOvlp);

      }

      from = nextFrom;
      to = from + subDiagSize;

   }
      
   // orthogonalize degenerated states
   PsiG psiI;
   PrecCoeffG scp;
   int i,j;
   SxOverlap S = hamPtr->getS ();
   for (i = 0; i < nStates; i++)  {
      psiI = psiAll.colRef(i);
      eps = epsAll(i);
      for (j = i-1; j>=0; j--)  {
         if (fabs(eps - epsAll(j)) > 1e-6) continue;
         if (fabs(eps - epsAll(j)) > 1e-2) break;
         scp = (psiAll.colRef(j) | S | psiI);
         if (scp.absSqr () > 1e-20)  {
            /*
            cout << "<" << i << '|' << j << "> = " << scp;
            cout << " (dEps =" << fabs(eps - epsAll(j));
            cout << ")" << endl;
            */
            psiI.plus_assign_ax (-scp, psiAll.colRef(j));
            S.normalize (&psiI);
         }
      }
   }
   VALIDATE_VECTOR (epsAll);
   return epsAll;
}

void SxHamSolver::propagateWaves (int ik, int jk)
{
   SX_NO_MPI;
   SxPWSet &waves = *wavesPtr;
   // --- copy lattice-periodic part from previous k-point
   SX_CHECK (ik >= 0 && ik < waves.getNk (), ik);
   SX_CHECK (jk >= 0 && jk < waves.getNk (), jk);

   SxArray<PsiGI> oldWavesNewBasis(nSpin);
   // important: jk must not be accessed before last access to
   // ik in case the waves are on disk
   const SxGBasis &Gi = waves.getGkBasis()(ik);
   const SxGBasis &Gj = waves.getGkBasis()(jk);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      oldWavesNewBasis(iSpin) = Gj.transferWaves (Gi | waves(iSpin,ik));
      oldWavesNewBasis(iSpin).handle->auxData.ik = jk;
   }

   // project to final basis (e.g. PAW) or copy data
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      waves(iSpin,jk) <<= (*waves(iSpin,jk).getBasisPtr ())
                          | oldWavesNewBasis(iSpin);
}

void SxHamSolver::diffEqtn (const SxSymbolTable *cmd, bool calc)
{
   SX_NO_MPI;

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SX_CLOCK (Timer::ElMinim);
   SX_CHECK (cmd);

   int maxItDiag = 5;
   try  {
      if (cmd->contains("dEnergy"))
         dEps = cmd->get("dEnergy")->toReal ();
      if (cmd->contains("maxStepsCCG"))
         maxItDiag   = cmd->get("maxStepsCCG")->toInt ();
           
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }



   // --- print input parameter
   cout << SX_SEPARATOR;
   cout << "| Band structure state-by-state conjugate gradient\n";
   cout << SX_SEPARATOR;
   cout << "|   conjugate gradient:\n";
   cout << "|      max. steps per state:   " << maxItDiag << endl;
   cout << "|      energy convergence      " << dEps << " H\n";
   cout << SX_SEPARATOR;

   if (!calc) return;

   int itDiag;
   SxDiracVec<TPrecCoeffG> psiI, hPsi, dPsi, dPsiOld, K, X, Xold, hX;
   SxDiracVec<TPrecCoeffG> psiGI, dPsiGI;
   SxDiracMat<TPrecCoeffG> h;
   Real8 eps, gamma = 0., theta, epsPrev=0.;
   PrecCoeffG tr, trOld = 0.;
   SxMatrix<TPrecCoeffG> M(2,2);

   cout.precision (10);
            
   // --- conjugate gradient ------
   psiI = waves(0,0,0);

   for (itDiag=0; /* empty */; itDiag++)  {

      // --- compute gradient: |dPsi> = H |psi>
      hPsi = *hamPtr | psiI;
      eps  = (psiI ^ hPsi).chop().re;
      cout << itDiag << ": eps: " << eps << ", dEps: "
         << fabs(eps - epsPrev) << endl;
      dPsi =  hPsi;

      // --- convergence on dEps
      if (fabs(eps - epsPrev) < dEps)  {
         if (itDiag > 0)  {
            break;
         }
      }

      // TODO: Preconditioner
      K = hamPtr->preconditioner(psiI) * dPsi;

      // --- get conjugate direction, ref1 (8.14b) / ref2 (10.6.7)
      tr = (dPsi ^ K).chop();
      if (itDiag > 0)  gamma = (tr - (dPsiOld^K).chop()) / trOld;
      dPsiOld.copy (dPsi);
      trOld = tr;

      // --- compute search direction, ref1 (8.14a)
      if (itDiag == 0) 
         X = K;
      else              X = K + gamma * Xold;
      Xold.resize (0); Xold.copy (X);
      X.normalize ();
      
      // --- get M matrix (in ref1: H matrix)
      hX = (*hamPtr | X);
      M(0,0) = eps;
      M(0,1) = (X ^ hPsi).chop();
      M(1,0) = M(0,1).conj();  M(1,1) = (X ^ hX).chop();

      // --- uniform transformation, ref1 (8.16)
      if ( M(0,0).re - M(1,1).re == 0. )
      {
         theta = 0.;
      }  else  {
         theta = 0.5 * atan (  (M(1,0).re + M(0,1).re) 
               / (M(0,0).re - M(1,1).re) );
         if (M(0,0).re > M(1,1).re)  theta -= PI/2.;
      }

      SX_CHECK_NUM (theta);
      if (fabs(theta) < 1e-9)  {
         break;
      }

      // --- compute new wavefunction, ref1 (8.15b)
      waves(0,0,0) <<= cos(theta) * psiI + sin(theta) * X;

      // --- convergence?
      if (itDiag+1 >= maxItDiag)  {
         break;
      }
      epsPrev = eps;

   } // itDiag

   // another call of hPsi necessary only to compute some output parameters
   // TODO: There should be a BETTER SOLUTION to do this
   cout << "final H|psi> call: " << endl;
   hPsi = *hamPtr | psiI;
   //double finalEnergy;
   //finalEnergy = hamPtr->getEnergy(waves, fermi);

   // --- print convergence statistics
   cout << SX_SEPARATOR;
   if ( fabs (eps - epsPrev) < dEps )  {
      sxprintf ("Convergence reached.\n");
   } else 

   if ( itDiag+1 >= maxItDiag )  {
      sxprintf ("WARNING: Maximum number of steps exceeded.\n");
      sxprintf ("         Convergence not yet reached.\n");
   }
   cout << "steps needed: " << itDiag << endl;
   cout << SX_SEPARATOR;
}

void SxHamSolver::hSqrCG (const SxSymbolTable *cmd, bool calc)
{
   SX_NO_MPI;
   SX_CLOCK (Timer::ElMinim);
   SX_CHECK (cmd);

   int maxStepsCCG = 5, printSteps = 20;
   int subDiagSize = -1, subDiagOvlp = -1;
   int nSloppy = 0;
   bool verbose = false, printRes = false;
   double dRes = -1.;
   bool autoSteps = false;
   int spareSteps;
   enum WavesGuess { None, Propagate } guess = None;
   bool keepOccFixed = false;
   bool noPreconditioning = false;
   try  {
      if (cmd->contains("dEps"))
         dEps = cmd->get("dEps")->toReal ();
      if (cmd->contains("dEnergy"))
         dEnergy   = cmd->get("dEnergy")->toReal ();
      else
         dEnergy   = dEps * nStates;
      if (!(cmd->contains("dEps")))
         dEps = dEnergy / nStates;
      if (cmd->contains("maxStepsCCG"))
         maxStepsCCG   = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("maxSteps"))
         maxSteps      = cmd->get("maxSteps")->toInt ();
      else
         maxSteps      = 20;
      if (cmd->contains("printSteps"))
         printSteps    = cmd->get("printSteps")->toInt();
      if (cmd->contains("autoSteps"))
         autoSteps = cmd->get("autoSteps")->toAttribute ();
      if (cmd->contains("verbose"))
         verbose = cmd->get("verbose")->toAttribute ();
      if (cmd->contains("noPreconditioning"))
         noPreconditioning = cmd->get("noPreconditioning")->toAttribute ();
      if (cmd->contains("printResidue"))
         printRes = cmd->get("printResidue")->toAttribute ();
      if (cmd->contains("dPsi"))
         dRes = cmd->get("dPsi")->toReal ();
      if (cmd->contains("nSloppy"))
         nSloppy = cmd->get("nSloppy")->toInt ();
      if (cmd->containsGroup ("subspaceDiag"))  {
         SxSymbolTable *subspaceGroup = cmd->getGroup ("subspaceDiag");
         subDiagSize = subspaceGroup->get("maxSize")->toInt ();
         subDiagOvlp = subspaceGroup->get("overlap")->toInt ();
         if (subDiagSize <= subDiagOvlp)  {
            cout << SxString(
                    "subspace block overlap must be (even considerably) "
                    "smaller than maximum block size."
                    ).wrap () << endl;
            SX_QUIT;
         }
      }
      if (cmd->contains("propagateWaves"))
         guess = cmd->get("propagateWaves")->toAttribute () ? Propagate : None;

      if (cmd->contains("keepOccFixed"))
         keepOccFixed = cmd->get("keepOccFixed")->toAttribute ();
           
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   // --- print input parameter
   cout << SX_SEPARATOR;
   cout << "| Band structure state-by-state conjugate gradient\n";
   cout << SX_SEPARATOR;
   cout << "|   max. steps:                " << maxSteps << endl;
   cout << "|   details printed every:     " << printSteps << " iteration(s)\n";
   cout << "|   dEnergy:                   " << dEnergy  << " H\n";
   cout << "|   sloppy (may not converge): " << nSloppy << " states\n";
   cout << "|   initial waves:             ";
   if (guess == Propagate)
      cout << "from previous k-point" << endl;
   else
      cout << "from initialGuess" << endl;
   cout << "|   conjugate gradient:\n";
   cout << "|      max. steps per state:   " << maxStepsCCG << endl;
   cout << "|      final stage tuning:     " << (autoSteps ? "on" : "off") 
                                             << endl;
   cout << "|      energy convergence      " << dEps << " H\n";
   cout << "|      residue convergence :   ";
   if (dRes > 0.)
      cout << dRes << endl;
   else
      cout << "off" << endl;
   cout << "|   subspace diagonalization:  ";
   if (subDiagSize > 0)  {
      cout << endl;
      cout << "|      max. blocksize:         " << subDiagSize << endl;
      cout << "|      block overlap:          " << subDiagOvlp << endl;
   } else {
      cout << "all states" << endl;
   }
   if (keepOccFixed)
      cout << "|   No occupation number update at the end." << endl;
   cout << SX_SEPARATOR;

   if (!calc) return;

   SX_CHECK (dynamic_cast<SxPW*>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW&>(*wavesPtr);
   int ik, nk         = waves.getNk ();
   int iSpin;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   int i;
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   int it, itDiag;
   int nConvRes;
   int maxItDiag;
   SxDiracVec<TPrecCoeffG> psiI, hPsi, hPsi2, hPsiSqr, dPsi, dPsi0, dPsiOld, K, X, Xold, hX;
   SxDiracVec<TPrecCoeffG> psiGI, dPsiGI;
   SxDiracMat<TPrecCoeffG> h;
   SxDiracSymMat<TPrecCoeffG> hSym;
   SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;
   Real8 eps, epsSqr, gamma = 0., theta, epsPrev=0., eTrial;
   PrecCoeffG scp, tr, trOld = 0.;
   SxMatrix<TPrecCoeffG> M(2,2), U(2,2);

   double eBand, lastEnergy;
   eTrial = hamPtr->getETrial();

   int convergeLimit, nUnconvDiag;

   cout.precision (10);
   int step = 0;
   for (ik=0; ik < nk; ik++)  {

      if (guess == Propagate && ik > 0)  {
         propagateWaves (ik-1, ik);
         for (iSpin = 0; iSpin < nSpin; ++iSpin) {
            fermi.eps(iSpin, ik) 
               = subspaceDiagonalization (&waves(iSpin,ik), subDiagSize,
                                          subDiagOvlp);
         }
      }
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         convergeLimit = nStates - 1 - nSloppy;
         lastEnergy = fermi.eps(iSpin,ik).sum (0, convergeLimit);
         for (it = 0; it < maxSteps; it++)  {
            nConvIt = 0;
            nConvAbs = 0;
            nConvRes = 0;
            nCosAcc = 0;
            nUnconvDiag = 0;
            spareSteps = 0;
            
            // --- conjugate gradient ------
            for (i=0; i < nStates; i++)  {
               psiI = waves(i,iSpin,ik);
               

               if (verbose)  {
                  sxprintf ("state: (%d,%d,%d) ", i, iSpin, ik); fflush(stdout);
               }

               maxItDiag = maxStepsCCG;
               if (autoSteps && spareSteps > 0)  {
                  if (spareSteps > maxStepsCCG)  {
                     maxItDiag += maxStepsCCG;
                     spareSteps -= maxStepsCCG;
                  } else  {
                     maxItDiag += spareSteps;
                     spareSteps = 0;
                  }
               }
               epsPrev = fermi.eps(i, iSpin, ik);
               for (itDiag=0; /* empty */; itDiag++)  {

                  // orthogonalize psiI to lower states
                  waves.setOrthogonal (&psiI, i, iSpin, ik,
                                       SxPW::NORMALIZE, 0.);
                  
                  // --- compute gradient: |dPsi> = (H-e) |psi>
                  hPsi = *hamPtr | psiI;
                  hPsi2 = *hamPtr | hPsi;
                  eps  = (psiI ^ hPsi).chop().re;
                  hPsiSqr = hPsi2 - (2. * eTrial) * hPsi
                     + (eTrial * eTrial) * psiI;
                  epsSqr = (psiI ^ hPsiSqr).chop().re;
                  dPsi =  // (H - eTrial)^2 Psi - eps^2 Psi
                     hPsiSqr - epsSqr * psiI;
                  fermi.eps(i,iSpin,ik) = eps;

                  if (printRes)  {
                     cout << "(" << i << "," << iSpin << "," << ik << ").";
                     cout << itDiag << ": res = " << dot(dPsi, dPsi).re;
                     cout << "; dEps = " << (eps - epsPrev) << endl;
                  }

                  // --- convergence on dEps
                  if (fabs(eps - epsPrev) < dEps)  {
                     if (itDiag > 0)  {
                        if (verbose)
                           sxprintf ("conv due to abs eps (%d of %d)\n", 
                                   itDiag+1, maxItDiag);
                        nConvAbs++;
                        break;
                     }
                     // --- for first iteration, residue provides a good
                     // indicator
                     if (dRes > 0. && dot(dPsi, dPsi).re < dRes)  {
                        if (verbose)
                           sxprintf ("conv due to residue (%d of %d) dEps=%g\n",
                                   itDiag+1, maxItDiag,
                                   eps - epsPrev);
                        nConvRes++;
                        break;
                     }
                  }

                  // --- orthogonalize preconditioned gradient
                  step++;
                  if (noPreconditioning)
                     K = dPsi;
                  else
                     K = hamPtr->preconditioner (psiI) * dPsi;
                  waves.setOrthogonal (&K, i+1, iSpin, ik, 
                                       SxPW::DONT_NORMALIZE, 0.l);

                  // --- get conjugate direction, ref1 (8.14b) / ref2 (10.6.7)
                  tr = (dPsi ^ K).chop();
                  if (itDiag > 0)  gamma = (tr - (dPsiOld^K).chop()) / trOld;
                  dPsiOld.copy (dPsi);
                  trOld = tr;

                  // --- compute search direction, ref1 (8.14a)
                  if (itDiag == 0)  X = K * 1.;
                  else              X = K + gamma * Xold;
                  Xold.resize (0); Xold.copy (X);
                  { SX_CLOCK (Timer::ortho);
                    scp = dot (psiI, X);
                    X -= psiI * scp;
                  }
                  X.normalize ();
                  // TODO: ugly
                  X.handle->auxData.i     = i;
                  X.handle->auxData.iSpin = iSpin;
                  X.handle->auxData.ik    = ik;
                  // --- get M matrix (in ref1: H matrix)
                  hX = (*hamPtr | X) - eTrial * X;
                  M(0,0) = eps * eps - 2. * eps * eTrial + eTrial * eTrial;
                  M(0,1) = (hX ^ (hPsi - eTrial * waves(i, iSpin, ik))).chop();
                  M(1,0) = M(0,1).conj();  M(1,1) = (hX ^ hX).chop();

                  // --- uniform transformation U, ref1 (8.16)
                  if ( M(0,0).re - M(1,1).re == 0. )  // '==' is correct!
                  {
                     /* This is a workaround for the numerical limit of
                        degenerate eigenstates. If an exact eigenstate has
                        been found, and X is a different eigenstate to the
                        same energy, the formula below doesn't work, as
                        atan(1/0) is undefined. Any theta is fine, so
                        we just take 0 to get out at the numerical limit
                        test. */
                     theta = 0.;
                  }  else  {
                     theta = 0.5 * atan (  (M(1,0).re + M(0,1).re) 
                                         / (M(0,0).re - M(1,1).re) );
                     // --- in case of being far off the quadr. regime theta
                     //     points to the maximum rather than to the minimum!  
                     if (M(0,0).re > M(1,1).re)  theta -= PI/2.;
                  }

                  SX_CHECK_NUM (theta);
                  if (fabs(theta) < 1e-9)  {
                     if (verbose)
                        sxprintf ("Numerical limit reached (%d of %d)\n",
                                itDiag+1, maxItDiag);
                     nCosAcc++;
                     break;
                  }
                  U(0,0) = cos (theta); U(0,1) = sin (theta);
                  U(1,0) = -U(0,1);     U(1,1) = U(0,0);

                  // --- compute new wavefunction, ref1 (8.15b)
                  waves(i,iSpin,ik) <<= U(0,0) * psiI + U(0,1) * X;

                  // --- convergence?
                  if (itDiag+1 >= maxItDiag)  {
                     if (verbose)  {
                        sxprintf (
                        "WARNING: Could not converge (%d of %d) dEps=%g\n",
                                itDiag+1, maxItDiag, eps-epsPrev);
                     }
                     nConvIt++;
                     break;
                  }
                  cout << step << ": Eps: " << eps*HA2EV
                     << ", dEps: " << eps - epsPrev << endl;
                  epsPrev = eps;
               } // itDiag


               if (autoSteps && (itDiag + 1 < maxItDiag))  {
                  // we spared some CCG steps with this state that we can use
                  // for higher states
                  spareSteps += maxItDiag - itDiag - 1;
               }

            } // i
                  
            // --- print convergence statistics
            cout << SX_SEPARATOR;
            cout << "| Convergence statistics:\n";
            cout << "|    absolute convergence:     " << nConvAbs << " times\n";
            cout << "|    small residue:            " << nConvRes << " times\n";
            cout << "|    numerical limit reached:  " << nCosAcc  << " times\n";
            cout << "|    exceeded max. iterations: " << nConvIt  << " times\n";
            cout << SX_SEPARATOR;

            if (subDiagSize > 0)
            {
               // --- subspace diagonalization ---
               SxDiracVec<TPrecEps> epsDiag;
               epsDiag = subspaceDiagonalization (&waves(iSpin,ik),
                                                  subDiagSize, subDiagOvlp);
               // --- copy energies and count number of unconverged states
               SxDiracVec<TPrecEps>::Iterator epsDiagIt, epsFermiIt;
               epsDiagIt = epsDiag.begin (),
               epsFermiIt = fermi.eps(iSpin,ik).begin ();
               for (i = 0; i <= convergeLimit; ++i)
                  if (fabs(*epsDiagIt - *epsFermiIt) > dEps) nUnconvDiag++;
               fermi.eps(iSpin, ik) = epsDiag;
            }


            // --- print results
            eBand  = fermi.eps(iSpin,ik).sum (0, convergeLimit);
            sxprintf ("ik=%d eBand(%d)=%15.12f\n", ik + 1, it+1, eBand);
            fflush (stdout);
            SxFileIO::appendToFile
            (    SxString(it)                                + "\t" // iteration
               + SxString(GETTIME(Timer::ElMinim), "%12.8f") + "\t" // time
               + SxString(eBand,  "%15.12f")                 + "\n" // E_band
            , "energy.dat");


            // --- convergence?
            cout << "eBand: " << eBand << endl;
            cout << "lastEnergy: " << lastEnergy << endl;
            cout << "eBand - lastEnergy: " << eBand - lastEnergy << endl;
            if ( fabs (eBand - lastEnergy) < dEnergy )  {
               sxprintf ("Convergence reached for (ik=%d, iSpin=%d).\n",
                       ik+1, iSpin+1);
               break;
            }
/*            if ( nUnconvDiag == 0 && nConvIt <= nSloppy)  {
               sxprintf ("Convergence reached for (ik=%d, iSpin=%d): "
                       "all eigenvalues converged to dEps.\n",
                       ik+1, iSpin+1);
               break;
            }*/
            if ( it+1 >= maxSteps )  {
               sxprintf ("WARNING: Maximum number of steps exceeded.\n");
               sxprintf ("         Convergence not yet reached.\n");
               break;
            }
            lastEnergy = eBand;
            hamPtr->printEnergies();
            cout << "eps; " << (psiI ^ hPsi).chop().re << endl;

            if ( !(it % printSteps) )  {
               fermi.epsSortList(iSpin, ik) 
                  = fermi.eps(iSpin, ik).getSortIdx ();
               cout << SX_SEPARATOR;
               fermi.printK (iSpin, ik, false);
               cout << SX_SEPARATOR;
            }
         }
         // --- print out
         fermi.epsSortList(iSpin, ik) = fermi.eps(iSpin, ik).getSortIdx ();
         cout << SX_SEPARATOR;
         fermi.printK (iSpin, ik, false, true);
         cout << SX_SEPARATOR;
      }
   }
      
   cout << "####################" << endl;
   cout << "steps needed: " << step << endl;
   cout << "####################" << endl;
}

void SxHamSolver::bandStructureCG (const SxSymbolTable *cmd, bool calc)
{
   SX_NO_MPI;
   SX_CLOCK (Timer::ElMinim);
   SX_CHECK (cmd);

   int maxStepsCCG = 5, printSteps = 20;
   int subDiagSize = -1, subDiagOvlp = -1;
   int nSloppy = 0;
   bool verbose = false, printRes = false;
   double dRes = -1.;
   bool autoSteps = false;
   int spareSteps;
   enum WavesGuess { None, Propagate } guess = None;
   bool keepOccFixed = false;
   try  {
      if (cmd->contains("dEps"))
         dEps = cmd->get("dEps")->toReal ();
      else
         dEps = 1e-4 / HA2EV;
      if (cmd->contains("dEnergy"))
         dEnergy   = cmd->get("dEnergy")->toReal ();
      else
         dEnergy   = dEps * nStates;
      if (cmd->contains("maxStepsCCG"))
         maxStepsCCG   = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("maxSteps"))
         maxSteps      = cmd->get("maxSteps")->toInt ();
      else
         maxSteps      = 20;
      if (cmd->contains("printSteps"))
         printSteps    = cmd->get("printSteps")->toInt();
      if (cmd->contains("autoSteps"))
         autoSteps = cmd->get("autoSteps")->toAttribute ();
      if (cmd->contains("verbose"))
         verbose = cmd->get("verbose")->toAttribute ();
      if (cmd->contains("printResidue"))
         printRes = cmd->get("printResidue")->toAttribute ();
      if (cmd->contains("dPsi"))
         dRes = cmd->get("dPsi")->toReal ();
      if (cmd->contains("nSloppy"))
         nSloppy = cmd->get("nSloppy")->toInt ();
      if (cmd->containsGroup ("subspaceDiag"))  {
         SxSymbolTable *subspaceGroup = cmd->getGroup ("subspaceDiag");
         subDiagSize = subspaceGroup->get("maxSize")->toInt ();
         subDiagOvlp = subspaceGroup->get("overlap")->toInt ();
         if (subDiagSize <= subDiagOvlp)  {
            cout << SxString(
                    "subspace block overlap must be (even considerably) "
                    "smaller than maximum block size."
                    ).wrap () << endl;
            SX_QUIT;
         }
      }
      if (cmd->contains("propagateWaves"))
         guess = cmd->get("propagateWaves")->toAttribute () ? Propagate : None;

      if (cmd->contains("keepOccFixed"))
         keepOccFixed = cmd->get("keepOccFixed")->toAttribute ();
           
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   // --- print input parameter
   cout << SX_SEPARATOR;
   cout << "| Band structure state-by-state conjugate gradient\n";
   cout << SX_SEPARATOR;
   cout << "|   max. steps:                " << maxSteps << endl;
   cout << "|   details printed every:     " << printSteps << " iteration(s)\n";
   cout << "|   dEnergy:                   " << dEnergy  << " H\n";
   cout << "|   sloppy (may not converge): " << nSloppy << " states\n";
   cout << "|   initial waves:             ";
   if (guess == Propagate)
      cout << "from previous k-point" << endl;
   else
      cout << "from initialGuess" << endl;
   cout << "|   conjugate gradient:\n";
   cout << "|      max. steps per state:   " << maxStepsCCG << endl;
   cout << "|      final stage tuning:     " << (autoSteps ? "on" : "off") 
                                             << endl;
   cout << "|      energy convergence      " << dEps << " H\n";
   cout << "|      residue convergence :   ";
   if (dRes > 0.)
      cout << dRes << endl;
   else
      cout << "off" << endl;
   cout << "|   subspace diagonalization:  ";
   if (subDiagSize > 0)  {
      cout << endl;
      cout << "|      max. blocksize:         " << subDiagSize << endl;
      cout << "|      block overlap:          " << subDiagOvlp << endl;
   } else {
      cout << "all states" << endl;
   }
   if (keepOccFixed)
      cout << "|   No occupation number update at the end." << endl;
   cout << SX_SEPARATOR;

   if (!calc) return;

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   int ik, nk         = waves.getNk ();
   int iSpin;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   int i;
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   int it, itDiag;
   int nConvRes;
   int maxItDiag;
   SxDiracVec<TPrecCoeffG> psiI, hPsi, dPsi, dPsiOld, K, X, Xold, hX;
   SxDiracVec<TPrecCoeffG> psiGI, dPsiGI;
   SxDiracMat<TPrecCoeffG> h;
   SxDiracSymMat<TPrecCoeffG> hSym;
   SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;
   Real8 eps, gamma = 0., theta, epsPrev=0.;
   PrecCoeffG scp, tr, trOld = 0.;
   SxMatrix<TPrecCoeffG> M(2,2), U(2,2);

   double eBand, lastEnergy;

   int convergeLimit, nUnconvDiag;

   SxOverlap S = hamPtr->getS ();
   std::streamsize oldPrec = cout.precision ();
   cout.precision (10);
   for (ik=0; ik < nk; ik++)  {
      if (guess == Propagate && ik > 0)  {
         propagateWaves (ik-1, ik);
         for (iSpin = 0; iSpin < nSpin; ++iSpin) {
            fermi.eps(iSpin, ik) 
               = subspaceDiagonalization (&waves(iSpin,ik), subDiagSize,
                                          subDiagOvlp);
         }
      }
      const SxGBasis &gk = waves.getGkBasis ()(ik);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         convergeLimit = nStates - 1 - nSloppy;
         lastEnergy = fermi.eps(iSpin,ik).sum (0, convergeLimit);
         for (it = 0; it < maxSteps; it++)  {
            nConvIt = 0;
            nConvAbs = 0;
            nConvRes = 0;
            nCosAcc = 0;
            nUnconvDiag = 0;
            spareSteps = 0;
            
            // --- conjugate gradient ------
            for (i=0; i < nStates; i++)  {
               psiI = waves(i,iSpin,ik);
               

               if (verbose)  {
                  sxprintf ("state: (%d,%d,%d) ", i, iSpin, ik); fflush(stdout);
               }

               maxItDiag = maxStepsCCG;
               if (autoSteps && spareSteps > 0)  {
                  if (spareSteps > maxStepsCCG)  {
                     maxItDiag += maxStepsCCG;
                     spareSteps -= maxStepsCCG;
                  } else  {
                     maxItDiag += spareSteps;
                     spareSteps = 0;
                  }
               }
               epsPrev = fermi.eps(i, iSpin, ik);
               for (itDiag=0; /* empty */; itDiag++)  {

                  // orthogonalize psiI to lower states
                  if (i > 0)
                     S.setOrthogonal (&psiI, waves.getBlock(0, i, iSpin, ik));
                  S.normalize (&psiI);
                  
                  // --- compute gradient: |dPsi> = (H-e) |psi>
                  hPsi = *hamPtr | psiI;
                  eps  = dot(gk | psiI, hPsi).re;
                  dPsi = hPsi - eps * (S * psiI);
                  fermi.eps(i,iSpin,ik) = eps;

                  if (printRes)  {
                     cout << "(" << i << "," << iSpin << "," << ik << ").";
                     cout << itDiag << ": res = " << dot(dPsi, dPsi).re;
                     cout << "; dEps = " << (eps - epsPrev) << endl;
                  }

                  // --- convergence on dEps
                  if (fabs(eps - epsPrev) < dEps)  {
                     if (itDiag > 0)  {
                        if (verbose)
                           sxprintf ("conv due to abs eps (%d of %d)\n", 
                                   itDiag+1, maxItDiag);
                        nConvAbs++;
                        break;
                     }
                     // --- for first iteration, residue provides a good
                     // indicator
                     if (dRes > 0. && dot(dPsi, dPsi).re < dRes)  {
                        if (verbose)
                           sxprintf ("conv due to residue (%d of %d) dEps=%g\n",
                                   itDiag+1, maxItDiag,
                                   eps - epsPrev);
                        nConvRes++;
                        break;
                     }
                  }

                  // --- orthogonalize preconditioned gradient
                  K = hamPtr->preconditioner (psiI) * dPsi;

                  // --- get conjugate direction, ref1 (8.14b) / ref2 (10.6.7)
                  tr = (dPsi ^ K).chop();
                  if (itDiag > 0)  gamma = (tr - (dPsiOld^K).chop()) / trOld;
   //             if (itDiag > 0)  gamma = tr / trOld;
                  dPsiOld.copy (dPsi);
                  trOld = tr;

                  // --- compute search direction, ref1 (8.14a)
                  if (itDiag == 0)  X = K * 1.;
                  else              X = K + gamma * Xold;
                  Xold.resize (0); Xold.copy (X);
                  X = (*psiI.getBasisPtr ()) | X;
                  S.setOrthogonal (&X, waves.getBlock (0, i+1, iSpin, ik));
                  { SX_CLOCK (Timer::ortho);
                    scp = (psiI | S | X);
   //               X -= psiI * scp;
                    X.plus_assign_ax (-scp, psiI);
                  }
                  S.normalize (&X);
                  // TODO: ugly
                  X.handle->auxData.i     = i;
                  X.handle->auxData.iSpin = iSpin;
                  X.handle->auxData.ik    = ik;

                  // --- get M matrix (in ref1: H matrix)
                  hX = *hamPtr | X;
                  M(0,0) = eps;            M(0,1) = dot(gk | X, hPsi);
                  M(1,0) = M(0,1).conj();  M(1,1) = dot(gk | X, hX);

                  // --- uniform transformation U, ref1 (8.16)
                  if ( M(0,0).re - M(1,1).re == 0. )  // '==' is correct!
                  {
                     /* This is a workaround for the numerical limit of
                        degenerate eigenstates. If an exact eigenstate has
                        been found, and X is a different eigenstate to the
                        same energy, the formula below doesn't work, as
                        atan(1/0) is undefined. Any theta is fine, so
                        we just take 0 to get out at the numerical limit
                        test. */
                     theta = 0.;
                  }  else  {
                     theta = 0.5 * atan (  (M(1,0).re + M(0,1).re) 
                                         / (M(0,0).re - M(1,1).re) );
                     // --- in case of being far off the quadr. regime theta
                     //     points to the maximum rather than to the minimum!  
                     if (M(0,0).re > M(1,1).re)  theta -= PI/2.;
                  }

                  SX_CHECK_NUM (theta);
                  if (fabs(theta) < 1e-8)  {
                     if (verbose)
                        sxprintf ("Numerical limit reached (%d of %d)\n",
                                itDiag+1, maxItDiag);
                     nCosAcc++;
                     break;
                  }
                  U(0,0) = cos (theta); U(0,1) = sin (theta);
                  U(1,0) = -U(0,1);     U(1,1) = U(0,0);

                  // --- compute new wavefunction, ref1 (8.15b)
                  waves(i,iSpin,ik) <<= U(0,0) * psiI + U(0,1) * X;

                  // --- convergence?
                  if (itDiag+1 >= maxItDiag)  {
                     if (verbose)  {
                        sxprintf (
                        "WARNING: Could not converge (%d of %d) dEps=%g\n",
                                itDiag+1, maxItDiag, eps-epsPrev);
                     }
                     nConvIt++;
                     break;
                  }
                  epsPrev = eps;
                 
               } // itDiag

               if (autoSteps && (itDiag + 1 < maxItDiag))  {
                  // we spared some CCG steps with this state that we can use
                  // for higher states
                  spareSteps += maxItDiag - itDiag - 1;
               }

            } // i
                  
            // --- print convergence statistics
            cout << SX_SEPARATOR;
            cout << "| Convergence statistics:\n";
            cout << "|    absolute convergence:     " << nConvAbs << " times\n";
            cout << "|    small residue:            " << nConvRes << " times\n";
            cout << "|    numerical limit reached:  " << nCosAcc  << " times\n";
            cout << "|    exceeded max. iterations: " << nConvIt  << " times\n";
            cout << SX_SEPARATOR;

            
            {
               // --- subspace diagonalization ---
               SxDiracVec<TPrecEps> epsDiag;
               epsDiag = subspaceDiagonalization (&waves(iSpin,ik),
                                                  subDiagSize, subDiagOvlp);
               // --- copy energies and count number of unconverged states
               SxDiracVec<TPrecEps>::Iterator epsDiagIt, epsFermiIt;
               epsDiagIt = epsDiag.begin (),
               epsFermiIt = fermi.eps(iSpin,ik).begin ();
               for (i = 0; i <= convergeLimit; ++i)
                  if (fabs(*epsDiagIt - *epsFermiIt) > dEps) nUnconvDiag++;
               fermi.eps(iSpin, ik) = epsDiag;
            }


            // --- print results
            eBand  = fermi.eps(iSpin,ik).sum (0, convergeLimit);
            sxprintf ("ik=%d eBand(%d)=%15.12f\n", ik + 1, it+1, eBand);
            fflush (stdout);
            SxFileIO::appendToFile
            (    SxString(it)                                + "\t" // iteration
               + SxString(GETTIME(Timer::ElMinim), "%12.8f") + "\t" // time
               + SxString(eBand,  "%15.12f")                 + "\n" // E_band
            , "energy.dat");


            // --- convergence?
            if ( fabs (eBand - lastEnergy) < dEnergy )  {
               sxprintf ("Convergence reached for (ik=%d, iSpin=%d).\n",
                       ik+1, iSpin+1);
               break;
            }
            if ( nUnconvDiag == 0 && nConvIt <= nSloppy)  {
               sxprintf ("Convergence reached for (ik=%d, iSpin=%d): "
                       "all eigenvalues converged to dEps.\n",
                       ik+1, iSpin+1);
               break;
            }
            if ( it+1 >= maxSteps )  {
               sxprintf ("WARNING: Maximum number of steps exceeded.\n");
               sxprintf ("         Convergence not yet reached.\n");
               break;
            }
            lastEnergy = eBand;

            if ( !(it % printSteps) )  {
               fermi.epsSortList(iSpin, ik) 
                  = fermi.eps(iSpin, ik).getSortIdx ();
               cout << SX_SEPARATOR;
               fermi.printK (iSpin, ik, false);
               cout << SX_SEPARATOR;
            }
         }
         // --- print out
         fermi.epsSortList(iSpin, ik) = fermi.eps(iSpin, ik).getSortIdx ();
         cout << SX_SEPARATOR;
         fermi.printK (iSpin, ik, false, true);
         cout << SX_SEPARATOR;
      }
   }
   cout.precision (oldPrec);

   if (!keepOccFixed) fermi.fermiDistribution (ekt);

}

void SxHamSolver::rmmDiis (const SxSymbolTable *cmd, bool sloppyConv, 
                           bool printOnly)
{
   SX_NO_MPI;
   SX_CLOCK (Timer::ElMinim);

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   int nk = waves.getNk ();
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());

   int maxItDiag = 4;
   bool verbose = false;
   dRelR = 0.3;
   try  {
      verbose = cmd->contains("verbose")
              ? cmd->get("verbose")->toAttribute ()
              : false;
      if (cmd->contains("maxStepsCCG"))
         maxItDiag = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("dRelR"))
         dRelR = cmd->get("dRelR")->toReal ();
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   if (printOnly)  {
      cout << "|   RMM-DIIS:\n";
      cout << "|      max. steps per state:   " << maxItDiag << endl;
      cout << "|      rel. change in residue: " << (dRelR) * 100. << "%\n";
      return;
   }

   nConvAbs = nCosAcc = nConvIt = nConvOcc = nConvRel = 0;
   nStepsAbs = nStepsCos = nStepsIt = nStepsOcc = nStepsRel = 0;
   SxOverlap S = hamPtr->getS ();
   for (int ik=0; ik < nk; ik++)  {
      for (int iSpin=0; iSpin < nSpin; iSpin++)  {

         // --- subspace diagonalization (ref 3, IIIC)
         fermi.eps(iSpin, ik) = subspaceDiagonalization (&waves(iSpin,ik));

         // --- DIIS for each band
         for (int i=0; i < nStates; i++)  {
            rmmDiis (i, iSpin, ik, maxItDiag, verbose, sloppyConv);
         }

         // --- orthonormalization
         S.orthonormalize (&waves(iSpin, ik));
         
      }
   }

   // print convergence statistics
   printConvergenceStatistics ();
   cout << SX_SEPARATOR;
}

enum RmmTimer { DiagRmm, Reduce};
SX_REGISTER_TIMERS(RmmTimer)
{
   regTimer (DiagRmm, "RMM-DIIS: diag");
   regTimer (Reduce , "RMM-DIIS: reduce");
}

SxDiracVec<TPrecCoeffG> 
SxHamSolver::solveRmmDiis (const SxDiracMat<TPrecCoeffG> &h2In,
                           const SxDiracMat<TPrecCoeffG> &hsshIn,
                           const SxDiracMat<TPrecCoeffG> &s2In,
                           const SxDiracMat<TPrecCoeffG> &s,
                           double *epsPtr,
                           int n)
{
   SX_CLOCK(Timer::RmmDiis);
   SX_CHECK(epsPtr);

   SxDiracMat<TPrecCoeffG> Rij, h2, hssh, s2;
   SxDiracMat<TPrecCoeffG> L = s.choleskyDecomposition ().inverse (),
                           U = L.adjoint ();
   SxDiracMat<TPrecCoeffG>::Eigensystem eig;
   SxDiracVec<TPrecCoeffG> alpha;
   h2   = L ^ h2In   ^ U;
   hssh = L ^ hsshIn ^ U;
   s2   = L ^ s2In   ^ U;

   int nPhi = (int)s2.nCols ();
   SxDiracVec<TPrecCoeffG> beta;
   beta.reformat (nPhi, n);

   for (int i = 0; i < n; ++i)  {
      double &eps = epsPtr[i];
      double dEps = -1., dEpsOld = -1., epsIn = eps, epsInOld;
      for (int it = 0; it < 30; ++it)  {
         Rij = h2 - eps * hssh + (eps*eps)*s2;
         {
            SX_CLOCK (DiagRmm);
            eig = Rij.eigensystem ();
         }
         epsInOld = epsIn;
         epsIn = eps;

         alpha = eig.vecs.colRef(0);
         eps = dot(alpha, hssh ^ alpha) / (2. * dot(alpha, s2 ^ alpha));
         dEpsOld = dEps;
         dEps = eps - epsIn;
         //cout << "rmmDIIS solver: dEps = " << dEps
         //     << "; R=" << eig.vals//(0) 
         //     << endl;
         if (fabs(dEps) < 1e-14 * fabs(eps)) break;
         if (fabs(epsInOld - epsIn) < 1e-14 * fabs(eps)) break;
         if (fabs(dEps) < 1e-16) break;

         if (it > 0)  {
            // predict optimal eps
            // assume that dEps is linear in eps
            double k = - (dEpsOld - dEps) / (epsInOld - epsIn);
            //cout << "k=" << k << endl;
            eps = epsIn + dEps / k;
         }
      }
      beta.colRef(i) <<= U ^ alpha;
      if (i + 1 >= n) break;
      // --- remove current direction from search space
      SX_CLOCK (Reduce);
      SxDiracMat<TPrecCoeffG> T, Ta;
      T = eig.vecs(SxIdx(nPhi, nPhi*nPhi-1));
      T.reshape (nPhi, nPhi-1);
      Ta = T.adjoint ();
      h2   = Ta ^ h2 ^ T;
      hssh = Ta ^ hssh ^ T;
      s2   = Ta ^ s2 ^ T;
      U = U ^ T;
      nPhi--;
   }
   return beta;
}

void SxHamSolver::rmmDiis (int i, int iSpin, int ik,
                           int maxItDiag, bool verbose, bool sloppyConv)
{
   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   SX_CHECK (ik >= 0 && ik < waves.getNk (), ik, waves.getNk ()); 
   SX_CHECK (iSpin >= 0 && iSpin < nSpin, iSpin, nSpin);
   SX_CHECK (i >= 0 && i < nStates, i, nStates);

   const SxGBasis &gk = waves.getGkBasis ()(ik);
   const SxBasis &waveBasis = waves.getBasis (ik);
   SxOverlap S = hamPtr->getS ();
   SxDiracVec<TPrecCoeffG> psiI, dPsi, hPsi, sPsi, X;

   psiI = waves(i,iSpin,ik);
   if (verbose)  {
      sxprintf ("state: (%d,%d,%d) ", i, iSpin, ik); fflush(stdout);
   }

   double eps, epsPrev = 0.;
   double resNorm = -1., resNorm0;

   PsiGI phiList, hList, sList;

   // --- compute gradient: |dPsi> = (H-e) |psi>
   hPsi = *hamPtr | psiI;
   eps  = dot (gk | psiI, hPsi).re;
   dPsi = hPsi - eps * (sPsi = S * psiI);
   fermi.eps(i,iSpin,ik) = eps;
   SxDiracAux auxData = dPsi.handle->auxData;

   phiList.copy (psiI);
   hList.copy (hPsi);
   sList.copy (sPsi);

   resNorm0 = dPsi.normSqr ();
   if (verbose)  {
      cout << "R=" << resNorm0 << endl;
   }
   for (int itDiag=0; /* empty */; itDiag++)  {

      if (itDiag >= maxItDiag)  {
         if (verbose)
            sxprintf ("WARNING: Could not converge (%d of %d) dEps=%g\n",
                    itDiag, maxItDiag, eps-epsPrev);
         nConvIt++;
         nStepsIt += itDiag;
         break;
      }

      // --- preconditioned gradient and search direction
      X = waveBasis | (hamPtr->preconditioner (gk | psiI) * dPsi);
      S.normalize (&X);

      // --- add to search space
      phiList = mergeCols (phiList, X);
      hList   = mergeCols (hList, *hamPtr | X);
      sList   = mergeCols (sList, S * X);

      epsPrev = eps;

      // --- solve non-linear DIIS problem in search space
      SxDiracVec<TPrecCoeffG> alpha;
      {
         SxDiracMat<TPrecCoeffG> h2, hssh, s2, s;
         s = sList.overlap(gk | phiList);
         if (s.eigenvalues ().real ().minval () < 1e-16)  {
            if (verbose)
               sxprintf ("numerical limit reached (%d of %d)\n", 
                       itDiag+1, maxItDiag);
            nCosAcc++;
            nStepsCos += itDiag + 1;
            break;
         }
         h2 = hList.overlap (hList);
         hssh = hList.overlap (sList);
         hssh += hssh.adjoint ();
         s2 = sList.overlap (sList);

         alpha = solveRmmDiis (h2, hssh, s2, s, &eps);
      }

      // --- update wave
      psiI(SxIdx(0,(int)psiI.getSize () - 1)) <<= phiList ^ alpha;
      fermi.eps(i,iSpin,ik) = eps;

      if (itDiag > 0 && fabs(eps - epsPrev) < dEps)  {
         if (verbose)
            sxprintf ("conv due to abs eps (%d of %d)\n", itDiag+1, maxItDiag);
         nConvAbs++;
         nStepsAbs += itDiag + 1;
         break;
      }

      // --- update gradient
      hPsi = hList ^ alpha;
      dPsi = hPsi - eps * (sList ^ alpha);
      dPsi.handle->auxData = auxData;

      if (verbose)  {
         double epsH = dot (gk | psiI, hPsi);
         cout << "eps(DIIS)=" << eps
              << " eps(H)=" << epsH
              << " delta=" << (eps-epsH)
              << endl;
      }
         
      resNorm = dPsi.normSqr ();
      if (verbose)  {
         cout << "R=" << resNorm << endl;
      }

      // --- convergence checks
      if ((resNorm < dRelR * resNorm0 && itDiag > 0)
          || resNorm < min(0.03, dRelR*dRelR) * resNorm0)  {
         if (verbose)
            sxprintf ("conv due to rel norm (%d of %d)\n",itDiag+1,maxItDiag);
         nConvRel++;
         nStepsRel += itDiag + 1;
         break;
      }

      if (sloppyConv && itDiag >= 1 && fabs(fermi.focc(i,iSpin,ik)) < 1e-8)
      {
         if (verbose)
            sxprintf ("conv. due to low occ. (%d of %d)\n", itDiag+1, maxItDiag);
         nConvOcc++;
         nStepsOcc+= itDiag + 1;
         break;
      }

   }

   if (verbose)  {
      cout << "R(final)=" << resNorm 
           << " factor=" << (resNorm/resNorm0)
           << endl;
   }

}

void SxHamSolver::blockRmmDiis (const SxSymbolTable *cmd, bool sloppyConv, 
                                bool printOnly)
{
   SX_NO_MPI;
   SX_CLOCK (Timer::ElMinim);

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SxOverlap S = hamPtr->getS ();
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());

   int maxItDiag = 4;
   bool verbose = false;
   int defBlockSize = 64;
   dRelR = 0.03;
   try  {
      verbose = cmd->contains("verbose")
              ? cmd->get("verbose")->toAttribute ()
              : false;
      if (cmd->contains("maxStepsCCG"))
         maxItDiag = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("blockSize"))
         defBlockSize = cmd->get("blockSize")->toInt ();
      if (cmd->contains("dRelR"))
         dRelR = cmd->get("dRelR")->toReal ();
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   if (printOnly)  {
      cout << "|   block RMM-DIIS:\n";
      cout << "|      block size:             " << defBlockSize << endl;
      cout << "|      max. steps per state:   " << maxItDiag << endl;
      cout << "|      rel. change in residue: " << (dRelR) * 100. << "%\n";
      return;
   }

   nConvAbs = nCosAcc = nConvIt = nConvOcc = nConvRel = 0;
   nConvCycle = 0;
   nStepsAbs = nStepsCos = nStepsIt = nStepsOcc = nStepsRel = 0;

   if (verbose)
      cout << "abs dEps for this step: " << dEps << endl;
   for (int ik=0; ik < waves.getNk (); ik++)  {
      for (int iSpin=0; iSpin < nSpin; iSpin++)  {

         // --- subspace diagonalization (ref 3, IIIC)
         fermi.eps(iSpin, ik) = subspaceDiagonalization (&waves(iSpin,ik));

         // --- DIIS for block of bands
         int blockSize;
         for (int i=0; i < nStates; i+= blockSize)  {
            blockSize = defBlockSize;
            if (i + blockSize >= nStates)  {
               blockSize = nStates-i;
            } else  {
               // try not to cut through degenerate states. Physically not
               // necessary, but improves stability wrt numerical noise
               while (   blockSize > 1 
                      && blockSize > defBlockSize - 3
                      && fabs( fermi.eps(i+blockSize-1, iSpin, ik) 
                              -fermi.eps(i+blockSize  , iSpin, ik)) < 1e-7)
                  blockSize--; 
            }
            blockRmmDiis (i, iSpin, ik, blockSize, maxItDiag, verbose, sloppyConv);
         }

         // orthonormalization
         S.orthonormalize (&waves(iSpin, ik));
      }
   }

   printConvergenceStatistics ();
   cout << SX_SEPARATOR;
}


void SxHamSolver::blockRmmDiis (int i, int iSpin, int ik, int blockSize, 
                                int maxItDiag, bool verbose, bool sloppyConv)
{
   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   SX_CHECK (ik >= 0 && ik < waves.getNk (), ik, waves.getNk ()); 
   SX_CHECK (iSpin >= 0 && iSpin < nSpin, iSpin, nSpin);

   const SxGBasis &gk = waves.getGkBasis ()(ik);
   const SxBasis &waveBasis = waves.getBasis (ik);
   SxOverlap S = hamPtr->getS ();
   SxDiracVec<TPrecCoeffG> psiI, dPsi, hPsi, sPsi, X;
   SxArray<bool> converged(blockSize);
   SxVector<Double> resNorm(blockSize), resNorm0(blockSize);
   SxDiracVec<Double> eps(blockSize), epsPrev(blockSize);
   converged.set (false);

   psiI = waves.getBlock(i, blockSize, iSpin,ik);
   eps  = fermi.eps(iSpin,ik)(SxIdx(i,i+blockSize-1));
   if (verbose)  {
      sxprintf ("states: (%d..%d,%d,%d)\n", i+1, i+blockSize, iSpin+1, ik+1);
      fflush(stdout);
   }
   

   PsiGI phiList, hList, sList;

   // --- compute gradient: |dPsi> = (H-e) |psi>
   hPsi = *hamPtr | psiI;
   sPsi = S * psiI;
   SxDiracAux auxData = hPsi.handle->auxData;
   {
      SX_CLOCK (Timer::SearchDir);
      dPsi.copy (hPsi);
      for (int iBlock = 0; iBlock < blockSize; ++iBlock)  {
         dPsi.colRef(iBlock).plus_assign_ax (-eps(iBlock), sPsi.colRef(iBlock));
      }
      dPsi.handle->auxData = auxData;
   }

   phiList = psiI.getCopy ();
   hList.copy (hPsi);
   sList.copy (sPsi);

   int nConvBlock = 0;

   // --- compute initial residue norm
   for (int iBlock = 0; iBlock < blockSize; ++iBlock)
      resNorm0(iBlock) = dPsi.colRef(iBlock).normSqr ();
   if (verbose)  {
      cout << "R=" << resNorm0 << endl;
   }
   resNorm <<= resNorm0;

   for (int itDiag=0; /* empty */; itDiag++)  {

      if (itDiag >= maxItDiag)  {
         for (int iBlock = 0; iBlock < blockSize; ++iBlock)  {
            if (converged(iBlock)) continue;
            if (verbose)  {
               sxprintf ("state %d: could not converge (%d of %d) R=%g\n",
                       i+iBlock+1, itDiag, maxItDiag, resNorm(iBlock));
            }
            nConvIt++;
            nConvBlock++;
            nStepsIt += itDiag;
         }
         break;
      }

      // count number of unnecessary cycles for converged states
      nConvCycle += nConvBlock;

      // --- preconditioned gradient and search direction
      {
         SX_CLOCK (Timer::SearchDir);
         SxDiracMat<TPrecCoeffG> K (gk.ng, blockSize);
         K.handle->auxData = dPsi.handle->auxData;
         for (int iBlock = 0; iBlock < blockSize; ++iBlock)  {
            K.colRef (iBlock) <<= 
               hamPtr->preconditioner (gk | psiI.colRef(iBlock))
               * dPsi.colRef (iBlock);
         }
         X = waveBasis | K;
         S.normalize (&X);
      }

      // --- add to search space
      phiList = mergeCols (phiList, X);
      hList   = mergeCols (hList, *hamPtr | X);
      sList   = mergeCols (sList, S * X);

      epsPrev.copy (eps);

      // --- solve non-linear DIIS problem in search space
      SxDiracVec<TPrecCoeffG> alpha;
      {
         SxDiracMat<TPrecCoeffG> h2, hssh, s2, s;
         s = sList.overlap(gk | phiList);
         SxDiracVec<Double> sev = s.eigenvalues ().real ();
         double cut = 1e-14 * sev.maxval ();
         if (sev.minval () < cut)  {
            // --- remove linear dependencies
            int nSingular, nPhi = (int)sev.getSize ();
            for (nSingular = 0; nSingular < nPhi; ++nSingular)
               if (sev(nSingular) > cut) break;
            
            cout << "Sij ev:" << sev << endl;
            if (nSingular >= blockSize)  {
               if (verbose)
                  sxprintf ("numerical limit reached (%d of %d)\n", 
                          itDiag+1, maxItDiag);
               nCosAcc+=blockSize;
               nStepsCos += (itDiag + 1) * blockSize;
               break;
            }
            if (verbose)  {
               cout << "Omitting " << nSingular 
                    << " linear dependent directions." << endl;
            }

            SxDiracMat<TPrecCoeffG>::Eigensystem eig = s.eigensystem ();
            SxDiracVec<TPrecCoeffG> U;
            U = eig.vecs (SxIdx(nPhi * nSingular, nPhi*nPhi - 1));
            U.reshape (nPhi, nPhi - nSingular);
            for (int iPhi = nSingular; iPhi < nPhi; ++iPhi)  {
               eig.vecs.colRef (iPhi) *= 1./sqrt(eig.vals(iPhi));
            }
            phiList = phiList ^ U;
            phiList.handle->auxData = psiI.handle->auxData;
            hList = hList ^ U;
            sList = sList ^ U;
            hList.handle->auxData = auxData;
            sList.handle->auxData = auxData;
            s = U.adjoint () ^ s ^ U;
         }
         // --- setup matrix representation
         SX_START_TIMER(Timer::SubspaceMatrix);
         h2 = hList.overlap (hList);
         hssh = hList.overlap (sList);
         hssh += hssh.adjoint ();
         s2 = sList.overlap (sList);
         SX_STOP_TIMER(Timer::SubspaceMatrix);

         alpha = solveRmmDiis (h2, hssh, s2, s, eps.elements, blockSize);
      }

      SX_START_TIMER(Timer::WaveUpdate);
      // --- update waves
      psiI(SxIdx(0,(int)psiI.getSize () - 1)) <<= phiList ^ alpha;

      // --- update gradient
      hPsi = hList ^ alpha;
      sPsi = PsiG ();
      sPsi = sList ^ alpha;
      dPsi.copy (hPsi);
      for (int iBlock = 0; iBlock < blockSize; ++iBlock)
         dPsi.colRef(iBlock).plus_assign_ax(-eps(iBlock), sPsi.colRef(iBlock));
      dPsi.handle->auxData = auxData;
      SX_STOP_TIMER(Timer::WaveUpdate);

      // --- compute initial residue norm
      for (int iBlock = 0; iBlock < blockSize; ++iBlock)
         resNorm(iBlock) = dPsi.colRef(iBlock).normSqr ();
      if (verbose)  {
         cout << "R=" << resNorm << endl;
         cout << "dEps=" << (eps-epsPrev) << endl;
      }

      // --- convergence checks
      for (int iBlock = 0; iBlock < blockSize; ++iBlock)  {
         if (converged(iBlock)) continue;

         if (fabs(eps(iBlock) - epsPrev(iBlock)) < dEps)  {
            if (verbose)
               sxprintf ("state %d: conv due to abs eps (%d of %d)\n", 
                       i+iBlock+1, itDiag+1, maxItDiag);
            converged(iBlock) = true;
            nConvBlock++;
            nConvAbs++;
            nStepsAbs += itDiag + 1;
            continue;
         }

         if (resNorm(iBlock) < dRelR * resNorm0(iBlock))  {
            if (verbose)  {
               sxprintf ("state %d: conv due to rel norm (%.3g, %d of %d)\n",
                       i+iBlock+1, resNorm(iBlock)/resNorm0(iBlock),
                       itDiag+1, maxItDiag);
            }

            converged(iBlock) = true;
            nConvBlock++;
            nConvRel++;
            nStepsRel += itDiag + 1;
            continue;
         }

         if (sloppyConv && itDiag >= 1
             && fabs(fermi.focc(i+iBlock,iSpin,ik))<1e-8)
         {
            if (verbose)
               sxprintf ("state %d: conv. due to low occ. (%d of %d)\n",
                       i+iBlock+1, itDiag+1, maxItDiag);
            converged(iBlock) = true;
            nConvBlock++;
            nConvOcc++;
            nStepsOcc+= itDiag + 1;
            //continue;
         }
      }

      if (nConvBlock == blockSize) { break; }

      if (itDiag + 1 >= maxItDiag) continue;
      if (blockSize - nConvBlock < min(8, blockSize / 4)
          || blockSize - nConvBlock < 2)
      {
         // only few states are not converged
         // use single-state algorithm for these
         for (int iBlock = 0; iBlock < blockSize; ++iBlock)
            if (!converged(iBlock))
               rmmDiis (i + iBlock, iSpin, ik, maxItDiag - itDiag, verbose,
                        sloppyConv);
         break;
      }
      if (nConvBlock > 0)  {
         // --- remove lowest converged states from block

         // find number of converged states at bottom
         int nConvBottom, top;
         // nConvBottom is the first unconverged state within block
         for (nConvBottom = 0; nConvBottom < blockSize; ++nConvBottom)
            if (!converged(nConvBottom)) break;
         // top is the lowest converged state from the top
         for (top = blockSize; top; )
            if (!converged(--top)) break;
         top++;

         if (nConvBottom > 0 || top < blockSize)  {
            int nRest = top - nConvBottom;
            SX_CHECK (nRest > 0, nRest);

            // --- reformat variables
            PsiG newBlock;
            int ng = gk.ng, nBasis = (int)X.nRows ();

            newBlock = SxDiracMat<TPrecCoeffG> (ng, nRest);
            newBlock.handle->auxData = dPsi.handle->auxData;
            newBlock <<= dPsi(SxIdx (nConvBottom * ng, ng * top - 1));
            dPsi = PsiG ();
            dPsi = newBlock;

            newBlock = SxDiracMat<TPrecCoeffG> (nBasis, nRest);
            newBlock.handle->auxData = X.handle->auxData;
            newBlock <<= X(SxIdx (nConvBottom * nBasis, nBasis * top - 1));
            X = PsiG ();
            X = newBlock;

            // transfer converged data
            for (int iRest = 0; iRest < nRest; ++iRest)  {
               converged(iRest) = converged(iRest + nConvBottom);
            }
            converged.resize (nRest, true);

            // --- transfer residue-norm data
            resNorm0 = resNorm0(SxIdx(nConvBottom, top - 1)).getCopy ();
            resNorm = resNorm(SxIdx(nConvBottom, top - 1)).getCopy ();

            // --- pick new state block
            psiI = PsiG ();
            eps = SxDiracVec<Double> ();
            i          += nConvBottom;
            nConvBlock -= blockSize - nRest;
            blockSize  = nRest;
            psiI = waves.getBlock(i, blockSize, iSpin,ik);
            eps  = fermi.eps(iSpin,ik)(SxIdx(i,i+blockSize-1));

            if (verbose)  {
               sxprintf ("states: (%d..%d,%d,%d)\n", i+1, i+blockSize, iSpin+1,
                       ik+1);
               fflush(stdout);
            }
         }
      }
   }
}

void SxHamSolver::blockRmmCG (const SxSymbolTable *cmd, bool sloppyConv, 
                              bool printOnly)
{
   SX_NO_MPI;
   SX_CLOCK (Timer::ElMinim);

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SxOverlap S = hamPtr->getS ();
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());

   int maxItDiag = 4;
   bool verbose = false;
   int defBlockSize = 64;
   dRelR = 0.03;
   try  {
      verbose = cmd->contains("verbose")
              ? cmd->get("verbose")->toAttribute ()
              : false;
      if (cmd->contains("maxStepsCCG"))
         maxItDiag = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("blockSize"))
         defBlockSize = cmd->get("blockSize")->toInt ();
      if (cmd->contains("dRelR"))
         dRelR = cmd->get("dRelR")->toReal ();
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   if (printOnly)  {
      cout << "|   block RMM-CG:\n";
      cout << "|      block size:             " << defBlockSize << endl;
      cout << "|      max. steps per state:   " << maxItDiag << endl;
      cout << "|      rel. change in residue: " << (dRelR) * 100. << "%\n";
      return;
   }

   nConvAbs = nCosAcc = nConvIt = nConvOcc = nConvRel = 0;
   nConvCycle = 0;
   nStepsAbs = nStepsCos = nStepsIt = nStepsOcc = nStepsRel = 0;

   if (verbose)
      cout << "abs dEps for this step: " << dEps << endl;
   for (int ik=0; ik < waves.getNk (); ik++)  {
      for (int iSpin=0; iSpin < nSpin; iSpin++)  {

         // --- subspace diagonalization (ref 3, IIIC)
         fermi.eps(iSpin, ik) = subspaceDiagonalization (&waves(iSpin,ik));

         // --- RMM for block of bands
         int blockSize;
         for (int i=0; i < nStates; i+= blockSize)  {
            blockSize = defBlockSize;
            if (i + blockSize >= nStates)  {
               blockSize = nStates-i;
            } else  {
               // try not to cut through degenerate states. Physically not
               // necessary, but improves stability wrt numerical noise
               while (   blockSize > 1 
                      && blockSize > defBlockSize - 3
                      && fabs( fermi.eps(i+blockSize-1, iSpin, ik) 
                              -fermi.eps(i+blockSize  , iSpin, ik)) < 1e-7)
                  blockSize--; 
            }
            blockRmmCG (i, iSpin, ik, blockSize, maxItDiag, verbose, sloppyConv);
         }

         // orthonormalization
         S.orthonormalize (&waves(iSpin, ik));

      }
   }

   printConvergenceStatistics ();
   cout << SX_SEPARATOR;
}

void SxHamSolver::blockRmmCG (int i, int iSpin, int ik, int blockSize, 
                                int maxItDiag, bool verbose, 
                                bool /* sloppyConv */)
{
   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (ik >= 0 && ik < waves.getNk (), ik, waves.getNk ()); 
   SX_CHECK (iSpin >= 0 && iSpin < nSpin, iSpin, nSpin);
   if (verbose)
      cout << "blockRmmCG for " << (i+1) << "->" << i+blockSize << endl;

   const SxGBasis &gk = waves.getGkBasis ()(ik);
   const SxBasis &waveBasis = waves.getBasis (ik);
   SxOverlap S = hamPtr->getS ();
   SxDiracVec<TPrecCoeffG> psiI, dPsi, hPsi, sPsi, X;
   double olddPsi2 = -1.;
   //SxArray<bool> converged(blockSize);
   SxVector<Double> resNorm(blockSize), resNorm0(blockSize);
   SxDiracVec<Double> eps(blockSize), epsPrev(blockSize);
   //converged.set (false);

   psiI = waves.getBlock(i, blockSize, iSpin,ik);
   eps  = fermi.eps(iSpin,ik)(SxIdx(i,i+blockSize-1));

   hPsi = (*hamPtr) | psiI;
   sPsi = S * psiI;

   double R0first = -1., R0last = -1.;

   for (int itDiag = 0; itDiag < maxItDiag; ++itDiag)  {
      
      if (verbose) cout << "it = " << (itDiag + 1) << endl;

      SxDiracMat<TPrecCoeffG> Lambda;
      Lambda = (gk | psiI).overlap (hPsi);
      PsiG Residuum = hPsi - (sPsi ^ Lambda);
      PsiG K = hamPtr->preconditioner (psiI);
      //K.set (1.);
      PsiG KR, KS, KKR;
      KKR.reformat (gk.ng, blockSize);
      KKR.handle->auxData = hPsi.handle->auxData;
      KR.copy (Residuum);
      KS.copy (sPsi);
      for (int j = 0; j < blockSize; ++j)  {
         KR.colRef (j) *= K;
         KS.colRef (j) *= K;
         KKR.colRef (j) <<= K * KR.colRef(j);
      }
      {
         KKR = waveBasis | (KKR);
         dPsi = (-2. * hPsi ^ (KS.overlap (KR)).real ())
              + (*hamPtr * KKR)
              - ((S * KKR) ^ Lambda);
      }
      // --- orthogonalize
      dPsi -= psiI ^ sPsi.overlap (dPsi);

      if (itDiag == 0)  {
         X = -dPsi;
         olddPsi2 = dPsi.normSqr ();
         R0first = KR.normSqr ();
      } else {
         double dPsi2 = dPsi.normSqr ();
         double gamma = dPsi2 / olddPsi2;
         olddPsi2 = dPsi2;
         X = -dPsi + gamma * X;
      }

      // --- orthogonalize
      X -= psiI ^ sPsi.overlap (X);

      // --- attach full basis
      X = waveBasis | X;

      // --- line search
      double R0 = KR.normSqr ();
      double dR = 2. * dot (dPsi, X).re;

      /*
      if (verbose)  {
         for (double lTrial = -0.5; lTrial < 0.5; lTrial+=0.05)  {
            PsiG psiTrial = psiI + lTrial * X;
            S.orthonormalize (&psiTrial);
            PsiG hPsiTrial = *hamPtr | psiTrial;
            PsiG resTrial = hPsiTrial 
                          - ((S * psiTrial) ^ (gk | psiTrial).overlap (hPsiTrial));
            for (int i = 0; i < blockSize; ++i)
               resTrial.colRef (i) *= K;
            cout << lTrial << " " << resTrial.normSqr () << endl;
         }
      }
      */
      double lTrial = 1;
      PsiG psiTrial = psiI + lTrial * X;
      S.orthonormalize (&psiTrial);
      PsiG hPsiTrial = *hamPtr | psiTrial;
      PsiG resTrial = hPsiTrial 
                    - ((S * psiTrial) ^ (gk | psiTrial).overlap (hPsiTrial));
      for (int j = 0; j < blockSize; ++j)
         resTrial.colRef (j) *= K;

      double curv;
      double lambdaMin = lineMinimization (R0, lTrial, resTrial.normSqr (),
                                           dR, &curv);
      R0last = R0 + lambdaMin * (dR + curv * lambdaMin);
      if (verbose)  {
         cout << "R0    = " << R0 << endl;
         cout << "dR    = " << dR << endl;
         cout << "RT    = " << resTrial.normSqr () << endl;
         cout << "lMin  = " << lambdaMin << endl;
         cout << "curv  = " << curv << endl;
         cout << "RPred = " << R0last << endl;
         cout << "|d|^2 = " << olddPsi2 << endl;
      }
      psiI.plus_assign_ax (lambdaMin, X);
      S.orthonormalize (&psiI);

      hPsi = (*hamPtr) * psiI;
      sPsi = PsiG ();
      sPsi = S * psiI;

   }
   cout << "R0(final)/R0(0) = " << R0last/R0first << endl;

   SxDiracMat<TPrecCoeffG> H = (gk | psiI).overlap (hPsi);
   SxDiracMat<TPrecCoeffG>::Eigensystem eig = H.eigensystem ();
   psiI.rotate (eig.vecs);
   eps <<= eig.vals;
}

void SxHamSolver::printConvStat (const SxString &name, int nConv, int nStep)
{
   SX_CHECK (name.getSize () <= 30, name.getSize ());
   cout << "|    " << name << ':';
   for (int i = int(name.getSize ()); i < 30; ++i) cout << ' ';
   cout << nConv << " times";
   if (nStep > 0 && nConv > 0)
      cout << " (" << (nStep / double(nConv)) << " steps)";
   cout << endl;
}

void SxHamSolver::printConvergenceStatistics () const
{
   // --- print convergence statistics
   cout << SX_SEPARATOR;
   cout << "| Convergence statistics:\n";
   
   printConvStat ("absolute convergence", nConvAbs, nStepsAbs);
   printConvStat ("relative convergence", nConvRel, nStepsRel);
   printConvStat ("low occupation", nConvOcc, nStepsOcc);
   printConvStat ("numerical cos limit reached", nCosAcc, nStepsCos);
   printConvStat ("exceeded max. iterations", nConvIt, nStepsIt);

   if (nConvCycle >= 0)
      cout << "|    number of unnecessary steps:   " << nConvCycle << endl;
   int stepSum = nStepsAbs + nStepsRel + nStepsOcc + nStepsCos + nStepsIt;
   if (nConvCycle > 0)
   stepSum += nConvCycle;
   if (stepSum > 0)
      cout << "|    number of steps (total):       " << stepSum << endl;
   // no SX_SEPARATOR, so that further stuff can be appended
   //cout << SX_SEPARATOR;
}

void SxHamSolver::printEnergyDatLine(int it, double eTot, double freeEnergy,
                                     double eBand, double entropy)
{
   SX_MPI_MASTER_ONLY {
      SxFileIO::appendToFile
      ( SxString(it)                                + "\t" // iteration
      + SxString(GETTIME(Timer::ElMinim), "%12.8f") + "\t" // time
      + SxString(eTot, "%15.12f")                   + "\t" // energy
      + SxString(freeEnergy, "%15.12f")             + "\t" // free energy
      + SxString(0.5 * (eTot + freeEnergy), "%15.12f") + "\t" // T->0 energy
      + SxString(eBand, "%15.12f")                  + "\t" // band energy
      + SxString(entropy, "%15.12f")                + "\n" // entropy
      , "energy.dat");
   }
}


void SxHamSolver::stateByStateCG (const SxSymbolTable *cmd, bool sloppyConv,
                                  bool printOnly)
{
   SX_CLOCK (Timer::ElMinim);

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   int nk = waves.getNk ();
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   nConvIt=nConvRel=nConvAbs=nConvOcc=nCosAcc = 0;
   int nEpsAcc=0;
   SxDiracVec<TPrecCoeffG> psiI, hPsi, dPsi, dPsiOld, K, X, Xold, hX;
   SxDiracVec<TPrecCoeffG> psiBackupG;
   Real8 eps, gamma = 0., theta, epsPrev=0., eps1st=0.;
   PrecCoeffG tr, trOld = 0.;
   SxMatrix<TPrecCoeffG> M(2,2), U(2,2);
   int maxItDiag = 10;
   bool verbose;

   try  {
      ekt     = cmd->contains("ekt")
              ? cmd->get("ekt")->toReal() / HA2EV 
              : ektDefault;
      dRelEps = cmd->contains("dRelEps")
              ?  cmd->get("dRelEps")->toReal()
              :  0.30;   // ref1, conv. conditions CG
      verbose = cmd->contains("verbose")
              ? cmd->get("verbose")->toAttribute ()
              : false;
      if (cmd->contains("maxStepsCCG"))
         maxItDiag = cmd->get("maxStepsCCG")->toInt ();
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   if (printOnly)  {
      cout << "|   conjugate gradient:\n";
      cout << "|      max. steps per state:   " << maxItDiag << endl;
      cout << "|      rel. change of eps:     " << 100. * dRelEps  << " %\n";
      return;
   }

   SxOverlap S = hamPtr->getS ();

   SX_NEW_LOOP (waves);
   for (int ik=0; ik < nk; ik++)
   {
      SX_MPI_LEVEL("waves-k");
      SX_ALLOC_CACHE;
      for (int iSpin=0; iSpin < nSpin; iSpin++)
      {
         if (SxLoopMPI::myWork(ik))    // LoopMPI
         {
            const SxGBasis &gk = waves.getGkBasis ()(ik);
            const SxBasis &waveBasis = waves.getBasis (ik);
            for (int i=0; i < nStates; i++)  {
               psiI = waves(i,iSpin,ik);
               if (verbose)  {
                  sxprintf ("state: (%d,%d,%d) ", i, iSpin, ik); fflush(stdout);
               }
               epsPrev = 1e50;  // just a large value
               for (int itDiag=0; /* empty */; itDiag++)  {

                  // orthogonalize psiI to lower states
                  if (i > 0)
                     S.setOrthogonal (&psiI, waves.getBlock (0, i, iSpin, ik));
                  S.normalize (&psiI);

                  // --- compute gradient: |dPsi> = (H-e) |psi>
                  hPsi = *hamPtr | psiI;
                  eps  = dot (gk | psiI, hPsi).re;
                  dPsi = hPsi - eps * (S * psiI);
                  fermi.eps(i,iSpin,ik) = eps;
                  if (itDiag == 1)  eps1st = fabs(eps - epsPrev);

                  if (itDiag > 0 && fabs(eps - epsPrev) < dEps)  {
                     if (verbose)
                        sxprintf ("conv due to abs eps (%d of %d)\n",
                              itDiag+1, maxItDiag);
                     nConvAbs++;
                     break;
                  }

                  if (eps > epsPrev)  {
                     if (verbose)
                        sxprintf ("numerical limit reached (%d of %d)\n",
                              itDiag+1, maxItDiag);
                     nEpsAcc++;
                     waves(i,iSpin,ik) <<= psiBackupG;
                     break;
                  }

                  // --- orthogonalize preconditioned gradient
                  K = waveBasis | (hamPtr->preconditioner (gk | psiI) * dPsi);
                  S.setOrthogonal (&K, waves.getBlock (0, i+1, iSpin, ik));

                  // --- get conjugate direction, ref1 (8.14b) / ref2 (10.6.7)
                  tr = (dPsi ^ (gk | K)).chop();
                  if (itDiag > 0)
                     gamma = (tr - (dPsiOld^(gk | K)).chop()) / trOld;
                  //             if (itDiag > 0)  gamma = tr / trOld;
                  dPsiOld.copy (dPsi);
                  trOld = tr;

                  // --- compute search direction, ref1 (8.14a)
                  if (itDiag == 0)  X = K * 1.;
                  else              X = K + gamma * Xold;
                  Xold.resize (0); Xold.copy (X);
                  {
                     SX_CLOCK (Timer::ortho);
                     PrecCoeffG scp = (psiI | S | X);
                     // X -= psiI * scp;
                     X.plus_assign_ax (-scp, psiI);
                  }
                  S.normalize (&X);
                  // TODO: ugly
                  X.handle->auxData.i     = i;
                  X.handle->auxData.iSpin = iSpin;
                  X.handle->auxData.ik    = ik;

                  // --- get M matrix (in ref1: H matrix)
                  hX = *hamPtr | X;
                  M(0,0) = eps;            M(0,1) = ((gk | X) ^ hPsi).chop();
                  M(1,0) = M(0,1).conj();  M(1,1) = ((gk | X) ^ hX).chop();

                  // --- uniform transformation U, ref1 (8.16)
                  if ( M(0,0).re == M(1,1).re )  // '==' is correct!
                  {
                     /* This is a workaround for the numerical limit of
                     degenerate eigenstates. If an exact eigenstate has
                     been found, and X is a different eigenstate to the
                     same energy, the formula below doesn't work, as
                     atan(1/0) is undefined. Any theta is fine, so
                     we just take 0 to get out at the numerical limit
                     test. */
                     theta = 0.;
                  }  else  {
                     theta = 0.5 * atan (  (M(1,0).re + M(0,1).re)
                           / (M(0,0).re - M(1,1).re) );
                     // --- in case of being far off the quadr. regime theta
                     //     points to the maximum rather than to the minimum!
                     if (M(0,0).re > M(1,1).re)  theta -= PI/2.;
                  }

                  SX_CHECK_NUM (theta);
                  if (fabs(theta) < 1e-8)  {
                     if (verbose)
                        sxprintf ("Numerical limit reached (%d of %d)\n",
                              itDiag+1, maxItDiag);
                     nCosAcc++;
                     break;
                  }
                  U(0,0) = cos (theta); U(0,1) = sin (theta);
                  U(1,0) = -U(0,1);     U(1,1) = U(0,0);

                  // --- compute new wavefunction, ref1 (8.15b)
                  psiBackupG = U(0,0) * psiI + U(0,1) * X;
                  waves(i,iSpin,ik) <<= psiBackupG;

                  // --- convergence?
                  if (sloppyConv && itDiag >= 2
                      && fabs(fermi.focc(i,iSpin,ik)) < 1e-8)
                  {
                     if (verbose)
                        sxprintf ("conv. due to low occ. (%d of %d)\n",
                              itDiag+1, maxItDiag);
                     nConvOcc++;
                     break;
                  }
                  if (itDiag > 0 && sloppyConv
                        && fabs(eps - epsPrev) < dRelEps * eps1st)
                  {
                     if (verbose)
                        sxprintf ("conv due to rel eps (%d of %d)\n",
                              itDiag+1, maxItDiag);
                     nConvRel++;
                     break;
                  }
                  if (itDiag+1 >= maxItDiag)  {
                     if (verbose)
                        sxprintf ("WARNING: Could not converge (%d of %d) dEps=%g\n",
                              itDiag+1, maxItDiag, eps-epsPrev);
                     nConvIt++;
                     break;
                  }
                  epsPrev = eps;
               }

            }
            // --- subspace diagonalization (for metalls only)
            fermi.eps(iSpin, ik) = subspaceDiagonalization (&waves(iSpin, ik));
         }   // LoopMPI
      }
   }
   fermi.eps.synMPI ();

   // --- print convergence statistics
   SX_MPI_SOURCE("waves-k", TaskGroupMaster);
   SX_MPI_TARGET(TopLevel, TaskGroupAll);
   cout << SX_SEPARATOR;
   cout << "| Convergence statistics:\n";
   cout << "|    absolute convergence:        " << SxLoopMPI::sum(nConvAbs) << " times\n";  // LoopMPI
   cout << "|    relative convergence:        " << SxLoopMPI::sum(nConvRel) << " times\n";  // LoopMPI
   cout << "|    low occupation:              " << SxLoopMPI::sum(nConvOcc) << " times\n";  // LoopMPI
   cout << "|    numerical cos limit reached: " << SxLoopMPI::sum(nCosAcc)  << " times\n";  // LoopMPI
   cout << "|    numerical eps limit reached: " << SxLoopMPI::sum(nEpsAcc)  << " times\n";  // LoopMPI
   cout << "|    exceeded max. iterations:    " << SxLoopMPI::sum(nConvIt)  << " times\n";  // LoopMPI
   cout << SX_SEPARATOR;

}

int SxHamSolver::blockStateByStateCG (const SxSymbolTable *cmd,
                                      MinimizerMode        runType)
{
   int maxItDiag = 5;
   int defBlockSize = (nStates > 64) ? 64 : (nStates / 2);
   if (defBlockSize < 2) defBlockSize = 2;
   int nSloppy = 0;
   bool verbose, numericalLimit = false;

   try  {
      dRelEps = cmd->contains("dRelEps")
              ?  cmd->get("dRelEps")->toReal()
              :  0.30;   // ref1, conv. conditions CG
      verbose = cmd->contains("verbose")
              ? cmd->get("verbose")->toAttribute ()
              : false;
      if (cmd->contains("blockSize"))
         defBlockSize = cmd->get("blockSize")->toInt ();
      if (cmd->contains("maxStepsCCG"))
         maxItDiag = cmd->get("maxStepsCCG")->toInt ();
      if (cmd->contains("numericalLimit"))
         numericalLimit = cmd->get("numericalLimit")->toAttribute ();
      if (cmd->contains("nSloppy"))
         nSloppy = cmd->get("nSloppy")->toInt ();
         
   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   return blockStateByStateCG (maxItDiag, defBlockSize, numericalLimit,
                               runType, nSloppy, verbose);
}

int SxHamSolver::blockStateByStateCG (int maxItDiag,
                                      int defBlockSize,
                                      bool numericalLimit,
                                      MinimizerMode runType,
                                      int nSloppy,
                                      bool verbose)
{
   if (defBlockSize > nStates) defBlockSize = nStates;
   if (runType == PrintOnly)  {
      cout << "|   blocked CG minimizer" << endl;
      if (nStates < 10)  {
         cout << "WARNING: blockCCG not recommended for very small systems."
              << endl;
      }
      cout << "|      block size:             " << defBlockSize << endl;
      cout << "|      max. CCG steps:         " << maxItDiag << endl;
      cout << "|      rel. change of eps:     " << 100. * dRelEps  << "%\n";
      cout << "|      numerical noise limit:  try "
           << (numericalLimit ? "hard to get it" : "to avoid it")
           << endl;
      if (nSloppy > 0)
         cout << "|      sloppy convergence:     " << nSloppy << " states\n";
      return -1;
   }

   SX_CLOCK (Timer::ElMinim);

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   int nk = waves.getNk ();
   int blockSize;
   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   SX_CHECK (nStates == waves.getNStates (), nStates, waves.getNStates ());
   
   nConvIt = nConvRel = nConvAbs = nConvOcc = nCosAcc = 0;
   int nEpsAcc=0, nAsSloppy=0;
   nStepsRel=nStepsAbs=nStepsOcc=nStepsCos=0;
   int nStepsEps=0;
   SxDiracVec<TPrecCoeffG> psiI, hPsi, X, Xold;
   SxDiracVec<TPrecEps> deltaEps, eps1st(defBlockSize);
   SxVector<Int> iter(defBlockSize);
   SxList<int> zeroX;
   Real8 gamma;
   double tr, trOld = 0.;
   int iDiagLoop;
   SxDiracSymMat<TPrecCoeffG> subspaceMatrix(nStates);
   SxOverlap S = hamPtr->getS ();

#  ifdef INDIVIDUAL_K_TIMING
      cout << SX_SEPARATOR;
      cout << "BEGIN_K_TIMING" << endl;
#  endif


   SX_START_TIMER (Timer::parKLoop);

   for (int ik=0; ik < nk; ik++)
   {

      iDiagLoop = 0;
      const SxGBasis &gk = waves.getGkBasis ()(ik);
      const SxBasis &waveBasis = waves.getBasis (ik);

      for (int iSpin=0; iSpin < nSpin; iSpin++)
      {

         SX_MPI_LEVEL("waves-k");
         SX_MPI_SOURCE(CurrentLevel, TaskGroupMaster);
         SX_MPI_TARGET(CurrentLevel, TaskGroupAll);
         if (SxLoopMPI::myWork(ik)) // LoopMPI
         {
            SX_ALLOC_CACHE;

#           ifdef INDIVIDUAL_K_TIMING
               double t0 = SxTime::getRealTime ();
#           endif

            iter.set (0);
            for (int i=0; i < nStates; /* done inside loop */)  {
               blockSize = defBlockSize;
               if (i + blockSize >= nStates)  {
                  blockSize = nStates-i;
               } else  {
                  // try not to cut through degenerate states. Physically not
                  // necessary, but improves stability wrt numerical noise
                  while (   blockSize > 1
                        && blockSize > defBlockSize - 3
                        && fabs( fermi.eps(i+blockSize-1, iSpin, ik)
                              -fermi.eps(i+blockSize  , iSpin, ik)) < 1e-7)
                     blockSize--;
#ifdef SX_LOOP_MPI
                  // synchronize block size
                  blockSize = SxLoopMPI::bcast (blockSize, 0);
#endif
               }
               psiI = PsiG ();
               psiI = waves.getBlock(i, blockSize, iSpin,ik);
               if (verbose)  {
                  sxprintf ("states: (%d..%d,%d,%d) ", i+1, i+blockSize, iSpin+1,
                        ik+1);
                  fflush(stdout);
               }
               for (int itDiag=0; /* empty */; itDiag++, iter+=1)  {

                  if (itDiag % 5 == 0)  {
                     // orthogonalize psiI to lower states
                     if (i > 0)
                        S.setOrthogonal (&psiI,waves.getBlock(0,i,iSpin,ik));
                     S.orthonormalize (&psiI, GramSchmidt);
#ifdef SX_LOOP_MPI
                     // synchronize psi
                     SxLoopMPI::bcast (psiI, 0);
#endif

                     // --- compute gradient
                     psiI.handle->auxData.iSpin = iSpin;
                     psiI.handle->auxData.ik = ik;
                     hPsi = *hamPtr | psiI;
                     if (itDiag == 0)  {
                        //eps1st = SxDiracVec<TPrecEps>(blockSize);
                        double eps;
                        for (int j = 0; j < blockSize; ++j)  {
                           eps = dot(gk | psiI.colRef(j), hPsi.colRef(j));
                           //eps1st(j) = fabs(eps - fermi.eps(i+j,iSpin,ik));
                           fermi.eps(i + j, iSpin, ik) = eps;
                        }
                     }
                  }

                  // --- get search direction
                  PsiGI Spsi = S * psiI;
                  {
                     SX_CLOCK(Timer::SearchDir);
                     SxDiracVec<TPrecCoeffG> dPsi;
                     SxDiracMat<TPrecCoeffG> K(gk.ng, blockSize);
                     // TODO: ugly
                     K.handle->auxData = psiI.handle->auxData;
                     K.setBasis (gk);
                     tr = 0.;
                     for (int j = 0; j < blockSize; ++j)  {
                        double eps = dot(gk | psiI.colRef(j),hPsi.colRef(j));
                        //dPsi = hPsi.colRef(j)
                        //              - eps * (Spsi.colRef(j));
                        //K.colRef(j) <<= hamPtr->preconditioner (gk | psiI.colRef(j))
                        //                * dPsi;
                        //tr += dPsi.normSqr ();
                        double tr_ = 0.;
                        SxDiracVec<Double> precond = hamPtr->preconditioner (gk | psiI.colRef(j));
                        ssize_t ng = gk.ng;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:tr_)
#endif
                        for (ssize_t ig = 0; ig < ng; ++ig)  {
                           SxComplex16 dPsi_ = hPsi(ig,j) - eps * Spsi(ig, j);
                           K(ig, j) = precond(ig) * dPsi_;
                           tr_ += dPsi_.absSqr ();
                        }
                        tr += tr_;
                     }
                     // change block to wave basis
                     K = waveBasis | K;
                     // --- orthogonalize preconditioned gradient
                     S.setOrthogonal (&K,
                           waves.getBlock (0, i+blockSize, iSpin, ik));

                     if (itDiag > 0)  {
                        // get conjugate direction, ref1 (8.14b) / ref2 (10.6.7)
                        gamma = tr / trOld;
                        // compute search direction, ref1 (8.14a)
                        X = K + gamma * Xold;
                     } else {
                        X = K;
                     }
                     trOld = tr;

                     Xold = PsiG (); Xold.copy (X);
                  }

                  {
                     // --- orthogonalize to psiI
                     SX_CLOCK (Timer::ortho);
                     //S.setOrthogonal (&X, psiI);
                     //X -= psiI ^ (Spsi.adjoint () ^ (gk | X));
                     X -= psiI ^ Spsi.overlap (X, gk.ng);

                     // --- orthonormalize (track zero search directions)
                     zeroX.resize (0);
                     SxDiracMat<TPrecCoeffG>::Eigensystem eig;
                     eig = S.getMatrix (X, X).eigensystem ();
                     if (eig.vals(0).re > 1e-12 * fabs(eig.vals(blockSize-1).re))  {
                        // Loewdin orthogonalisation
                        SxDiracMat<TPrecCoeffG> Ieps;
                        Ieps = Ieps.identity (1. / sqrt(eig.vals.real ()));
                        X.rotate (eig.vecs ^ Ieps ^ eig.vecs.adjoint ());
                     } else {
                        cout << "near-linear dependencies in search space"
                              << endl;
                        // --- orthonormalize (track zero search directions)
                        for (int j = 0; j < blockSize; ++j)  {
                           SxDiracVec<Complex16> Xj = X.colRef(j);
                           // --- Gram-Schmidt orthogonalization
                           if (j > 0)  {
                              int ng = (int)X.nRows ();
                              SxDiracVec<Complex16> Xlow = X(SxIdx(0, j * ng - 1));
                              Xlow.reshape (ng, j);
                              S.setOrthogonal (&Xj, Xlow);
                           }
                           double norm2 = (Xj | S | Xj).re;
                           if (norm2 > 1e-26)  {
                              // Xj /= sqrt(norm2);
                              // workaround for 1e-10 check in SxVec::operator/=
                              Xj *= 1./sqrt(norm2);
                           } else  {
                              Xj.set (0.);
                              zeroX.append(j);
                              if (verbose) cout << "X(" << (i+j) << ")=0" << endl;
                           }
                        }
                     }
                  }
                  Spsi = SxDiracVec<TPrecCoeffG> (); // free memory
                  // TODO: ugly
                  X.handle->auxData.i     = -1;
                  X.handle->auxData.iSpin = iSpin;
                  X.handle->auxData.ik    = ik;
                  X.setBasis (waveBasis);

                  // --- compute new wavefunctions
                  SxArray<bool> numCosLim(2*blockSize);
                  {
                     PsiGI hX = *hamPtr | X;
                     SxDiracSymMat<TPrecCoeffG> ham(2 * blockSize);
                     if (blockSize > 1) {
                        SX_CLOCK(Timer::SubspaceMatrix);
                        PsiG psiIadj = (gk | psiI).adjoint ();
                        PsiG psiHpsi = psiIadj ^ hPsi;
                        PsiG psiHX   = psiIadj ^ hX;
                        psiIadj = PsiG (); // free memory
                        PsiG XhX     = (gk | X).adjoint () ^ hX;

                        if (zeroX.getSize () > 0)  {
                           // set diagonal of zero X to huge number
                           SxList<int>::Iterator it;
                           double bigH = psiHpsi.norm ();
                           for (it = zeroX.begin (); it != zeroX.end (); ++it)  {
                              XhX(*it,*it) = bigH;
                           }
                        }

                        int iL, iR, iLx, iRx;
                        for (iL = 0; iL < blockSize; ++iL)  {
                           iLx = iL + blockSize;
                           for (iR = 0; iR < blockSize; ++iR)  {
                              iRx = iR + blockSize;
                              if (iL <= iR)  {
                                 ham(iL, iR)   = 0.5 * (  psiHpsi(iL,iR)
                                       + psiHpsi(iR,iL).conj ());
                                 ham(iLx, iRx) = 0.5 * (  XhX(iL,iR)
                                       + XhX(iR,iL).conj ());
                              }
                              ham(iL,iRx) = psiHX(iL,iR);
                           }
                        }
                     } else {
                        SX_CLOCK(Timer::SubspaceMatrix);
                        ham(0,0) = dot(gk | psiI, hPsi);
                        ham(0,1) = dot(gk | psiI, hX);
                        ham(1,1) = dot(gk | X   , hX);
                     }
                     for (int iL = 0; iL < blockSize; ++iL)  {
                        int iR = iL+blockSize;
                        double deltaE = ham(iL,iL).re - ham(iR,iR).re;
                        numCosLim(iL) = ham(iL,iR).re < 0.5e-9 * fabs(deltaE)
                                               || (deltaE == 0.);
                     }

                     // --- minimize band energy by computing eigensystem
                     SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;
                     {
                        SX_CLOCK (Timer::SubspaceDiag);
                        eig = ham.eigensystem ();
                     }

                     // --- update eigenenergies and their change
                     deltaEps = eig.vals(SxIdx(0,blockSize-1))
                                       - fermi.eps(iSpin, ik)(SxIdx(i,i+blockSize-1));
                     for (int j = 0; j < blockSize; ++j)
                        if (iter(j) == 0) eps1st(j) = fabs(deltaEps(j));
                     if (verbose)  {
                        /*
                     cout << endl << "eigenvals=" << eig.vals << endl
                          << "fermi.eps="
                          << fermi.eps(iSpin, ik)(SxIdx(i,i+blockSize-1));
                         */
                        cout << endl << "deltaEps = " << deltaEps << endl;
                     }
                     fermi.eps(iSpin, ik)(SxIdx(i,i+blockSize-1))
                                 <<= eig.vals(SxIdx(0,blockSize-1));

                     // --- update waves and their derivatives
                     if (blockSize > 1)  {
                        SX_CLOCK(Timer::WaveUpdate);
                        PsiG eigVec;
                        SxDiracMat<TPrecCoeffG> eigVecI(blockSize,blockSize),
                              eigVecX(blockSize,blockSize);
                        for (int j = 0; j < blockSize; ++j)  {
                           eigVec = eig.vecs.colRef(j);
                           // I-dependent part of lower eigenvectors
                           eigVecI.colRef(j) <<= eigVec(SxIdx(0,blockSize-1));
                           // X-dependent part of lower eigenvectors
                           eigVecX.colRef(j) <<= eigVec(SxIdx(blockSize,
                                 2*blockSize-1));
                           // higher eigenvectors are not used
                        }
                        /*
                     // more elegant, but needs more memory
                     psiI <<= (psiI ^ eigVecI) + (X ^ eigVecX);
                     hPsi   = (hPsi ^ eigVecI) + (hX ^ eigVecX);
                         */
                        // same as above, but needs less memory
                        psiI.rotate (eigVecI);
                        psiI += X ^ eigVecX;
                        hPsi.rotate (eigVecI);
                        hPsi += hX ^ eigVecX;

                        // rotate search direction for next step accordingly
                        Xold = Xold ^ eigVecI;
                     } else {
                        SX_CLOCK(Timer::WaveUpdate);
                        // vector class confused for single psi
                        psiI <<= eig.vecs(0,0) * psiI + eig.vecs(1,0) * X;
                        hPsi   = eig.vecs(0,0) * hPsi + eig.vecs(1,0) * hX;
                        Xold *= eig.vecs(0,0);
                     }
                  }
#ifdef SX_LOOP_MPI
                  // synchronize psi
                  SxLoopMPI::bcast (psiI, 0);
#endif

                  // --- check for convergence
                  int j;
                  for (j = 0; j < blockSize; ++j)  {
                     if (numCosLim(j) && !numericalLimit)  {
                        if (verbose)
                           sxprintf ("numerical cos limit reached (%d of %d)\n",
                                 iter(j)+1, maxItDiag);
                        nCosAcc++;
                        nStepsCos += iter(j);
                     } else if (   (fabs(deltaEps(j)) < dEps
                           && (iter(j) > 0 || runType != SCFrun))
                           || (fabs(deltaEps(j)) < 1e-2 * dEps
                                 && !numericalLimit)
                                 || (fabs(deltaEps(j)) < 1e-4
                                       && runType == PropagateRun))  {
                        if (verbose)
                           sxprintf ("conv due to abs eps (%d of %d)\n",
                                 iter(j)+1, maxItDiag);
                        nConvAbs++;
                        nStepsAbs += iter(j);
                     } else if (iter(j) > 0 && runType == SCFrun
                           && fabs(deltaEps(j)) < dRelEps * eps1st(j)) {
                        if (verbose)
                           sxprintf ("conv due to rel eps (%d of %d)\n",
                                 iter(j)+1, maxItDiag);
                        nConvRel++;
                        nStepsRel += iter(j);
                     } else if (runType == SCFrun && iter(j) >= 2
                           && fermi.focc(i+j,iSpin,ik) < 1e-4)  {
                        if (verbose)
                           sxprintf ("conv. due to low occ. (%d of %d)\n",
                                 iter(j)+1, maxItDiag);
                        nConvOcc++;
                        nStepsOcc += iter(j);
                     } else if (iter(j)+1 >= maxItDiag)  {
                        if (verbose)
                           sxprintf ("WARNING: Could not converge (%d of %d) dEps=%g\n",
                                 iter(j)+1, maxItDiag, deltaEps(j));
                        nConvIt++;
                        if (runType != SCFrun && i+j+nSloppy >= nStates)
                           nAsSloppy++;
                     } else {
                        break;
                     }
                  }
#ifdef SX_LOOP_MPI
                  // synchronize number of converged states
                  j = SxLoopMPI::bcast (j, 0);
#endif
                  if (j > 0)  {
                     // --- setup (parts of) subspace matrix
                     int ng = (int)hPsi.nRows ();
                     SxDiracMat<TPrecCoeffG> psiLow, hPsiConv;
                     psiLow = (gk | waves.getBlock(0, i+j, iSpin, ik));
                     hPsiConv = hPsi(SxIdx(0, ng*j - 1));
                     hPsiConv.reshape(ng, j);
                     PsiG psiConvHPsiLow = hPsiConv.adjoint () ^ psiLow;
                     int iL, iR;
                     for (iR = 0; iR < i+j; ++iR)
                        for (iL = 0; iL < j; ++iL)
                           if (iR >= i + iL)
                              subspaceMatrix(i+iL,iR) = psiConvHPsiLow(iL,iR);
                           else
                              subspaceMatrix(iR,i+iL) = psiConvHPsiLow(iL,iR)
                              .conj ();

                     // --- update convergence meta data
                     i+=j; // j many states are converged
                     int k;
                     for (k = 0; k < blockSize - j; ++k) {
                        iter(k) = iter(k+j) + 1; // inc iter for unconverged states
                        eps1st(j) = deltaEps(k+j);
                     }
                     // reset iterations for next states
                     for (; k < defBlockSize; ++k) iter(k) = 0;

                     break;  // new block
                  }
               }
            }
            {
               // --- diagonalize subspace matrix and update waves
               SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;
               {
                  SX_CLOCK(Timer::SubspaceDiag);
                  eig = subspaceMatrix.eigensystem ();
               }
#ifdef SX_LOOP_MPI
                // synchronize eigenvalues
                SxLoopMPI::bcast (eig.vals, 0);
#endif
               fermi.eps(iSpin,ik) = eig.vals;

               SX_CLOCK(Timer::WaveUpdate);
#ifdef SX_LOOP_MPI
                // synchronize subspace rotation
                SxLoopMPI::bcast (eig.vecs, 0);
#endif
               // waves(iSpin,ik) <<= waves(iSpin,ik) ^ eig.vecs;
               waves(iSpin, ik).rotate (eig.vecs);

            }
            if (runType == PropagateRun)  {
               if (iDiagLoop >= 5 || nConvIt == nAsSloppy)  {
                  iDiagLoop = 0;
                  nConvIt = nAsSloppy = 0;
               } else {
                  iDiagLoop++;
                  // do this (iSpin, ik) again
                  cout << "Propagation loop " << iDiagLoop
                        << " for ik=" << (ik+1) << " iSpin=" << (iSpin+1) << endl;
                  iSpin--;
               }
            }

#           ifdef INDIVIDUAL_K_TIMING
               double t1 = SxTime::getRealTime ();
               // cout << "k=" << ik << " took " << (t1-t0) << "s" << endl;
               cout << ik << "   " << (t1-t0) << endl;
#           endif

         } // LoopMPI
      } // iSpin
      // khr:  the following block is not part of the test example
      if (runType == PropagateRun && (ik+1 < waves.getNk ()))  {
         propagateWaves (ik, ik+1);
         // --- the next thing is quite memory consuming:
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)
            fermi.eps(iSpin, ik) = subspaceDiagonalization (&waves(iSpin,ik+1));
      }
   } // ik


   SX_STOP_TIMER (Timer::parKLoop);
   fermi.eps.synMPI ();



#  ifdef INDIVIDUAL_K_TIMING
      cout << "END_K_TIMING" << endl;
      cout << SX_SEPARATOR;
#  endif


   if (runType == PropagateRun) return 1; // no convergence statistics


   // do the sums using a single MPI call for performance reasons
   {
      SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
      SX_MPI_TARGET (TopLevel, TaskGroupAll);
      // SX_MPI_SUM_VARS (nConvAbs, nStepsAbs, nConvRel, nStepsRel,
      //                  nConvOcc, nStepsOcc, nCosAcc, nStepsCos,
      //                  nEpsAcc, nStepsEps, nConvIt, nAsSloppy);
      SxVector<Int> sumTmp(12);
      sumTmp(0)  = nConvAbs;   sumTmp(1)  = nStepsAbs;  sumTmp(2)  = nConvRel;
      sumTmp(3)  = nStepsRel;  sumTmp(4)  = nConvOcc;   sumTmp(5)  = nStepsOcc;
      sumTmp(6)  = nCosAcc;    sumTmp(7)  = nStepsCos;  sumTmp(8)  = nEpsAcc;
      sumTmp(9)  = nStepsEps;  sumTmp(10) = nConvIt;    sumTmp(11) = nAsSloppy;
      //
      SxLoopMPI::sum (sumTmp);
      //
      nConvAbs   = sumTmp(0);  nStepsAbs  = sumTmp(1);  nConvRel   = sumTmp(2);
      nStepsRel  = sumTmp(3);  nConvOcc   = sumTmp(4);  nStepsOcc  = sumTmp(5);
      nCosAcc    = sumTmp(6);  nStepsCos  = sumTmp(7);  nEpsAcc    = sumTmp(8);
      nStepsEps  = sumTmp(9);  nConvIt    = sumTmp(10); nAsSloppy  = sumTmp(11);
   }


   // --- print convergence statistics
   cout << SX_SEPARATOR;
   cout << "| Convergence statistics:\n";
   printConvStat ("absolute convergence", nConvAbs, nStepsAbs);
   printConvStat ("relative convergence", nConvRel, nStepsRel);
   printConvStat ("low occupation", nConvOcc, nStepsOcc);
   printConvStat ("numerical cos limit reached", nCosAcc, nStepsCos);
   printConvStat ("numerical eps limit reached", nEpsAcc, nStepsEps);
   cout << "|    exceeded max. iterations:      " << nConvIt << " times\n";
   if (runType == BandStructureRun && nSloppy > 0)
   cout << "|       thereof rated as sloppy     " << nAsSloppy << " times\n";
   cout << SX_SEPARATOR;
   return nConvIt - nAsSloppy;
}


void SxHamSolver::directDiagonalization (const SxSymbolTable *cmd, bool calc)
{
   SX_NO_MPI;
   // TODO: use hamPtr here...
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ())); 
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr); 
   SX_CLOCK (Timer::ElMinim);
   SX_CHECK (cmd);

   cout << SX_SEPARATOR;
   cout << "| Direct Diagonalisation" << endl;

   // --- get input parameters
   int                   printSteps;
   bool                  exx;
   RhoR                  newRhoR;

   try  {
      maxSteps     = (cmd->contains("maxSteps"))
                   ?  cmd->get("maxSteps")->toInt()
                   :  1;
      printSteps   = (cmd->contains("printSteps"))
                   ?  cmd->get("printSteps")->toInt()
                   :  1;
      dEnergy      = (cmd->contains("dEnergy"))
                   ?  cmd->get("dEnergy")->toReal()
                   :  -1.;
      exx          = (cmd->contains("exx"))
                   ?  cmd->get("exx")->toAttribute()
                   :  false;
      keepOcc      = (cmd->contains("keepOccFixed"))
                   ?  cmd->get("keepOccFixed")->toAttribute()
                   :  false;
      keepRho      = (cmd->contains("keepRhoFixed"))
                   ?  cmd->get("keepRhoFixed")->toAttribute()
                   :  false;

   }  catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   // --- dump input parameters
   if (exx)  {
      sxprintf ("|   diag. within EXX calculation\n");
   }  else  {
      sxprintf ("|   max. Steps:           %d\n", maxSteps);
      sxprintf ("|   details printed every %d iteration(s)\n", printSteps);
   }
   if (keepOcc)  sxprintf ("|   occupation:           fixed\n");
   if (keepRho)  sxprintf ("|   Bandstructure calculation\n");
   else          sxprintf ("|   SCF electronic minimization\n");
   sxprintf ("|\n");

   SxRhoMixer mixer (cmd, nSpin == 2);
   mixer.print ();

   cout << SX_SEPARATOR;  cout.flush ();

   if (!calc)  return;

   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;
   SxGkBasis &Gk = waves.getGkBasis();
   int         ik, nk = waves.getNk();
   SxList<int> nPerK;

   SX_CHECK (Gk.getNk() == waves.getNk(), Gk.getNk(), waves.getNk());

   // --- setup n123inv in G
   if (G.n123inv.getSize() == 0)  G.n123invSetup ();

   // --- use full, within our numerics complete, basis
   for (ik = 0; ik < nk; ik++)  nPerK << Gk(ik).ng;

   // --- diagonalization loop
   int                         it;
   int                         i, n, ng, iSpin;
   PrecEnergy                  eTot = 0., eTotOld = 0.;
   SxDiracSymMat<TPrecCoeffG>  HComplete;
   SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;

   for (it = 0; it < maxSteps; it++)  {
      mixer.addRhoIn (H.rho);

      // --- direct diagonalization H |psi> = e |psi>
      for (ik = 0; ik < nk; ik++)  {
         ng = nPerK(ik);

         for (iSpin = 0; iSpin < nSpin; iSpin++)  {
            {
               SX_CLOCK(Timer::SubspaceMatrix);
               HComplete = H.getMatrix (Gk, iSpin, ik, ng);
            }
            { 
               SX_CLOCK(Timer::SubspaceDiag);
               eig = HComplete.eigensystem ();  // TODO: 'true' to sort values?
            }

            // --- get energies
            n = (int)fermi.eps(iSpin,ik).getSize ();
            fermi.eps(iSpin,ik) <<= eig.vals( SxIdx(0,n-1) );
            SX_CHECK (waves(iSpin,ik).nCols() == n,
                      waves(iSpin,ik).nCols(), n);

            // --- store waves, '<<' fills trailing part with zeroes
            for (i = 0; i < n; i++)  {
               waves(iSpin,ik).colRef(i) << eig.vecs.colRef(i);
            }  // :i
         }  // :iSpin
      }  // :ik

      // --- fermi occupation
      if (!keepOcc)  fermi.fermiDistribution (ekt);

      if (!keepRho)  {
         // --- update rho
         H.computeRho (fermi.focc, waves);

         // --- mixing
         mixer.addRhoOut (H.rho);
         H.rho = mixer.getMixedRho ();
      }

      // --- total energy
      H.update (fermi);
      eTot = H.eTotal;
      sxprintf ("| Total energy (%d): %15.12g\n", it+1, eTot);  cout.flush ();

      // --- exx
      if (exx)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| Direct diagonalisation: done.\n");
         break;
      }

      // --- convergence?
      if (fabs(eTotOld - eTot) < dEnergy)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| Direct diagonalisation: convergence reached.\n");
         break;
      }

      if (it >= maxSteps-1)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| WARNING: Maximum number of steps exceeded.\n");
         sxprintf ("|          Convergence not yet reached.\n");
         break;
      }

      // --- dump detailed info
      if ( !(it % printSteps) )  {
         H.printEnergies ();
         fermi.printOccupation (true);
      }

      eTotOld = eTot;
   }  // :it
   energy = eTot - ekt * fermi.getEntropy ();
}

void SxHamSolver::linGradTest (const SxSymbolTable *cmd, bool calc)
{
   SX_NO_MPI;
   // TODO: use hamPtr here...
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ())); 
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr); 
   SX_CHECK (dynamic_cast<SxPW*>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW&>(*wavesPtr);
   SxGkBasis Gk = waves.getGkBasis();
   int i, iSpin, ik, nSteps;
   double lambdaMin, lambdaMax;
   try  {
      i         = cmd->get("i")->toInt ();
      iSpin     = cmd->get("iSpin")->toInt ();
      ik        = cmd->get("ik")->toInt ();
      lambdaMin = cmd->get("lambdaMin")->toReal ();
      lambdaMax = cmd->get("lambdaMax")->toReal ();
      nSteps    = cmd->get("nSteps")->toInt ();
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }

   cout << SX_SEPARATOR;
   cout << "| Linearized gradient test\n";
   cout << SX_SEPARATOR;
   cout << "|    State (" << i << "," << iSpin << "," << ik << ")\n";
   cout << "|    Range: " << lambdaMin << " - " << lambdaMax << endl;
   cout << SX_SEPARATOR;
   cout.flush ();


   if (!calc)  return;

   if (!(H.contrib & SxPWHamiltonian::CALC_RHO))  {
      sxprintf ("ERROR: Make sure that CALC_RHO is part of hContrib.\n");
      SX_QUIT;
   }

   // --- create header for output file
   try { SxFSAction::rm ("lingrad.dat"); } catch (const SxException&) { /* empty */ }
   SxFileIO::appendToFile ("# lambdaT  eTot[H]   eLin[H]\n", "lingrad.dat");

   double eLin, eLin0, eGrd0, eps0, dLambda = (lambdaMax - lambdaMin)/nSteps;
   PsiG psi0, g, X;

   // --- remove all states other than (i,iSpin,ik) from waves
   SX_CHECK(&Gk == &(waves.getGkBasis()));
   H.rho.nElectrons = 1.;
   int j, jSpin, jk;
   for (jk=0; jk < waves.getNk (); ++jk)  
      for (jSpin=0; jSpin < waves.getNSpin(); ++jSpin)  
         for (j=0; j < waves.getNStates(); ++j) 
               fermi.focc(j,jSpin,jk) = 0.;

   fermi.focc(i,iSpin,ik) = 1.;
   Gk.weights(ik) = 1.;
   H.computeRho (fermi.focc, waves);
   H.update (fermi);
   

   // --- get gradient vector and y-offset at origin
   psi0.copy (waves(i,iSpin,ik));
   g = H | psi0;
   eps0  = (psi0 ^ g).chop ().re;
   eLin0 = 2. * eps0;
   eGrd0 = H.eTotal;
   sxprintf ("lin0 = %g, grd0 = %g\n", eLin0, eGrd0);

   // --- loop over lambda
   for (double lambda = lambdaMin; lambda <= lambdaMax; lambda += dLambda)  {

      // --- search direction along the steepest descent gradient
      X  = g - psi0*eps0;

      waves(i,iSpin,ik) <<= psi0 - lambda * X;
      H.computeRho (fermi.focc, waves.orthonormalize ());
      // TODO: should be compute (...)
      H.getEnergy (waves, fermi);
      eLin = 2. * (g ^ waves(i,iSpin,ik)).chop().re;

      // --- write to stdout and output file
      cout << SX_CLEARLINE << "eTotal(l=" << lambda << ") = " << H.eTotal << endl;
      cout.flush ();
      SxFileIO::appendToFile
      (  SxString(lambda,           "%15.12f") + "\t"
       + SxString(eLin - eLin0,     "%15.12f") + "\t" 
       + SxString(H.eTotal - eGrd0, "%15.12f") + "\n"
      , "lingrad.dat");

   }
   
}

const SxPW &SxHamSolver::getWaves () const
{
   SX_CHECK (dynamic_cast<const SxPW *> (wavesPtr.getPtr ()));
   return dynamic_cast<const SxPW &> (*wavesPtr);
}

void SxHamSolver::extrapolateWaves (double alpha,  const SxPW &oldWaves)
{ 
   SX_NO_MPI;
   SX_CHECK (dynamic_cast<SxPW *> (wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &> (*wavesPtr);
   int iStates, iSpin, ik, nk;

   SX_CHECK (nSpin == waves.getNSpin (), nSpin, waves.getNSpin ());
   nk      = waves.getNk ();
   SX_CHECK (nStates = waves.getNStates (), nStates, waves.getNStates ());

   for (ik = 0; ik < nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (iStates=0; iStates < nStates; iStates++)  {
            waves(iStates, iSpin, ik).set ( 
                  (1. + alpha) *    waves (iStates, iSpin, ik)
                      - alpha  * oldWaves(iStates, iSpin, ik)  );
         }      
      }
   }

   waves.orthonormalize (Loewdin);

}

SxAtomicStructure SxHamSolver::getStructure () const
{
   return structure;
}


void SxHamSolver::write (const SxString &filename) const
{
   // --- write waves to disk
#if defined USE_LOOPMPI and !defined USE_PARALLEL_NETCDF4
   if (SxLoopMPI::nr () > 1)  {
      SX_MPI_MASTER_ONLY {
         cout << "Writing waves is omitted since Sx was compiled without support for parallel NetCDF4 IO." << endl;
      }
   } else
#endif
   wavesPtr->writeWavesFile(filename, fermi, structure, true);
}

void SxHamSolver::printParameters (const SxString &prefix,
                                   const SxString &title) const
{
   cout << SX_SEPARATOR;
   cout << prefix << title << endl;
   cout << SX_SEPARATOR;
   if ( fabs(deltaT) > 1e-10)
      cout << prefix << "   time step:    " << deltaT << endl;
   if ( fabs(ekt) > 1e-10)
      cout << prefix << "   elect. temp.: " << ekt*HA2EV << " eV" << endl;
   if ( fabs(rhoMixing) > 1e-10)  
      cout << prefix << "   rho mixing:   " << rhoMixing  << endl;
   if ( fabs(foccMixing) > 1e-10)
      cout << prefix << "   Fermi mixing: " << foccMixing  << endl;
   if ( fabs(dEnergy) > 1e-10)
      cout << prefix << "   eTot conv.:   " << dEnergy << " H" << endl;
   if ( fabs(dEps) > 1e-10)
      cout << prefix << "   eps conv.:    " << dEps*HA2EV << " eV" << endl;
   if ( fabs(dPsiConv) > 1e-10)
      cout << prefix << "   psi conv.:    " << dPsiConv << endl;
   if ( maxSteps > 0)
      cout << prefix << "   max. steps:   " << maxSteps  << endl;
   cout << SX_SEPARATOR;
}
