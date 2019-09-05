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

#include <SxTBHamSolver.h>
#include <SxNeighbors.h>
#include <SxFileIO.h>

// ref[1]: Phys. Stat. Sol. (b) 217, 41 (2000)
// ref[2]: Phys. Rev. B 58, 11 (1998)
// ref[3]: J. Comput. Chem. 24, 565-581 (2003)

SxTBHamSolver::SxTBHamSolver (const SxAtomicStructure &str,
                              const SxSymbolTable *table)
   : structure(str,SxAtomicStructure::Copy),
     data (table),
     kPoints (str.cell, table)
{
   // periodic ?
   periodic = structure.isPeriodic ();

   //if (periodic)
   //   kPoints = SxKPoints (str.cell, table);
   //else {
   //   cout << endl;
   //   cout << " NON PERIODIC STRUCTURE: GAMMA K-POINT CALCULATIONS " << endl;
   //   kPoints = SxKPoints ());
   //}

   str.print (data);
   firstTime = true;

   // clean up files produced by SxTBHamSolver
   SxBinIO::deleteFile ("energy-tb.dat");
}

SxTBHamSolver::~SxTBHamSolver ()
{
   // empty
}

bool SxTBHamSolver::isRegistered (const SxSymbolTable *cmd) const
{
   if (!cmd) return false;
   SxString str = cmd->getName();
   return (  str == "tbInitialGuess" || str == "DIAG" );
}

void SxTBHamSolver::execute (const SxSymbolTable *cmd, bool calc)
{
   SxString str = cmd->getName();
   SxString indent = " "; // = cmd->level * ' ';
   SxString prefix = "|" + indent;

   if      (str == "tbInitialGuess")   tbInitialGuess (cmd, calc);
   else if (str == "DIAG")             diagonalize  (cmd, calc);
   else {
      SX_EXIT;
   }

   if (calc)  {  // write data to disc
      cout << SX_SEPARATOR;
      cout << "| Saving data\n";
      cout << SX_SEPARATOR;
      cout << "|   Eigenspectrum...    "; cout.flush();
      fermi.writeSpectrum ("eps-tb", "dat"); cout << "done\n";
      cout << SX_SEPARATOR;
      if (str != "tbInitialGuess") 
         printTiming (/* restart= */true,true); 
      else printTiming (/* restart= */true,false);
   }
}


void SxTBHamSolver::tbInitialGuess (const SxSymbolTable *cmd, bool calc)
{
   SX_CLOCK (Timer::initialization);
   cout << SX_SEPARATOR;
   sxprintf ("| This is SxTBHamSolver::tbInitialGuess\n");

   // --- read input data from cmd
   try  {
      withoutERepulsive = (cmd->contains("withoutERepulsive"))
                        ? cmd->get("withoutERepulsive")->toAttribute ()
                        : false;
   }  catch (SxException e)  {
      e.print();
      SX_EXIT;
   }

   // --- print input data
   if (withoutERepulsive)  {
      cout << SX_SEPARATOR;
      cout << "| WARNING :\n";
      cout << "| withoutERepulsive flag will cause that the\n";
      cout << "| Repulsive Energy coefficients will not be read from files\n"
           << "| and Repulsive Energy will be set to ZERO everywhere \n"; 
   }

   if (!calc)  return;
   
   // --- perform calculations

   // tbInitialGuess reads slater-koster files and 
   // initializes tbAtomAtom objects
   
   nSpecies = structure.getNSpecies();
   nTlAtoms = structure.getNAtoms();
   nAtoms.resize (nSpecies);
   int iS, jS; 
   
   for (iS = 0; iS < nSpecies; iS++)  { 
      nAtoms(iS) = structure.getNAtoms (iS);
   }
   
   tbAtomAtom.resize (nSpecies);
      for (iS = 0; iS < nSpecies; iS++)
      tbAtomAtom(iS).resize (nSpecies);
   
   cout << SX_SEPARATOR;

   // --- constructing slater-koster files name using skElementName from
   // --- TBSpecies data, and then initializing tbAtomAtom objects
   for (iS = 0; iS < nSpecies; iS++)  {
      for (jS = 0; jS < nSpecies; jS++)  {
         tbAtomAtom(iS)(jS).init (iS, jS, data, withoutERepulsive); 
      }
   }
   cout << SX_SEPARATOR;
   
   // initialize indices, and find lOrbMax
   initIdx (structure);
   
   // --- number of electrons
   nElect.resize (nSpecies);
   nTlElect = 0;  

   for (iS = 0; iS < nSpecies; iS++)  { 
      nElect(iS) = tbAtomAtom(iS)(iS).getNElect();
      nTlElect  += nAtoms(iS)*nElect(iS);
      if (tbAtomAtom(iS)(0).lMax1 != data.lMax(iS)) {
         cout << "| WARNING: lMax in species data for " << data.skElementName(iS) << "\n"
              << "|          is different from the value provided in the SK-files.\n" 
              << "|          lmax = " << tbAtomAtom(iS)(0).lMax1 << " versus " 
                                      << data.lMax(iS) << endl;
         cout << "|          Only the value in the SK-files is considered." << endl;
      }
   }

   // --- Gamma k-point ?
   isGammaPoint = (kPoints.nk == 1) && (kPoints.getK(0).norm() < 1e-8); 

   // --- other initialization 
   eTotal     = 0.; 
   eBand      = 0.;
   eRepulsive = 0.;
   // --- warning! only spin compensated: nSpin = 1. TODO
   fermi      = SxFermi (nTlElect, lOrbMax, 1, kPoints);

   // --- find the cutoff radius for neighbours 
   double radius;
   cutoff = 0.;
   for (iS = 0; iS < nSpecies; iS++)  { 
      for (jS = 0; jS < nSpecies; jS++)  {
         radius = tbAtomAtom(iS)(jS).getNeighborsCutoff();
         if (radius > cutoff) cutoff = radius;
      }
   }

   // R & G spheres not initialized yet 
   spheresInit = false;
}

void SxTBHamSolver::diagonalize (const SxSymbolTable *cmd, bool calc)
{
   SX_CLOCK (Timer::tbElMinim);
   cout << SX_SEPARATOR;
   sxprintf ("| This is SxTBHamSolver::diagonalize\n");

   SX_CHECK (calc == true);

   isSCF = false; 

   // --- read input data from cmd
   try  {
      ekt = cmd->get("ekt")->toReal();
   }  catch (SxException e)  {
      e.print();
      SX_EXIT;
   }
   SxSymbolTable *scfGroup = NULL;
   SxString mixingStr;
   if (cmd->containsGroup("SCF"))  {
      isSCF = true;
      scfGroup = cmd->getGroup("SCF");
      int mixingMethod;
      try  {
         maxSteps     = scfGroup->get("maxSteps"  )->toInt();
         printSteps   = scfGroup->get("printSteps")->toInt();
         dCharge      = scfGroup->get("dCharge"   )->toReal();
         mixingMethod = scfGroup->get("mixingMethod")->toInt();
         nMixingSteps = scfGroup->get("nMixingSteps")->toInt();
         rhoMixing    = scfGroup->get("rhoMixing" )->toReal();
         switch (mixingMethod)  {
         case 0 : mixerType = SxRhoMixer::Linear;
                  mixingStr = "linear mixer";
                  break;
         case 1 : mixerType = SxRhoMixer::Pulay;
                  mixingStr = "Pulay mixer";
                  cout << " Only linear mixing can be used for now !! " << endl;
                  SX_EXIT;
                  break;
         default: sxprintf ("Unknown mixing method: %d\n", mixingMethod);
                  SX_QUIT;
         }
      }  catch (SxException e)  {
         e.print();
         SX_EXIT;
      }
   } 

   // --- print input data
   cout << SX_SEPARATOR;
   cout << "| Diagonalization of the Tight-Binding Hamiltonian\n"; 
   if (isGammaPoint) cout << "|   Gamma-point calculations\n"; 
   if (isSCF)   {
      cout << "|   self-consistent tight-binding\n";
      cout << "|     details printed every    : " << printSteps 
                                                  << " iteration(s)\n";
      cout << "|     maximum iterations steps : " << maxSteps 
                                                  << " iteration(s)\n";
      cout << "|     mixer                    : " << mixingStr << endl;
      cout << "|     mixing factor            : " << rhoMixing << endl;
      cout << "|     mixing steps             : " << nMixingSteps << endl;
      cout << "|     delta charge convergence : " << dCharge   << endl;
   }  else  {
      cout << "|   non-self-consistent tight-binding\n";
   }
   cout << SX_SEPARATOR;

   if (!calc)  return;

   // --- perform calculation
      
   if (isSCF && periodic && !spheresInit)  {
      initRadii ();
      initRSphere(); 
      initGSphere(); 
      spheresInit = true;
      cout << SX_SEPARATOR;
   }

   // set up the grid for calculating neighbors
   grid = SxGrid (structure, 20);

   // --- diagonalizer
   if (isGammaPoint) diagonalizeGamma();
   else              diagonalizeCmplx(); 

   // clean the grid ?

}

void SxTBHamSolver::diagonalizeGamma ()
{
   SxArray<SxVector<Double> > rho(nTlAtoms);
   SxArray<SxVector<Double> > rhoNew(nTlAtoms), rhoOld(nTlAtoms);
   SxMatrix<Double>           HG(lOrbMax, lOrbMax);
   SxVector<Double>           rhoVec;
   double residue;
   bool   converged = false;
   int it, iA; 

   // calculate the non-self-consistent H and S
   calcHG0AndSG(); 
    
   // --- diagonalize and set fermi distribution
   eigG     = diag (HG0, SG);
   setFermi (eigG.vals);

   eRepulsive = getRepulsiveEnergy(); // doesn't change during iterations 
   eBand      = fermi.getEBand (SxFermi::UseFocc);
   eTotal     = eBand + eRepulsive;   // non-scf tot. energy
   
   if (isSCF)  {
      sxprintf("| NON-SCF eBand  = %12.9f\n", eBand);
      cout << "| TODO: Pulay mixing " << endl;
      
      // --- calculate initial mulliken and atomic charges 
      atomicRho   = getAtomicRho ();
      if (firstTime) {
         mullikenRho = getMullikenRhoG();
         rho         = SxArray<SxVector<Double> > (mullikenRho);
         firstTime   = false;
      } else {
         rho = SxArray<SxVector<Double> > (rhoMemory);
      }
      // --- calculate gammas
      if (periodic) gammaMatrix = getGammaMatPeriodic ();
      else          gammaMatrix = getGammaMatNonPeriodic ();

      // --- define the charge mixer
      //SxRhoMixer mixer (mixerType, rhoMixing, nMixingSteps);
      //mixer.setNorm(nTlElect,1.0);
      
      cout << "|\n| start iterations :\n" << endl; 
      // --- start iterations
      for (it = 0; it < maxSteps; it++)  {
         //mixer.addRhoIn (putRhoInVec (rho));

         calcHG1 (rho);           // Hamiltonian elements modification
         HG   = HG0 + HG1;      
         eigG = diag (HG, SG);    // modefy eiigensystem 
         setFermi (eigG.vals);
         
         mullikenRho = getMullikenRhoG();
         rhoOld      = SxArray<SxVector<Double> > (rho); 
         rhoNew      = SxArray<SxVector<Double> > (mullikenRho);
         
         eCorrection = getECorrection();
         eBand       = fermi.getEBand (SxFermi::UseFocc) - eCorrection;
         eTotal      = eBand  + eRepulsive;
         
         residue = 0;
         for (iA = 0; iA < nTlAtoms; iA++)
            residue += fabs(  rhoNew(iA).absSqr().sum() 
                            - rhoOld(iA).absSqr().sum() );
         
         // --- print out, and to file energy-tb.dat 
         cout << "iteration = " << it+1 << " eTotal = " << eTotal
              << " eBand = " << eBand << " rho residue = " << residue 
              << endl << endl;
         SxFileIO::appendToFile
         (  SxString(it+1)                               + "\t" // iteration
          + SxString(GETTIME(Timer::tbElMinim),"%12.8f") + "\t" // time
          + SxString(eTotal,"%15.12f")                   + "\t" // E_tot  [H]
          + SxString(eBand ,"%15.12f")                   + "\n" // E_band [H]
         , "energy-tb.dat");
         if ( !((it+1) % printSteps) )  {
            cout << SX_SEPARATOR;
            sxprintf("| Band structure energy = %12.9f\n", eBand);
            sxprintf("| Repulsive energy      = %12.9f\n", eRepulsive);
            sxprintf("| TOTAL ENERGY          = %12.9f\n", eTotal);
            cout << SX_SEPARATOR;
            fermi.printOccupation ();
         }
         
         // --- check convergence 
         if (residue < dCharge)  {
            converged = true;
            cout << SX_SEPARATOR;
            cout << "| convergence reached\n";
            cout << SX_SEPARATOR;
            break; 
         }
         
         // --- mix and get mixed rho
         //mixer.addRhoOut (putRhoInVec(rhoNew));
         //rhoVec = toVector (mixer.getMixedRho()(0));
         //rho    = putRhoInArrayOfVec (rhoVec); 
         rhoVec = (1.0 - rhoMixing)*putRhoInVec(rhoOld) 
                + rhoMixing*putRhoInVec(rhoNew);
         rho = putRhoInArrayOfVec (rhoVec); 
      } // iterations   

      if (!converged)  { 
         cout << " WARNING:\n";
         cout << "   the maximum number of iterations has been reached\n";
         cout << "   without achieving convergence !!! \n" ;
      }
      rhoMemory  = SxArray<SxVector<Double> > (rho);

      // --- anyway, print Mulliken charges to a file 
      printMullikenRho ();

   } // isSCF

   cout << SX_SEPARATOR;
   sxprintf("| Band structure energy = %12.9f\n", eBand);
   sxprintf("| Repulsive energy      = %12.9f\n", eRepulsive);
   sxprintf("| TOTAL ENERGY          = %12.9f\n", eTotal);
   cout << SX_SEPARATOR;
   fermi.printOccupation();

}

void SxTBHamSolver::calcHG0AndSG () 
{
   SX_CLOCK (Timer::hamOvl);

   int iS, jS, iA, jA, iN, iOrb, jOrb; 
   int nOrb1, nOrb2, nNeighbors;
   SxVector3<Double> distVec;

   HG0.reformat (lOrbMax, lOrbMax);
   SG.reformat  (lOrbMax, lOrbMax);
   HG0.set(0);
   SG.set(0);

   SxNeighbors neighbors;

   // --- calculating Hamiltonian (H) and Ovelrap (S) matrix elements
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS         = structure.getISpecies (iA);
      neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                         ,SxNeighbors::StoreRel);
      nNeighbors = neighbors.getSize();
      for (iN = -1; iN < nNeighbors; iN++)  { 
         if (iN == -1)  {   // the atom itself 
            jA = iA;
            jS = iS;
            tbAtomAtom(iS)(jS).setDiagElements();
            nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
            nOrb2 = nOrb1;
         }  else { 
            jA = neighbors.idx(iN);
            jS = structure.getISpecies (jA);
            distVec = neighbors.relPositions(iN);
            tbAtomAtom(iS)(jS).setDistVec (distVec);
            tbAtomAtom(iS)(jS).setElements();
            nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
            nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
         }
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  { 
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
               HG0 (idx(iA)+iOrb, idx(jA)+jOrb) 
                  += tbAtomAtom(iS)(jS).getHam (iOrb, jOrb); 
               SG  (idx(iA)+iOrb, idx(jA)+jOrb)
                  += tbAtomAtom(iS)(jS).getOvl (iOrb, jOrb); 
            }
         }
      }
   }

}
void SxTBHamSolver::calcHG1 (SxArray<SxVector<Double> > &rhoIn) 
{
   SX_CLOCK (Timer::hamMod);

   int iA, jA, iS, jS, iOrb, jOrb;
   int nOrb1, nOrb2;
   int kA, kOrb, kS, nOrb3;
   double sum;

   HG1.reformat (lOrbMax, lOrbMax);
   HG1.set(0);
   HG1Part1.reformat (lOrbMax, lOrbMax);
   HG1Part1.set(0);

   // rhoIn is mulliken charges. Rho here should be :
   SxArray<SxVector<Double> > rho = rhoIn - atomicRho;

   // see ref[1], 2nd term of eq(34) we call H1_{mu,nu}, this routin 
   // stores also part of HG1 , for convenience 
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS = structure.getISpecies (jA);
         nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
               sum = 0;
               for (kA = 0; kA < nTlAtoms; kA++)  { 
                  kS    = structure.getISpecies (kA);
                  nOrb3 = tbAtomAtom(kS)(0).getNOrb1();
                  for (kOrb = 0; kOrb < nOrb3; kOrb++)  {  
                     sum += rho(kA)(kOrb) * 
                           (  gammaMatrix (idx(iA)+iOrb, idx(kA)+kOrb)
                            + gammaMatrix (idx(jA)+jOrb, idx(kA)+kOrb) );
                  }
               }
               HG1Part1 (idx(iA)+iOrb, idx(jA)+jOrb) += sum;
            }
         }
      }
   }

   HG1Part1 *= 0.5;
   HG1 = HG1Part1 * SG;
}

SxSymMatrix<Double>::Eigensystem SxTBHamSolver::
                     diag (SxMatrix<Double> &HIn, SxMatrix<Double> &SIn) const
{
   SX_CLOCK (Timer::diagonalization);
   int iOrb, jOrb; 
   SxSymMatrix<Double>::Eigensystem eigIn;

   // --- Diagonalization det(H - S*Ei) = 0, this is done via 
   // --- det ((L^dagger . H . L) - Ei * 1) = 0, where S^-1 = L . L^dagger
   SxMatrix<Double>    HInNew, SInInv, L;
   SxSymMatrix<Double> HInSym(lOrbMax);
   SInInv = SIn.inverse();
   L      = SInInv.choleskyDecomposition();//S^-1 = L . L^dagger
   HInNew = (L.adjoint() ^ HIn ^ L);
   
   for (iOrb = 0; iOrb < lOrbMax; iOrb++)
      for (jOrb = iOrb; jOrb < lOrbMax; jOrb++)
         HInSym(iOrb,jOrb) = HInNew(iOrb,jOrb);
   eigIn = HInSym.eigensystem();
   
   eigIn.vecs = L ^ eigIn.vecs;
   
   return eigIn;
}

void SxTBHamSolver::setFermi (SxVector<Double> &eigenvalues) 
{
   int iOrb; 
   for (iOrb = 0; iOrb < lOrbMax; iOrb++) {
      fermi.eps(iOrb, 0, 0) = eigenvalues(iOrb);
   }
   //fermi.fermiDistribution (ekt, fermi.NO_DAMPING);
   fermi.fermiDistribution (ekt);
}

void SxTBHamSolver::setFermi () 
{
   int iOrb, ik;
   for (ik = 0; ik < fermi.getNk(); ik++)  {
      for (iOrb = 0; iOrb < lOrbMax; iOrb++) {
         fermi.eps(iOrb, 0, ik) = eig(ik).vals(iOrb);
      }
   }
   //fermi.fermiDistribution (ekt, fermi.NO_DAMPING);
   fermi.fermiDistribution (ekt);
}


void SxTBHamSolver::diagonalizeCmplx ()
{
   SxArray<SxVector<Double> > rho(nTlAtoms);
   SxArray<SxVector<Double> > rhoNew(nTlAtoms), rhoOld(nTlAtoms);
   
   SxArray<SxMatrix<Complex16> > H(fermi.getNk());

   int    ik, it, iA; 
   bool   converged = false;
   double residue;
   SxVector<Double> rhoVec;

   eig.resize (fermi.getNk());

   // calculate the non-SCF H and S
   calcH0AndS ();
   
   // --- diagonalize and set Fermi distribution
   for (ik = 0; ik < fermi.getNk(); ik++)  {
      eig(ik) = diag (H0(ik),S(ik));
   }
   setFermi();
   fermi.printOccupation();  

   eRepulsive = getRepulsiveEnergy(); // doesn't change during iterations 
   eBand      = fermi.getEBand (SxFermi::UseFocc);
   eTotal     = eBand + eRepulsive;   // non-scf tot. energy

   if (isSCF)  {
      sxprintf("| NON-SCF eBand  = %12.9f\n", eBand);
      cout << "| TODO: Pulay mixing " << endl;
      
      // --- calculate initial mulliken and atomic charges 
      atomicRho   = getAtomicRho();
      if (firstTime) {
         mullikenRho = getMullikenRho();
         rho         = SxArray<SxVector<Double> > (mullikenRho);
         firstTime   = false;
      } else {
         rho = SxArray<SxVector<Double> > (rhoMemory);
      }

      // --- calculate gammas
      if (periodic) gammaMatrix = getGammaMatPeriodic ();
      else          gammaMatrix = getGammaMatNonPeriodic ();
      
      // --- define the charge mixer
      //SxRhoMixer mixer (mixerType, rhoMixing, nMixingSteps);
      //mixer.setNorm(nTlElect,1.0);
      
      for (ik = 0; ik < fermi.getNk(); ik++) 
         H(ik).reformat (lOrbMax,lOrbMax);
      
      cout << "|\n| start iterations :\n" << endl; 
      // --- start iterations
      for (it = 0; it < maxSteps; it++)  {
         //mixer.addRhoIn (putRhoInVec (rho));
         
         calcH1 (rho);           // Hamiltonian elements modification 
         for (ik = 0; ik < fermi.getNk(); ik++)  {
            H(ik)   = H0(ik) + H1(ik);      
            eig(ik) = diag (H(ik),S(ik));      // modefy eiigensystem
         }
         setFermi ();
         
         mullikenRho = getMullikenRho();
         rhoOld      = SxArray<SxVector<Double> > (rho); 
         rhoNew      = SxArray<SxVector<Double> > (mullikenRho);
         
         eCorrection = getECorrection();
         eBand       = fermi.getEBand (SxFermi::UseFocc) - eCorrection;
         eTotal      = eBand  + eRepulsive;
         
         residue = 0;
         for (iA = 0; iA < nTlAtoms; iA++)
            residue += fabs(  rhoNew(iA).absSqr().sum() 
                            - rhoOld(iA).absSqr().sum() );
         
         // --- print out, and to file energy-tb.dat 
         cout << "iteration = " << it+1 << " eTotal = " << eTotal
              << " eBand = " << eBand << " rho residue = " << residue 
              << endl << endl;
         SxFileIO::appendToFile
         (  SxString(it+1)                               + "\t" // iteration
          + SxString(GETTIME(Timer::tbElMinim),"%12.8f") + "\t" // time
          + SxString(eTotal,"%15.12f")                   + "\t" // E_tot  [H]
          + SxString(eBand ,"%15.12f")                   + "\n" // E_band [H]
         , "energy-tb.dat");
         if ( !((it+1) % printSteps) )  {
            cout << SX_SEPARATOR;
            sxprintf("| Band structure energy = %12.9f\n", eBand);
            sxprintf("| Repulsive energy      = %12.9f\n", eRepulsive);
            sxprintf("| TOTAL ENERGY          = %12.9f\n", eTotal);
            cout << SX_SEPARATOR;
            fermi.printOccupation ();
         }
         
         // --- check convergence 
         if (residue < dCharge)  {
            converged = true;
            cout << SX_SEPARATOR;
            cout << "| convergence reached\n";
            cout << SX_SEPARATOR;
            break; 
         }
         
         // --- mix and get mixed rho
         //mixer.addRhoOut (putRhoInVec(rhoNew));
         //rhoVec = toVector (mixer.getMixedRho()(0));
         //rho    = putRhoInArrayOfVec (rhoVec); 
         rhoVec = (1.0 - rhoMixing)*putRhoInVec(rhoOld) 
                + rhoMixing*putRhoInVec(rhoNew);
         rho = putRhoInArrayOfVec (rhoVec); 
      } // iterations   

      if (!converged)  { 
         cout << " WARNING:\n";
         cout << "   the maximum number of iterations has been reached\n";
         cout << "   without achieving convergence !!! \n" ;
      }
      rhoMemory  = SxArray<SxVector<Double> > (rho);
      
      // --- anyway, print Mulliken charges to a file 
      printMullikenRho ();
      
   } // isSCF
    
   cout << SX_SEPARATOR;
   sxprintf("| Band structure energy = %12.9f\n", eBand);
   sxprintf("| Repulsive energy      = %12.9f\n", eRepulsive);
   sxprintf("| TOTAL ENERGY          = %12.9f\n", eTotal);
   cout << SX_SEPARATOR;
   fermi.printOccupation();

}

SxAtomicStructure SxTBHamSolver::getForces (const SxAtomicStructure &tau,
                                            const SxSymbolTable *cmd)
{
   SX_CLOCK (Timer::forces);

   cout << SX_SEPARATOR;
   sxprintf ("| This is SxTBHamSolver::getForces\n");

   structure = tau;
   execute (cmd);

   fTotal = structure.getNewStr ();
   fTotal.set(Coord(0.0,0.0,0.0));

   // --- resize and initialize forces 
   hamForces.resize (nTlAtoms);
   repulsiveForces.resize (nTlAtoms);
   for (int iA = 0; iA < nTlAtoms; iA++)  { 
      hamForces(iA).set(0.);
      repulsiveForces(iA).set(0.);
   }

   // set up the grid for calculating neighbors
   grid = SxGrid (structure, 20);

   // --- do the calculation  
   if (isGammaPoint) calcForcesGamma ();
   else              calcForcesCmplx ();
   
   return fTotal;
}

void SxTBHamSolver::calcForcesGamma ()
{
   int iS, jS, iA, jA, iN; 
   int i, idx1, idx2;
   int nOrb1, nOrb2, nNeighbors;
   SxVector3<Double> deltaVec, distVec0, fVec, distVec;
   SxMatrix<Double>  ham1, ham2, ovl1, ovl2;
   SxMatrix<Double>  subPDens, subPergDens, subHG1Part1;
   SxMatrix<Double>  hamDiff, ovlDiff;
   double            gradient;

   SxMatrix<Double> pDens    = getPDensMat();
   SxMatrix<Double> pErgDens = getPErgDensMat();
     
   SxNeighbors neighbors;

   // --- calculating the forces
   deltaVec.set (0.01);
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS         = structure.getISpecies (iA);
      neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                         ,SxNeighbors::StoreRel);
      nNeighbors = neighbors.getSize();

      for (iN = 0; iN < nNeighbors; iN++)  { 
         jA    = neighbors.idx(iN);
         jS    = structure.getISpecies (jA);
         nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
         idx1  = idx(iA);
         idx2  = idx(jA);

         distVec0 = neighbors.relPositions(iN);
         tbAtomAtom(iS)(jS).setDistVec (distVec0);
            
         // --- repulsive forces 
         fVec    =(distVec0/distVec0.norm())
                 * tbAtomAtom(iS)(jS).getRepulsiveEGrad (distVec0.norm());
         repulsiveForces(iA) += fVec; 

         // --- Hamiltonian forces 
         for (i = 0; i < 3; i++)  {
            distVec = SxVector3<Double> (distVec0); 
            distVec(i) -= deltaVec(i); 
            tbAtomAtom(iS)(jS).setDistVec (distVec);
            tbAtomAtom(iS)(jS).setElements();
            ham1  = tbAtomAtom(iS)(jS).getAtomAtomHam();
            ovl1  = tbAtomAtom(iS)(jS).getAtomAtomOvl();

            distVec = SxVector3<Double> (distVec0); 
            distVec(i) += deltaVec(i); 
            tbAtomAtom(iS)(jS).setDistVec (distVec);
            tbAtomAtom(iS)(jS).setElements();
            ham2  = tbAtomAtom(iS)(jS).getAtomAtomHam();
            ovl2  = tbAtomAtom(iS)(jS).getAtomAtomOvl();
            
            subPDens    = extractMat (pDens,idx1,idx2,nOrb1,nOrb2);
            subPergDens = extractMat (pErgDens,idx1,idx2,nOrb1,nOrb2);
            hamDiff = ham2 - ham1;
            ovlDiff = ovl2 - ovl1;

            if (isSCF)  {
               subHG1Part1 = extractMat (HG1Part1,idx1,idx2,nOrb1,nOrb2);
               gradient    = ( ( (hamDiff*subPDens - ovlDiff*(subPergDens
                                - subPDens*subHG1Part1)) 
                             ) /(2.*deltaVec(i))) .sum(); 
               hamForces(iA)(i) += 2. * gradient;
            }  else  {
               gradient = ( (hamDiff*subPDens - ovlDiff*subPergDens)
                          /(2.*deltaVec(i)) ).sum(); 
               hamForces(iA)(i) += 2. * gradient;
            }
         }
      }
   }
   // --- contributions of scf
   if (isSCF)  {
      SxArray<SxVector3<Double> > deltaRhoForces;
      if (!periodic)  deltaRhoForces = getDRhoForcesNonPeriodic ();
      else            deltaRhoForces = getDRhoForcesPeriodic ();
      for (iA = 0; iA < nTlAtoms; iA++)
         hamForces(iA) += deltaRhoForces(iA);
   }

   // --- sum contributions
   for (iA = 0; iA < nTlAtoms; iA++) 
      fTotal.ref(iA) = repulsiveForces(iA) + hamForces(iA);
}


void SxTBHamSolver::calcForcesCmplx ()
{
   int iS, jS, iA, jA, iN, ik; 
   int i, idx1, idx2;
   int nOrb1, nOrb2, nNeighbors;
   SxVector3<Double>    deltaVec, distVec0, fVec, distVec, kVec;
   SxMatrix<Complex16>  ham1, ham2, ovl1, ovl2;
   SxMatrix<Complex16>  subPDens, subPergDens, subH1Part1;
   SxMatrix<Complex16>  hamDiff, ovlDiff;
   SxMatrix<Complex16>  pDens, pErgDens;
   double               gradient, arg;
   SxComplex16          phase;
   SxNeighbors          neighbors;

   // --- calculating the forces
   // --- 1) repulsive forces 
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS         = structure.getISpecies (iA);
      neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                         ,SxNeighbors::StoreRel);
      nNeighbors = neighbors.getSize();
      
      for (iN = 0; iN < nNeighbors; iN++)  { 
         jA       = neighbors.idx(iN);
         jS       = structure.getISpecies (jA);
         distVec0 = neighbors.relPositions(iN);
         tbAtomAtom(iS)(jS).setDistVec (distVec0);
         
         fVec =(distVec0/distVec0.norm())
              * tbAtomAtom(iS)(jS).getRepulsiveEGrad (distVec0.norm());
         repulsiveForces(iA) += fVec; 
      }
   } 
   
   // --- 2) Hamiltonian forces 
   deltaVec.set (0.01);
   for (ik = 0; ik < fermi.getNk(); ik++)  {
      kVec     = kPoints.getK (ik);
      pDens    = getPDensMat (ik);
      pErgDens = getPErgDensMat (ik);

      for (iA = 0; iA < nTlAtoms; iA++)  { 
         iS         = structure.getISpecies (iA);
         neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                            ,SxNeighbors::StoreRel);
         nNeighbors = neighbors.getSize();

         for (iN = 0; iN < nNeighbors; iN++)  { 
            jA    = neighbors.idx(iN);
            jS    = structure.getISpecies (jA);
            nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
            nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
            idx1  = idx(iA);
            idx2  = idx(jA);
            
            distVec0 = neighbors.relPositions(iN);
            
            arg   = (distVec0*kVec).sum(); 
            phase = SxComplex16 (cos (arg), -sin (arg));  // exp {-ik.R} !

            for (i = 0; i < 3; i++)  {
               distVec = SxVector3<Double> (distVec0); 
               distVec(i) -= deltaVec(i); 
               tbAtomAtom(iS)(jS).setDistVec (distVec);
               tbAtomAtom(iS)(jS).setElements();
               ham1  = tbAtomAtom(iS)(jS).getAtomAtomHam();
               ovl1  = tbAtomAtom(iS)(jS).getAtomAtomOvl();
               
               distVec = SxVector3<Double> (distVec0); 
               distVec(i) += deltaVec(i); 
               tbAtomAtom(iS)(jS).setDistVec (distVec);
               tbAtomAtom(iS)(jS).setElements();
               ham2  = tbAtomAtom(iS)(jS).getAtomAtomHam();
               ovl2  = tbAtomAtom(iS)(jS).getAtomAtomOvl();
               
               subPDens    = extractMat (pDens,idx1,idx2,nOrb1,nOrb2);
               subPergDens = extractMat (pErgDens,idx1,idx2,nOrb1,nOrb2);
               hamDiff     = ham2 - ham1;
               ovlDiff     = ovl2 - ovl1;

               if (isSCF)  {
                  subH1Part1 = extractMat (H1Part1,idx1,idx2,nOrb1,nOrb2);
                  gradient    = ( ( (hamDiff*subPDens - ovlDiff*(subPergDens
                                     - subPDens*subH1Part1)) 
                                ) /(2.*deltaVec(i))) .sum() * phase; 

                  hamForces(iA)(i) += 2. * kPoints.weights(ik) * gradient;
               }  else  {
                  gradient = ( (hamDiff*subPDens - ovlDiff*subPergDens)
                                /(2.*deltaVec(i)) ).sum() * phase; 
                  
                  hamForces(iA)(i) += 2. * kPoints.weights(ik) * gradient;
               }
            }
         }
      }
   } 
   // --- contributions of scf
   if (isSCF)  {
      SxArray<SxVector3<Double> > deltaRhoForces; 
      if (!periodic)  deltaRhoForces = getDRhoForcesNonPeriodic ();
      else            deltaRhoForces = getDRhoForcesPeriodic ();
      for (iA = 0; iA < nTlAtoms; iA++)
         hamForces(iA) += deltaRhoForces(iA);
   } 

   // --- sum contributions
   for (iA = 0; iA < nTlAtoms; iA++) 
         fTotal.ref(iA) = repulsiveForces(iA) + hamForces(iA);

}

SxArray<SxVector3<Double> > SxTBHamSolver::getDRhoForcesNonPeriodic ()
{
   int iS, jS, iA, jA;
   int iOrb, jOrb, nOrb1, nOrb2;

   SxVector3<Double> fVec, distVec;

   SxArray<SxVector3<Double> > deltaRhoF (nTlAtoms);
   for (iA = 0; iA < nTlAtoms; iA++) deltaRhoF(iA).set(0.);

   // --- calculate the derivative of gamma (analytical)
   SxMatrix<Double>  dGammaMatrix = getDGammaMatNonPeriodic ();

   // calculate deltaRho
   SxArray<SxVector<Double> > deltaRho = mullikenRho - atomicRho;

   // --- do the calculation 
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  {  
         if (jA != iA)  {
            jS      = structure.getISpecies (jA);
            distVec = structure.getAtom (jA) - structure.getAtom(iA); 
            nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
            nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
            fVec.set(0);
            for (iOrb = 0; iOrb < nOrb1; iOrb++)  {
               for (jOrb = 0; jOrb < nOrb2; jOrb++)  {
                  fVec += (distVec/(distVec.norm())) * deltaRho(iA)(iOrb)
                     *(  dGammaMatrix (idx(iA)+iOrb, idx(jA)+jOrb)
                           *deltaRho(jA)(jOrb) );
               }
            }
            deltaRhoF(iA) += fVec;
         }
      }
   }	
   
   return deltaRhoF;
}

SxArray<SxVector3<Double> > SxTBHamSolver::getDRhoForcesPeriodic ()
{
   int iS, jS, iA, jA;
   int iOrb, jOrb, nOrb1, nOrb2;

   SxVector3<Double> fVec, distVec;

   SxArray<SxVector3<Double> > deltaRhoF (nTlAtoms);
   for (iA = 0; iA < nTlAtoms; iA++) deltaRhoF(iA).set(0.);

   // --- calculate the derivative of gamma (analytical)
   SxArray<SxMatrix<Double> >  dGammaMatrix = getDGammaMatPeriodic();

   // calculate deltaRho
   SxArray<SxVector<Double> > deltaRho = mullikenRho - atomicRho;

   // --- do the calculation 
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  {  
         if (jA != iA)  {
            jS      = structure.getISpecies (jA);
            distVec = structure.getAtom (jA) - structure.getAtom(iA); 
            nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
            nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
            fVec.set(0);
            for (iOrb = 0; iOrb < nOrb1; iOrb++)  {
               for (jOrb = 0; jOrb < nOrb2; jOrb++)  {
                  for (int i=0; i < 3 ; i++)   { 
                     fVec(i) += deltaRho(iA)(iOrb)
                             *( dGammaMatrix(i) (idx(iA)+iOrb, idx(jA)+jOrb)
                                *deltaRho(jA)(jOrb) );
                  }
               }
            }
            deltaRhoF(iA) += fVec;
         }
      }
   }	
   
   return deltaRhoF;
}

void SxTBHamSolver::calcH0AndS ()
{
   SX_CLOCK (Timer::hamOvl);

   int iS, jS, iA, jA, iN, iOrb, jOrb, ik; 
   int nOrb1, nOrb2, nNeighbors;

   H0.resize (fermi.getNk());
   S.resize (fermi.getNk());

   double arg;
   SxVector3<Double> kVec, distVec;
   SxComplex16       phase;

   SxNeighbors neighbors;

   // --- calculating Hamiltonian (H) and Ovelrap (S) matrix elements
   for (ik = 0; ik < fermi.getNk(); ik++)  {
      kVec = kPoints.getK (ik);
      H0(ik).reformat (lOrbMax, lOrbMax);
      S(ik).reformat (lOrbMax, lOrbMax);
      H0(ik).set(0);
      S(ik).set(0);
      for (iA = 0; iA < nTlAtoms; iA++)  { 
         iS         = structure.getISpecies (iA);
         neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                            ,SxNeighbors::StoreRel);
         nNeighbors = neighbors.getSize();
         for (iN = -1; iN < nNeighbors; iN++)  { 
            if (iN == -1)  {   // the atom itself 
               jA = iA;
               jS = iS;
               tbAtomAtom(iS)(jS).setDiagElements();
               phase = SxComplex16 (1., 0);
               nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
               nOrb2 = nOrb1;
            }  else { 
               jA = neighbors.idx(iN);
               jS = structure.getISpecies (jA);
               distVec = neighbors.relPositions(iN);
               tbAtomAtom(iS)(jS).setDistVec (distVec);
               tbAtomAtom(iS)(jS).setElements();
               arg   = (distVec*kVec).sum(); 
               phase = SxComplex16 (cos (arg), sin (arg));  // exp {ik.r}
               nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
               nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
            }
            for (iOrb = 0; iOrb < nOrb1; iOrb++)  { 
               for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
                  H0(ik) (idx(iA)+iOrb, idx(jA)+jOrb) 
                     += tbAtomAtom(iS)(jS).getHam (iOrb, jOrb) * phase; 
                  S (ik) (idx(iA)+iOrb, idx(jA)+jOrb)
                     += tbAtomAtom(iS)(jS).getOvl (iOrb, jOrb) * phase; 
               }
            }
         }
      } 
   }

}

void SxTBHamSolver::calcH1 (SxArray<SxVector<Double> > &rhoIn) 
{
   SX_CLOCK (Timer::hamMod);

   int ik, iA, jA, iS, jS, iOrb, jOrb;
   int nOrb1, nOrb2;
   int kA, kOrb, kS, nOrb3;
   double sum;

   H1Part1.reformat (lOrbMax, lOrbMax);
   H1Part1.set(0);

   // rhoIn is mulliken charges. Rho here should be :
   SxArray<SxVector<Double> > rho = rhoIn - atomicRho;

   // see ref[1], 2nd term of eq(34) we call H1_{mu,nu}, this routin 
   // stores also part of H1, for convenience 
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS = structure.getISpecies (jA);
         nOrb1 = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2 = tbAtomAtom(iS)(jS).getNOrb2();
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
               sum = 0;
               for (kA = 0; kA < nTlAtoms; kA++)  { 
                  kS    = structure.getISpecies (kA);
                  nOrb3 = tbAtomAtom(kS)(0).getNOrb1();
                  for (kOrb = 0; kOrb < nOrb3; kOrb++)  {  
                     sum += rho(kA)(kOrb) * 
                        (  gammaMatrix (idx(iA)+iOrb, idx(kA)+kOrb)
                         + gammaMatrix (idx(jA)+jOrb, idx(kA)+kOrb) );
                  }
               }
               H1Part1 (idx(iA)+iOrb, idx(jA)+jOrb) += sum;
            }
         }
      }
   }

   H1Part1 *= 0.5;
   H1.resize(fermi.getNk());
   for (ik = 0; ik < fermi.getNk(); ik++)  {
      H1(ik).reformat (lOrbMax, lOrbMax);
      H1(ik) = H1Part1 * S(ik);
      isHermitian(H1(ik));
   }

}

SxSymMatrix<Complex16>::Eigensystem SxTBHamSolver::
            diag (SxMatrix<Complex16> &HIn, SxMatrix<Complex16> &SIn) const
{
   SX_CLOCK (Timer::diagonalization);
   int iOrb, jOrb; 
   SxSymMatrix<Complex16>::Eigensystem eigIn;

   // --- Diagonalization det(H - S*Ei) = 0, this is done via 
   // --- det ((L^dagger . H . L) - Ei * 1) = 0, where S^-1 = L . L^dagger
   SxMatrix<Complex16>    HInNew, SInInv, L;
   SxSymMatrix<Complex16> HInSym(lOrbMax);
   SInInv = SIn.inverse();
   L      = SInInv.choleskyDecomposition();//S^-1 = L . L^dagger
   HInNew = (L.adjoint() ^ HIn ^ L);
   
   for (iOrb = 0; iOrb < lOrbMax; iOrb++)
      for (jOrb = iOrb; jOrb < lOrbMax; jOrb++)
         HInSym(iOrb,jOrb) = HInNew(iOrb,jOrb);
   eigIn = HInSym.eigensystem();
   
   eigIn.vecs = L ^ eigIn.vecs;
   
   return eigIn;
}

double SxTBHamSolver::getRepulsiveEnergy ()
{
   int iS, jS, iA, jA, iN, nNeighbors; 
   SxVector3<Double> distVec;
   double repulsive = 0.0;
  
   SxNeighbors neighbors;

   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      neighbors.compute (grid,structure,structure.getAtom(iA),cutoff
                         ,SxNeighbors::StoreRel);
      nNeighbors = neighbors.getSize();
      for (iN = 0; iN < nNeighbors; iN++)  { 
         jA = neighbors.idx(iN);
         jS = structure.getISpecies (jA);
         distVec = neighbors.relPositions(iN);
         tbAtomAtom(iS)(jS).setDistVec (distVec);
         repulsive += tbAtomAtom(iS)(jS).getERepulsive (distVec.norm());
      }
   }
   repulsive /= 2.; // to eliminate double counting 

   return repulsive;
}

double SxTBHamSolver::getECorrection () const
{
   double eCorr = 0.;
   int iA, jA, iOrb, jOrb;
   int iS, jS, nOrb1, nOrb2;
   SxMatrix<Double>  subGammaMat;
   SxVector3<Double> distVec;

   SxArray<SxVector<Double> > deltaRho = mullikenRho - atomicRho;

   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS      = structure.getISpecies (jA);
         distVec = structure.getAtom(jA) - structure.getAtom(iA); 
         nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
         subGammaMat = extractMat (gammaMatrix,idx(iA),idx(jA),nOrb1,nOrb2);
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  {
               eCorr += 0.5*subGammaMat(iOrb,jOrb)*deltaRho(iA)(iOrb)
                      *( mullikenRho(jA)(jOrb) + atomicRho(jA)(jOrb) );
            }
         }
      }
   }
   return eCorr;
}


void SxTBHamSolver::initIdx (const SxAtomicStructure &str) 
{
   idx.resize (nTlAtoms);
   idx.set(0);
   int iS, iA, nOrb1;
   int index = 0;
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS      = str.getISpecies (iA);
      idx(iA) = index;
      nOrb1   = tbAtomAtom(iS)(iS).getNOrb1();
      index  += nOrb1; 
   }
   lOrbMax = index;
}


SxSpeciesData SxTBHamSolver::getSpeciesData () const
{
   return data;
}

PrecEnergy SxTBHamSolver::getEnergy () const
{
   return eTotal;
}

SxArray<SxVector<Double> > SxTBHamSolver::getAtomicRho () //:iAtom, iOrb
{
   int    iOrb, lMaxIA, l, ml;
   int    nOrbIA, iA, iS;
   double nElectIA;
   SxArray<SxVector<Double> > atomicRhoIn(nTlAtoms);

   // --- calaculate isolated atoms charges (atomicRho), per orbital
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS       = structure.getISpecies (iA);
      lMaxIA   = tbAtomAtom(iS)(0).lMax1;      // ugly !
      nOrbIA   = tbAtomAtom(iS)(0).getNOrb1(); // ugly !
      nElectIA = nElect(iS);

      atomicRhoIn(iA).resize(nOrbIA);
      atomicRhoIn(iA).set(0.);
      iOrb = 0;
      for (l = 0; l <= lMaxIA; l++)  {
         for (ml = 0; ml < (2*l + 1); ml++, iOrb++)  {
            if (nElectIA >= 2.0*(2.0*l + 1.0) )  { 
               atomicRhoIn(iA)(iOrb) += 2.0;
            } else  {
               atomicRhoIn(iA)(iOrb) += nElectIA/(2.0*l + 1.0);
            }
         }
         nElectIA -= 2.0*(2.0*l + 1.0);
      }
   }
  
   return atomicRhoIn;
}

SxArray<SxVector<Double> > SxTBHamSolver::getMullikenRhoG () //:iAtom, iOrb
{
   SX_CLOCK (Timer::mullikenRho);

   int    nOrbIA, iOrb, iA, idx1, iS;

   SxArray<SxVector<Double> > mullikenRhoIn(nTlAtoms);

   // calculate the density matrix
   SxMatrix<Double> densityMat = getPDensMat ();

   // --- Mulliken charges are the diagonal elements of PSMat
   SxMatrix<Double> PSMat = densityMat ^ SG; 

   // --- atoms charges (q) (Mulliken charges) per orbital. See ref[1] eq(33)
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      idx1   = idx(iA);
      iS     = structure.getISpecies (iA);
      nOrbIA = tbAtomAtom(iS)(0).getNOrb1(); // ugly !
      mullikenRhoIn(iA).resize(nOrbIA);
      mullikenRhoIn(iA).set(0.);
      for (iOrb = 0; iOrb < nOrbIA; iOrb++)  {
         mullikenRhoIn(iA)(iOrb) = PSMat(idx1+iOrb,idx1+iOrb);
      }
   }
      
   return mullikenRhoIn;
}

SxArray<SxVector<Double> > SxTBHamSolver::getMullikenRho () //:iAtom, iOrb
{
   SX_CLOCK (Timer::mullikenRho);

   int    nOrbIA, iOrb, iA, idx1, iS, ik;

   SxArray<SxVector<Double> >     mullikenRhoIn(nTlAtoms);
   SxArray<SxMatrix<Complex16> >  densityMat(fermi.getNk ());
   SxArray<SxMatrix<Complex16> >  PSMat(fermi.getNk ());

   for (ik = 0; ik < fermi.getNk(); ik++)  {
      densityMat (ik) = getPDensMat (ik);
      PSMat(ik) = densityMat (ik) ^ S(ik);
   }
   
   // --- atoms charges (q) (Mulliken charges) per orbital. See ref[1] eq(33)
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      idx1   = idx(iA);
      iS     = structure.getISpecies (iA);
      nOrbIA = tbAtomAtom(iS)(0).getNOrb1(); // ugly !
      mullikenRhoIn(iA).resize(nOrbIA);
      mullikenRhoIn(iA).set(0.);
      for (iOrb = 0; iOrb < nOrbIA; iOrb++)  {
         for (ik = 0; ik < fermi.getNk(); ik++)  {
            mullikenRhoIn(iA)(iOrb) += double (kPoints.weights(ik)*
                                       PSMat (ik)(idx1+iOrb,idx1+iOrb));
         }
      }
   }

   /*
   // old code
   // --- the long expression contains a second inner loop over atoms
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      idx1   = idx(iA);
      iS     = structure.getISpecies (iA);
      nOrbIA = tbAtomAtom(iS)(0).getNOrb1(); // ugly !
      mullikenRhoIn(iA).resize(nOrbIA);
      mullikenRhoIn(iA).set(0.);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         idx2     = idx(jA);
         jS       = structure.getISpecies (jA);
         nOrb1    = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2    = tbAtomAtom(iS)(jS).getNOrb2();
         for (ik = 0; ik < fermi.getNk(); ik++)  {
            // calculate the density matrix
            densityMat = getPDensMat (ik); // ugly! 

            subPDens = extractMat (densityMat,idx1,idx2,nOrb1,nOrb2);
            subS     = extractMat (S(ik),idx1,idx2,nOrb1,nOrb2);

            for (iOrb = 0; iOrb < nOrb1; iOrb++)  {
               for (jOrb = 0; jOrb < nOrb2; jOrb++)  {
                  //mullikenRhoIn(iA)(iOrb) += double (0.5*kPoints.weights(ik)*
                        //(  subPDens(iOrb,jOrb)*subS(iOrb,jOrb).conj()
                        //+ subPDens(iOrb,jOrb).conj()*subS(iOrb,jOrb) ));
                  mullikenRhoIn(iA)(iOrb) += double (kPoints.weights(ik)*
                          subPDens(iOrb,jOrb)*subS(iOrb,jOrb).conj());
               }
            }
         }
      }
   }
   */

   return mullikenRhoIn;
}

SxMatrix<Double> SxTBHamSolver::getPDensMat () const
{
   SxMatrix<Double> tempMat, densMat;
   SxVector<Double> focc = toVector (fermi.focc(0,0));// iSpin=0, ik=0 

   tempMat  = tempMat.identity (focc);
   densMat  = (eigG.vecs ^ tempMat ^ eigG.vecs.transpose()); 

   return densMat;
}

SxMatrix<Double> SxTBHamSolver::getPErgDensMat () const
{
   SxMatrix<Double> tempMat, densMat;
   SxVector<Double> focc = toVector (fermi.focc(0,0));// iSpin=0, ik=0
   
   tempMat  = tempMat.identity (focc*eigG.vals);
   densMat  = (eigG.vecs ^ tempMat ^ eigG.vecs.transpose()); 

   return densMat;
}

SxMatrix<Complex16> SxTBHamSolver::getPDensMat (int ik) const
{
   SxMatrix<Complex16> densMat;
   SxMatrix<Double>    tempMat;
   SxVector<Double>    focc = toVector (fermi.focc(0,ik));// iSpin=0

   tempMat  = tempMat.identity (focc);
   densMat  = (eig(ik).vecs ^ tempMat ^ eig(ik).vecs.adjoint()); 

   return densMat;
}

SxMatrix<Complex16> SxTBHamSolver::getPErgDensMat (int ik) const
{
   SxMatrix<Complex16> densMat;
   SxMatrix<Double>    tempMat;
   SxVector<Double>    focc = toVector (fermi.focc(0,ik));// iSpin=0
   
   tempMat  = tempMat.identity (focc*eig(ik).vals);
   densMat  = (eig(ik).vecs ^ tempMat ^ eig(ik).vecs.adjoint()); 

   return densMat;
}

//----------------------------------------------------------------------------
// gamma related parts
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// a) gamma for non-periodic systems (e.g. clusters)
//----------------------------------------------------------------------------

SxMatrix<Double> SxTBHamSolver::getGammaMatNonPeriodic () const
{
   SX_CLOCK (Timer::gammaMatrix);

   SxMatrix<Double>  gammaMat(lOrbMax,lOrbMax);
   SxVector3<Double> distVec;
   int iA, jA, iS, jS, nOrb1, nOrb2, iOrb, jOrb;

   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS      = structure.getISpecies (jA);
         nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
         distVec = structure.getAtom(jA) - structure.getAtom(iA); 
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
               gammaMat(idx(iA)+iOrb, idx(jA)+jOrb) 
                  = gammaSCF (iA,jA,iOrb,jOrb,distVec.norm());
            }
         }
      }
   }

   return gammaMat;
}

double SxTBHamSolver::gammaSCF (int iA, int jA, int iOrb, int jOrb,
                                double distance) const 
{
   double gamma, tauIA, tauJA;

   int iS = structure.getISpecies (iA);
   int jS = structure.getISpecies (jA);
   tauIA  = (16./5.) * tbAtomAtom(iS)(iS).getHubbardU(iOrb);
   tauJA  = (16./5.) * tbAtomAtom(jS)(jS).getHubbardU(jOrb);

   // see ref[1] eq(30), ref[3] eq(6)
   if (iA == jA)  {       // i.e. distance = 0 
      if (iOrb == jOrb) 
         gamma = tbAtomAtom(iS)(jS).getHubbardU(iOrb);  
      else 
         gamma = (tauIA*tauJA*(tauIA*tauIA + 3.*tauIA*tauJA + tauJA*tauJA))
                 /(2.*pow(tauIA + tauJA,3)); 
   }  else  {
      gamma = 1.0/distance - shortRange (tauIA,tauJA,distance);  
   }

   return gamma;
}

double SxTBHamSolver::shortRange (double tA, double tB, double distance) const
{
   SX_CHECK (distance > 1e-8, distance);

   double shortRangePart;

   if (fabs(tA - tB) < 1e-8) // if same species, or tA = tB  somehow!
      //This function is the limit of the 2nd one (in the else part) as tB->tA
      shortRangePart = (48. + distance*tA*(33. + distance*tA*(9.+distance*tA)))
                      /(48.*distance*exp(distance*tA)); 
   else 
      shortRangePart = exp(-tA*distance)*K (tA,tB,distance) 
                     + exp(-tB*distance)*K (tB,tA,distance);  

   return shortRangePart;
}

double SxTBHamSolver::K (double tA, double tB, double distance) const 
{
   SX_CHECK (fabs(tA - tB) > 1e-8, tA, tB);

   double k;
   // a function, see ref[3] eq(7)  
   k =  pow(tB,4)*tA/(2.*pow((tA*tA - tB*tB),2))
     - (pow(tB,6) - 3*pow(tB,4)*tA*tA)/(pow((tA*tA - tB*tB),3)*distance); 

   return k;
}

//----------------------------------------------------------------------------
// b) gamma for periodic systems (e.g. solids)
//----------------------------------------------------------------------------

SxMatrix<Double> SxTBHamSolver::getGammaMatPeriodic () const
{
   SX_CLOCK (Timer::gammaMatrix);

   SxMatrix<Double>  gammaMat(lOrbMax,lOrbMax);
   SxVector3<Double> distVec;
   int iA, jA, iS, jS, nOrb1, nOrb2, iOrb, jOrb;
   double shortRangePart, longRangePart;

   gammaMat.set(0);
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS      = structure.getISpecies (jA);
         nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
         distVec = structure.getAtom(jA) - structure.getAtom(iA); 

         // SUM_R { 1/|r - R| }
         longRangePart = longRangeSum (iA,jA,distVec);  

         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 

               // SUM_R {short range function}
               shortRangePart = shortRangeSum (iA,jA,iOrb,jOrb,distVec);

               gammaMat(idx(iA)+iOrb,idx(jA)+jOrb) 
                  = longRangePart - shortRangePart;
            }
         }
      }
   }

   return gammaMat;
}

double SxTBHamSolver::shortRangeSum (int iA, int jA, int iOrb, int jOrb,
                                     SxVector3<Double> distVec) const 
{
   int iS = structure.getISpecies (iA);
   int jS = structure.getISpecies (jA);
   double tauIA  = (16./5.) * tbAtomAtom(iS)(iS).getHubbardU(iOrb);
   double tauJA  = (16./5.) * tbAtomAtom(jS)(jS).getHubbardU(jOrb);

   int n1, n2, n3;
   int nMax = 50;
   int nMin = 3;
   int nR;
   SxVector3<Double> latVec;
   double lastTerm = 0;
   double rSum  = 0.;
   double r;
   
   nR = 0;
   
   do  {
      lastTerm = 0;
      for (n1 = -nR; n1 <= nR; n1++)  {
         for (n2 = -nR; n2 <= nR; n2++)  {
            for (n3 = -nR; n3 <= nR; n3++)  {
               // only R's for the new shell contributes. Previous R's 
               // contributions were added in the previous do loop, and so on
               if (nR == abs(n1) || nR == abs(n2) || nR == abs(n3) )  {
                  
                  latVec  = structure.cell ^ SxVector3<Int> (n1,n2,n3);

                  if (iA == jA) {        
                     r = latVec.norm();
                     if (r > 1e-8)  {    // exclud rLat = 0
                        lastTerm =  shortRange (tauIA,tauJA,r);
                        rSum    +=  lastTerm;
                     }  else  {         
                        // this will add : 
                        // limit_r->0 [ ShortRange(r) - 1/r ] 
                        // = - limit_r->0 [ gamma(r) ]
                        lastTerm = -gammaSCF (iA,jA,iOrb,jOrb,r); 
                        rSum    += lastTerm; 
                     }
                  }  else  {
                     r = (latVec + distVec).norm();
                     lastTerm =  shortRange (tauIA,tauJA,r);
                     rSum    +=  lastTerm;
                  }     
                     
               }
            }
         }
      }
      nR ++; 
   }  while (nR <= nMax && (nR <= nMin || fabs(lastTerm) > 1e-25) ); 

   return rSum;
}

double SxTBHamSolver::longRangeSum (int iA, int jA,
                                    SxVector3<Double> distVec) const 
{
   SX_CLOCK (Timer::ewald); // TODO
   int iLat; 
   int nRLat = (int)rLat.getSize (); 
   int nGLat = (int)gLat.getSize ();
   double rSum  = 0., gSum = 0;
   double ewald = 0.; 
   double r;
   double g2, arg;
   double factor = 4.*PI/structure.cell.volume;
   SxVector3<Double> gVec;
   SxComplex16       phase;

   // --- do the real space sum
   if (iA == jA) {
      for (iLat = 0; iLat < nRLat; iLat++ )  { 
         r = rLat(iLat).norm();  
         if (r > 1e-8)      // exclud rLat = 0
            rSum += derfc(eta*r)/r;
      }
   } else {
      for (iLat = 0; iLat < nRLat; iLat++ )  { 
         r     = (distVec + rLat(iLat)).norm();  
         rSum += derfc(eta*r)/r;
      }
   }

   // --- do the reciprocal space sum
   for (iLat = 0; iLat < nGLat; iLat++ )  { // G = 0 was excluded 
      gVec  = gLat(iLat);
      g2    = gVec.absSqr().sum();
      arg   = distVec ^ gVec;
      phase = SxComplex16 (cos(arg), sin(arg));
      gSum += (double)( (exp(-g2/(4.*eta*eta)) * phase)/g2 );
   }
   gSum *= factor;

   // --- add different contributions
   ewald = rSum + gSum - factor/(4.*eta*eta);
   if (iA == jA) ewald -= 2.*eta/sqrt(PI);
   
   return ewald; 
}

void SxTBHamSolver::initRadii()
{
   double tol = 1e-5;
   int    i;
   bool   tolRFound = false, tolGFound = false;
   double R, G, G2;
   
   // initial R & G that are related to the dimentions of the cell
   double initR = pow(structure.cell.volume,1./3.);
   double initG = 2*PI/initR;

   // smaller eta -> increase # of R vectors needede, and vise versa 
   eta = getEta ();

   // --- find largest R & G vectors to get a good convergence
   // --- in both real & reciprocal space sums of the ewald sum 
   // --- (which is calculated in longRangeSum ( )). 
   for (i=1; i<1000; i++)  {
      R = i * initR;
      if (fabs (derfc(eta*R)/R) <= tol)  {
         rEwald = R;
         tolRFound = true;
         break;
      }
   }
   for (i=1; i<1000; i++)  {
      G  = i * initG;
      G2 = G * G;
      if (fabs (exp(-G2/(4.*eta*eta))/G2) <= tol)  {
         gEwald = G;
         tolGFound = true;
         break;
      }
   }
   
   if (!tolRFound)  { 
      cout << " tolerance not acheived in R-space part of ewald sum !\n"; 
      SX_EXIT;
   }  else if (!tolGFound)  { 
      cout << " tolerance not acheived in G-space part of ewald sum !\n"; 
      SX_EXIT;
   }

   cout << "|  rEwald = " << rEwald << endl;
   cout << "|  gEwald = " << gEwald << endl;

}

double SxTBHamSolver::getEta ()
{
   double etaOut, shortR, shortG, difference;
   double eta1   = 1e-5;
   double eta2   = 10;
   int    nLoops = 0;

   const CellMat b (structure.cell.getReciprocalCell());

   shortR = minimum (structure.cell(0).norm(), structure.cell(1).norm(),
                     structure.cell(2).norm());
   shortG = minimum (b(0).norm(), b(1).norm(), b(2).norm());

   do {
      eta1 *= 2.0;
      // difference is negative for relatively small eta !
      difference = diffReciprocalReal (shortR, shortG, eta1);
   } while ( difference < 0. );

   do {
      eta2 /= 2.0;
      // difference is positive for relatively larg eta !
      difference = diffReciprocalReal (shortR, shortG, eta2);
   } while ( difference > 0. );

   etaOut = (eta1 + eta2)/2.0;

   do {
      if (diffReciprocalReal (shortR, shortG, etaOut) < 0.0) {
         eta2 = etaOut;
      }
      if (diffReciprocalReal (shortR, shortG, etaOut) > 0.0) {
         eta1 = etaOut;
      }
      etaOut = (eta1 + eta2)/2.0;
      nLoops += 1;
   } while (fabs(diffReciprocalReal (shortR, shortG, etaOut)) > 1e-5 && 
            nLoops < 30);

   if (nLoops == 30)  {
      cout << " WARNING: no optimal eta found !\n" ;
      cout << " eta1 = " << eta1 << " eta2 = " << eta2 << endl;
      cout << " eta = " << etaOut << endl;
   }
   
   return etaOut;
}

double SxTBHamSolver::diffReciprocalReal (double sR, double sG, double etaVal)
{
   double diffReal, diffReciprocal, diff;
   double sG2 = sG * sG;

   // --- both diifReal & diffReciprocal should be positive
   diffReal       = derfc(etaVal*2.*sR)/(2.*sR) 
                  - derfc(etaVal*3.*sR)/(3.*sR); 

   diffReciprocal = exp(-(2.*sG2)/(4.*etaVal*etaVal))/sG2
                  - exp(-(3.*sG2)/(4.*etaVal*etaVal))/sG2; 
      
   diff =  diffReciprocal - diffReal; 

   return diff;
}

void SxTBHamSolver::initRSphere()
{
   const Real8 rEwald2 = rEwald * rEwald;
   double latVec2;
   int n1, n2, n3;
   SxVector3<TReal8> latVec;

   int nC1 = (int)(rEwald/structure.cell(0).norm() + 0.5) + 2;
   int nC2 = (int)(rEwald/structure.cell(1).norm() + 0.5) + 2;
   int nC3 = (int)(rEwald/structure.cell(2).norm() + 0.5) + 2;

   // --- get R vectors that are inside cut-off sphere
   rLat.resize (0);
   for (n1 = -nC1; n1 <= nC1; n1++)  {
      for (n2 = -nC2; n2 <= nC2; n2++)  {
         for (n3 = -nC3; n3 <= nC3; n3++)  {
            latVec  = structure.cell ^ SxVector3<Int> (n1,n2,n3);
            latVec2 = latVec^latVec;
            if (latVec2  < rEwald2 ) // don't exclude R = 0  
               rLat.append (SxVector3<TReal8> (latVec));
         }
      }
   }
   cout << "|  number of R vectors   = " << rLat.getSize() << endl;
}

void SxTBHamSolver::initGSphere()
{
   const Real8 gEwald2 = gEwald * gEwald;
   double g2;
   int    n1, n2, n3;
   SxVector3<Double> gLatVec;
   
   const CellMat b (structure.cell.getReciprocalCell());

   int nC1 = (int)(gEwald/b(0).norm() + 0.5) + 2;
   int nC2 = (int)(gEwald/b(1).norm() + 0.5) + 2;
   int nC3 = (int)(gEwald/b(2).norm() + 0.5) + 2;

   // --- get G vectors that are inside cut-off sphere
   gLat.resize (0);
   for (n1 = -nC1; n1 <= nC1; n1++)  {
      for (n2 = -nC2; n2 <= nC2; n2++)  {
         for (n3 = -nC3; n3 <= nC3; n3++)  {
            gLatVec = b(0)*n1 + b(1)*n2 + b(2)*n3;
            g2      = gLatVec.absSqr().sum();
            if ( g2 < gEwald2 && g2 > 1e-10) // excluding G = 0 
               gLat.append (SxVector3<TReal8> (gLatVec));
         }
      }
   }
   cout << "|  number of G vectors   = " << gLat.getSize() << endl;
}
//----------------------------------------------------------------------------
// gamma derivative
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// a) for non-periodic systems (e.g. clusters)
//----------------------------------------------------------------------------

SxMatrix<Double> SxTBHamSolver::getDGammaMatNonPeriodic () const
{
   SxMatrix<Double>  dGammaMat(lOrbMax,lOrbMax);
   SxVector3<Double> distVec;
   int iA, jA, iS, jS, nOrb1, nOrb2, iOrb, jOrb;

   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS      = structure.getISpecies (jA);
         nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
         distVec = structure.getAtom(jA) - structure.getAtom(iA); 
         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 
               dGammaMat(idx(iA)+iOrb, idx(jA)+jOrb) 
                  = dGammaSCF (iA,jA,iOrb,jOrb,distVec.norm());
            }
         }
      }
   }

   return dGammaMat;
}

double SxTBHamSolver::dGammaSCF (int iA, int jA, int iOrb, int jOrb, 
                                 double distance) const 
{
   double dGamma, tauIA, tauJA;

   int iS = structure.getISpecies (iA);
   int jS = structure.getISpecies (jA);
   tauIA  = (16./5.) * tbAtomAtom(iS)(iS).getHubbardU(iOrb);
   tauJA  = (16./5.) * tbAtomAtom(jS)(jS).getHubbardU(jOrb);

   if (iA == jA)  
      // for non-periodic str. means r = 0., don't use it for periodic str.
      dGamma = 0.;
   else 
      dGamma = -1.0/(distance*distance) - dShortRange (tauIA,tauJA,distance);
  
   return dGamma;
}

double SxTBHamSolver::dShortRange (double tA, double tB, double distance) const
{
   double dShort;

   if (fabs(tA - tB) < 1e-10) // if same species, or tA = tB  somehow!
      //This function is the limit of the 2nd one (in the else part) as tB->tA
      dShort =( -48 - distance*tA*(4 + distance*tA)
               *(12 + distance*tA*(3 + distance*tA) )
              )/(48.*distance*distance*exp(distance*tA) ); 
   else 
      dShort = - tA*exp(-tA*distance)*K (tA,tB,distance)
               +    exp(-tA*distance)*dK(tA,tB,distance)
               - tB*exp(-tB*distance)*K (tB,tA,distance)
               +    exp(-tB*distance)*dK(tB,tA,distance);
              
   return dShort;
}

double SxTBHamSolver::dK (double tA, double tB, double distance) const 
{
   SX_CHECK (fabs(tA - tB) > 1e-8, tA, tB);

   double diffK;

   diffK = (-3*tA*tA*pow(tB,4) + pow(tB,6))
           /(distance*distance*pow(tA*tA - tB*tB,3)); 

   return diffK;
}

//----------------------------------------------------------------------------
// b) for periodic systems (e.g. solids)
//----------------------------------------------------------------------------

SxArray<SxMatrix<Double> > SxTBHamSolver::getDGammaMatPeriodic () const
{
   SxArray<SxMatrix<Double> > dGammaMat(3);
   SxMatrix<Double> matX(lOrbMax,lOrbMax);
   SxMatrix<Double> matY(lOrbMax,lOrbMax);
   SxMatrix<Double> matZ(lOrbMax,lOrbMax);

   SxVector3<Double> distVec;
   SxVector3<Double> dShortRangePart, dLongRangePart;
   int iA, jA, iS, jS, nOrb1, nOrb2, iOrb, jOrb;

   for (int i=0; i<3; i++)  {
      dGammaMat(i).reformat (lOrbMax,lOrbMax);
      dGammaMat(i).set(0.);
   }
   for (iA = 0; iA < nTlAtoms; iA++)  { 
      iS = structure.getISpecies (iA);
      for (jA = 0; jA < nTlAtoms; jA++)  { 
         jS      = structure.getISpecies (jA);
         nOrb1   = tbAtomAtom(iS)(jS).getNOrb1();
         nOrb2   = tbAtomAtom(iS)(jS).getNOrb2();
         distVec = structure.getAtom(jA) - structure.getAtom(iA); 

         dLongRangePart = dLongRangeSum (iA,jA,distVec); 

         for (iOrb = 0; iOrb < nOrb1; iOrb++)  {     
            for (jOrb = 0; jOrb < nOrb2; jOrb++)  { 

               dShortRangePart = dShortRangeSum (iA,jA,iOrb,jOrb,distVec);
               
               matX(idx(iA)+iOrb,idx(jA)+jOrb)
                  = dLongRangePart(0) - dShortRangePart(0);
               matY(idx(iA)+iOrb,idx(jA)+jOrb)
                  = dLongRangePart(1) - dShortRangePart(1);
               matZ(idx(iA)+iOrb,idx(jA)+jOrb)
                  = dLongRangePart(2) - dShortRangePart(2);
            }
         }
      }
   }
   dGammaMat(0) = SxMatrix<Double> (matX);
   dGammaMat(1) = SxMatrix<Double> (matY);
   dGammaMat(2) = SxMatrix<Double> (matZ);

   return dGammaMat;
}

SxVector3<Double> SxTBHamSolver::dShortRangeSum 
                                (int iA, int jA, int iOrb, int jOrb,
                                 SxVector3<Double> distVec) const 
{
   int iS = structure.getISpecies (iA);
   int jS = structure.getISpecies (jA);
   double tauIA  = (16./5.) * tbAtomAtom(iS)(iS).getHubbardU(iOrb);
   double tauJA  = (16./5.) * tbAtomAtom(jS)(jS).getHubbardU(jOrb);
   
   int n1, n2, n3, nR;
   int nMax = 50;
   int nMin = 3;
   double lastTerm = 0;
   double r;
   SxVector3<Double> latVec, rVec;
   SxVector3<Double> rSumVec; 

   rSumVec.set(0.);
   nR = 0;
   
   do  {
      lastTerm = 0;
      for (n1 = -nR; n1 <= nR; n1++)  {
         for (n2 = -nR; n2 <= nR; n2++)  {
            for (n3 = -nR; n3 <= nR; n3++)  {
               // only R's for the new shell contributes. Previous R's 
               // contributions were added in the previous do loop, and so on
               if (nR == abs(n1) || nR == abs(n2) || nR == abs(n3) )  {
                  
                  latVec  = structure.cell ^ SxVector3<Int> (n1,n2,n3);
                  
                  if (iA != jA) {     
                     rVec  = distVec - latVec;
                     r     = rVec.norm();  
                     lastTerm = dShortRange (tauIA,tauJA,r);
                     rSumVec += lastTerm * rVec/r;
                  }
                  // if (iA == jA) no derivative needed for the forces, 
                  // so, keep it zero.
               }
            }
         }
      }
      nR ++; 
   }  while (nR <= nMax && (nR <= nMin || fabs(lastTerm) > 1e-25) ); 
   
   return rSumVec;
}

SxVector3<Double> SxTBHamSolver::dLongRangeSum (int iA, int jA,
                                    SxVector3<Double> distVec) const 
{
   int iLat; 
   int nRLat = (int)rLat.getSize (); 
   int nGLat = (int)gLat.getSize ();
   SxVector3<Double> rSum, gSum;
   SxVector3<Double> dEwald; 
   SxVector3<Double> gVec, rVec;
   double r;
   double g2, arg;
   double factor = 4.*PI/structure.cell.volume;

   dEwald.set(0); rSum.set(0); gSum.set(0);
   if (iA != jA )  {
      // --- real space sum
      for (iLat = 0; iLat < nRLat; iLat++ )  { 
         rVec  = distVec - rLat(iLat);  
         r     = rVec.norm();  
         rSum += (-2./sqrt(PI)*exp(-eta*eta*r*r)*eta*r - derfc(eta*r))/(r*r*r)
               * rVec ;
      }
      // --- reciprocal space sum
      for (iLat = 0; iLat < nGLat; iLat++ )  { // G = 0 was excluded 
         gVec  = gLat(iLat);
         g2    = gVec.absSqr().sum();
         arg   = distVec ^ gVec;
         gSum += (exp(-g2/(4.*eta*eta))*sin(arg))/g2 * -gVec;
      }
      gSum *= factor;
   }
   // if (iA == jA) no derivative needed for the forces,
   // so, keep it zero
   
   // add contributions
   dEwald = rSum + gSum;
   
   return dEwald;
}


//----------------------------------------------------------------------------
// service routines
//----------------------------------------------------------------------------

SxMatrix<Double> SxTBHamSolver::extractMat (const SxMatrix<Double> &pMat,
                                            int index1, int index2, 
                                            int nOrb1,  int nOrb2) const
{
   SxMatrix<Double> subPMat;
   subPMat.reformat (nOrb1,nOrb2);
   int i, j;

   for (i = 0; i < nOrb1; i++)  {
      for (j = 0; j < nOrb2; j++)  {
         subPMat (i, j) = pMat(index1+i, index2+j);
      }
   }
   return subPMat;
}

SxMatrix<Complex16> SxTBHamSolver::extractMat (const SxMatrix<Complex16> &pMat,
                                               int index1, int index2, 
                                               int nOrb1,  int nOrb2) const
{
   SxMatrix<Complex16> subPMat;
   subPMat.reformat (nOrb1,nOrb2);
   int i, j;

   for (i = 0; i < nOrb1; i++)  {
      for (j = 0; j < nOrb2; j++)  {
         subPMat (i, j) = pMat(index1+i, index2+j);
      }
   }
   return subPMat;
}
SxVector<Double> SxTBHamSolver::putRhoInVec 
                             (SxArray<SxVector<Double> > &rhoArray) const
{
   SxVector<Double> rhoVec(lOrbMax);
   int iA, iS, iOrb, nOrb;
   
   for (iA = 0; iA < nTlAtoms; iA++)  {
      iS    = structure.getISpecies (iA);
      nOrb  = tbAtomAtom(iS)(0).getNOrb1();
      for (iOrb = 0; iOrb < nOrb; iOrb++)  { 
         rhoVec (idx(iA)+iOrb) = rhoArray(iA)(iOrb);
      }
   }

   return rhoVec;
}

SxArray<SxVector<Double> > SxTBHamSolver::putRhoInArrayOfVec 
                             (SxVector<Double> &rhoVec) const
{
   SxArray<SxVector<Double> > rhoArray(nTlAtoms);

   int iA, iS, iOrb, nOrb;
   
   for (iA = 0; iA < nTlAtoms; iA++)  {
      iS    = structure.getISpecies (iA);
      nOrb  = tbAtomAtom(iS)(0).getNOrb1();
      rhoArray(iA).resize(nOrb);
      for (iOrb = 0; iOrb < nOrb; iOrb++)  { 
         rhoArray(iA)(iOrb) = rhoVec (idx(iA)+iOrb) ;
      }
   }

   return rhoArray;
}

void SxTBHamSolver::printMullikenRho () const 
{ 
   int iA, iOrb;
   FILE *fp = NULL;

   fp = fopen ("mullikenRho.dat", "w");
   if (fp)  {
      fprintf (fp,"#Atom, total charge, charge per orbital\n");
      for (iA = 0; iA < nTlAtoms; iA++)  {
         fprintf (fp,"  %d   %12.9f : ", iA+1, mullikenRho(iA).sum()); 
         for (iOrb = 0; iOrb < mullikenRho(iA).getSize(); iOrb++)
            fprintf (fp,"%12.9f",mullikenRho(iA)(iOrb));
         fprintf (fp,"\n");
      }
   }
   fclose (fp);
}

bool SxTBHamSolver::isHermitian (const SxMatrix<Complex16> &matIn) const
{
   int iC, iR;
   bool isHerm = true;

   for (iR = 0 ; iR < matIn.nRows(); iR++)  {
      for (iC = 0 ; iC < matIn.nCols(); iC++)  {
        if (fabs(matIn(iR,iC).re - matIn(iC,iR).re) > 1e-8)  {
           isHerm = false;
           cout << " in real parts    = " << matIn(iR,iC) << " and " 
                                          << matIn(iC,iR) << endl;
        }
        if (fabs(matIn(iR,iC).im - matIn(iC,iR).conj().im) > 1e-8)  {
           isHerm = false;
           cout << " in complex parts = " << matIn(iR,iC) << " and " 
                                          << matIn(iC,iR) << endl;
        }

      }
   }
   if (!isHerm) { 
      cout << " matrix is not Hermiitian !\n";
      SX_EXIT;
   } 

   return isHerm;
}

