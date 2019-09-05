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

#include <SxInitPAW.h>
#include <SxHamSolver.h>
#include <SxFileParser.h>
#include <SxHubbardMO.h>

SxInitPAW::SxInitPAW (const SxSymbolTable *cmd,
                      const SxPtr<SxSpeciesData> &potPtrIn,
                      const SxPtr<SxGkBasis> &gkPtr,
                      const SxHamSolver *hamSolver)
   : potPtr (potPtrIn)
{
   SX_CHECK (hamSolver);
   const SxAtomicStructure &structure = hamSolver->structure;
   SxGBasis &G = const_cast<SxGBasis&>(hamSolver->G);
   const SxRBasis &R = hamSolver->R;
   int nSpin = hamSolver->nSpin;

   const SxSymbolTable *topLvl = cmd->topLevel ();
   SxSymbolTable *exchangeGroup =  cmd->containsGroup ("exchange")
                                ? cmd->getGroup("exchange")
                                : NULL;
   SxSymbolTable *wavesGroup =  cmd->containsGroup ("waves")
                                ? cmd->getGroup("waves")
                                : NULL;
   pawHam = SxPtr<SxPAWHamiltonian>::create ();
   pawHam->structure = structure;
   pawHam->potPtr    = SxPtr<SxPAWPot>(potPtr);
   pawHam->rPtr      = const_cast<SxRBasis*>(&R);

   pawHam->setupProj (*gkPtr);
   pawHam->setupPSRefG (G);
   pawHam->xcPtr = SxPtr<SxXC>::create (nSpin);
   try  {
      const SxSymbolTable *hamGroup = topLvl->getGroup ("PAWHamiltonian");
      pawHam->xcPtr->read (hamGroup, structure.cell);
      pawHam->read (hamGroup);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   if (pawHam->hubbardU)  pawHam->setupForHubbard (*gkPtr, nSpin);

   if (!pawHam->xcPtr->xcBasisPtr)  {
      unsigned oldMode = 0;
      SxFFT::quickFFTPlanner (SxFFT::Estimate, &oldMode);
      SxMesh3D xcMesh = 2 * R.fft3d.mesh;
      pawHam->xcPtr->xcBasisPtr=SxPtr<SxRBasis>::create(xcMesh,structure.cell);
      pawHam->xcPtr->xcBasisPtr->registerBasis (G);
      SxFFT::restorePlannerMode (oldMode);
      cout << "xc will be computed on " << xcMesh << " mesh.\n";
   }
   pawHam->xcPtr->origXcBasis = pawHam->xcPtr->xcBasisPtr;

   if (SxXC::isHybrid (pawHam->xcPtr->xcFunctional))  {
      SxPtr<SxPAWExchange> xPtr = xPtr.create ();
      pawHam->exchangePtr = xPtr;

      xPtr->potPtr = pawHam->potPtr;
      pawHam->potPtr->setupQijL ();
      pawHam->potPtr->computeCoreX (); 
      pawHam->potPtr->computeXKernel ();

      SxPtr<SxBasis> GX;
      // --- create new G basis for exchange
      // double gCutX = G.g2(G.ng-1) * 0.5;
      // GX = GX.create (mesh, structure, gCutX, Coord(0.,0.,0.), true,
      //                 SxGBasis::SaveTime);

      // --- use normal G basis for exchange (but precompute phase factors)
      GX = G.getThis ();
      G.memMode = SxGBasis::SaveTime;
      G.cleanPhaseFactors ();
      SX_LOOP(iR) G.fft3d(iR).createArrays (SxFFT::InArrayZero);
      G.changeTau (structure);

      // --- bases
      xPtr->pBasis = pawHam->pBasis;
      xPtr->rPtr   = R.getThis ();
      xPtr->gPtr   = GX;
      xPtr->pBasisX = pawHam->pBasis;

      // --- Ylm rotation group
      xPtr->ylmRot = pawHam->ylmRot;

      SxString xFile;
      if (exchangeGroup && exchangeGroup->contains ("file")) {
         xFile = exchangeGroup->get("file")->toString ();
      } else if (wavesGroup && wavesGroup->contains("file")) {
         xFile = wavesGroup->get("file")->toString ();
         cout << "WARNING: duplicating waves in PAW init!" << endl;
         // TODO: if waves are read from file, we should use
         // SxHamSolver's wavesPtr...
      } else  {
         cout << "Exchange waves has to be initialized from file!" << endl;
         SX_EXIT;
      }
      SxBinIO io (xFile, SxBinIO::BINARY_READ_ONLY);
      SxPtr<SxPAWSet> xWavesPtr = 
         SxPtr<SxPAWSet>::create(pawHam->potPtr, pawHam->structure, io);
      SxPtr<SxFermi> xFermiPtr = SxPtr<SxFermi>::create(io);
      xPtr->wavesPtr = xWavesPtr;
      xPtr->setFocc (xFermiPtr->focc);
      xPtr->pBasisX = xWavesPtr->getPBasis ();
      xPtr->computeXij (SxPAWHamiltonian::computeDij (*xWavesPtr, 
                        xFermiPtr->focc, pawHam->potPtr, 
                        *xPtr->pBasisX, pawHam->ylmRot ));
   }
   pawHam->pawRho.potPtr = potPtr;
   pawHam->pawRho.pwRho = SxRho(R, nSpin, hamSolver->nElectrons);
   pawHam->pawRho.Dij.resize (structure.atomInfo, nSpin);
}

void SxInitPAW::randomRho ()
{
   const SxRBasis &R = *pawHam->rPtr;
   const SxGBasis &G = R.getGBasis ();
   int nSpin = pawHam->pawRho.getNSpin ();
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   pawHam->pawRho.pwRho;
   pawHam->pawRho.pwRho.randomize ();
   pawHam->pawRho.pwRho.nElectrons = -1.;
   pawHam->pawRho.Dij.resize (structure.atomInfo, *pawHam->potPtr, nSpin);
   pawHam->pawRho.Dij.set (0.);
   PsiG coresG(G);
   coresG.set (0.);
   for (int is = 0; is < structure.getNSpecies (); ++is)
      coresG += G.structureFactors(is)
                * pawHam->rhoCorePSG(is);
   coresG /= nSpin;
   SxMeshR coresR = coresG.to (R);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      pawHam->pawRho.pwRho(iSpin) += coresR;
   }
}

static void spinMismatch (const char *where)
{
   cout << SX_SEPARATOR << "ERROR:" << endl;
   cout << "Initial guess has " << where
        << ", but PAWHamiltonian is not spin-polarized." << endl;
   cout << "You must set the spinPolarized flag in PAWHamiltonian!" << endl;
   SX_QUIT;
}

static void checkNoTotalSpin (const SxSymbolTable* rhoGroup)
{
   if (rhoGroup->contains("spinMoment")) {
      cout << SX_SEPARATOR << "ERROR:" << endl;
      cout << "You want to set up single spins on atoms and a total spin "
           << " for the whole system?" << endl;
      cout << "That does not work." << endl;
      SX_QUIT;
   }
}

void SxInitPAW::atomicRho (const SxSymbolTable *rhoGroup,
                           const SxAtomicStructure &structure)
{
   const SxRBasis &R = *pawHam->rPtr;
   int nSpin = pawHam->pawRho.getNSpin ();
   double nElectrons = pawHam->pawRho.pwRho.nElectrons;
   pawHam->pawRho.pwRho.nElectrons = -1.; // suppress normalization
   if (rhoGroup->containsGroup ("atomicSpin")) {
      if (nSpin != 2) spinMismatch ("atomicSpin");
      checkNoTotalSpin (rhoGroup);

      SxArray<SxArray<double> > atomSpin (structure.getNSpecies ());
      for (int is=0; is < structure.getNSpecies (); is++)  {
         atomSpin(is).resize (structure.getNAtoms (is));
         atomSpin(is).set (0.);
      }

      const SxSymbolTable *spinGroup;
      for (spinGroup = rhoGroup->getGroup("atomicSpin");
           spinGroup != NULL;
           spinGroup = spinGroup->nextSibling("atomicSpin"))
      {
         if (spinGroup->contains ("file"))  {
            SxFileParser fp(spinGroup->get("file")->toString ());
            SxVector<Double> spins = fp.getVector (structure.getNAtoms());
            SX_LOOP2(is,ia)  {
               atomSpin(is)(ia) = spins(structure.getIAtom ((int)is,(int)ia));
            }
         } else {
            const SxArray<SxString> &labels = structure.getLabels ();
            // get label and spin for this label
            const SxString spinLabel = spinGroup->get("label")->toString();
            double spin = spinGroup->get("spin")->toReal();
            cout << "The spin for the label " << spinLabel << " is " << spin << endl;
            // set atomSpin for atoms that have our label
            for (int is=0; is < structure.getNSpecies(); is++){
               for (int ia=0; ia < structure.getNAtoms(is); ia++) {
                  int iTlAtom = structure.getIAtom (is,ia);
                  if (spinLabel == labels(iTlAtom)) {
                     atomSpin(is)(ia) = spin;
                     cout << "atomSpin("<<is<<")("<<ia<<")= :" << atomSpin(is)(ia) << endl;
                  }
               }
            }
         }
      }
      pawHam->pawRho.atomicChargeDensity (structure, potPtr, R, atomSpin);
   } else if (!rhoGroup || !rhoGroup->containsGroup ("atom"))  {
      pawHam->pawRho.atomicChargeDensity (structure, potPtr, R, nSpin);
      if (rhoGroup && rhoGroup->contains("spinMoment")) {
         if (nSpin != 2) spinMismatch ("atoms with spin");
         double spinMoment = rhoGroup->get("spinMoment")->toReal();
         SxMeshR corePS = 0.5 * pawHam->computeCorePS ();
         SX_LOOP(iSpin) pawHam->pawRho.pwRho(iSpin) -= corePS;
         pawHam->pawRho.pwRho(0) *= 1. + spinMoment/nElectrons;
         pawHam->pawRho.pwRho(1) *= 1. - spinMoment/nElectrons;
         SX_LOOP(iSpin) pawHam->pawRho.pwRho(iSpin) += corePS;
         SX_LOOP2(is, ia)  {
            pawHam->pawRho.Dij(0,is,ia) *= 1. + spinMoment/nElectrons;
            pawHam->pawRho.Dij(1,is,ia) *= 1. - spinMoment/nElectrons;
         }
         cout << "The spin is " << pawHam->pawRho.getSpin() << endl;
      }
   } else {
      // --- use label-specific occupation numbers
      const SxArray<SxString> &labels = structure.getLabels ();
      SxArray2<SxVector<Double> > focc(nSpin, structure.getNAtoms ());
      // --- initialize focc with potential values
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         for (int is=0, iTlAtom=0; is < structure.getNSpecies (); ++is)
            for (int ia=0; ia < structure.getNAtoms(is); ++ia,++iTlAtom)
               focc(iSpin, iTlAtom) = pawHam->potPtr->foccInit(is)
                                    / double(nSpin);

      // --- run through atom types
      const SxSymbolTable *grp;
      for (grp = rhoGroup->getGroup ("atom"); grp != NULL;
           grp = grp->nextSibling ("atom"))
      {
         // --- read symbol table
         SxString label;
         SxVector<Double> atomFocc;
         double spin = 0.;
         try {
            label = grp->get ("label")->toString ();
            if (grp->contains ("focc"))  {
               atomFocc = grp->get ("focc")->toList ();
            }
            if (nSpin == 2 && grp->contains ("spin"))  {
               spin = grp->get ("spin")->toReal ();
cout << "WARNING: using atom { label=... spin=... } is deprecated." <<endl
  << "         Reason: atom {} sets up the density from potentially"
  << endl
  << "         truncated atomic wavefunctions." << endl
  << "         Use atomicSpin {} instead" << endl;
            }
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }

         // --- set up occupations
         for (int ia = 0; ia < structure.getNAtoms (); ++ia)  {
            if (labels(ia) != label) continue;
            int nf =  (int)focc(0, ia).getSize ();
            if (atomFocc.getSize () % nf != 0 
                || atomFocc.getSize () / nf > nSpin) {
               cout << "Label = " << label << ": ";
               cout << "atomic occupation has " << atomFocc.getSize () 
                    << " components, but potential has "
                    << focc(0, ia).getSize ()
                    << endl;
            }
            if (atomFocc.getSize () == 0 && nSpin == 2)  {
               // --- spin initialization
               SxVector<Double> oldFocc = focc(0, ia);
               double n = 2. * oldFocc.sum ();
               focc(0, ia) = (n + spin)/n * oldFocc;
               focc(1, ia) = (n - spin)/n * oldFocc;
            } else if (nSpin == 2)  {
               // --- focc initialization
               if (atomFocc.getSize () == nf) {
                  focc(0, ia) = atomFocc / 2.;
                  focc(1, ia) = atomFocc / 2.;
               } else {
                  focc(0, ia).copy (atomFocc(SxIdx (0, nf-1)));
                  focc(1, ia).copy (atomFocc(SxIdx (nf, 2 * nf - 1)));
               }
            } else {
               if (atomFocc.getSize () > 0) focc(0, ia) = atomFocc;
            }
         }
      }
      // --- calculate density
      pawHam->pawRho.atomicChargeDensity (structure, potPtr, R, focc);
   }
}

void SxInitPAW::readRho (const SxString &file)
{
   pawHam->pawRho.pwRho.nElectrons = -1.;
   pawHam->pawRho.readRho (file);
   const SxRBasis &R = *pawHam->rPtr;
   const SxGBasis &G = R.getGBasis (); 
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   int nSpin = pawHam->pawRho.getNSpin ();

   SxAtomicStructure tauFile;
   try {
      SxBinIO io(file, SxBinIO::BINARY_READ_ONLY);
      if (io.contains ("tau")) tauFile.read (io);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   if (tauFile.getNAtoms () > 0)  {
      if (tauFile.getNSpecies () == structure.getNSpecies ()  
          && tauFile.getNAtoms () == tauFile.getNAtoms ()
          && (  tauFile  .atomInfo->nAtoms 
              - structure.atomInfo->nAtoms).sqr ().sum () == 0)
      {
         tauFile.replaceInfo (structure.atomInfo);
         if ((tauFile - structure).absSqr ().sum ()
               > 1e-6 * structure.getNAtoms ())
         {
            // --- structure has changed: move atomic density
            cout << SX_SEPARATOR;
            cout << "| WARNING: initial density from displaced atoms."
                 << endl << SX_SEPARATOR;
            SxVector<Double> coordsNow = structure.coordRef ().getCopy ();
            const_cast<SxGBasis&>(G).changeTau (tauFile);
            SxPAWRho atomRho;
            atomRho.atomicChargeDensity (tauFile, potPtr, R, nSpin);
            pawHam->pawRho -= atomRho;

            structure.coordRef () <<= coordsNow;

            const_cast<SxGBasis&>(G).changeTau (structure);
            atomRho.atomicChargeDensity (structure, potPtr, R, nSpin);
            pawHam->pawRho += atomRho;
         }
      } else {
         cout << "WARNING: initial density from different structure"
              << endl;
      }
   }
}

void SxInitPAW::addCharge (const SxMeshR &chargeR)
{
   int nSpin = pawHam->pawRho.getNSpin ();
   for (int iSpin = 0.; iSpin < nSpin; ++iSpin)
      pawHam->pawRho.pwRho(iSpin).plus_assign_ax(-1./nSpin, chargeR);
}

SxPtr<SxPWSet> SxInitPAW::setupWaves (const SxSymbolTable *wavesGroup,
                                      const SxPtr<SxGkBasis> &gkPtr,
                                      int nStates, int nSpin)
{
   SxString tmpDir = "";
   bool pawBasis = true;
   if (wavesGroup)  {
      if (wavesGroup->contains ("pawBasis"))
         pawBasis = wavesGroup->get("pawBasis")->toAttribute ();
      if (wavesGroup->contains("keepWavesOnDisk")
          && wavesGroup->get("keepWavesOnDisk")->toAttribute ())
         tmpDir = ".";
   }
      
   if (pawBasis)
      return SxPtr<SxPAWSet>::create (gkPtr, pawHam->pBasis, 
                                      nStates, nSpin, tmpDir);
   return SxPtr<SxPW>::create (nStates, nSpin, gkPtr, tmpDir);
}

void SxInitPAW::randomize (const SxPtr<SxPWSet> &wavesPtr)
{
   const SxGkBasis &gk = wavesPtr->getGkBasis ();
   SxOverlap S = pawHam->getS ();
   SX_NEW_LOOP (*wavesPtr);
   for (int ik = 0; ik < wavesPtr->getGkBasis ().getNk (); ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork (ik)) {
         for (int iSpin = 0; iSpin < wavesPtr->getNSpin (); ++iSpin)  {
            // randomize plane-waves & projections
            (*wavesPtr)(iSpin,ik).randomize ();
            if (dynamic_cast<const SxPAWBasis*>(&wavesPtr->getBasis (ik)))  {
               // compute projections consistent with plane-waves
               (*wavesPtr)(iSpin, ik) 
                  <<= (wavesPtr->getBasis (ik) | (gk(ik) | (*wavesPtr)(iSpin, ik)));
            }
            S.orthonormalize (&(*wavesPtr)(iSpin,ik), GramSchmidt);
         }
      }
   }
}

void SxInitPAW::rhoFromWaves (const SxPtr<SxPWSet> &wavesPtr,
                              const Focc &focc)
{
   pawHam->pawRho.pwRho.nElectrons = -1.;
   pawHam->computeRho (focc, *wavesPtr);
}

SxPtr<SxHamiltonian> SxInitPAW::setupHam (const SxPtr<SxPWSet> &,
                                          const SxPtr<SxSpeciesData> &,
                                          const SxAtomicStructure &)
{
   return pawHam;
}
