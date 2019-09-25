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


#include <SxPAWHamiltonian.h>
#include <SxRadBasis.h>
#include <SxProjector.h>
#include <SxNeighbors.h>
#include <SxYlm.h>
#include <SxSphereGrid.h>
#include <SxGaussIJ.h>
#include <SxPAWOverlap.h>
#include <SxForceSym.h>
#include <SxRadialAtom.h>
#include <SxPartialWaveBasis.h>
#include <SxPAWBasis.h>
#include <SxLoopMPI.h>
#include <SxTaskGroup.h>
#include <SxParallelHierarchy.h>
#include <SxPAWSet.h>
#include <SxHubbardMO.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif


double SxPAWHamiltonian::rSoft = 1.0;

/*
References:
 1. P. E. Bloechl, Phys. Rev. B 50, 17953 (1994).
 2. C. Freysoldt, PAW implementation notes
 3. C. Freysoldt, Hubbard U implementation notes
*/

SxPAWHamiltonian::SxPAWHamiltonian ()
   : SxHamiltonian (),
     hContrib(CalcAll),
     eTotal (0.),
     rPtr (NULL)
{
   // empty
}


SxPAWHamiltonian::SxPAWHamiltonian (const SxGBasis  &G,
                      const SxPtr<SxPAWPot> &potPtrIn,
                      const SxGkBasis &Gk,
                      const SxString  &rhoFile,
                      const SxSymbolTable *table)
/*
                      //const SxFermi &fermi,
                      //const SxPtr<SxPWSet> &wavesPtr)
*/
   : hContrib(CalcAll),
     eTotal (0.),
     potPtr(potPtrIn),
     rPtr (NULL),
     pawRho (potPtrIn)
{
   SX_CHECK (G.structPtr);
   structure = *G.structPtr;
   const SxRBasis &R = G.getRBasis ();
   int nSpin = SxHamiltonian::getNSpin (table);
   pBasis    = SxPtr<SxPartialWaveBasis>::create (potPtr, structure);
   rPtr      = const_cast<SxRBasis*>(&R);
   setupProj (Gk);
   setupPSRefG (G);
   read (table);
   if (hubbardU) setupForHubbard (Gk, nSpin);
   xcPtr = SxPtr<SxXC>::create (nSpin);
   try  {
      xcPtr->read (table->getGroup ("PAWHamiltonian"), structure.cell);
   } catch (const SxException &e)  {
      e.print ();
      SX_EXIT;
   }
   if (!xcPtr->xcBasisPtr)  {
      unsigned oldMode = 0;
      SxFFT::quickFFTPlanner (SxFFT::Estimate, &oldMode);
      SxMesh3D mesh2 = 2 * R.fft3d.mesh;
      xcPtr->xcBasisPtr = SxPtr<SxRBasis>::create (mesh2, structure.cell);
      xcPtr->xcBasisPtr->registerBasis (G);
      SxFFT::restorePlannerMode (oldMode);
      cout << "xc will be computed on " << mesh2 << " mesh.\n";
   }
   xcPtr->origXcBasis = xcPtr->xcBasisPtr;

   if (SxXC::isHybrid (xcPtr->xcFunctional))  {
      // Look at initialGuess and think about how it is done correctly
      SX_EXIT;
   }
   pawRho.pwRho = SxRho(R, nSpin, -1.);
   pawRho.readRho (rhoFile);
   /*
   if (pawHamPtr->pawRho.Dij(0).getNSpecies () != structure.getNSpecies ())  {
      cout << SX_SEPARATOR;
      cout << "| WARNING " << endl;
      cout << "| Species in rho and in structure differ, rho is now calculated from waves!" << endl;
      cout << SX_SEPARATOR;
      pawHamPtr->computeRho (fermi.focc, *wavesPtr);
   }
   */

   computeRhoTerms ();

}

SxPAWHamiltonian::~SxPAWHamiltonian ()
{
   // empty
}

void SxPAWHamiltonian::compute (bool, bool)
{
   SX_EXIT; // not yet implemented
}

void SxPAWHamiltonian::printEnergies () const
{
   //SX_EXIT;
}
void SxPAWHamiltonian::read (const SxSymbolTable *table)
{
   const SxSymbolTable *hamGroup = getHamiltonianGroup (table);
   if (hamGroup->contains ("rSoft"))  {
      rSoft = hamGroup->get ("rSoft")->toReal ();
   }
   // --- dipole correction
   if (withDipoleCorrection (hamGroup))  {
      dipoleCorr = SxPtr<SxDipoleCorrZ>::create ();
      try {
         if (hamGroup->contains("fixedDipoleZ"))  {
            double dipolePos = hamGroup->get("fixedDipoleZ")->toReal();
            dipoleCorr->fixedDipoleZ = true;
            dipoleCorr->dipolePos = dipolePos;
         }
         if (hamGroup->contains("zField"))
            dipoleCorr->extraField = hamGroup->get("zField")->toReal() / HA2EV;
      } catch (const SxException &e) {
         e.print ();
         SX_EXIT;
      }
   }

   // --- vdW correction
   if (hamGroup->containsGroup ("vdwCorrection")) {
      applyVDWCorrection = true;
      vdwCorrection = SxVDW(structure, hamGroup);
      SxArray<double> freeVolume(structure.nTlAtoms);
      for (int i = 0; i<structure.nTlAtoms; i++) {
         freeVolume(i) = 1.;
      }
   }

   if (hamGroup->containsGroup ("HubbardU"))  {
      hubbardU = SxPtr<SxHubbardU>::create ();
      hubbardU->read (hamGroup->getGroup ("HubbardU"), potPtr, structure);
      if (hubbardU->getSize () == 0) {
         cout << "WARNING: HubbardU produces no sites" << endl;
         hubbardU = SxPtr<SxHubbardU> ();
      }
   }

   // --- external potential
   if (hamGroup->containsGroup("vExt"))  {
      SxString vExtFile = hamGroup->getGroup("vExt")->get("file")->toString();

      try  {
         SxBinIO io (vExtFile, SxBinIO::BINARY_READ_ONLY);
         SX_CHECK (rPtr);
         const SxRBasis &R = *rPtr;
         int nSpinFile = io.getDimension ("nMeshes");
         if (nSpinFile != 1)  {
            cout << "Inconsistency when reading external potential '"
                 << vExtFile << "'.\n";
            cout << "There should be one channel, but found " 
                 << nSpinFile << " in file." << endl;
            SX_QUIT;
         }
         SxMatrix3<Double> cellFile;
         SxVector3<Int>    dim;
         vExt = io.readMesh (&cellFile, &dim)(0);
         if ((cellFile - R.cell).absSqr ().sum () > 1e-7)  {
            cout << "cell inconsistency in external potential." << endl;
            cout << "Expected " << R.cell << endl;
            cout << "Found " << cellFile << endl;
            SX_QUIT;
         }
         if (!(dim == R.fft3d.mesh))  {
            cout << "Mesh inconsistency in external potential." << endl;
            cout << "Expected " << R.fft3d.mesh << ", found " << dim << endl;
            SX_QUIT;
         }
         vExt.setBasis (&R);
      } catch (const SxException &e) {
         e.print ();
         SX_EXIT;
      }

      //Check if vExt is symmetry consistent
      cout << SX_SEPARATOR;
      cout << "Check for symmetry consistent vExt ... ";
      SxDiracVec<Double> vExtSym = rPtr->symmetrize(vExt);
      if((vExtSym - vExt).normSqr() > 1e-6)   {
         cout << endl;
         cout << SX_SEPARATOR;
         cout << "vExt seems to break symmetry!\nS/PHI/nX quits here!" << endl;
         cout << SX_SEPARATOR;
         SX_QUIT;
      }
      cout << "ok" << endl;
      cout << SX_SEPARATOR;
   }
   if (hamGroup->contains ("hContrib"))  {
      hContrib = lround(hamGroup->get ("hContrib")->toReal ());
      cout << SX_SEPARATOR;
      cout << "| hContrib:" << endl;
      if (hContrib & CalcNuc)
         cout << "|    nuclear potential" << endl;
      if (hContrib & CalcBar)
         cout << "|    local pseudopotential vBar" << endl;
      if (hContrib & CalcXcRad)
         cout << "|    radial xc" << endl;
      if (hContrib & CalcXcPS)
         cout << "|    pseudo xc" << endl;
      if (hContrib & CalcKinRad)
         cout << "|    kinetic energy on-site" << endl;
      if (hContrib & CalcKinPS)
         cout << "|    kinetic energy" << endl;
      if (hContrib & CalcNHatRad)
         cout << "|    compensation charge (radial)" << endl;
      if (hContrib & CalcNHatGHard)
         cout << "|    compensation charge hard" << endl;
      if (hContrib & CalcNHatGSoft)
         cout << "|    compensation charge soft" << endl;
      if (hContrib & CalcCore)
         cout << "|    core energy" << endl;
      if (hContrib & CalcHartreeRad)
         cout << "|    Hartree (radial)" << endl;
      if (hContrib & CalcHartreePS)
         cout << "|    Hartree (pseudo)" << endl;
      if (hContrib & CalcU)
         cout << "|    \"U\" (compensation hard/soft)" << endl;
      if (hContrib & CalcCorValExch)
         cout << "|    core-valence exchange (if any)" << endl;
      if (hContrib & CalcValExch)
         cout << "|    valence exchange (if any)" << endl;
      if ((hContrib & CalcHubbardU) && hubbardU)
         cout << "|    Hubbard U" << endl;
      cout << SX_SEPARATOR;
   }
}

bool SxPAWHamiltonian::rereadTable (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SX_CHECK (xcPtr);
   bool changed;
   // xc mesh
   changed = xcPtr->rereadTable (table);

   // --- dipole correction
   if (table->contains ("dipoleCorrection"))  {
      if (!dipoleCorr)  {
         cout << "WARNING: no dipole correction in PAWHamiltonian"
              << endl;
         cout << "WARNING: From here on, there will be a dipole correction."
              << endl;
         dipoleCorr = SxPtr<SxDipoleCorrZ>::create ();
         hContrib &= ~CalcDipoleC; // until now, it was suppressed
      }
      changed = true;
      bool dipole = table->get ("dipoleCorrection")->toAttribute ();
      cout << "| Dipole correction is now switched "
           << (dipole ? "on" : "off") << "." << endl;

      bool recompute = (dipole ? 0 : 1) != ((hContrib & CalcDipoleC) ? 0 : 1);
      if (dipole)
         hContrib |= CalcDipoleC; // switch on
      else
         hContrib &= ~CalcDipoleC; // switch off

      // is dipole correction being switched?
      if (recompute) {
         cout << "| Recomputing potentials ..." << endl;
         computeRhoTerms ();
      }
   }
   return changed;
}

void SxPAWHamiltonian::backToDefault (const SxSymbolTable *)
{
   SX_CHECK (xcPtr);
   // xc basis
   xcPtr->rereadTable (NULL);

   // --- dipole correction
   if (dipoleCorr && !(hContrib & CalcDipoleC))  {
      cout << "| Dipole correction is now switched on." << endl;
      cout << "| Recomputing potentials ..." << endl;
      hContrib |= CalcDipoleC;
      computeRhoTerms ();
   }
   hContrib |= CalcDipoleC;

}

const SxSymbolTable *
SxPAWHamiltonian::getHamiltonianGroup (const SxSymbolTable *table)
{
   SxSymbolTable *hGroup = NULL;
   if (table->name == "PAWHamiltonian") return table;
   try  {
      hGroup = table->getGroup("PAWHamiltonian");
   }  catch (const SxException&) {  }
   return hGroup;
}

void SxPAWHamiltonian::setupForHubbard (const SxGkBasis &gk, int nSpin)
{
   SX_CHECK(hubbardU);
   SxArray<SxPtr<SxHubbardMO> > &moSite = hubbardU->moSite;
   pawRho.blockAO.resize (moSite.getSize ());
   for (int iMoSite = 0; iMoSite < moSite.getSize (); iMoSite++)  {
      // setup projectors
      moSite(iMoSite)->setupProjGk (gk);
      // cleanup setup-specific stuff
      moSite(iMoSite)->finalize ();
      // --- create entry in paw density, ...
      pawRho.blockAO(iMoSite) = SxPtr<SxBlockDensityMatrix>::create ();
      // ... resize density and "potential", ...
      SxBlockDensityMatrix &blockRho = *pawRho.blockAO(iMoSite),
                           &blockHam = moSite(iMoSite)->hamProj;
      int nSite = moSite(iMoSite)->getNSite ();
      int nAO = moSite(iMoSite)->getNAOperSite ();
      blockRho.resize (nSpin, nSite);
      blockHam.resize (nSpin, nSite);
      SX_LOOP2(iSpin,iSite)  {
         blockRho(iSpin,iSite).reformat (nAO, nAO);
         blockHam(iSpin,iSite).reformat (nAO, nAO);
         // ... and set to zero
         blockRho(iSpin,iSite).set (0.);
         blockHam(iSpin,iSite).set (0.);
      }
   }
}

void SxPAWHamiltonian::setupProj (const SxGkBasis &gk)
{
   SX_CHECK (gk.getNk() > 0);
   SX_CHECK (potPtr);

   SxArray<SxDiracMat<Double> > projRad(potPtr->getNSpecies ());
   for (int is = 0; is < potPtr->getNSpecies (); ++is)  {
      if (    potPtr->pPsFine.getSize () > 0
           && potPtr->pPsFine(is).getSize () > 0)
      {
         projRad(is) = potPtr->pPsFine(is);
      } else {
         projRad(is) = potPtr->pPS(is);
      }
   }
   projBasis = SxPtr<SxAOBasis>::create (gk, projRad, potPtr->lPhi);

   // --- setup partial wave basis
   pBasis = SxPtr<SxPartialWaveBasis>::create (potPtr, structure);
   pBasis->projectors = projBasis;

   // --- Setup symmetry group
   int lmax = potPtr->getLMax ();
   ylmRot = SxYlm::computeYlmRotMatrices 
      (structure.cell.symGroupPtr->getSymmorphic (), lmax);
}

SxRadMat SxPAWHamiltonian::computeDij (const SxPWSet &waves,
                                                const Focc &focc) const
{
   return computeDij (waves, focc, potPtr, *pBasis, ylmRot);
}

SxRadMat
SxPAWHamiltonian::computeDij (const SxPWSet &waves,
                              const Focc &focc,
                              const SxConstPtr<SxPAWPot> &potPtr,
                              const SxPartialWaveBasis &pBasis,
                              const SxYlmRotGroup &ylmRot)
{
   SX_CLOCK (Timer::Dij);

   const SxVector<Double> &weights = waves.getGkBasis ().weights;
   const SxAtomicStructure &structure = waves.getGkBasis ().getTau ();
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;
   int nSpin = waves.getNSpin ();
   SxRadMat res;
   res.resize (structure.atomInfo, pot, nSpin);
   res.set (0.);

   SX_NEW_LOOP (waves);
   for (int ik = 0; ik < waves.getNk (); ++ik) {
      for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
         SX_MPI_LEVEL("waves-k");
         if (SxLoopMPI::myWork(ik))
         {  // LoopMPI
            double weight = weights(ik);
            SxDiracVec<Complex16> p = pBasis | waves(iSpin,ik);
            int offset = 0;
            for (int is = 0; is < structure.getNSpecies (); ++is)  {
               int nProjLocal = pot.getNProj (is);

               // --- set up i^l
               SxVector<Complex16> il(nProjLocal); // i^(l)  (:ip)
               for (int ipt = 0, ipl = 0; ipt < pot.lPhi(is).getSize (); ++ipt)  {
                  int l = pot.lPhi(is)(ipt);
                  SxComplex16 c = (l & 1) ? I : SxComplex16(1.);
                  if (l & 2) c=-c;
                  for (int m = -l; m <=l; ++m, ++ipl)  {
                     il(ipl) = c;
                  }
               }

               for (int ia = 0; ia < structure.getNAtoms(is); ++ia, offset+=nProjLocal)  {
                  SxMatrix<Double> &D = res(iSpin, is, ia);
                  for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
                     for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
                        SxComplex16 phase = il(ipl) * il(jpl).conj ();
                        for (int iState = 0; iState < p.nCols (); ++iState)  {
                           // ref. 1, Eq. (52)
                           // + projector i^l phase-correction (Ref. 2)
                           // + use k-dependent weights
                           // + exploit k/-k symmetry (=> .re)
                           D(ipl, jpl) += weight * focc(iState,iSpin,ik)
                                           * (  p(offset+ipl,iState).conj ()
                                                 * p(offset+jpl,iState) * phase).re;
                        }
                     }
                  }
               }
            }
         }   // LoopMPI
      }
   }

   // --- symmetrize
   SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
   SX_MPI_TARGET (TopLevel, TaskGroupAll);
   res.sumMPI (); // LoopMPI
   res.symmetrize (structure, pot, ylmRot);

   return res;
}

enum LocTimer { Setup, TimeYlmGl, Loop1, Loop2, Loop3, PsStuff, NHatG,
                PsStuff1, PsStuff2, PsStuff3, PsStuff4, PsStuff5, PsStuff6 };
SX_REGISTER_TIMERS (LocTimer)
{
   regTimer (Setup, "setup");
   regTimer (TimeYlmGl, "YlmGl");
   regTimer (Loop1, "loop 1");
   regTimer (NHatG, "nHatG");
   regTimer (Loop2, "loop 2");
   regTimer (Loop3, "loop 3");
   regTimer (PsStuff, "ps stuff");
   regTimer (PsStuff1, "ps stuff 1");
   regTimer (PsStuff2, "ps stuff 2");
   regTimer (PsStuff3, "ps stuff 3");
   regTimer (PsStuff4, "ps stuff 4");
   regTimer (PsStuff5, "ps stuff 5");
   regTimer (PsStuff6, "ps stuff 6");
}


void SxPAWHamiltonian::computeRhoTerms ()
{
   SX_CLOCK (Timer::ComputeHamRho);
   const SxRadMat &Dij = pawRho.Dij;
   const RhoR &rhoPSR = pawRho.pwRho.rhoR;
   SX_START_TIMER (Setup);
   const double SQRT_4PI = sqrt(FOUR_PI);
   double rNorm = FOUR_PI / sqrt(structure.cell.volume);
   SX_CHECK_NUM (rNorm);
   SX_ALLOC_CACHE;

   SX_CHECK (potPtr);
   SX_CHECK (xcPtr);
   SX_CHECK (rhoPSR.getSize () > 0, rhoPSR.getSize ());
   const SxPAWPot &pawpot = *potPtr;

   // --- get Dirac bases
   const SxRadBasis &rad  = pawpot.getRadBasis ();
   const SxRBasis &rBasis = rhoPSR(0).getBasis<SxRBasis> ();
   const SxGBasis &gBasis = rBasis.getGBasis ();
   // indices without G=0 component
   SxIdx noZero(1,gBasis.ng - 1);

   int nSpin = Dij.getNSpin ();
   SX_CHECK  (nSpin > 0 && nSpin <= 2, nSpin);
   SX_CHECK (rhoPSR.getSize () == nSpin, rhoPSR.getSize (), nSpin);

   // compute total rhoPS in G-space
   PsiG rhoPSG = gBasis | ((nSpin == 1) ? rhoPSR(0)
                                        : rhoPSR(0) + rhoPSR(1));

   double eXcRadial      = 0.;
   double eHartreeRadial = 0.;
   double eBarRadial     = 0.;
   double eR0            = 0.;
   double vHatGZero = 0.;
   eDoubleCounting = 0.;

   if (hubbardU) hubbardU->setZero ();

   //double eVnc = 0.;
   Vij.resize (Dij);
   Vij.set (0.);

   // Ref 1, Eq. 35
   SxArray<SxArray<SxVector<Double> > > dEdQrl(structure.getNSpecies ());
   Qrl.resize (structure.getNSpecies ());

   // \hat n, \hat n' in G-space
   PsiG nHatG(gBasis), nHatpG(gBasis);
   nHatG.set (0.);     nHatpG.set (0.);

   // |G|^l * Ylm(G)
   SX_START_TIMER (TimeYlmGl);
   int lMaxRho = pawpot.lMaxRho.maxval ();
   SxDiracMat<Double> YlmGl(gBasis.ng, sqr(lMaxRho + 1));
   YlmGl.setBasis (gBasis);
   for (int lm = 0, l = 0; l <= lMaxRho; ++l)  {
      SxDiracVec<Double> gl = pow (gBasis.g2, 0.5 * l);
      for (int m = -l; m <= l; ++m,++lm)
         YlmGl.colRef(lm) <<= gBasis.getYlm (l,m) * gl;
   }
   SX_STOP_TIMER (TimeYlmGl);

   // Hartree potential from rhoPSG (G != 0)
   PsiG vRhoPSG = FOUR_PI / gBasis.g2(noZero) * rhoPSG(noZero);
   if (!(hContrib & CalcHartreePS)) vRhoPSG.set (0.);

   SX_STOP_TIMER (Setup);

   SxVector<Double> MSpin(structure.getNAtoms()) ;
   MSpin.set(0.);
   SX_CHECK(potPtr);
   const SxPAWPot &pawPot = *potPtr;

   // --- loop over atoms
   SX_START_TIMER (Loop1);
   for (int is = 0; is < structure.getNSpecies (); ++is)
   {

      SxRadialMesh v1AE,v1PS;
      SxArray<SxRadialMesh> vXcAE,vXcPS;
      SxArray<SxDiracVec<Double> > gRl;

      double logDr = rad.logDr(is);
      int lMax = pawpot.lMaxRho (is);
      double r0cub = pow(rad.radFunc(is)(0), 3.);

      // --- resize objects
      Qrl(is).resize (structure.getNAtoms (is));
      dEdQrl(is).resize (structure.getNAtoms (is));
      v1AE.resize (rad, is, lMax);
      v1PS.resize (rad, is, lMax);

      // --- set up generalized Gaussians
      gRl.resize (lMax+1);
      for (int l = 0; l <= lMax; ++l)  {
         gRl(l) = pawpot.getGrl (is, l);
      }
      // --- Gaussians in G space
      double rc2 = sqr(pawpot.rc(is)), rcp2 = sqr(rSoft);
      SxDiracVec<Double> gaussHard = exp ( -(0.25 * rc2 ) * gBasis.g2),
            gaussSoft = exp ( -(0.25 * rcp2) * gBasis.g2);

      SxRadialAtom radial(xcPtr->xcFunctional, pawpot.aGridType(is));

      // ---
      /*
      SxDiracVec<Double> VncAE
         = -(SQRT_4PI * pawpot.nuclearCharge(is))/rad.radFunc(is)
         + radial.getHartreePotential (pawpot.rhoCoreAE(is));
      SxDiracVec<Double> VncPS;
      {
         SxDiracVec<Double> dRhoCore = pawpot.rhoCorePS(is).getCopy ();
         dRhoCore -= gRl(0) * (tr (pawpot.rhoCorePS(is)) + SQRT_1_4PI * pawpot.valenceCharge(is));
         VncPS = radial.getHartreePotential (dRhoCore);
      }
      */
      
      // --- setup omegaij(ipl,jpl) for spin constraints
      SxMatrix<Double> omegaij(pawPot.getNProj(is),pawPot.getNProj(is));
      omegaij.set (0.);
      SX_LOOP2(ipt, jpt)  {
         int li = pawPot.lPhi(is)(ipt);
         int lj = pawPot.lPhi(is)(jpt);
         if (li == lj)  {
            int iOff = pawPot.offset(is)(ipt),
                jOff = pawPot.offset(is)(jpt);
            for (int im = 0; im < 2 * li + 1; ++im)
               omegaij(im+iOff,im+jOff) = pawPot.omegaPAW(is)(ipt,jpt);
         }
      }

      for (int ia = 0; ia < structure.getNAtoms (is); ++ia)
      {
         int iTlAtom = structure.getIAtom (is, ia);
         if (SxLoopMPI::myWork(iTlAtom)) // SxParallelHierarchy
         {

            // Compute 1-center density
            SxArray<SxRadialMesh> rhoRadPS(nSpin), rhoRadAE(nSpin);
            for (int iSpin = 0; iSpin <  nSpin; ++iSpin)  {
               rhoRadPS(iSpin) = pawpot.computeRhoPS (Dij(iSpin,is,ia),is,nSpin);
               rhoRadAE(iSpin) = pawpot.computeRhoAE (Dij(iSpin,is,ia),is,nSpin);
            }

            // --- xc part
            if (hContrib & CalcXcRad)  {
               if (   pawpot.kjXcShape.getSize () > 0
                     && pawpot.kjXcShape(is).getSize () > 0)  {
                  // --- Kresse-Joubert variant: include compensation charge
                  // in pseudo-xc
                  SxArray<SxRadialMesh> rhoXcKj(nSpin);
                  double Q = 0.;
                  for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                     Q += ((rhoRadAE(iSpin)(0,0) - rhoRadPS(iSpin)(0,0))
                           * rad.radFunc(is).cub ()).integrate (logDr);
                  }
                  // remove core charges (not included in Kresse-Joubert xc)
                  Q -= ((pawpot.rhoCoreAE(is) - pawpot.rhoCorePS(is))
                        * rad.radFunc(is).cub ()).integrate (logDr);
                  // cout << "Q=" << Q << endl;
                  for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                     rhoXcKj(iSpin).resize (rhoRadPS(iSpin));
                     rhoXcKj(iSpin).meshData.copy (rhoRadPS(iSpin).meshData);
                     rhoXcKj(iSpin)(0,0) += Q/nSpin * pawpot.kjXcShape(is);
                  }

                  vXcPS = radial.computeXC (rhoXcKj);
                  eXcRadial -= radial.eXc;
               } else {
                  vXcPS = radial.computeXC (rhoRadPS);
                  eXcRadial -= radial.eXc;
               }
               // --- double counting correction PS xc
               eDoubleCounting -= radial.eXc;
               for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                  for (int l = 0; l <= lMax; ++l)  {
                     for (int m = -l; m <= l; ++m)  {
                        eDoubleCounting +=  tr(  rhoRadPS(iSpin)(l,m)
                                               * vXcPS(iSpin)(l,m));
                     }
                  }
                  // rhoRadPS contains pseudo-core electrons, which are not 
                  // double-counted
                  if (pawpot.rhoCorePS(is).getSize () > 0)
                     eDoubleCounting -=  tr(pawpot.rhoCorePS(is)
                                            * vXcPS(iSpin)(0,0)) / nSpin;
               }
               // --- all electron xc
               vXcAE = radial.computeXC (rhoRadAE);
               eXcRadial += radial.eXc;

               // --- double counting correction AE xc
               eDoubleCounting += radial.eXc;
               for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                  for (int l = 0; l <= lMax; ++l)  {
                    for (int m = -l; m <= l; ++m)  {
                      eDoubleCounting -=  tr(rhoRadAE(iSpin)(l,m)
                                             * vXcAE(iSpin)(l,m));
                    }
                  }
                  // rhoRadAE contains core electrons, which are not 
                  // double-counted
                  if (pawpot.rhoCoreAE(is).getSize () > 0)
                     eDoubleCounting +=  tr(pawpot.rhoCoreAE(is)
                                            * vXcAE(iSpin)(0,0)) / nSpin;
               }
            }

            if (nSpin == 2)  {
               // switch to total densities, save in iSpin=0
               rhoRadPS(0) += rhoRadPS(1);
               rhoRadPS(1) = SxRadialMesh ();
               rhoRadAE(0) += rhoRadAE(1);
               rhoRadAE(1) = SxRadialMesh ();
            }
            // --- vBar part (Ref. 1, Eq. 21) without ps core
            if (hContrib & CalcBar)  {
               bool hasCore = pawpot.rhoCorePS(is).getSize () > 0;
               SxDiracVec<Double> rhoVal 
                  = hasCore ? rhoRadPS(0)(0,0) - pawpot.rhoCorePS(is)
                            : rhoRadPS(0)(0,0);
               eBarRadial += tr (rhoVal * pawpot.vBar(is));
               eR0 += rhoVal(0) * pawpot.vBar(is)(0) * r0cub / 3.;
               //cout << "eBarR0=" << rhoVal(0) * pawpot.vBar(is)(0) * r0cub / 3. << endl;
            }

            // --- compute QRL
            SxRadialMesh rhoPAWcorr = rhoRadAE(0) - rhoRadPS(0);
            SxVector<Double> &Q = Qrl(is)(ia);
            Q.resize(rhoPAWcorr.meshData.nCols ());
            for (int l = 0; l <= rhoPAWcorr.lmax; ++l)  {
               SxDiracVec<Double> rl = pow(rad.radFunc(is),double(l+3));
               for (int m = -l; m <= l; ++m)  {
                  int lm = SxYlm::combineLm(l,m);
                  // Ref. 1, Eq. (24)
                  Q(lm) = (rhoPAWcorr(l,m) * rl).integrate (rad.logDr(is));
                  if (l == 0) Q(lm) -= SQRT_1_4PI * pawpot.nuclearCharge(is);
                  //cout << "Q(l=" << l << "; m=" << m << ")=" << Q(lm) << endl;

                  // add compensation charge to pseudo-density
                  if (hContrib & CalcNHatRad)
                     rhoRadPS(0)(l,m).plus_assign_ax(Q(lm), gRl(l));
               }
            }

            // --- Hartree energy and potentials
            dEdQrl(is)(ia).resize(rhoPAWcorr.meshData.nCols ());
            if (hContrib & CalcHartreeRad)  {
               for (int l = 0; l <= rhoPAWcorr.lmax; ++l)  {
                  for (int m = -l; m <= l; ++m)  {
                     v1AE(l,m) <<= radial.getHartreePotential (rhoRadAE(0)(l,m));
                     v1PS(l,m) <<= radial.getHartreePotential (rhoRadPS(0)(l,m));
                     // Ref. 1, Eq. 20, Hartree term without n^Z
                     double E =  0.5 * tr(rhoRadAE(0)(l,m) * v1AE(l,m));
                     eHartreeRadial += E;
                     eDoubleCounting -= E;
                     // Ref. 1, Eq. 21, Hartree term
                     E = 0.5 * tr(rhoRadPS(0)(l,m) * v1PS(l,m));
                     eHartreeRadial -= E;
                     eDoubleCounting += E;
                     if (hContrib & CalcNHatRad)  {
                        int lm = SxYlm::combineLm(l,m);
                        // Ref. 1, Eq. 36, 3rd term
                        dEdQrl(is)(ia)(lm) = -tr(gRl(l) * v1PS(l,m));
                     }
                  }
               }
               if (pawpot.rhoCoreAE(is).getSize () > 0)
                  eDoubleCounting += tr(pawpot.rhoCoreAE(is) * v1AE(0,0));
               if (pawpot.rhoCorePS(is).getSize () > 0)
                  eDoubleCounting -= tr(pawpot.rhoCorePS(is) * v1PS(0,0));
               eR0 += rhoRadAE(0)(0,0,0) * v1AE(0,0,0) * r0cub / 3.;
               //cout << "H r0 AE="
               //     << (rhoRadAE(0)(0,0,0) * v1AE(0,0,0) * r0cub / 3.) << endl;
               eR0 -= rhoRadPS(0)(0,0,0) * v1PS(0,0,0) * r0cub / 3.;
               //cout << "H r0 PS="
               //     << (-rhoRadPS(0)(0,0,0) * v1PS(0,0,0) * r0cub / 3.) << endl;
               if (!(hContrib & CalcNHatRad)) dEdQrl(is)(ia).set (0.);
            } else {
               v1AE.set (0.);
               v1PS.set (0.);
               dEdQrl(is)(ia).set (0.);
            }

            if (hContrib & CalcNuc)  {
               // Ref. 1, Eq. 20, Hartree n^Z parts
               // \int r^2 dr rho_AE(r) * (-Z)/r = -Z \int d(log r) r^2 rho_AE(r)
               eHartreeRadial -= SQRT_4PI * pawpot.nuclearCharge(is)
                                 * (rhoRadAE(0)(0,0)  * rad.radFunc(is).sqr ())
                                 .integrate (logDr);
               eDoubleCounting -= SQRT_4PI * pawpot.nuclearCharge(is)
                                 * (pawpot.rhoCoreAE(is) * rad.radFunc(is).sqr ())
                                 .integrate (logDr);
               // Ref. 1, Eq. (38): n^Z part
               v1AE(0,0) -= (SQRT_4PI * pawpot.nuclearCharge(is))
                               / rad.radFunc(is);

               // --- Estimate integral 0..r0 from rho(r0) and rho'(r0)
               double r0  = rad.radFunc(is)(0),
                     r1  = rad.radFunc(is)(1),
                     r01 = r0 - r1,
                     w0  = logDr / 3.,      // includes Simpson weight
                     w1  = 4. * logDr / 3., // includes Simpson weight
                     rho0 = rhoRadAE(0)(0,0,0),
                     rho1 = rhoRadAE(0)(1,0,0),
                     Z    = pawpot.nuclearCharge(is);
               eHartreeRadial -= SQRT_4PI * Z * (rho0 * 0.5 * r0 * r0
                     - (rho0  -rho1 ) / r01 * r0cub/6.);
               eDoubleCounting -= SQRT_4PI * Z *
                                  (pawpot.rhoCoreAE(is)(0) * 0.5 * r0 * r0
                     - (pawpot.rhoCoreAE(is)(0)  -pawpot.rhoCoreAE(is)(1) ) / r01 * r0cub/6.);
               v1AE(0,0,0) -= SQRT_4PI * Z *
                     (0.5 / r0 - 1./(6. * r01)) / w0;
               v1AE(1,0,0) -= SQRT_4PI * Z * 1./(6. * r01)
                                 *r0cub / (r1*r1*r1) / w1;
            }
            //eVnc += tr (VncAE * (rhoRadAE(0)(0,0) - pawpot.rhoCoreAE(is)));
            //eVnc -= tr (VncPS * (rhoRadPS(0)(0,0) - pawpot.rhoCorePS(is)));

            // Missing part vBar in Ref. 1, Eq. (39) (cf. Eq. 21)
            if (hContrib & CalcBar)
               v1PS(0,0) += pawpot.vBar(is);

            if (hContrib & CalcXcRad)  {
               // Ref. 1, Eq. 38/39: Hartree + xc part
               for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
                  // Ref 1, Eq. 40 with Eq. 38/39 without vR0 (see below)
                  Vij(iSpin,is,ia)
                        = pawpot.getVMatEl (v1AE + vXcAE(iSpin), pawpot.phiAE(is))
                        - pawpot.getVMatEl (v1PS + vXcPS(iSpin), pawpot.phiPS(is));
               }
            } else  {
               // Ref. 1, Eq. 38/39: Hartree, but not xc part
               // Ref 1, Eq. 40 with Eq. 38/39 without vR0 (see below)
               Vij(0,is,ia) = pawpot.getVMatEl (v1AE, pawpot.phiAE(is))
                            - pawpot.getVMatEl (v1PS, pawpot.phiPS(is));
               if (nSpin == 2)
                  Vij(1,is,ia) <<= Vij(0,is,ia);
            }

            bool spinConstraint = nuA.getSize () > 0; // TODO: fix condition
            if (nSpin == 2) {
               SX_LOOP2(i,j)
                  MSpin(iTlAtom) += omegaij(i,j) * (  Dij(0,is,ia)(i,j) 
                                                    - Dij(1,is,ia)(i,j));
            }
            if (spinConstraint)  {
               Vij(0,is,ia) -= nuA(iTlAtom) * omegaij;
               Vij(1,is,ia) += nuA(iTlAtom) * omegaij;
               eDoubleCounting += nuA(iTlAtom) * MSpin(iTlAtom);
            }

            if (hubbardU && (hContrib & CalcHubbardU))  {
               SX_LOOP(iSpin)
                  Vij(iSpin,is,ia) 
                     += hubbardU->computeAtom (iTlAtom, Dij(iSpin,is,ia),
                                               (nSpin == 1 ? 2. : 1.));
            }

            // --- compute compensation charges in G-space
            if (hContrib & CalcNHatG)  {
               SX_CLOCK (NHatG);
               PsiG T = gBasis.getPhaseFactors (is,ia);
               int nlm = sqr(lMax+1);
               // --- compute prefactors
               SxArray<SxComplex16> weightLm(nlm), normLm(nlm);
               double prefac = 1.;
               for (int l = 0; l <= lMax; ++l)  {
                  prefac /= 2. * l + 1.;
                  for (int m = -l; m <= l; ++m)  {
                     int lm = SxYlm::combineLm(l,m);
                     double N = rNorm * SxYlm::getYlmNormFactor(l,m) * prefac;
                     // i^l factor, see also Ref. 2, Eq. 3
                     SxComplex16 il = 1.;
                     if ((l & 3) == 1) il = I;
                     else if ((l & 3) == 2) il = -1.;
                     else if ((l & 3) == 3) il = -I;
                     weightLm(lm) = N * Qrl(is)(ia)(lm) * il;
                     normLm(lm) = N * il;
                  }
               }
#ifdef USE_OPENMP
               SxArray<SxArray<SxComplex16> > vPsGlmPerThread;
               int nThread = omp_get_max_threads ();
               vPsGlmPerThread.resize (nThread);
#pragma omp parallel
               {
                  SxArray<SxComplex16>& vPsGlm 
                     = vPsGlmPerThread(omp_get_thread_num ());
                  vPsGlm.resize (nlm);
#else
                  SxArray<SxComplex16> vPsGlm(nlm);
#endif
                  vPsGlm.set (0.);
#ifdef USE_OPENMP
#pragma omp for
#endif
                  for (int ig = 0; ig < gBasis.ng; ++ig)  {
                     // radial shape & atomic phase factor
                     SxComplex16 cHard = gaussHard(ig) * T(ig);
                     SxComplex16 cSoft = gaussSoft(ig) * T(ig);
                     SxComplex16 vGaussT;
                     // radial shape & atomic phase factor * potential
                     if (ig > 0) vGaussT = vRhoPSG(ig-1).conj () * cHard;
                     else        vGaussT = 0.;
                     // --- now multiply by YlmGl(ig,lm) and collect
                     for (int lm = 0; lm < nlm; ++lm)  {
                        double ylmgl = YlmGl(ig, lm);
                        if (hContrib & CalcNHatGHard) {
                           nHatG(ig)  += weightLm(lm) * ylmgl * cHard;
                           vPsGlm(lm) += ylmgl * vGaussT;
                        }
                        if (hContrib & CalcNHatGSoft)
                           nHatpG(ig) += weightLm(lm) * ylmgl * cSoft;
                     }
                  }
#ifdef USE_OPENMP
               }
               // --- sum over threads
               for (int iThread = 0; iThread < nThread; ++iThread)  {
                  const SxArray<SxComplex16>& vPsGlm = vPsGlmPerThread(iThread);
                  if (vPsGlm.getSize () == 0) continue;
                  for (int lm = 0; lm < nlm; ++lm)
                     dEdQrl(is)(ia)(lm) += (vPsGlm(lm) * normLm(lm)).re;
               }
#else
               for (int lm = 0; lm < nlm; ++lm)
                  dEdQrl(is)(ia)(lm) += (vPsGlm(lm) * normLm(lm)).re;
#endif
               prefac = 1.;

               if ((hContrib & CalcNHatGHard) && (hContrib & CalcNHatGSoft))  {
                  // absolute alignment of vHat (depends on rc2/rcp2)
                  // 4 PI/G^2 * Qrl * ( grl - gprl)
                  // G -> 0: 4 PI/G^2 * Qrl * N * G^L (1 - rc2/4 G^2 - 1 + rcp2/4 G^2)
                  //         = PI * N * Qrl * G^l (rcp2 - rc2)
                  // non-zero limit for l==0
                  prefac = PI * rNorm * SxYlm::getYlmNormFactor(0,0) * (rcp2 - rc2);
                  // for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {
                  vHatGZero += prefac * Qrl(is)(ia)(0);
                  if (hContrib & CalcHartreePS)
                     dEdQrl(is)(ia)(0) += prefac * rhoPSG(0).re;
                  // }
               }
            }
         } // LoopMPI
      } // ia loop
   } // is loop

   if (hubbardU)  {
      eDoubleCounting += hubbardU->eDoubleCounting;
      hubbardU->energy = SxLoopMPI::sum(hubbardU->energy);
   }

   {
      // --- MPI summations
      eXcRadial       = SxLoopMPI::sum (eXcRadial);
      eHartreeRadial  = SxLoopMPI::sum (eHartreeRadial);
      eBarRadial      = SxLoopMPI::sum (eBarRadial);
      eR0             = SxLoopMPI::sum (eR0);
      vHatGZero       = SxLoopMPI::sum (vHatGZero);
      eDoubleCounting = SxLoopMPI::sum (eDoubleCounting);
      SxLoopMPI::sum (nHatG);
      SxLoopMPI::sum (nHatpG);
      SxLoopMPI::sum (MSpin);
   }
   // ---
   SX_STOP_TIMER(Loop1);

   /// No MPI parallelization
   if (hubbardU && (hContrib & CalcHubbardU))  {
      hubbardU->eDoubleCounting = 0.;
      hubbardU->computeMO (pawRho.blockAO, structure);
      eDoubleCounting += hubbardU->eDoubleCounting;
   }

   if (nuA.getSize ()) SX_LOOP(ia) {
      cout << "nu(" << ia << ") = " << nuA(ia) << endl;
      cout << "Spin of atom " << ia << " = " << MSpin(ia) << endl;
   }

//   eVnc = SxLoopMPI::sum (eVnc);
//   sxprintf ("eVnc = %.12f\n", eVnc);

#ifdef USE_LOOPMPI
   // sync' Qrl between MPI processes
   if (SxLoopMPI::nr () > 1)
   {
      for (int is = 0; is < structure.getNSpecies (); ++is) {
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia) {
            int size = 0;
            int iatom = structure.getIAtom (is, ia);
            int sender = SxLoopMPI::whoseWork (iatom);
            if (SxLoopMPI::myWork (iatom)) {
               size = (int)Qrl(is)(ia).getSize ();
            }
            size = SxLoopMPI::bcast (size, sender);
            if ((!SxLoopMPI::myWork (iatom)) && (size != Qrl(is)(ia).getSize ())) {
               Qrl(is)(ia).resize (size);
            }
            SxLoopMPI::bcast (Qrl(is)(ia), sender);
         }
      }
   }
#endif

   SX_START_TIMER (PsStuff);

   // --- the pseudo-density related parts
   PsiG rhoHartreeHard = rhoPSG + nHatG,
        rhoHartreeSoft = rhoPSG + nHatpG;

   double sqrtVol = sqrt(structure.cell.volume);

   if (rhoHartreeHard(0).re * sqrtVol > 1e-8)  {
      cout << "Hartree charge = " << (rhoHartreeHard(0).re * sqrtVol) << endl;
   }

   PsiG rhoCoreG;
   if (hContrib & (CalcHartreePS | CalcBar | CalcNHatG))
      rhoCoreG = gBasis | computeCorePS ();

   double eHartreePS = 0.;

   PsiG vHartreeHard(gBasis);
   // Ref. 1, Eq. 34, first & second term (cf. Eq. 27)
   if (hContrib & CalcHartreePS)  {
      vHartreeHard(0) = vHatGZero;
      if (hContrib & CalcNHatGHard)
         vHartreeHard(noZero) = FOUR_PI * rhoHartreeHard(noZero)
                              / gBasis.g2(noZero);
      else
         vHartreeHard(noZero) = FOUR_PI * rhoPSG(noZero)
                              / gBasis.g2(noZero);
   } else {
      vHartreeHard.set (0.);
   }

   // --- vBar terms (Ref. 1, Eq. 34, third term)
   PsiG vBarTotal(gBasis);
   vBarTotal.set (0.);
   if (hContrib & CalcBar)  {
      for (int is = 0; is < structure.getNSpecies (); ++is)
         vBarTotal += vBarG(is) * gBasis.structureFactors(is);

      // --- Ref. 1, Eq. 19/21 vBar term (without ps core)
      eBar = dot(rhoPSG - rhoCoreG, vBarTotal) - eBarRadial;
   } else {
      eBar = 0.;
   }

   // Hartree terms in Ref. 1, Eq. 34
   SxMeshR vPSR = (vHartreeHard + vBarTotal).to (rBasis);

   // Initialization of dipole potential   
   double eDipole = 0.;
   SxMeshR vDipole;
   double vAlign = 0.;
   if (dipoleCorr && (hContrib & CalcDipoleC))  {
      vDipole = SxMeshR(rBasis);
      vDipole.set (0.);
      dipoleCorr->charge = rhoHartreeSoft(0) * sqrtVol;
      if (fabs(dipoleCorr->charge) < 1e-6)
         dipoleCorr->charge = 0.;
#ifdef USE_LOOPMPI
      dipoleCorr->charge = SxLoopMPI::bcast (dipoleCorr->charge, 0);
#endif
      SxMeshR rhoHartreeR = rhoHartreeSoft.to (rBasis);
      dipoleCorr->update (rhoPSG.to (rBasis), rhoHartreeR);
      dipoleCorr->correctPotential (&vDipole);
      // correction to potential
      vPSR += vDipole;

      // correction to energy
      eDipole =  0.5 * tr(rhoHartreeSoft.to (rBasis) * vDipole);
      eDoubleCounting -= eDipole;
      eDoubleCounting += tr(rhoCoreG.to (rBasis) * vDipole);

      if (fabs(dipoleCorr->charge) > 1e-6)  {
         // --- align to right virtual electrode at z=0
         vAlign = dipoleCorr->getVRZero (vPSR);
         cout << "vAlign=" << vAlign * HA2EV << " eV" << endl;
         vPSR -= vAlign;
         vDipole -= vAlign;
         eDipole         -= 0.5 * rhoHartreeSoft(0).re * sqrtVol * vAlign;
         eDoubleCounting += 0.5 * rhoHartreeSoft(0).re * sqrtVol * vAlign;
         eDoubleCounting -= rhoCoreG(0).re       * sqrtVol * vAlign;
         double electrodeEnergy = dipoleCorr->getElectrodeEnergy (rBasis.cell);
         eDipole += electrodeEnergy;
         eDoubleCounting += electrodeEnergy;
      } else {
         double eField = 0.5 * dipoleCorr->extraField * dipoleCorr->dipole;
         eDipole += eField;
         eDoubleCounting += eField;

      }
      sxprintf ("eDipole         = %18.12f H\n", eDipole);
   }


   SX_STOP_TIMER(PsStuff);

   if (hContrib & CalcNHatGSoft)
   {
      SX_CLOCK (Loop2);
      // --- compute dE/dQRl terms that depend on nHatp
      PsiG vHatpG(gBasis);
      vHatpG(0) = 0.;
      vHatpG(noZero) = FOUR_PI / gBasis.g2(noZero) * nHatpG(noZero);
      if (vExt.getSize () > 0) {
         vHatpG += gBasis | vExt;
      }
      if (dipoleCorr && (hContrib & CalcDipoleC))  {
         vHatpG += gBasis | vDipole;
      }
      double rcp2 = sqr(rSoft);
      SxDiracVec<Double> gaussSoft = exp ( -(0.25 * rcp2) * gBasis.g2);
      for (int is = 0; is < structure.getNSpecies (); ++is) {
         int lMax = pawpot.lMaxRho (is);
         int nlm = sqr(lMax+1);
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia) {
            if (SxLoopMPI::myWork(structure.getIAtom (is,ia))) { // SxParallelHierarchy
               PsiG T = gBasis.getPhaseFactors (is,ia);
#ifdef USE_OPENMP
               SxArray<SxArray<SxComplex16> > part;
               part.resize (omp_get_max_threads ());
#pragma omp parallel
               {
                  SxArray<SxComplex16>& myPart = part(omp_get_thread_num ());
                  myPart.resize (nlm);
#else
                  SxArray<SxComplex16> myPart(nlm);
#endif
                  myPart.set (0.);
#ifdef USE_OPENMP
#pragma omp for
#endif
                  for (int ig = 0; ig < gBasis.ng; ++ig)  {
                     SxComplex16 vGaussT = vHatpG(ig).conj ()
                                         * gaussSoft(ig)
                                         * T(ig);
                     for (int lm = 0; lm < nlm; ++lm)  {
                        // Ref. 1 Eq. 36, first line, 2nd term
                        myPart(lm) += YlmGl(ig + gBasis.ng * lm) * vGaussT;
                     }
                  }
#ifdef USE_OPENMP
               }
               for (int iThread = 0; iThread < part.getSize (); ++iThread)  {
                  const SxArray<SxComplex16>& myPart = part(iThread);
                  if (myPart.getSize () == 0) continue;
#endif
               double prefac = 1.;
               for (int l = 0; l <= lMax; ++l)  {
                  prefac /= 2 * l + 1;
                  // i^l factor
                  SxComplex16 il = 1.;
                       if ((l & 3) == 1) il = I;
                  else if ((l & 3) == 2) il = -1.;
                  else if ((l & 3) == 3) il = -I;
                  for (int m = -l; m <= l; ++m)  {
                     double N = rNorm * SxYlm::getYlmNormFactor(l,m) * prefac;
                     int lm = SxYlm::combineLm(l,m);
                     // Ref. 1 Eq. 36, first line, 2nd term
                     dEdQrl(is)(ia)(lm) += N * (il * myPart(lm)).re;
                  }
               }
#ifdef USE_OPENMP
               }
#endif
            } // LoopMPI
         }
      }
   }


   // compute U and U-related part of dE/dQrl
   double eHartreeU = 0.;
   if (hContrib & CalcU)  {
      eHartreeU = computeU (&dEdQrl);
      //cout << "U = " << eHartreeU << endl;
   }


   /// --- Compute vR0 contribution
   SX_START_TIMER (Loop3);
   for (int is = 0; is < structure.getNSpecies (); ++is)
   {
      double eDoubleCountingPrime = 0.;

      {
         double Zps = -pawpot.nuclearCharge(is) / sqrt(FOUR_PI);
         if (pawpot.rhoCoreAE(is).getSize () > 0)
            Zps += tr(pawpot.rhoCoreAE(is));
         if (pawpot.rhoCorePS(is).getSize () > 0)
            Zps -= tr(pawpot.rhoCorePS(is));
         
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia) {
            if(SxLoopMPI::myWork(structure.getIAtom (is,ia))) { // SxParallelHierarchy
               eDoubleCountingPrime += dEdQrl(is)(ia)(0) * Zps;
            }
         }
      }

      int lMax = pawpot.lMaxRho (is);
      // int np = pawpot.getNProj (is);
      int nProjType = pawpot.getNProjType (is);

      SxRadialMesh vR0;
      vR0.resize (rad, is, lMax);

      for (int ia = 0; ia < structure.getNAtoms (is); ++ia) {
         if (SxLoopMPI::myWork(structure.getIAtom (is,ia))) { // SxParallelHierarchy

            for (int lm = 0, l = 0; l <= lMax; ++l)  {
               SxDiracVec<Double> rl = pow(rad.radFunc(is), double(l));

               // Ref. 1, Eq. 37
               for (int m = -l; m <= l; ++m, ++lm) {
                  vR0(l,m) <<= dEdQrl(is)(ia)(lm) * rl;
               }
            }


            // Ref 1, Eq. 40 with Eq. 38/39: vR0 term and kinetic energy
            SxMatrix<Double> vR0ij =   pawpot.getVMatEl (vR0, pawpot.phiAE(is))
                                     - pawpot.getVMatEl (vR0, pawpot.phiPS(is));
            for (int iSpin = 0; iSpin < nSpin; ++iSpin)   {
               SxMatrix<Double> &VijR = Vij(iSpin,is,ia);
               // add VR0 term
               VijR += vR0ij;
               // --- add kinetic Energy
               if (hContrib & CalcKinRad)  {
                  for (int ipt = 0; ipt < nProjType; ++ipt)  {
                     int l = pawpot.lPhi(is)(ipt), nm = 2 * l + 1;
                     int offsetI = pawpot.offset(is)(ipt);
                     for (int jpt = 0; jpt < nProjType; ++jpt)  {
                        if (l == pawpot.lPhi(is)(jpt))  {
                           int offsetJ = pawpot.offset(is)(jpt);
                           double Tij = pawpot.deltaKin(is)(ipt,jpt);
                           for (int m = 0; m < nm; ++m)
                              VijR(offsetI + m, offsetJ + m) += Tij;
                        }
                     }
                  }
               }
               // --- add exact core-valence exchange
               if ((hContrib & CalcCorValExch) && exchangePtr) {
                  for (int ipt = 0; ipt < nProjType; ++ipt)  {
                     int l = pawpot.lPhi(is)(ipt), nm = 2 * l + 1;
                     int offsetI = pawpot.offset(is)(ipt);
                     for (int jpt = 0; jpt < nProjType; ++jpt)  {
                        if (l == pawpot.lPhi(is)(jpt))  {
                           int offsetJ = pawpot.offset(is)(jpt);
                           double Xij = SxXCFunctional::alphaHybrid
                                      * pawpot.coreX(is)(ipt,jpt);
                           for (int m = 0; m < nm; ++m)
                              VijR(offsetI + m, offsetJ + m) += Xij;
                        }
                     }
                  }
               }
            }
         } // LoopMPI
      }

#ifdef USE_LOOPMPI
      if (SxLoopMPI::nr () > 1) {
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia) {
            for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
               SxMatrix<Double> &VijR = Vij(iSpin,is,ia);
               int sender = SxLoopMPI::whoseWork (structure.getIAtom (is, ia));
               SxLoopMPI::bcast (VijR, sender);
            }
         }
         eDoubleCountingPrime = SxLoopMPI::sum (eDoubleCountingPrime);
      }
#endif

      eDoubleCounting += eDoubleCountingPrime;
   } // end is loop
   SX_STOP_TIMER (Loop3);


   SX_CLOCK (PsStuff);

   if (hContrib & CalcHartreePS)  {
      // Ref. 1, Eq. 25, first term via Ref. 1, Eq. 26
      eHartreePS = TWO_PI * (rhoHartreeSoft.absSqr ()(noZero)
                             / gBasis.g2(noZero)).sum ();
      eDoubleCounting -= eHartreePS;
      eDoubleCounting += FOUR_PI * dot(rhoCoreG(noZero), rhoHartreeSoft(noZero)
                                                       / gBasis.g2(noZero)).re;
   } else if (hContrib & CalcNHatGSoft)  {
      eHartreePS = TWO_PI * (nHatpG.absSqr ()(noZero)
                             / gBasis.g2(noZero)).sum ();
      eDoubleCounting -= eHartreePS;
      eDoubleCounting += FOUR_PI * dot(rhoCoreG(noZero), nHatpG(noZero)
                                                       / gBasis.g2(noZero)).re;
   }
   //sxprintf ("eHartreePS      = %18.12f H\n", eHartreePS);
   
   // Ref. 1, Eq. 27 in G-space
   PsiG vHatG(gBasis);
   vHatG(0) = vHatGZero; // see above
   vHatG(noZero) = FOUR_PI/gBasis.g2(noZero) * (nHatG - nHatpG)(noZero);

   if (hContrib & CalcHartreePS)  {
      // Ref. 1, Eq. 25, second term
      double eVHat = dot(rhoPSG, vHatG).re;
      //sxprintf ("eVHat           = %18.12f H\n", eVHat);
      eHartreePS += eVHat;
      eDoubleCounting -= eVHat;
      eDoubleCounting += dot(rhoCoreG, vHatG).re;
   }

   if (hContrib & CalcU)  {
      // Ref. 1, Eq. 25, third term
      //sxprintf ("eHartreeU       = %18.12f H\n", eHartreeU);
      eHartreePS += eHartreeU;
   }


   // compose total Hartree energy
   //sxprintf ("eHartreeRadial  = %18.12f H\n", eHartreeRadial);
   eHartree = eHartreePS + eHartreeRadial + eDipole;

   // SX_START_TIMER (PsStuff4);
   // doesn't take much time

   if (vExt.getSize () > 0)  {
      // --- external potential
      eExt = tr (rhoHartreeSoft.to (rBasis) * vExt);
      eDoubleCounting += tr(rhoCoreG.to (rBasis) * vExt);
      vPSR += vExt;
   } else {
      eExt = 0.;
   }

#ifdef USE_LOOPMPI
   SX_MPI_MASTER_ONLY  {
      vPSR = rBasis.symmetrize (vPSR);
   }
   SxLoopMPI::bcast (vPSR, 0);
#else
   vPSR = rBasis.symmetrize (vPSR);
#endif

   // SX_STOP_TIMER (PsStuff4);
   // doesn't take much time


   // SX_START_TIMER (PsStuff5);

   SX_MPI_MASTER_ONLY
   {
      RhoR vElStat(1);
      vElStat(0) = vPSR * HA2EV;
      SxRho(vElStat).writeRho ("vElStat-eV.sxb");
   }

   // SX_STOP_TIMER (PsStuff5);
   // doesn't take much time



   SX_START_TIMER (PsStuff1);

   vPS.resize (nSpin);
   // --- xc part of pseudo-density
   if (hContrib & CalcXcPS)  {
      if (pawpot.kjXcShape.getSize () > 0)  {

         SX_START_TIMER (PsStuff2);

         PsiG rhoXcG(gBasis);
         rhoXcG.set (0.);

         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            if (pawpot.kjXcShape(is).getSize () > 0)  {
               cout << "Using Kresse-Joubert formulation of xc for "
                    << pawpot.prettyName(is) << endl;
               double Qcore = -pawpot.nuclearCharge(is) / sqrt(FOUR_PI)
                            +  ((pawpot.rhoCoreAE(is) - pawpot.rhoCorePS(is))
                                * rad.radFunc(is).cub ()).integrate (rad.logDr(is));
               PsiG kjxcG = gBasis | pawpot.kjXcShape(is);
               for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                  double Q = Qrl(is)(ia)(0) - Qcore;
                  //cout << "Q=" << Q << endl;
                  PsiG T = gBasis.getPhaseFactors(is,ia);
                  rhoXcG.plus_assign_ax (Q, T * kjxcG);
               }
            }
         }
         rhoXcG /= nSpin;
         RhoR rhoXc(nSpin);
         rhoXc(0) = rBasis | rhoXcG;
         //cout << "nXcKj=" << structure.cell.volume * rhoXc(0).sum ()
         //        / rhoXc(0).getSize () << endl;
         if (nSpin == 2)
            rhoXc(1) = rhoXc(0) + rhoPSR(1);
         rhoXc(0) += rhoPSR(0);
         //cout << "nPs=" << structure.cell.volume * rhoPSR(0).sum ()
         //        / rhoPSR(0).getSize () << endl;
         cout << "nXc=" << structure.cell.volume * rhoXc(0).sum ()
                 / (PrecTauR)rhoXc(0).getSize () << endl;
         xcPtr->updateXC (rhoXc);

         SX_STOP_TIMER (PsStuff2);

      } else {
         SX_START_TIMER (PsStuff3);
         xcPtr->updateXC (rhoPSR);
         SX_STOP_TIMER (PsStuff3);
      }
      {
         SxMeshR rhoCoreR = computeCorePS ();
         if (nSpin == 2) rhoCoreR *= 0.5;
         SX_LOOP(iSpin)
            eDoubleCounting -= dot(rhoPSR(iSpin) - rhoCoreR,
                                   xcPtr->vXc(iSpin)) * rBasis.dOmega;
      }
      eDoubleCounting += xcPtr->eXc;

      SX_START_TIMER (PsStuff4);
      eXc = eXcRadial + xcPtr->eXc;
      SX_STOP_TIMER (PsStuff4);

      SX_START_TIMER (PsStuff5);
      // Ref. 1, Eq. 34, last term
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         vPS(iSpin) = vPSR + xcPtr->vXc(iSpin);
      SX_STOP_TIMER (PsStuff5);

   } else {
      SX_START_TIMER (PsStuff6);
      eXc = eXcRadial;
      for (int iSpin = 0; iSpin < nSpin; ++iSpin) vPS(iSpin) = vPSR;
      SX_STOP_TIMER (PsStuff6);
   }

   SX_STOP_TIMER (PsStuff1);

   eCore = 0.;
   if (hContrib & CalcCore)
      eCore = (structure.atomInfo->nAtoms * potPtr->coreEnergy).sum ();
   eDoubleCounting -= eCore;

#ifdef USE_LOOPMPI
   {
      double eDCmaster = SxLoopMPI::bcast (eDoubleCounting, 0);
      if (eDCmaster != eDoubleCounting)  {
         cout << "eDC mismatch: " << (eDoubleCounting - eDCmaster) << endl;
      }
      eDoubleCounting = eDCmaster;
   }
#endif

   //cout << "eXc = " << eXc << endl;
   //cout << "eExt = " << eExt << endl;
   //cout << "eHartree = " << eHartree << endl;
   //cout << "eR0 = " << eR0 << endl;
}


PrecEnergy SxPAWHamiltonian::getEnergy (const SxPsiSet &pw,
                                        const SxFermi &fermi)
{
   // TODO: rPtr => SxPtr
   SX_CHECK (rPtr);
   SX_CHECK(dynamic_cast<const SxPWSet *>(&pw));
   compute (dynamic_cast<const SxPWSet &>(pw), fermi, false);
   return eTotal;
}

double SxPAWHamiltonian::getHarrisEnergy (const SxFermi &fermi)
{
   double eHarrisFoulkes = eDoubleCounting;
   //cout << "eDoubleCounting = " << eDoubleCounting << endl;
   SX_LOOP2(iSpin,ik)
      eHarrisFoulkes += dot(fermi.eps(iSpin,ik),
                            fermi.focc(iSpin,ik)) * fermi.kpPtr->weights(ik);
   //cout << "eBand = " << (eHarrisFoulkes - eDoubleCounting) << endl;
   return eHarrisFoulkes;
}

void SxPAWHamiltonian::compute (const SxPWSet &waves, const SxFermi &fermi,
                                bool computeNewRho)
{
   SX_CLOCK(Timer::ComputeHam);
   SX_CHECK (potPtr);
   SX_CHECK (xcPtr);

   // compute new density
   if (computeNewRho) computeRho (fermi.focc, waves);

   //sxprintf ("eHarrisFoulkes = %20.16f\n", getHarrisEnergy (fermi));

   // --- Hartree & xc parts
   computeRhoTerms ();


   if (exchangePtr && (hContrib & CalcValExch))  {
      double eX = exchangePtr->compute (waves.getThis (), fermi.focc,
                                        pawRho.Dij);
      //cout << "eX = " << eX << endl;
      eXc += SxXCFunctional::alphaHybrid * eX;
      eDoubleCounting -= SxXCFunctional::alphaHybrid * eX;
   }

   // kinetic energy
   double eKin = 0.;
   // --- pseudo-kinetic enery
   int nSpin = waves.getNSpin ();

   if (hContrib & CalcKinPS)  {
      SX_NEW_LOOP (waves);
      SxLaplacian L;
      for (int ik = 0; ik < waves.getNk (); ++ik)  {
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            // Ref. 1, Eq. 19, kinetic energy
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik)) {
               double weight = fermi.kpPtr->weights(ik);
               const SxGBasis &gk = waves.getGkBasis ()(ik);
               for (int iState = 0; iState < waves.getNStates(ik); ++iState)  {
                  double focc = fermi.focc(iState,iSpin,ik);
                  if (fabs(focc) < 1e-16) continue;
                  PsiG psiG = gk | waves(iState,iSpin,ik);
                  eKin += 0.5 * weight * focc * (psiG | L | psiG);
               }
            } // LoopMPI
         }
      }
      SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
      SX_MPI_TARGET (TopLevel, TaskGroupAll);
      eKin = SxLoopMPI::sum(eKin);
   }

   bool corValEx = exchangePtr && (hContrib & CalcCorValExch);
   bool kin = (hContrib & CalcKinRad );

   if (kin || corValEx)  {
      // --- 1-center correction of kinetic energy
      const SxPAWPot &pot = *potPtr;
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         // Ref. 1, Eq. 20/21 with Eq. 52
         // \sum_ij Dij * ( <phiAE_j| -1/2 nabla^2 | phiAE_i>
         //                -<phiPS_j| -1/2 nabla^2 | phiPS_i>)
         int nProjType = pot.getNProjType (is);
         for (int ipt = 0; ipt < nProjType; ++ipt)  {
            int l = pot.lPhi(is)(ipt);
            int offsetI = pot.offset(is)(ipt) + l;
            for (int jpt = 0; jpt < nProjType; ++jpt)  {
               if (l == pot.lPhi(is)(jpt))  {
                  int offsetJ = pot.offset(is)(jpt) + l;
                  double mSumD = 0.;
                  for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
                     for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                        for (int m = -l; m <= l; ++m)  {
                           mSumD += pawRho.Dij(iSpin,is,ia)(offsetI + m,
                                                            offsetJ + m);
                        }
                     }
                  }
                  if (kin)
                     eKin += mSumD * pot.deltaKin(is)(jpt, ipt);
                  if (corValEx)
                     eXc += SxXCFunctional::alphaHybrid
                          * mSumD * pot.coreX(is)(jpt, ipt);
               }
            }
         }
      }
   }


   sxprintf ("eKin      = % 19.12f H\n", eKin);
   sxprintf ("eHartree  = % 19.12f H\n", eHartree);
   sxprintf ("eXc       = % 19.12f H\n", eXc);
   sxprintf ("eBar      = % 19.12f H\n", eBar);
   sxprintf ("eCore     = % 19.12f H\n", eCore);
   if (applyVDWCorrection) {
      SxArray<SxVector3<Double> > newCoords (structure.nTlAtoms);
      SxArray <double> effectiveVolume = getEffectiveVolume(vdwCorrection.correctionType);
      for (int i=0; i<structure.nTlAtoms; i++) {
          newCoords(i) = structure.ref(i);
      }
      vdwCorrection.update(newCoords, effectiveVolume);
      vdwCorrection.compute();
      eVDW = vdwCorrection.totalEnergy;
      sxprintf ("eVDW     = % 19.12f H\n", eVDW);
   }
   if (vExt.getSize () > 0)
      sxprintf ("eExt      = %19.12f H\n", eExt);
   eTotal = eKin + eHartree + eXc + eBar - eCore + eExt;
   if (applyVDWCorrection) {
      eTotal = eTotal + eVDW;
      eDoubleCounting = eDoubleCounting + eVDW;
   }
   if (hubbardU) {
      SX_CHECK (fabs(fermi.fFull * nSpin - 2.) < 1e-12, fermi.fFull, nSpin);
      sxprintf ("eHubbard  = % 19.12f H\n", hubbardU->energy);
      eTotal += hubbardU->energy;
   }

   sxprintf ("eTot(Val) = % 19.12f H\n", eTotal);

}

SxArray<Double> SxPAWHamiltonian::getEffectiveVolume (SxString correctiontype)
{
   SxArray<Double> effectivevolume(structure.nTlAtoms);
   if (correctiontype == SxString("TS")) {
      SxArray<double> hirshfeldVolume = pawRho.getHirshfeldVolume ();
      for (int i=0; i<structure.nTlAtoms; i++) {
         effectivevolume(i) = hirshfeldVolume(i) / freeVolume(structure.elements(i));
         //effectivevolume(i) = 1.;
      }

   }
   else {
      for (int i=0; i<structure.nTlAtoms; i++) {
         effectivevolume(i) = 1.;
      }
   }

   return effectivevolume;
}

double SxPAWHamiltonian::computeU (SxArray<SxArray<SxVector<Double> > > *dEdQrl)
{
   SX_CLOCK(Timer::ComputeU);
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;
   
   // --- for neighbor search
   SxGrid grid (structure, 10);
   SxNeighbors neighbors;
   int neighborParams = SxNeighbors::StoreRel
                      | SxNeighbors::IncludeZeroDistance;
   // note: for rcut=13, 1 bohr wide Gaussians behave like point multipoles
   // within 1e-16
   // we assume that multipole Gaussians are always below 1 bohr wide
   double rcut = 13. * rSoft;

   double U = 0.;
   const SxYlm::SxClebschTable &cg = pot.clebschGordan;
   SxGaussIJ gaussIJ;

   // --- loop over atoms (index i => R in Bloechl notation)
   for (int is = 0; is < structure.getNSpecies (); ++is)  {
      int lmaxI = pot.lMaxRho(is);
      for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
         if(SxLoopMPI::myWork(structure.getIAtom (is, ia)))
         {
            // compute neighbors within rcut
            neighbors.compute (grid, structure, structure.getAtom(is,ia),
                  rcut, neighborParams);

            // --- loop over neighbors j (R' in Bloechl notation)
            for (int js=0; js < structure.getNSpecies (); js++)  { // --- species j
               int lmaxJ = pot.lMaxRho(js);
               double a2Hard = sqr(pot.rc(is)) + sqr(pot.rc(js));
               double a2Soft = 2. * sqr(rSoft); // sqr(1.) + sqr(1.);
               for (int ja = 0; ja < neighbors.relPositions.getNAtoms(js); ++ja)
               {
                  // i => R
                  // j => R'
                  // r = R'-R
                  Coord r = neighbors.relPositions.getAtom(js, ja);
                  // Ref. 1 Eq. (28): kernel
                  gaussIJ.setDelta (r, a2Hard, a2Soft, lmaxI + lmaxJ);
                  SxMatrix<Double> delta = gaussIJ.compute (lmaxI, lmaxJ, cg);
                  // sum over moments of j=R'
                  int jaGlobal = neighbors.relPositions.getIAtom(js,ja);
                  int jRefAtom = neighbors.idx(jaGlobal)
                               - structure.atomInfo->offset(js);
                  SxVector<Double> sumJ = delta ^ Qrl(js)(jRefAtom);

                  // Ref. 36, 2nd term (dU/d Qrl)
                  if (dEdQrl) (*dEdQrl)(is)(ia) += sumJ;
                  // Ref. 1 Eq. (28)
                  U += 0.5 * (Qrl(is)(ia) * sumJ).sum ();
               }
            }
         } // LoopMPI
      }
   }

//   eDoubleCounting -= U;
//   U = SxLoopMPI::sum(U);
   
   U = SxLoopMPI::sum(U);
   eDoubleCounting -= U;
   
   return U;
}

void SxPAWHamiltonian::setupPSRefG(const SxGBasis &gBasis)
{
   SX_CHECK(potPtr);
   // pot must have a radial basis
   SX_CHECK(potPtr->getBasisPtr ());
   int nSpecies = potPtr->getNSpecies ();
   vBarG.resize (nSpecies);
   rhoCorePSG.resize (nSpecies);
   for (int is = 0; is < nSpecies; ++is)  {
      vBarG(is)      = (gBasis | potPtr->vBar(is));
      rhoCorePSG(is) = (gBasis | potPtr->rhoCorePS(is));
   }
}

PsiG SxPAWHamiltonian::apply (const PsiG &psi) const
{
   SX_CLOCK (Timer::PAW_H_psi);
   SX_CHECK (psi.handle);
   int iSpin = psi.handle->auxData.iSpin;
   SX_CHECK (psi.nCols () > 0, psi.nCols ());
   int nStates = (int)psi.nCols ();
   const SxGBasis *gkPtr = dynamic_cast<const SxGBasis *> (psi.getBasisPtr ());
   if (!gkPtr)  {
      const SxPAWBasis *paw 
         = dynamic_cast<const SxPAWBasis *> (psi.getBasisPtr ());
      SX_CHECK (paw);
      gkPtr = paw->gBasis.getPtr ();
   }
   SX_CHECK (gkPtr);
   const SxGBasis &gk = *gkPtr;
   SX_CHECK (rPtr);
   const SxRBasis &R = *rPtr;
   SX_CHECK (pBasis);
   const SxPartialWaveBasis &p = *pBasis;

   // Ref. 1, Eq. 40, |pi>(...)<pj term
   SX_START_TIMER(Timer::PAW_H_nl);
   PsiG res = (gk | p) * (p | (Vij ^ (p | psi) ));
   SX_STOP_TIMER(Timer::PAW_H_nl);

   if (hubbardU && hubbardU->moSite.getSize () > 0)  {
      SX_CLOCK (Timer::PAW_H_hubMO);
      SX_LOOP(i)  {
         SxBlockDensityMatrix hMo = hubbardU->moSite(i)->hamProj;
         SxAOBasis pMo = *hubbardU->moSite(i)->aoProj;
         res += gk | (hMo ^ (pMo | psi) );
      }
   }

#ifdef USE_FFT2d1d
   bool direct = true;
#else
   bool direct = false;
#endif
   if (direct && (hContrib & CalcVPS))  {
      SX_CLOCK(Timer::PAW_H_loc);
      res += gk.convolute (psi, vPS(iSpin));
   }

   for (int iState = 0; iState < nStates; ++iState)  {
      PsiG psiG = gk | psi.colRef(iState);
      if (hContrib & CalcKinPS)  {
         // pseudo-kinetic energy (Ref. 1, Eq. 40, first term)
         res.colRef(iState).plus_assign_ax(0.5, gk.g2 * psiG);
      }
      if ((hContrib & CalcVPS) && !direct)  {
         SX_CLOCK(Timer::PAW_H_loc);
         // local potential       (Ref. 1, Eq. 40, second term)
         res.colRef(iState) += gk | (vPS(iSpin) * (R | psiG));
      }

   }
   if (exchangePtr && (hContrib & CalcValExch))  {
      PsiG xPsi = exchangePtr->apply (psi);
      res.plus_assign_ax (SxXCFunctional::alphaHybrid, xPsi);
      // --- synchronize across MPI tasks for this k-point
      SX_MPI_LEVEL("waves-k");
      SX_MPI_SOURCE(CurrentLevel,TaskGroupMaster);
      SX_MPI_TARGET(CurrentLevel,TaskGroupAll);
      SxLoopMPI::bcast (res, 0);
   }

   res.handle->auxData  = psi.handle->auxData;
   res.setBasis (&gk);
   return res;
}

SxDiracVec<TPrecCoeffG::TReal> SxPAWHamiltonian::preconditioner (const PsiG &dPsi) const
{
   const SxGBasis *Gk = dynamic_cast<const SxGBasis *>(dPsi.getBasisPtr());
   if (!Gk)  {
      const SxPAWBasis* pawBasis
         = dynamic_cast<const SxPAWBasis *>(dPsi.getBasisPtr ());
      SX_CHECK(pawBasis);
      Gk = pawBasis->gBasis.getPtr ();
   }
   SX_CHECK (Gk);

   double kin = (*Gk | dPsi).laplacian();
   /* --- beautiful and slow
   SxDiracVec<TPrecCoeffG::TReal> x, x2, x3, x4, n, K;
   x  = Gk->g2 / kin;
   x2 = x.sqr();
   x3 = x.cub();
   x4 = x2.sqr();
   n  = 27. + 18.*x + 12.*x2 + 8.*x3;
   K  = n / (n + 16.*x4);
   K /= kin; // make K approach 1/(G^2 - kin)
   */
   // --- less beautiful, but much faster
   double kinInv = 1. / kin;
   SxDiracVec<TPrecCoeffG::TReal> K;
   K.resize (Gk->ng);
   K.setBasis (Gk);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
   for (ssize_t ig = 0; ig < Gk->ng ; ++ig)   {
      double x_ = Gk->g2(ig) * kinInv;
      double x2_ = x_ * x_;
      double x3_ = x_ * x2_;
      double x4_ = x2_ * x2_;
      double n_ = 27. + 18. * x_ + 12. * x2_ + 8. * x3_;
      K(ig) = n_ / (kin * (n_ + 16. * x4_));
   }
   return K;
}

SxDiracMat<Complex16>
SxPAWHamiltonian::getHnm (const PsiGI &waveBlock) const
{
   ssize_t nStates = waveBlock.nCols ();
   int bSize = 64;
   if (waveBlock.getSize () < 1000000)  {
      PsiGI HPsi = apply(waveBlock);
      return waveBlock.overlap (HPsi, HPsi.nRows ());
   }
   SxDiracMat<Complex16> Hnm(nStates, nStates);
   ssize_t ngPAW = waveBlock.nRows ();
   for (ssize_t i = 0; i < nStates; i += bSize)  {
      ssize_t ni = i + bSize < nStates ? bSize : (nStates - i);
      SxDiracMat<Complex16> blockI = waveBlock.getBlock (0, i, ngPAW, ni);
      blockI.setBasis (waveBlock.getBasisPtr ());
      PsiGI HPsi = apply(blockI);
      Hnm.setBlock (waveBlock.overlap (HPsi, HPsi.nRows ()), 0, i);
   }
   return Hnm;
}

void SxPAWHamiltonian::computeForces (const SxPWSet &waves, SxFermi &fermi)
{
   SX_CLOCK(Timer::Force);
   bool printForceParts = false;
   SX_CHECK(rPtr);
   SX_CHECK(potPtr);

   forces = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   const SxGBasis &gBasis      = rPtr->getGBasis();
   const SxGkBasis &gk         = waves.getGkBasis();

   SxIdx noZero(1,gBasis.ng - 1);
   const SxDiracMat<Double> G  = SxDiracMat<Double>(gBasis.gVec);
   const SxDiracVec<Double> &G2 = gBasis.g2;

   int nSpin = waves.getNSpin();
   const SxRho &rhoPSR = pawRho.pwRho;

   PsiG rhoPSG = gBasis | ((nSpin == 1) ? rhoPSR(0)
                                        : rhoPSR(0) + rhoPSR(1));

   // |G|^l * Ylm(G)
   SX_START_TIMER (TimeYlmGl);
   int lMaxRho = potPtr->lMaxRho.maxval ();
   SxDiracMat<Double> YlmGl(gBasis.ng, sqr(lMaxRho + 1));
   YlmGl.setBasis (gBasis);
   for (int lm = 0, l = 0; l <= lMaxRho; ++l) {
      SxDiracVec<Double> gl = pow (gBasis.g2, 0.5 * l);
      for (int m = -l; m <= l; ++m,++lm)
         YlmGl.colRef(lm) <<= gBasis.getYlm (l,m) * gl;
   }
   SX_STOP_TIMER (TimeYlmGl);

   double rcp2 = sqr(rSoft);

   // khr:  the calculations of f11...f14 do not need to be parallelized

   SxAtomicStructure f11 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   SxAtomicStructure f12 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   SxAtomicStructure f13 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));

   SX_START_TIMER (Timer::F11);
   // --- calculate vHartreeSoft
   PsiG rhoHatPG (gBasis);
   rhoHatPG.set(0.0);
   if (hContrib & CalcNHatGHard) {
      for(int is = 0; is < structure.getNSpecies (); is++)
         for(int ia = 0; ia < structure.getNAtoms (is); ia++)
            rhoHatPG += computeRhoHatG (is,ia,YlmGl, rcp2);
   }

   PsiG vHartreeSoft(gBasis);
   vHartreeSoft(0)        = 0.;
   vHartreeSoft(noZero) <<= FOUR_PI * rhoHatPG(noZero) / G2(noZero);
   if (hContrib & CalcHartree)
      vHartreeSoft(noZero) += FOUR_PI * rhoPSG(noZero) / G2(noZero);

   if (dipoleCorr && (hContrib & CalcDipoleC))  {
      // --- add dipole field to potential
      SxMeshR vDipole(*rPtr);
      vDipole.set (0.);
      dipoleCorr->correctPotential (&vDipole);
      vHartreeSoft += vDipole.to (gBasis);
   }

   if (vExt.getSize () > 0)  {
      // add external field to potential
      vHartreeSoft += vExt.to (gBasis);
   }
   SX_STOP_TIMER (Timer::F11);

   SxMeshG vPSG = (gBasis | ((nSpin == 1) ? vPS(0)
                                          : vPS(0) + vPS(1)));

   // --- calculate f11 f12 f13 (Ref. 1, Eq. (50))
   SX_LOOP2(is,ia) {
      PsiG rhoHatRPG = computeRhoHatG ((int)is,(int)ia,YlmGl, rcp2);
      PsiG T = gBasis.getPhaseFactors ((int)is,(int)ia);
      // --- Ref.1 Eq (50) 1st term,
      {  SX_CLOCK (Timer::F11);
         f11.ref(is,ia) -= ((I * vHartreeSoft.conj() * rhoHatRPG) ^ G)
                           .real ().toVector3 ();
      }
      // --- Ref.1 Eq (50) 2nd term,
      {  SX_CLOCK (Timer::F12);
         PsiG rhoHatRG(gBasis), v(gBasis);
         double rc2 = sqr(potPtr->rc(is));
         // NHatG part
         if (hContrib & CalcNHatGHard)
            rhoHatRG = computeRhoHatG ((int)is, (int)ia, YlmGl, rc2);
         else
            rhoHatRG.set(0.0);
         if (!(hContrib & CalcNHatGSoft)) rhoHatRPG.set (0.);

         v(0)        = 0.;
         v(noZero) <<= FOUR_PI * (rhoHatRG(noZero)-rhoHatRPG(noZero)) / G2(noZero);
         // vBar part
         if (hContrib & CalcBar) v += vBarG(is) * T;
         f12.ref(is,ia) -= ((-I * rhoPSG * v.conj()) ^ G).real ().toVector3 ();
      }
      // Ref.1 Eq (50) 3rd term, contribution due to pseudopot and core density
      {  SX_CLOCK (Timer::F13);
         f13.ref(is,ia) += 1./nSpin * ((rhoCorePSG(is) * T * vPSG.conj()) ^ G)
                                      .imag ().toVector3 ();
      }
   }
   if (printForceParts) cout << "F11" << endl << f11 << endl;
   if (printForceParts) cout << "F12" << endl << f12 << endl;
   if (printForceParts) cout << "F13" << endl << f13 << endl;

   // Ref.1 Eq (50) 4th term, U-Term
   SX_START_TIMER (Timer::F14);
   SxAtomicStructure f14 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   if (hContrib & CalcU)  f14 -= computeGradU ();
   SX_STOP_TIMER (Timer::F14);
   if (printForceParts) cout << "F14" << endl << f14 << endl;

   forces = f11 + f12 + f13 + f14;

   // F2 force computed below
   // F2 forces needs gradDij, which is computed together with F3
   SxArray<SxRadMat> gradDij(3);
   for (int iDir = 0; iDir < 3; iDir++) {
      gradDij(iDir).resize (structure.atomInfo, *potPtr, nSpin);
      gradDij(iDir).set (0.0);
   }

   // F2-like forces for the Hubbard part, cf. ref. 3, Eq. (38b-40)
   // need analog to gradDij
   SxArray2<SxBlockDensityMatrix> gradPij;
   if (hubbardU && hubbardU->moSite.getSize () > 0)  {
      gradPij.reformat (hubbardU->moSite.getSize (), nSpin);
      SX_LOOP2(iMO, iSpin)  {
         gradPij(iMO, iSpin).resize (3, hubbardU->moSite(iMO)->getNSite ());
         int nAO = hubbardU->moSite(iMO)->getNAOperSite ();
         SX_LOOP2(iDir,iSite)  {
            gradPij(iMO, iSpin)(iDir, iSite).reformat (nAO, nAO);
            gradPij(iMO, iSpin)(iDir, iSite).set (0.);
         }
      }
   }

   // Ref.1 Eq (58)
   // khr:  parallelized over the k points
   SX_START_TIMER (Timer::F3);
   SxAtomicStructure f3 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   int nStates = waves.getNStates ();

   // k loop, parallelize
   for (int ik = 0; ik < waves.getNk (); ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < waves.getNSpin (); ++iSpin) {
         SxDiracVec<Complex16> P = *pBasis | waves(iSpin,ik);
         SxArray<SxDiracMat<Complex16> > gradP
         = projBasis->gradProject (gk(ik) | waves(iSpin,ik));

         // --- set up gradDij for F2
         SxArray<SxRadMat> gradDijBAK =
               computeUnSymmGradDij (fermi, iSpin, ik, P, gradP);
         SX_LOOP3(iDir,is,ia) {
            gradDij(iDir)(iSpin,is,ia) += gradDijBAK(iDir)(0,is,ia);
         }
         // --- set up gradPij for F2-like forces for Hubbard
         if (hubbardU)  {
            for (int iMO = 0; iMO < hubbardU->moSite.getSize (); ++iMO)  {
               const SxHubbardMO &mo = *hubbardU->moSite(iMO);
               // get <p|psi>
               SxDiracVec<Complex16> PMo = (*mo.aoProj) | waves(iSpin, ik);
               // get d/dR <p|psi>
               SxArray<SxDiracMat<Complex16> > gradPMo
                  = mo.aoProj->gradProject (waves(iSpin, ik), ik);
               gradPij(iMO, iSpin) += mo.computeGradPij (fermi, PMo, gradPMo);
            }
         }

         // --- setup Hmn
         SxDiracMat<Complex16> Hnm;
         if (   dynamic_cast<const SxPW*>(&waves) 
             || dynamic_cast<const SxPAWSet*>(&waves))
         {
            Hnm = getHnm (waves(iSpin, ik));
         } else {
            Hnm.reformat (nStates,nStates);
            // Hmn is a hermitian matrix
            for (int jState = 0; jState < nStates; ++jState) {
               PsiG HPsiJ = apply(waves(jState,iSpin,ik));
               for (int iState = jState; iState < nStates; ++iState) {
                  Hnm(iState,jState) = dot(gk(ik)|waves(iState,iSpin,ik),
                                           HPsiJ);
                  Hnm(jState,iState) = Hnm(iState,jState).conj();
               }
            }
         }

         SX_LOOP3(iDir, is, ia) {
            double sum = 0;
            SX_LOOP2(iState, jState) {
               double focc = 0.5 * (fermi.focc((int)iState,iSpin,ik) 
                                    + fermi.focc((int)jState,iSpin,ik));
               sum += (focc * Hnm(iState, jState)
                      * pawGradDot ((int)is, (int)ia, (int)jState, (int)iState,
                                    P, gradP(iDir))).re;
            }
            f3.ref(is,ia)(iDir) -= gk.weights(ik) * sum;
         }
      }
   }

   if (SxLoopMPI::nr () > 1) {
      SX_MPI_SOURCE("waves-k", TaskGroupMaster);
      SX_MPI_TARGET(TopLevel, TaskGroupAll);
      // --- gradDij
      for (int iDir = 0; iDir < 3; iDir++) {
         gradDij(iDir).sumMPI ();
      }
      // gradPij
      if (gradPij.getSize () > 0)
         SX_LOOP2(iMO, iSpin) gradPij(iMO, iSpin).sumMPI ();
      // f3
      SxLoopMPI::sum (f3);
   }

   SX_STOP_TIMER (Timer::F3);
   if (printForceParts) cout << "F3" << endl << f3 << endl;

   forces += f3;

   // Ref.1 Eq (57)
   // khr:  does not need to be parallelized
   SX_START_TIMER (Timer::F2);
   gradDij = symmetrizeGradDij(gradDij);
   SxAtomicStructure f2 = structure.getNewStr ().set(Coord(0.0,0.0,0.0));
   SX_LOOP4 (iDir, iSpin, is, ia)  {
      double sum = 0.0;
      SX_LOOP2 (jpl,ipl) {
         sum += gradDij(iDir)(iSpin,is,ia)(ipl,jpl)
              * Vij(iSpin,is,ia)(ipl,jpl);
      }
      f2.ref(is,ia)(iDir) += sum;
   }

   if (gradPij.getSize () > 0)  {
      SX_LOOP2(iMO, iSpin)  {
         const SxHubbardMO &mo = *hubbardU->moSite(iMO);
         mo.symmetrizeGradPij (&gradPij(iMO,iSpin), structure, ylmRot);
         f2 += mo.getForce (gradPij(iMO,iSpin),iSpin);
      }
   }
   SX_STOP_TIMER (Timer::F2);

   if (hubbardU)  {
      // AO->MO transformation forces
      // they arise because the MO sites depend on the relative
      // position of the atoms
      // cf. ref. 3, Eq. (38a), (41)
      for (int iMO = 0; iMO < hubbardU->moSite.getSize (); ++iMO)  {
         forces += hubbardU->moSite(iMO)->trafoForce;
      }
   }

   forces += f2;
   if (printForceParts) cout << "F2" << endl << f2 << endl;


   // force is - dE/dR
   forces *= -1.0;

   // --- final symmetrize forces
   SxForceSym forceSymmetrizer;
   forceSymmetrizer.setup(structure);
   forces = forceSymmetrizer * forces;
   if (SxLoopMPI::nr () > 1)  {
      // make sure that forces are consistent
      SxLoopMPI::bcast (forces.coords, 0);
   }

   cout << SX_SEPARATOR;
   cout << "Total Forces:" << endl;
   cout << SX_SEPARATOR;
   SX_LOOP2(is,ia) {
      cout << "Species: " << is << "\t"
           << "Atom: "<< ia << "\t"
           << forces.getAtom(is,ia) << "\t"
           << forces.getAtom(is,ia).norm() << endl;
   }
}


SxAtomicStructure SxPAWHamiltonian::computeGradU ()   
{
   // --- computes last summand Ref 1. eq (50)
   SX_CLOCK(Timer::ComputeGradU);
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;

   // --- for neighbor search
   SxGrid grid (structure, 10);
   SxNeighbors neighbors;
   // note: for rcut=13, 1 bohr wide Gaussians behave like point multipoles
   // within 1e-16
   // we assume that multipole Gaussians are always below 1 bohr wide
   double rcut = 13.0 * rSoft;

   SxAtomicStructure GradU = structure.getNewStr ().set(Coord(0.0,0.0,0.0));

   const SxYlm::SxClebschTable &cg = pot.clebschGordan;
   SxGaussIJ gaussIJ;

   // --- loop over atoms (index i => R in Bloechl notation)
   for (int is = 0; is < structure.getNSpecies (); is++)   {
      int lmaxI = pot.lMaxRho(is);
      for (int ia = 0; ia < structure.getNAtoms(is); ia++)   {
         // compute neighbors within rcut
         neighbors.compute (grid, structure, structure.getAtom(is,ia),
                            rcut, SxNeighbors::StoreRel);
         // --- loop over neighbors j (R' in Bloechl notation)
         for (int js = 0; js < structure.getNSpecies (); js++)   {
            int lmaxJ = pot.lMaxRho(js);
            double a2Hard = sqr(pot.rc(is)) + sqr(pot.rc(js));
            double a2Soft = 2.0 * sqr(rSoft);
            for (int ja = 0; ja < neighbors.relPositions.getNAtoms(js); ja++)
            {
               // i => R
               // j => R'
               // r = R' - R
               Coord r = neighbors.relPositions.getAtom(js, ja);
               int jaGlobal = neighbors.relPositions.getIAtom(js,ja);
               int jRefAtom = neighbors.idx(jaGlobal)
                            - structure.atomInfo->offset(js);
               gaussIJ.setDelta (r, a2Hard, a2Soft, lmaxI + lmaxJ + 1);
               gaussIJ.compute (lmaxI, lmaxJ, cg, 1); // prepare forces

               Coord sum(0., 0., 0.);
               SX_LOOP2 (lm2, lm1) {
                  sum += Qrl(is)(ia)(lm1)
                       * Qrl(js)(jRefAtom)(lm2)
                       * gaussIJ.getForce ((int)lm1, (int)lm2, cg);
               }
               sum *= 0.5;
               GradU.ref(is,ia) += sum;
               GradU.ref(js,jRefAtom) -=sum;
            }
         }
      }
   }
   return GradU;
}

PsiG SxPAWHamiltonian::computeRhoHatG (int is, int ia, SxDiracMat<Double> &YlmGl, double rc2)
{  
   double rNorm = FOUR_PI / sqrt(structure.cell.volume);
   SX_CHECK_NUM (rNorm);

   SX_CHECK(potPtr);
  
   const SxPAWPot &pawpot = *potPtr; 
   const SxGBasis &gBasis = YlmGl.getBasis<SxGBasis> ();

   PsiG nHatG (gBasis);
   nHatG.set(0.0);

   int lMax = pawpot.lMaxRho (is);
   SxDiracVec<Double> gauss = exp ( -(0.25 * rc2) * gBasis.g2);
  
   double dFac = 1.;
   PsiG T = gBasis.getPhaseFactors (is,ia);
   for (int l = 0; l <= lMax; l++, dFac *= 2 * l + 1)   {
      for (int m = -l; m <= l; m++)   {
         int lm = SxYlm::combineLm(l,m);
         double N = rNorm * SxYlm::getYlmNormFactor(l,m) / dFac;
         double weight = N * Qrl(is)(ia)(lm);
         // i^l factor
         SxComplex16 il = 1.;
                       if ((l & 3) == 1) il = I; 
                  else if ((l & 3) == 2) il = -1.;
                  else if ((l & 3) == 3) il = -I;
         PsiG grl = (YlmGl.colRef(lm) * gauss) * T;
         nHatG.plus_assign_ax(il * weight,grl);
      }
   }

   return nHatG;
}

SxArray<SxRadMat> SxPAWHamiltonian::computeUnSymmGradDij (const SxFermi &fermi,
                                                          const int iSpin,
                                                          const int ik,
                                                          SxDiracVec<Complex16> &P,
                                                          SxArray<SxDiracMat<Complex16> > &gradP)
{
   SX_CLOCK (Timer::GradDij);
   SX_CHECK (fermi.kpPtr);
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;
   int nStates = fermi.getNStates ();
   SxArray<SxRadMat> res(3);
   for (int iDir = 0; iDir < 3; iDir++)   {
      res(iDir).resize (structure.atomInfo, pot, 1);
      res(iDir).set (0.0);
   }

   double weight = fermi.kpPtr->weights(ik);
   for (int iDir = 0; iDir < 3; iDir++)   {
      int offset = 0;
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         int nProjLocal = pot.getNProj (is);

         // --- set up i^l 
         SxVector<Complex16> il(nProjLocal); // i^(-l)  (:ipl)
         for (int ipt = 0, ipl = 0; ipt < pot.lPhi(is).getSize (); ++ipt)  {
            int l = pot.lPhi(is)(ipt);
            SxComplex16 c = (l & 1) ? I : SxComplex16(1.);
            if (!(l & 2)) c=-c;
            for (int m = -l; m <=l; ++m, ++ipl)  {
               il(ipl) = c;
            }
         }

         for (int ia = 0; ia < structure.getNAtoms(is); ++ia, offset+=nProjLocal)  {
            SxMatrix<Double> &D = res(iDir)(0,is,ia);
            for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
               for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
                  SxComplex16 phase = il(ipl) * il(jpl).conj ();
                  for (int iState = 0; iState < nStates; ++iState)  {
                     // ref. 1, Eq. (51)
                     // + projector i^l phase-correction (Ref. 2)
                     // + use k-dependent weights
                     // + exploit k/-k symmetry (=> .re)

                     // Somehow a -1 is in the derivative of Dij  
                     D(ipl, jpl) += -1.0 * weight * fermi.focc(iState,iSpin,ik)
                        * (( gradP(iDir)(offset+ipl,iState).conj() * P(offset+jpl,iState)
                                 +   P(offset+ipl,iState).conj() * gradP(iDir)(offset+jpl,iState)) * phase
                          ).re;
                  }
               }
            }
         }
      }
   }
   return res;
}

SxComplex16
SxPAWHamiltonian::pawGradDot (int iSpecies, int iAtom, int i, int j,
                              const SxDiracMat<Complex16> &P,
                              const SxDiracMat<Complex16> &gradP) const
{
   // SX_CLOCK (Timer::pawGradDot); // don't use: clock takes 50% of time
   SX_CHECK (potPtr);
   const SxArray<int> &lPhi = potPtr->lPhi(iSpecies);
   const SxArray<int> &offset = potPtr->offset(iSpecies);
   const SxDiracMat<Double> &dS = potPtr->deltaS(iSpecies);
   SxComplex16 res = 0.;
   int offsetA = 0;
   for (int js = 0; js < iSpecies; ++js)  {
      offsetA += int(projBasis->refOrbMap(js).getSize ())
                 * structure.getNAtoms (js);
   }
   offsetA += int(projBasis->refOrbMap(iSpecies).getSize ()) * iAtom;
   for (int ipt = 0; ipt < lPhi.getSize (); ++ipt)   {
      int l = lPhi(ipt), nm = 2 * l + 1;
      for(int jpt = 0; jpt < lPhi.getSize (); jpt++)   {
         if (lPhi(jpt) != l) continue;
         // find position in global projector counting 
         int ipg = offset(ipt) + offsetA;
         int jpg = offset(jpt) + offsetA;
         for (int m = 0; m < nm; ++m, ++ipg, ++jpg)  {
            res -= (  gradP(ipg,i).conj() *     P(jpg,j)
                    +     P(ipg,i).conj() * gradP(jpg,j))     
                 * dS(ipt,jpt);
         }
      }
   }
   return res;
}

SxArray<SxRadMat> SxPAWHamiltonian::symmetrizeGradDij (SxArray<SxRadMat> &GDij)
{
   SX_CLOCK (Timer::GDijSym);
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;

   SX_CHECK (structure.cell.symGroupPtr);
   const SxSymGroup & syms = *structure.cell.symGroupPtr;
   int nSyms = syms.getNSymmorphic ();
   SxVector<Int> symAtomIdx(nSyms);
   SxGrid grid (structure,10);
   int nSpin = GDij(0).getNSpin ();

   SxArray<SxRadMat> res(3);
   for (int iDir = 0; iDir < 3; iDir++)   {
      res(iDir).resize (structure.atomInfo, pot, nSpin);
      res(iDir).set (0.0);
   }


   //cout << "GradDij symmetrization has to be checked" << endl;
   // Concept like SxRadMat symmetrization

   for (int iSpin = 0; iSpin < nSpin; iSpin++)   { 
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         int nProjTypes = pot.getNProjType (is);
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  { 
            // --- find equivalent atoms
            for (int iSym = 0; iSym < nSyms; ++iSym)  {
               symAtomIdx(iSym) = structure.find (syms.getSymmorphic(iSym) 
                     ^ structure.getAtom(is,ia), grid)
                     -structure.atomInfo->offset(is);
            }
            // check if symmetrization has been performed already
            if (symAtomIdx.minval () < ia) continue;

            // ---loop over partial wave types
            for (int ipt = 0; ipt < nProjTypes; ++ipt)  {
               int l1 = pot.lPhi(is)(ipt), nm1 = 2 * l1 + 1;
               int offset1 = pot.offset(is)(ipt);
               for (int ipt2 = 0; ipt2 < nProjTypes; ++ipt2)  {
                  int l2 = pot.lPhi(is)(ipt2), nm2 = 2 * l2 + 1;
                  int offset2 = pot.offset(is)(ipt2);
                  // --- construct symmetrized Dij (type1,type2)
                  SxArray<SxMatrix<Double> > sym1 (3);
                  SxMatrix<Double> rotGDijDir (nm1,nm2); 
                  for(int iDir = 0; iDir < 3; iDir++)   {
                     sym1(iDir).reformat(nm1,nm2);
                     sym1(iDir).set (0.0);
                  }
                  for (int iSym = 0; iSym < nSyms; ++iSym)  {
                     SymMat S = syms.getSymmorphic(iSym);
                     for (int iDir = 0; iDir < 3 ; iDir++)   {
                        // --- collect rotgradDij components
                        SxMatrix<Double> &GDijDir = GDij(iDir)(iSpin,is,symAtomIdx(iSym));
                        const SxMatrix<Double> &Dl1 = ylmRot(iSym)(l1);
                        const SxMatrix<Double> &Dl2 = ylmRot(iSym)(l2);
                        for (int m1 = 0; m1 < nm1; ++m1)   {
                           for (int m2 = 0; m2 < nm2; ++m2)   {
                              rotGDijDir(m1,m2) = GDijDir(offset1 + m1, offset2 + m2);
                           }
                        }
                        // debug
                        /*
                        if (ipt == 0 && ipt2 == 0 && iDir == 2)  {
                           cout << iSym << '\t' << rotGDijDir;
                           cout << '\t' << (S(2,2) * rotGDijDir) << endl;
                        }
                        */
                        rotGDijDir = Dl1.transpose () ^ rotGDijDir ^ Dl2;
                        for(int jDir = 0; jDir < 3; jDir++)   {
                           sym1(jDir) += S(iDir, jDir) * rotGDijDir;
                           //sym1(jDir) += S(jDir, iDir) * rotGDijDir;
                        }
                     }
                  }

                  for (int iDir = 0; iDir < 3 ; iDir++)   {
                     sym1(iDir) /= double(nSyms);
                  }

                  // --- distribute symmetrized gradient
                  for (int iSym = 0; iSym < nSyms; ++iSym)  {
                     SymMat S = syms.getSymmorphic(iSym);
                     for (int iDir = 0; iDir < 3 ; iDir++)   {
                        // --- collect rotgradDij component
                        const SxMatrix<Double> &Dl1 = ylmRot(iSym)(l1);
                        const SxMatrix<Double> &Dl2 = ylmRot(iSym)(l2);
                        rotGDijDir.set(0.0);
                        for(int jDir = 0; jDir < 3; jDir++)   {
                           rotGDijDir += S(iDir,jDir) * sym1(jDir);
                           //rotGDijDir += S(jDir,iDir) * sym1(jDir);

                        }
                        rotGDijDir = Dl1 ^ rotGDijDir ^ Dl2.transpose ();

                        // distribute
                        for (int m1 = 0; m1 < nm1; ++m1)   {
                           for (int m2 = 0; m2 < nm2; ++m2)   {
                              res(iDir)(iSpin,is,symAtomIdx(iSym))(offset1 + m1, offset2 + m2) = rotGDijDir(m1,m2);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }         

   return res;

}

SxConstPtr<SxOverlapBase> SxPAWHamiltonian::getS () const
{
   return SxPtr<SxPAWOverlap>::create (pBasis, potPtr);
}

SxMeshR SxPAWHamiltonian::computeCorePS () const
{
   SX_CHECK (rPtr);
   const SxGBasis &gBasis = rPtr->getGBasis ();
   PsiG coresG(gBasis);
   coresG.set (0.);
   for (int is = 0; is < structure.getNSpecies (); ++is)
      coresG += rhoCorePSG(is) * gBasis.structureFactors(is);
   VALIDATE_VECTOR (coresG);
   SxMeshR res = coresG.to (*rPtr);
   res = rPtr->symmetrize (res);
   return res;
}

void SxPAWHamiltonian::computeRho (const Focc &focc, const SxPsiSet &psi)
{
   SX_CHECK (rPtr);
   const SxPWSet &waves = dynamic_cast<const SxPWSet &> (psi);
   // compute 1-center density matrix
   pawRho.Dij = computeDij (waves, focc);

   // --- compute pseudo-density
   int nSpin = waves.getNSpin ();
   if (pawRho.pwRho.rhoR.getSize () == 0)  {
      // initialize pwRho
      pawRho.pwRho = SxRho (*rPtr, nSpin, 0. ); // 0. = no renormalization
      // set paw potential
      pawRho.potPtr = potPtr;
   }
   pawRho.pwRho.computeRho (focc, waves);

   // --- add pseudo-core densities
   SxMeshR coresR = computeCorePS ();
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      pawRho.pwRho(iSpin).plus_assign_ax (1. / nSpin, coresR);
   //pawRho.psCore = coresR;

   if (hubbardU && hubbardU->moSite.getSize () > 0)  {
      // --- compute MO sites' density matrix
      SX_LOOP3(i,iSpin,iSite) (*pawRho.blockAO(i))(iSpin,iSite).set (0.);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         for (int ik = 0; ik < waves.getNk (); ++ik)  {
            if (SxLoopMPI::myWork (ik))  {
               double weight = waves.getGkBasis ().weights(ik);
               SX_LOOP(i) {
                  hubbardU->moSite(i)->addToRho (pawRho.blockAO(i).getPtr (),
                                                 waves(iSpin,ik),
                                                 weight, focc(iSpin,ik));
               }
            }
         }
      }
      // sum over MPI nodes and symmetrize in the end
      SX_LOOP(i)  {
         pawRho.blockAO(i)->sumMPI ();
         hubbardU->moSite(i)->symmetrize(pawRho.blockAO(i).getPtr (),
                                         structure, ylmRot);
      }
   }

}

