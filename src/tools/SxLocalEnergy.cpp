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
#include <SxProjector.h>
#include <SxAOBasis.h>
#include <SxPAWOverlap.h>
#include <SxPWOverlap.h>
#include <SxFermi.h>
#include <SxPW.h>
#include <SxHamSolver.h>
#include <SxPWHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxXCFunctional.h>
#include <SxRadialAtom.h>
#include <SxNeighbors.h>
#include <SxGaussIJ.h>


typedef SxArray<SxArray<SxDiracMat<Complex16> > > SxExpCoeffs;
typedef SxNArray<double, 3> SxClebschTable;

namespace Timer {
   enum LocalEnergyTimer {
      Startup, AOBasis, CalcCoeff, ComputeLocalEnergies, FinalPrint, 
      ExpansionCoeffs, ProjectWaves, RhoProjQ, ComputeEtot, 
      ComputeEkinG, ComputeEkinR,ComputeExc, ComputeEbar, 
      ComputeEHartree};
}

SX_REGISTER_TIMERS(Timer::LocalEnergyTimer)
{
   using namespace Timer;
   regTimer (Startup, "startup");
   regTimer (AOBasis, "set up AOBasis");
   regTimer (CalcCoeff, "calculate C");
   regTimer (ComputeLocalEnergies, "compute ELoc");
   regTimer (FinalPrint, "final result print");
   regTimer (ExpansionCoeffs, "C = S^-1 <mu|Psi>");
   regTimer (ProjectWaves, "P|Psi>");
   regTimer (RhoProjQ, "rhoProjQ");
   regTimer (ComputeEtot, "Etot");
   regTimer (ComputeEkinG, "EkinG");
   regTimer (ComputeEkinR, "EkinR");
   regTimer (ComputeExc, "Exc");
   regTimer (ComputeEbar, "Ebar");
   regTimer (ComputeEHartree, "EHartree");
}

bool keepWavesOnDisk = false;
bool saveMemory = false;
bool useSymmetries = false;
bool langeFreysoldt = false;
double norm = 0.0;

SxPtr<SxDipoleCorrZ> dipoleCorr;

void matPrint (SxMatrix<Double> &mat)
{
   ssize_t nRows = mat.nRows ();
   ssize_t nCols = mat.nCols ();

   for (ssize_t iRow = 0; iRow < nRows; iRow++)  {
      sxprintf ("[");
      for (ssize_t iCol = 0; iCol < nCols; iCol++)  {
         if (fabs(mat(iRow,iCol)> 1e-6))
            sxprintf ("%8.6e\t", mat(iRow,iCol));
         else
             sxprintf ("0.0\t\t");
         if ((iCol+1)%8 == 0) sxprintf ("\n");
      }
      sxprintf ("]\n");
   }
}

SxDiracVec<Double> createMask (SxPAWRho &projRho, 
                               SxPAWRho &fullRho, 
                               const SxArray<SxVector<Int> > &atomicSpace, 
                               const SxVector<Int> &map)
{
   ssize_t dim = projRho.pwRho(0).getSize ();
   int nSpin = projRho.getNSpin ();
   
   SxDiracVec<Double> result (projRho.pwRho(0).getSize ());
   result.set(0.0);

   SxDiracVec<Double> rhoProjR = (nSpin == 1) ? 1.0 * projRho.pwRho.rhoR(0)
                               : projRho.pwRho.rhoR(0) + projRho.pwRho.rhoR(1);
   SxDiracVec<Double> rhoFullR = (nSpin == 1) ? 1.0 * fullRho.pwRho.rhoR(0)
                               : fullRho.pwRho.rhoR(0) + fullRho.pwRho.rhoR(1);

   for(ssize_t i = 0; i < dim; i++)  {
      double weight = 0.0;
      if (fabs (rhoFullR(i)) < 1e-6)  {
         ssize_t nPts = atomicSpace(i).getSize ();
         double value = 1.0/double(nPts);
         for (ssize_t j = 0; j < nPts; j++)  {
            for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
              if ( atomicSpace(i)(j) == map(iAtomX) )  {
                 weight += value;
              }
            }
         }
      } else weight = rhoProjR(i) / rhoFullR(i);
      result(i) = weight;
   }

   return result;
}

SxArray<SxVector<Int> > atomicSpaceRegion (const SxRBasis *rPtr, const SxAtomicStructure &structure) 
{
   const SxRBasis &R = *rPtr;
   SxVector3<Int> dims = R.getMesh ();
   SxMesh3D mesh (dims);
   int nDim = dims.product ();

   // grid for neighbor search
   SxGrid grid (structure, 10);

   SxArray<SxVector<Int> > result(nDim);

   for (int i = 0; i < nDim; i++)  {
      // get Coord
      Coord abs = 
         structure.cell.relToCar (mesh.getMeshVec(i,SxMesh3D::Positive)/Coord(mesh)); 

      // find Neighbours
      SxNeighbors neighbors;
      double rCut = 6.0;
      ssize_t nNeighbors = 0;
      do  {
         int mode = SxNeighbors::StoreIdx |
            SxNeighbors::StoreDistSqr |
            SxNeighbors::IncludeZeroDistance;
         neighbors.compute(grid, structure, abs, rCut, mode);
         rCut *= 2.0;
         nNeighbors = neighbors.idx.getSize ();
      } while (nNeighbors == 0);

      int minIdx = -1;
      neighbors.distSqr.minval(&minIdx);
      SxStack<int> atomIndices;
      for (ssize_t iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++)  {
         if (fabs(neighbors.distSqr(iNeighbor) - neighbors.distSqr(minIdx)) < 1e-4) 
           atomIndices << neighbors.idx(iNeighbor);
      } 
      result(i) = SxVector<Int>(atomIndices);
   }

#ifndef NDEBUG
   for (int i = 0; i < nDim; i++)  {
      if (result(i).getSize () < 1) {
         cout << "Spacepoint has no atomic neighbor" << endl;
         SX_EXIT
      }
   }
#endif

   return result;
}

double getSpaceNorm (const SxPW &proj, const SxPW &ref, SxPtr<SxPAWOverlap> SPtr)
{
   double result = 0.0;

   int nStates = ref.getNStates ();
   int nk      = ref.getNk();
   int nSpin   = ref.getNSpin ();

   const SxGkBasis &gk = *ref.getGkBasisPtr ();

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         SxDiracMat<Complex16> SRefSet = SPtr->apply(ref(iSpin, ik));
         for (int iState = 0; iState < nStates; iState++)  {
            result += gk.weights(ik) 
               * dot(proj(iState,iSpin,ik), SRefSet.colRef(iState)).re;
         }
      }
   }

   result = SxLoopMPI::sum(result);

   return result;
}

double getSpaceNorm (const SxPW &proj, 
                     const SxPW &ref, 
                     SxPtr<SxPAWOverlap> SPtr, 
                     SxDiracVec<Double> &filter)
{
   double result = 0.0;

   int nStates = ref.getNStates ();
   int nk      = ref.getNk();
   int nSpin   = ref.getNSpin ();

   const SxGkBasis &gk = *ref.getGkBasisPtr ();
   const SxRBasis &R = filter.getBasis<SxRBasis> ();
   double dOmega = R.dOmega;

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         SxDiracMat<Complex16> SRefSet = SPtr->apply(ref(iSpin, ik));
         for (int iState = 0; iState < nStates; iState++)  {
            PsiG psiGProj = (gk(ik) | proj(iState,iSpin,ik));
            PsiG psiGRef  = (gk(ik) | SRefSet.colRef(iState));
            PsiR psiRProj = filter * (R | psiGProj );
            PsiR psiRRef  = (R | psiGRef);
            result += gk.weights(ik) 
               * (psiRProj.conj() * psiRRef).sum ().re * dOmega;
         }
      }
   }

   SxLoopMPI::sum (result);

   return result;
}

SxPW wavesDiff(const SxPW &waves1,const SxPW &waves2)
{

   SX_CHECK(waves1.getNStates () == waves2.getNStates ());
   SX_CHECK(waves1.getNSpin () == waves2.getNSpin ());
   SX_CHECK(waves1.getNk () == waves2.getNk ());

   int nk = waves1.getNk ();
   int nSpin = waves1.getNSpin ();
   int nStates = waves1.getNStates ();

   SxPtr<SxGkBasis> gkBasisPtr = waves1.getGkBasisPtr ();

   SxString tmpDir = "";
   if (keepWavesOnDisk) tmpDir = ".";

   SxPW result(nStates, nSpin, gkBasisPtr, tmpDir);

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         result(iSpin,ik) = waves1(iSpin,ik) - waves2(iSpin,ik);
      }
   }
   
   result.setGkBasisPtr (gkBasisPtr);

   return result;
}

SxExpCoeffs computeExpCoeffs (const SxAOBasis &aoBasis, const SxPW &waves)
{
SX_CLOCK(Timer::ExpansionCoeffs);

   int nk      = waves.getNk();
   int nSpin   = waves.getNSpin ();
   int nStates = waves.getNStates ();
   int nOrb = aoBasis.getNOrb ();

   SxPtr<SxGkBasis> gkBasisPtr = waves.getGkBasisPtr ();

   SxExpCoeffs result (nk);
   for (int ik = 0; ik < nk; ik++)  {
      result(ik).resize(nSpin);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         result(ik)(iSpin).reformat(nOrb,nStates);
         result(ik)(iSpin).set(0.0);
      }
   }

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      const SxDiracMat<Complex16> &Sinv = aoBasis.getInverseOverlap(ik);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         const PsiG psiG = waves(iSpin,ik);
         SxDiracMat<Complex16> muPsiG = (aoBasis | psiG);
         result(ik)(iSpin) = (Sinv ^ muPsiG);
      }
   }

   for (int ik = 0; ik < nk; ik++)  {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         SxLoopMPI::sum (result(ik)(iSpin));
      }
   }
   
   return result;
}

SxPW projectWaves (const SxAOBasis &aoBasis, 
                   const SxPW &waves,
                   const SxExpCoeffs C, 
                   const SxAtomicStructure &structure, 
                   const SxVector<Int> &map)
{
   SX_CLOCK(Timer::ProjectWaves);

   int nStates = waves.getNStates ();
   int nk      = waves.getNk();
   int nSpin   = waves.getNSpin ();
   int nOrb = aoBasis.getNOrb ();

   SxPtr<SxGkBasis> gkBasisPtr = waves.getGkBasisPtr ();

   SxString tmpDir = "";
   if (keepWavesOnDisk) tmpDir = ".";

   SxPW projWaves (nStates, nSpin, gkBasisPtr,tmpDir);
   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         projWaves(iSpin,ik).set(0.);
         for (int iOrb = 0; iOrb < nOrb; iOrb++)  {
            for (int iAtomX = 0; 
                     iAtomX < map.getSize (); 
                     iAtomX++) {
               int iAtom;
               int iSpecies = structure.getISpecies(map(iAtomX), &iAtom);
               if ( iSpecies == aoBasis.orbitalMap(iOrb).is 
                    && iAtom == aoBasis.orbitalMap(iOrb).ia)  {
                  PsiG nu = aoBasis.getAOinG (ik, iOrb);
                  for (int iState = 0; iState < nStates; iState++)  {
                     projWaves(iState, iSpin, ik)
                        .plus_assign_ax (C(ik)(iSpin)(iOrb,iState),nu);
                  }
               }
            }
         }
      }
   }

   projWaves.setGkBasisPtr (gkBasisPtr);
   
   return projWaves;
}

SxPAWRho computeRhoProjS (const SxPAWRho &projRhoQ, 
                          const SxPAWRho &rho, 
                          const SxDiracVec<Double> &filter)
{
   SxPAWRho result;
   result = rho - projRhoQ;

   int nSpin = rho.pwRho.getNSpin ();

   for (int iSpin = 0; iSpin < nSpin; iSpin++)
      result.pwRho(iSpin) *= filter;
   result.Dij.set(0.0);

   return result;
}

SxPAWRho computeRhoProjQ (const Focc &focc, 
                     const SxPW &projPsi, 
                     const SxPW &fullPsi,
                     const SxPAWRho refRho,
                     const SxRBasis *rPtr,
                     SxPtr<SxPAWPot> pawPotPtr, 
                     const SxAtomicStructure &structure, 
                     const SxVector<Int> &map)
{
   SX_CLOCK(Timer::RhoProjQ);

   int nk      = fullPsi.getNk ();
   int nSpin   = fullPsi.getNSpin ();
   int nStates = fullPsi.getNStates ();

   const SxRBasis &R = *rPtr;
   const SxGBasis &G = R.getGBasis ();
   const SxGkBasis &gk = fullPsi.getGkBasis ();

   SxPAWRho result(pawPotPtr);

   // compute 1-center density matrix Dij
   result.Dij.resize (structure.atomInfo, *pawPotPtr, nSpin);
   result.Dij.set(0.0);
   SX_LOOP3(iSpin,iSpecies,iAtom)  {
      int iAtomTl = structure.getIAtom(iSpecies,iAtom);
      for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  
         if ( iAtomTl == map(iAtomX) )
            result.Dij(iSpin,iSpecies,iAtom) 
                   <<= refRho.Dij(iSpin,iSpecies,iAtom);
   }

   // initialize pwRho
   result.pwRho = SxRho (R, nSpin, 0. ); // 0. = no renormalization

   // --- initialize rhoR
   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      result.pwRho(iSpin).resize(R.getMeshSize ());
      result.pwRho(iSpin).set(0.0);
      result.pwRho(iSpin).setBasis(rPtr);
   }

   // --- compute rhoR
   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int iState = 0; iState < nStates; iState++)  {
            double f = gk.weights(ik) * focc(iState,iSpin,ik);
            if (f > 1e-12)  {
               PsiR psiFullR = (R | gk(ik)) * (gk(ik) | fullPsi(iState,iSpin,ik));
               PsiR psiProjR = (R | gk(ik)) * (gk(ik) | projPsi(iState,iSpin,ik));
               result.pwRho(iSpin)
                  .plus_assign_ax (f, (psiProjR.conj() * psiFullR).real ());
            }
         }
      }
   }

   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      SxLoopMPI::sum (result.pwRho(iSpin));
   }

   // --- add pseudo-core densities
   int nSpecies = structure.getNSpecies ();
   SxArray<PsiG> rhoCorePSG(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      rhoCorePSG(is) = (G | pawPotPtr->rhoCorePS(is));
   }
   PsiG coresG(G);
   coresG.set (0.);
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
      int ia;
      int is = structure.getISpecies(map(iAtomX), &ia);
      coresG += rhoCorePSG(is) * G.getPhaseFactors (is,ia);
   }
   SxMeshR coresR = coresG.to (R);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      result.pwRho(iSpin).plus_assign_ax (1. / nSpin, coresR);
   }

   // --- symmetrize
   if (useSymmetries) {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  { 
         result.pwRho(iSpin) = R.symmetrize (result.pwRho.rhoR(iSpin));
      }
   }

   return result;
}

double computeU (const SxArray<SxArray<SxVector<Double> > > &Qrl,
                 const SxPAWPot &pawPot, 
                 const SxAtomicStructure &structure, 
                 const SxVector<Int> &map)
{
   // --- for neighbor search
   SxGrid grid (structure, 10);
   int neighborParams = SxNeighbors::StoreRel
                      | SxNeighbors::IncludeZeroDistance;
   // note: for rcut=13, 1 bohr wide Gaussians behave like point multipoles
   // within 1e-16
   // we assume that multipole Gaussians are always below 1 bohr wide
   double rSoft = SxPAWHamiltonian::rSoft;
   double rcut = 13. * rSoft;

   const SxYlm::SxClebschTable &cg = pawPot.clebschGordan;
   SxGaussIJ gaussIJ;

   // set up full U matrix
   int nAtoms = structure.getNAtoms ();
   SxArray<SxNeighbors> neighbors (nAtoms);
   int maxNeighbors = 0;
   // --- loop over atoms (index i => R in Bloechl notation)
   for (int iAtomI = 0; iAtomI < nAtoms; iAtomI++)  {
      int ia;
      int is = structure.getISpecies(iAtomI, &ia);
      // compute neighbors within rcut
      neighbors(iAtomI).compute (grid, structure, structure.getAtom(is,ia),
                         rcut, neighborParams);
      int nNeighbors = neighbors(iAtomI).relPositions.getNAtoms();
      if (maxNeighbors < nNeighbors) maxNeighbors = nNeighbors;
   }
   
   SxMatrix<Double> UMat(nAtoms, maxNeighbors);
   UMat.set(0.0);
   for (int iAtomI = 0; iAtomI < nAtoms; iAtomI++)  {
      int ia;
      int is = structure.getISpecies(iAtomI, &ia);
      int lmaxI = pawPot.lMaxRho(is);
      int nNeighbors = neighbors(iAtomI).relPositions.getNAtoms();
      // --- loop over neighbors j (R' in Bloechl notation)
      for (int iAtomJ = 0; iAtomJ < nNeighbors; iAtomJ++)  { 
         int ja;
         int js = neighbors(iAtomI).relPositions.getISpecies(iAtomJ, &ja);
         int lmaxJ = pawPot.lMaxRho(js);
         // i => R
         // j => R'
         // r = R'-R
         Coord r = neighbors(iAtomI).relPositions.getAtom(js, ja);
         // Ref. 1 Eq. (28): kernel
         SxMatrix<Double> delta(sqr(lmaxI+1),sqr(lmaxJ+1));
         double a2Hard = sqr(pawPot.rc(is)) + sqr(pawPot.rc(js));
         double a2Soft = sqr(rSoft) + sqr(rSoft);
         double a2SoftHard = sqr(rSoft) + sqr(pawPot.rc(js));
         double a2HardSoft = sqr(pawPot.rc(is)) + sqr(rSoft);
         

         if (langeFreysoldt)  {
            // --- Lange-Freysoldt variant
            gaussIJ.setDelta (r, a2Hard, a2HardSoft, lmaxI + lmaxJ);
            gaussIJ.compute (lmaxI, lmaxJ, cg, &delta);
            gaussIJ.setDelta (r, a2SoftHard, a2Soft, lmaxI + lmaxJ);
            delta -= gaussIJ.compute (lmaxI, lmaxJ, cg);
         } else {
            // --- Bloechl variant
            gaussIJ.setDelta (r, a2Hard, a2Soft, lmaxI + lmaxJ);
            gaussIJ.compute (lmaxI, lmaxJ, cg, &delta);
         }
         
         // sum over moments of j=R'
         int jaGlobal = neighbors(iAtomI).relPositions.getIAtom(js,ja);
         int jRefAtom = neighbors(iAtomI).idx(jaGlobal)
            - structure.atomInfo->offset(js);
         SxVector<Double> sumJ = delta ^ Qrl(js)(jRefAtom);

         // Ref. 1 Eq. (28)
         UMat(iAtomI,iAtomJ) += 0.5 * (Qrl(is)(ia) * sumJ).sum ();
      }
   }
   //matPrint(UMat);

   double U = 0.;
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
      U += UMat.row(map(iAtomX)).sum ();
   }

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf("U = %18.12f \n", U / norm);
      cout << SX_SEPARATOR;
   }

   return U;
}

double computeEHartree (const SxPAWRho &densityAtom,
                    const SxPAWRho &densityTot,
                    const SxPAWPot &pawPot, 
                    const SxAtomicStructure &structure, 
                    const SxVector<Int> &map,
                    double charge)
{
   SX_CLOCK(Timer::ComputeEHartree);
   int nSpecies = structure.getNSpecies ();
   int nSpin = densityAtom.getNSpin ();
   
   const RhoR &rhoAtomR = densityAtom.pwRho.rhoR;
   const RhoR &rhoTotR = densityTot.pwRho.rhoR;
   const SxRBasis &rBasis = rhoAtomR(0).getBasis<SxRBasis> ();
   const SxGBasis &gBasis = rBasis.getGBasis ();
   const SxRadBasis &rad = pawPot.getRadBasis ();

   double omega = structure.cell.volume;

   // get YlmGl
   int lMaxAll =  pawPot.lMaxRho.maxval ();
   SxDiracMat<Double> YlmGl(gBasis.ng, sqr(lMaxAll + 1));
   YlmGl.setBasis (gBasis);
   for (int lm = 0, l = 0; l <= lMaxAll; ++l)  {
      SxDiracVec<Double> gl = pow (gBasis.g2, 0.5 * l);
      for (int m = -l; m <= l; ++m,++lm)
         YlmGl.colRef(lm) <<= gBasis.getYlm (l,m) * gl;
   }

   // get gRl and YlmGl
   SxArray<SxArray<SxDiracVec<Double> > > gRl (nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = pawPot.lMaxRho(is);
      gRl(is).resize(lMax+1);
      for (int l = 0; l <= lMax; l++)  {
         gRl(is)(l) = pawPot.getGrl (is, l);
      }
   }

   // setup Qrl
   SxArray<SxArray<SxVector<Double> > > Qrl (nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = pawPot.lMaxRho(is);
      Qrl(is).resize (structure.getNAtoms (is));
      for (int ia = 0; ia < structure.getNAtoms (is); ia++)  {
         SxVector<Double> &Q = Qrl(is)(ia);
         SxRadialMesh rhoRadPS, rhoRadAE;
         rhoRadPS.resize(rad,is,lMax);
         rhoRadAE.resize(rad,is,lMax);
         rhoRadPS.set(0.0);
         rhoRadAE.set(0.0);
         for (int iSpin = 0; iSpin <  nSpin; ++iSpin)  {
            rhoRadPS += 
               pawPot.computeRhoPS (densityTot.Dij(iSpin,is,ia),is,nSpin);
            rhoRadAE += 
               pawPot.computeRhoAE (densityTot.Dij(iSpin,is,ia),is,nSpin);
         }
         SxRadialMesh rhoPAWcorr = rhoRadAE - rhoRadPS;
         Q.resize(sqr(lMax+1));
         for (int l = 0; l <= rhoPAWcorr.lmax; ++l)  {
            SxDiracVec<Double> rl = pow(rad.radFunc(is),double(l+3));
            for (int m = -l; m <= l; ++m)  {
               int lm = SxYlm::combineLm(l,m);
               Q(lm) = (rhoPAWcorr(l,m) * rl).integrate (rad.logDr(is));
               if (l == 0) Q(lm) -= SQRT_1_4PI * pawPot.nuclearCharge(is);
            }
         }
      }
   }
   
   // setup nHatG and nHatPG
   PsiG nHatTotG(gBasis), nHatpTotG(gBasis), 
        nHatAtomG(gBasis), nHatpAtomG(gBasis);
   nHatTotG.set(0.0); nHatpTotG.set(0.0);
   nHatAtomG.set(0.0); nHatpAtomG.set(0.0);
   double rNorm = FOUR_PI / sqrt(omega);
   double rSoft = SxPAWHamiltonian::rSoft;
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = pawPot.lMaxRho(is);
      int nAtoms = structure.getNAtoms(is);
      double rc2 = sqr(pawPot.rc(is)), rcp2 = sqr(rSoft);
      SxDiracVec<Double> gaussHard = exp ( -(0.25 * rc2 ) * gBasis.g2);
      SxDiracVec<Double> gaussSoft = exp ( -(0.25 * rcp2) * gBasis.g2);
      for (int ia = 0; ia < nAtoms; ia++)  {
         PsiG T = gBasis.getPhaseFactors (is,ia);
         bool isChoice = false;
         int iAtomTl = structure.getIAtom(is,ia);
         for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  
            if ( iAtomTl == map(iAtomX) ) isChoice = true;
         double prefac = 1.0;
         for (int l = 0; l <= lMax; ++l)  {
            prefac /= (2.0 * l + 1.0);
            for (int m = -l; m <= l; ++m)  {
               int lm = SxYlm::combineLm(l,m);
               double N = rNorm * SxYlm::getYlmNormFactor(l,m) * prefac;
               double weight = N * Qrl(is)(ia)(lm);
               // i^l factor
               SxComplex16 il = 1.;
               if ((l & 3) == 1) il = I;
               else if ((l & 3) == 2) il = -1.;
               else if ((l & 3) == 3) il = -I;
               PsiG grl  = (YlmGl.colRef(lm) * gaussHard) * T;
               PsiG gprl = (YlmGl.colRef(lm) * gaussSoft) * T;
               nHatTotG.plus_assign_ax(il * weight, grl);
               nHatpTotG.plus_assign_ax(il * weight, gprl);
               if (isChoice)  {
                  nHatAtomG.plus_assign_ax(il * weight, grl);
                  nHatpAtomG.plus_assign_ax(il * weight, gprl);
               }
            }
         }
      }
   }

   double eHartreeRadial = 0.0;
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
      int ia;
      int is = structure.getISpecies(map(iAtomX), &ia);
      int lMax = pawPot.lMaxRho(is); 
      SxRadialMesh rhoRadPS, rhoRadAE, v1AE, v1PS;
      rhoRadPS.resize(rad,is,lMax);
      rhoRadAE.resize(rad,is,lMax);
      v1AE.resize(rad,is,lMax);
      v1PS.resize(rad,is,lMax);
      rhoRadPS.set(0.0);
      rhoRadAE.set(0.0);
      SxRadialAtom radial(SxXC::Unknown, pawPot.aGridType(is));
      for (int iSpin = 0; iSpin <  nSpin; ++iSpin)  {
         rhoRadPS += 
            pawPot.computeRhoPS (densityTot.Dij(iSpin,is,ia),is,nSpin);
         rhoRadAE += 
            pawPot.computeRhoAE (densityTot.Dij(iSpin,is,ia),is,nSpin);
      }
      for (int l = 0; l <= lMax; l++)  {
         for (int m = -l; m <= l; m++ )  {
            int lm = SxYlm::combineLm(l,m);
            SxVector<Double> &Q = Qrl(is)(ia);
            rhoRadPS(l,m).plus_assign_ax(Q(lm), gRl(is)(l));
            v1AE(l,m) <<= radial.getHartreePotential (rhoRadAE(l,m));
            v1PS(l,m) <<= radial.getHartreePotential (rhoRadPS(l,m));
            eHartreeRadial += 0.5 * tr(rhoRadAE(l,m) * v1AE(l,m));
            eHartreeRadial -= 0.5 * tr(rhoRadPS(l,m) * v1PS(l,m));
         }
      }
      double Z = pawPot.nuclearCharge(is),
             r0 = rad.radFunc(is)(0),
             r1 = rad.radFunc(is)(1),
             r01 = r0 - r1,
             rho0 = rhoRadAE(0,0,0),
             rho1 = rhoRadAE(1,0,0),
             logDr = rad.logDr(is),
             SQRT_4PI = sqrt(FOUR_PI),
             r0cub = r0 * r0 * r0;
      eHartreeRadial -= SQRT_4PI * Z
         * (rhoRadAE(0,0)  * rad.radFunc(is).sqr ())
         .integrate (logDr);
      eHartreeRadial -= SQRT_4PI * Z * (rho0 * 0.5 * r0 * r0
            - (rho0  -rho1 ) / r01 * r0cub/6.);
   }

   const PsiG &g2 = gBasis.g2;
   SxIdx noZero (1,gBasis.ng - 1);
   
   // vHatGAtom & vHatGTot
   PsiG vHatAtomG(gBasis);
   PsiG vHatTotG(gBasis);
   double vHatAtomGZero = 0.0;
   double vHatTotGZero = 0.0;
   for (int is = 0; is < nSpecies; is++)  {
      double rc2 = sqr(pawPot.rc(is)), rcp2 = sqr(rSoft);
      int nAtoms = structure.getNAtoms(is);
      for (int ia = 0; ia < nAtoms; ia++)  {
         int iTlAtom = structure.getIAtom(is,ia);
         bool isChoice = false;
         for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)
            if (iTlAtom == map(iAtomX)) isChoice = true;
         double prefac = PI * rNorm * SxYlm::getYlmNormFactor(0,0) * (rcp2 - rc2);
         vHatTotGZero += prefac * Qrl(is)(ia)(0);
         if (isChoice) vHatAtomGZero -= prefac * Qrl(is)(ia)(0);
      }
   }
   vHatAtomG(0) = vHatAtomGZero;
   vHatTotG(0) = vHatTotGZero;
   vHatAtomG(noZero) = FOUR_PI/g2(noZero) * (nHatAtomG - nHatpAtomG)(noZero);
   vHatTotG(noZero) = FOUR_PI/g2(noZero) * (nHatTotG - nHatpTotG)(noZero);

   PsiG rhoAtomG = gBasis | ((nSpin == 1) ? 1.0 * rhoAtomR(0)
                                          : rhoAtomR(0) + rhoAtomR(1));
   PsiG rhoTotG = gBasis | ((nSpin == 1) ? 1.0 * rhoTotR(0)
                                          : rhoTotR(0) + rhoTotR(1));

   // pure Soft = <g/g^2 rhoAtom(g), g/g^2 rhoTot(g)> = <rhoAtom(g),rhoTot(g)/g^2>' ' excludes zero
   PsiG rhoHartreeSoftTot = rhoTotG + nHatpTotG;
   PsiG rhoHartreeSoftAtom = rhoAtomG + nHatpAtomG;
   PsiG vHartreeTotPS(gBasis);
   vHartreeTotPS(noZero) = FOUR_PI * rhoHartreeSoftTot(noZero) / g2(noZero);
   vHartreeTotPS(0) = 0.0; // cell is neutral
   PsiG vHartreeAtomPS(gBasis);
   vHartreeAtomPS(noZero) = FOUR_PI * rhoHartreeSoftAtom(noZero) / g2(noZero);
   // atom might be charged, is this right?
   vHartreeAtomPS(0) = FOUR_PI * charge / sqrt(omega); 
   
   /* Symmetric form
   double eHartreePS = 0.25 * dot(rhoHartreeSoftAtom,vHartreeTotPS).re
                     + 0.25 * dot(rhoHartreeSoftTot,vHartreeAtomPS).re;
   */

   // rhoAtom in Potential of rhoTot
   double eHartreePS = 0.5 * dot(rhoHartreeSoftAtom,vHartreeTotPS).re;

   SX_MPI_MASTER_ONLY { 
      cout << SX_SEPARATOR;
      sxprintf ("rhoTotG(0)            = %18.12f\n", rhoTotG(0).re);
      sxprintf ("rhoAtomG(0)           = %18.12f\n", rhoAtomG(0).re / norm);
      sxprintf ("vHatTotG(0)           = %18.12f\n", vHatTotG(0).re);
      sxprintf ("vHatAtomG(0)          = %18.12f\n", vHatAtomG(0).re / norm);
      sxprintf ("rhoHartreeSoftAtom(0) = %18.12f\n", rhoHartreeSoftAtom(0).re / norm);
      sxprintf ("rhoHartreeSoftTot(0)  = %18.12f\n", rhoHartreeSoftTot(0).re);
      sxprintf ("vHartreeTotPS(0)      = %18.12f\n", vHartreeTotPS(0).re);
      sxprintf ("vHartreeAtomPS(0)     = %18.12f\n", vHartreeAtomPS(0).re / norm);
      cout << SX_SEPARATOR;
   }
   
   double eVHat = 0.0;
   if (langeFreysoldt)  {
      // --- Lange-Freysoldt variant
      eVHat += 0.5 * dot(nHatAtomG - nHatpAtomG,vHartreeTotPS).re;
   
      eVHat += 0.5 * dot(rhoHartreeSoftAtom, vHatTotG).re;
   } else {
      // Bloechl variant
      eVHat += dot(rhoAtomG, vHatTotG).re;
   }

   double eHartreeU = 0.0;
   eHartreeU += computeU (Qrl, pawPot, structure, map);

   double eDipole = 0.0;
   SxMeshR vDipole(rBasis);
   vDipole.set (0.);
   if (dipoleCorr)  {
      SxMeshR rhoHartreeRTot = rhoHartreeSoftTot.to (rBasis);
      SxMeshR rhoHartreeRAtom = rhoHartreeSoftAtom.to (rBasis);
      dipoleCorr->update (rhoTotG.to (rBasis), rhoHartreeRTot);
      dipoleCorr->correctPotential (&vDipole);
      // correction to energy
      eDipole = 0.5 * tr(rhoHartreeRAtom * vDipole);
   }

   double eHartree = eHartreePS + eVHat + eHartreeU + eHartreeRadial + eDipole;

   if ((map.getSize() == structure.getNAtoms()))  {
      PsiG rhoHartreeHard = rhoTotG + nHatTotG;
      PsiG vHartreeHard(gBasis);
      vHartreeHard(0) = vHatTotGZero;
      vHartreeHard(noZero) = FOUR_PI * rhoHartreeHard(noZero)
         / gBasis.g2(noZero);
      RhoR vElStat(1);
      vElStat(0) = vHartreeHard.to (rBasis) + vDipole;
      SxRho(vElStat).writeRho ("vHartreeHard-H.sxb");
   }

   //double eCorr = charge * (vHatTotGZero) / sqrt(omega);

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("eHartreePS      = %18.12f H\n", eHartreePS / norm);
      sxprintf ("eVHat           = %18.12f H\n", eVHat / norm);
      sxprintf ("eHartreeU       = %18.12f H\n", eHartreeU / norm);
      sxprintf ("eHartreeRadial  = %18.12f H\n", eHartreeRadial / norm);
      if (dipoleCorr)
         sxprintf ("eDipole         = %18.12f H\n", eDipole);
      //sxprintf ("Q               = %18.12f H\n", charge / norm); 
      //sxprintf ("dV(n - n')      = %18.12f H\n", vHatTotGZero);
      //sxprintf ("Q * dV(n - n') /sqrt(Omega) = %18.12f H\n", eCorr / norm);
      sxprintf ("Omega  = %18.12f H\n", omega); 
      cout << SX_SEPARATOR;
   }
   
   //return eHartree - eCorr;
   //return eHartree - 0.5 * eCorr;
   return eHartree;
}

double computeEBar (const SxPAWRho &densityAtom, 
                    const SxPAWPot &pawPot, 
                    const SxAtomicStructure &structure, 
                    const SxVector<Int> &map)
{
   SX_CLOCK(Timer::ComputeEbar);
   int nSpecies = structure.getNSpecies ();
   int nSpin = densityAtom.getNSpin ();
   
   RhoR rhoAtom(nSpin);
   for (int iSpin = 0; iSpin < nSpin; iSpin++)
      rhoAtom(iSpin) = 1.0 * densityAtom.pwRho.rhoR(iSpin); 
   const SxRBasis &rBasis = rhoAtom(0).getBasis<SxRBasis> ();
   const SxGBasis &gBasis = rBasis.getGBasis ();
   const SxRadBasis &rad = pawPot.getRadBasis ();


   double eBarRadial = 0.0;
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
      int ia;
      int is = structure.getISpecies(map(iAtomX), &ia);
      ssize_t nr = rad.radFunc(is).getSize ();
      int lMax = pawPot.lMaxRho(is); 
      SxRadialMesh rhoRadPS(nr,lMax,0.0);
      for (int iSpin = 0; iSpin <  nSpin; ++iSpin)  {
         rhoRadPS += 
            pawPot.computeRhoPS (densityAtom.Dij(iSpin,is,ia),is,nSpin);
      }
      eBarRadial += tr (rhoRadPS(0,0) * pawPot.vBar(is));
   }

   PsiG rhoPSAtomG = gBasis | ((nSpin == 1) ? rhoAtom(0)
                                        : rhoAtom(0) + rhoAtom(1));
   PsiG vBarTotal(gBasis);
   vBarTotal.set(0.0);
   for (int is = 0; is < nSpecies; is++)
      vBarTotal += (gBasis | pawPot.vBar(is)) * gBasis.structureFactors(is);

   double eBarPS = dot(rhoPSAtomG, vBarTotal);
   
   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("eBarPS       = %18.12f H\n", eBarPS / norm);
      sxprintf ("eBarRadial   = %18.12f H\n", eBarRadial / norm);
      cout << SX_SEPARATOR;
   }

   if ((map.getSize() == structure.getNAtoms()))  {
      RhoR vElStat(1);
      vElStat(0) = vBarTotal.to (rBasis);
      SxRho(vElStat).writeRho ("vBarTotal-H.sxb");
   }

   return eBarPS - eBarRadial;
}

double computeXC (const SxPAWRho &densityAtom, 
                  const SxPAWRho &densityTot, 
                  SxPtr<SxXC> &xcPtr, 
                  const SxPAWPot &pawPot, 
                  const SxAtomicStructure &structure, 
                  const SxVector<Int> &map)
{
   SX_CLOCK(Timer::ComputeExc);
   int nSpin = densityTot.getNSpin ();
   const SxRBasis &rBasis = densityTot.pwRho.rhoR(0).getBasis<SxRBasis> ();
   RhoR rhoXcRAtom(nSpin);
   RhoR rhoXcRTot(nSpin);

   // define rhoXcR
   for (int iSpin=0; iSpin < nSpin; iSpin++)  {
      if (xcPtr->nlcc)  {
         cout << "nlcc not supportet" << endl;
         SX_QUIT;
      } else {
         rhoXcRAtom(iSpin) = 1.0 * densityAtom.pwRho.rhoR(iSpin);
         rhoXcRTot(iSpin) = 1.0 * densityTot.pwRho.rhoR(iSpin);
      }
   }

   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>(rhoXcRTot(0).getBasisPtr());
   double dOmega = rBasisPtr->dOmega;
   
   if (xcPtr->xcBasisPtr 
         && (xcPtr->xcBasisPtr.getPtr () != rBasisPtr))  {
      const SxGBasis &G = rBasisPtr->getGBasis ();
      xcPtr->xcBasisPtr->registerGBasis (G);
      dOmega = xcPtr->xcBasisPtr->dOmega;
      for (int iSpin=0; iSpin < nSpin; iSpin++)  {
         // Fourier interpolation
         rhoXcRAtom(iSpin) = rhoXcRAtom(iSpin).to(G).to(*xcPtr->xcBasisPtr);
         rhoXcRTot(iSpin) = rhoXcRTot(iSpin).to(G).to(*xcPtr->xcBasisPtr);
      }
   }

   SxXCFunctional xc (nSpin);
   size_t dim = rhoXcRTot(0).getSize ();
   SxMeshR epsXC (dim);
   SxArray<SxArray<SxMeshR> > gradRhoSpin(nSpin);

   if (xcPtr->xcFunctional == SxXC::PBE) {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  
         gradRhoSpin(iSpin) = xcPtr->gradient(rhoXcRTot(iSpin));
   }
   
   for (size_t i = 0; i < dim; i++)  {
      if (xcPtr->xcFunctional == SxXC::LDA)  {
         if (nSpin == 1) xc.computeLDA  (rhoXcRTot(0)(i));
         else  xc.computeLSDA (rhoXcRTot(0)(i), rhoXcRTot(1)(i));
      }
      else if (xcPtr->xcFunctional == SxXC::LDA_PW)
         if (nSpin == 1) {
            double rhoPart = 0.5 * rhoXcRTot(0)(i);
            xc.computeLSDA_PW  (rhoPart, rhoPart);
         }
         else xc.computeLSDA_PW (rhoXcRTot(0)(i),rhoXcRTot(1)(i));
      else if (xcPtr->xcFunctional == SxXC::PBE) {
         int up = 0, down = nSpin-1;
         SxVector3<TPrecRhoR> gradRhoUp, gradRhoDown;
         SxVector3<TPrecRhoR> gradRhoTl;
         for (int iDir = 0; iDir < 3; ++iDir)  {
            gradRhoUp(iDir)   = gradRhoSpin(up)(iDir)(i);
            gradRhoDown(iDir) = gradRhoSpin(down)(iDir)(i);
         }
         if (nSpin == 1)
            gradRhoTl = gradRhoUp;
         else
            gradRhoTl = gradRhoUp + gradRhoDown;

         double gradRhoNormUp   = gradRhoUp.norm (),
                gradRhoNormDown = gradRhoDown.norm (),
                gradRhoNormTl   = gradRhoTl.norm ();
         xc.computePBE(rhoXcRTot(up)(i),gradRhoNormUp,
                       rhoXcRTot(down)(i),gradRhoNormDown,
                       gradRhoNormTl);
      }
      else  {
         sxprintf ("Unknown XC functional.\n");
         SX_EXIT;
      }
      epsXC(i) = xc.epsXc;
   }

   SxMeshR rhoAtomTl = 1.0 * rhoXcRAtom(0);
   if (nSpin == 2) rhoAtomTl += rhoXcRAtom(1);
   double eXcPS = (rhoAtomTl * epsXC).sum () * dOmega;
   
   // PAW - 1 center contributions
   double eXcRadial = 0.0;
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++)  {
      int ia;
      int is = structure.getISpecies(map(iAtomX), &ia);
      SxRadialAtom radial(xcPtr->xcFunctional, pawPot.aGridType(is));
      SxArray<SxRadialMesh> rhoRadPS(nSpin), rhoRadAE(nSpin);
      for (int iSpin = 0; iSpin <  nSpin; ++iSpin)  {
         rhoRadPS(iSpin) 
            = pawPot.computeRhoPS (densityTot.Dij(iSpin,is,ia),is,nSpin);
         rhoRadAE(iSpin) 
            = pawPot.computeRhoAE (densityTot.Dij(iSpin,is,ia),is,nSpin);
      }
      if (pawPot.kjXcShape.getSize () > 0 
            && pawPot.kjXcShape(is).getSize () > 0)  {
         // TODO
         // --- Kresse-Joubert variant: include compensation charge
         // in pseudo-xc
         cout << SX_SEPARATOR;
         cout << " Kresse-Joubert variant" << endl;
         cout << SX_SEPARATOR;
         SX_EXIT;
      } else {
         radial.computeXC(rhoRadPS);
         eXcRadial -= radial.eXc;
      }
      radial.computeXC(rhoRadAE);
      eXcRadial += radial.eXc;
   }

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("eXcPS      = %18.12f H\n", eXcPS / norm);
      sxprintf ("eXcRadial  = %18.12f H\n", eXcRadial / norm);
      cout << SX_SEPARATOR;
   }

   if((map.getSize() == structure.getNAtoms())) {
      RhoR vXC(1);
      epsXC.setBasis(rBasis);
      vXC(0) = epsXC;
      SxRho(vXC).writeRho ("vXC-H.sxb");
   }
   
   return (eXcPS + eXcRadial);
}

PsiR getEkinSpillageDensity(const SxPW &fullWavesQ,
                            const SxPW &fullWaves,
                            const SxFermi &fermi, 
                            const SxPAWRho &rho)
{
   SX_CLOCK(Timer::ComputeEkinR);
   int nk       = fullWaves.getNk ();
   int nSpin    = fullWaves.getNSpin ();
   int nStates  = fullWaves.getNStates ();

   SxPW fullWavesS = wavesDiff(fullWaves,fullWavesQ); 
   
   const SxRBasis &R = rho.pwRho.rhoR(0).getBasis<SxRBasis> ();
   double dOmega = R.dOmega;
   PsiR eKinSpillageDensity(R);
   eKinSpillageDensity.set(0.0);

   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      const SxGBasis &g = fullWaves.getGkBasis ()(ik);
      double weight = fullWaves.getGkBasis ().weights(ik);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int iState = 0; iState < nStates; iState++)  {
            double occ = fermi.focc(iState,iSpin,ik);
            if (occ > 1e-12)  {
               PsiG psiGProj = (g | fullWavesS(iState,iSpin,ik));
               PsiG psiGFull = (g | fullWaves(iState,iSpin,ik));
               PsiR psiRProj = (R | psiGProj );
               PsiR D2PsiRFull = (R | g.g2 * psiGFull);
               eKinSpillageDensity += 0.5 * occ * weight * dOmega 
                                    * psiRProj.conj() * D2PsiRFull;
            }
         }
      }
   }

   SxLoopMPI::sum (eKinSpillageDensity);

   return eKinSpillageDensity;
}

double computeEKinR(const PsiR &eKinSpillageDensity,
                    const SxDiracVec<Double> &filter)
{
   SX_CLOCK(Timer::ComputeEkinR);
   double eKin = (filter * eKinSpillageDensity).sum().re;
   
   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("eKinSpillage    = %18.12f H\n", eKin / norm);
      cout << SX_SEPARATOR;
   }

   return eKin;
}

double computeEKinG(const SxPW &projWaves,
                    const SxPW &fullWaves,
                    const SxFermi &fermi, 
                    const SxPAWRho &rho,
                    const SxPAWPot &pawPot, 
                    const SxAtomicStructure &structure)
{
   SX_CLOCK(Timer::ComputeEkinG);
   double eKin = 0.0;
   
   int nk       = fullWaves.getNk ();
   int nSpin    = fullWaves.getNSpin ();
   int nStates  = fullWaves.getNStates ();
   int nSpecies = structure.getNSpecies ();
  
   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      const SxGBasis &g = fullWaves.getGkBasis ()(ik);
      double weight = fullWaves.getGkBasis ().weights(ik);
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int iState = 0; iState < nStates; iState++)  {
            double occ = fermi.focc(iState,iSpin,ik);
            if (occ > 1e-12)  {
               PsiG psiGProj = (g | projWaves(iState,iSpin,ik));
               PsiG psiGFull = (g | fullWaves(iState,iSpin,ik));
               PsiG D2PsiGFull = (g | g.g2 * psiGFull);
               eKin += 0.5 * occ * weight * dot(psiGProj,D2PsiGFull).re;
            }
         }
      }
   }

   SxLoopMPI::sum(eKin);

   double eKinPAWCorr = 0.0;
   // --- 1-center correction of kinetic energy
   for (int iSpin = 0; iSpin < nSpin; ++iSpin) {
      for (int is = 0; is < nSpecies; is++) {
         int nAtoms = structure.getNAtoms(is);
         for (int ia = 0; ia < nAtoms; ia++) {
            // Ref. 1, Eq. 20/21 with Eq. 52
            // \sum_ij Dij * ( <phiAE_j| -1/2 nabla^2 | phiAE_i>
            //                -<phiPS_j| -1/2 nabla^2 | phiPS_i>)
            int nProjType = pawPot.getNProjType (is);
            for (int ipt = 0; ipt < nProjType; ipt++)  {
               int l = pawPot.lPhi(is)(ipt);
               int offsetI = pawPot.offset(is)(ipt) + l;
               for (int jpt = 0; jpt < nProjType; jpt++)  {
                  if (l == pawPot.lPhi(is)(jpt))  {
                     int offsetJ = pawPot.offset(is)(jpt) + l;
                     for (int m = -l; m <= l; m++)  {
                        eKinPAWCorr += 
                           rho.Dij(iSpin,is,ia)(offsetI + m,offsetJ + m)
                           * pawPot.deltaKin(is)(jpt, ipt);

                     }
                  }
               }
            }
         }
      }
   }

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("eKin       = %18.12f H\n", eKin / norm);
      sxprintf ("eKinPAW    = %18.12f H\n", eKinPAWCorr / norm);
      cout << SX_SEPARATOR;
   }

   return eKin + eKinPAWCorr;
}

double computeEKin (const SxPW &projWavesQ,
                    const SxPW &fullWaves,
                    const SxFermi &fermi, 
                    const SxPAWRho &rho,
                    const SxPAWPot &pawPot, 
                    const SxAtomicStructure &structure,
                    const PsiR &eKinSpillageDensity,
                    const SxDiracVec<Double> &filter)
{
   double eKin = computeEKinG(projWavesQ, fullWaves, fermi, rho, pawPot, structure);
   eKin += computeEKinR(eKinSpillageDensity, filter);

   return eKin;
}

SxPAWRho combineRho(const SxPAWRho &projRhoQ,
                    const SxPAWRho &projRhoS)
{
   SxPAWRho result;

   result.potPtr = projRhoQ.potPtr;

   int nSpin = projRhoQ.getNSpin ();
   SX_CHECK(nSpin == projRhoS.getNSpin ());

   result.pwRho.rBasisPtr = projRhoQ.pwRho.rBasisPtr;
   result.pwRho.rhoR.resize(nSpin);
   //Spillage Dij is zero
   result.Dij = projRhoQ.Dij;
   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      // plane wave part
      result.pwRho(iSpin) = projRhoQ.pwRho(iSpin) 
                          + projRhoS.pwRho(iSpin);
   }
   
   return result;
}

double computeEtot (const SxPW &projWavesQ,
                    const SxPW &fullWaves,
                    const SxPAWRho &projRho, 
                    const SxPAWRho &fullRho,
                    const SxFermi &fermi,
                    SxPtr<SxPAWHamiltonian> HPtr, 
                    const SxAtomicStructure &structure, 
                    const SxVector<Int> &map,
                    const SxDiracVec<Double> &filter,
                    const PsiR &eKinSpillageDensity,
                    double charge)
{
   SX_CLOCK(Timer::ComputeEtot);
   double result = 0.0;
  
   //--- define potential
   SxPAWPot &pawPot = *(HPtr->potPtr);
   SxPtr<SxXC> xcPtr = HPtr->xcPtr;

   //--- eKin
   double eKin = computeEKin (projWavesQ, fullWaves, fermi, 
                              projRho, pawPot, structure,
                              eKinSpillageDensity, filter); 
   // --- eHartree
   double eHartree = computeEHartree (projRho, fullRho, pawPot, structure, map, charge);

   // --- eXc
   double eXc = computeXC (projRho, fullRho, xcPtr, pawPot, structure, map);

   // --- eBar
   double eBar = computeEBar (projRho, pawPot, structure, map);

   // core energy
   double eCore = 0.0;
   for (int iAtomX = 0; iAtomX < map.getSize (); iAtomX++) {
      int iSpecies = structure.getISpecies(map(iAtomX));
      eCore += pawPot.coreEnergy(iSpecies);
   }
      
   // --- output
   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("Local eKin    = %18.12f H\n", eKin / norm);
      sxprintf ("Local Hartree = %18.12f H\n", eHartree / norm);
      sxprintf ("Local eXc     = %18.12f H\n", eXc / norm);
      sxprintf ("Local eBar    = %18.12f H\n", eBar / norm);
      sxprintf ("Local eCore   = %18.12f H\n", eCore / norm);
      result = eKin + eHartree + eXc + eBar - eCore;
   
      sxprintf ("Local eTot    = %18.12f H\n",result / norm); 
      cout << SX_SEPARATOR;
   }

   return result;   
}

int main (int argc, char **argv)
{

   SX_START_TIMER (Timer::Startup);
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxLoopMPI::init (argc, argv);  // LoopMPI

   SxCLI cli (argc,argv);
   cli.authors = "B. Lange";

   SxString inputFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   SxString waveFile = cli.option ("-w|--waves", "file", "S/PHI/nX waves file")
                     .toString ("waves.sxb");
   
   SxString rhoFile = cli.option ("-r|--rho", "file", "S/PHI/nX rho file")
                     .toString ("rho.sxb");

   SxString basisFile = cli.option ("-b|--basis", "file", "S/PHI/nX basis file")
                     .toString ("basis.sx");
   double vAlign = cli.option ("--vAlign", "number", "constant Alignment potential")
                     .toDouble (0.0)/HA2EV;
   double rSoft = cli.option ("--rSoft", "number", "rSoft")
                     .toDouble (SxPAWHamiltonian::rSoft);
   keepWavesOnDisk = cli.option ("--keepWavesOnDisk", "keep Waves on disk").toBool ();
   saveMemory = cli.option ("--saveMemory", "save memory in G+k basis").toBool ();

   useSymmetries = cli.option ("--useSyms", "use Symmetries").toBool ();

   bool noSpillage = cli.option ("--noSpillage", "neglect Spillage contribution").toBool ();
   langeFreysoldt = cli.option ("--langeFreysoldt", "langeFreysoldt variant of EVHat").toBool ();

   SxArray<int> atomList = cli.option ("--atoms", "vector", 
                   "atoms of interest").toIntList();
   bool printRho = cli.option ("--printRho", "print electronic density in for each atom").toBool ();

   cli.finalize ();

   SxFFT::quickFFTPlanner ();

   SxPAWHamiltonian::rSoft = rSoft;

   SxParser parser;
   SxParser::Table table = parser.read (inputFile);
   if (table->containsGroup("pawPot"))
      cout << "PAW Hamiltonian found" << endl;  
   else {
      cout << "Tool works only for PAWHamiltonian" << endl;
      SX_QUIT;
   }
   SxString tmpDir = "";;
   if (keepWavesOnDisk) tmpDir = ".";
   SxHamSolver hamSolver (waveFile, rhoFile, table, tmpDir, saveMemory);

   const SxPW &waves = dynamic_cast<const SxPW&>(*hamSolver.wavesPtr);

   SxPtr<SxGkBasis> gkBasisPtr = waves.getGkBasisPtr ();
   SxFermi fermi = hamSolver.fermi;
 
   SxAtomicStructure &structure = hamSolver.structure;
   
   gkBasisPtr->changeTau(structure);

   SxPtr<SxPAWHamiltonian> HPtr 
         = SxPtr<SxPAWHamiltonian>(hamSolver.hamPtr);
   SxPtr<SxPartialWaveBasis> pBasis 
      = SxPtr<SxPartialWaveBasis>::create (HPtr->potPtr, structure);
   pBasis->createProjBasis (*gkBasisPtr);
   SxPtr<SxPAWOverlap> SPtr 
      = SxPtr<SxPAWOverlap>::create (pBasis, HPtr->potPtr);
   
   // get atomic space distrubution
   SxArray<SxVector<Int> > atomicSpace = atomicSpaceRegion (HPtr->rPtr, structure);
   
   SX_START_TIMER(Timer::AOBasis);
   SxPtr<SxAOBasis> aoPtr;
   try {
      SxParser aoParser;
      SxConstPtr<SxSymbolTable> aoTable = aoParser.read(basisFile);
      SxSymbolTable *aoGroup = aoTable->getGroup("AOBasis");
      SxAtomicOrbitals orbs;
      orbs.setup(aoGroup);
      aoPtr = SxPtr<SxAOBasis>::create (*gkBasisPtr, orbs, SPtr);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   SxAOBasis &aoBasis = *aoPtr;
   SX_MPI_MASTER_ONLY { cout << "Setting up ao basis ... done" << endl;}
   SX_STOP_TIMER(Timer::AOBasis);

   if (aoBasis.getNSpecies () != structure.getNSpecies ())  {
      cout << "Incompatible structure information in basis file and structure file" << endl;
      SX_QUIT;
   }

   // initial computation of PW-Energy
   HPtr->compute(waves, fermi);
   SxPAWRho rho = HPtr->pawRho;
   dipoleCorr = HPtr->dipoleCorr;
   double eTotPW = HPtr->eTotal;

   int nStates = waves.getNStates ();
   int nk      = waves.getNk();
   int nSpin   = waves.getNSpin ();
   aoBasis.setOverlapCaching (SxAOBasis::CacheCurrentK);
   aoBasis.setInvOverlapCaching (SxAOBasis::CacheCurrentK);
   SX_STOP_TIMER(Timer::Startup);

   SxVector<Int> mapAll = structure.match(SxGrid(structure,3),structure)->parentMap;
   SX_MPI_MASTER_ONLY { cout << "Structure matching ... done" << endl;}
   SxExpCoeffs C = computeExpCoeffs (aoBasis,waves);
   SxPW fullWavesQ = projectWaves(aoBasis, waves, C, structure, mapAll);
   SX_MPI_MASTER_ONLY { cout << "Quamol projection of waves ... done" << endl;}
   SxPAWRho fullRhoQ = computeRhoProjQ (fermi.focc, fullWavesQ, waves, rho,
                                      HPtr->rPtr, HPtr->potPtr, structure, mapAll);
   SX_MPI_MASTER_ONLY { cout << "Quamol projection of density ... done" << endl;}

   // --Print out Spillage for analysis
   SxVector<Double> stateSpillage(nStates);
   double spaceNorm = 0.0;
   stateSpillage.set(0.0);
   for (int ik = 0; ik < nk; ik++)  {
      if (!SxLoopMPI::myWork(ik)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         const SxDiracMat<Complex16> &PsiPWSet = waves(iSpin, ik);
         const SxDiracMat<Complex16> &PsiQuamolSet = fullWavesQ(iSpin, ik);
         SxDiracMat<Complex16> SPsiPWSet = SPtr->apply(waves(iSpin, ik));
         for (int iState = 0; iState < nStates; iState++)  {
            spaceNorm += fermi.kpPtr->weights(ik)
               * dot(PsiPWSet.colRef(iState),SPsiPWSet.colRef(iState)).re;
            stateSpillage(iState) += fermi.kpPtr->weights(ik)
               * (dot(PsiPWSet.colRef(iState),SPsiPWSet.colRef(iState)) 
                     - dot(PsiQuamolSet.colRef(iState),SPsiPWSet.colRef(iState))).re;
         }
      }
   }

   spaceNorm = SxLoopMPI::sum(spaceNorm);
   SxLoopMPI::sum (stateSpillage);
   stateSpillage /= spaceNorm;

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      cout << "Weightingfactors of each state equals 1!" << endl;
      sxprintf("Full Spillage is %18.12f\n",stateSpillage.sum());
      for (int iState = 0; iState < nStates; iState++)  {
         sxprintf("Spillage iState = %i is %18.12f\n",iState,stateSpillage(iState));
      }
      cout << SX_SEPARATOR;
   }

   // Get all local energies
   SX_START_TIMER(Timer::ComputeLocalEnergies);
   int nAtoms = structure.getNAtoms();
   SxVector<Double> eLocTot (nAtoms+1);
   SxVector<Double> atomCharge (nAtoms+1);

   // set atomic subset
   PsiR eKinSpillageDensity;
   if (noSpillage) 
      eKinSpillageDensity = getEkinSpillageDensity(fullWavesQ, fullWavesQ, fermi, rho);
   else
      eKinSpillageDensity = getEkinSpillageDensity(fullWavesQ, waves, fermi, rho);

   SxAtomicStructure choice;
   SxVector<Int> map;
   double eAlignTot = 0.0;

   SxArray<int> equivalentAtoms(nAtoms);
   equivalentAtoms.set(-1);
   structure.getInequivalentAtoms(&equivalentAtoms);

   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      sxprintf ("Atom %i equivalent to Atom %i\n", iAtom, equivalentAtoms(iAtom));
   }

   for (int iAtom = -1; iAtom < nAtoms; iAtom++)  {
      bool compute = true;
      norm = 1.0;

      SX_MPI_MASTER_ONLY {
         cout << SX_SEPARATOR;
         if (iAtom == -1) sxprintf ("All Atoms\n");
         else if (iAtom >= 0)  {
            Coord pos = structure(iAtom);
            sxprintf ("Atom %i @ (%.4f,%.4f,%.4f)\n", iAtom, pos(0), pos(1), pos(2));
         }
         cout << SX_SEPARATOR;
      }
      
      if (iAtom == -1) {
         choice = SxAtomicStructure(structure);
         map = structure.match(SxGrid(structure,3),choice)->parentMap;
      } else {
         SxList<SxList<Coord> > currentAtom;
         currentAtom.resize(structure.getNSpecies ());
         if (useSymmetries) {
            int equivalentAtom = equivalentAtoms(iAtom);
            int iSpecies = structure.getISpecies(equivalentAtom);
            for (int jAtom = 0; jAtom < nAtoms; jAtom++)  { 
               if (equivalentAtoms(jAtom) == equivalentAtom)
                  currentAtom(iSpecies).append(structure(jAtom));
            }
            if (iAtom == equivalentAtom) {
               cout << currentAtom(iSpecies).getSize () 
                    << " Atoms taken into account." << endl;
            }
            if (iAtom > equivalentAtom) {
               compute = false;
               eLocTot(iAtom + 1) = eLocTot(equivalentAtom + 1);
               atomCharge(iAtom+1) = atomCharge(equivalentAtom + 1);
               cout << "equivalent to Atom " << equivalentAtom << endl;
            }
         } else {
            int iSpecies = structure.getISpecies(iAtom); 
            currentAtom(iSpecies).resize(1);
            currentAtom(iSpecies)(0) = structure(iAtom);
         }
         if (atomList.getSize () > 0)  {
            compute = false;
            for (int iElem = 0; iElem < atomList.getSize (); iElem++)  {
               int jAtom = atomList(iElem);
               if (jAtom == iAtom) compute = true;
            }
         }

         choice = SxAtomicStructure(structure.cell, currentAtom);
         map = structure.match(SxGrid(structure,3),choice)->parentMap;
         if (useSymmetries) norm *= (double)map.getSize ();
      }

      if (compute)  {
      
         SxPW projWavesQ = projectWaves(aoBasis, waves, C, structure, map);
         SxPAWRho projRhoQ = computeRhoProjQ (fermi.focc, projWavesQ, waves, rho,
               HPtr->rPtr, HPtr->potPtr, structure, map);
         SxDiracVec<Double> filter = createMask (projRhoQ, fullRhoQ, atomicSpace, map);
         filter.setBasis(&*HPtr->rPtr);

         SX_MPI_MASTER_ONLY {
            cout << SX_SEPARATOR;
            cout << "Atom has a weight of " << filter.sum () / norm << endl;
            cout << SX_SEPARATOR;
         }
      
         //--- Set up densities
         SxPAWRho projRhoS;
         if (noSpillage)
            projRhoS = computeRhoProjS (fullRhoQ, fullRhoQ, filter);
         else 
            projRhoS = computeRhoProjS (fullRhoQ, rho, filter);

         // compute pseudo-core densities
         double projCore = 0.0;
         int nSpecies = structure.getNSpecies ();
         SxArray<double> pseudoCoreElectrons (nSpecies);
         const SxRadBasis &rad = HPtr->potPtr->getRadBasis ();
         for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
            SxDiracVec<Double> r = rad.radFunc(iSpecies);
            pseudoCoreElectrons(iSpecies) 
               = sqrt(FOUR_PI) 
               * (HPtr->potPtr->rhoCorePS(iSpecies) * r.cub ())
               .integrate(rad.logDr(iSpecies));
         }

         for (int iAtomTl = 0; iAtomTl < structure.getNAtoms (); iAtomTl++)  {
            int iSpecies = structure.getISpecies(iAtomTl);
            for (int iAtomX = 0; iAtomX < map.getSize(); iAtomX++)
               if (iAtomTl == map(iAtomX)) 
                  projCore += pseudoCoreElectrons (iSpecies);
         }
      
         SxPAWRho projRho = combineRho(projRhoQ,projRhoS);
         if (printRho) {
            SxString fileName = "projRho-" + SxString(iAtom) +".sxb";
            projRho.writeRho (fileName);
         }
         // compute valence charges
         double projRhoNormQ = projRhoQ.getNorm () - projCore;
         double projRhoNormS = projRhoS.getNorm ();
         double projRhoNorm  = projRho.getNorm () - projCore;

         // compute eAlign
         double chargeQ = projRhoNormQ;
         double chargeS = projRhoNormS;
         double charge  = projRhoNorm;
         for (int iAtomX = 0; iAtomX < map.getSize(); iAtomX++)  {
            int iSpecies = structure.getISpecies(map(iAtomX));
            chargeQ -= HPtr->potPtr->valenceCharge(iSpecies);
            charge  -= HPtr->potPtr->valenceCharge(iSpecies);
         }
         double eAlignQ = chargeQ * vAlign;
         double eAlignS = chargeS * vAlign;
         double eAlign  = charge * vAlign;
         if (iAtom == -1) eAlignTot = eAlign;
         SX_MPI_MASTER_ONLY {
            cout << SX_SEPARATOR;
            sxprintf ("chargeQ   = %18.12f\n", chargeQ / norm);
            sxprintf ("chargeS   = %18.12f\n", chargeS / norm);
            sxprintf ("charge    = %18.12f\n", charge / norm);
            sxprintf ("vAlign    = %18.12f\n", vAlign);
            sxprintf ("eAlignQ   = %18.12f\n", eAlignQ / norm);
            sxprintf ("eAlignS   = %18.12f\n", eAlignS / norm);
            sxprintf ("eAlign    = %18.12f\n", eAlign / norm);
            cout << SX_SEPARATOR;
         }

         atomCharge(iAtom + 1) = charge / norm;

         eLocTot(iAtom + 1) = computeEtot (projWavesQ, waves, 
               projRho, rho, fermi, HPtr, 
               structure, map, filter, 
               eKinSpillageDensity, charge);
      
         eLocTot(iAtom + 1) += eAlign;

         eLocTot(iAtom + 1) /= norm;

         SX_MPI_MASTER_ONLY {
            sxprintf ("Local eTot    = %18.12f H\n",eLocTot(iAtom + 1)); 
            cout << SX_SEPARATOR;
            sxprintf ("rhoSubSetQ valence    = %18.12f\n", projRhoNormQ / norm);
            sxprintf ("rhoSubSetS valence    = %18.12f\n", projRhoNormS / norm);
            sxprintf ("rhoSubSet valence     = %18.12f\n", projRhoNorm / norm);
            cout << SX_SEPARATOR;
         }
      }
   }
      
   SX_STOP_TIMER(Timer::ComputeLocalEnergies);

   SX_MPI_MASTER_ONLY {
      cout << SX_SEPARATOR;
      sxprintf ("ETot PW       = %18.12f H\n", eTotPW);
      sxprintf ("Full Spillage = %18.12f\n",stateSpillage.sum());
      sxprintf ("All Atoms     = %18.12f H\n", eLocTot(0));
      sxprintf ("eAlign        = %18.12f H\n", eAlignTot);
      sxprintf ("PW Diff       = %18.12f H\n", eLocTot(0) - eTotPW - eAlignTot);

      cout << SX_SEPARATOR;
      for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
         Coord pos = structure(iAtom);
         sxprintf ("Atom %i @ (%.4f,%.4f,%.4f) = %18.12f H, charge = %18.12f e-\n", 
               iAtom, pos(0), pos(1), pos(2), eLocTot(iAtom+1), atomCharge(iAtom+1));
      }
   
      SxIdx idx (1,nAtoms);
      cout << SX_SEPARATOR;
      sxprintf ("DIFF Etot = %18.12f H, charge = %18.12f e-\n", 
            eLocTot(idx).sum() - eLocTot(0), 
            atomCharge(idx).sum() - atomCharge(0));
      cout << SX_SEPARATOR;

      SxTimer::getGlobalTimer().print ();
      cout << "SxLocalEnergy completes successfully" << endl;
   }

   return 0;
}
