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

#include <SxPAWExchange.h>
#include <SxProjector.h>
#include <SxGrid.h>
#include <SxPAWHamiltonian.h>
#include <SxNeighbors.h>
#include <SxSpeciesRef.h>
#include <SxPAWSet.h>
#include <SxGaussIJ.h>
#include <SxTaskGroup.h>
#include <SxParallelHierarchy.h>
#include <SxXCFunctional.h>

#include <SxLoopMPI.h>

// Ref 1: J. Paier, R. Hirschl, M. Marsmann, G. Kresse,
//        J. Chem. Phys. 122, 234102 (2005)
// Ref 2: HSE implementation notes

#define FOCC_MIN 1e-12

namespace Timer {
   enum PAWExTimer { ComputePAWX, ApplyPAWX, ProjX, ProjPsi, PrepareRotq,
                     NHat, NHatPhase, Qrl, UKernelX, ComputeVkq, PrepareRot,
                     RotateU, RotatePx, Xij, OnSite, GetUrS,
                     VxToR, CompExchange, GRlExchange,
                     Qlm, Gkq2};
}

SX_REGISTER_TIMERS (Timer::PAWExTimer)
{
   using namespace Timer;
   regTimer (ComputePAWX, "Compute X");
   regTimer (ApplyPAWX, "Apply X");
   regTimer (ProjX, "<p|psi_x>");
   regTimer (ProjPsi, "<p|psi>");
   regTimer (PrepareRotq, "prepare rot q");
   regTimer (NHat, "comp. charge");
   regTimer (NHatPhase, "nHat phase");
   regTimer (Qrl, "U-type part");
   regTimer (GRlExchange, "gRlExchange");
   regTimer (UKernelX, "U kernel");
   regTimer (ComputeVkq, "Compute Vkq");
   regTimer (PrepareRot, "prepare rot");
   regTimer (RotateU, "rot. <r|psi_x>");
   regTimer (RotatePx, "rot. <p|psi_x>");
   regTimer (Xij, "Xij");
   regTimer (OnSite, "on-site");
   regTimer (GetUrS, "u(G)->u(r)");
   regTimer (VxToR, "Vij(G)->Vij(r)");
   regTimer (CompExchange, "comp. exchange");
   regTimer (Qlm, "QLM");
   regTimer (Gkq2, "|G+k-q|^2");
}

SxPAWExchange::SxPAWExchange ()
{
   // empty
}

void SxPAWExchange::setFocc (const Focc &foccIn)
{
   if (wavesPtr)  {
      SX_CHECK(foccIn.getNk () == wavesPtr->getNk (),
               foccIn.getNk (), wavesPtr->getNk ());
      SX_CHECK(foccIn.getNSpin () == wavesPtr->getNSpin (),
               foccIn.getNSpin (), wavesPtr->getNSpin ());
      SX_CHECK(foccIn.getNStates () == wavesPtr->getNStates (),
               foccIn.getNStates (), wavesPtr->getNStates ());
   }
   focc = foccIn;
}

double SxPAWExchange::computeV0 (const Coord &kVec) const
{
   SX_CHECK (wavesPtr);
   const SxGBasis &G = *gPtr;
   SX_CHECK (G.structPtr);
   SxSymGroup &syms = *G.structPtr->cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();
   int nSym2 = 2 * nSym; // include -k symmetry

   double qMax = 0.;
   const SxKPoints &kp = wavesPtr->getGkBasis ();
   for (int iq = 0; iq < kp.getNk (); ++iq)  {
      double q = kp.getK (iq).norm ();
      if (q > qMax) qMax = q;
   }
   double omega = sqrt(min(30., G.g2(G.ng-1)) - qMax) / 6.5;
   double a2 = 1. / sqr(omega);
   // cout << sqrt(a2) * 0.5 << endl;
   double gMax2 = sqr(6.5 * omega + qMax);
   
   // 4pi \int dg exp(-g^2/omega^2)/g^2
   double V0 = TWO_PI * SQRT_PI * omega;

   // cout << V0 << endl;
   double dq0 = 0.;
   // --- now substract numerical integral (excluding G+k-q=0)
   for (int iq = 0; iq < kp.getNk (); ++iq)  {
      double dq = kp.weights(iq) * TWO_PI * TWO_PI * TWO_PI 
               / G.structPtr->cell.volume / double(nSym2);
      for (int iSym = 0; iSym < nSym2; ++iSym)  {
         Coord qVec = syms.getSymmorphic (iSym % nSym) ^ kp.getK (iq);
         if (iSym >= nSym) qVec = -qVec;
         for (int ig = 0; ig < G.ng && G.g2(ig) < gMax2; ++ig)  {
            double Gkq2 = (G.getG(ig) + kVec - qVec).normSqr ();
            if (Gkq2 > 1e-12) {
               V0 -= dq * exp(-a2 * Gkq2) / Gkq2;
            } else  {
               V0 += dq * a2; // Gkq2->0 = 1 - a2 * Gkq
               dq0 += dq;
            }
         }
      }
   }
   // V0 now contains the dq-weighted correction term
   // multiply with 4pi since V(G) = 4pi/G^2
   return (dq0 == 0.) ? 0. : (FOUR_PI * V0 / dq0);

}

void SxPAWExchange::prepareRotData(int iq, Coord &kVec) const
{
   SX_CLOCK (Timer::PrepareRotq);
   SX_CHECK (wavesPtr);
   SX_CHECK (gPtr);
   SX_CHECK (rPtr);
   SX_CHECK (potPtr);
   const SxGkBasis &gkBasis = wavesPtr->getGkBasis ();
   const SxRBasis &R = *rPtr;
   const SxGBasis &G = *gPtr;
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   const SxMesh3D mesh = R.fft3d.mesh;
   const SxCell   &cell = R.cell;
   SxCell meshCell (cell.basis(0) / mesh(0),
                    cell.basis(1) / mesh(1),
                    cell.basis(2) / mesh(2));
   int nSym = cell.symGroupPtr->getNSymmorphic ();
   int nSym2 = 2 * nSym;

   // --- value for v(G=0)
   double vG0 = 0.;
   bool screened = (SxXCFunctional::omegaHSE > 1e-10);
   bool timeReversal = false;
   if (screened) {
      vG0 = PI / sqr(SxXCFunctional::omegaHSE);
   } else {
      vG0 = computeV0 (kVec);
   }

   vkq.resize (nSym2);
   symRot.resize (nSym);
   YlmGl.resize (nSym2);
   gHard.resize (nSym2);
   gSoft.resize (nSym2);
   if (screened)  {
      smallG.resize (nSym2);
      smallG2.resize (nSym2);
      smallV.resize (nSym2);
   }
   if (pwProjComp) {
      pwProjComp->cacheRefOrb = SxAOBasis::CacheAll;
      pwProjComp->SPtr = SxPtr<SxPWOverlap>::create ();
      pwProjComp->refOrbitals.resize (nSym2);
      ssize_t nQlm = ((potPtr->lMaxRho + 1).sqr () * structure.atomInfo->nAtoms).sum ();
      pwProjComp->orbitalMap.resize (nQlm);
      for (int is = 0, io = 0; is < structure.getNSpecies (); ++is)
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia)
            for (int l = 0, lm = 0; l <= potPtr->lMaxRho(is); ++l)
               for (int m = -l; m <= l; ++m, ++lm)
                  pwProjComp->orbitalMap (io++) = SxAOBasis::OrbitalIndex (is, ia, lm);

   }
   idxGkq0.resize (nSym2);
   idxGkq0.set (-1);
   for (int iSym = 0; iSym < nSym2; ++iSym)  {
#if defined USE_LOOPMPI && ! defined NDEBUG
      if (SxLoopMPI::getLevel ()->getName ().tokenize (':').last () != "iSym")  {
         cout << "Wrong MPI level " << SxLoopMPI::getLevel ()->getName ()
              << endl;
         cout << "Should be ...:iSym" << endl;
         SX_EXIT;
      }
#endif
      if (!SxLoopMPI::myWork(iSym)) continue;
      // S rotates ks->(+-)q
      SxMatrix3<Double> S = cell.symGroupPtr->getSymmorphic (iSym % nSym);
      if (iSym >= nSym) timeReversal = true;
      Coord qVec = S ^ gkBasis.getK(iq);
                      
      if (timeReversal) qVec *= -1.;
      Coord dkq = kVec - qVec;

      //bool useFolding =    (gkBasis.foldingMP(0) > 0)
      //                  && (gkBasis.foldingMP(1) > 0)
      //                  && (gkBasis.foldingMP(2) > 0);
      bool useFolding = false;
      {
         SxStack<ssize_t> whichSmallG;
         SxStack<double> absG;
         SxStack<double> vLR;

         // --- setup v = 4pi/|G+q-k|^2
         //SX_CLOCK (Timer::ComputeVkq);
         SxDiracVec<Double> &vkqRot = vkq(iSym);
         vkqRot.resize (G.ng);

         SxCell recCell = cell.getReciprocalCell ();
         for (int ig = 0; ig < G.ng; ++ig)  {
            Coord Gkq = G.getG (ig) + dkq;
            double gkq2 = Gkq.normSqr ();
            if (gkq2 < 1e-16)  {
               SX_CHECK (idxGkq0(iSym) == -1); // can't happen twice
               vkqRot(ig) = vG0;
               idxGkq0(iSym) = ig;
               continue;
            }
            Coord GkqRel = recCell.carToRel (Gkq) * gkBasis.foldingMP;
            if (useFolding
                && fabs(GkqRel(0)) < 0.4999
                && fabs(GkqRel(1)) < 0.4999
                && fabs(GkqRel(2)) < 0.4999)
            {
               int N = 5;
               double N3 = double(2 * N + 1); N3 = (N3 * N3 * N3);
               double vAvg = vG0;
               SxVector3<Int> i;
               Coord foldInv = 1./Coord(gkBasis.foldingMP);
               double minG2 = 1e20;
               foldInv /= double(2 * N + 1);
               for (i(0) = -N; i(0) <= N; ++i(0))  {
                  for (i(1) = -N; i(1) <= N; ++i(1))  {
                     for (i(2) = -N; i(2) <= N; ++i(2))  {
                        Coord dq = recCell.relToCar (i * foldInv);
                        double g2 = (Gkq + dq).normSqr ();
                        if (g2 < minG2)  {
                           vAvg += FOUR_PI / minG2;
                           minG2 = g2;
                        } else {
                           vAvg += FOUR_PI / g2;
                        }
                     }
                  }
               }
               vkqRot(ig) = vAvg / N3;
            } else {
               vkqRot(ig) = FOUR_PI / gkq2;
            }
            if (screened && gkq2 < 169. * sqr(SxXCFunctional::omegaHSE))  {
               double omega2 = sqr(SxXCFunctional::omegaHSE);
               double gaussHSE = exp(-0.25/omega2 * gkq2 );
               whichSmallG << ig;
               absG        << gkq2;
               // Ref 2, Eq. (3b)
               vLR         << gaussHSE * vkqRot(ig);
               // Ref 2, Eq. (2)
               vkqRot(ig) *= 1. - gaussHSE;
            }
         }
         if (screened)  {
            if (whichSmallG.getSize () == 0)  {
               cout << "small G iSym =" << iSym << ": size is zero" << endl;
               SX_EXIT;
            }
            smallG(iSym) = whichSmallG;
            smallG2(iSym) = absG;
            smallV(iSym) = vLR;
         }
      }

      //int ng = gkBasis(0).ng;
      //vkqRot(SxIdx(ng, vkqRot.getSize () -1)).set (0.); // Fourier filter

      // --- set up rotation index map
      if (iSym < nSym) {
         //SX_CLOCK (Timer::PrepareRot);
         // get associated real-space rotation (transpose=inverse)
         // cf. Ref 2, Eq. (36)
         SxMatrix3<Int> rot = meshCell.carToRel (S.transpose ());

         SxArray<ssize_t> &idxmap = symRot(iSym);
         idxmap.resize (R.getNElements ());
         ssize_t ir, irRot;
         SxVector3<Int> rVec;

         for (rVec(0)=0; rVec(0) < mesh(0); rVec(0)++)  {
            for (rVec(1)=0; rVec(1) < mesh(1); rVec(1)++)  {
               for (rVec(2)=0; rVec(2) < mesh(2); rVec(2)++)  {
                  ir = mesh.getMeshIdx(rVec, SxMesh3D::Positive);
                  // compute index after rotation
                  irRot = mesh.getMeshIdx(rot^rVec, SxMesh3D::Unknown);
                  idxmap(ir) = irRot;
               }
            }
         }
      }

      SxDiracVec<Double> gkq2 = getGkq2 (dkq);
      
      // --- |G+k-q|^l * Ylm(G+k-q)
      {
         SxDiracVec<Double> gkqAbs = sqrt(gkq2), gkqL(G);
         gkqL.copy (gkqAbs);
         int lMaxRho = potPtr->lMaxRho.maxval ();
         YlmGl(iSym) = SxRadBasis::realYlm(lMaxRho, G, dkq);
         for (int lm = 1, l = 1; l <= lMaxRho; ++l)  {
            // Ylm * g^l
            for (int m = -l; m <= l; ++m,++lm)  {
               //YlmGl(iSym).colRef(lm) *= pow (gkq2, 0.5 * l);
               YlmGl(iSym).colRef(lm) *= gkqL;
            }
            if (l < lMaxRho) gkqL *= gkqAbs;
         }
      }


      // soft compensation charges
      gSoft(iSym) = exp(-(0.25 * sqr(SxPAWHamiltonian::rSoft)) * gkq2);

      // --- hard compensation charge shapes
      gHard(iSym) = SxDiracMat<Double> (G.ng, structure.getNSpecies ());
      {
         if (pwProjComp) pwProjComp->refOrbitals(iSym).resize (structure.getNSpecies ());
         for (int is = 0; is < structure.getNSpecies (); ++is)  {
            int nProjRef = sqr(potPtr->lMaxRho(is) + 1);
            if (pwProjComp) {
               pwProjComp->refOrbitals(iSym)(is).reformat (G.ng, nProjRef);
               pwProjComp->refOrbitals(iSym)(is).setBasis (G);
            }
            double rc2 = 0.25 * sqr(potPtr->rc(is));
            SxDiracVec<Double> gauss = exp( -rc2 * gkq2);
            gHard(iSym).colRef(is) <<= gauss;
            if (pwProjComp)  {
               double prefactor = FOUR_PI / sqrt(structure.cell.volume);
               for (int l = 0, lm = 0; l <= potPtr->lMaxRho(is); ++l)  {
                  
                  prefactor /= 2 * l + 1;

                  for (int m = -l; m <= l; ++m, ++lm)  {
                     PsiG projRef = gauss * YlmGl(iSym).colRef (lm);
                     projRef *= prefactor * SxYlm::getYlmNormFactor(l, m);
                     pwProjComp->refOrbitals(iSym)(is).colRef(lm) <<= projRef;
                  }
               }
            }
         }
      }
   }
}

void SxPAWExchange::rotCleanup () const
{
   vkq.resize (0);
   YlmGl.resize (0);
   symRot.resize (0);
   pwProjComp = SxPtr<SxAOBasis> ();
}

enum CTimer { kLoop, Summation, Get_nG, vG };
SX_REGISTER_TIMERS(CTimer)
{
   regTimer (kLoop, "k loop");
   regTimer (Summation, "sum");
   regTimer (Get_nG, "nG");
   regTimer (vG, "v in G (hard)");
}

enum XContrib {
   XijTerms  = 0x00100000, // small
   PSTerms   = 0x00200000,
   CompTerms = 0x00400000, // small
   HardSoft  = 0x00800000, // small
   Comp0     = 0x01000000,
   OnSiteLR  = 0x02000000,
   AllXTerms = 0x0ff00000

};
const int xContrib = AllXTerms & (~OnSiteLR);

double SxPAWExchange::compute (const SxConstPtr<SxPWSet> &wavesIn,
                               const Focc &foccIn,
                               const SxRadMat &Dij)
{
   cout.precision (15);
   SX_CLOCK (Timer::ComputePAWX);
   SX_CHECK (wavesIn);
   if (wavesPtr.getPtr () != NULL && (&wavesPtr->getGkBasis () != &wavesIn->getGkBasis ())) 
      cout << "WARNING: SxPAWExchange::compute changes pointer to wave functions!" << endl;
   wavesPtr = wavesIn;
   focc = foccIn;
   SX_CHECK (gPtr);
   SX_CHECK (rPtr);
   SX_CHECK (potPtr);
   SX_CHECK (potPtr->QijL.getSize () > 0);
   SX_CHECK (pBasis);
   {
      // --- check if pBasisX needs to be updated
      const SxPAWSet* pawSet=dynamic_cast<const SxPAWSet*>(wavesPtr.getPtr ());
      if (pawSet)  {
         SxConstPtr<SxPartialWaveBasis> pBasisW
            = dynamic_cast<const SxPAWBasis&>(pawSet->getBasis (0)).pBasis;
         if (pBasisW.getPtr () != pBasisX.getPtr ())  {
            cout << "Warning: using exchange waves with new basis" << endl;
            pBasisX = pBasisW;
         }
      }
   }
   const SxPWSet &waves = *wavesPtr;
   const SxGkBasis &gkBasis = waves.getGkBasis ();
   const SxRBasis &R = *rPtr;
   const SxGBasis &G = *gPtr;
   const int jBlockSize = 16;
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   SX_CHECK (structure.cell.symGroupPtr);
   int nSym = structure.cell.symGroupPtr->getNSymmorphic ();
   int nSym2 = 2 * nSym;

   bool screened = (SxXCFunctional::omegaHSE > 1e-10);
   if (screened && (xContrib & OnSiteLR))  {
      cout << SX_SEPARATOR;
      cout << "| WARNING: on-site LR corrections are switched on" << endl
           << "|          for the screened exchange energy."
           << endl;
      cout << "|          The contribution to the Hamiltonian is missing."
           << endl;
      cout << SX_SEPARATOR;
   }

   int nk = waves.getNk ();
   int nSpin = waves.getNSpin ();

   SxComplex16 energy = 0.;
   double energy58b = 0.;

   {
      int nStates = waves.getNStates ();
      SxAutoLevel ()
         .append ("pawxc:comp:ik", nk)
         .append ("pawxc:comp:iq", nk)
         .append ("pawxc:comp:iSpin", nSpin)
         .append ("pawxc:comp:i", nStates)
         .append ("pawxc:comp:iSym", nSym2);
   }
      

   if (!pwProjComp)
      pwProjComp = pwProjComp.create ();

   // --- prepare Xij and U kernel
   computeXij (Dij);

   // --- compute (pseudo+compensation) exchange energy
   for (int ik = 0; ik < nk; ++ik)  { // nk:1
      SX_MPI_LEVEL("pawxc:comp:ik");
      if (! SxLoopMPI::myWork(ik)) continue;

      SX_CLOCK (kLoop);
      Coord kVec = gkBasis.getK (ik);

      SxDiracVec<TPrecCoeffG> pX, pXrot, pPsi;
      SxDiracVec<TPrecCoeffG> uR, uRs;

      for (int iq = 0; iq < nk; ++iq)  { // nk:1
         SX_MPI_LEVEL("pawxc:comp:iq");
         if (! SxLoopMPI::myWork(iq)) continue;

         double dq = gkBasis.weights(iq) * gkBasis.weights(ik) / double(nSym2);
         if (dq < FOCC_MIN) continue; // skip k-points with no weight

         // --- compute for all S ^ q
         // * vkq         v(k-q)
         // * symRot     real-space rotation map
         // * YlmGl       |G+k-q|^l Ylm(G+k-q)
         // * pwProjComp  hard compensation projectors
         { 
            SX_MPI_LEVEL("pawxc:comp:iSpin");
            SX_MPI_LEVEL("pawxc:comp:i");
            SX_MPI_LEVEL("pawxc:comp:iSym");
            prepareRotData (iq, kVec);
         }

         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {  // 1
            SX_MPI_LEVEL("pawxc:comp:iSpin");
            if (! SxLoopMPI::myWork(iSpin)) continue;

            {
               SX_CLOCK (Timer::ProjPsi);
               pPsi = *pBasisX | waves(iSpin,ik);
            }
            {
               SX_CLOCK (Timer::ProjX);
               pX = *pBasisX | waves(iSpin,iq);
            }

            for (ssize_t j0 = 0; j0 < waves.getNStates (ik); j0 += jBlockSize)  {
               int jBlock = min(jBlockSize, int(waves.getNStates (ik) - j0));
               for ( ; jBlock > 0 && focc(j0 + jBlock - 1, iSpin, ik) < FOCC_MIN ; jBlock--) {}
               if (jBlock <= 0)  {
                  // no more occupied states in current block -> run up with j0, probably until j0=nStates
                  for (++j0 ; j0 < waves.getNStates (ik)
                              && focc(j0, iSpin, ik) < FOCC_MIN; ++j0) {}
                  continue;
               }

               SxArray<PsiR> psiR(jBlock);
               SX_LOOP(j)
                  psiR(j) = R | waves(j0 + j, iSpin, ik);
               //double fJ = focc(j, iSpin, ik);
               //if (fJ < FOCC_MIN) continue;

               SxDiracVec<TPrecCoeffG> pPsiJ;
               pPsiJ = pPsi(SxIdx(int(pPsi.nRows () * j0), 
                                  int(pPsi.nRows () * (j0 + jBlock) - 1)));
               pPsiJ.reshape (pPsi.nRows (), jBlock);

               // --- loop over exchange states
               for (int i = 0; i < waves.getNStates (iq); ++i)  { // Si-8 : 8 states
                  SX_MPI_LEVEL("pawxc:comp:i");
                  if (! SxLoopMPI::myWork(i)) continue;

                  // --- get occupation number
                  double fOccI = focc(i, iSpin, iq);
                  if (nSpin == 1) fOccI *= 0.5; // only one spin channel!
                  // ignore unoccupied states
                  if (fabs(fOccI * focc(j0,iSpin,ik)) < FOCC_MIN) continue;

                  // get u in real space at special q
                  { 
                     SX_CLOCK (Timer::GetUrS);
                     uRs = (R | waves(i, iSpin, iq));
                  }

                  for (int iSym = 0; iSym < nSym2; ++iSym)  { // Si8:96
                     SX_MPI_LEVEL("pawxc:comp:iSym");
                     if (! SxLoopMPI::myWork(iSym)) continue;

                     bool conj = (iSym >= nSym); // time reversal k->(-k)
                     SymMat S = structure.cell.symGroupPtr->getSymmorphic (iSym % nSym);
                     Coord qVec = S ^ gkBasis.getK (iq);
                     if (iSym >= nSym) qVec *= -1.;
                     Coord dkq = kVec - qVec;

                     {
                        // --- get rotated u 
                        SX_CLOCK (Timer::RotateU);
                        uR = uRs.getSorted (symRot(iSym % nSym));
                        if (conj) uR = uR.conj (); // Ref 2, Eq. (37)
                     }

                     // rotated pX
                     pXrot = getRotated(pX.colRef(i), gkBasis.getK (iq), iSym);
                     // get compensation moments
                     SxQlm qlm = getQlm (pXrot, pPsiJ);
                     // get compensation charge phases
                     SxArray<PsiG> nHatPhase;
                     if (xContrib & CompTerms)  {
                        nHatPhase = getNHatUnshaped (qlm,iSym,dkq);
                     } else {
                        nHatPhase.resize (structure.getNSpecies ());
                        PsiG zero;
                        zero.reformat (G.ng, jBlock);
                        zero.setBasis (G);
                        zero.set (0.);
                        zero.handle->auxData.ik = iSym; // cf. getNHatUnshaped
                        for (int is = 0; is < nHatPhase.getSize (); ++is)  {
                           nHatPhase(is) = zero;
                        }
                     }

                     PsiG vGSoft;
                     vGSoft.reformat (G.ng, jBlock);
                     vGSoft.setBasis (G);
                     vGSoft.handle->auxData.ik = iSym; // needed for pwProjComp
                     for (int j1 = 0; j1 < jBlock; ++j1)  {
                        // n*(r') = u*(r') psi(r') in G-space
                        // Ref. 1, Eq. 8 ; also Ref 2, Eq. (31)
                        PsiG nG;
                        if (xContrib & PSTerms)  {
                           SX_CLOCK (Get_nG);
                           nG = (G | (uR.conj () * psiR(j1)));
                        } else {
                           nG = PsiG (G);
                           nG.set (0.);
                        }

                        PsiG nHatHard = getNHat (nHatPhase, j1);
                        PsiG nHatSoft = getNHatSoft (nHatPhase,j1);

                        // two-orbital exchange potential - Ref. 1, Eq. (25)
                        SX_START_TIMER (vG);
                        PsiG vGHard = vkq(iSym) * (nG + nHatHard);
                        SX_STOP_TIMER (vG);
                        vGSoft.colRef(j1) <<= vkq(iSym) * (nG + nHatSoft);

                        if (idxGkq0(iSym) >= 0 && (xContrib & Comp0)
                            && !screened)
                        {
                           // --- special treatment for G+k-q = 0
                           // include 2nd order term
                           // vkq = 4pi/|G+k-q|^2
                           // nHat = nHatPhase * (1 - rc^2/4 * |G+k-q|^2)
                           // => extra -PI nHatPhase rc^2
                           // cf. Hartree energy in PAW Hamiltonian
                           double rSoft2 = sqr(SxPAWHamiltonian::rSoft);
                           for (int is = 0; is < nHatPhase.getSize (); ++is)  {
                              double rHard2 = sqr(potPtr->rc(is));
                              SxComplex16 nHat0  = nHatPhase(is)(idxGkq0(iSym), j1);
                              vGSoft(idxGkq0(iSym),j1) -= PI * rSoft2 * nHat0;
                              vGHard(idxGkq0(iSym)) -= PI * rHard2 * nHat0;
                           }
                        }

                        double fOcc = fOccI * focc(j0+j1, iSpin, ik);
                        {
                           SX_CLOCK(Summation);
                           // Ref. 2, Eq. (40a)
                           if (xContrib & PSTerms)  {
                              SxComplex16 dE = dq * fOcc * dot(nG, vGHard);
                              energy -= dE;
                              //exij -= dE.re;
                           }

                           // Ref. 2, Eq. (40b)
                           if (xContrib & CompTerms)  {
                              SxComplex16 dE = dq * fOcc
                                             * dot (nHatHard, vGSoft.colRef (j1));
                              energy -= dE;
                           }

                           if (idxGkq0(iSym) >= 0 && (xContrib & Comp0)
                               && !screened)
                           {
                              // --- special treatment for G+k-q = 0
                              // include 2nd order term
                              // vkq = 4pi/|G+k-q|^2
                              // nHat = nHatPhase * (1 - rc^2/4 * |G+k-q|^2 + ...)
                              // => extra -PI nHatPhase rc^2
                              // cf. Hartree energy in PAW Hamiltonian
                              // note: nHatSoft(idxGkq0) = nHatHard(idxGkq0)
                              SxComplex16 n0 = nG(idxGkq0(iSym)) 
                                             + nHatHard(idxGkq0(iSym));
                              SxComplex16 dE = 0.;
                              for (int is = 0; is < nHatPhase.getSize (); ++is) {
                                 double rHard2 = sqr(potPtr->rc(is));
                                 SxComplex16 nHat0 
                                    = nHatPhase(is)(idxGkq0(iSym), j1);
                                 if (!(xContrib & CompTerms))  {
                                   nHat0 = 0.;
                                   for (int ia=0; ia<structure.getNAtoms(is); ++ia)
                                     nHat0 += qlm(is)(ia)(0, j1);
                                   nHat0 *= sqrt(FOUR_PI / structure.cell.volume);
                                 }
                                 dE -= PI * rHard2 * n0 * nHat0.conj ();
                              }
                              dE *= dq * fOcc;
                              energy -= dE;
                              //exij -= dE.re;
                           }
                        }
                        
                        if (xContrib & HardSoft && screened)  {
                           SxComplex16 dE = 0.;
                           // Ref. 2, Eq. (70)
                           SX_LOOP(igSmall)  {
                              ssize_t ig = smallG(iSym)(igSmall);
                              SxComplex16 nH = nHatHard(ig),
                                          nS = nHatSoft(ig);
                              dE -= (nH.conj () * (nH - nS))
                                 * smallV(iSym)(igSmall);
                           }
                           if (idxGkq0(iSym) >= 0)  {
                              // G+k-q = 0 case
                              SxComplex16 nH0 = 0.;
                              double r2Soft = sqr(SxPAWHamiltonian::rSoft);
                              SX_LOOP(is)
                                 nH0 += nHatPhase(is)(idxGkq0(iSym), j1)
                                      * (sqr(potPtr->rc(is)) - r2Soft);
                              dE += nHatHard(idxGkq0(iSym)).conj () * PI * nH0;
                           }
                           dE *= dq * fOcc;
                           energy -= dE;
                        }

                        if (screened && (xContrib & OnSiteLR))  {
                           SxComplex16 M0;
                           PsiG deltaN = computePAWCorrG (pXrot, pPsiJ.colRef(j1),
                                                          iSym, dkq, &M0);
                           PsiG vDeltaN = deltaN * smallV(iSym);
                           double dE = 0.;
                           SX_LOOP(ig)  {
                              ssize_t idxg = smallG(iSym)(ig);
                              SxComplex16 nPSComp = nG(idxg) + nHatHard(idxg);
                              // Ref. 2, Eq. 58b
                              dE += ((2. * nPSComp + deltaN(ig)).conj ()
                                     * vDeltaN(ig)                      ).re;
                           }
                           if (idxGkq0(iSym) >= 0)  {
                              ssize_t ig0 = idxGkq0(iSym);
                              // special treatment for G=0
                              // vDeltaN = (Y_00 * M0 * G^2) * 4pi/G^2
                              dE += 2. * sqrt(FOUR_PI)
                                  * ((nG(ig0) + nHatHard(ig0)).conj () * M0).re; 
                           }

                           energy += dq * fOcc * dE;
                           energy58b += dq * fOcc * dE;
                        }
                     }
                     if (xContrib & HardSoft)  {
                        // Ref. 2, Eq. (40c)
                        SxDiracVec<TPrecFocc> foccJ
                           = focc(iSpin,ik)(SxIdx((int)j0, int(j0 + jBlock -1)));
                        energy -= dq * fOccI * qRlExchange (dkq, qlm, NULL, foccJ);
                     }
                  }
                  // cout << endl;
                  /*
                  cout << "(" << i << "," << iq
                       << "|X(" << iSpin << ")|"
                       << j << "," << ik << ")="
                       << exij << endl;
                  */
               } // i loop over iq states
            } // j loop over ik states
         } // iSpin loop
      } // iq loop
   } // ik loop

   // clear cached data
   rotCleanup ();

   energy.re = SxLoopMPI::sum(energy.re);
   energy.im = SxLoopMPI::sum(energy.im);

   cout << "ps+comp energy=" << (0.5 * energy) << endl;
   if (screened && (xContrib & OnSiteLR))
      cout << "energy58b = " << 0.5 * energy58b << endl;
   // --- on-site terms
   if (xContrib & XijTerms) {
      SX_CLOCK (Timer::OnSite);
      energy += tr(Dij, Xij);
   }

   energy *= 0.5; // double counting
   cout << "energy=" << energy << endl;
   return energy;
}

PsiG SxPAWExchange::apply (const PsiG &psi) const
{
   SX_CLOCK (Timer::ApplyPAWX);
   SX_CHECK (psi.handle);
   SX_CHECK (wavesPtr);
   SX_CHECK (gPtr);
   SX_CHECK (rPtr);
   SX_CHECK (potPtr);
   SX_CHECK (potPtr->QijL.getSize () > 0);
   SX_CHECK (pBasis);
   // --- get the G basis
   const SxGBasis *gBasisRes 
      = dynamic_cast<const SxGBasis*>(psi.getBasisPtr ());
   if (!gBasisRes)  {
      const SxPAWBasis *pawBasis
        = dynamic_cast<const SxPAWBasis*>(psi.getBasisPtr ());
      SX_CHECK (pawBasis);
      gBasisRes = pawBasis->gBasis.getPtr ();
   }
   // get partial-wave and radial contributions in G-space
   const SxGBasis &resG = *gBasisRes;

   const int jBlockSize = 16;
   int nStates = (int)psi.nCols ();

   if (nStates > jBlockSize)  {
      PsiG res;
      res.reformat (gBasisRes->ng, nStates);
      res.setBasis (gBasisRes);
      for (int j0 = 0; j0 < nStates; j0 += jBlockSize)  {
         int block = min(jBlockSize, nStates - j0);
         SxIdx iRes(j0 * resG.ng, (j0 + block) * resG.ng - 1),
               iPsi(j0 * (int)psi.nRows (), (j0 + block) * (int)psi.nRows () - 1);
         PsiG psiBlock = psi(iPsi);
         psiBlock.reshape (psi.nRows (), block);
         res (iRes) <<= apply(psiBlock);
      }
      return res;
   }
   const SxPWSet &waves = *wavesPtr;
   const SxGkBasis &gkBasis = waves.getGkBasis ();
   const SxRBasis &R = *rPtr;
   const SxGBasis &G = *gPtr;
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   SX_CHECK (structure.cell.symGroupPtr);
   int nSym = structure.cell.symGroupPtr->getNSymmorphic ();
   int nSym2 = 2 * nSym;
   bool screened = (SxXCFunctional::omegaHSE > 1e-10);

   Coord kVec = resG.getK ();

   int nSpin = waves.getNSpin ();
   int nk = waves.getNk ();
   int iSpin = psi.handle->auxData.iSpin;

   SxAutoLevel ()
      .append ("pawxc:apply:iq", nk)
      .append ("pawxc:apply:i", waves.getNStates ())
      .append ("pawxc:apply:iSym", nSym2);

   // initialize hard compensation projectors
   if (!pwProjComp)
      pwProjComp = pwProjComp.create ();

   SxDiracVec<TPrecCoeffG> pPsi = *pBasis | psi;
   SxArray<PsiR> psiR(nStates);
   SX_LOOP(j) psiR(j) = R | (resG | psi.colRef(j));

   // --- compute on-site exchange terms
   SxDiracMat<Complex16> Xin;
   SxDiracMat<Complex16> Xin0;
   if (xContrib & XijTerms)  {
      Xin0 = Xij ^ pPsi;
   } else {
      Xin0 = 0. * pPsi;
   }
   Xin = 0. * pPsi;

   // --- compute (pseudo+compensation) exchange energy
   SxDiracVec<TPrecCoeffG> pX, pXrot;
   SxDiracVec<TPrecCoeffG> uR, uRs;
   SxArray<SxDiracVec<TPrecCoeffG> > xPsi(nStates);
   SX_LOOP (j)  {
      xPsi(j) = SxDiracVec<TPrecCoeffG>(R);
      xPsi(j).set (0.);
   }

   for (int iq = 0; iq < nk; ++iq)  { // Si8:nk:1
      SX_MPI_LEVEL("pawxc:apply:iq");
      if (! SxLoopMPI::myWork(iq)) continue;

      double dq = gkBasis.weights(iq) / double(nSym2);
      if (dq < FOCC_MIN) continue; // skip k-points with no weight

      // --- compute for all S ^ q
      // * vkq         v(k-q)
      // * symRot      real-space rotation map
      // * YlmGl       |G+k-q|^l Ylm(G+k-q)
      // * pwProjComp  hard compensation projectors
      { 
         SX_MPI_LEVEL("pawxc:apply:i");
         SX_MPI_LEVEL("pawxc:apply:iSym");
         prepareRotData (iq, kVec);
      }

      SX_START_TIMER (Timer::ProjX);
      pX = *pBasisX |waves(iSpin,iq);
      SX_STOP_TIMER (Timer::ProjX);

      // --- loop over exchange states
      for (int i = 0; i < waves.getNStates (iq); ++i)  {
         SX_MPI_LEVEL("pawxc:apply:i");
         if (! SxLoopMPI::myWork(i)) continue;


         // --- get occupation number
         double fOcc = focc(i, iSpin, iq);
         if (nSpin == 1) fOcc *= 0.5; // only one spin channel!
         if (fOcc < FOCC_MIN) continue; // ignore unoccupied states

         // get u in real space at special q
         { 
            SX_CLOCK (Timer::GetUrS);
            uRs = (R | waves(i, iSpin, iq));
         }

         for (int iSym = 0; iSym < nSym2; ++iSym)  { // Si8 : nSym2==96
            SX_MPI_LEVEL("pawxc:apply:iSym");
            if (! SxLoopMPI::myWork(iSym)) continue;

            //(cout << '.').flush ();

            bool conj = (iSym >= nSym); // time reversal k->(-k)
            SymMat S = structure.cell.symGroupPtr->getSymmorphic (iSym % nSym);
            Coord qVec = S ^ gkBasis.getK (iq);
            if (iSym >= nSym) qVec *= -1.;
            Coord dkq = kVec - qVec;

            {
               // --- get rotated u 
               SX_CLOCK (Timer::RotateU);
               uR = uRs.getSorted (symRot(iSym % nSym));
               if (conj) uR = uR.conj (); // Ref 2, Eq. (37)
            }

            // rotated pX
            pXrot = getRotated(pX.colRef(i), gkBasis.getK (iq), iSym);

            // get compensation moments
            SxQlm qlm = getQlm (pXrot, pPsi);
            // get compensation charge phases
            SxArray<PsiG> nHatPhase;
            if (xContrib & CompTerms)  {
               nHatPhase = getNHatUnshaped (qlm,iSym,dkq);
            } else {
               nHatPhase.resize (structure.getNSpecies ());
               PsiG zero;
               zero.reformat (G.ng, nStates);
               zero.setBasis (G);
               zero.set (0.);
               zero.handle->auxData.ik = iSym; // cf. getNHatUnshaped
               for (int is = 0; is < nHatPhase.getSize (); ++is)  {
                  nHatPhase(is) = zero;
               }
            }

            PsiG vGSoft;
            vGSoft.reformat (G.ng, nStates);
            vGSoft.setBasis (G);
            vGSoft.handle->auxData.ik = iSym; // needed for pwProjComp
            SxArray<SxComplex16> n0(nStates);
            n0.set (0.);
            for (int j = 0; j < nStates; ++j)  {
               // n*(r') = u*(r') psi(r') in G-space
               // Ref. 1, Eq. 8 ; also Ref 2, Eq. (31)
               PsiG nG;
               if (xContrib & PSTerms)  {
                  SX_CLOCK (Get_nG);
                  nG = (G | (uR.conj () * psiR(j)));
               } else {
                  nG = PsiG (G);
                  nG.set (0.);
               }

               // two-orbital exchange potential - Ref. 1, Eq. (25)
               PsiG nHatHard = nG + getNHat    (nHatPhase, j);
               PsiG nHatSoft = nG + getNHatSoft(nHatPhase, j);
               PsiG vGHard = vkq(iSym) * nHatHard;
               vGSoft.colRef(j) <<= vkq(iSym) * nHatSoft;

               if (idxGkq0(iSym) >= 0 && (xContrib & Comp0) && !screened)  {
                  // --- special treatment for G+k-q = 0
                  // include 2nd order term
                  // vkq = 4pi/|G+k-q|^2
                  // nHat = nHatPhase * (1 - rc^2/4 * |G+k-q|^2)
                  // => extra -PI nHatPhase rc^2
                  // cf. Hartree energy in PAW Hamiltonian
                  double rSoft2 = sqr(SxPAWHamiltonian::rSoft);
                  for (int is = 0; is < nHatPhase.getSize (); ++is)  {
                     double rHard2 = sqr(potPtr->rc(is));
                     SxComplex16 nHat0  = nHatPhase(is)(idxGkq0(iSym), j);
                     vGSoft(idxGkq0(iSym),j) -= PI * rSoft2 * nHat0;
                     vGHard(idxGkq0(iSym)) -= PI * rHard2 * nHat0;
                  }
               }
               if (screened)  {
                  // gradient of Ref. 2, Eq. (70)
                  SX_LOOP(igSmall)  {
                     ssize_t ig = smallG(iSym)(igSmall);
                     SxComplex16 nH = nHatHard(ig),
                                 nS = nHatSoft(ig);
                     vGSoft(ig,j) -= (nH - nS) * smallV(iSym)(igSmall);
                  }
                  if (idxGkq0(iSym) >= 0)  {
                     // G+k-q = 0 case
                     SxComplex16 nH0 = 0.;
                     double r2Soft = sqr(SxPAWHamiltonian::rSoft);
                     SX_LOOP(is)
                        nH0 += nHatPhase(is)(idxGkq0(iSym), j)
                             * (sqr(potPtr->rc(is)) - r2Soft);
                     vGSoft(idxGkq0(iSym),j) += PI * nH0;
                  }
               }

               // Ref. 1, left term of Eq. 26 = Eq. (27), Ref. 2, Eq. (42a)
               if (xContrib & PSTerms)  {
                  SX_CLOCK (Timer::VxToR);
                  xPsi(j).plus_assign_ax (-dq * fOcc, uR * (R | vGHard));
               }

               if (idxGkq0(iSym) >= 0 && (xContrib & Comp0) && !screened)  {
                  n0(j) = nG(idxGkq0(iSym));
                  for (int is = 0; is < nHatPhase.getSize (); ++is)
                    n0(j) += nHatPhase(is)(idxGkq0(iSym), j);
               }
               // cout << endl;
            }
            // Ref. 1, right term of Eq. 26 = Eq. (28) (without p_i)
            // Ref 2, Eq. (42b) without |p_i>
            Xin.plus_assign_ax (-dq * fOcc,
                                compensationExchange (vGSoft,
                                               dkq, qlm, pXrot, n0));
         }
      }
   }
   // clear cached data
   rotCleanup ();
   (cout << Xin.getSize () << endl).flush ();
   
   SX_MPI_SOURCE(CurrentLevel,TaskGroupAll);
   SX_MPI_TARGET(CurrentLevel,TaskGroupAll);
   SxLoopMPI::sum(Xin);
   for (int j = 0; j < nStates; ++j) SxLoopMPI::sum(xPsi(j));

   Xin += Xin0;
   
   VALIDATE_VECTOR (Xin);
#ifndef NDEBUG
   SX_LOOP(j) VALIDATE_VECTOR (xPsi(j));
#endif
   //cout << "Xin: " << Xin.normSqr () << endl;
   //cout << "xPsi: " << xPsi.normSqr () << endl;
   PsiG res = (resG | Xin);
   SX_LOOP(j)  {
      res.colRef (j) += (resG | xPsi(j));
   }
   return  res;
}

SxPAWExchange::SxQlm SxPAWExchange::gRlExchange (const PsiG &vG,
                                                 const Coord &dkq,
                                                 const SxQlm &qlm) const
{
   SX_CLOCK (Timer::GRlExchange);
   SX_CHECK (potPtr);
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   int nSpecies = str.getNSpecies ();
   ssize_t nStates = vG.nCols ();
   
   SxQlm res(nSpecies);
   // Ref 2, Eq. (43a)
   SxDiracVec<TPrecCoeffG> gRlV;
   if (xContrib & CompTerms)  {
      gRlV = pwProjComp->fromPWBasis (vG);
   } else {
      gRlV.reformat (pwProjComp->orbitalMap.getSize (), nStates);
      gRlV.set (0.);
   }

   int ipg = 0;
   for (int is = 0; is < nSpecies; ++is)  {
      int nlm = sqr(potPtr->lMaxRho(is) + 1);
      int nAtoms = str.getNAtoms (is);
      res(is).resize (nAtoms);

      // --- extract <gRL | vSoft> from gRlV and add missing (k-q) phase factors
      //     and i^l factor 
      for (int ia = 0; ia < nAtoms; ++ia)  {
         res(is)(ia).reformat (nlm, nStates);
         // e^{+i (k-q) r} (so far only e^{-i G r} in gRL(G) )
         SxComplex16 kPhase = exp (I * (dkq ^ str.getAtom (is, ia)));

         for (int iState = 0; iState < nStates; ++iState)  {
            for (int l = 0, lm = 0; l <= potPtr->lMaxRho(is); ++l)  {
               SxComplex16 il = (l & 1) ? -I : SxComplex16(1.);
               if (l & 2) il=-il;
               for (int m =-l; m <= l; ++m, ++lm)  {
                  res(is)(ia)(lm,iState) = gRlV(ipg + lm, iState) * kPhase * il;
               }
            }
         }
         ipg += nlm;
         VALIDATE_VECTOR (res(is)(ia));
      }
   }

   if (xContrib & HardSoft)  {
      // --- hard/soft correction (cf SxPAWHamiltonian::computeU)
      // so far we have computed   (gRl | v[nPs + nHat'])
      // correction is             (gRl | v[nHat - nHat'])
      qRlExchange (dkq, qlm, &res);
   }
   return res;
}

void SxPAWExchange::computeUKernel ()
{
   SX_CLOCK (Timer::UKernelX);
   SX_CHECK (potPtr);
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   const SxPAWPot &pot = *potPtr;
   int nSpecies = str.getNSpecies ();

   // --- for neighbor search
   SxGrid grid (str, 10);
   SxNeighbors neighbors;
   int neighborParams = SxNeighbors::StoreRel
                      | SxNeighbors::IncludeZeroDistance;
   // note: for rcut=13, 1 bohr wide Gaussians behave like point multipoles
   // within 1e-16
   double rcut = 13. * SxPAWHamiltonian::rSoft;

   const SxYlm::SxClebschTable &cg = pot.clebschGordan;
   SxMatrix<Double> delta;
   
   ssize_t deltaMem = 0;
   neighborDist.resize (str.getSize ());
   kernelU.resize (str.getSize ());

   SxGaussIJ gaussIJ;
   for (int is = 0, iTlAtom=0; is < nSpecies; ++is)  {
      // correction is q^*(RL) g(RL) | v | q(R'L') [g(R'L') - g'(R'L')]
      int lmaxI = pot.lMaxRho(is);
      
      // --- loop over atoms (index i => gRl)
      for (int ia = 0; ia < str.getNAtoms(is); ++ia, iTlAtom++)  {
         // compute neighbors within rcut
         neighbors.compute (grid, str, str.getAtom(is,ia), 
                            rcut, neighborParams);
         neighborDist(iTlAtom).copy (neighbors.relPositions);
         kernelU(iTlAtom).resize (neighbors.relPositions.getSize ());

         // --- loop over neighbors j (=> nHat)
         for (int js=0, jTlAtom = 0; js < str.getNSpecies (); js++)  {
            int lmaxJ = pot.lMaxRho(js);
            double a2Hard = sqr(pot.rc(is)) + sqr(pot.rc(js));
            double a2Soft = sqr(pot.rc(is)) + sqr(SxPAWHamiltonian::rSoft);

            for (int ja = 0; ja < neighbors.relPositions.getNAtoms(js); ++ja)
            {
               // i => R
               // j => R'
               // r = R'-R
               Coord r = neighbors.relPositions.getAtom(js, ja);
               gaussIJ.setDelta (r, a2Hard, a2Soft, lmaxI + lmaxJ);
               delta = gaussIJ.compute (lmaxI, lmaxJ, cg).getCopy ();
               deltaMem += sizeof(double) * delta.getSize ();
               kernelU(iTlAtom)(jTlAtom++) = delta;
            }
         }
      }
   }
   cout << "kernelU mem=" << deltaMem << endl;
}

SxComplex16
SxPAWExchange::qRlExchange (const Coord &dkq,
                            const SxQlm &qlm,
                            SxQlm *grl,
                            const SxDiracVec<TPrecFocc> &foccJ) const
{
   SX_CLOCK (Timer::Qrl);
   SX_CHECK (potPtr);
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   int nSpecies = str.getNSpecies ();

   SxComplex16 res = 0.; // energy (all atoms)
#ifdef USE_OPENMP
#  pragma omp parallel
#endif
   {
      SxComplex16 energy = 0.; // energy (1 thread)
#ifdef USE_OPENMP
#     pragma omp for
#endif
      for (int iTlAtom = 0; iTlAtom < str.getNAtoms (); ++iTlAtom)  {
         int ia, is = str.getISpecies (iTlAtom, &ia);
         const SxVector<Int> &neighborIdx 
            = neighborDist(iTlAtom).atomInfo->parentMap;
         // --- loop over neighbors j (=> nHat)
         for (int js=0; js < nSpecies; js++)  { // --- species j
            for (int ja = 0; ja < neighborDist(iTlAtom).getNAtoms(js); ++ja)
            {
               int jTlAtom = neighborDist(iTlAtom).getIAtom (js, ja);
               // i => R
               // j => R'
               // r = R'-R
               Coord r = neighborDist(iTlAtom).getAtom(jTlAtom);
               const SxMatrix<Double> &delta = kernelU(iTlAtom)(jTlAtom);

               // sum over moments of j=R'
               int jRefAtom = neighborIdx(jTlAtom)
                            - str.atomInfo->offset(js);

               Coord R = r + str(is,ia) - str(js, jRefAtom);
               SxComplex16 kPhase = exp (I * (dkq ^ R));
               // --- the following part is time-critical, so we
               //     do not use a matrix notation: matrices/vectors
               //     are tiny, and intermediate memory allocation
               //     becomes the bottleneck (note: delta is real,
               //     qlm is complex, so delta must be promoted first)
               if (grl)  {
                        SxVector<Complex16> &gRL  = (*grl)(is)(ia);
                  const SxVector<Complex16> &QLMj = qlm(js)(jRefAtom);
                  ssize_t nStates = QLMj.nCols ();
                  for (ssize_t iState = 0; iState < nStates; ++iState)  {
                     for (int ilm = 0; ilm < delta.nRows (); ++ilm)  {
                        SxComplex16 sum = 0.;
                        for (int jlm = 0; jlm < delta.nCols (); ++jlm)
                           sum += delta(ilm, jlm) * QLMj(jlm, iState);
                        gRL(ilm, iState) += kPhase * sum;
                     }
                  }
               } else {
                  const SxVector<Complex16> &QLMi = qlm(is)(ia),
                                            &QLMj = qlm(js)(jRefAtom);
                  SxComplex16 sumI = 0.;
                  ssize_t nStates = QLMj.nCols ();
                  SX_CHECK (foccJ.getSize () == nStates);
                  for (ssize_t iState = 0; iState < nStates; ++iState)  {
                     for (int ilm = 0; ilm < delta.nRows (); ++ilm)  {
                        SxComplex16 sumJ = 0.;
                        for (int jlm = 0; jlm < delta.nCols (); ++jlm)
                           sumJ += delta(ilm, jlm) * QLMj(jlm, iState);
                        sumI += foccJ(iState) * QLMi(ilm, iState).conj () * sumJ;
                     }
                  }
                  energy += kPhase * sumI;
               }
            }
         }
#ifndef NDEBUG
         if (grl) { VALIDATE_VECTOR ((*grl)(is)(ia)); }
#endif
      }
      // --- get energy from all threads
#ifdef USE_OPENMP
#     pragma omp atomic
      res.re += energy.re;
#     pragma omp atomic
      res.im += energy.im;
#else
      res = energy;
#endif
   }
   return res;
}


SxDiracVec<TPrecCoeffG> 
SxPAWExchange::compensationExchange (const PsiG &vG,
                                     const Coord &dkq,
                                     const SxQlm &qlm,
                                     const SxDiracVec<TPrecCoeffG> &pX,
                                     const SxArray<SxComplex16> &n0) const
{
   SX_CLOCK (Timer::CompExchange);
   SX_CHECK (pX.nCols () <= 1, pX.nCols ());
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;
   ssize_t nStates = vG.nCols ();
   SX_CHECK (n0.getSize () == nStates, n0.getSize (), nStates);


   // Compute exchange with normalized compensation charges
   // Ref 2, Eq. (43)
   SxQlm xHatRl = gRlExchange (vG, dkq, qlm);

   bool screened = SxXCFunctional::omegaHSE > 1e-10;
   if ((xContrib & Comp0) && !screened)  {
      for (ssize_t iState = 0; iState < nStates; ++iState)  {
         if (n0(iState).absSqr () > 1e-28 )  {
            // --- special treatment for G+k-q = 0
            // include 2nd order term
            // vkq = 4pi/|G+k-q|^2
            // nHat = nHatPhase * (1 - rc^2/4 * |G+k-q|^2 + ...)
            // => extra -PI nHatPhase rc^2
            // cf. Hartree energy in PAW Hamiltonian
            SX_LOOP(is)  {
               double rHard2 = sqr(potPtr->rc(is));
               SxComplex16 x0 = -PI * rHard2 * n0(iState) 
                              * sqrt(FOUR_PI / str.cell.volume);
               SX_LOOP(ia) xHatRl(is)(ia)(0,iState) += x0;
            }
         }
      }
   }

   SxDiracVec<TPrecCoeffG> res;
   res.reformat (pX.getSize (), nStates);
   res.set (0.);
   // TODO: ref
   int offset = 0;
   for (int is = 0; is < str.getNSpecies (); ++is)  {
      int npt = potPtr->getNProjType (is);
      int npl = pot.getNProj (is);

      // --- set up i^l 
      SxVector<Complex16> il(npt); // i^l
      for (int ipt = 0; ipt < npt; ++ipt)  {
         int l = pot.lPhi(is)(ipt);
         SxComplex16 c = (l & 1) ? I : SxComplex16(1.);
         if (l & 2) c=-c;
         il(ipt) = c;
      }

      SX_LOOP (ia)  {
         SX_LOOP2 (ipt,jpt)  {
            int l1 = pot.lPhi(is)(ipt), off1 = pot.offset(is)(ipt) + l1;
            int l2 = pot.lPhi(is)(jpt), off2 = pot.offset(is)(jpt) + l2;

            const SxVector<Double> &Qij = pot.QijL(is)(ipt,jpt);

            int maxL = min(l2 + l1, pot.lMaxRho(is));

            for (ssize_t iState = 0; iState < nStates; ++iState) {
               // now the angular momentum sums
               for (int m1 = -l1; m1 <= l1; ++m1)  {
                  int lm1 = SxYlm::combineLm(l1,m1);
                  for (int m2 = -l2; m2 <= l2; ++m2)  {
                     int lm2 = SxYlm::combineLm(l2,m2);

                     SxComplex16 sum = 0.;
                     // attach i^l phase
                     SxComplex16 pB = pX(offset + off1 + m1) * il(ipt).conj ();

                     // Ref ... with precomputed radial integrals
                     // Qij(l) = \int r^2 dr phi_i(r) phi_j(r) r^l
                     // loop over radial mesh's lm's
                     // possible l from triangular condition & parity
                     for (int l = abs(l2-l1); l <= maxL; l+=2)  {
                        for (int m = -l; m <= l; ++m)  {
                           int lm = SxYlm::combineLm(l,m);
                           double cg = pot.clebschGordan(lm1, lm2, lm);
                           // no i^l factors for PAW radial ylm
                           // but is l1+l2-l (from cg definition)
                           if ((l1 + l2 - l) & 2) cg = -cg;
                           // Ref. 2, Eq. ???
                        double QijL = cg * Qij(l);
                        sum += pB * QijL * xHatRl(is)(ia)(lm, iState);
                        }
                     }
                     res(offset + off2 + m2, iState) += il(jpt) * sum;
                  }
               }
            }
         }
         offset += npl;
      }
   }
   VALIDATE_VECTOR (res);
   return res;
}




SxPAWExchange::SxQlm
SxPAWExchange::getQlm (const SxDiracVec<TPrecCoeffG> pA,
                       const SxDiracVec<TPrecCoeffG> pB) const
{
   SX_CLOCK (Timer::Qlm);
   SX_CHECK (potPtr);
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   const SxPAWPot &pot = *potPtr;
   int nSpecies = str.getNSpecies ();
   SxQlm qlm(nSpecies);
   int offset = 0;
   int nB = (int)pB.nCols ();
   for (int is = 0; is < nSpecies; ++is)  {
      int nAtoms = str.getNAtoms (is);
      int npl = pot.getNProj (is);
      int npt = pot.getNProjType (is);

      // --- set up i^l 
      SxVector<Complex16> il(npt); // i^l
      for (int ipt = 0; ipt < npt; ++ipt)  {
         int l = pot.lPhi(is)(ipt);
         SxComplex16 c = (l & 1) ? I : SxComplex16(1.);
         if (l & 2) c=-c;
         il(ipt) = c;
      }

      qlm(is).resize (nAtoms);
      for (int ia = 0; ia < nAtoms; ++ia)  {
         SxVector<Complex16> &qRlm = qlm(is)(ia);
         qRlm.reformat(sqr(pot.lMaxRho(is) + 1), nB);
         qRlm.set (0.);
         // --- loop over partial wave types
         for (int ipt = 0; ipt < npt; ++ipt)  {
            int l1 = pot.lPhi(is)(ipt), off1 = pot.offset(is)(ipt) + l1;
            for (int jpt = 0; jpt < npt; ++jpt)  {
               int l2 = pot.lPhi(is)(jpt), off2 = pot.offset(is)(jpt) + l2;

               const SxVector<Double> &Qij = pot.QijL(is)(ipt,jpt);

               int maxL = min(l2 + l1, pot.lMaxRho(is));

               // now the angular momentum sums
               for (int m1 = -l1; m1 <= l1; ++m1)  {
                  int lm1 = SxYlm::combineLm(l1,m1);
                  for (int m2 = -l2; m2 <= l2; ++m2)  {
                     int lm2 = SxYlm::combineLm(l2,m2);
                     
                     // Ref. 2, Eq. (34): projection A + i^l phases
                     SxComplex16 ilpA = pA(offset + off1 + m1).conj ()
                                      * il(ipt) * il(jpt).conj ();

                     // Ref. 2, Eq. (35) with precomputed radial integrals
                     // Qij(l) = \int r^2 dr phi_i(r) phi_j(r) r^l
                     // loop over radial mesh's lm's
                     // possible l from triangular condition & parity
                     for (int l = abs(l2-l1); l <= maxL; l+=2)  {
                        for (int m = -l; m <= l; ++m)  {
                           int lm = SxYlm::combineLm(l,m);
                           double cg = pot.clebschGordan(lm1, lm2, lm);
                           // no i^l factors for PAW radial ylm
                           // but is l1+l2-l (from cg definition)
                           if ((l1 + l2 - l) & 2) cg = -cg;
                           // Ref. 2, Eq. (35)
                           double QijL = cg * Qij(l);
                           for (int ib = 0; ib < nB; ++ib)  {
                              // Ref. 2, Eq. (34)
                              qRlm(lm, ib) += QijL * ilpA
                                            * pB(offset + off2 + m2, ib);
                           }
                        }
                     }
                  }
               }
            }
         }
         VALIDATE_VECTOR (qRlm);
         offset += npl;
      }
   }
   return qlm;
}

SxArray<PsiG> 
SxPAWExchange::getNHatUnshaped (const SxQlm &qlm,
                                int iSym,
                                const Coord &dkq) const
{
   SX_CLOCK (Timer::NHatPhase);
   SX_CHECK (wavesPtr);
   const SxGBasis &G = *gPtr;
   SX_CHECK (G.structPtr);
   const SxAtomicStructure &structure = *G.structPtr;
   int nSpecies = structure.getNSpecies ();
   SxArray<PsiG> res(nSpecies);

   /*
   const SxDiracMat<Double> &Ylmgl = YlmGl(iSym);
   for (int is = 0; is < nSpecies; ++is)  {
      int nlm = sqr(potPtr->lMaxRho(is) + 1);
      SxDiracVec<Complex16> weights(nlm);
      res(is) = PsiG (G);
      res(is).set (0.);
      res(is).handle->auxData.ik = iSym;
      for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {
         // e^{iGr}
         PsiG T = G.getPhaseFactors (is,ia);
         //PsiG sum(G.ng);
         //sum.set (0.);
         // e^{-i(k-q)r}
         SxComplex16 kPhase = exp(-I * (dkq ^ structure.getAtom(is, ia)));

         double prefac = FOUR_PI / sqrt(structure.cell.volume);
         for (int l = 0; l <= potPtr->lMaxRho(is); ++l)  {
            prefac /= 2 * l + 1;
            // i^l factor
            SxComplex16 il = 1.;
                 if ((l & 3) == 1) il = I; 
            else if ((l & 3) == 2) il = -1.;
            else if ((l & 3) == 3) il = -I;
            for (int m = -l; m <= l; ++m)  {
               int lm = SxYlm::combineLm(l,m);
               double N = SxYlm::getYlmNormFactor(l,m) * prefac;
               weights(lm) = N * qlm(is)(ia)(lm) * il * kPhase;
               // Ref 2, Eq. (33) (without gRl radial shape, though)
               //sum.plus_assign_ax (weights(lm), Ylmgl.colRef(lm));
            }
         }
         PsiG sum = Ylmgl ^ weights;
         {
            // --- variant 1
            // res(is) += sum * T;

            // --- variant 2 (should be 25% faster than variant 1)
            PsiG &resIs = res(is);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (ssize_t ig = 0; ig < G.ng; ++ig)
               resIs(ig) += sum(ig) * T(ig);
         }
      }
      VALIDATE_VECTOR (res(is));
   }
   */
   
   for (int is = 0; is < nSpecies; ++is)  {
      int nlm = sqr(potPtr->lMaxRho(is) + 1);
      int nAtoms = structure.getNAtoms (is);

      // --- collect qlms for this species with all prefactors
      int nPairs = (nAtoms > 0) ? (int)qlm(is)(0).nCols () : 0;
      SxDiracVec<Complex16> qlms;
      qlms.reformat (nlm * nAtoms, nPairs);
      for (int ia = 0, ipg = 0; ia < nAtoms; ++ia)  {
         SxComplex16 kPhase = exp(-I * (dkq ^ structure.getAtom(is, ia)));
         double prefac = FOUR_PI / sqrt(structure.cell.volume);
         for (int l = 0; l <= potPtr->lMaxRho(is); ++l)  {
            prefac /= 2 * l + 1;
            // i^l factor
            SxComplex16 il = 1.;
                 if ((l & 3) == 1) il = I; 
            else if ((l & 3) == 2) il = -1.;
            else if ((l & 3) == 3) il = -I;
            for (int m = -l; m <= l; ++m)  {
               int lm = SxYlm::combineLm(l,m);
               double N = SxYlm::getYlmNormFactor(l,m) * prefac;
               for (int iPair = 0; iPair < nPairs; ++iPair)
                  qlms(ipg, iPair) = N * qlm(is)(ia)(lm, iPair) * il * kPhase;
               ipg++;
            }
         }
      }
      VALIDATE_VECTOR (qlms);

      // --- set up fast projector
      SxDiracVec<Complex16> myYlm = YlmGl(iSym)(SxIdx(0, nlm * G.ng-1));
      myYlm.reshape (G.ng, nlm);
      SxAOBasis phaseProj(myYlm, is);
      res(is) = phaseProj.toPWBasis (&G, qlms);
      res(is).setBasis (G);
      res(is).handle->auxData.ik = iSym;
   }
   return res;
}

SxDiracVec<Double> SxPAWExchange::getGkq2(const Coord &dkq) const
{
   SX_CLOCK (Timer::Gkq2);
   SxDiracVec<Double> gkq2(*gPtr);
   gkq2.set (0.);
   for (int i = 0; i < 3; ++i)
      gkq2 += (gPtr->gVec.colRef(i) + dkq(i)).sqr ();
   return gkq2;
}

PsiG SxPAWExchange::getNHat (const SxArray<PsiG> &nHatUnshaped, ssize_t j) const
{
   SX_CLOCK (Timer::NHat);
   int iSym = nHatUnshaped(0).handle->auxData.ik;
   SX_CHECK (gPtr);
   PsiG res(*gPtr);
   res.set (0.);

   // use "hard" Gaussian shape from paw potential
   for (int is = 0; is < nHatUnshaped.getSize (); ++is)  {
      res += gHard(iSym).colRef(is) * nHatUnshaped(is).colRef (j);
   }
   return res;
}

PsiG
SxPAWExchange::getNHatSoft (const SxArray<PsiG> &nHatUnshaped, ssize_t j) const
{
   SX_CLOCK (Timer::NHat);
   int iSym = nHatUnshaped(0).handle->auxData.ik;
   SX_CHECK (gPtr);
   PsiG res(*gPtr);
   res.set (0.);
   // sum phases (all same shape)
   for (int is = 0; is < nHatUnshaped.getSize (); ++is)  {
      res += nHatUnshaped(is).colRef (j);
   }
   // use "soft" Gaussian shape from paw potential
   res *= gSoft(iSym);
   return res;
}

SxDiracVec<TPrecCoeffG> 
SxPAWExchange::getRotated (const SxDiracVec<TPrecCoeffG> &proj,
                           const Coord &kVec, int iSym) const
{
   SX_CLOCK (Timer::RotatePx);
   SX_CHECK (potPtr);
   const SxPAWPot &pot = *potPtr;
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   SX_CHECK (str.cell.symGroupPtr);
   const SxSymGroup &syms = *str.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();
   SxGrid grid (str, 10);
   bool timeReversal = (iSym >= nSym);
   if (timeReversal) iSym -= nSym;

   SxDiracVec<TPrecCoeffG> res (proj.getSize ());
   // --- symmetry operation
   // <p(lm,tau)|psi_S(k)>
   // = <p(lm,tau)|S(G+k)><S(G+k)|psi_S(k)>
   // = p(|G+k|) Y_lm[S(G+k)]^* exp(i[S(G+k)][tau]) <G+k|psi_k>
   // = p(|G+k|) Dmm'(S) Y_lm'(G+k)^* exp(i(G+k)[S^{-1}tau] <G+k|psi_k>
   // = D_mm'(S) <p(lm',S^{-1}tau)|psi_k> 

   // --- time reversal
   // <p(lm,tau)|psi_(-k)>
   // = <p(lm,tau)|r><r|psi_k>*
   // = [<p(lm,tau)|psi_k>]* if <p|r> is real
   // in our case, <p|r> is real for even l
   //              <p|r> is purely imaginary for odd l

   int offsetIs = 0;
   for (int is = 0; is < str.getNSpecies (); ++is)  {
      int npl = pot.getNProj (is);
      for (int ia = 0; ia < str.getNAtoms (is); ++ia)  {
         // --- find equivalent atom
         Coord rotTau = syms.getSymmorphic(iSym).transpose ()
                      ^ str.getAtom(is,ia);
         int symAtomIdx = str.find (rotTau, grid)
                    - str.atomInfo->offset(is);
         SxComplex16 phase = exp(I * (kVec ^ (rotTau - str(is,symAtomIdx))) );
         if (timeReversal) phase = phase.conj ();

         // --- loop of partial wave types
         for (int ipt = 0; ipt < pot.getNProjType (is); ++ipt)  {
            int l = pot.lPhi(is)(ipt), nm = 2 * l + 1;
            int offset = pot.offset(is)(ipt) + offsetIs;
            int offsetProj = offset + npl * symAtomIdx; // Ref. 2, Eq. (38)
            int offsetRes  = offset + npl * ia;
            const SxMatrix<Double> &Dl  = ylmRot(iSym)(l);
            
            // --- rotate Ylm
            for (int m = 0; m < nm; ++m)  {
               SxComplex16 sum = 0.;
               for (int m2 = 0; m2 < nm; ++m2)
                  sum += Dl(m,m2) * proj(offsetProj + m2); // Ref. 2, Eq. (38)
               if (timeReversal)  {
                  // Ref. 2, Eq. (39)
                  sum = sum.conj ();
                  if (l & 1) sum = -sum; // (i)^l / (i^l)*
               }
               res(offsetRes + m) = sum * phase;
            }
         }
      }
      offsetIs += npl * str.getNAtoms (is);
   }
   VALIDATE_VECTOR (res);
   return res;
}

void SxPAWExchange::computeXij (const SxRadMat &Dij)
{
   // Compute U kernel (done here by convenience)
   computeUKernel ();
   SX_CLOCK (Timer::Xij);
   SX_CHECK (potPtr);
   SX_CHECK (gPtr->structPtr);
   const SxAtomicStructure &str = *gPtr->structPtr;
   const SxPAWPot &pot = *potPtr;
   int nSpecies = str.getNSpecies ();
   int nSpin = Dij.getNSpin ();
   Xij.resize (Dij);
   double spinFactor = (nSpin == 1) ? 0.5 : 1.;
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      for (int is = 0; is < nSpecies; ++is)  {
         int nAtoms = str.getNAtoms (is);
         int npt = pot.getNProjType (is);
         const SxNArray<double, 2> &K = pot.xKernel(is);
         //SX_CHECK (pot.coreX(is).getSize () > 0);

         for (int ia = 0; ia < nAtoms; ++ia)  {
            SxMatrix<Double> &XijR = Xij(iSpin, is, ia);
            const SxMatrix<Double> &DijR = Dij(iSpin, is, ia);
            // --- loop over outer partial waves
            for (int jpt = 0; jpt < npt; ++jpt)  {
               int l2 = pot.lPhi(is)(jpt), off2 = pot.offset(is)(jpt) + l2;
               for (int m2 = -l2; m2 <= l2; ++m2)  {
                  int lm2 = SxYlm::combineLm(l2,m2);
                  for (int ipt = 0; ipt < npt; ++ipt)  {
                     int l1 = pot.lPhi(is)(ipt), off1 = pot.offset(is)(ipt)+l1;
                     for (int m1 = -l1; m1 <= l1; ++m1)  {
                        int lm1 = SxYlm::combineLm(l1,m1);
                        // compute Ref. 2, Eq. (52) for this Xij
                        XijR(off1 + m1, off2 + m2) = computeXijSum (is,DijR,K,
                                                                    ipt,l1,lm1,
                                                                    jpt,l2,lm2)
                                                   * spinFactor
                                                   ;//+ pot.coreX(is)(ipt, jpt);
                     }
                  }
               }
            }
         }
      }
   }
}

double SxPAWExchange::computeXijSum (int is,
                                     const SxMatrix<Double> &DijR,
                                     const SxNArray<double, 2> &K,
                                     int ipt, int l1, int lm1,
                                     int lpt, int l4, int lm4)
{
   const SxPAWPot &pot = *potPtr;
   int npt = pot.getNProjType (is);
   double res = 0.;
   // --- loop over inner partial wave types
   for (int jpt = 0; jpt < npt; ++jpt)  {
      int l2 = pot.lPhi(is)(jpt), off2 = pot.offset(is)(jpt) + l2;
      for (int kpt = 0; kpt < npt; ++kpt)  {
         int l3 = pot.lPhi(is)(kpt), off3 = pot.offset(is)(kpt) + l3;

         // --- possible l from triangular condition & parity

         // sum of all l must be even (exchange vanishes otherwise)
         if ((l1 + l2 + l3 + l4) & 1) continue;

         int minL = max(abs(l1-l2), abs(l3-l4));
         int maxL = min (l2 + l1, l3 + l4);
         maxL = min(maxL, pot.lMaxRho(is));
         if (maxL < minL) continue;

         int ijkl = pot.get4Idx (ipt, jpt, kpt, lpt, npt);

         // --- now the inner index angular momentum sums
         for (int m2 = -l2; m2 <= l2; ++m2)  {
            int lm2 = SxYlm::combineLm(l2,m2);
            for (int m3 = -l3; m3 <= l3; ++m3)  {
               int lm3 = SxYlm::combineLm(l3,m3);
               
               // loop over radial mesh's lm's
               for (int l = minL; l <= maxL; l+=2)  {
                  for (int m = -l; m <= l; ++m)  {
                     int lm = SxYlm::combineLm(l,m);
                     double cg = pot.clebschGordan(lm1, lm2, lm);
                     if (fabs(cg) < 1e-12) continue;
                     cg *= pot.clebschGordan(lm3, lm4, lm);
                     if (fabs(cg) < 1e-12) continue;
                     // i^l factors cancel out
                     // Ref. 2, Eq. (52) 
                     res -= cg * DijR(off2 + m2, off3 + m3) * K(l, ijkl);
                  }
               }
            }
         }
      }
   }
   return res;
}

PsiG SxPAWExchange::computePAWCorrG (const SxDiracVec<TPrecCoeffG> &pa,
                                     const SxDiracVec<TPrecCoeffG> &pb,
                                     int iSym,
                                     const Coord &dkq,
                                     SxComplex16 *Mab0)
{
   *Mab0 = 0.;
   PsiG res(smallG(iSym).getSize ());
   res.set (0.);
   const SxPAWPot &pot = *potPtr;
   ssize_t offset = 0;
   SxArray2<SxComplex16> Mab;
   SX_LOOP2(is,ia)  {
      int maxM = (int)pot.MijL(is).getDim (0);
      Mab.reformat (sqr(pot.lMaxRho(is) + 1), maxM);
      Mab.set (0.);
      SX_LOOP2(ipt,jpt)  {
         int l1 = pot.lPhi(is)(ipt);
         int l2 = pot.lPhi(is)(jpt);
         int ij = pot.get2Idx ((int)ipt, (int)jpt, pot.getNProjType (is));
         for (int m1 = -l1; m1 <= l1; ++m1)  {
            int lm1 = SxYlm::combineLm(l1,m1);
            for (int m2 = -l2; m2 <= l2; ++m2)  {
               int lm2 = SxYlm::combineLm(l2,m2);
               // Ref. 2, Eq. (69) <psi_a|p_i><p_j|psi_b> part
               SxComplex16 pij 
                  = pa(offset + pot.offset(is)(ipt) + l1 + m1).conj ()
                  * pb(offset + pot.offset(is)(jpt) + l2 + m2);

               for (int m = 0; m < maxM; ++m)  {
                  for (int L = abs(l1-l2); L <= pot.lMaxRho(is); L+=2) {
                     for (int M = -L; M <= L; ++M)  {
                        int LM = SxYlm::combineLm(L,M);
                        // Ref. 2, Eq. (69)
                        Mab(LM,m) += pij 
                           * pot.clebschGordan(lm1,lm2,LM)
                           * pot.MijL(is)(m,L,ij);
                     }
                  }
               }
            }
         }
      }
      offset += pot.getNProj ((int)is);
      SxComplex16 kPhase = exp(I * (dkq ^ gPtr->structPtr->getAtom (is, ia)));
      // Ref. 2, Eq. (67), sum over LM and m = 2,4,... 
      SX_LOOP(ig)  {
         ssize_t idxg = smallG(iSym)(ig);
         SxComplex16 phase = gPtr->getPhaseFactors((int)is,(int)ia,idxg)
                           * kPhase;
         double g2 = smallG2(iSym)(ig);
         SX_LOOP(LM)  {
            SxComplex16 Mabm = 0.;
            double g2m = g2;
            for (int m = 0; m < maxM; ++m) {
               Mabm += Mab(LM,m) * g2m;
               g2m *= g2;
            }
            res(ig) += Mabm * YlmGl(iSym)(idxg,ssize_t(LM)) * phase;
         }
      }
      *Mab0 += Mab(0,0); // needed for 2nd-order term in V
   }
   return res;
}

#include <SxCLI.h>
#include <SxFermi.h>
#include <SxPAWHamiltonian.h>

void SxPAWExchange::test (int argc, char **argv)
{
   cout.precision (3);
   SxCLI cli (argc, argv);
   SxString inputFile = cli.option ("-i", "input file", "input file")
                        .toString ("input.sx");
   SxString wavesFile = cli.option ("-w", "waves file", "waves file")
                     .toString ("waves.sxb");
   cli.finalize ();

   SxParser parser;
   SxParser::Table table = parser.read (inputFile);
   SxPtr<SxPAWPot> pawPot = SxPtr<SxPAWPot>::create (&*table);
   potPtr = pawPot;
   pawPot->setupQijL ();

   SxFermi fermi;
   SxAtomicStructure structure;
   SxPW waves;
   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      waves.read (io);
      fermi.read (io);
      //fermi.focc.set (0.);
      //fermi.focc(0,0,0) = 2.;
      structure.read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxGkBasis &gkBasis = waves.getGkBasis();
   wavesPtr = waves.getThis ();
   SxMesh3D mesh = gkBasis(0).fft3d(0).mesh;
   SxRBasis R(mesh, structure.cell);
   SxGBasis G(mesh, structure,
              SxGBasis::getGCut(SxGBasis::getECut(&*table)));
   G.registerRBasis (R);
   R.registerGBasis (G);
   rPtr = R.getThis ();
   gPtr = G.getThis ();
   gkBasis.changeTau (structure);

   focc = fermi.focc;
   SxXCFunctional::omegaHSE = 0.;

   SxPtr<SxRadBasis> radPtr 
      = SxPtr<SxRadBasis>::create(potPtr->rad, potPtr->logDr);
   pawPot->setBasis (radPtr);
   pawPot->computeCoreX ();
   pawPot->computeXKernel ();

   SxPtr<SxPartialWaveBasis> thePBasis = thePBasis.create(pawPot, structure);
   thePBasis->createProjBasis (gkBasis);
   pBasisX = pBasis = thePBasis;

   ylmRot = SxYlm::computeYlmRotMatrices 
      (structure.cell.symGroupPtr->getSymmorphic (), potPtr->getLMax ());

   SxRadMat Dij = SxPAWHamiltonian::computeDij (waves, focc, potPtr,
                                                *pBasis, ylmRot);
   /*
   computeXij (Dij);
   cout << Dij(0)(0,0) << endl;
   cout << Xij(0)(0,0) << endl;

   for (int i = 0; i < waves.getNStates (0); ++i)  {
      PsiG res = apply (waves(i,0,0));
      for (int j = 0; j < waves.getNStates (0); ++j)  {
         cout << "<" << j << "00|X|" << i << "00> = " 
              << dot(waves(j,0,0), res) << endl;
      }
   }
   */

   compute (waves.getThis (), focc, Dij);

   printTiming ();

}

