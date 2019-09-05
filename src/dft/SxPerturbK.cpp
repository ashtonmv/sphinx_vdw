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

#include <SxPerturbK.h>
#include <SxRadBasis.h>
#include <SxProjector.h>
#include <SxProjMatrix.h>
#include <SxYlm.h>

enum {PHASEFACTORS, PHI_NL, PHI_NLR, EKB,
      KINETIC, PSIPHI, SUM_PHI,
      TOTAL_TIME, N_TIMER};

SxPerturbK::SxPerturbK ()
{
   timer.init (N_TIMER);
   timer.setName (PHI_NL, "phiNl setup");
   timer.setName (PHI_NLR, "phiNlR setup");
   timer.setName (EKB, "eKB setup");
   timer.setName (KINETIC, "kinetic part");
   timer.setName (PHASEFACTORS, "phase factors");
   timer.setName (PSIPHI, "<psi|phi>");
   timer.setName (SUM_PHI, "sum_phi");
   timer.setName (TOTAL_TIME, "total");
   nonLocal = false;
   blockSize = 64;
}

void SxPerturbK::set (const SxPseudoPot &psPot, const SxGkBasis &gk,
                   const SxAtomicStructure &structure)
{
   timer.start (TOTAL_TIME);
   (cout << "Setting up (phi) and (ir * phi)..." << endl).flush ();
   int nk = gk.getNk ();
   gBases.resize (nk);
   phiNl.resize (nk);
   phiNlR.resize (nk);

   SxRadBasis r(psPot.rad, psPot.logDr);

   int nSpecies = structure.getNSpecies ();
   SX_CHECK (nSpecies == psPot.getNSpecies (),
             nSpecies, psPot.getNSpecies ());

   // --- determine nOrb, lmax
   int nOrb = 0;
   int nAllOrb = 0;
   eKB.resize(nSpecies);
   int is, iOrb, l, idir, ik;
   int lmax = -1;
   for (is = 0; is < nSpecies; is++)  {
      eKB(is).resize (psPot.lMax(is) + 1);
      for (l = 0; l <= psPot.lMax(is); l++)
         if (!(psPot.lLoc(is) == l))  {
            eKB(is)(l).resize (/* n = */1 );
            nOrb += 2*l + 1;
            nAllOrb += (2*l + 1) * structure.getNAtoms(is);
            if (lmax < l) lmax = l;
         }
   }

   if (nOrb == 0)  {
      nonLocal = false;
      return;
   }

   // --- set up orbInfo
   orbInfo.resize (nAllOrb);
   {
      int iTlOrb = 0, iOrbStart;
      SxQuantumNumbers orb;
      orb.ik = orb.iSpin = -1;
      orb.n = 0; // multiple projectors to be implemented
      orb.i = 0; // reference orbital
      for (orb.iSpecies = 0; orb.iSpecies < nSpecies; orb.iSpecies++)  {
         iOrbStart = orb.i;
         for (orb.iAtom = 0; orb.iAtom < structure.getNAtoms(orb.iSpecies); 
              orb.iAtom++)
         {
            orb.i = iOrbStart;
            for (orb.l = 0; orb.l <= psPot.lMax(orb.iSpecies); orb.l++)  {
               if (psPot.lLoc(orb.iSpecies) == orb.l) continue;
               for (orb.m = -orb.l; orb.m <= orb.l; orb.m++, orb.i++)
                  orbInfo(iTlOrb++) = orb;
            }
         }
      }
      SX_CHECK (iTlOrb == nAllOrb);
      SX_CHECK (orb.i = nOrb);
   }

   // resize phiNl, phiNlR, setup gBases
   int ng;
   for (ik = 0; ik < nk; ik++)  {
      gBases(ik) = &gk(ik);
      ng = gk(ik).ng;
      phiNl(ik).reformat (ng, nOrb);
      phiNlR(ik).resize(3);
      for (idir = 0; idir < 3; idir++)  {
         phiNlR(ik)(idir).reformat(ng, nOrb);
         phiNlR(ik)(idir).set (0.);
      }
   }

   // vars for calc of phi_nl
   SxDiracVec<Double> psi_pot, rad, psi, pot; // :r
   double ldr, clebsch;

   SxYlm::SxClebschTable ClebschGordan = SxYlm::getClebschGordan(lmax, 1, lmax + 1);

   PsiG psi_potG;
   int m, l2, m2, offset;
   int i1,i2[3],i3;
   i2[0] = SxRadBasis::IPX;
   i2[1] = SxRadBasis::IPY;
   i2[2] = SxRadBasis::IPZ;

   iOrb = 0;
   for (is = 0; is < nSpecies; is++)  {
      for (l = 0; l <= psPot.lMax(is); l++)  {
         if (psPot.lLoc(is) == l) continue;
         
         // --- retrieve pseudopotential information
         rad = psPot.rad(is);
         psi = psPot.pseudoPsi(is)(l);
         pot = psPot.pseudoPot(is)(l);

         // TODO: ugly!!!
         psi_pot = psi * pot;
         psi_pot.setBasis (&r);
         psi_pot.handle->auxData.is = is;
         psi_pot.handle->auxData.l  = l;
         ldr = psPot.logDr(is);
         
         // Kleinman-Bylander energies
         timer.start (EKB);
         eKB(is)(l)(0 /*n*/) = (psi.sqr () * rad.cub () * pot).integrate (ldr);
         timer.stop (EKB);

         offset = iOrb + l; // offset for m = 0
         
         timer.start (PHI_NL);
         for (m=-l; m <=l; m++, iOrb++)  {
            psi_pot.handle->auxData.m  = m;
            for (ik = 0; ik < nk; ik++)  {
               (cout << '.').flush ();
               phiNl(ik).colRef(iOrb) <<= gk(ik) | psi_pot;
            }
         } // m
         timer.stop (PHI_NL);

         psi_pot *= rad; // must ensure that l2 <> l
         /*
         psi_pot = psi * pot * rad;
         psi_pot.setBasis (&r);
         psi_pot.handle->auxData.is = is;
         */
         psi_pot *= sqrt(FOUR_PI / 3.);  // normalization factor of Y1m

         timer.start (PHI_NLR);
         for (l2 = l-1; l2 <= l+1; l2 += 2)  {
            if (l2 < 0) continue;
            psi_pot.handle->auxData.l  = l2;
            for (m2 = -l2; m2 <= l2; m2++)  {
               i3 = SxYlm::combineLm (l2,m2);
               psi_pot.handle->auxData.m  = m2;
               for (ik = 0; ik < nk; ik++)  {
                  (cout << '.').flush ();
                  psi_potG = gk(ik) | psi_pot;
                  for (m = -l; m <= l; m++)  {
                     i1 = SxYlm::combineLm(l,m);
                     for (idir = 0; idir < 3; idir++)  {
                        clebsch = ClebschGordan(i1, i2[idir], i3);
                        if (fabs(clebsch) > 1e-10)
                           phiNlR(ik)(idir).colRef(offset + m)
                              .plus_assign_ax(clebsch, psi_potG);
                     } // idir
                  } // m
               } // ik
            } // m2
         } // l2
         timer.stop (PHI_NLR);
      }  // l
   } // is
   cout << endl;

   nonLocal = true;
   timer.stop (TOTAL_TIME);

}

void SxPerturbK::printTimer ()
{
   timer.print (TOTAL_TIME);
}

SxDiracMat<TPrecCoeffG> 
SxPerturbK::getMatrixElements (const SxGBasis::TPsi &psiL, 
                            const SxGBasis::TPsi &psiR) const
{
   timer.start(TOTAL_TIME);
   const SxGBasis *gBasis 
      = dynamic_cast <const SxGBasis *> (psiL.getBasisPtr ());
   SX_CHECK (gBasis);
   SX_CHECK (gBasis == dynamic_cast <const SxGBasis *> (psiR.getBasisPtr ()));

   int ik, idir;

   int nL = (int)psiL.nCols ();
   SX_CHECK (nL > 0, nL);
   int nR = (int)psiR.nCols ();
   SX_CHECK (nR > 0, nR);

   SxDiracMat<TPrecCoeffG> result(nL * nR, 3);
   SxDiracVec<TPrecCoeffG> resCol;
   // --- kinetic part
   {
      timer.start (KINETIC);
      SxDiracVec<TPrecCoeffG> gPsi;
      SxDiracVec<TPrecCoeffG>::Iterator gPsiIt, psiIt;
      SxDiracVec<Double>::Iterator gVecIt;
      int iR, ig, ng = (int)psiR.nRows ();
      for (idir = 0; idir < 3; idir++)  {
         // get result matrix
         resCol = result.colRef (idir);
         resCol.reshape (nL, nR);

         /*
         for (int iR = 0; iR < nR; iR++)  {
            gPsi = gBasis->gVec.colRef(idir) * psiR.colRef(iR);
            // this is inefficient, but will improve a lot
            // with the new XPress library
            // resCol.colRef(iR) <<= (psiL.adjoint () ^ gPsi);
            if (nL == 1)
               resCol(iR) = (psiL ^ gPsi).chop ();
            else
               resCol.colRef(iR) <<= (gPsi ^ psiL.conj ());
         }
         */
         // gPsi = PsiG (); gPsi.copy(psiR);
         gPsi.reformat (ng, nR);
         gPsiIt = gPsi.begin ();
         psiIt = psiR.begin ();
         for (iR = 0; iR < nR; iR++)  {
            // gPsi.colRef(iR) *= gBasis->gVec.colRef(idir);
            gVecIt = gBasis->gVec.colRef(idir).begin ();
            for (ig = 0; ig < ng; ig++, ++psiIt, ++gPsiIt, ++gVecIt)
               *gPsiIt = (*gVecIt) * (*psiIt);
         }
         VALIDATE_VECTOR (gPsi);
         VALIDATE_VECTOR (psiL);
         resCol <<= psiL.adjoint () ^ gPsi;


      } // idir
      timer.stop (KINETIC);
   } // kinetic part brace

   if (!nonLocal)  {
      timer.stop (TOTAL_TIME);
      return result;
   }

   // --- determine ik
   int nk = int(gBases.getSize ());
   for (ik = 0; ik < nk; ik++)
      if (gBases(ik) == gBasis) break;
   SX_CHECK (ik < nk); // unknown basis if ik == nk

   // --- nonlocal part
   /*   The projections onto the nonlocal pseudopotential projectors
        and the derivative projectors uses a block-algorithm (SxProjMatrix).
        Here we implement the projector setup.
    */
   class SxNonLocalPart : public SxProjMatrix<TPrecCoeffG>::SaveProjections
   {
      private:
         // reference (unshifted) projectors
         const PsiGI             &phi;     // :(ig:iOrb)
         const SxArray<PsiGI>    &phiR;    // :idir:(ig:iOrb)
         // the G-basis (provides phase factors)
         const SxGBasis          &g;
         // cached phase factor
         mutable PsiG T;
         mutable SxQuantumNumbers lastOrb;
         // reference to SxPerturbK class
         const SxPerturbK &kp;
      public:
         // --- the usual blabla: constructor & destructor
         virtual ~SxNonLocalPart () {}
         SxNonLocalPart (const PsiGI &phi_, const SxArray<PsiG> &phiR_,
                         const SxGBasis &g_,
                         const SxPerturbK *kp_) 
            : SxProjMatrix<TPrecCoeffG>::SaveProjections (
                  4 * int(kp_->orbInfo.getSize ()),
                  kp_->blockSize,
                  (int)phi_.nRows ()),
              phi(phi_), phiR(phiR_), g(g_), kp(*kp_)
         { }

         // --- pseudo-define unnecessary virtual functions
         HAS_TARGET_GETPROJECTOR;
         NO_GETFACTOR;

         // --- the real stuff: get a projector. We have (3 + 1) * nOrb
         //     projectors, for each orbital 3 derivative projectors and
         //     one normal. In order to make use of cached phase factors
         //     we make iOrb the slow index and (3+1) the fast.
         virtual void getProjector(int i, SxDiracVec<TPrecCoeffG> *res) const
         {
            int iPhi = i / 4; // slow index: projector id
            int dir  = i % 4; // fast index: direction (0..2) and normal (3)
            SxQuantumNumbers orb = kp.orbInfo(iPhi);
            // --- get new phase factor (if necessary)
            if (   lastOrb.iAtom    != orb.iAtom
                || lastOrb.iSpecies != orb.iSpecies)
            {
               kp.timer.start(PHASEFACTORS);
               T = g.getPhaseFactors(orb.iSpecies, orb.iAtom);
               kp.timer.stop(PHASEFACTORS);
            }
            // --- multiply reference orbital with phase factor
            PsiG::Iterator resIt = res->begin (), 
                           itT = T.begin (), 
                           phiIt;
            if (dir == 3)  {
               double invEKB = 1. / kp.eKB(orb.iSpecies)(orb.l)(orb.n);
               phiIt = phi.colRef(orb.i).begin ();
               for (int ig = 0; ig < nElements; ++ig)
                  *resIt++ = *itT++ * *phiIt++ * invEKB;
            } else {
               phiIt = phiR(dir).colRef(orb.i).begin ();
               for (int ig = 0; ig < nElements; ++ig)
                  *resIt++ = *itT++ * *phiIt++;
            }
            lastOrb = orb;
         }

   /** Don't get confused. This here is an embedded class with a single
       instantiation, so we can declare the variable nlpp together with
       the class SxNonLocalPart.
    */
   } nlpp(phiNl(ik), phiNlR(ik), *gBasis, this);
   

   // --- get projections
   int nOrb = int(orbInfo.getSize ());
   SxDiracVec<TPrecCoeffG> phiPsiL, phiPsiR, lSide, rSide;

   timer.start(PSIPHI);
   
   phiPsiL = nlpp.getProjection (psiL); // (4 x nOrb, nL)
   phiPsiL.reshape (4, nOrb * nL); // (4, nOrb x nL)
   phiPsiL = phiPsiL.transpose (); // -> (nOrb x nL, 4)

   phiPsiR = nlpp.getProjection (psiR); // (4 x nOrb, nR)
   phiPsiR.reshape (4, nOrb * nR); // (4, nOrb x nR)
   phiPsiR = phiPsiR.transpose (); // -> (nOrb x nR, 4)

   // --- precalculate adjoint for phiPsiL(3)
   lSide = phiPsiL.colRef(3);
   lSide.reshape (nOrb, nL);
   lSide <<= lSide.adjoint ();

   timer.stop(PSIPHI);

   // --- sum over projectors
   timer.start (SUM_PHI);
   for (idir = 0; idir < 3; idir++)  {
      // get result matrix
      resCol = result.colRef(idir);
      resCol.reshape (nL, nR);

      // --- <psiL|ir * phi><phi|psiR>
      lSide = phiPsiL.colRef (idir);
      lSide.reshape (nOrb, nL);
      rSide = phiPsiR.colRef (3);
      rSide.reshape (nOrb,nR);

      resCol -= lSide.adjoint () ^ rSide;

      // --- <psiL|phi><ir * phi|psiR>
      rSide = phiPsiR.colRef (idir);
      rSide.reshape (nOrb, nR);
      lSide = phiPsiL.colRef (3);
      // lSide.reshape (nOrb, nL);
      // resCol -= lSide.adjoint () ^ rSide;
      lSide.reshape (nL, nOrb); // adjoint done above
      resCol -= lSide ^ rSide;
   }
   timer.stop (SUM_PHI);
         
   timer.stop(TOTAL_TIME);

   return result;
}

