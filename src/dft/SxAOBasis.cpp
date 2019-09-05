#include <SxDFTConfig.h>
#include <SxAOBasis.h>
#include <SxProjector.h>
#include <SxProjMatrix.h>
#include <SxPWOverlap.h>
#ifdef USE_SXGEMMM
#include <SxGemmm.h>
#endif
#include <SxPAWBasis.h>

class SxAOBasisProj
   : public SxProjMatrix<TGBasisType>::SaveProjections
{
   public:
      const SxArray<SxAOBasis::OrbitalIndex> &orbitalMap;
      const SxArray<SxDiracMat<TPrecCoeffG> > &refOrbitals;
      mutable PsiG T;
      SxAOBasisProj(const SxAOBasis &ao, int ik)
         : SxProjMatrix<TGBasisType>::SaveProjections (
               ao.getNOrb (),
               ao.blockSize,
               (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).nRows ()),
           orbitalMap(ao.orbitalMap),
           refOrbitals(ao.refOrbitals(ik))
      { /* empty */ }

      virtual ~SxAOBasisProj () {}

      virtual void getProjector(int iOrb, SxDiracVec<TGBasisType> *target) const
      {
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxAOBasis::OrbitalIndex &idx = orbitalMap(iOrb);
         const SxDiracMat<TGBasisType> &refOrb
            = refOrbitals(idx.is).colRef (idx.io);

         if (!T.handle
             || T.handle->auxData.is != idx.is
             || T.handle->auxData.ia != idx.ia)
         {
            T = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
         }
         SX_CHECK (T.handle->auxData.is == idx.is,
                   T.handle->auxData.is, idx.is);
         SX_CHECK (T.handle->auxData.ia == idx.ia,
                   T.handle->auxData.ia, idx.ia);

#ifdef USE_OPENMP
#pragma omp parallel for if (nElements > sxChunkSize)
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) = T(ig) * refOrb(ig);
      }

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
};

SxAOBasis::SxAOBasis ()
{
   cacheRefOrb = Unknown;
   init ();
}

namespace Timer {
   // timer ids
   enum AOTimers{ AoRefOrbitalSetup, AoOverlap, AoSInversion, AoProjection,
                  AoGradient, AoGradProject, AoTotalTime} ;
}

SX_REGISTER_TIMERS(Timer::AOTimers)
{
   using namespace Timer;
   regTimer (AoRefOrbitalSetup,"Ref. orbital setup");
   regTimer (AoOverlap,"Overlap setup");
   regTimer (AoSInversion,"Overlap inversion");
   regTimer (AoProjection,"AO projection");
   regTimer (AoGradient,"AO gradient");
   regTimer (AoGradProject,"AO projection d/dR");
   regTimer (AoTotalTime, "AOBasis total");
}

void SxAOBasis::init ()  {
   cacheOverlap = cacheInverse = Unknown;
   refOrbCachedK = overlapCachedK = invOverlapCachedK = -1;
   blockSize = 64;
}

SxAOBasis::~SxAOBasis ()
{
   deregisterAll ();
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxRadBasis &rad,
                      const SxArray<SxArray<SxDiracVec<TReal8> > > &psiRad,
                      const SxConstPtr<SxOverlapBase> SPtrIN)
{
   init ();
   set (gk, rad, psiRad, SPtrIN);
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxAtomicOrbitals &ao,
                      const SxConstPtr<SxOverlapBase> SPtrIN)
{
   init ();
   set (gk, *ao.getRadBasisPtr (), ao.getMuSet (), SPtrIN);
}

SxAOBasis::SxAOBasis (const SxDiracVec<Complex16> &YlmGl, int iSpecies)
{
   init ();
   SX_CHECK (YlmGl.getBasis<SxGBasis> ().structPtr);
   const SxAtomicStructure &structure = *(YlmGl.getBasis<SxGBasis> ().structPtr);
   int nRef = (int)YlmGl.nCols ();
   int nOrb = nRef * structure.getNAtoms (iSpecies);

   orbitalMap.resize (nOrb);
   for (int ia = 0; ia < structure.getNAtoms (iSpecies); ++ia)  {
      for (int io = 0; io < nRef; ++io)  {
         orbitalMap(io + nRef * ia) = OrbitalIndex(iSpecies, ia, io);
      }
   }

   refOrbitals.resize (1);
   refOrbitals(0).resize (structure.getNSpecies ());
   refOrbitals(0)(iSpecies) = YlmGl;
   cacheRefOrb = CacheAll;

}

void SxAOBasis::set (const SxGkBasis &gk,
                     const SxRadBasis &rad,
                     const SxArray<SxArray<SxDiracVec<TReal8> > > &psiRad,
                     const SxConstPtr<SxOverlapBase> SPtrIN)
{
   SX_CLOCK (Timer::AoTotalTime);
   SPtr = SPtrIN;
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   int nSpecies = int(psiRad.getSize ());
   SX_CHECK (gk.getNk () > 0, gk.getNk ());
   SX_CHECK (nSpecies > 0, nSpecies);
   cacheRefOrb = CacheAll;

   // --- setup orbital map io -> (n,l,m) for each species
   refOrbMap.resize (nSpecies);
   SxVector<Int> nAtoms(nSpecies), nOrbital(nSpecies);

   int nOrb = 0;
   {
      for (int is = 0; is < nSpecies; ++is)  {
         int nOrbType = int(psiRad(is).getSize ());
         int nOrbSpecies = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType)
            nOrbSpecies += 2 * psiRad(is)(iOrbType).handle->auxData.l + 1;
         refOrbMap(is).resize (nOrbSpecies);
         int io = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType) {
            int l = psiRad(is)(iOrbType).handle->auxData.l;
            for (int m = -l; m <= l; ++m, ++io)  {
               refOrbMap(is)(io) = AoIndex(iOrbType, l, m);
            }
         }
         nAtoms(is) = gk.getTau ().getNAtoms(is);
         nOrbital(is) = io;
         nOrb += io * nAtoms(is);
      }
   }

   // -- setup orbital map iOrb -> (is, ia, io)
   orbitalMap.resize (nOrb);
   {
      int iOrb = 0, ia, is, io;
      for (is = 0; is < nSpecies; ++is)
         for (ia = 0; ia < nAtoms(is); ++ia)
            for (io = 0; io < nOrbital(is); ++io, ++iOrb)
               orbitalMap(iOrb) = OrbitalIndex (is, ia, io);
   }

   // --- setup reference orbitals
   computeRefOrbitals(gk, rad, psiRad);
   // clean overlap cache
   overlap.resize (0);
   invOverlap.resize (0);
   overlapCachedK = invOverlapCachedK = -1;
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxArray<SxDiracMat<Double> > &psi,
                      const SxArray<SxArray<int> > &lPsi)
{
   init ();
   SX_CLOCK (Timer::AoTotalTime);
   //SPtr = SPtrIN;
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   int nSpecies = int(psi.getSize ());
   int nk = gk.nk;
   SX_CHECK (nk > 0, nk);
   SX_CHECK (nSpecies > 0, nSpecies);
   cacheRefOrb = CacheAll;
   // clean overlap cache
   overlap.resize (0);
   invOverlap.resize (0);
   overlapCachedK = invOverlapCachedK = -1;

   // --- setup orbital map io -> (n,l,m) for each species
   refOrbMap.resize (nSpecies);
   SxVector<Int> nAtoms(nSpecies), nOrbital(nSpecies);

   int nOrb = 0;
   {
      for (int is = 0; is < nSpecies; ++is)  {
         int nOrbType = (int)psi(is).nCols ();
         int nOrbSpecies = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType)
            nOrbSpecies += 2 * lPsi(is)(iOrbType) + 1;
         refOrbMap(is).resize (nOrbSpecies);
         int io = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType) {
            int l = lPsi(is)(iOrbType);
            for (int m = -l; m <= l; ++m, ++io)  {
               refOrbMap(is)(io) = AoIndex(iOrbType, l, m);
            }
         }
         nAtoms(is) = gk.getTau ().getNAtoms(is);
         nOrbital(is) = io;
         nOrb += io * nAtoms(is);
      }
   }

   // -- setup orbital map iOrb -> (is, ia, io)
   orbitalMap.resize (nOrb);
   {
      int iOrb = 0, ia, is, io;
      for (is = 0; is < nSpecies; ++is)
         for (ia = 0; ia < nAtoms(is); ++ia)
            for (io = 0; io < nOrbital(is); ++io, ++iOrb)
               orbitalMap(iOrb) = OrbitalIndex (is, ia, io);
   }

   // --- setup reference orbitals
   refOrbitals.resize (nk);
   double gMax = 0.;
   for (int ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))  {
         refOrbitals(ik).resize (nSpecies);
         for (int is = 0; is < nSpecies; ++is)  {
            refOrbitals(ik)(is).reformat (gk(ik).ng, nOrbital(is));
            refOrbitals(ik)(is).setBasis (gk(ik));
            refOrbitals(ik)(is).handle->auxData.ik = ik;
         }
         gMax = max(gMax, sqrt(gk(ik).g2(gk(ik).ng - 1)));
      }
   }
   if (gMax == 0.) return; // no k-points

   gMax *= 1.1;
   double cutVol = FOUR_PI / 3. * gMax * gMax * gMax;
   double dg = 0.003;
   int ng1 = int(cutVol / gk.getTau ().cell.getReciprocalCell ().volume);
   int ng2 = int(gMax/dg) + 1;
   bool viaRadG =  (ng1 * nk > ng2);
   SxRadGBasis radG;
   if (viaRadG) radG.set (0., gMax, ng2, SxRadGBasis::Linear);

   for (int is = 0; is < nSpecies; ++is)  {

      const SxRadGBasis *myRadG 
         = dynamic_cast<const SxRadGBasis*>(psi(is).getBasisPtr ());
      for (int n = 0, io = 0; n < lPsi(is).getSize (); ++n)  {
         int l = lPsi(is)(n);
         SxDiracVec<Double> proj = psi(is).colRef (n);
         proj.handle->auxData.is = is;
         proj.handle->auxData.ia = -1;
         proj.handle->auxData.n  = n;
         proj.handle->auxData.l  = l;
         proj.setBasis (psi(is).getBasisPtr ());
         if (myRadG)  {
            proj = myRadG->toSpline (proj);
         } else if (viaRadG)  {
            proj = radG.toSpline (radG | proj);
         }
         for (int m = -l; m <= l; ++m, ++io)  {
            proj.handle->auxData.m = m;
            for (int ik = 0; ik < nk; ++ik)  {
               if (refOrbitals(ik).getSize () > 0) {
                  refOrbitals(ik)(is).colRef (io) <<= gk(ik) | proj;
               }
            }
         }
      }
   }

}

void SxAOBasis::computeRefOrbitals
   (const SxGkBasis                              &gk,
    const SxRadBasis                             &rad,
    const SxArray<SxArray<SxDiracVec<TReal8> > > &psiRad)
{
   SX_CLOCK (Timer::AoRefOrbitalSetup);
   SxDiracVec<TReal8> psiRef;
   ssize_t nSpecies = psiRad.getSize ();
   ssize_t nk = gk.getNk ();

   double gMax = 0.;
   refOrbitals.resize (nk);
   for (int ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))  {
         refOrbitals(ik).resize(nSpecies);
         SX_LOOP(is)  {
            ssize_t nOrbital = refOrbMap(is).getSize ();
            if (nOrbital > 0)  {
               refOrbitals(ik)(is).reformat (gk(ik).ng, nOrbital);
               refOrbitals(ik)(is).setBasis (&gk(ik));
            }
         }
         gMax = max(gMax, sqrt(gk(ik).g2(gk(ik).ng - 1)));
      }
   }
   if (gMax == 0.) return; // no k-points

   double dg = 0.003;
   int nPoints = int(gMax/dg) + 1;
   SxRadGBasis gRad(0.0, gMax, nPoints, SxRadGBasis::Linear);

   SX_LOOP(is)  {
      SxDiracVec<Double> psiGRef;
      if (refOrbMap(is).getSize () == 0) continue;
      SX_LOOP(io) {
         AoIndex nlm = refOrbMap(is)(io);
         psiRef = psiRad(is)(nlm.n);
         psiRef.setBasis (&rad);
         psiRef.handle->auxData.is = (int)is;
         psiRef.handle->auxData.n = nlm.n;
         psiRef.handle->auxData.l = nlm.l;
         psiRef.handle->auxData.m = nlm.m;
         if (psiGRef.getSize () == 0
             || psiGRef.handle->auxData.n != nlm.n
             || psiGRef.handle->auxData.l != nlm.l)
         {
            psiGRef = gRad.toSpline (gRad | psiRef);
         }
         psiGRef.handle->auxData.m = nlm.m;
         // store <G+k|mu> = sum(r) <G+k|r><r|mu>
         //(*refOrbPtr)(is).colRef(io) <<= ( g | psiRef );
         // store <G+k|mu> = sum(r,gRad) <G+k|gRad><gRad|r><r|mu>
         SX_LOOP(ik)
            if (refOrbitals(ik).getSize () > 0)
               refOrbitals(ik)(is).colRef(io) <<= gk(ik) | psiGRef;
      }
   }
}

SxDiracMat<TAOBasisType> SxAOBasis::calculateOverlap (int ik) const
{
   int iOrb;
#ifndef NDEBUG
   if (cacheRefOrb == CacheAll)  {
      int nk = int(refOrbitals.getSize ());
      SX_CHECK (nk > 0);
      SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   } else if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
   }
#endif
   int nOrb = int(orbitalMap.getSize ());

   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoOverlap);
   SxDiracMat<TAOBasisType> ovlp(nOrb, nOrb);
   TPsi aoMu, aoNu;
   TAOBasisType::Type SMuNu;
   SxOverlap S(SPtr);
   //SxAOBasisProj aop(*this, ik);
   // --- high-mem
   /*
   const SxGBasis *gBasis = dynamic_cast<const SxGBasis*>
                           (refOrbitals(ik)(0).getBasisPtr ());
   SX_CHECK (gBasis);
   int ng = gBasis->ng;
   aoMu.reformat(ng, nOrb);
   for (iOrb = 0; iOrb < nOrb; ++iOrb)
      aoMu.colRef(iOrb) <<= S | getAOinG(ik, iOrb);
   ovlp = aop.getProjection(aoMu);
   */

   // ---medium-mem
   ovlp.reformat(nOrb, nOrb);
   // todo: blocking on iOrb
   for (iOrb = 0; iOrb < nOrb; ++iOrb)  {
      //ovlp.colRef(iOrb) <<= aop.getProjection(S | getAOinG(ik, iOrb));
      ovlp.colRef(iOrb) <<= fromPWBasis (getAOinG(ik, iOrb));
   }

   /*
   // --- low-mem
   for (iOrb = 0; iOrb < nOrb; ++iOrb)  {
      aoMu = getAOinG (ik, iOrb);
      for (jOrb = iOrb; jOrb < nOrb; ++jOrb)  {
         aoNu = S | getAOinG(ik, jOrb);
         SMuNu = dot(aoMu,aoNu);
         ovlp(iOrb,jOrb) = SMuNu;
         ovlp(jOrb,iOrb) = SMuNu.conj ();
      }
   }
   */
   return ovlp;
}


SxDiracVec<SxAOBasis::TBasisType>
SxAOBasis::identity  (const SxAOBasis *basis,
                      const SxDiracVec<TBasisType> &psiAO) const
{
   SX_CHECK (basis == this);
   return psiAO;
}

SxDiracVec<TGBasisType>
SxAOBasis::toPWBasis (const SxGBasis *gPtr,
                      const SxDiracVec<TBasisType> &psiAO) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradient);
   // no projection of empty vectors
   SX_CHECK (psiAO.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiAO.handle->auxData.ik;
   int nk = int(refOrbitals.getSize ());
   if (ik >= 0 && ik < nk)  {
#ifndef NDEBUG
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (gPtr == refOrbitals(ik)(is).getBasisPtr ());
#endif
   } else if (cacheRefOrb == CacheAll)  {
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik).getSize () == 0) continue; // different MPI task
         if (refOrbitals(ik)(is).getBasisPtr () == gPtr)
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }


   int nPsi = (int)psiAO.nCols ();
#ifndef USE_SXGEMMM
   /* This is the BLAS algorithm
      It collects the projectors and uses zgemm calls to do
      the data reduction.
      sxgemm3m is an alternative direct algorithm, see below.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   if (nPsi > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);
      SxDiracVec<TBasisType> psiAOphase;
      if (extraPhase.getSize () > 0)  {
         const SxDiracVec<Complex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         psiAOphase.reformat (getNElements (), nPsi);
         for (ssize_t iState = 0, i = 0; iState < nPsi; ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++i)
               psiAOphase(i) = psiAO(i) * phase(iOrb);
      } else {
         psiAOphase = psiAO;
      }
      SxDiracMat<TPrecCoeffG> res = aop.gradient (psiAOphase);
      res.handle->auxData = psiAO.handle->auxData;
      res.setBasis (gPtr);
      res.handle->auxData.ik = ik;
      return res;
   } else if (nPsi > 1)  {
      // use algorithm below for each psi
      SxDiracMat<Complex16> res (gPtr->ng, nPsi);
      res.handle->auxData = psiAO.handle->auxData;
      res.setBasis (gPtr);
      res.handle->auxData.ik = ik;

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= toPWBasis (gPtr, psiAO.colRef(i));
      return res;
   }
#endif

   /* --- Use atom-blocking
      For a single psi, we follow this strategy:

      For
      \sum_{is,ia,ipl} T(is,ia) |p(is,ipl>  A(is, ipl, ia)
      we first compute
      |pA(is,ia)> = \sum_{ipl} |p(is,ipl)> A(is, ipl, ia)
      as a matrix operation.

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of pA.

      --- Direct algorithm (sxgemmm algorithm)
      We also have a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia).conj () * psiAO(ipl, ia, iPsi)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
   */
   int nOrb = getNOrb ();
   const SxArray<SxDiracMat<TPrecCoeffG> > &proj = refOrbitals(ik);

   PsiGI res;
   res.reformat (gPtr->ng, nPsi);
   res.set (0.);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      // --- find number of projectors for this species
      int npl = (int)proj(is).nCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }
      iOrb = iOrbStart;

#ifndef USE_SXGEMMM
      // --- compute gradient (with blocking on atoms)
      if (npl >= nAtoms)  {
         SxDiracMat<TPrecCoeffG> pA;
         for (int ia = 0, iab = 0; ia < nAtoms; ++ia)  {
            if (iab == pA.nCols ())  {
               iab = 0;
               // discard old block
               pA = SxDiracMat<TPrecCoeffG> ();
               // --- compute new pA block
               int nAtomBlock = min(nAtoms-ia, blockSize);
               SX_CHECK (nAtomBlock > 0, nAtomBlock);
               SxDiracVec<Complex16> psiAOphase;
               if (extraPhase.getSize () > 0)  {
                  const SxDiracVec<Complex16> &phase = extraPhase(ik);
                  psiAOphase.reformat (npl, nAtomBlock);
                  for (ssize_t i = 0; i < psiAOphase.getSize (); ++i)
                     psiAOphase(i) = psiAO(iOrb + i) * phase(iOrb+i);
               } else {
                  SxIdx blockIdx(iOrb, iOrb + nAtomBlock * npl - 1);
                  psiAOphase = psiAO(blockIdx).reshape (npl, nAtomBlock);
               }
               pA = proj(is) ^ psiAOphase;
               iOrb += nAtomBlock * npl;
            }
            // apply translation
            //res += gk.getPhaseFactors (is, ia) * pA.colRef(iab++);
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            PsiG T = gPtr->getPhaseFactors (is, iAtom),
                 pAcol = pA.colRef(iab++);
#ifdef USE_OPENMP
#pragma omp parallel for if (res.getSize () > sxChunkSize)
#endif
            for (int ig = 0; ig < res.getSize (); ++ig)
               res(ig) += T(ig) * pAcol(ig);
         }
         SX_CHECK (iOrbStart + nAtoms * npl == iOrb,
                   iOrb - iOrbStart, nAtoms, npl);
      } else {
         SxDiracMat<TPrecCoeffG> phaseNl(gPtr->ng, npl);
         phaseNl.set (0.);
#else /* USE_SXGEMMM */
      {
#endif
         SxDiracMat<TPrecCoeffG> T;
         for (int ia = 0, iab = 0; /* done inside */; ++ia)  {
            if (iab == T.nCols ())  {
               if (iab > 0)  {
                  SxDiracVec<Complex16> psiAOphase;
                  if (extraPhase.getSize () > 0)  {
                     const SxDiracVec<Complex16> &phase = extraPhase(ik);
#ifdef USE_SXGEMMM
                     // --- multiply psiAO by orbital-dependent phase
                     ssize_t nOrbBlock = iab * npl;
                     psiAOphase.reformat (nOrbBlock, nPsi);
                     for (ssize_t iPsi = 0, iPh = 0; iPsi < nPsi; ++iPsi)  {
                        ssize_t iAO = iOrb + iPsi * nOrb;
                        for (ssize_t jOrb = 0; jOrb < nOrbBlock; ++jOrb)  {
                           psiAOphase(iPh++) = psiAO(iAO++) * phase(iOrb+jOrb);
                        }
                     }
                     // do sxgemm summation
                     sxpgemm3m(nPsi, iab, npl, gPtr->ng,
                               T.elements, proj(is).real ().elements,
                               psiAOphase.elements, npl, nOrbBlock,
                               res.elements);
#else
                     // --- multiply psiAO by orbital-dependent phase
                     //     and reorder to (iAtom, ipl)
                     psiAOphase.reformat (iab, npl);
                     for (ssize_t jab = 0; jab < iab; ++jab)
                        for (ssize_t ipl = 0; ipl < npl; ++ipl)
                           psiAOphase(jab, ipl) = psiAO(iOrb + ipl + jab * npl)
                                                * phase(iOrb + ipl + jab * npl);
                     phaseNl += T ^ psiAOphase;
#endif
                  } else {
#ifdef USE_SXGEMMM
                     sxpgemm3m(nPsi, iab, npl, gPtr->ng,
                               T.elements, proj(is).real ().elements,
                               psiAO.elements + iOrb, npl, nOrb,
                               res.elements);
#else
                     // sum over atoms
                     SxIdx blockIdx(iOrb, iOrb + iab * npl - 1);
                     phaseNl += T ^ psiAO(blockIdx).reshape (npl, iab)
                                                   .transpose ();
#endif
                  }
                  iOrb += iab * npl;
               }
               if (ia == nAtoms) break;
               iab = 0;
               int nAtomBlock = min(nAtoms-ia, blockSize);
               SX_CHECK (nAtomBlock > 0, nAtomBlock);
               // discard old block
               T = SxDiracMat<TPrecCoeffG> ();
               T = SxDiracMat<TPrecCoeffG> (gPtr->ng, nAtomBlock);
            }
            // collect phase factors
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            T.colRef (iab++) <<= gPtr->getPhaseFactors(is, iAtom);
         }

#ifndef USE_SXGEMMM
         // --- sum over local projectors
         int ng = gPtr->ng;
         /*
         for (int ipl = 0; ipl < npl; ++ipl)
            res += proj(is).colRef(ipl) * phaseNl.colRef(ipl);
         */
         for (int ipl = 0; ipl < npl; ++ipl)  {
            PsiG phiRef = proj(is).colRef(ipl);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
            for (int ig = 0; ig < ng; ig++)
               res(ig) += phiRef(ig) * phaseNl(ig + ng * ipl);
         }
         /*
         int gBlock = 64;
         int ig = 0, ig0;
         for (int igstop = 0; ; igstop += gBlock)  {
            if (igstop > ng) igstop = ng;
            ig0 = ig;
            for (int ipl = 0; ipl < npl; ++ipl)  {
               const PsiG &phiRef = proj(is)(ipl);
               PrecCoeffG *srcPtr = &phaseNl(ig0, ipl);
               for (ig = ig0; ig < igstop; ++ig)
                  //res(ig) += phiRef(ig) * phaseNl(ig, ipl);
                  res(ig) += phiRef(ig) * *srcPtr++;
            }
            if (ig == ng) break;
         }
         */
#endif
      }
   }

   VALIDATE_VECTOR(res);

   res.handle->auxData = psiAO.handle->auxData;
   res.handle->auxData.ik = ik;
   res.setBasis (gPtr);
   return res;
}

SxAOBasis::TPsi
SxAOBasis::fromPWBasis (const SxDiracVec<TBasisType> &psiG) const
{
   // no projection of empty vectors
   SX_CHECK (psiG.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiG.handle->auxData.ik;
   int nk = int(refOrbitals.getSize ());
   int is = orbitalMap(0).is; // first species, may be > 0
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiG.getBasisPtr () == refOrbitals(ik)(is).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(is).getBasisPtr () == psiG.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   PsiG psiS = SPtr ? SPtr->apply (psiG) : psiG;
   return fromPWBasis(psiS, ik);
}

SxAOBasis::TPsi
SxAOBasis::fromPWBasis (const PsiG &psiS, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoProjection);

#ifndef USE_SXGEMMM
   /* This is the BLAS algorithm
      It collects the projectors and uses zgemm calls to do
      the data reduction.
      sxgemmm is an alternative direct algorithm, see below.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   if (psiS.nCols () > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);

      TPsi result = aop.getProjectionFromExtended(psiS);

      VALIDATE_VECTOR (result);
      result.handle->auxData = psiS.handle->auxData;
      result.setBasis (this);
      if (extraPhase.getSize () > 0)  {
         const SxDiracVec<Complex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         for (ssize_t iState = 0, iRes = 0; iState < result.nCols (); ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
               result(iRes) *= phase(iOrb).conj ();
      }
      return result;
   } else if (psiS.nCols () > 1)  {
      int nPsi = (int)psiS.nCols ();
      // use algorithm below for each psi
      SxDiracMat<Complex16> res (getNOrb (), nPsi);

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= fromPWBasis (psiS.colRef(i));
      res.handle->auxData = psiS.handle->auxData;
      res.setBasis (this);
      return res;
   }
#endif

   /* --- Use atom-blocking (BLAS algorithm)
      For a single psi, we follow this strategy:

      First note that
      <p(is,ipl) * T(is,ia)|psi> = <p(is,ipl) | T*(is,ia) * psi>

      We now perform the right-hand formula as matrix operation for
      each species:
      <T * p|psi> (npl * nAtoms) = p(npl x ng) ^ psiT(ng x nAtoms)

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of psiT.

      --- Direct algorithm (sxgemmm algorithm)
      We also have a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia).conj () * psi(ig, iState)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
      
   */
   const SxArray<SxDiracMat<TPrecCoeffG> > &proj = refOrbitals(ik);
   const SxGBasis *gk = dynamic_cast<const SxGBasis*>(psiS.getBasisPtr ());
   // get G basis from PAW basis
   if (!gk) gk = psiS.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gk);
   int nOrb = getNOrb ();
   ssize_t nStates = psiS.nCols ();
   TPsi res(nOrb * nStates);
   res.reshape (nOrb, nStates);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;
      SX_CHECK (gk == proj(is).getBasisPtr ());

      // number of projectors for this species
      int npl = (int)proj(is).nCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }
      iOrb = iOrbStart;

#ifndef USE_SXGEMMM
      // --- shortcut for nAtoms = 1
      if (nAtoms == 1)  {
         PsiG psiTrans = gk->getPhaseFactors(is, orbitalMap(iOrb).ia).conj()
                       * psiS;
         res (SxIdx(iOrb,iOrb+npl-1)) <<= proj(is).overlap (psiTrans);
         iOrb += npl;
         continue;
      }
#endif

      // --- perform projections (with blocking on atoms)
      SxDiracMat<TPrecCoeffG> atomBlock;
      int ng = gk->ng;
      for (int ia = 0, iab = 0; ia < nAtoms; ++ia)  {
         if (iab == 0)  {
            // start new (psi * T) block
            atomBlock.reformat (gk->ng, min(nAtoms-ia, blockSize));
         }
         {
            // collect into atomBlock:
            // - for BLAS-based algorithm: psi * T.conj ()
            // - for sxgemmm-based case: T.conj ();
            //psiTrans.colRef (iab++)<<= psi * gk.getPhaseFactors(is, ia).conj ();
            // get the real atom index
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            PsiG bcol = atomBlock.colRef(iab), T = gk->getPhaseFactors(is, iAtom);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
            for (int ig = 0; ig < ng; ig++)
#ifdef USE_SXGEMMM
               bcol(ig) = T(ig).conj ();
#else
               bcol(ig) = psiS(ig) * T(ig).conj ();
#endif
            iab++;
         }

         if (iab == atomBlock.nCols ())  {
            // project current block
            SX_CHECK (proj(is).imag ().normSqr () < 1e-20);
#ifdef USE_SXGEMMM
            sxgemmm(nStates, iab, npl, ng,
                    psiS.elements, psiS.nRows (), atomBlock.elements,
                    proj(is).real ().elements,
                    res.elements + iOrb, npl, nOrb);
#else
            SxIdx currentOrbs(iOrb, iOrb + iab*npl - 1);
            res(currentOrbs) <<= proj(is).overlap (atomBlock);
#endif
            iOrb += iab * npl;
            iab = 0;
            atomBlock = SxDiracMat<TPrecCoeffG> ();
         }
      }
      SX_CHECK (iOrbStart + nAtoms * npl == iOrb,
                iOrb - iOrbStart, nAtoms, npl);
   }
   VALIDATE_VECTOR(res);
   res.handle->auxData = psiS.handle->auxData;
   res.setBasis (this);
   if (extraPhase.getSize () > 0)  {
      const SxDiracVec<Complex16> &phase = extraPhase(ik);
      SX_CHECK (phase.getSize () == getNElements (),
                phase.getSize (), getNElements ());
      for (ssize_t iState = 0, iRes = 0; iState < res.nCols (); ++iState)
         for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
            res(iRes) *= phase(iOrb).conj ();
   }
   return res;
}

int SxAOBasis::getLMax () const
{
   int result = 0;
   int nOrbs = getNOrb ();
   for (int iOrb = 0; iOrb < nOrbs; iOrb++)  {
      int is = orbitalMap(iOrb).is;
      int io = orbitalMap(iOrb).io;
      int l = refOrbMap(is)(io).l;
      if (result < l) result = l;
   }

   return result;
}

SxDiracMat<TGBasisType> SxAOBasis::getAOinG (int ik) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   const SxGBasis *gPtr 
      = dynamic_cast<const SxGBasis *> (refOrbitals(ik)(0).getBasisPtr());
   ssize_t ng = gPtr->g2.getSize ();
   ssize_t nOrbitals = orbitalMap.getSize ();
   SxDiracMat<TGBasisType> result(ng,nOrbitals);
   for (int iOrb = 0; iOrb < nOrbitals; iOrb++)  {
     result.colRef(iOrb) << getAOinG (ik,iOrb);
   }
   result.setBasis(gPtr);
   result.handle->auxData.ik = ik;

   return result; 
}

SxDiracVec<TGBasisType> SxAOBasis::getAOinG (int ik, int iOrb) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   int trueK = ik;
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
             iOrb, orbitalMap.getSize ());
   const OrbitalIndex &idx = orbitalMap(iOrb);
   const SxDiracMat<TGBasisType> &refOrb = refOrbitals(ik)(idx.is);

   const SxGBasis *gPtr
      = dynamic_cast<const SxGBasis *> (refOrb.getBasisPtr ());
   SX_CHECK (gPtr);
   SxDiracVec<TGBasisType> res;
   res = gPtr->getPhaseFactors(idx.is,idx.ia) * refOrb.colRef(idx.io);
   if (extraPhase.getSize () > 0) res *= extraPhase(ik)(iOrb);
   const AoIndex &nlm = refOrbMap(idx.is)(idx.io);
   res.handle->auxData.ik = trueK;
   res.handle->auxData.is = idx.is;
   res.handle->auxData.ia = idx.ia;
   res.handle->auxData.n = nlm.n;
   res.handle->auxData.l = nlm.l;
   res.handle->auxData.m = nlm.m;

   VALIDATE_VECTOR (res);
   return res;
}

SxDiracMat<TAOBasisType> SxAOBasis::getOverlap (int ik) const
{
   // --- single k caching
   if (cacheOverlap == CacheCurrentK)  {
      if (ik != overlapCachedK)  {
         overlapCachedK = ik;
         overlap(0) = SxDiracMat<TAOBasisType> ();
         overlap(0) = calculateOverlap(ik);
      }
      return overlap(0);
   }
   // --- all k caching
   if (cacheOverlap == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < overlap.getSize (),
                ik, overlap.getSize ());
      if (overlap(ik).getSize () == 0)
         overlap(ik) = calculateOverlap(ik);
      return overlap(ik);
   }
   // no caching
   SX_CHECK (cacheOverlap == Unknown || cacheOverlap == Recompute);

   if (cacheOverlap == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   return calculateOverlap(ik);
}

SxDiracMat<TAOBasisType> SxAOBasis::getInverseOverlap (int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   // --- single k caching
   if (cacheInverse == CacheCurrentK)  {
      if (ik != invOverlapCachedK)  {
         invOverlapCachedK = ik;
         invOverlap(0) = SxDiracMat<TAOBasisType> ();
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(0) = getOverlap(ik).inverse ();
      }
      return invOverlap(0);
   }
   // --- all k caching
   if (cacheInverse == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < invOverlap.getSize (),
                ik, overlap.getSize ());
      if (invOverlap(ik).getSize () == 0)  {
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(ik) = getOverlap(ik).inverse ();
      }
      return invOverlap(ik);
   }
   // no caching
   SX_CHECK (cacheInverse == Unknown || cacheInverse == Recompute);

   if (cacheInverse == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   SxDiracMat<TAOBasisType> res = getOverlap(ik);
   SX_CLOCK (Timer::AoSInversion);
   res = res.inverse ();

   return res;
}

void SxAOBasis::setOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      overlap.resize (0);
      overlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      overlap.resize (1);
      overlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheOverlap == CacheCurrentK && overlapCachedK != -1)  {
         SxDiracMat<TAOBasisType> ovlp = overlap(0);
         overlap.resize (nk);
         overlap(overlapCachedK) = ovlp;
      }
      overlapCachedK = -1;
   }
   cacheOverlap = mode;
}

void SxAOBasis::setInvOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      invOverlap.resize (0);
      invOverlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      invOverlap.resize (1);
      invOverlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheInverse == CacheCurrentK && invOverlapCachedK != -1)  {
         SxDiracMat<TAOBasisType> ovlp = overlap(0);
         overlap.resize (nk);
         overlap(invOverlapCachedK) = ovlp;
      } else {
         overlap.resize (nk);
      }
      invOverlapCachedK = -1;
   }
   cacheInverse = mode;
}

class SxAOBasisProjGrad
: public SxProjMatrix<TPrecCoeffG>::SaveProjections
{
   public:
      const SxArray<SxAOBasis::OrbitalIndex> &orbitalMap;
      const SxArray<SxDiracMat<TPrecCoeffG> > &refOrbitals;
      mutable PsiG T; // cached phase factors ("atom translation")

      const SxGBasis &gk;

      SxAOBasisProjGrad (const SxAOBasis &ao, int ik)
         : SxProjMatrix<TPrecCoeffG>::SaveProjections (
              3 * ao.getNOrb (),
              ao.blockSize,
              (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).nRows ()),
         orbitalMap(ao.orbitalMap),
         refOrbitals(ao.refOrbitals(ik)),
         gk(refOrbitals(ao.orbitalMap(0).is).getBasis<SxGBasis> ())
      {
         // empty
      }
      virtual ~SxAOBasisProjGrad () {}

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
      virtual void getProjector(int i, SxDiracVec<TGBasisType> *target) const
      {
         int iOrb = i / 3;
         int iDir = i - 3 * iOrb;
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxAOBasis::OrbitalIndex &idx = orbitalMap(iOrb);
         const SxDiracMat<TGBasisType> &refOrb
            = refOrbitals(idx.is).colRef (idx.io);

         if (!T.handle
             || T.handle->auxData.is != idx.is
             || T.handle->auxData.ia != idx.ia)
         {
            PsiG phi = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
            T = SxDiracMat<TPrecCoeffG> (nElements, 3);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (int jDir = 0; jDir < 3; ++jDir)  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
               for (int ig = 0; ig < nElements; ++ig)  {
                  T(ig + jDir * nElements) = I * phi(ig) 
                                           * gk.gVec(ig + jDir * nElements);
               }
            }
            T.handle->auxData.is = idx.is;
            T.handle->auxData.ia = idx.ia;
         }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) =  T(ig + iDir * nElements) * refOrb(ig);
      }
};

SxArray<SxDiracMat<Complex16> >
SxAOBasis::gradProject (const PsiG &psi) const
{
   //SX_CLOCK (Timer::GradProject);
   SX_CHECK(psi.getSize () >0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psi.handle->auxData.ik;
   int nk = int(refOrbitals.getSize ());
   PsiG psiS = SPtr->apply (psi);
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiS.getBasisPtr () == refOrbitals(ik)(0).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(0).getBasisPtr () == psiS.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   return gradProject(psiS, ik);
}

SxArray<SxDiracMat<Complex16> >
SxAOBasis::gradProject (const PsiG &psiS, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradProject);
   int nPsi = (int)psiS.nCols ();
   if (nPsi==0) nPsi = 1;
   int nOrb = getNOrb ();

   SxArray<SxDiracMat<Complex16> > result(3);
   SxDiracVec<Complex16> allProj;
   allProj = SxAOBasisProjGrad(*this, ik).getProjectionFromExtended (psiS);
   allProj.reshape (3, nOrb * nPsi);
   for (int iDir = 0; iDir < 3; ++iDir)  {
      result(iDir) = allProj.row (iDir); // performs a copy!
      result(iDir).reshape (nOrb, nPsi);
      if (extraPhase.getSize () > 0)  {
         const SxDiracVec<Complex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         for (ssize_t iState = 0, iRes = 0; iState < nPsi; ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
               result(iDir)(iRes) *= phase(iOrb).conj ();
      }
   }
   return result;
}



