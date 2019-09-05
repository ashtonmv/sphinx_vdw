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


#include <SxMuPW.h>
#include <SxProjector.h>
#include <SxLoopMPI.h>  // LoopMPI

SxMuPW::SxMuPW (const SxArray<SxArray<SxRadBasis::TPsi> > &in,
                SxConstPtr<SxRadBasis> radBasisPtrIn,
                SxPtr<SxGkBasis>  gkBasisPtrIn,
                bool  initOrbitals)
   : SxPWSet (gkBasisPtrIn),
     mu(in, radBasisPtrIn),
     radBasisPtr(radBasisPtrIn.getConstPtr()),
     isFinalized (false)
{
   SX_CHECK (gkBasisPtr);
   SX_CHECK (radBasisPtrIn);

   int nk = gkBasisPtr->nk;
   muPhi.resize (nk);
   muPhiList.resize (nk);

   if (initOrbitals)  {
      for (int iSpecies=0; iSpecies < radBasisPtr->radFunc.getSize(); iSpecies++)  {
         for (int iot=0; iot < mu.getNOrbTypes(iSpecies); iot++)
            addOrbital (iSpecies, iot);
      }
      finalize ();
   }
   registerMemoryObservers ();
}

SxMuPW::SxMuPW (const SxAtomicOrbitals &in,
              SxPtr<SxGkBasis> gkBasisPtrIn,
              bool initOrbitals)
    : SxPWSet (gkBasisPtrIn),
      isFinalized (false)
{

   mu = in;
   radBasisPtr = in.getRadBasisPtr ().getConstPtr ();

   int nk = gkBasisPtr->nk;
   muPhi.resize (nk);
   muPhiList.resize (nk);

   if (initOrbitals)  {
      for (int iSpecies=0; iSpecies < mu.getNSpecies (); iSpecies++)  {
         for (int iot=0; iot < mu.getNOrbTypes(iSpecies); iot++)
            addOrbital (iSpecies, iot);
      }
      finalize ();
   }
   registerMemoryObservers ();
}

SxMuPW::SxMuPW (const SxArray<SxArray<SxRadBasis::TPsi> > &in,
                SxConstPtr<SxRadBasis> radBasisPtrIn,
                const SxArray<SxArray<bool> > &occMu,
                SxPtr<SxGkBasis>  gkBasisPtrIn)
   : SxPWSet (gkBasisPtrIn),
     mu(in, radBasisPtrIn),
     radBasisPtr(radBasisPtrIn.getConstPtr ()),
     isFinalized (false)
{
   SX_CHECK (gkBasisPtr);
   SX_CHECK (radBasisPtr);

   int nk = gkBasisPtr->nk;
   muPhi.resize (nk);
   muPhiList.resize (nk);

   muPhiList.resize(nk);

   for (int iSpecies=0; iSpecies < radBasisPtr->radFunc.getSize(); iSpecies++)  
      for (int iot=0; iot < mu.getNOrbTypes (iSpecies); iot++)  
         if (occMu(iSpecies)(iot))  addOrbital (iSpecies, iot);
   
   finalize ();
   registerMemoryObservers ();
}

SxMuPW::SxMuPW (const SxAtomicOrbitals &in,
                const SxArray<SxArray<bool> > &occMu,
                SxPtr<SxGkBasis> gkBasisPtrIn)
   : SxPWSet (gkBasisPtrIn),
     isFinalized (false)
{
   SX_CHECK (gkBasisPtr);

   mu = in;
   radBasisPtr = in.getRadBasisPtr ().getConstPtr ();

   int nk = gkBasisPtr->nk;
   muPhi.resize (nk);
   muPhiList.resize (nk);

   muPhiList.resize(nk);

   for (int iSpecies=0; iSpecies < mu.getNSpecies (); iSpecies++)  
      for (int iot=0; iot < mu.getNOrbTypes (iSpecies); iot++)  
         if (occMu(iSpecies)(iot))  addOrbital (iSpecies, iot);
   
   finalize ();
   registerMemoryObservers ();
}


SxMuPW::~SxMuPW ()
{
   // empty
}


void SxMuPW::addOrbital (int iSpecies, int orbitalType)
{
   SX_MPI_LEVEL("waves-k");
   SX_CHECK (!isFinalized);
   const SxGkBasis &Gk = *gkBasisPtr;
   
   SxQuantumNumbers qN = mu.getQuantumNumbers(iSpecies, orbitalType);
   int n = qN.n;
   int l = qN.l;

   int gkIndex = 0;
#ifdef USE_LOOPMPI
   // --- find some k-point that is ours
   for (gkIndex = 0; gkIndex < Gk.nk; ++gkIndex)
      if (SxLoopMPI::myWork (gkIndex)) break;
   if (gkIndex == Gk.nk) return; // nothing to do
#endif

   {
      // --- update both orbitals and refOrbitals
      SX_CHECK (Gk(gkIndex).structPtr);

      int nAtoms = Gk(gkIndex).structPtr->getNAtoms(iSpecies);

      int refOrbStart = int(muPhiList(gkIndex).getSize ());
      for (int iAtom=0; iAtom < nAtoms; iAtom++)  {
         int refOrb = refOrbStart;
         for (int m=-l; m <= l; m++, refOrb++)  {
            psiOrbList.append (SxQuantumNumbers (iSpecies, iAtom, n, l, m));
            phiIdxList.append (refOrb);
         }
      }
   }

   // --- prepare projection via radial interpolation
   // ---get max. |G| 
   double gMax = 0.;
   for (int ik = 0; ik < Gk.nk; ++ik)  {
      if (SxLoopMPI::myWork(ik))  {
         gMax = max(gMax, sqrt(Gk(ik).g2(Gk(ik).ng - 1)));
      }
   }
   gMax *= 1.1;
   // --- estimate number of |G+k| points to compute
   double cutVol = FOUR_PI / 3. * gMax * gMax * gMax;
   double dg = 0.003;
   int ng1 = int(cutVol / Gk(gkIndex).structPtr->cell.getReciprocalCell ().volume);
   int ng2 = int(gMax/dg) + 1;
   // --- setup radial G basis if necessary
   bool viaRadG =  (ng1 * Gk.nk > ng2);
   SxRadGBasis radG;
   SxDiracVec<Double> orbital;
   if (viaRadG) {
      radG.set (0., gMax, ng2, SxRadGBasis::Linear);
      orbital = radG.toSpline (radG | mu(iSpecies,-1,n,l,0));
   }
   
   // --- <G+k|R Y> projection to plane wave basis
   for (int ik=0; ik < Gk.nk; ik++)  {
      if (SxLoopMPI::myWork(ik)) {
         for (int m=-l; m <= l; m++)  {
            if (viaRadG)  {
               orbital.handle->auxData.m = m;
               muPhiList(ik) << ( Gk(ik) | orbital );
            } else {
               muPhiList(ik) << ( Gk(ik) | mu(iSpecies,-1,n,l,m) );
            }
            SxDiracAux &auxData = muPhiList(ik).last ().handle->auxData;
            auxData.is = iSpecies;
            auxData.n = n;
            auxData.l = l;
            auxData.m = m;
         }
      }
   }
   UPDATE_MEMORY (muPhiList);
}

const SxGBasis::TPsi SxMuPW::operator() (int i, int iSpin, int ik) const
{
   SX_CHECK (i >= 0 && i < psiOrb.getSize(), i, psiOrb.getSize());
   SX_CHECK (isFinalized);
   const SxGkBasis &Gk = *gkBasisPtr;

   const SxQuantumNumbers &orb = psiOrb(i);

   int refIdx = phiIdx(i);

   PsiG psiG;
   // --- move projected orbital to the atom's position
   psiG = muPhi(ik)(refIdx) * Gk(ik).getPhaseFactors(orb.iSpecies, orb.iAtom);

   // TODO: ugly !!!
   psiG.handle->auxData.i     = i;
   psiG.handle->auxData.iSpin = iSpin;
   psiG.handle->auxData.ik    = ik;

   // set Basis
   psiG.setBasis(&Gk(ik));

   return psiG;

}

void SxMuPW::finalize ()
{
   SX_CHECK (!isFinalized);
   SX_CHECK (gkBasisPtr);

   phiIdx = phiIdxList;
   psiOrb = psiOrbList;
   for (int ik=0; ik < gkBasisPtr->nk; ik++)
   {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))
      {
         muPhi(ik) = muPhiList(ik);
      }  // LoopMPI
   }

# ifndef NDEBUG
   // Check for unique Quantumnumbers
   for (int i = 0; i < muPhi(0).getSize () - 1; i++)   {
      SxDiracAux &auxDataI = muPhi(0)(i).handle->auxData;
      for (int j = i+1; j < muPhi(0).getSize (); j++)   {
         SxDiracAux &auxDataJ = muPhi(0)(j).handle->auxData;
         if ( (auxDataI.is == auxDataJ.is) && (auxDataI.l == auxDataJ.l) 
               && (auxDataI.m == auxDataJ.m))   {
            SX_CHECK (auxDataI.n != auxDataJ.n, auxDataI.n);
         }
      }
   }
# endif

   phiIdxList.removeAll ();
   psiOrbList.removeAll ();
   muPhiList.resize (0);

   isFinalized = true;
}


int SxMuPW::getNStates (int) const
{
   return int(psiOrb.getSize());
}

int SxMuPW::getNSpin () const
{
   return 1;
}

int SxMuPW::getNk () const
{
   return int(muPhi.getSize());
}

void SxMuPW::setGkBasisPtr (SxPtr<SxGkBasis> gkBasisPtrIn)
{
   gkBasisPtr = gkBasisPtrIn;
   int nk = getNk ();
   // Gk basis must fit unless this is uninitialized
   SX_CHECK (nk <= 0 || gkBasisPtr->getNk () == nk, nk, gkBasisPtr->getNk ());

   // set basis in muPhi
   int ik, iOrb;
   for (ik = 0; ik < nk; ik++)
   {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))
      {
         for (iOrb = 0; iOrb < getNStates (); iOrb++)
         {
            muPhi(ik)(iOrb).setBasis( &(*gkBasisPtr)(ik) );
         }
      }  // LoopMPI
   }
   
}

void SxMuPW::registerMemoryObservers ()
{
   TRACK_MEMORY (psiAtom);
   TRACK_MEMORY (muPhi);
   TRACK_MEMORY (muPhiList);
}
