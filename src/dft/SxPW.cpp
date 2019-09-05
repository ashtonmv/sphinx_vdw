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

#include <SxPW.h>
#include <SxSymMatrix.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <SxMathLib.h>
//#include <SFHIngX.h>
#include <SxLoopMPI.h>  // LoopMPI
#include <SxAllocCache.h>

#ifndef WIN32
#   include <unistd.h>
#endif /* WIN32 */

//Ref1: Rev. Mod. Phys. (64), 1045-1097 (1992)

SxPW::SxPW () : SxPWSet ()
{
   nStates = nSpin = nkPoints = 0;
   loadedK = -1;
   loadedSpin = -1;
   stored = Unknown;
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (int nStatesIn, int nSpinIn, SxPtr<SxGkBasis> gkBasisPtrIn,
            const SxString &tmpDirIn)
   : SxPWSet ()
{
   SX_CHECK (gkBasisPtrIn.getPtr () != NULL);
   SX_CHECK (gkBasisPtrIn->gBasisList(0));
   nStates    = nStatesIn;
   nSpin      = nSpinIn;
   tmpDir     = tmpDirIn;
   gkBasisPtr = gkBasisPtrIn;
   int ng, ik, nk, iSpin;
   nk = gkBasisPtr->nk;
   int nComp = gkBasisPtr->gBasisList(0)->getNComp();
   SxGkBasis &gkBasis = *gkBasisPtr;

   SX_CHECK  (nk > 0, nk);

   nkPoints   = nk;
   if ( tmpDir.getSize () == 0 )  {
      // --- keep waves in memory
      stored = InMemory;
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)
      {
         ng = gkBasis(ik).ng;
         waves(ik).resize (nSpin);
         nStatesPerK.append (nStates);
         nGIPerK.append (nStates*ng);
         for (iSpin=0; iSpin < nSpin; iSpin++)
         {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik))
            {  // LoopMPI
               waves(ik)(iSpin).resize (nStates * ng * nComp);
               waves(ik)(iSpin).reshape (ng * nComp, nStates);
               waves(ik)(iSpin).setBasis (&gkBasis(ik));
            }  // LoopMPI
         }
      }
   }  else  {
      // --- keep waves on disk
      stored = KeepOnDisk;
      createScratchFile ();
      loadWaves (0,0);
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}

SxPW::SxPW (const SxPW &in)
   : SxPWSet ()
{
   this->gkBasisPtr = in.gkBasisPtr;
   SxGkBasis &gkBasis = *gkBasisPtr;
   int ik, nk, iSpin, nGI, nRows, nCols;
   nk = in.getNk();
   nSpin = in.getNSpin ();
   stored = in.stored;
   tmpDir     = in.tmpDir;
   if (stored != InMemory)  {
      SX_EXIT;
      stored = KeepOnDisk;
      createScratchFile ();
   }

   nkPoints   = nk;
   nStates = in.nStates;
   if ( stored == InMemory )  {
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)  {
         waves(ik).resize (nSpin);
         nStatesPerK.append (in.getNStates (ik));
         if (nStates >= 0)  {
            if (nStates != nStatesPerK.last ())
               nStates = -1;
         }
         nGIPerK.append ((int)in.waves(ik)(0).getSize());
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            nRows = (int)in.waves(ik)(iSpin).nRows();
            nCols = (int)in.waves(ik)(iSpin).nCols();
            nGI   = (int)in.waves(ik)(iSpin).getSize();  // {ng x i}
            waves(ik)(iSpin).resize (nGI);
            waves(ik)(iSpin).reshape (nRows, nCols);
            waves(ik)(iSpin) <<= in.waves(ik)(iSpin);
            waves(ik)(iSpin).setBasis( &gkBasis(ik) );
         }
      }
   }  else  {  // keep waves on disk
      SX_EXIT;
   }

   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (int nSpinIn, const SxList<int> &ngPerK,
            const SxString &tmpDirIn)
   : SxPWSet ()
{
   SX_EXIT;
   gkBasisPtr = SxPtr<SxGkBasis>::create();  // ???
   tmpDir = tmpDirIn;
   int ng, ik, nk, iSpin;
   nk = int(ngPerK.getSize());
   bool keepOnDisc = (tmpDir != "");
   nkPoints   = nk;

   if (keepOnDisc)  {
      SX_EXIT;
   }

   SX_CHECK  (nk > 0, nk);
   SX_CHECK  (nSpin == 1 || nSpin == 2, nSpin);
   nSpin = nSpinIn;


   if (!keepOnDisc)  {
      stored = InMemory;
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)  {
         ng = ngPerK(ik);
         waves(ik).resize (nSpin);
         nStatesPerK.append (nStates);
         nGIPerK.append (nStates*ng);
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            waves(ik)(iSpin).resize (nStates * ng);
            waves(ik)(iSpin).reshape (ng, nStates);
         }
      }
   }  else  {
      // --- keep waves on disk
      stored = KeepOnDisk;
      createScratchFile ();
      loadWaves (0,0);
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (const SxString filename, enum StorageModel how, const SxString &tmpDirIn)
   : SxPWSet (),
     nkPoints(-1),
     nStates(-1),
     nSpin(-1),
     stored(how)
{
   SX_CHECK (how == InMemory || how == ReadOnDemand ||
             how == ReadOneByOne || how == KeepOnDisk);
   try {
      wavesFile.open (filename, SxBinIO::BINARY_READ_ONLY);
      if (how == InMemory)  {
         read (wavesFile);
         wavesFile.close ();
      } else if (how == KeepOnDisk) {
         if (tmpDirIn.isEmpty ()) tmpDir =".";
         else tmpDir = tmpDirIn;
         stored = KeepOnDisk;
         createScratchFile ();
         read (wavesFile);
         wavesFile.close ();
      } else {
         nkPoints = wavesFile.getDimension ("nk");
         nSpin    = wavesFile.getDimension ("nSpin");
         nStatesPerK.resize (nkPoints);
         wavesFile.read ("nPerK", &nStatesPerK, nkPoints);

         // get nGIPerK
         SxDiracVec<Int> nGk(nkPoints);
         wavesFile.read ("nGk", &nGk, nkPoints);
         nGIPerK.resize (nkPoints);
         for (int ik = 0; ik < nkPoints; ++ik)
            nGIPerK(ik) = nGk(ik) * nStatesPerK(ik);

         loadedK = -1;
         loadedSpin = -1;
         waves.resize(1);
         waves(0).resize (1);
         nStates = nStatesPerK(0); // ad hoc
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}

SxPW::~SxPW ()
{
   if (stored == KeepOnDisk)  {
      SxString tmpFilename = wavesFile.filename;
      try {
         wavesFile.close ();
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
      unlink (tmpFilename.ascii());
   }
   if (stored == ReadOnDemand || stored == ReadOneByOne)  {
      try {
         wavesFile.close ();
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
   }
}

void SxPW::operator= (const SxPW &in)
{
   SxPWSet::operator= (in);
   //cout << "call to SxPW::operator=" << endl; // identify hidden calls
   SX_CHECK (stored != KeepOnDisk); // would need closing of current file
   if (stored == ReadOnDemand || stored == ReadOneByOne)
      wavesFile.close ();
   this->gkBasisPtr = in.gkBasisPtr;
   nStates     = in.nStates;
   nSpin       = in.nSpin;
   nkPoints    = in.getNk();
   nStatesPerK = in.nStatesPerK;
   nGIPerK     = in.nGIPerK;

   stored     = in.stored;
   if (stored == KeepOnDisk)  {
      // take over file
      tmpDir      = in.tmpDir;
      loadedK     = in.loadedK;
      loadedSpin = in.loadedSpin;
      wavesFile.fp          = in.wavesFile.fp;
      in.wavesFile.fp       = NULL;
      wavesFile.mode        = SxBinIO::SCRATCH_READ_WRITE;
      in.wavesFile.mode     = SxBinIO::UNKNOWN;
      wavesFile.filename    = in.wavesFile.filename;
      in.wavesFile.filename = "";
      wavesFile.isOpen      = true;
      in.wavesFile.isOpen   = false;
      in.stored             = Unknown;
   }
   if (stored == ReadOnDemand)  {
      loadedK = in.loadedK;
      loadedSpin = in.loadedSpin;
      try  {
         wavesFile.open (in.wavesFile.filename, SxBinIO::BINARY_READ_ONLY);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }
   waves = in.waves;
   UPDATE_MEMORY (waves);
}

SxPtr<SxPWSet> SxPW::getNew () const
{
   return SxPtr<SxPW>::create (getNStates (), getNSpin (), gkBasisPtr,
                               tmpDir);
}

SxGkBasis::TPsi& SxPW::operator() (int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand);
   // no ReadOneByOne (since one by one implies only one state at a time)
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }
      return waves(0)(0);
   }  else  {
      waves(ik)(iSpin).handle->auxData.i     = -1;
      waves(ik)(iSpin).handle->auxData.iSpin = iSpin;
      waves(ik)(iSpin).handle->auxData.ik    = ik;
      /*
      SX_CHECK(waves(ik)(iSpin).handle->auxData.i     == -1);
      SX_CHECK(waves(ik)(iSpin).handle->auxData.iSpin == iSpin);
      SX_CHECK(waves(ik)(iSpin).handle->auxData.ik    == ik);
      */
      return waves(ik)(iSpin);
   }
}


const SxGkBasis::TPsi& SxPW::operator() (int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand);
   SxGkBasis::TPsi psiOut;
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         if (loadedK >= 0)  {
            SX_CHECK  (iSpin >= 0, iSpin);
            flushWaves ();
         }
         loadWaves (ik,iSpin);
      }
      return waves(0)(0);
   }  else  {
      waves(ik)(iSpin).handle->auxData.i     = -1;
      waves(ik)(iSpin).handle->auxData.iSpin = iSpin;
      waves(ik)(iSpin).handle->auxData.ik    = ik;
      /*
      SX_CHECK(waves(ik)(iSpin).handle->auxData.i     == -1);
      SX_CHECK(waves(ik)(iSpin).handle->auxData.iSpin == iSpin);
      SX_CHECK(waves(ik)(iSpin).handle->auxData.ik    == ik);
      */

      return waves(ik)(iSpin);
   }
}

SxGkBasis::TPsi SxPW::operator() (const SxIdx &idx, int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (idx.start >= 0, idx.start);
   SX_CHECK (idx.end <= getNStates(ik), idx.end);
   SxGkBasis::TPsi &allStates = operator() (iSpin, ik), res;
   int ng = (int)allStates.nRows ();
   res = allStates(SxIdx (idx.start * ng, idx.end * ng + ng - 1));
   res.reshape (ng, idx.end - idx.start + 1);
   res.handle->auxData = allStates.handle->auxData;
   return res;
}

const SxGkBasis::TPsi
SxPW::operator() (const SxIdx &idx, int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (idx.start >= 0, idx.start);
   SX_CHECK (idx.end <= getNStates(ik), idx.end);
   const SxGkBasis::TPsi &allStates = operator() (iSpin, ik);
   SxGkBasis::TPsi res;
   int ng = (int)allStates.nRows ();
   res = allStates(SxIdx (idx.start * ng, idx.end * ng + ng - 1));
   res.reshape (ng, idx.end - idx.start + 1);
   res.handle->auxData = allStates.handle->auxData;
   return res;
}

SxGkBasis::TPsi SxPW::operator() (int i, int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand || stored == ReadOneByOne);
   SxGkBasis::TPsi psiOut;
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }

      psiOut = waves(0)(0).colRef(i);
   } else if (stored == ReadOneByOne) {
      if (   waves(0)(0).getSize () == 0
          || waves(0)(0).handle->auxData.i != i
          || waves(0)(0).handle->auxData.ik != ik)
         loadWaves(i, iSpin, ik);
      psiOut = waves(0)(0);
   }  else  {
      psiOut = waves(ik)(iSpin).colRef(i);
   }
   psiOut.handle->auxData.i     = i;
   psiOut.handle->auxData.iSpin = iSpin;
   psiOut.handle->auxData.ik    = ik;
   return psiOut;
}


const SxGkBasis::TPsi SxPW::operator() (int i, int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk
             || stored == ReadOnDemand || stored == ReadOneByOne);
   SxGkBasis::TPsi psiOut;
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }
      psiOut = waves(0)(0).colRef(i);
   } else if (stored == ReadOneByOne) {
      if (   waves(0)(0).getSize () == 0
          || waves(0)(0).handle->auxData.i != i
          || waves(0)(0).handle->auxData.ik != ik)
         loadWaves(i, iSpin, ik);
      psiOut = waves(0)(0);
   }  else if (stored == InMemory)  {
      psiOut = waves(ik)(iSpin).colRef(i);
   }
   psiOut.handle->auxData.i     = i;
   psiOut.handle->auxData.iSpin = iSpin;
   psiOut.handle->auxData.ik    = ik;
   return psiOut;
}


int SxPW::getNStates (int ik) const
{
   if (ik == -1) return nStates;
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   if (stored == InMemory && waves(ik).getSize () > 0
       && waves(ik)(0).getSize () > 0)
   {
      SX_CHECK (nStatesPerK(ik) == waves(ik)(0).nCols (),
                nStatesPerK(ik), waves(ik)(0).nCols ());
      return (int)waves(ik)(0).nCols ();
   }
   return nStatesPerK(ik);
}


void SxPW::changeNStates (int newN)
{
   SX_CHECK (newN > getNStates() || stored == ReadOnDemand ||
             stored == ReadOneByOne, newN, getNStates());

   int i, nOrig, iSpin, ng, ik, nk = getNk ();
   nkPoints = nk;
   nSpin = getNSpin();

   if ( stored == KeepOnDisk )  {
      // --- resize temporary file
      SX_EXIT;
      SxString tmpFilename = wavesFile.filename;
      wavesFile.close ();
      unlink (tmpFilename.ascii());  // remove old file
      createScratchFile ();
   }
   if (stored == ReadOnDemand || stored == ReadOneByOne)  {
      for (ik = 0; ik < nk; ++ik)  {
         if (nStatesPerK(ik) < newN)  {
            cout << "Too few states in '" << wavesFile.filename;
            cout << "' for SxPW::changeNStates to " << newN << endl;
            SX_EXIT;
         }
         nGIPerK(ik) /= nStatesPerK(ik);
         nStatesPerK(ik) = newN;
         nGIPerK(ik) *= newN;
      }
      loadedK = -1;
      loadedSpin = -1;
      return;
   }

   SxPW origWaves ( *this );

   for (ik=0; ik < nk; ik++)  {
      //waves.resize (nk);
      nOrig = origWaves.getNStates (ik);
      //waves(ik).resize (nSpin);
      nStatesPerK(ik) = newN;
      ng          = (int)origWaves.waves(ik)(0).nRows();
      nGIPerK(ik) = (int)origWaves.waves(ik)(0).nRows() * newN;
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         waves(ik)(iSpin).resize  (ng * newN);
         waves(ik)(iSpin).reshape (ng, newN);
         for (i=0; i < nOrig; i++)
            waves(ik)(iSpin).colRef(i) <<= origWaves(i,iSpin,ik);
         for (i=nOrig; i < newN; i++)
            waves(ik)(iSpin).colRef(i).randomize();
      }
   }
   orthonormalize ();
}



void SxPW::randomize ()
{
   int ik, nk = getNk(), iSpin;
   nkPoints = nk;
   SX_MPI_LEVEL("waves-k");
   for (ik=0; ik < nk; ik++)  {
      if (SxLoopMPI::myWork(ik))  {
         for (iSpin=0; iSpin < nSpin; iSpin++)
            (*this)(iSpin,ik).randomize();
      }
   }

//   orthonormalize ();
}



SxPW &SxPW::normalize ()
{
   int iSpin, ik, nk = int(waves.getSize());
   nkPoints = nk;
   for (ik=0; ik < nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         waves(ik)(iSpin).normalize ();
   return *this;
}


SxPW &SxPW::orthogonalize (enum SxOrthoMethod)
{
   SX_EXIT;
   return *this;
}


SxPW &SxPW::orthogonalize (SxPW &)
{
   SX_EXIT;
   return *this;
}

void SxPW::setOrthogonal (SxDiracVec<TPrecCoeffG> *psiPtr,
                          int firstN, int iSpin, int ik,
                          enum NormMethod normal,
                          double threshold)
{
   // by C. Freysoldt
   // This routine is heavily commented on technical details
   // This is because I'm implementing it now, but the XPress library
   // will need changes here, but maybe not by me
   SX_CLOCK (Timer::ortho);

   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (firstN >= 0 && firstN <= getNStates (ik),
             firstN, getNStates (ik));
   SX_CHECK (threshold >= 0., threshold);
   SX_CHECK (   normal == DONT_NORMALIZE
             || normal == NORMALIZE
             || normal == RENORMALIZE);

   SX_CHECK (psiPtr);
   if (firstN == 0)  {
      if (normal == NORMALIZE) psiPtr->normalize ();
      return;
   }
   SxDiracVec<TPrecCoeffG> &psi = *psiPtr;
   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;

   // psi must be normalized
   SX_CHECK (! (normal == RENORMALIZE)
             || fabs (1. - dot(psi,psi).re) < 1e-7,
             1. - dot(psi,psi).re        );

   int j;
   SxDiracVec<TPrecCoeffG> psiJ;
   PrecCoeffG scp;
   if (firstN > 150 || nPsi > 1)  {
      // for low n, the direct (convential) BLAS1 approach is faster
      // because there is a certain overhead, and BLAS2 routines may
      // not be fastet for BLAS1-like tasks
      // the cross-over point should be tested
      // this should be the zdotc / gemv crossover for nPsi == 1
      // this should be the gemv / gemm crossover for nPsi > 1
      SX_CHECK (stored == InMemory || loadedK == ik, loadedK, ik);
      SX_CHECK  (loadedK >= 0, loadedK);
      SX_CHECK  (loadedSpin >= 0, loadedSpin);

      SxDiracVec<TPrecCoeffG> *allStates
         = (stored == InMemory) ? &(waves(ik)(iSpin))
                                : &(waves(0)(0));

      // --- get matrix |j> for j < firstN
      int ng = (int)allStates->nRows ();
      SxDiracVec<TPrecCoeffG> lowerStates
         = (*allStates)(SxIdx(0, ng * firstN - 1));
      lowerStates.reshape(ng, firstN);
      // --- get <j|psi> for all 0 <= j < firstN

      SxDiracVec<TPrecCoeffG> scps;
      if (nPsi > firstN)
         // --- variant 1 (proper for XPress)
         scps = lowerStates.adjoint () ^ psi;
      else
         // --- variant 2
         scps = (psi.adjoint () ^ lowerStates).adjoint ();

      // --- subtract |j><j|psi>
      // bool loopSubtr = (firstN < 500); // axpy / gemv crossover, to be tested
      bool loopSubtr = false;
      /*
      if (! loopSubtr && nPsi == 1)  {
         // count negligible scps
         int scpNegligible = firstN / 2; // if fewer negligible scps, use BLAS2
         SxDiracVec<TPrecCoeffG>::Iterator it = scps.begin ();
         for (j = 0; (j < firstN) && (scpNegligible > 0); j++, ++it)
            if ((*it).absSqr () > threshold) scpNegligible--;
         loopSubtr = (scpNegligible > 0);
      }
      */
      if (nPsi == 1 && loopSubtr)  {
         // BLAS1 version
         SxDiracVec<TPrecCoeffG>::Iterator it = scps.begin ();
         for (j = 0; j < firstN; j++, ++it)
            if ((*it).absSqr () > threshold)
               psi.plus_assign_ax (-(*it), lowerStates.colRef(j));
      } else {
         // BLAS2 (nPsi == 1) or BLAS3 version
         psi -= (lowerStates ^ scps);
      }

      // --- renormalize
      if (normal == RENORMALIZE)  {
         SX_CHECK (nPsi == 1, nPsi);
         double deltaNorm = dot (scps, scps).re;
         if (deltaNorm > threshold)  {
            psi /= sqrt(1. - deltaNorm);
         }
      }

   } else {
      // --- conventional BLAS1 approach
      double norm = 1.;
      double scp2;
      for (j=0; j < firstN; j++)  {
         psiJ = operator() (j,iSpin,ik);
         scp = dot (psiJ, psi); // <j|psi>
         scp2 = scp.absSqr ();
         if (normal == RENORMALIZE) norm -= scp2;
//       psi -= psiJ * scp;
         if (scp2 > threshold)
            psi.plus_assign_ax (-scp, psiJ); // psi -= |j><j|psi>
      }

      // --- renormalization
      if (normal == RENORMALIZE)  {
         if (fabs(1. - norm) > threshold)  {
            psi /= sqrt(norm);
         }
      }
   }

   if (normal == NORMALIZE)  {
      if (nPsi == 1)
         psi.normalize ();
      else
         for (int i = 0; i < nPsi; ++i) psi.colRef(i).normalize ();
   }

   /*
   // --- check normalization
   if (normal == RENORMALIZE)  {
      double x = 1. - dot(psi, psi).re;
      if (fabs(x) > threshold)  {
        cout << x << endl;
        psi.normalize ();
      }
   }
   */
}

SxPW
&SxPW::orthonormalize (enum SxOrthoMethod method,
                       SxArray<SxArray<SxDiracMat<TPrecCoeffG> > > *uPtr)
{
   SxPW &psiSet = *this;
//   timer->start (TIME_ORTHO);

   if (method == GramSchmidt)  {
      //   SX_CLOCK (Timer::ortho); // done in setOrthogonal
      // --- Gram/Schmidt orthogonalization
      SxGkBasis::TPsi psiI, psiJ;
      PrecCoeffG scp;
      int i, nI, ik, iSpin, nk = getNk();
      SX_MPI_LEVEL("waves-k");
      for (ik=0; ik < nk; ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         nI = nStatesPerK(ik);
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            for (i=0; i < nI; i++)  {
               psiI = psiSet(i,iSpin,ik);
               setOrthogonal (&psiI, i, iSpin, ik, NORMALIZE, 0.);
            }
         }
      }
   }  else  {
   // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SX_CLOCK (Timer::ortho);
      PsiGI psiGI;
      int i, j, iSpin, ik, nk = getNk();

      SxDiracMat<TPrecCoeffG> U, L, Ieps;      // S being the overlap matrix,
                                               // U is U^-1/2, and
      SxDiracSymMat<TPrecCoeffG> S;
      SxDiracSymMat<TPrecCoeffG>::Eigensystem eig;

      SX_MPI_LEVEL("waves-k");
      for (ik=0; ik < nk; ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            S = SxDiracSymMat<TPrecCoeffG> (nStatesPerK(ik));
            {
               SX_CLOCK (Timer::loewdinS);
               for (i=0; i < getNStates(ik); i++)  {
                  for (j=i; j < getNStates(ik); j++)  {
                     S(i,j)  = (psiSet(i,iSpin,ik) ^ psiSet(j,iSpin,ik)).chop();
                  }
               }
            }

            eig   = S.eigensystem();
            Ieps  = Ieps.identity(1. / sqrt(eig.vals) );
            U     = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

            // --- store U, if required
            if (uPtr)  (*uPtr)(ik)(iSpin) = U;

//          psiSet(iSpin,ik) <<= (U ^ psiSet(iSpin,ik).adjoint()).adjoint();
            psiSet(iSpin,ik) <<= (psiSet(iSpin,ik) ^ U);

         }
      }
   }
//   timer->stop (TIME_ORTHO);

   return *this;
}


bool SxPW::isOrthogonal ()
{
   SxPW &psiSet = *this;
   SxGkBasis::TPsi psiI, psiJ;
   int i, j, ik, iSpin, nk = getNk();
   PrecCoeffG scp;
   for (ik=0; ik < nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (i=0; i < getNStates(ik); i++)  {
            psiI = psiSet(i,iSpin,ik);
            for (j=0; j < getNStates(ik); j++)  {
               scp = (psiI ^ psiSet(j,iSpin,ik)).chop();
               if ( i == j )  {
                  if ( fabs(scp.re /* - 1.*/ ) < 1e-8 )  {
                     cout << i << "/" << j << ": (a) " << scp << endl;
                     return false;
                  }
               }  else  {
                  if ( fabs(scp.re) > 1e-8 )  {
                     cout << i << "/" << j << ": (b) " << scp << endl;
                     return false;
                  }
               }
               if ( fabs(scp.im)    > 1e-8 )  {
                  cout << i << "/" << j << ": (c) " << scp << endl;
                  return false;
               }
            }
         }
      }
   }
   return true;
}



bool SxPW::isOrthonormal ()
{
   SxPW &psiSet = *this;
   SxGkBasis::TPsi psiI, psiJ;
   int i, j, ik, iSpin, nk = getNk();
   PrecCoeffG scp;
   for (ik=0; ik < nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         for (i=0; i < nStatesPerK(ik); i++)  {
            psiI = psiSet(i,iSpin,ik);
            for (j=0; j < nStatesPerK(ik); j++)  {
               scp = (psiI ^ psiSet(j,iSpin,ik)).chop();
               if ( i == j )  {
                  if ( fabs(scp.re - 1. ) > 1e-8 )  {
                     sxprintf ("(%d,%d,%d): ", i, iSpin, ik);
                     cout << i << "/" << j << ": (a) " << scp << endl;
                     return false;
                  }
                  if ( fabs(scp.im      ) > 1e-8 )  {
                     sxprintf ("(%d,%d,%d): ", i, iSpin, ik);
                     cout << i << "/" << j << ": (a) " << scp << endl;
                     return false;
                  }
               }  else  {
                  if ( fabs(scp.re) > 1e-8 )  {
                     sxprintf ("(%d,%d,%d): ", i, iSpin, ik);
                     cout << i << "/" << j << ": (b) " << scp << endl;
                     return false;
                  }
               }
               if ( fabs(scp.im)    > 1e-8 )  {
                  sxprintf ("(%d,%d,%d): ", i, iSpin, ik);
                  cout << i << "/" << j << ": (c) " << scp << endl;
                  return false;
               }
            }
         }
      }
   }
   return true;
}

void SxPW::setGkBasisPtr (SxPtr<SxGkBasis> gkBasisPtrIn)
{
   gkBasisPtr = gkBasisPtrIn;
   SxGkBasis &gkBasis = *gkBasisPtr;
   int nk = getNk ();
   // Gk basis must fit unless this is uninitialized
   SX_CHECK (nk <= 0 || gkBasisPtr->getNk () == nk, nk, gkBasisPtr->getNk ());

   // set basis in waves
   if (stored == InMemory)  {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int ik = 0; ik < nk; ik++)  {
           waves(ik)(iSpin).setBasis( &gkBasis(ik) );
         }
      }
   }   else if (loadedK >= 0)  {
      waves(0)(0).setBasis( &gkBasis(loadedK) );
   }
}


void SxPW::setNStates (const SxArray<int> &nPerK)
{
   SX_CHECK (waves.getSize() == nPerK.getSize(),
             waves.getSize(),   nPerK.getSize());
   SX_CHECK      (gkBasisPtr);

   int i, n, iSpin, ik, nGk, sizeOld, nOld;
   int nk = getNk ();
   const SxGkBasis &Gk = *gkBasisPtr;

   nStates = nPerK(0);
   if (stored == InMemory)  {

      for (ik=0; ik < nk; ik++)  {
         n = nPerK(ik);
         if (n < nStates) nStates = n; // get min. no. of states
         nGk = Gk(ik).ng;
         if ( n > nGk)  n = nGk;
//         SX_CHECK (n == nGk, n, nGk);
         nGIPerK(ik) = (int)waves(ik)(0).nRows() * nPerK(ik);
         for (iSpin=0;  iSpin < nSpin;  iSpin++)  {
            sizeOld = (int)waves(ik)(iSpin).getSize ();
            nOld = (int)waves(ik)(iSpin).nCols();
            waves(ik)(iSpin).reshape (sizeOld, 1);
            waves(ik)(iSpin).resize  (nGk*n, true);
            waves(ik)(iSpin).reshape (nGk, n);

            // --- fill rest of wavefunctions with random numbers
            if (nOld < n)  {
               for (i = nOld; i < n; i++)  {
                  waves(ik)(iSpin).colRef(i).randomize ();
//                  waves(ik)(iSpin).normalize ();  // maybe, that gets necess.
               }
            }
         }
      }
   }  else  {
      // TODO: setNStates(nPerK) for case 'keepOnDisk' not yet implemented
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);
}


void SxPW::createScratchFile ()
{
   SX_CHECK (stored == KeepOnDisk);
   SX_CHECK (!wavesFile.isOpen);
   SX_CHECK (gkBasisPtr);
   const SxGkBasis &gkBasis = *gkBasisPtr;
   try  {
      wavesFile.open (SxBinIO::createTempName(tmpDir), SxBinIO::SCRATCH_READ_WRITE);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   waves.resize (1);  // only 1 k-point is kept in memory
   waves(0).resize (1); // resize by 1 spin
   size_t nElem, byteLen = sizeof (waves(0)(0)(0));
   int ik, iSpin;
   int nk = nkPoints, ng;
   nStatesPerK.resize (0);
   nGIPerK.resize (0);
   for (ik=0; ik < nk; ik++)  {
      ng = gkBasis(ik).ng;
      nElem = nStates*ng;
      nStatesPerK.append (nStates);
      nGIPerK.append ((int)nElem);
      // --- fill file buffer
      waves(0)(0).reformat(ng, nStates);
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         fwrite (waves(0)(0).elements, byteLen, nElem, wavesFile.fp);
      }
   }
   loadedK    = -1;
   loadedSpin = -1;
}

void SxPW::flushWaves () const
{
   if (stored == ReadOnDemand || stored == ReadOneByOne) return;
   SX_CHECK (stored == KeepOnDisk);
   if (loadedK == -1) return;
   SX_CHECK (loadedK <= nkPoints, loadedK, nkPoints);
   int i;
   size_t pos, byteLen = sizeof (waves(0)(0)(0));
#  ifndef NDEBUG
   cout << "Flushing k-point " << (loadedK+1) << endl;
#  endif /* NDEBUG */
   for (pos=0, i=0; i < loadedK; i++)  pos += nGIPerK (i)*nSpin;
   pos += nGIPerK (loadedK)*loadedSpin;
   fseek (wavesFile.fp, pos*byteLen, SEEK_SET);
   fwrite (waves(0)(0).elements, byteLen, nGIPerK(loadedK),
              wavesFile.fp);
}


void SxPW::loadWaves (int ik, int iSpin) const
{
   SX_CHECK (stored == KeepOnDisk || stored == ReadOnDemand);
   int i, ng;
   size_t pos=0, byteLen = sizeof (PrecCoeffG);
   SxDiracVec<Int> nPerK;
   /* FFT mapping not in use: assume G-basis fits, because it can be read
      from file.  This allows to use ReadOnDemand with modified FFT meshes.

   SxDiracVec<Int> fftIdx;

   */
#  ifndef NDEBUG
   cout << "Loading k-point " << (ik+1) << endl;
#  endif /* NDEBUG */

   bool readKwise = false;

   // --- get starting point
   if (stored == KeepOnDisk)  {
      for (i=0; i < ik; i++) pos += nGIPerK (i)*nSpin;
      pos += nGIPerK (ik)*iSpin;
      fseek (wavesFile.fp, pos*byteLen, SEEK_SET);
      ng = nGIPerK(ik) / nStatesPerK(ik);
   } else {
      nPerK.resize (nkPoints);
      wavesFile.read ("nPerK", &nPerK, nkPoints);
      SxDiracVec<Int> nGk(nkPoints);
      wavesFile.read ("nGk", &nGk, nkPoints);
      int fftOffset = 0;
      if (wavesFile.containsDim("nCoeff"))  {
         readKwise = false;
         for (i=0; i < ik; i++)  {
            pos += nPerK(i) * nGk(i) * nSpin;
            fftOffset += nGk(i);
         }
      } else {
         readKwise = true;
         pos = 0;
      }
      ng = nGk(ik);
      pos += nPerK(ik) * nGk(ik) * iSpin;
      /* FFT mapping not in use, see above
      if (gkBasisPtr)  {
         fftIdx.resize (ng);
         wavesFile.read ("fftIdx", &fftIdx, ng, fftOffset);
      }
      */
   }

   PsiG psi, psiTmp;
   size_t offset = pos; // for ReadOnDemand
   if ((stored == ReadOnDemand) && (!gkBasisPtr))
      psiTmp.resize (ng);

   SxString psiVarName = readKwise ? SxString("psi-") + SxString(ik+1)
                                   : SxString("psi");

   waves(0)(0) = PsiGI ();
   waves(0)(0).reformat(ng, nStatesPerK(ik));
   if (stored == KeepOnDisk)  {
      size_t nRead
         = fread (waves(0)(0).elements, byteLen, nGIPerK(ik), wavesFile.fp);
      if (nRead != byteLen * nGIPerK(ik)) { SX_EXIT; }
   } else  {
        for (i = 0; i < nStatesPerK(ik); ++i)  {
            psi = waves(0)(0).colRef (i);
            /* FFT mapping not in use, see above

            if (gkBasisPtr)  {
               wavesFile.read (psiVarName, &psiTmp, ng, offset);
               psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx);
            } else
            */
            {
               wavesFile.read (psiVarName, &psi, ng, (int)offset);
            }
            offset += ng;
         }
         offset += ng * (nPerK(ik) - nStatesPerK(ik));
      }
   if (gkBasisPtr)
   waves(0)(0).setBasis (& (*gkBasisPtr)(ik));
   waves(0)(0).handle->auxData.i     = -1;
   waves(0)(0).handle->auxData.iSpin = iSpin;
   waves(0)(0).handle->auxData.ik    = ik;
   loadedK = ik;
   loadedSpin = iSpin;

}

void SxPW::loadWaves (int i, int iSpin, int ik) const
{
   SX_CHECK (stored == ReadOneByOne);
   waves(0)(0) = PsiGI ();
   waves(0)(0) = readPsi(wavesFile, i, iSpin, ik);
   if (gkBasisPtr)
      waves(0)(0).setBasis (& (*gkBasisPtr)(ik));
}

PsiG SxPW::readPsi (const SxBinIO &io, int i, int iSpin, int ik)
{
   // SX_CHECK (nSpin == 1,nSpin);

   int ng;
   size_t pos=0;
   /* FFT mapping not in use: assume G-basis fits, because it can be read
      from file.  This allows to use ReadOnDemand with modified FFT meshes.

   SxDiracVec<Int> fftIdx;

   */

   // --- get starting point
   int nk, nSpin;
   try {
      nk    = io.getDimension ("nk");
      nSpin = io.getDimension ("nSpin");
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   SxDiracVec<Int> nPerK(nk);
   SxDiracVec<Int> nGk(nk);
   try {
      io.read ("nPerK", &nPerK, nk);
      io.read ("nGk", &nGk, nk);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   size_t fftOffset = 0;
   SxString psiVarName;
   if (io.containsDim("nCoeff"))  {
      psiVarName = "psi";
      for (i=0; i < ik; i++)  {
         pos += nPerK(i) * nGk(i);
         fftOffset += nGk(i);
      }
   } else {
      psiVarName = "psi-" + SxString(ik+1);
      pos = 0;
   }
   ng = nGk(ik);
   /* FFT mapping not in use, see above
   if (gkBasisPtr)  {
      fftIdx.resize (ng);
      io.read ("fftIdx", &fftIdx, ng, fftOffset);
   }
   */

   size_t offset = nSpin * pos + size_t(i) * ng;

   PsiG psi(ng);

   try {
      /* FFT mapping not in use, see above
      if (gkBasisPtr)  {
         PsiG psiTmp(ng);
         io.read (psiVarName, &psiTmp, ng, offset);
         psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx);
      } else
      */
      io.read (psiVarName, &psi, ng, (int)offset);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   SxDiracAux &auxData = psi.handle->auxData;
   auxData.i     = i;
   auxData.iSpin = iSpin;
   auxData.ik    = ik;
   return psi;
}

void SxPW::memMinimize ()
{
   if (stored == InMemory) return; // nothing to do
   // write to disk (if necessary)
   if (stored == KeepOnDisk) flushWaves ();
   // unload waves
   loadedK = -1;
   loadedSpin = -1;
   waves(0)(0) = SxDiracVec<TPrecCoeffG> ();
}

void SxPW::read (const SxBinIO &io, int mode)
{
   SX_CHECK (stored != ReadOnDemand && stored != ReadOneByOne);
   if (stored == Unknown) stored = InMemory;
   SX_CHECK (stored == InMemory || mode & KeepNStates); // to be implemented

   bool hasBasis = (gkBasisPtr.getPtr() != NULL);
   bool needsReorder = hasBasis;

   cout << "SxPW::read: |G+k> basis ";
   if (!hasBasis) cout << "not ";
   cout << "present." << endl;

   // --- variable for k-wise read (allows larger files)
   bool readKwise;
   SxString psiVarName;

   try  {
      int ik, nk, iSpin, nSpinFile, ng, i, iOffset;
#ifndef NDEBUG
      int nCoeff = -1;
#endif
      int nComp = 1;
      if(io.containsDim("nComp"))
         nComp = io.getDimension("nComp");
      nk        = io.getDimension ("nk");
      nSpinFile = io.getDimension ("nSpin");
      if (nSpinFile != nSpin && nSpin > 0)  {
         // can happen in non-DEBUG cases, so SX_CHECK not sufficient
         cout << "Inconsistency: this->nSpin = " << this->nSpin;
         cout << "; nSpin = " << nSpinFile << " in " << io.filename << endl;
         SX_EXIT;
      } else {
         nSpin = nSpinFile;
      }
      readKwise = !io.containsDim("nCoeff");
      if (readKwise)  {
         cout << "Reading waves k-wise" << endl;
      } else {
#ifndef NDEBUG
         nCoeff   = io.getDimension ("nCoeff");
#endif
         psiVarName = "psi";
      }
      sxprintf ("Reading wavefunction:\n");
      sxprintf ("---------------------\n");
      sxprintf ("nk: %d, nSpin: %d\n", nk, nSpin);

      SxString name;
      int offset = 0;
      int nAllGk = io.getDimension ("nAllGk");
      SxDiracVec<Int> nPerK (nk), nGk (nk), fftIdxIn(nAllGk);
      SxVector<TPrecFFTIdx> fftIdx;

      io.read ("nPerK",    &nPerK,  nk);
      io.read ("nGk",      &nGk,    nk);
      io.read ("fftIdx",   &fftIdxIn, nAllGk);
      // --- check if #k, k and #G(k) are correct
      if (hasBasis && nk != gkBasisPtr->nk)  {
         sxprintf ("Error: Input file is corrupt.\n");
         sxprintf ("       Number of k-points: %d (should be %d)\n", nk,
                 gkBasisPtr->nk);
         SX_EXIT;
      }
      if (hasBasis)  {
         SxMatrix<TPrecG> kVecs(nk,3); // need matrix for reading
         Coord kVec;
         io.read ("kVec", &kVecs, nk, 3);
         for (ik = 0; ik < nk; ik++)  {
            kVec = kVecs.row(ik).toVector3 ();
            if ( (kVec - gkBasisPtr->getK(ik)).normSqr () > 1e-10)  {
               cout << endl;
               cout << "Consistency check failed for k-point " << (ik+1);
               cout << ".\n Expected " << gkBasisPtr->getK(ik) << " but found ";
               cout << kVec<< " in '" << io.filename << "'." << endl;
               SX_EXIT; // may indicate changes in the k-point setup!
               // if this occurs, we should remap the k-points
            }
         }
      }
      for (ik=0; ik < nk; ik++)  {
         if (hasBasis && nGk(ik) != (*gkBasisPtr)(ik).ng)  {
            sxprintf ("Error: Input file is corrupt.\n");
            sxprintf ("       Number of |G+k> vectors: %d (should be %d)\n",
                     nGk(ik), (*gkBasisPtr)(ik).ng);
            SX_EXIT;
         }
      }

      SxMesh3D mesh;
      io.read("meshDim", &mesh);

      // --- resize to kpoints
      if (stored == InMemory)  {
         waves.resize (nk);
      }
      if (!(mode & KeepNStates))  nStatesPerK.resize (nk);
      SX_CHECK (nStatesPerK.getSize () == nk, nStatesPerK.getSize (), nk);

      nGIPerK.resize (nk);

      PsiG psi, psiTmp;

      iOffset = 0;
      for (ik=0; ik < nk; ik++)  {
         if (readKwise)  {
            offset = 0;
            psiVarName = "psi-" + SxString(ik+1);
         }
         if (stored == InMemory)
            waves(ik).resize (nSpin);
         if (mode & KeepNStates)
            nStates = nStatesPerK(ik);
         else
            nStatesPerK(ik) = nStates = nPerK(ik);
         ng = nGk(ik);
         int n = nStatesPerK(ik);
         nGIPerK(ik) = n * ng;
         fftIdx = SxVector<TPrecFFTIdx> ();
         fftIdx = toVector(fftIdxIn(SxIdx(iOffset, iOffset+ng*nComp-1)));
         iOffset += ng*nComp;
         if (hasBasis)  {
            if (mesh == (*gkBasisPtr)(ik).fft3d(0).mesh)  {
               // --- check n123 vs fftIdx
               SxVector<TPrecFFTIdx>::Iterator n123It, fftIdxIt;
               n123It   = (*gkBasisPtr)(ik).n123(0).begin ();
               fftIdxIt = fftIdx.begin ();
               int ig;
               for (ig = 0; ig < ng; ++ig)
                  if (*fftIdxIt++ != *n123It++) break;
               needsReorder = (ig < ng);
            } else {
               needsReorder = true;
               SxVector3<Int> vec;
               const SxFFT3d &fft = (*gkBasisPtr)(ik).fft3d(0);
               // --- check fftIdx fits into mesh
               for (int ig = 0; ig < ng*nComp; ++ig)  {
                  vec = mesh.getMeshVec (fftIdx(ig), SxMesh3D::Origin);
                  if (   2 * abs(vec(0)) > fft.mesh(0)
                      || 2 * abs(vec(1)) > fft.mesh(1)
                      || 2 * abs(vec(2)) > fft.mesh(2))  {
                     cout << "Mesh incompatibility between file and basis."
                          << endl;
                     SX_EXIT;
                  }
               }
            }
         }

         for (iSpin=0; iSpin < nSpin; iSpin++)  {

            sxprintf ("reading state (%d,%d)...%dx%d elements\n",
                  iSpin, ik, ng * nComp, n < nPerK(ik) ? n : nPerK(ik));
            fflush (stdout);

            // resize waves
            if (stored == InMemory)  {
               waves(ik)(iSpin).reformat (ng * nComp, n);
               if (gkBasisPtr.getPtr() != NULL)
                  waves(ik)(iSpin).setBasis (&(*gkBasisPtr)(ik));  // UGLY
            } else  {
               // stored == KeepOnDisk
               waves(0)(0).reformat (ng*nComp, n);
               loadedK = ik;
               loadedSpin = iSpin;
            }
            for (i=0; i < nPerK(ik); i++)  {
               if (i >= n)  {
                  // read fewer states than there are in file
                  // jump other states
                  offset += ng * (nPerK(ik) - n);
                  break;
               }

               psi   = (*this)(i,iSpin,ik);
               if (hasBasis && needsReorder)  {
                  psiTmp.resize (ng*nComp);
                  io.read (psiVarName, &psiTmp, ng*nComp, offset);
                  VALIDATE_VECTOR (psiTmp);
                  psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx, mesh);
               } else {
                  io.read (psiVarName, &psi, ng*nComp, offset);
                  VALIDATE_VECTOR (psi);
               }
               offset += ng*nComp;
            }

            // randomize rest of waves
            for (i = nPerK(ik); i < n; i++)  {
               (*this)(i, iSpin, ik).randomize ();
               (*this)(i, iSpin, ik).normalize ();
            }
            if (stored == KeepOnDisk)
               flushWaves ();

         } // iSpin

         // at the end, all coefficients must be read
         SX_CHECK ( (!readKwise) || offset == ng * nPerK(ik) * nSpin,
                   offset, ng * nComp*nPerK(ik));
      }
      // at the end, all coefficients must be read
      SX_CHECK (readKwise || offset == nCoeff, offset, nCoeff);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);

   //--- if GkBasis should be read from file create it from file and set
   if(!(mode & KeepGkBasis)) {
      //TODO CHECK MACRO for debugging only, overwrite GK basis is a bad idea but who knows when it is needed.
      SX_CHECK(gkBasisPtr.getPtr() == NULL);
      if (gkBasisPtr.getPtr() != NULL)  {
         cout << SX_SEPARATOR;
         cout << "WARNING: Overwrite GkBasisPtr in SxPW!" << endl;
         cout << SX_SEPARATOR;
      }
      if (mode & SaveMemory) gkBasisPtr = SxPtr<SxGkBasis>::create(io, true, true);
      else gkBasisPtr = SxPtr<SxGkBasis>::create(io, true, false);
      setGkBasisPtr(gkBasisPtr);
   }
}

/* TODO
void SxPW::readHDF5 (const SxString &file, int mode)
{
   SX_CHECK (stored != ReadOnDemand && stored != ReadOneByOne);
   if (stored == Unknown) stored = InMemory;
   SX_CHECK (stored == InMemory || mode & KeepNStates); // to be implemented

   bool hasBasis = (gkBasisPtr.getPtr() != NULL);
   bool needsReorder = hasBasis;

   cout << "SxPW::read: |G+k> basis ";
   if (!hasBasis) cout << "not ";
   cout << "present." << endl;

   // --- variable for k-wise read (allows larger files)
   bool readKwise;
   SxString psiVarName;

   //read HDF5 file
   hid_t waveFile = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

   try  {
      int ik, nk, iSpin, nSpinFile, ng, i, iOffset, nCoeff = -1;
      int nComp = 1;
      if(io.containsDim("nComp"))
         nComp = io.getDimension("nComp");
      nk        = io.getDimension ("nk");
      nSpinFile = io.getDimension ("nSpin");
      if (nSpinFile != nSpin && nSpin > 0)  {
         // can happen in non-DEBUG cases, so SX_CHECK not sufficient
         cout << "Inconsistency: this->nSpin = " << this->nSpin;
         cout << "; nSpin = " << nSpinFile << " in " << io.filename << endl;
         SX_EXIT;
      } else {
         nSpin = nSpinFile;
      }
      readKwise = !io.containsDim("nCoeff");
      if (readKwise)  {
         cout << "Reading waves k-wise" << endl;
      } else {
         nCoeff   = io.getDimension ("nCoeff");
         psiVarName = "psi";
      }
      sxprintf ("Reading wavefunction:\n");
      sxprintf ("---------------------\n");
      sxprintf ("nk: %d, nSpin: %d\n", nk, nSpin);

      SxString name;
      int offset = 0;
      int nAllGk = io.getDimension ("nAllGk");
      SxDiracVec<Int> nPerK (nk), nGk (nk), fftIdxIn(nAllGk), fftIdx;

      io.read ("nPerK",    &nPerK,  nk);
      io.read ("nGk",      &nGk,    nk);
      io.read ("fftIdx",   &fftIdxIn, nAllGk);
      // --- check if #k, k and #G(k) are correct
      if (hasBasis && nk != gkBasisPtr->nk)  {
         sxprintf ("Error: Input file is corrupt.\n");
         sxprintf ("       Number of k-points: %d (should be %d)\n", nk,
                 gkBasisPtr->nk);
         SX_EXIT;
      }
      if (hasBasis)  {
         SxMatrix<TPrecG> kVecs(nk,3); // need matrix for reading
         Coord kVec;
         io.read ("kVec", &kVecs, nk, 3);
         for (ik = 0; ik < nk; ik++)  {
            kVec = kVecs.row(ik).toVector3 ();
            if ( (kVec - gkBasisPtr->getK(ik)).normSqr () > 1e-10)  {
               cout << endl;
               cout << "Consistency check failed for k-point " << (ik+1);
               cout << ".\n Expected " << gkBasisPtr->getK(ik) << " but found ";
               cout << kVec<< " in '" << io.filename << "'." << endl;
               SX_EXIT; // may indicate changes in the k-point setup!
               // if this occurs, we should remap the k-points
            }
         }
      }
      for (ik=0; ik < nk; ik++)  {
         if (hasBasis && nGk(ik) != (*gkBasisPtr)(ik).ng)  {
            sxprintf ("Error: Input file is corrupt.\n");
            sxprintf ("       Number of |G+k> vectors: %d (should be %d)\n",
                     nGk(ik), (*gkBasisPtr)(ik).ng);
            SX_EXIT;
         }
      }

      SxMesh3D mesh;
      io.read("meshDim", &mesh);

      // --- resize to kpoints
      if (stored == InMemory)  {
         waves.resize (nk);
      }
      if (!(mode & KeepNStates))  nStatesPerK.resize (nk);
      SX_CHECK (nStatesPerK.getSize () == nk, nStatesPerK.getSize (), nk);

      nGIPerK.resize (nk);

      PsiG psi, psiTmp;

      iOffset = 0;
      for (ik=0; ik < nk; ik++)  {
         if (readKwise)  {
            offset = 0;
            psiVarName = "psi-" + SxString(ik+1);
         }
         if (stored == InMemory)
            waves(ik).resize (nSpin);
         if (mode & KeepNStates)
            nStates = nStatesPerK(ik);
         else
            nStatesPerK(ik) = nPerK(ik);
         ng = nGk(ik);
         int n = nStatesPerK(ik);
         nGIPerK(ik) = n * ng;
         fftIdx = SxDiracVec<Int> ();
         fftIdx = fftIdxIn(SxIdx(iOffset, iOffset+ng*nComp-1));
         iOffset += ng*nComp;
         if (hasBasis)  {
            if (mesh == (*gkBasisPtr)(ik).fft3d(0).mesh)  {
               // --- check n123 vs fftIdx
               SxDiracVec<Int>::Iterator n123It, fftIdxIt;
               n123It   = (*gkBasisPtr)(ik).n123(0).begin ();
               fftIdxIt = fftIdx.begin ();
               int ig;
               for (ig = 0; ig < ng; ++ig)
                  if (*fftIdxIt++ != *n123It++) break;
               needsReorder = (ig < ng);
            } else {
               needsReorder = true;
               SxVector3<Int> vec;
               const SxFFT3d &fft = (*gkBasisPtr)(ik).fft3d(0);
               // --- check fftIdx fits into mesh
               for (int ig = 0; ig < ng*nComp; ++ig)  {
                  vec = mesh.getMeshVec (fftIdx(ig), SxMesh3D::Origin);
                  if (   2 * abs(vec(0)) > fft.mesh(0)
                      || 2 * abs(vec(1)) > fft.mesh(1)
                      || 2 * abs(vec(2)) > fft.mesh(2))  {
                     cout << "Mesh incompatibility between file and basis."
                          << endl;
                     SX_EXIT;
                  }
               }
            }
         }

         for (iSpin=0; iSpin < nSpin; iSpin++)  {

            sxprintf ("reading state (%d,%d)...%dx%d elements\n",
                  iSpin, ik, ng * nComp, n < nPerK(ik) ? n : nPerK(ik));
            fflush (stdout);

            // resize waves
            if (stored == InMemory)  {
               waves(ik)(iSpin).reformat (ng * nComp, n);
               if (gkBasisPtr.getPtr() != NULL)
                  waves(ik)(iSpin).setBasis (&(*gkBasisPtr)(ik));  // UGLY
            } else  {
               // stored == KeepOnDisk
               waves(0)(0).reformat (ng*nComp, n);
               loadedK = ik;
               loadedSpin = iSpin;
            }
            for (i=0; i < nPerK(ik); i++)  {
               if (i >= n)  {
                  // read fewer states than there are in file
                  // jump other states
                  offset += ng * (nPerK(ik) - n);
                  break;
               }

               psi   = (*this)(i,iSpin,ik);
               if (hasBasis && needsReorder)  {
                  psiTmp.resize (ng*nComp);
                  io.read (psiVarName, &psiTmp, ng*nComp, offset);
                  VALIDATE_VECTOR (psiTmp);
                  psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx, mesh);
               } else {
                  io.read (psiVarName, &psi, ng*nComp, offset);
                  VALIDATE_VECTOR (psi);
               }
               offset += ng*nComp;
            }

            // randomize rest of waves
            for (i = nPerK(ik); i < n; i++)  {
               (*this)(i, iSpin, ik).randomize ();
               (*this)(i, iSpin, ik).normalize ();
            }
            if (stored == KeepOnDisk)
               flushWaves ();

         } // iSpin

         // at the end, all coefficients must be read
         SX_CHECK ( (!readKwise) || offset == ng * nPerK(ik) * nSpin,
                   offset, ng * nComp*nPerK(ik));
      }
      // at the end, all coefficients must be read
      SX_CHECK (readKwise || offset == nCoeff, offset, nCoeff);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);

   //--- if GkBasis should be read from file create it from file and set
   if(!(mode & KeepGkBasis)) {
      //TODO CHECK MACRO for debugging only, overwrite GK basis is a bad idea but who knows when it is needed.
      SX_CHECK(gkBasisPtr.getPtr() == NULL);
      if (gkBasisPtr.getPtr() != NULL)  {
         cout << SX_SEPARATOR;
         cout << "WARNING: Overwrite GkBasisPtr in SxPW!" << endl;
         cout << SX_SEPARATOR;
      }
      gkBasisPtr = SxPtr<SxGkBasis>::create(io);
      setGkBasisPtr(gkBasisPtr);
   }

   // Close HDF5 file
   H5Fclose(waveFile);
}
*/

void SxPW::write (SxBinIO &io) const
{
   int nk = getNk ();
   // for large numbers of kPoints reaching ncdump limit
   // set writekwise kPoint dependent
   bool writeKwise = nk < 1000;
   if (stored == KeepOnDisk && io.ncMode == SxBinIO::WRITE_DATA)
      flushWaves ();
   const SxGkBasis &Gk = *gkBasisPtr;
   int nComp = Gk.gBasis ? Gk.gBasis->getNComp() : 1;
   try  {
      int iSpin;
      int ik;
      int nCoeff = 0, i, n, ng;

      SxVector<Int> nPerK(nk), nGk(nk);
      for (ik=0; ik < nk; ik++)  {
         nPerK(ik) = n  = getNStates (ik);
         nGk(ik)   = ng = nGIPerK(ik) / n;
         nCoeff += n * ng * nSpin * nComp;
      }
      int offset = 0;

      // --- create dimensions
      io.addDimension ("nk",        nk);
      if (writeKwise) {
         for (ik = 0; ik < nk; ++ik)
            io.addDimension ("nCoeff-"+SxString(ik+1), nGIPerK(ik) * nSpin * nComp);
      } else {
         io.addDimension ("nCoeff",    nCoeff);
      }
      io.addDimension ("nSpin",     nSpin);
      io.addDimension ("nComp",     nComp);

      // --- write data
      if (!io.contains ("nGk") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nGk",   nGk,   "nk");
      if (!io.contains ("nPerK") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nPerK", nPerK, "nk");
      PsiG psi;
      SxString psiVarName = "psi",
               coeffDimName = "nCoeff";
      SX_NEW_LOOP (*this);
      for (ik=0; ik < nk; ik++)  {
         if (writeKwise)  {
            offset = 0;
            SxString sK(ik+1);
            psiVarName = "psi-" + sK;
            coeffDimName = "nCoeff-" + sK;
         }
         for (iSpin=0; iSpin < nSpin; iSpin++)  {

            if (   stored != InMemory && stored != Unknown
                && io.ncMode == SxBinIO::WRITE_DATA)
            {
               loadWaves (ik,iSpin);
            }

            for (i=0; i < nPerK (ik); i++)  {
               if (io.ncMode == SxBinIO::WRITE_DATA)  {
                  psi = (*this)(i, iSpin, ik);
               }
               io.write (psiVarName, psi, coeffDimName, offset);
               SX_CHECK (   io.ncMode != SxBinIO::WRITE_DATA
                         || psi.getSize () == nGk(ik) * nComp,
                         psi.getSize (), nGk(ik) * nComp);
               offset += nGk(ik) * nComp;
            }
         }
         // at the end, all coefficients must be written
         SX_CHECK ((!writeKwise) || offset == nGIPerK(ik) * nSpin * nComp,
                   offset, nGIPerK(ik) * nSpin * nComp);
      }

      // at the end, all coefficients must be written
      SX_CHECK (writeKwise || offset == nCoeff, offset, nCoeff);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


size_t SxPW::getNBytes () const
{
   return ::getNBytes (waves)
        + ::getNBytes (nStates)
        + ::getNBytes (nGIPerK);
}

void SxPW::registerMemoryObservers ()
{
   TRACK_MEMORY (waves);
   TRACK_MEMORY (nStates);
   TRACK_MEMORY (nGIPerK);
   TRACK_MALLOC (*this, 1);
}

void SxPW::setZero()
{
   SX_CHECK(stored == InMemory);
   for (int ik = 0; ik < getNk(); ik++)
      for (int iSpin = 0; iSpin < getNSpin(); iSpin++)
            waves(ik)(iSpin).set(0.);

}

#ifdef SXPW_FINDLOOP
void SxPW::findLoopUpdate (int iSpin, int ik) const
{
   if (iSpin == findLoopLastIk && ik == findLoopLastIk) return;
   if (iSpin == 0 && ik == 0)  {
      cout << "NEW waves loop!" << endl;
   }
   findLoopLastIk = ik; findLoopLastIspin = iSpin;
}

void SxPW::findLoopReset () const
{
   findLoopLastIk = findLoopLastIspin = 0;
}
#endif

