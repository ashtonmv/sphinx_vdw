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

#include <SxPAWSet.h>
#include <SxLoopMPI.h> // LoopMPI
#include <SxProjector.h>

SxPAWSet::SxPAWSet (SxPtr<SxGkBasis> gkPtrIn,
                    const SxPtr<SxPartialWaveBasis> &pBasis,
                    int n, int nSpinIn, const SxString &tmpDir)
: SxPWSet (gkPtrIn),
  nStates(n), nSpin (nSpinIn),
  keepWavesOnDisk (tmpDir.getSize () > 0),
  loadedK(-1), loadedSpin(-1)
{
   SxGkBasis &gk = *gkPtrIn;
   int nk = gk.getNk ();
   pawBasis.resize (nk);
   if (keepWavesOnDisk)  {
      waves.reformat (1, 1);
   } else {
      waves.reformat (nk, nSpin);
   }
   for (int ik = 0; ik < nk; ++ik)  {
      pawBasis(ik) = SxPtr<SxPAWBasis>::create (gk.gBasisList(ik), pBasis);
      if (!keepWavesOnDisk)  {
         int N = int(pawBasis(ik)->getNElements ());
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik))
            {  // LoopMPI
               waves(ik, iSpin) = SxDiracMat<Complex16> (N, nStates);
               waves(ik, iSpin).handle->auxData.ik    = ik;
               waves(ik, iSpin).handle->auxData.iSpin = iSpin;
               waves(ik, iSpin).setBasis (pawBasis(ik).getPtr ());
            }  // LoopMPI
         }
      }
   }
   if (keepWavesOnDisk) createScratchFile (tmpDir);
}

SxPAWSet::SxPAWSet (SxPtr<SxGkBasis> gkPtrIn,
                    const SxArray<SxPtr<SxPAWBasis> > pawBasisIn,
                    int n, int nSpinIn, const SxString &tmpDir)
 : SxPWSet (gkPtrIn),
   nStates(n), nSpin (nSpinIn),
   keepWavesOnDisk (tmpDir.getSize () > 0),
   loadedK(-1), loadedSpin(-1)
{
   SxGkBasis &gk = *gkPtrIn;
   int nk = gk.getNk ();
   if (keepWavesOnDisk)  {
      waves.reformat (1,1);
   } else {
      waves.reformat (nk, nSpin);
   }
   pawBasis.resize (nk);
   for (int ik = 0; ik < nk; ++ik)  {
      pawBasis(ik) = pawBasisIn(ik);
      SX_MPI_LEVEL("waves-k");
      if (!keepWavesOnDisk && SxLoopMPI::myWork (ik))  {
         int N = int(pawBasis(ik)->getNElements ());
         for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
            waves(ik, iSpin) = SxDiracMat<Complex16> (N, nStates);
            waves(ik, iSpin).handle->auxData.ik    = ik;
            waves(ik, iSpin).handle->auxData.iSpin = iSpin;
            waves(ik, iSpin).setBasis (pawBasis(ik).getPtr ());
         }
      }
   }
   if (keepWavesOnDisk) createScratchFile (tmpDir);
}

SxPAWSet::SxPAWSet (SxConstPtr<SxPAWPot> potPtr, 
                    const SxAtomicStructure &structure, 
                    const SxBinIO &io)
 : nStates(-1), 
   nSpin (-1),
   keepWavesOnDisk (false),
   loadedK(-1), 
   loadedSpin(-1),
   wavesFile(io)
{
   readPAWBasis(wavesFile, potPtr, structure);
   
   read(wavesFile);
}

SxPAWSet::~SxPAWSet ()
{
   if (keepWavesOnDisk)  {
      SxString tmpFilename = wavesFile.filename;
      try {
         wavesFile.close ();
      } catch (SxException e) {
         e.print ();
         SX_EXIT;
      }
      unlink (tmpFilename.ascii());
   }
}


void SxPAWSet::setGkBasisPtr (SxPtr<SxGkBasis> inPtr)
{
   gkBasisPtr = inPtr;
}

SxPtr<SxPWSet> SxPAWSet::getNew () const
{
   SxString tmpDir;
   if (keepWavesOnDisk)
      tmpDir = SxFileInfo (wavesFile.filename).getPath ();
   return SxPtr<SxPAWSet>::create (gkBasisPtr, pawBasis, nStates, 
                                   getNSpin (), tmpDir);
}

const SxGBasis::TPsi SxPAWSet::operator() (int i, int iSpin, int ik) const
{
   SX_CHECK (i >= 0 && i < nStates, i, nStates);
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   PsiG res;
   if (keepWavesOnDisk)  {
      if (loadedK != ik || loadedSpin != iSpin)  {
         flushWaves ();
         loadWaves(iSpin, ik);
      }
      res = waves(0,0).colRef (i);
   } else {
      res = waves(ik, iSpin).colRef (i);
   }
   res.handle->auxData.i     = i;
   res.handle->auxData.iSpin = iSpin;
   res.handle->auxData.ik    = ik;
   res.setBasis (&getBasis(ik));
   return res;
}

SxGBasis::TPsi SxPAWSet::operator() (int i, int iSpin, int ik)
{
   SX_CHECK (i >= 0 && i < nStates, i, nStates);
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   PsiG res;
   if (keepWavesOnDisk)  {
      if (loadedK != ik || loadedSpin != iSpin)  {
         flushWaves ();
         loadWaves(iSpin, ik);
      }
      res = waves(0,0).colRef (i);
   } else {
      res = waves(ik, iSpin).colRef (i);
   }
   res.handle->auxData.i     = i;
   res.handle->auxData.iSpin = iSpin;
   res.handle->auxData.ik    = ik;
   res.setBasis (&getBasis(ik));
   return res;
}

const SxGBasis::TPsi& SxPAWSet::operator() (int iSpin, int ik) const
{
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   if (keepWavesOnDisk)  {
      if (loadedK != ik || loadedSpin != iSpin)  {
         flushWaves ();
         loadWaves(iSpin, ik);
      }
      return waves(0,0);
   }
   return waves(ik, iSpin);
}

SxGBasis::TPsi& SxPAWSet::operator() (int iSpin, int ik)
{
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   if (keepWavesOnDisk)  {
      if (loadedK != ik || loadedSpin != iSpin)  {
         flushWaves ();
         loadWaves(iSpin, ik);
      }
      return waves(0,0);
   }
   return waves(ik, iSpin);
}

PsiG SxPAWSet::getBlock (int iStart, int nBlock, int iSpin, int ik)
{
   SX_CHECK (nBlock > 0, nBlock);
   SX_CHECK (iStart >= 0 && iStart + nBlock <= nStates,
             iStart, nBlock, nStates);
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   PsiGI &waveRef = (*this)(iSpin, ik); // note: might be from disk
   int N = (int)waveRef.nRows ();
   PsiG res = waveRef(SxIdx (iStart * N, (iStart+nBlock)*N - 1));
   res.reshape (N, nBlock);
   res.handle->auxData.ik = ik;
   res.handle->auxData.iSpin = iSpin;
   res.setBasis (&getBasis(ik));
   return res;
}

void SxPAWSet::createScratchFile (const SxString &tmpDir)
{
   // --- create scratch file
   const SxString fileName = SxBinIO::createTempName(tmpDir);
   wavesFile.open (fileName, SxBinIO::SCRATCH_READ_WRITE);
   // --- fill file buffer
   for (int ik = 0; ik < getNk (); ++ik)  {
      size_t size = pawBasis(ik)->getNElements () * nStates;
      SxVector<Complex16> psi(size);
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         fwrite (psi.elements, sizeof(SxComplex16), size, wavesFile.fp);
      }
   }
}

void SxPAWSet::loadWaves (int iSpin, int ik) const
{
   SX_CHECK (keepWavesOnDisk);
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   //cout << "Loading iSpin=" << iSpin << ", ik=" << ik << endl;

   // --- prepare memory
   waves(0,0) = SxDiracMat<Complex16> (); // release memory (might be buffered)
   waves(0,0) = SxDiracMat<Complex16> (pawBasis(ik)->getNElements (),
                                         nStates);
   SxDiracAux &aux = waves(0,0).handle->auxData;
   aux.iSpin = iSpin;
   aux.ik    = ik;
   waves(0,0).setBasis (&getBasis(ik));
   loadedSpin = iSpin;
   loadedK    = ik;

   // --- read waves
   size_t pos = 0, size_k = pawBasis(ik)->getNElements () * nStates;
   for (int jk = 0; jk < ik; ++jk) 
      pos += pawBasis(jk)->getNElements () * nStates * getNSpin ();
   pos += size_k * iSpin;
   fseek (wavesFile.fp, pos * sizeof(SxComplex16), SEEK_SET);
   SX_CHECK (waves(0,0).getSize () == ssize_t(size_k),
             waves(0,0).getSize (), size_k);
   size_t nRead = fread (waves(0,0).elements, sizeof(SxComplex16), size_k,
                          wavesFile.fp);
   if (nRead != sizeof(SxComplex16) * size_k) { SX_EXIT; }
}

void SxPAWSet::flushWaves () const
{
   SX_CHECK (keepWavesOnDisk);
   if (loadedK < 0 || loadedSpin < 0) return;
   SX_CHECK (loadedK >= 0 && loadedK < getNk (), loadedK, getNk ());
   SX_CHECK (loadedSpin >= 0 && loadedSpin < getNSpin (),
             loadedSpin, getNSpin ());
   //cout << "Flushing iSpin=" << loadedSpin << ", ik=" << loadedK << endl;

   // --- write waves
   size_t pos = 0, size_k = pawBasis(loadedK)->getNElements () * nStates;
   for (int jk = 0; jk < loadedK; ++jk) 
      pos += pawBasis(jk)->getNElements () * nStates * getNSpin ();
   pos += size_k * loadedSpin;
   fseek (wavesFile.fp, pos * sizeof(SxComplex16), SEEK_SET);
   SX_CHECK (waves(0,0).getSize () == ssize_t(size_k),
             waves(0,0).getSize (), size_k);
   fwrite (waves(0,0).elements, sizeof(SxComplex16), size_k, wavesFile.fp);
}

void SxPAWSet::memMinimize ()
{
   if (keepWavesOnDisk)  {
      flushWaves ();
      loadedK = loadedSpin = -1;
      waves(0,0) = SxDiracMat<TPrecCoeffG> ();
   }
}

// TODO   Check support for parallel NetCDF4 / SxParallelHierarchy
void SxPAWSet::write (SxBinIO &io) const
{
   int nk = getNk ();

   bool writeKwise;
   if (io.mode == SxBinIO::BINARY_WRITE_PARALLEL)
      writeKwise = true;
   else
      writeKwise = nk < 200;

   if (keepWavesOnDisk && io.ncMode == SxBinIO::WRITE_DATA)
      flushWaves ();

   const SxGkBasis &Gk = *gkBasisPtr;
   SX_CHECK (Gk.gBasis->getNComp() == 1, Gk.gBasis->getNComp ());
   int nComp = 1;
   int nProj = (int)pawBasis(0)->pBasis->getNElements ();
   try  {
      int nCoeff = 0, nProjAll = 0;

      SxVector<Int> nPerK(nk), nGk(nk), nGIPerK(nk);
      for (int ik=0; ik < nk; ik++)  {
         int n, ng;
         nPerK(ik) = n  = getNStates (ik);
         nGk(ik)   = ng = Gk(ik).ng;
         nGIPerK(ik) = n * ng;
         nCoeff += n * ng * nSpin;
         nProjAll += n * nProj * nSpin;
      }

      // --- create dimensions
      io.addDimension ("nk",        nk);
      if (writeKwise) {
         for (int ik = 0; ik < nk; ++ik)
            io.addDimension ("nCoeff-"+SxString(ik+1), nGIPerK(ik) * nSpin);
      } else {
         io.addDimension ("nCoeff",    nCoeff);
      }
      io.addDimension ("nSpin",     nSpin);
      io.addDimension ("nComp",     nComp);
      io.addDimension ("nProj",     nProj);
      io.addDimension ("nProjAll",  nProjAll);

      // --- write data
      int offset = 0, offsetP = 0;
      if (!io.contains ("nGk") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nGk",   nGk,   "nk");
      if (!io.contains ("nPerK") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nPerK", nPerK, "nk");
      PsiG psi;
      SxString psiVarName = "psi",
               coeffDimName = "nCoeff";
      //SX_NEW_LOOP (*this);
      SX_MPI_LEVEL("waves-k");
      for (int ik=0; ik < nk; ik++)  {
         if (writeKwise)  {
            offset = 0;
            SxString sK(ik+1);
            psiVarName = "psi-" + sK;
            coeffDimName = "nCoeff-" + sK;
         }
         SxIdx gkPart(0, nGk(ik) - 1),
               pPart (nGk(ik), (int)pawBasis(ik)->getNElements () - 1);
         for (int iSpin=0; iSpin < nSpin; iSpin++)  {

            if (keepWavesOnDisk && io.ncMode == SxBinIO::WRITE_DATA)
               loadWaves (iSpin, ik);
         
            for (int i=0; i < nPerK (ik); i++)  {
               if (io.ncMode == SxBinIO::WRITE_DATA)  {
                  psi = (*this)(i, iSpin, ik);
               } else {
                  psi.resize (pawBasis(ik)->getNElements ());
               }
               if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::myWork(ik)))
               {
                  io.write ("projPsi", psi(pPart), "nProjAll", offsetP);
                  io.write (psiVarName, psi(gkPart), coeffDimName, offset);
               }
               offset  += nGk(ik);
               offsetP += nProj;
            }
         }
         // at the end, all coefficients must be written
         SX_CHECK ((!writeKwise) || offset == nGIPerK(ik) * nSpin,
                   offset, nGIPerK(ik) * nSpin);
      }
      
      // at the end, all coefficients must be written
      SX_CHECK (writeKwise || offset == nCoeff, offset, nCoeff);
      SX_CHECK (offsetP == nProjAll, offsetP, nProjAll);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxPAWSet::read (const SxBinIO &io, int mode)
{
   SX_CHECK (gkBasisPtr.getPtr() != NULL);
   SX_CHECK (pawBasis.getSize () > 0);

   bool needsReorder = true;

#ifndef NDEBUG
   SX_CHECK (gkBasisPtr->getNk () == pawBasis.getSize (),
             gkBasisPtr->getNk (), pawBasis.getSize ());
   for (int ik = 0; ik < pawBasis.getSize (); ++ik)  {
      SX_CHECK (&(*gkBasisPtr)(ik) == pawBasis(ik)->gBasis.getPtr ());
   }
#endif

   //--- if GkBasis should be read from file create it from file and set
   if(!(mode & KeepGkBasis)) {
      cout << SX_SEPARATOR;
      cout << "WARNING: Overwrite GkBasisPtr in SxPW!" << endl;
      cout << SX_SEPARATOR;
      gkBasisPtr = SxPtr<SxGkBasis>::create(io);
   }

   // --- variable for k-wise read (allows larger files)
   bool readKwise;
   SxString psiVarName;
   bool hasProj = false;
   
   try  {
      int nk, nSpinFile, ng;
#ifndef NDEBUG
      int nCoeff = -1, nProjAll = -1;
#endif
      int nComp = 1;
      if(io.containsDim("nComp"))
         nComp = io.getDimension("nComp");
      if (nComp != 1)  {
         cout << "No support for nComp != 1 in PAW basis" << endl;
         SX_EXIT;
      }
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
      if ((readKwise = !io.containsDim("nCoeff")))  {
         cout << "Reading waves k-wise" << endl;
      } else {
#ifndef NDEBUG
         nCoeff   = io.getDimension ("nCoeff");
#endif
         psiVarName = "psi";
      }

      hasProj = io.containsDim ("nProj");

      int nProj = (int)pawBasis(0)->pBasis->getNElements ();
      if (hasProj)  {
         int nProjFile = io.getDimension ("nProj");
         if (nProjFile != nProj)  {
            cout << "Inconsistent number of projectors" << endl;
            cout << io.filename << " has " << nProjFile << " projectors." << endl;
            cout << "current SxPartialWaveBasis has " << nProj << endl;
            SX_QUIT;
         }
#ifndef NDEBUG
         nProjAll = io.getDimension ("nProjAll");
#endif
      }
      sxprintf ("Reading wavefunction:\n");
      sxprintf ("---------------------\n");
      sxprintf ("nk: %d, nSpin: %d\n", nk, nSpin);
      
      SxString name;
      int nAllGk = io.getDimension ("nAllGk");
      SxDiracVec<Int> nPerK (nk), nGk (nk), fftIdxIn(nAllGk);
      SxVector<TPrecFFTIdx> fftIdx;

      io.read ("nPerK",    &nPerK,  nk);
      io.read ("nGk",      &nGk,    nk);
      io.read ("fftIdx",   &fftIdxIn, nAllGk);
      // --- check if #k, k and #G(k) are correct
      if (nk != gkBasisPtr->nk)  {
         sxprintf ("Error: Input file is corrupt.\n");
         sxprintf ("       Number of k-points: %d (should be %d)\n", nk,
                 gkBasisPtr->nk);
         SX_EXIT;
      }
      SxMatrix<TPrecG> kVecs(nk,3); // need matrix for reading
      Coord kVec;
      io.read ("kVec", &kVecs, nk, 3);
      for (int ik = 0; ik < nk; ik++)  {
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
      for (int ik=0; ik < nk; ik++)  {
         SX_MPI_LEVEL("waves-k");
         if (!SxLoopMPI::myWork(ik)) continue; // SxParallelHierarchy
         if (nGk(ik) != (*gkBasisPtr)(ik).ng)  {
            sxprintf ("Error: Input file is corrupt.\n");
            sxprintf ("       Number of |G+k> vectors: %d (should be %d)\n",
                     nGk(ik), (*gkBasisPtr)(ik).ng);
            SX_EXIT;
         }
      }

      SxMesh3D mesh;
      io.read("meshDim", &mesh);

      // --- resize to kpoints
      if (!keepWavesOnDisk)  {
         waves.reformat (nk, nSpin);
      }

      PsiG psi, psiTmp;
      
      int iOffset = 0;
      int offset = 0, offsetP = 0;

      // set nStates from first k-point 
      if (!(mode & KeepNStates)) nStates = nPerK(0);

      for (int ik=0; ik < nk; ik++)  {
         SX_MPI_LEVEL("waves-k");
         // --- check number of states
         if (!(mode & KeepNStates))  {
            if (nStates != nPerK(ik))  {
               cout << "PAW set does not support varying number of states per "
                    << "k-point" << endl;
               cout << "Incompatability while reading waves file "
                    << io.filename << endl; 
               cout << "k-point 1 has " << nStates << " states" << endl;
               cout << "k-point " << (ik+1) << " has " << nPerK(ik) 
                    << " states " << endl;
               SX_EXIT;
            }
         }
         // --- move on if this k-point is not for this MPI task
         if (!SxLoopMPI::myWork(ik))  {
            iOffset += nGk(ik)*nComp;
            if (!readKwise) offset += nGk(ik) * nSpin * nPerK(ik);
            if (hasProj) offsetP += nProj * nSpin * nPerK(ik);
            continue;
         }
         if (readKwise)  {
            offset = 0;
            psiVarName = "psi-" + SxString(ik+1);
         }
         ng = nGk(ik);
         fftIdx = SxVector<TPrecFFTIdx> ();
         fftIdx = toVector(fftIdxIn(SxIdx(iOffset, iOffset+ng*nComp-1)));
         iOffset += ng*nComp;
         if (mesh == (*gkBasisPtr)(ik).fft3d(0).mesh)  {
            // --- check n123 vs fftIdx
            SxVector<Int>::Iterator n123It, fftIdxIt;
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

         for (int iSpin=0; iSpin < nSpin; iSpin++)  {
            
            if (hasProj)
               sxprintf ("reading state (%d,%d)...(%d+%d)x%d elements\n", iSpin, ik,
                     ng * nComp, int(pawBasis(ik)->pBasis->getNElements ()),
                     nStates < nPerK(ik) ? nStates : nPerK(ik));
            else 
               sxprintf ("reading state (%d,%d)...%dx%d elements\n", iSpin, ik,
                     ng * nComp, nStates < nPerK(ik) ? nStates : nPerK(ik));
            fflush (stdout);

            // resize waves
            if (keepWavesOnDisk)  {
               // stored == KeepOnDisk
               waves(0,0).reformat (ng*nComp, nStates);
               loadedK = ik;
               loadedSpin = iSpin;
            } else  {
               waves(ik, iSpin).reformat (pawBasis(ik)->getNElements (), nStates);
               waves(ik, iSpin).setBasis (*pawBasis(ik));
               waves(ik, iSpin).handle->auxData.ik = ik;
               waves(ik, iSpin).handle->auxData.iSpin = iSpin;
            }
            for (int i=0; i < nPerK(ik); i++)  {
               if (i >= nStates)  {
                  // read fewer states than there are in file
                  // jump other states
                  offset += ng * (nPerK(ik) - nStates);
                  break;
               }
                  
               psi   = (*this)(i,iSpin,ik);
               psiTmp.resize (ng*nComp);
               io.read (psiVarName, &psiTmp, ng*nComp, offset);
               offset += ng*nComp;
               if (needsReorder)  {
                  VALIDATE_VECTOR (psiTmp);
                  psiTmp <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx, mesh);
               }
               if (hasProj)  {
                  PsiG proj(*pawBasis(ik)->pBasis);
                  io.read ("projPsi", &proj, (int)proj.getSize (), offsetP);
                  offsetP += (int)proj.getSize ();
                  psi(SxIdx(0, ng-1)) <<= psiTmp;
                  psi(SxIdx(ng, (int)psi.getSize () - 1)) <<= proj;
               } else {
                  psiTmp.setBasis ((*gkBasisPtr)(ik));
                  psiTmp.handle->auxData.ik = ik;
                  psi <<= *pawBasis(ik) | psiTmp; 
               }
               VALIDATE_VECTOR (psi);
            }

            // randomize rest of waves
            for (int i = nPerK(ik); i < nStates; i++)  {
               (*this)(i, iSpin, ik).randomize ();
               (*this)(i, iSpin, ik).normalize ();
            }
            if (keepWavesOnDisk)
               flushWaves ();

         } // iSpin

         // at the end, all coefficients must be read
         SX_CHECK ((!readKwise) || offset == ng * nPerK(ik) * nSpin,
                   offset, ng * nComp*nPerK(ik));
      }
      // at the end, all coefficients must be read
      SX_CHECK (readKwise || offset == nCoeff, offset, nCoeff);
      SX_CHECK (!hasProj || offsetP == nProjAll, offsetP, nProjAll);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);
}

void SxPAWSet::readPAWBasis (const SxBinIO &io, SxConstPtr<SxPAWPot> potPtr, const SxAtomicStructure &structure)
{
   // 
   SX_CHECK (gkBasisPtr.getPtr () == NULL);
   SX_CHECK (pawBasis.getSize () == 0);

   // Initialize Gk Basis + set structure
   gkBasisPtr = SxPtr<SxGkBasis>::create (io);
   gkBasisPtr-> changeTau (structure);

   // Initialize partial wave basis + generate projectors
   SxPtr<SxPartialWaveBasis> pBasis = 
      SxPtr<SxPartialWaveBasis>::create (potPtr, structure);
   pBasis->projectors = pBasis->createProjBasis (*gkBasisPtr, *potPtr);

   // Initialize PAW basis
   int nk = gkBasisPtr->getNk ();
   pawBasis.resize (nk);
   for (int ik = 0; ik < nk; ++ik)  {  
      pawBasis(ik) = SxPtr<SxPAWBasis>::create (gkBasisPtr->gBasisList(ik), pBasis);
   }
}

