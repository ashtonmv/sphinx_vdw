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

#include <SxGBasis.h>
#include <SxConstants.h>
#include <SxError.h>
//#include <main.h>
#include <math.h>
#include <SxBinIO.h>
#include <SxRBasis.h>
#include <SxAOBasis.h>
#include <SxPartialWaveBasis.h>
#include <SxPAWBasis.h>

namespace Timer {
   enum GBasisTimer {
      GtoRad,
      GBessel,
      Convolute,
      ConvoluteFFT,
      AddToRho,
      AddToRhoFFT
   };
}

SX_REGISTER_TIMERS (Timer::GBasisTimer)
{
   using namespace Timer;
   regTimer (GtoRad,       "GtoRad");
   regTimer (GBessel,      "GBessel");
   regTimer (Convolute,    "convolute");
   regTimer (ConvoluteFFT, "convolute FFT");
   regTimer (AddToRho, "rho(2+1D)");
   regTimer (AddToRhoFFT, "rho(2+1D) FFT");
}


SxGBasis::SxGBasis (enum MemoryUsage memModeIn)
   : rBasisPtr(NULL)
{
   memMode = memModeIn;
   ng = -1;
   structPtr = NULL;
   useTimeRevSym = false;
   gCut = 0.;
   registerMemoryObservers ();
   nComp = 1;
}
SxGBasis::SxGBasis (const SxVector3<Int> &mesh,
                    const SxAtomicStructure &structIn,
                    PrecEnergy gCut_,
                    const SxVector3<TPrecG> &dG_,
                    bool useTimeRevSymIn,
                    enum MemoryUsage memModeIn)
{
   set(mesh, structIn, gCut_, dG_, useTimeRevSymIn, memModeIn);
}

SxGBasis::SxGBasis (const SxGBasis &origBasis, int nCompNew)
{
ng = origBasis.ng;
memMode = origBasis.memMode;
structPtr = origBasis.structPtr;
gCut = origBasis.gCut;
dG = origBasis.dG;
useTimeRevSym = origBasis.useTimeRevSym;
nComp = nCompNew;
rBasisPtr = NULL;

g2.copy (origBasis.g2);
g2.setBasis (this);
gVec.copy (origBasis.gVec);
gVec.setBasis (this);
structureFactors = origBasis.structureFactors;

fft3d.resize (1);
fft3d(0) = origBasis.fft3d(0);
n123.resize (1);
n123(0).copy (origBasis.n123(0));
}

void SxGBasis::set (const SxVector3<Int> &mesh,
                    const SxAtomicStructure &structIn,
                    PrecEnergy gCut_,
                    const SxVector3<TPrecG> &dG_,
                    bool useTimeRevSymIn,
                    enum MemoryUsage memModeIn)
{
   memMode = memModeIn;
   rBasisPtr = NULL;
   structPtr = &structIn;
   if (fft3d.getSize () == 0) fft3d.resize (1);
   SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
   fft3d(0) = SxFFT3d (SxFFT3d::Forward, mesh(0), mesh(1), mesh(2),
                       structIn.cell.volume, true),
   fft3d(0).autoscale = false;
   if (memMode == SaveTime)
      fft3d(0).clean ();
   else
      fft3d(0).destroyArrays ();
   gCut = gCut_;
   dG   = dG_;
   if (useTimeRevSymIn && (dG^dG) < 1e-10 )  useTimeRevSym = true;
   else                                      useTimeRevSym = false;

   compute ();
   nComp = 1;
   registerMemoryObservers ();
}

void SxGBasis::set (const SxVector3<Int> &mesh,
                    const SxCell &cell,
                    PrecEnergy gCut_,
                    const SxVector3<TPrecG> &dG_,
                    bool useTimeRevSymIn,
                    MemoryUsage memModeIn)
{
   memMode = memModeIn;
   SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
   if (fft3d.getSize () == 0) fft3d.resize (1);
   fft3d(0).symmetric = true;
   fft3d(0).dir = SxFFT3d::Forward;
   fft3d(0).setMesh (mesh(0), mesh(1), mesh(2), cell.volume);
   if (memMode == SaveTime)
      fft3d(0).clean ();
   else
      fft3d(0).destroyArrays ();
   fft3d(0).autoscale = false;
   gCut = gCut_;
   dG = dG_;
   useTimeRevSym = useTimeRevSymIn && ((dG^dG) < 1e-10);

   structPtr = NULL;
   phaseFactors.resize     (0);
   structureFactors.resize (0);
   update (cell);

   nComp = 1;
   registerMemoryObservers ();
}

SxGBasis::~SxGBasis ()
{
   deregisterAll ();
}

void SxGBasis::setNComp (int comp)
{
   nComp = comp;
}

int SxGBasis::getNComp () const
{
   return nComp;
}

Real8 SxGBasis::getECut (const SxSymbolTable *table)
{
   Real8 eCut = 0.0;

   SX_CHECK (table);
   try  {
      eCut = table->getGroup("basis")->get("eCut")->toReal();
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return eCut;
}

Real8 SxGBasis::getGCut (const SxSymbolTable *table)
{
   Real8 gCut = 0.0;

   SX_CHECK (table);
   try  {
      const SxSymbolTable *basisGrp = table->getGroup("basis");
      if (basisGrp->contains("gCut"))  {
         gCut = basisGrp->get("gCut")->toReal ();
      } else {
         gCut = getGCut(basisGrp->get("eCut")->toReal());
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return gCut;
}

SxVector3<Int> SxGBasis::getMesh (const SxSymbolTable *table)
{
   SxCell cell (table);
   if (!cell.symGroupPtr)  {
      cell.symGroupPtr 
         = SxPtr<SxSymGroup>::create(cell.getLatticeSymmetries ());
   }
   return getMesh(table, cell);
}


SxVector3<Int> SxGBasis::getMesh (const SxSymbolTable *table,
                                  const SxCell &cell)
{
   SxVector3<Int> fftMeshSize;

   Real8 eCut = 0.;
   Real8 meshAccuracy = 1.0;
   try {
      eCut = SxGBasis::getECut (table);

      SxSymbolTable *basis = table->getGroup("basis");
      if (basis->contains("meshAccuracy"))
         meshAccuracy = basis->get("meshAccuracy")->toReal();

      if (basis->contains("mesh"))  {
         fftMeshSize = SxVector3<Int> (basis->get("mesh")->toIntList());
         if (!isCommensurableMesh (fftMeshSize, cell, true))  {
            sxprintf ("ERROR: The provided mesh is not commensurable!\n");
            getCommensurableMesh (eCut, cell, meshAccuracy);
            SX_QUIT;
         }
         return fftMeshSize;
      }


   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return getCommensurableMesh (eCut, cell, meshAccuracy);

}

SxVector3<Int> SxGBasis::getMeshSize (Real8 eCut,
                                      const SxMatrix3<TPrecTauR> &aMat,
                                      Real8 meshAccuracy)
{
   SxVector3<Int> fftMeshSize;
   int size;
   for (int d=0; d < 3; d++)  {
      size = (int)ceil(2. * sqrt(eCut * aMat(d).absSqr().sum()) / PI
                       * meshAccuracy);
      fftMeshSize(d) = SxFFT::getNextMeshSize(size);
   }
   return fftMeshSize;
}


SxVector3<Int> SxGBasis::getCommensurableMesh (Real8 eCut,
                                               const SxCell &cell,
                                               double meshAccuracy)
{
   return getCommensurableMesh (cell, getMeshSize (eCut, cell, meshAccuracy));
}

SxVector3<Int> SxGBasis::getCommensurableMesh (const SxCell &cell,
                                               SxVector3<Int> meshIn)
{
   SxVector3<Int> mesh = meshIn;
   SxVector3<Int> minMesh(0,0,0);
   //cout << "FFT-optimal mesh = " << mesh << endl;
   int i, j, k;
   const int maxTrials = 10;
   if ( !isCommensurableMesh (mesh, cell, true) )  {
      cout << "Commensurable mesh search..." << endl;
      for (i=0; i < maxTrials; ++i)  {
         for (j=0; j < maxTrials; ++j)  {
            for (k=0; k < maxTrials; ++k)  {
               if ( isCommensurableMesh(mesh, cell) )  {
                  cout << "Possible mesh: " << mesh << " ("
                       << mesh.product() << " elements)\n";
                  if (minMesh(0) == 0 || minMesh.product () > mesh.product ())
                    minMesh = mesh;
               }
               mesh(2) = SxFFT::getNextMeshSize(mesh(2)+1);
            }
            mesh(1) = SxFFT::getNextMeshSize(mesh(1)+1);
            mesh(2) = meshIn(2);
         }
         mesh(0) = SxFFT::getNextMeshSize(mesh(0)+1);
         mesh(1) = meshIn(1);
      }
      if (minMesh(0) == 0)  {
         cout << "No evenly-spaced commensurable mesh found." << endl;
         mesh = meshIn;
         for (int d = 0; d < 3; ++d)  {
            cout << "Tried mesh(" << d << ")=";
            for (i=0; i < maxTrials; ++i)  {
               cout << " " << mesh(d);
               mesh(d) = SxFFT::getNextMeshSize(mesh(d)+1);
            }
            cout << endl;
         }
         minMesh = max(max(meshIn(0),meshIn(1)), meshIn(2));
         if (!isCommensurableMesh(minMesh, cell))  {
            // scaling a cell cannot change its symmetries!
            cout << "Symmetry inconsistency error" << endl;
            SX_EXIT;
         }
         cout << "Fall-back to regular mesh " << minMesh << endl;
      }
   } else {
      minMesh = mesh;
   }

   return minMesh;
}

bool SxGBasis::isCommensurableMesh (const SxVector3<Int> &mesh,
                                    const SxCell         &cell,
                                    bool  verbose_)
{
   if (!cell.symGroupPtr) return true;
   SxSymGroup &S = *cell.symGroupPtr;
   int iOp, nOp = S.getNSymmorphic ();
   Coord x;
   SxCell meshCell (cell(0)/double(mesh(0)),
                    cell(1)/double(mesh(1)),
                    cell(2)/double(mesh(2)));

   bool symOk = true;
   for (iOp=0; iOp < nOp; iOp++)  {
      if (!meshCell.isLatticeSymmetry (S.getSymmorphic(iOp)))  {
         if (!verbose_)  return false;
         else  {
            cout << "Incompatible symmetry " << S.getSymmorphic(iOp) << endl;
            cout << meshCell.carToRel (S.getSymmorphic(iOp)) << endl;
            symOk = false;
         }
      }
   }
   return symOk;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxGBasis::compute ()
{
   SX_CLOCK (Timer::GBasis);
   SX_CHECK (structPtr);
   const SxAtomicStructure &str = *structPtr;
   update ();

   changeTau (str);
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxGBasis::update ()
{
   //useTimeRevSym = false;
   SX_CHECK (structPtr);
   update (structPtr->cell);
}

void SxGBasis::update (const SxCell &cell)
{
   SxVector3<Int> i;
   SxMesh3D mesh = fft3d(0).mesh;
   SxVector3<TPrecG> g;
   PrecEnergy gAbs;
   int ig = 0;
   int i0, j0, k0;
   SxList<PrecG> g2List;
   SxList<PrecFFTIdx>      n123List;
   SxList<int >  iList, jList, kList;
   SxList<PrecG> gVecListX, gVecListY, gVecListZ;
   SxMatrix3<Double> bMat = TWO_PI * cell.inv;
   for (int idx = 0; idx < 3; idx++)  {
     if (gCut > 4.*sqr(PI * mesh(idx) / cell.basis(idx).norm ()))  { // factor 4 correct here??
        cout << "Mesh " << mesh << " too small for gCut = " << gCut
             << " with cell = " << cell << "!" << endl;
        SxMesh3D minMesh;
        for (int ib = 0; ib < 3; ++ib)
           minMesh(ib) = (int)ceil(sqrt(gCut/4.)/PI*cell.basis(ib).norm ());
        cout << "Minimum size (cutoff = " << gCut << " Ry) is "
             << minMesh << endl;
        //SX_EXIT;

     }
   }


   // use only upper half space for 'G = -G' symmetry
   i0 = (useTimeRevSym) ? 0 : -mesh(0)/2;
   for (i(0) = i0; i(0) <= mesh(0)/2; i(0)++)  {
      j0 = (useTimeRevSym)
         ? ( (!i(0))  ?  0  :  -mesh(1)/2 )
         : -mesh(1)/2;
      for (i(1) = j0; i(1) <= mesh(1)/2; i(1)++)  {
         k0 = (useTimeRevSym)
            ? ( (!i(0) && !i(1))  ?  1  :  -mesh(2)/2 )
            : -mesh(2)/2;
         for (i(2) = k0; i(2) <= mesh(2)/2; i(2)++)  {

            g = (bMat.transpose() ^ i) + dG;

            gAbs = g ^ g;
            if (gAbs < gCut)  {
               g2List.append (gAbs);
               gVecListX.append (g(0));
               gVecListY.append (g(1));
               gVecListZ.append (g(2));
               // setup index array for quasi 1-dimensional fft working mesh
               n123List.append ((int)mesh.getMeshIdx (i, SxMesh3D::Origin));
               if (useTimeRevSym)  {
                  iList.append (i(0));
                  jList.append (i(1));
                  kList.append (i(2));
               }
               ig++;
            }
         }
      }
   }

   ng = ig;
   if (ng == 0)  {
      cout << SX_SEPARATOR;
      cout << "| FATAL ERROR:" << endl;
      cout << SX_SEPARATOR;
      cout << "  G+k basis does not have any plane waves!" << endl;
      cout << "  cutoff is " << gCut << " Ry" << endl;
      cout << "  k = " << dG << "; |k|^2 = " << dG.normSqr () << endl;
      cout << "  volume of reciprocal cell is 4pi/3|G|^3 with |G|^2="
           << pow(SxCell(bMat).volume / FOUR_PI * 3., 2./3.)
           << endl;
      cout << "  Possible mistakes: " << endl
           << "  * too small cutoff" << endl
           << "  * too small cell" << endl
           << "  * k-point far outside first Brillouin zone" << endl;
      cout << "  Please check cell volume, cutoff, and/or k-points" << endl;
      SX_QUIT;
   }

   n123.resize (1);
   g2      = SxDiracVec<TPrecG>      (g2List);
   n123(0) = SxVector<TPrecFFTIdx> (n123List);
   gVec.resize (ng*3); gVec.reshape (ng,3);
   gVec.colRef(0) <<= SxDiracVec<TPrecG>  (gVecListX);
   gVec.colRef(1) <<= SxDiracVec<TPrecG>  (gVecListY);
   gVec.colRef(2) <<= SxDiracVec<TPrecG>  (gVecListZ);

   // sort helper arrays according to g2
   SxArray<ssize_t> sortIdx = g2.getSortIdx();
   g2.sortByIdx (sortIdx);
   gVec.colRef(0).sortByIdx (sortIdx);
   gVec.colRef(1).sortByIdx (sortIdx);
   gVec.colRef(2).sortByIdx (sortIdx);
   n123(0).sortByIdx (sortIdx);


   // --- expand remaining half space
   if (useTimeRevSym)  {
      int fg, fullSx = 2*ng + 1;                   // ... + 1:  G=0 component
      SxDiracVec<Int>  iVec (iList), jVec (jList), kVec (kList);
      iVec.sortByIdx(sortIdx); jVec.sortByIdx(sortIdx); kVec.sortByIdx(sortIdx);
      g2.resize  (fullSx, true);
      // --- resize gVec
      SxDiracMat<TPrecG> gVecTmp (fullSx,3);
      gVecTmp.reshape (fullSx, 3);
      gVecTmp.colRef(0,0,ng-1) <<= gVec.colRef(0);
      gVecTmp.colRef(1,0,ng-1) <<= gVec.colRef(1);
      gVecTmp.colRef(2,0,ng-1) <<= gVec.colRef(2);
      gVec = gVecTmp;
      gVec.setBasis (this);

      n123(0).resize    (fullSx, true);
      for (ig=ng-1, fg=fullSx-1;
           ig >= 0;
           ig--, fg--)
      {
         // --- (+G) component
         g2(fg)      = g2(ig);
         gVec(fg,0) = gVec(ig,0);
         gVec(fg,1) = gVec(ig,1);
         gVec(fg,2) = gVec(ig,2);

         n123(0)(fg) = (int)mesh.getMeshIdx (iVec(ig), jVec(ig), kVec(ig),
                                             SxMesh3D::Origin);
         // --- (-G) component
         fg--;
         g2(fg)      = g2(ig);
         gVec(fg,0) = -gVec(ig,0);
         gVec(fg,1) = -gVec(ig,1);
         gVec(fg,2) = -gVec(ig,2);

         n123(0)(fg) = (int)mesh.getMeshIdx (-iVec(ig),-jVec(ig),-kVec(ig),
                                             SxMesh3D::Origin);
      }
      // --- G=0 component
      g2(0)     = 0.;
      gVec(0,0) = 0.; gVec(0,1) = 0.; gVec(0,2) = 0.;
      n123(0)(fg)  = (int)mesh.getMeshIdx (0, 0, 0, SxMesh3D::Positive);

      ng = fullSx;
   }

   UPDATE_MEMORY (structureFactors);
   UPDATE_MEMORY (phaseFactors);
}


void SxGBasis::changeTau (const SxAtomicStructure &tauList)
{
   structPtr = &tauList;
   int is, ia;
   int nSpecies = tauList.getNSpecies ();
   // find non-empty species
   for (is = 0; tauList.getNAtoms (is) == 0; is++) { /* empty */ }

   if (phaseFactors.getSize () != nSpecies
       || phaseFactors(is)(0).getSize () == 0)  {
      // --- allocate memory for both phase and structure factors
      phaseFactors.resize (nSpecies);
      structureFactors.resize (nSpecies);
      for (is=0; is < nSpecies; is++)  {
         structureFactors(is).resize (ng);
         phaseFactors(is).resize (tauList.getNAtoms(is));
         if (memMode == SaveTime)  {
            for (ia=0; ia < tauList.getNAtoms(is); ia++)  {
               phaseFactors(is)(ia).resize (ng);
               phaseFactors(is)(ia).handle->auxData.is = is;
               phaseFactors(is)(ia).handle->auxData.ia = ia;
            }
         }
      }
   }
   SX_CHECK (structureFactors.getSize () == nSpecies,
             structureFactors.getSize (), nSpecies);

   // --- compute phase and structure factors
   SxVector3<TPrecTauR> tau;
   int ig;
   PrecPhase  phase;
   PrecCoeffG c;
   SxDiracVec<TPrecCoeffG>::Iterator phFacIt, strFacIt;
   for (is=0; is < nSpecies; is++)  {
      structureFactors(is).set (0.);
      for (ia=0; ia < tauList.getNAtoms(is); ia++)  {
         tau = tauList(is,ia);
         strFacIt = structureFactors(is).begin ();
         if (memMode == SaveTime) {
            phFacIt = phaseFactors(is)(ia).begin ();
            for (ig=0; ig < ng; ++ig, ++phFacIt, ++strFacIt)  {
               phase = -(getG(ig) ^ tau);
               c     = PrecCoeffG (cos((double)phase),
                                   sin((double)phase));
               // phaseFactors(is)(ia)(ig)  = c;
               *phFacIt = c;
               // structureFactors(is)(ig) += c;
               *strFacIt += c;
            }
         } else {
            // delete old partial phase factors
            bool hasPhase = (phaseFactors(is)(ia).getSize () > 0);
            phaseFactors(is)(ia) = SxDiracVec<TPrecCoeffG> ();
            // --- compute structure factors
            structureFactors(is) += getPhaseFactors(is,ia);
            if (!hasPhase)
               phaseFactors(is)(ia).resize (0);
         }
      }
   }
   UPDATE_MEMORY (phaseFactors);
   UPDATE_MEMORY (structureFactors);
}

SxDiracVec<TPrecPhase> SxGBasis::getPhaseFactors(int is, int ia) const
{
   SX_CHECK (is >= 0 && is < phaseFactors.getSize (),
             is, phaseFactors.getSize () );
   SX_CHECK (ia >= 0 && ia < phaseFactors(is).getSize (),
             ia, phaseFactors(is).getSize () );
   if (memMode == SaveTime)  {
      return phaseFactors(is)(ia);
   }

   // --- calculate phase factors

   setupPhase1D (is, ia);

   SxDiracVec<TPrecPhase> res = composePhase (phaseFactors(is)(ia));
   res.handle->auxData.is = is;
   res.handle->auxData.ia = ia;
   return res;
}

SxComplex16 SxGBasis::getPhaseFactors(int is, int ia, ssize_t ig) const
{
   SX_CHECK (is >= 0 && is < phaseFactors.getSize (),
             is, phaseFactors.getSize () );
   SX_CHECK (ia >= 0 && ia < phaseFactors(is).getSize (),
             ia, phaseFactors(is).getSize () );

   if (memMode == SaveTime) return phaseFactors(is)(ia)(ig);

   // --- calculate 1D phase factors if necessary
   setupPhase1D (is, ia);
   SxDiracVec<TPrecPhase> &phFac = phaseFactors(is)(ia);
   return  phFac(packedGrel(3*ig))
         * phFac(packedGrel(3*ig+1))
         * phFac(packedGrel(3*ig+2));
}

SxDiracVec<Complex16>
SxGBasis::composePhase (const SxDiracVec<Complex16> &phFac) const
{
   SX_CHECK (phFac.getSize () == fft3d(0).mesh.sum (),
             phFac.getSize (), fft3d(0).mesh.sum ());
   SxDiracVec<Complex16> res(ng);

   if (packedGrel.getSize () != 3*ng) setupPackedG ();
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
   for (int ig = 0; ig < ng; ++ig)  {
      res(ig) = phFac(packedGrel(3*ig))
              * phFac(packedGrel(3*ig+1))
              * phFac(packedGrel(3*ig+2));
   }
#else
   typename SxDiracVec<Complex16>::Iterator resIt = res.begin ();
   ssize_t offset = 0;
   for (int ig = ng; ig; --ig, offset+=3)  {
      *resIt++ = phFac(packedGrel(offset  ))
               * phFac(packedGrel(offset+1))
               * phFac(packedGrel(offset+2));
   }
#endif
   res.setBasis (this);
   return res;
}

void SxGBasis::applyComposedPhase (const SxDiracVec<Complex16> &phFac,
                                  SxDiracVec<Complex16> *resPtr) const
{
   SX_CHECK (resPtr);
   SxDiracVec<Complex16> &res = *resPtr;
   SX_CHECK (res.getBasisPtr () == this);
   SX_CHECK (res.getSize () == ng, res.getSize (), ng);
   SX_CHECK (packedGrel.getSize () == 3*ng, packedGrel.getSize (), ng);
   SX_CHECK (phFac.getSize () == fft3d(0).mesh.sum (),
             phFac.getSize (), fft3d(0).mesh.sum ());

   SxDiracVec<Complex16>::Iterator resIt= res.begin ();
   ssize_t offset = 0;
   for (int ig = ng; ig; --ig, offset+=3)  {
      *resIt++ *= phFac(packedGrel(offset  ))
                * phFac(packedGrel(offset+1))
                * phFac(packedGrel(offset+2));
   }
}

SxAtomicStructure 
SxGBasis::get1DPackedVecs () const
{
   SX_CHECK (structPtr);
   SxCell recCell = structPtr->cell.getReciprocalCell ();
   SxMesh3D &mesh = fft3d(0).mesh;
   SxAtomicStructure res(recCell, mesh.sum ());
   int offset = 0;
   for (int dim = 0; dim < 3; ++dim)  {
      int n = mesh(dim);
      Coord gRel (0.,0.,0.);
      for (int i = 0; i < n; ++i, ++offset)  {
         gRel(dim) = (i+i > n) ? (i - n) : i;
         res.ref(offset) = dG + recCell.relToCar (gRel);
      }
   }

   if (packedGrel.getSize () != 3*ng)  {
      // --- maps ig vectors to 1D phases
      packedGrel.resize (3*ng);
      SxVector3<Int> Grel;
      SxVector<TPrecFFTIdx>::Iterator n123It = n123(0).begin ();
      // x,y,z are just names. We work in relative coordinates here.
      ssize_t nx = fft3d(0).mesh(0), ny = fft3d(0).mesh(1);
      for (int ig = 0; ig < ng; ++ig, ++n123It)  {
         Grel = fft3d(0).mesh.getMeshVec(*n123It, SxMesh3D::Positive);
         packedGrel(3*ig)   = Grel(0);
         packedGrel(3*ig+1) = Grel(1) + nx;
         packedGrel(3*ig+2) = Grel(2) + nx + ny;
      }
   }

   return res;
}

void SxGBasis::setupPhase1D (int is, int ia) const
{
   SX_CHECK (memMode == SaveMemory);
   SX_CHECK (structPtr);
   const SxAtomicStructure &str = *structPtr;
   SX_CHECK (is >= 0 && is < str.getNSpecies (), is, str.getNSpecies ());
   SX_CHECK (ia >= 0 && ia < str.getNAtoms (is), ia, str.getNAtoms (is));

   SxDiracVec<TPrecPhase> &phFac = phaseFactors(is)(ia);

   if (phFac.getSize () != fft3d(0).mesh.sum ())
      phFac = setupPhase1D(str.constRef(is,ia));
   if (packedGrel.getSize () != 3*ng) setupPackedG ();
}

void SxGBasis::setupPackedG () const
{
   if (packedGrel.getSize () != 3*ng)  {
      // --- maps ig vectors to 1D phases
      packedGrel.resize (3*ng);
      SxVector3<Int> Grel;
      SxVector<TPrecFFTIdx>::Iterator n123It = n123(0).begin ();
      // x,y,z are just names. We work in relative coordinates here.
      ssize_t nx = fft3d(0).mesh(0), ny = fft3d(0).mesh(1);
      for (int ig = 0; ig < ng; ++ig, ++n123It)  {
         Grel = fft3d(0).mesh.getMeshVec(*n123It, SxMesh3D::Positive);
         packedGrel(3*ig)   = Grel(0);
         packedGrel(3*ig+1) = Grel(1) + nx;
         packedGrel(3*ig+2) = Grel(2) + nx + ny;
      }
   }
}

void SxGBasis::registerUnknown(const SxBasis &basis) const
{
   const SxRBasis *rPtr = dynamic_cast<const SxRBasis*>(&basis);
   if (rPtr) registerBasis (*rPtr); // special case
   else      registerBasis (basis); // generic case
}



namespace Timer {
   enum registerBasisTimer {
      regbas1, regbas2, regbas3, regbas4, regbas5, regbas6, regbas7, regbas8
   };
}
SX_REGISTER_TIMERS (Timer::registerBasisTimer)
{
   using namespace Timer;
   regTimer (regbas1, "regbas1");
   regTimer (regbas2, "regbas2");
   regTimer (regbas3, "regbas3");
   regTimer (regbas4, "regbas4");
   regTimer (regbas5, "regbas5");
   regTimer (regbas6, "regbas6");
   regTimer (regbas7, "regbas7");
   regTimer (regbas8, "regbas8");
}

void SxGBasis::registerBasis (const SxRBasis &rBasis) const
{
   if (!SxBasis::registerBasis (rBasis)) return;
   int idx = getBasisId (&rBasis);
   if (idx == 0)  {
      // --- first real-space basis
      if ((rBasis.fft3d.mesh - fft3d(0).mesh).absSqr ().sum () != 0)  {

         SX_START_TIMER (Timer::regbas1);
         SxFFT3d newFFT (SxFFT3d::Forward,
                       rBasis.fft3d.mesh(0),
                       rBasis.fft3d.mesh(1),
                       rBasis.fft3d.mesh(2),
                       rBasis.cell.volume, true);
         SX_STOP_TIMER (Timer::regbas1);

         SX_START_TIMER(Timer::regbas2);
         const_cast<SxGBasis*>(this)->replaceMesh (newFFT);
         SX_STOP_TIMER (Timer::regbas2);

      } else {
         return;
      }
   } else {

      SX_START_TIMER(Timer::regbas3);
      SxFFT3d newFFT (SxFFT3d::Forward,
                    rBasis.fft3d.mesh(0),
                    rBasis.fft3d.mesh(1),
                    rBasis.fft3d.mesh(2),
                    rBasis.cell.volume, true);
      SX_STOP_TIMER(Timer::regbas3);

      // --- add one more fft
      SX_CHECK (idx == fft3d.getSize (), idx, fft3d.getSize ());
      SX_CHECK (idx == n123.getSize (), idx, n123.getSize ());

      {
         SX_CLOCK(Timer::regbas4);

         fft3d.resize (idx + 1, true);
         n123.resize  (idx + 1, true);
         fft3d(idx) = newFFT;
         n123(idx).resize (ng);
      }

      // --- set up new n123
      const SxFFT3d &origFFT = fft3d(0);
      SxVector<TPrecFFTIdx>::Iterator n123OrigIt = n123(0).begin (),
                                        n123NewIt = n123(idx).begin ();
      SxVector3<Int> vec;

      {
         SX_CLOCK(Timer::regbas5);

         for (int ig = 0; ig < ng; ig++, ++n123OrigIt, ++n123NewIt)  {
            vec = origFFT.mesh.getMeshVec(*n123OrigIt, SxMesh3D::Origin);
            if (   (2*abs(vec(0)) > newFFT.mesh(0))
                || (2*abs(vec(1)) > newFFT.mesh(1))
                || (2*abs(vec(2)) > newFFT.mesh(2)))
            {
               cout << "new FFT mesh too small for some G vectors." << endl;
               SX_EXIT;
            }
            *n123NewIt = (int)newFFT.mesh.getMeshIdx(vec, SxMesh3D::Origin);
         }
      }
   }

   {
      SX_CLOCK(Timer::regbas6);

      fft3d(idx).autoscale = false;
      SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
      if (memMode == SaveTime)
         for (int i = 0; i < fft3d.getSize (); ++i) fft3d(i).clean ();
      else
         for (int i = 0; i < fft3d.getSize (); ++i) fft3d(i).destroyArrays ();
   }

}

void SxGBasis::deregister (const SxBasis* basis) const
{
   const SxRBasis *rPtr = dynamic_cast<const SxRBasis*>(basis);
   int nKnownR = int(fft3d.getSize ());
   SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
   if (rPtr && nKnownR > 1) {
      for (int i = getBasisId (rPtr) + 1; i < nKnownR; ++i)  {
         fft3d(i-1) = fft3d(i);
         n123(i-1) = n123(i);
      }
      fft3d.resize(nKnownR - 1, true);
      if (memMode == SaveTime)
         for (int i = 0; i < fft3d.getSize (); ++i) fft3d(i).clean ();
      else
         for (int i = 0; i < fft3d.getSize (); ++i) fft3d(i).destroyArrays ();
      n123.resize (nKnownR - 1, true);
   }
   if (rPtr == rBasisPtr) rBasisPtr = NULL;
   deregisterBasis (basis);
}


void SxGBasis::setupMixedFFT (const SxRBasis &rBasis)
{
   registerBasis (rBasis);
   if (mixedFFT && mixedFFT->realMesh == rBasis.fft3d.mesh)
      return;
   int idx = getBasisId (&rBasis);
   int nonZero = 0;
   SxVector<TPrecFFTIdx>::Iterator n123It = n123(idx).begin ();
   const SxMesh3D &mesh = rBasis.fft3d.mesh;
   for ( ; n123It != n123(idx).end (); n123It++)  {
      SxVector3<Int> vec = mesh.getMeshVec (*n123It, SxMesh3D::Origin);
      int ix = (vec(2) < 0) ? (-vec(2)) : (vec(2) + 1);
      if (ix > nonZero) nonZero = ix;
   }
   mixedFFT = SxPtr<SxFFT2d1d>::create (mesh(0), mesh(1), mesh(2), nonZero);
   n231 = mixedFFT->getN231 (n123(idx));
}

PsiGI SxGBasis::convolute (const PsiGI &allPsi,
                           const SxDiracVec<Double> &V) const
{
   const SxRBasis &R = V.getBasis<SxRBasis> ();
   SX_CHECK (V.getSize () == R.fft3d.mesh.getSize (),
             V.getSize (), R.fft3d.mesh.getSize ());

   // --- check that mixedFFT is for this R basis
   if (mixedFFT)  {
      if (R.fft3d.mesh != mixedFFT->realMesh)  {
         cout << "Direct convolution: mesh mismatch!" << endl;
         cout << "Potential has " << R.fft3d.mesh << endl;
         cout << "2+1D FFT has " << mixedFFT->realMesh << endl;
         SX_EXIT;
      }
   } else {
      const_cast<SxGBasis*>(this)->setupMixedFFT (R);
   }
   SX_CHECK (n231.getSize () == ng, n231.getSize (), ng);

   SX_CLOCK (Timer::Convolute);
   ssize_t nStates = allPsi.nCols ();
   double scale = 1./double(R.fft3d.mesh.getSize ());
   PsiGI res;
   res.reformat (ng, nStates);

   SxComplex16* mesh = mixedFFT->getMesh ();
   for (ssize_t iState = 0; iState < nStates; ++iState)  {
      // clean the mesh
      mixedFFT->clean ();
      PsiG psi = allPsi.colRef (iState);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      // put data on mesh
      for (ssize_t ig = 0; ig < ng; ++ig)
         mesh[n231(ig)] = scale * psi(ig);

      // do the direct convolution
      {
         SX_CLOCK (Timer::ConvoluteFFT);
         mixedFFT->convolute (V.elements);
      }

      PsiG Vpsi = res.colRef (iState);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
      // fetch data from mesh
      for (ssize_t ig = 0; ig < n231.getSize (); ++ig)
         Vpsi(ig) = mesh [n231(ig)];

   }
   mixedFFT->freeMesh ();
   res.setBasis (this);
   return res;
}

void SxGBasis::addToRho (const PsiGI &allPsi,
                         const SxDiracVec<Double> &focc,
                         SxDiracVec<Double> *rho) const
{
   SX_CHECK (rho);
   const SxRBasis &R = rho->getBasis<SxRBasis> ();
   SX_CHECK (rho->getSize () == R.fft3d.mesh.getSize (),
             rho->getSize (), R.fft3d.mesh.getSize ());

   // --- check that mixedFFT is for this R basis
   if (mixedFFT)  {
      if (R.fft3d.mesh != mixedFFT->realMesh)  {
         cout << "Direct rho: mesh mismatch!" << endl;
         cout << "Potential has " << R.fft3d.mesh << endl;
         cout << "2+1D FFT has " << mixedFFT->realMesh << endl;
         SX_EXIT;
      }
   } else {
      const_cast<SxGBasis*>(this)->setupMixedFFT (R);
   }
   SX_CHECK (n231.getSize () == ng, n231.getSize (), ng);

   SX_CLOCK (Timer::AddToRho);
   ssize_t nStates = allPsi.nCols ();
   int idx = getBasisId (&R);
   double scale = fft3d(idx).scaleForFFT;

   SxComplex16* mesh = mixedFFT->getMesh ();
   for (ssize_t iState = 0; iState < nStates; ++iState)  {
      if (fabs(focc(iState)) <= 1e-12) continue;
      // clean the mesh
      mixedFFT->clean ();
      PsiG psi = allPsi.colRef (iState);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      // put data on mesh
      for (ssize_t ig = 0; ig < ng; ++ig)
         mesh[n231(ig)] = psi(ig);

      // do the direct convolution
      {
         SX_CLOCK (Timer::AddToRhoFFT);
         mixedFFT->addToRho ((scale*scale) * focc(iState), rho->elements);
      }
   }
   mixedFFT->freeMesh ();
}


void SxGBasis::cleanPhaseFactors () const
{
   SX_CHECK (structPtr);
   const SxAtomicStructure &str = *structPtr;
   int is, nSpecies = str.getNSpecies ();
   int ia, nAtoms;
   for (is = 0; is < nSpecies; ++is)  {
      nAtoms = str.getNAtoms(is);
      for (ia = 0; ia < nAtoms; ++ia)
         phaseFactors(is)(ia).resize (0);
   }
}


SxDiracVec<Complex16> SxGBasis::setupPhase1D (const Coord &tau) const
{
   SX_CHECK (structPtr);
   const SxAtomicStructure &str = *structPtr;
   SxDiracVec<Complex16> phase1(fft3d(0).mesh.sum ());
   // relative coordinates, 2pi factor
   Coord tauRel = str.cell.carToRel (tau);
   Coord dGrel = str.cell.getReciprocalCell ().carToRel (dG);
   // --- setup 1D phase
   for (int ip = 0, idir = 0; idir < 3; idir++)  {
      double r = TWO_PI * tauRel(idir);
      double dg = dGrel(idir);
      int n = fft3d(0).mesh(idir);
      for (int i = 0; i < n; i++, ip++)  {
         double g = dg + ((2 * i > n) ? (i - n) : i);
         SxComplex16 &res = phase1(ip);
         sincos (-g * r, &res.im, &res.re);
      }
   }
   return phase1;
}

SxArray<SxDiracMat<TPrecPhase> >
SxGBasis::getPhaseFactors (const SxAtomicStructure &structure) const
{
   int nSpecies = structure.getNSpecies ();
   SX_CHECK (nSpecies > 0, nSpecies);
   SX_CHECK (structure.getNAtoms () > 0, structure.getNAtoms ());

   SxArray<SxDiracMat<TPrecPhase> > result(nSpecies);

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      int nAtoms = structure.getNAtoms (iSpecies);
      result(iSpecies).reformat (ng, nAtoms);
      for (int iAtom = 0; iAtom < nAtoms; ++iAtom)
         result(iSpecies).colRef (iAtom)
            <<= getPhaseFactors(structure(iSpecies, iAtom));
      VALIDATE_VECTOR (result(iSpecies));
   }
   return result;
}

void SxGBasis::replaceMesh (const SxFFT3d &newFFT)
{
   SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
   SX_CHECK (newFFT.meshSize > 0, newFFT.meshSize);
   SX_CHECK (newFFT.dir == SxFFT3d::Forward || newFFT.dir == SxFFT3d::Both);
   SxVector3<Int> vec;
   for (int ig = 0; ig < ng; ig++)  {
      vec = fft3d(0).mesh.getMeshVec(n123(0)(ig), SxMesh3D::Origin);
      if (   (2*abs(vec(0)) > newFFT.mesh(0))
          || (2*abs(vec(1)) > newFFT.mesh(1))
          || (2*abs(vec(2)) > newFFT.mesh(2)))
      {
         cout << "Invalid change of FFT meshsize: old FFT is too large.";
         cout << endl;
         SX_EXIT;
      }
      n123(0)(ig) = (int)newFFT.mesh.getMeshIdx(vec, SxMesh3D::Origin);
   } // ig
   fft3d(0) = newFFT;
   fft3d(0).autoscale = false;
   if (memMode == SaveTime)
      fft3d(0).clean ();
   else
      fft3d(0).destroyArrays ();
   // n123Inv is no longer up-to-date
   n123inv.resize (0);
   if (memMode == SaveMemory)  {
      packedGrel.resize (0);
      if (structPtr) cleanPhaseFactors ();
   }
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int SxGBasis::getMaxNG ()
{
   SX_EXIT;
   return -1;
//   double maxNG = 8./(3. ) * geom->omega
//                * (double)pow ((double)control->eCut, (double)1.5) * 1.2;
//   return (int)maxNG;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxVector3<TPrecG> SxGBasis::getG (int ig) const
{
   return SxVector3<TPrecG> (gVec(ig,0), gVec(ig,1), gVec(ig,2));
}

SxVector3<TPrecG> SxGBasis::getK () const
{
   SxVector<TPrecFFTIdx> &idx = n123(0);
   // note: (0,0,0) should be one of the first elements (sorted by |G+k|)
   for (int ig = 0; ig < ng; ++ig)
      if (idx(ig) == 0) return getG(ig);
   SX_EXIT;
   return SxVector3<TPrecG> ();
}

SxDiracVec<SxGBasis::TBasisType>
SxGBasis::identity (const SxGBasis *basis,
                    const SxDiracVec<TBasisType> &vec) const
{
   // --- verify that this is really the identity projector, i.e.
   //     that the provided basis is the same as the vector's basis
   SX_CHECK (vec.getBasisPtr() == this);
   SX_CHECK (basis == this);
   return vec;
}

SxDiracVec<SxGBasis::TBasisType>
SxGBasis::toRealSpace (const SxRBasis *rBasis_,
                       const SxDiracVec<TGBasisType> &psiG) const
{
   SX_CLOCK (Timer::GtoR);
   // check that wavefunctions match to this basis
   SX_CHECK (psiG.getSize() == ng*nComp, psiG.getSize(), ng*nComp);
   SX_CHECK (rBasis_->nComp == nComp, rBasis_->nComp, nComp);

   SX_CHECK (rBasis_);
   // "handshake"
   registerBasis(*rBasis_); // does nothing if rBasis_ is already registered

   // --- set up fft mesh
   int ig, ngk   = (int)psiG.getSize ()/nComp;
   int idx = getBasisId (rBasis_);
   SX_CHECK (idx >= 0 && idx < fft3d.getSize (), idx, fft3d.getSize ());

   PrecFFTIdx *n123Ptr;
   SxFFT3d &fft = fft3d(idx);
   SX_CHECK (fft.mesh(0) == rBasis_->fft3d.mesh(0),
             fft.mesh(0), rBasis_->fft3d.mesh(0));
   SX_CHECK (fft.mesh(1) == rBasis_->fft3d.mesh(1),
             fft.mesh(1), rBasis_->fft3d.mesh(1));
   SX_CHECK (fft.mesh(2) == rBasis_->fft3d.mesh(2),
             fft.mesh(2), rBasis_->fft3d.mesh(2));

   SxDiracVec<Complex16> psiOut (fft.meshSize*nComp);
   psiOut.handle->auxData = psiG.handle->auxData;

   if (hasMixedFFT(rBasis_->fft3d.mesh))  {
      SxComplex16* mesh = mixedFFT->getMesh ();

      for (int c = 0; c < nComp; c++)  {
         mixedFFT->clean ();
         SX_FFT_REAL scale = fft.scaleForFFT;
#ifdef USE_OPENMP
#pragma omp parallel for if (ngk > sxChunkSize)
#endif
         // put data on mesh
         for (ssize_t ig = 0; ig < ngk; ++ig)
            mesh[n231(ig)] = scale * psiG(ig + c * ngk);

         // do the FT
         SX_CLOCK (Timer::GtoR_FFT);
         mixedFFT->fftForward (psiOut.elements + c * fft.meshSize);
      }
      mixedFFT->freeMesh ();
      psiOut.setBasis (rBasis_);
      VALIDATE_VECTOR (psiOut);
      return psiOut;
   }

   bool needWorkArray = (fft.inArray == NULL);
   if (needWorkArray) fft.createArrays(SxFFT::InArrayZero);

   int size = fft.meshSize;
   SxDiracVec<TGBasisType> psiGC;

   // --- loop over components
   for (int c = 0; c < nComp; c++)  {
      psiGC = SxDiracVec<TGBasisType> (); // release previous psiGC reference!
      psiGC = psiG(SxIdx(c*ngk, (c+1)*ngk-1));

      const SxDiracVec<SX_T_FFT_COMPLEX> psiIn  (psiGC);
      SX_FFT_COMPLEX *psiGPtr = psiIn.elements;
      n123Ptr = n123(idx).elements;
      SX_FFT_REAL s = fft.scaleForFFT;
      fft.clean();
      // make sure that fft is cleaned: mesh/2 should always be outside cutoff
      SX_CHECK ((double&)fft.inArray[fft.mesh.getMeshIdx(fft.mesh/2,
                SxMesh3D::Unknown)] == 0.);
      //   fft.clean ();

#  ifdef USE_OPENMP
      if (ngk > sxChunkSize)  {
#        pragma omp parallel for shared(n123Ptr,psiGPtr,ngk) \
         private(ig) schedule(static)
         for (ig=0; ig < ngk; ig++)
            fft.setElement (n123Ptr[ig], s * psiGPtr[ig]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ig=0; ig < ngk; ig++)
            fft.setElement (*n123Ptr++, s * *psiGPtr++);
      }

      {  SX_CLOCK (Timer::GtoR_FFT);
      fft.fftForward (fft.meshSize, fft.inArray,
            psiOut(SxIdx(c * size, (c+1)*size -1)).elements);
      }

   }
   if (needWorkArray) fft.destroyArrays ();

   // --- Remember! The following statement performs the copy only in case
   //     types of psiOut and psiR. See description of reference counting
   //     in SxVec<T>
   SxDiracVec<TBasisType> psiR (psiOut);
   psiR.setBasis (rBasis_);

   VALIDATE_VECTOR (psiR);

   // cout << "--> toRealSpace called!  nComp = " << nComp << endl;

   return psiR;
}

SxDiracVec<SxGBasis::TBasisType>
SxGBasis::toRadialSpace ( const SxRadBasis *radPtr,
                          const SxDiracVec<TGBasisType> &psiG) const
{
   SX_CLOCK(Timer::GtoRad);
   // check that wavefunctions math to this basis
   SX_CHECK (psiG.getSize() == ng*nComp, psiG.getSize(), ng*nComp);

   SX_CHECK (radPtr);
   SX_CHECK (structPtr);
   SX_CHECK (structPtr->cell.volume > 0.);
   double rNorm = FOUR_PI / sqrt(structPtr->cell.volume);
   const SxRadBasis &radBasis = *radPtr;
   // Save Quantumnumbers
   int is = psiG.handle->auxData.is;
   SX_CHECK (is >= 0, is);
   int ia = psiG.handle->auxData.ia;
   int n = psiG.handle->auxData.n;
   int l = psiG.handle->auxData.l; 
   SX_CHECK (l >= 0, l);
   int m = psiG.handle->auxData.m;

   SxDiracVec<TGBasisType> work = psiG.getCopy ();
   
   // Transform for a special ia ?
   if (ia != -1)   {
      SxDiracVec<Complex16> phase = getPhaseFactors(is,ia).conj();
      work *= phase;
   }
   // Transform for a special m ? 
   if (m != NONE_M)   {
      SX_CHECK (abs(m) <= l,abs(m),l);
      SxDiracVec<Double> YlmVec = getYlm(l,m);
      work *= rNorm * SxYlm::getYlmNormFactor(l,m) * YlmVec;
   }

   // Get jl
   size_t gDim = g2.getSize ();
   size_t rDim = radBasis.radFunc(is).getSize();
   SxDiracVec<Double> result(rDim);
   result.set(0.0);
   SxDiracVec<Double> rAbs = radBasis.radFunc(is);
   SxDiracVec<Double> jl(gDim);
   for(size_t ir = 0; ir < rDim; ir++)   {
      if ( rAbs(ir) < radBasis.cutoff )   {
         jl.set(0.0);
         double gLast = 0.0;
         SX_CLOCK (Timer::GBessel);
         for (size_t ig = 0; ig < gDim; ig++)   {
            if ((ig == 0) ||(sqrt(g2(ig)) != gLast))   {
               jl(ig) = SxYlm::jsb(l,rAbs(ir)*sqrt(g2(ig)));
               gLast = sqrt(g2(ig));
            }
            else   {
               jl(ig) = jl(ig-1);
            }
         }
         result(ir) = (work * jl).sum();
      }
   }
   result.handle->auxData.is = is; 
   result.handle->auxData.ia = ia; 
   result.handle->auxData.n  = n; 
   result.handle->auxData.l  = l; 
   result.handle->auxData.m  = m;
   result.setBasis(&radBasis);  
   
   // cout << "--> toRadialSpace called!" << endl;

   return result;
}

SxDiracVec<TAOBasisType>
SxGBasis::toAO (const SxAOBasis *aoBasis_,
                const SxDiracVec<TGBasisType> &psiG) const
{
   SX_CHECK (aoBasis_);
   SX_CHECK (psiG.getBasisPtr () == this);
   return aoBasis_->fromPWBasis (psiG);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Real8 SxGBasis::laplacian (const void *in) const
{
   SX_CHECK (in);
   const SxDiracVec<TPrecCoeffG> &coeff = *(const SxDiracVec<TPrecCoeffG> *)in;

   SX_CHECK (g2.getSize() * nComp == coeff.getSize(),
             g2.getSize(), nComp, coeff.getSize());
   SX_CHECK (coeff.getBasisPtr() == this);

   Real8 res = 0.;
   for (int iComp = 0; iComp < nComp; ++iComp)  {
      ssize_t compStart = iComp * ng;
#pragma omp parallel for if (ng > sxChunkSize) reduction(+:res)
      for (ssize_t ig = 0; ig < ng; ++ig)  {
         res += g2(ig) * coeff(ig + compStart).absSqr ();
      }
   }

   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxVector<TPrecCoeffR>
SxGBasis::GtoR111 (const SxDiracVec<TPrecCoeffG> &vec,
                   int nr, int shift, const SxCell &cell) const
{
   SX_CHECK (vec.getSize() > 0, vec.getSize());
   SX_CHECK (nr > 0, nr);
   SX_CHECK (shift >= 0, shift);

   int ir, ig;
   PrecPhase  phase;
   PrecCoeffG expph, sum;
   SxVector3<TPrecG> r;
   SxDiracVec<TPrecCoeffR> res(nr, 0.);
   PrecG a = exp( log (4. * cell.volume) / 3. );  // for FCC only!!!

   for (ir = 0; ir < nr; ir++)  {
      r = SxVector3<TPrecG> (1., 1., 1.) * (ir-shift) * a * 0.01; // Haeh??
      sum = 0.;
      for (ig = 0; ig < ng; ig++)  {
         phase  = (getG(ig) ^ r);
         expph  = PrecCoeffG (cos ((double)phase), sin ((double)phase));
         sum   += vec(ig) * (PrecCoeffG)expph;
      }
      res(ir) = sum;
   }
   return toVector(res) * HA2EV / cell.volume;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxGBasis::n123invSetup () const
{
   int ig;
   int fftSize = fft3d(0).mesh.product();

   if (n123inv.getSize () == fftSize) return;
   n123inv.resize (fftSize);
   n123inv = -1;

   for (ig = 0; ig < ng; ig++)  n123inv( n123(0)(ig) ) = ig;
}


//------------------------------------------------------------------------------
// scalar product
//------------------------------------------------------------------------------
SxComplex16 SxGBasis::scalarProduct (const SxGBasis &aBasis,
                                     const SxDiracVec<Complex16> &a,
                                     const SxDiracVec<Complex16> &b)
{
   SX_CHECK (a.getBasisPtr() == &aBasis);
   SX_CHECK (b.getBasisPtr() == this);

   int         igA, igB;
   SxComplex16 out (0.);

   for (igB = 0; igB < ng; igB++)  {
      igA = conveyIdx (igB, &aBasis);
      if (igA > -1)  out += a(igA).conj() * b(igB);
   }

   return out;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxDiracVec<Complex16>
SxGBasis::shiftK (const SxDiracVec<Complex16> &f,
                  const SxVector3<TPrecG>     &deltaGCart,
                  const SxArray<Coord>        &rMesh)
{
   SX_CHECK (fft3d(0).mesh.product() == rMesh.getSize(),
             fft3d(0).mesh.product(), rMesh.getSize());
   SX_CHECK      (rBasisPtr);
   SX_CHECK      (f.getBasisPtr() == this);

   int                    meshSize = fft3d(0).mesh.product();
   double                 phase;
   SxDiracVec<Complex16>  shiftFactor(meshSize);
   SxDiracVec<Complex16>  fShifted;

   for (int rIdx = 0; rIdx < meshSize; rIdx++)  {
      phase             = -(deltaGCart ^ rMesh(rIdx));
      shiftFactor(rIdx) = PrecCoeffR (cos(phase), sin(phase));
   }

   fShifted = shiftFactor * toRealSpace (rBasisPtr, f);

   return rBasisPtr->toGSpace (this, fShifted);
}

SxDiracVec<Complex16>
SxGBasis::shiftK (const SxDiracVec<Complex16> &f, RelVec &deltaGRel)
{
   SX_CHECK (f.getBasisPtr() == this);
   int                    idx;
   SxDiracVec<Complex16>  fShifted(ng);
   fShifted.set (0.);
   fShifted.setBasis (this);

   for (int ig = 0; ig < ng; ig++)  {
      idx = getIdxGSum (ig, deltaGRel);
      if (idx > -1)  fShifted(idx) = f(ig);
   }

   return fShifted;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
PrecEnergy SxGBasis::getGCut (PrecEnergy eCut)
{
   return 4. * eCut;
}


void SxGBasis::write (SxBinIO &/*io*/)
{
   SX_EXIT;
}


PsiG SxGBasis::mapToFFT (const PsiG &psiIn,
                         const SxVector<TPrecFFTIdx> &fftIdxIn,
                         const SxMesh3D        &mesh) const
{
   SX_CHECK (psiIn.getSize() == ng, psiIn.getSize(), ng);

   PsiG psiOut (ng);

   PsiG psi3d (mesh.getSize ());

   SxVector<TPrecFFTIdx>::Iterator n123It;

   n123It = fftIdxIn.begin ();
   for (int ig=0; ig < ng; ig++, ++n123It)  psi3d(*n123It) = psiIn(ig);

   n123It = n123(0).begin ();
   if (mesh == fft3d(0).mesh)  {
      for (int ig = 0; ig < ng; ig++, ++n123It)
         psiOut(ig) = psi3d(*n123It);
   } else  {
      const SxMesh3D oldMesh = fft3d(0).mesh;
      int idx;
      for (int ig = 0; ig < ng; ++ig, ++n123It)  {
         idx = (int)mesh.getMeshIdx (oldMesh.getMeshVec(*n123It,
                                     SxMesh3D::Origin),       SxMesh3D::Origin);
         psiOut(ig) = psi3d(idx);
      }
   }

   VALIDATE_VECTOR (psiOut);

   return psiOut;
}


void SxGBasis::registerMemoryObservers ()
{
   TRACK_MEMORY (gVec);
   TRACK_MEMORY (g2);
   TRACK_MEMORY (structureFactors);
   TRACK_MEMORY (phaseFactors);
   TRACK_MEMORY (packedGrel);
   TRACK_MEMORY (n123);
   TRACK_MEMORY (n123inv);


}


void SxGBasis::print () const
{
   sxprintf ("SxGBasis:\n");
   for (int i = 0; i < fft3d.getSize (); ++i)  {
      sxprintf ("   FFT mesh size:     %d x %d x %d\n",
              fft3d(i).mesh(0), fft3d(i).mesh(1), fft3d(i).mesh(2));
   }
   sxprintf ("   gCut:              %g\n", gCut);
   sxprintf ("   dG:                (%g,%g,%g)\n", dG(0), dG(1), dG(2));
   sxprintf ("   use E(G) = E(-G):  %d\n", useTimeRevSym);
   sxprintf ("   # G-vectors:       %d\n", (int)g2.getSize());
}

SxDiracVec<TPrecCoeffG>
SxGBasis::transferWaves (const SxDiracVec<TPrecCoeffG> &in) const
{
   // --- get in's G basis and dimensions
   SX_CHECK (in.getSize () > 0);
   const SxGBasis *oldBasis = dynamic_cast<const SxGBasis*>(in.getBasisPtr ());
   SX_CHECK (oldBasis);
   int ngOld = oldBasis->ng, nStates = (int)in.nCols ();
   SX_CHECK (ngOld == in.nRows (), ngOld, in.nRows ());
   if (nStates == 0 && in.getSize () == ngOld) nStates = 1;

   SX_CHECK (oldBasis->fft3d(0).meshSize == fft3d(0).meshSize,
             oldBasis->fft3d(0).meshSize, fft3d(0).meshSize);

   // --- prepare output
   SxDiracVec<TPrecCoeffG> res;
   res.reformat(ng, nStates);
   res.handle->auxData = in.handle->auxData;
   res.setBasis (this);

   // --- transfer by shuffling data on fft mesh
   PsiG psiMesh(fft3d(0).meshSize);
   psiMesh.set (0.);
   SxVector<TPrecFFTIdx>::Iterator n123It;
   PsiG::Iterator inIt = in.begin (), resIt = res.begin ();

   int ig, i;
   for (i = 0; i < nStates; ++i)  {
      n123It = oldBasis->n123(0).begin ();
      for (ig = 0; ig < ngOld; ++ig) psiMesh(*n123It++) = *inIt++;
      n123It = n123(0).begin ();
      for (ig = 0; ig < ng; ++ig) *resIt++ = psiMesh(*n123It++);
   }
   SX_CHECK (resIt == res.end ());
   SX_CHECK (inIt == in.end ());

   return res;

}

void SxGBasis::read (const SxBinIO &io, int ngIn, int offset,
                     const SxVector3<Int> &mesh)
{
   ng = ngIn;
   structPtr = NULL;
   n123.resize (1);
   // --- handle FFT
   if (mesh.product () != 0)  {
      fft3d.resize (1);
      fft3d(0).symmetric = true;
      fft3d(0).dir = SxFFT3d::Forward;
      fft3d(0).setMesh(mesh(0),mesh(1),mesh(2), SxCell(io).volume);
      SX_CHECK (memMode == SaveTime || memMode == SaveMemory);
      if (memMode == SaveTime)
         fft3d(0).clean ();
      else
         fft3d(0).destroyArrays ();
      fft3d(0).autoscale = false;
   } else {
      fft3d.resize (0);
   }
   SxDiracMat<Double> gVecMatrixRef;
   try {

      // read gVec
      gVec.reformat(ng,3);
      gVecMatrixRef = gVec; // ugly
      io.read ("gkVec", &gVecMatrixRef, ng, 3, offset);

      // read n123
      n123(0).resize (ng);
      io.read ("fftIdx", &n123(0), ng, offset);

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   // calculate g2
   g2.resize (ng);
   for (int ig = 0; ig < ng; ig++)
      g2(ig) = gVec.row(ig).absSqr ().sum ();
}

SxDiracVec<Complex16>
SxGBasis::toPartials (const SxPartialWaveBasis *pBasis,
                      const SxDiracVec<TGBasisType> &psi) const
{
   SX_CHECK (pBasis);
   SX_CHECK (pBasis->projectors);
   SX_CHECK (psi.handle);
   SxDiracVec<TPrecCoeffG> res = pBasis->projectors->fromPWBasis (psi);
   res.setBasis (pBasis);
   res.handle->auxData.ik = psi.handle->auxData.ik;
   return res;
}

SxDiracVec<Complex16> 
SxGBasis::toPAW (const SxPAWBasis *pawBasis,
                 const SxDiracVec<TGBasisType> &psi) const
{
   SX_CHECK (pawBasis);
   SX_CHECK (pawBasis->pBasis);
   SX_CHECK (pawBasis->pBasis->projectors);
   SX_CHECK (psi.handle);
   SxDiracVec<TPrecCoeffG> proj;
   proj = pawBasis->pBasis->projectors->fromPWBasis (psi);

   int nPsi = (int)psi.nCols ();
   if (nPsi == 0) nPsi = 1;
   SX_CHECK (psi.getSize () > 0);
   int N = int(pawBasis->getNElements ());
   SxDiracVec<Complex16> res;
   res.reformat (N, nPsi);
   res.handle->auxData = psi.handle->auxData;
   res.setBasis (pawBasis);
   // --- assemble data: G-vector coefficients, partial wave coefficients
   for (int n = 0; n < nPsi; ++n)  {
      res( SxIdx(N * n     , N * n + ng - 1) ) <<= psi.colRef (n);
      res( SxIdx(N * n + ng, N * (n+1) - 1) )  <<= proj.colRef (n);
   }
   return res;
}

SxDiracVec<Complex16> SxGBasis::toAny (const SxBasis *basis,
                                       const SxDiracVec<Complex16> &in) const
{
   // --- identity ?
   if (dynamic_cast<const SxGBasis*> (basis))  {
      SX_CHECK (basis == this);
      return in;
   }
   // --- PAW basis?
   if (const SxPAWBasis* pawBasis = dynamic_cast<const SxPAWBasis*> (basis))
      return toPAW (pawBasis, in);

   SX_EXIT; return in;
}

void SxGBasis::setupRealYlm (int lmax) const
{
   SX_CHECK (lmax >= 0, lmax);
   int nl = lmax + 1;
   int nLm = nl*nl;
   SX_CHECK (ng >= 0, ng);
   SxVector<Double> Ylm(nLm);

   realYlm.reformat(ng, nLm);

   // set up Ylm for all G vectors
   int lm, ig = 0;
   if (g2(0) < 1e-12)  {
      // exclude G+k = 0
      realYlm(0,0) = 1.;
      for (lm = 1; lm < nLm; lm++)
         realYlm(0,lm) = 0.;
      ig++;
   }
   while (ig < ng) {
      SxYlm::getYlmArray(lmax, gVec(ig,0), gVec(ig,1), 
                         gVec(ig,2), &Ylm);
      for (lm = 0; lm < nLm; lm++)
         realYlm(ig,lm) = Ylm(lm);
      ig++;
   }
   realYlm.setBasis (this);
}

SxDiracVec<Double> SxGBasis::getYlm(int l, int m) const
{
   int lmMax = (int)realYlm.nCols ();
   int lm = SxYlm::combineLm(l,m);
   if (lm >= lmMax) setupRealYlm (l);

   return realYlm.colRef(lm);
}
