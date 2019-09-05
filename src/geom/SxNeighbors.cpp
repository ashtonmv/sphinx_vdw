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

#include <SxNeighbors.h>
#include <SxStack.h>

inline int floorDivide(int a, int b)
{
   int res = a / b;
   if (a < res * b) --res;
   return res;
}

void SxNeighbors::compute (const SxGrid            &grid,
                           const SxAtomicStructure &str,
                           const Coord             &r0,
                           double                  rcut,
                           int                     mode)
{

   SxVector3<Int> from, to;
   {
      // circumscribing box
      Coord f = grid.gridCell.getBoundingBox (rcut);
      f += 1.01; // safety range for numerical noise at grid cell boundaries
      // atoms at the boundaries may be registered one grid cell away
      
      // range of grid cells (mesh vector)
      Coord r0Rel = grid.gridCell.carToRel(r0);

      from = floor(r0Rel-f);
      to   = floor(r0Rel+f);
   }
   
   SxVector3<Int> meshVec, tRel;
   SxVector<Int>::Iterator idxIt, idxEnd;
   double rcut2 = rcut * rcut, d2;
   int iGrid;
   Coord delta, t, r0shifted;

   bool storePos   = (mode & StoreAbsNoInfo),
        storeDelta = (mode & StoreRelNoInfo),
        storeD2    = (mode & StoreDistSqr),
        storeIdx   = (mode & StoreIdx);
   SxArray<SxStack<Coord> >  posList, deltaList;
   SxArray<SxStack<double> > distSqrList;
   SxArray<SxStack<int> >    idxList;
   
   int nSpecies = (str.atomInfo && (mode & StoreAbsAndRel)) 
                  ? str.getNSpecies () 
                  : 1;
   if (storePos)   posList.resize (nSpecies);
   if (storeDelta) deltaList.resize (nSpecies);
   if (storeD2)    distSqrList.resize (nSpecies);
   if (storeIdx)   idxList.resize (nSpecies);

   int nFound = 0;
   double zeroLimit 
      = (mode & IncludeZeroDistance) ? -1. 
                                     : 3.*str.epsEqual*str.epsEqual;

   // --- find neighbors
   for (meshVec(0) = from(0); meshVec(0) <= to(0); meshVec(0)++) {
      tRel(0) = floorDivide (meshVec(0), grid.gridMesh(0));
      for (meshVec(1) = from(1); meshVec(1) <= to(1); meshVec(1)++) {
         tRel(1) = floorDivide (meshVec(1), grid.gridMesh(1));
         for (meshVec(2) = from(2); meshVec(2) <= to(2); meshVec(2)++) {
            tRel(2) = floorDivide (meshVec(2), grid.gridMesh(2));
            
            // compute translation vector in absolute coordinates
            t = grid.regCell.relToCar(tRel);
            
            // absolute position of possible neighbors is x+t
            // instead of shifting x, we shift r0 in the opposite direction
            r0shifted = r0 - t;

            // get grid cell
            iGrid = (int)grid.gridMesh.getMeshIdx(meshVec, SxMesh3D::Unknown);

            // --- loop over atoms in grid cell
            idxIt  = grid.getIdxIt (iGrid);
            idxEnd = grid.computeItEnd (iGrid);
            for ( ; idxIt != idxEnd; ++idxIt)  {
               // --- check distance
               delta = grid.mappedStr(*idxIt) - r0shifted;
               d2 = delta.normSqr ();
               if ( d2  <= rcut2 && d2 > zeroLimit) {
                  nFound++;
                  int is = (nSpecies > 1) ? str.getISpecies (*idxIt) : 0;
                  if (storePos)   posList(is)     << (grid.mappedStr(*idxIt)+t);
                  if (storeDelta) deltaList(is)   << delta;
                  if (storeD2)    distSqrList(is) << d2;
                  if (storeIdx)   idxList(is)     << *idxIt;
               }
            }
         }
      }
   }
   
   // --- copy collected data into member variables

   if (storeIdx)  {
      // --- atom indices
      idx = SxVector<Int> (nFound);
      int nNeighbors, offset = 0;
      for (int is = 0; is < nSpecies; ++is)  {
         nNeighbors = int(idxList(is).getSize ());
         idx.set (idxList(is), nNeighbors, offset);
         offset += nNeighbors;
      }
   }

   if (storeD2)  {
      // --- square of distances
      distSqr = SxVector<TPrecTauR> (nFound);
      int offset = 0, nNeighbors;
      for (int is = 0; is < nSpecies; ++is)  {
         nNeighbors = int(distSqrList(is).getSize ());
         distSqr.set (distSqrList(is), nNeighbors, offset);
         offset += nNeighbors;
      }
   }

   if ( (!storePos) && (!storeDelta)) return;
   
   // --- setup info for SxAtomicStructure
   SxAtomInfo::Ptr info;
   if (str.atomInfo)
      info = SxAtomInfo::derive (str.atomInfo);
   else
      info = SxAtomInfo::Ptr::create ();
   
   info->nSpecies = nSpecies;
   info->nAtoms.resize (nSpecies);
   if (storeIdx) info->parentMap = idx;
   
   if (storePos)  {
      // --- absolute neighbor positions
      for (int is = 0; is < nSpecies; ++is) 
         info->nAtoms(is) = int(posList(is).getSize ());
      info->setupOffset ();
      
      absPositions.resize (nFound);
      absPositions.atomInfo = info;
      for (int is = 0; is < nSpecies; ++is)
         absPositions.importStack(posList(is), is);
      
      //absPositions.cell     = str.cell;
      absPositions.epsEqual = str.epsEqual;
   }
   
   if (storeDelta)  {
      // --- neighbor positions relative to r0
      if (!storePos)  {
         for (int is = 0; is < nSpecies; ++is) 
            info->nAtoms(is) = int(deltaList(is).getSize ());
         info->setupOffset ();
      } 
         
      relPositions.resize (nFound);
      relPositions.atomInfo = info;
      for (int is = 0; is < nSpecies; ++is)
         relPositions.importStack (deltaList(is), is);

      //relPositions.cell     = str.cell;
      relPositions.epsEqual = str.epsEqual;
   }
}

int SxNeighbors::getSize () const
{
   ssize_t size;
   if ((size = idx.getSize ()) > 0) return (int)size;
   if ((size = distSqr.getSize ()) > 0) return (int)size;
   if ((size = absPositions.getNAtoms ()) > 0) return (int)size;
   return (int)relPositions.getNAtoms ();
}



