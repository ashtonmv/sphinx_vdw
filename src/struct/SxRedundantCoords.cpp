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

#include <SxRedundantCoords.h>
#include <SxGrid.h>
#include <SxNeighbors.h>
#include <SxSortedList.h>
#include <SxElemDB.h>

SxRedundantCoords::SxRedundantCoords ()
   : nBond (-1),
     nAngle (0)
{
   primarySetLimit = 0.05;
   rmsThreshold = 3.;
   planeCutLimit = 0.95;
   maxDist = 10.;

   verbose = false;
}

namespace {
   // An auxiliary class for setting up the nearest neighbors
   // in SxRedundantCoords::setup. It is used to store neighbor
   // candidate information
   class SxNearN {
      public:
         ssize_t jTl;
         Coord relPos;
         int type;
         bool done;
         /// Constructor
         SxNearN (ssize_t jTlIn, const Coord &relPosIn, int typeIn)
         : jTl(jTlIn), relPos(relPosIn), type(typeIn), done(false)
         { }
         SxNearN () {}

         // --- make this class sortable
         bool operator> (const SxNearN &in) const {
            return relPos.normSqr () > in.relPos.normSqr ();
         }
         bool operator< (const SxNearN &in) const {
            return relPos.normSqr () < in.relPos.normSqr ();
         }
   };

   // combine unsorted (i,j) into a unique index
   inline ssize_t getPairIdx (ssize_t i, ssize_t j)
   {
      if (i < j) { ssize_t k = i; i = j; j = k; }
      return ((i * (i + 1)) >> 1) + j;
   }
}


SxArray<ssize_t>
SxRedundantCoords::classifyBonds (const SxArray<double> &pairDist,
                                  SxList<ssize_t> *setStartPtr) const
{
   SX_CHECK (setStartPtr);
   SxList<ssize_t> &setStart = *setStartPtr;
   setStart.removeAll ();
   // the cycle set stores sets as cyclic index loops
   // cycleSet(i) contains the index of the next member
   // the last member refers back to the first member
   // this makes it easy to join sets
   SxArray<ssize_t> cycleSet(pairDist.getSize ());
   SX_LOOP(i) cycleSet(i) = i;

   // --- put all distances that differ little (as defined by
   //     primarySetLimit) into one set
   for (ssize_t i = 1; i < pairDist.getSize (); i++)  {
      // distances differ by less than ~5% ?
      if (pairDist(i) - pairDist(i-1) < 2. * primarySetLimit)  {
         // join sets
         ssize_t j = cycleSet(i);
         cycleSet(i) = cycleSet(i-1);
         cycleSet(i-1) = j;
      }
   }

   // --- acquire basic statistics on primary sets
   SxList<double> dSum, d2Sum;
   SxList<ssize_t> nSet;
   {
      for (ssize_t start = 0; start < pairDist.getSize (); start++)  {
         double dSum_ = 0., d2Sum_ = 0.;
         ssize_t nSet_ = 0;
         ssize_t idx = start;
         // run through cycleSet until we hit the start
         // or any part below the start (which we have already
         // encountered)
         do {
            dSum_ += pairDist(idx);
            d2Sum_ += sqr(pairDist(idx));
            nSet_++;
            idx = cycleSet(idx);
         } while (idx > start);
         if (idx == start)  {
            // --- found a new set
            dSum  << dSum_;
            d2Sum << d2Sum_;
            nSet  << nSet_;
            setStart << start;
            // --- output
            if (verbose)  {
               cout << "Primary set " << (nSet.getSize ())
                    << ": d=" << exp(0.5 * dSum_ / double(nSet_));
               double deltaSum = fabs(d2Sum_ - sqr(dSum_) / double(nSet_));
               if (nSet_ > 1)
                  cout << " rms rel="
                       << 0.5 * sqrt(deltaSum / double(nSet_ - 1)) ;
               cout << " n=" << nSet_ << endl;
            }
         }
      }
   }

   // --- join primary sets
   for (ssize_t i = 1; i < nSet.getSize (); i++)  {
      double msNow = d2Sum(i-1) - sqr(dSum(i-1))/double(nSet(i-1));
      if (nSet(i-1) > 1) msNow /= double(nSet(i-1) - 1);
      double msNext = d2Sum(i) - sqr(dSum(i))/double(nSet(i));
      if (nSet(i) > 1) msNext /= double(nSet(i) - 1);
      // check distance of centers compared to joined rms
      if (sqr(dSum(i) - dSum(i-1)) < sqr(rmsThreshold) * (msNext + msNow))  {
         // join these sets
         {
            ssize_t start = setStart(i-1),
                    startNext = setStart(i);
            ssize_t j = cycleSet(start);
            cycleSet(start) = cycleSet(startNext);
            cycleSet(startNext) = j;
         }

         // update statistics
         dSum(i-1) += dSum(i);
         d2Sum(i-1) += d2Sum(i);
         nSet(i-1) += nSet(i);

         // remove joined set's statistics
         dSum.remove (i);
         d2Sum.remove (i);
         nSet.remove (i);
         setStart.remove (i);

         if (verbose) cout << "Joined set " << i << " with next one" << endl;
         // go back 2 steps, maybe we now join with previous set
         if (i > 1) i -= 2; else i = 0;
      }
   }
   return cycleSet;
}

bool SxRedundantCoords::verifyClasses (const SxAtomicStructure &str) const
{
   SxArray<SxStack<double> > pairDist_(nBond);
   // --- acquire all bond lengthes of currently known bonds
   SX_LOOP(iTl)  {
      if (idAtom(iTl).getSize () == 0) continue;
      SX_LOOP(iN)  {
         ssize_t jTl = idAtom(iTl)(iN);
         Coord rel = str.constRef(jTl) + latShift(iTl)(iN) - str.constRef(iTl);
         pairDist_(ricType(iTl)(iN)) << log(rel.normSqr ());
      }
   }

   // switch of verbose output
   bool verboseSafe = verbose;
   verbose = false;

   // --- rerun classify with current bond type
   bool ok = true;
   SX_LOOP(iType)  {
      SxArray<double> pairDist = pairDist_(iType);
      pairDist.sort ();
      SxList<ssize_t> startIdx;
      classifyBonds(pairDist, &startIdx);
      // if the current bond type splits into several classes, the
      // redundant coordinate search should be restarted
      ok = (startIdx.getSize () == 1);
      if (!ok) break;
   }

   // restore verbosity
   verbose = verboseSafe;

   return ok;
}

SxArray<SxSortedList<SxNearN> >
SxRedundantCoords::getCandidates(const SxAtomicStructure &structure)
{
   const SxArray<SxString> &elem
      = structure.atomInfo->meta.get (SxAtomicStructure::Elements);

   int nSpecies = structure.getNSpecies ();
   ssize_t nAtoms = structure.getNAtoms ();
   int nElemPairs = nSpecies * (nSpecies + 1) / 2;

   // temporary neighbor lists, sorted by element combinations
   SxArray<SxStack<long> > pairIJ_(nElemPairs);
   SxArray<SxStack<Coord> > relPos_(nElemPairs);
   SxArray<SxStack<double> > pairDist_(nElemPairs);

   // --- collect all neighbor relations
   {
      SxGrid grid(structure, 20);
      SxNeighbors nn;
      int nnParam = SxNeighbors::StoreDistSqr | SxNeighbors::StoreIdx
                  | SxNeighbors::StoreRel;
      SX_LOOP2(is,ia)  {
         ssize_t iTl = structure.getIAtom((int)is,(int)ia);
         nn.compute (grid, structure, structure(is,ia), maxDist, nnParam);
         if (nn.getSize () == 0)  {
            cout << "Warning: atom " << elem(is) << (ia+1)
                 << " @ " << structure(is,ia) << " has no neighbors." << endl;
            continue;
         }
         SX_LOOP(jN)  {
            int jTl = nn.idx(jN);
            if (jTl < iTl) continue;
            ssize_t js = structure.getISpecies (jTl);
            ssize_t iElemPair = getPairIdx (is, js);
            pairIJ_(iElemPair) << ( (long(iTl)<<32) + jTl);
            pairDist_(iElemPair) << log(nn.distSqr(jN));
            relPos_(iElemPair) << nn.relPositions(jN);
         }
      }
   }

   // --- sort interatomic distances for each combination of elements
   //     into sets ("types")
   //     and setup up neighbor candidate list for all atoms
   SxArray<SxSortedList<SxNearN> > candidates(nAtoms);
   nBond = 0;
   SX_LOOP(iElemPair)  {
      if (pairIJ_(iElemPair).getSize () == 0) continue;
      SxArray<long> pairIJ = pairIJ_(iElemPair);
      SxArray<double> pairDist = pairDist_(iElemPair);
      SxArray<Coord> relPos = relPos_(iElemPair);
      if (verbose) {
         int iTl = int(pairIJ(0) >> 32);
         int jTl = int(pairIJ(0) & (~int(0)));
         int is = structure.getISpecies (iTl);
         int js = structure.getISpecies (jTl);
         cout << elem(is) << "-" << elem(js) << endl;
      }

      // --- sort by distance
      if (pairDist.getSize () > 1) {
         SxArray<ssize_t> sortIdx = pairDist.getSortIdx ();
         pairIJ.sortByIdx (sortIdx);
         pairDist.sortByIdx (sortIdx);
         relPos.sortByIdx (sortIdx);
      }

      SxList<ssize_t> setStart;
      // the cycle set stores sets as cyclic index loops
      // cycleSet(i) contains the index of the next member
      // the last member refers back to the first member
      SxArray<ssize_t> cycleSet = classifyBonds (pairDist, &setStart);

      // --- now set up the candidates for all atoms with type information
      for (ssize_t i = 0; i < setStart.getSize (); i++)  {
         ssize_t start = setStart(i);
         ssize_t idx = start;
         // --- run through set
         do {
            int iTl = int(pairIJ(idx) >> 32);
            int jTl = int(pairIJ(idx) & (~int(0)));
            candidates(iTl) << SxNearN (jTl, relPos(idx), nBond);
            if (iTl != jTl)
               candidates(jTl) << SxNearN (iTl, -relPos(idx), nBond);
            idx = cycleSet(idx);
         } while (idx != start);
         nBond++;
      }
   }
   return candidates;
}

void SxRedundantCoords::setup (const SxAtomicStructure &structure)
{
   const SxArray<SxString> &elem
      = structure.atomInfo->meta.get (SxAtomicStructure::Elements);

   // --- resize
   {
      ssize_t nAtoms = structure.getNAtoms ();
      idAtom.resize (nAtoms);
      ricType.resize (nAtoms);
      latShift.resize (nAtoms);
   }

   SxArray<SxSortedList<SxNearN> > candidates = getCandidates (structure);

   // --- map preliminary types to types finally used
   SxArray<int> finalType(nBond);
   finalType.set (-1);
   int nTypeFinal = 0;

   // --- now run over candidates
   SX_LOOP(iTl)  {
      SxList<ssize_t> nearestNeighbors;
      if (verbose)
         cout << "iTl=" << iTl << ": " << candidates(iTl).getSize ()
              << " candidates" << endl;
      for (ssize_t iC = 0; iC < candidates(iTl).getSize (); ++iC)  {
         SxNearN &pair = candidates(iTl)(iC);
         if (pair.done) continue;
         //cout << pair.relPos << " (" << pair.type << ")";

         // --- check if this candidate is inside the polyhedron
         //     defined by the current neighbors
         bool keep = true;
         for (ssize_t jC = 0; jC < nearestNeighbors.getSize (); ++jC)  {
            Coord relPosJ = candidates(iTl)(nearestNeighbors(jC)).relPos;
            if ((pair.relPos ^ relPosJ)/relPosJ.normSqr () > planeCutLimit)  {
               // this candididate lies beyond the half-plane
               // defined by jC => skip it
               //cout << " <=> " << relPosJ;
               keep=false;
               break;
            }
         }

         //cout << (keep ? " ok" : " skipped" ) << endl;

         // --- if so, keep this and all its siblings
         if (keep)  {
            pair.done = true;
            nearestNeighbors << iC;
            // also keep all others of the same type
            for (ssize_t jC = 0; jC < candidates(iTl).getSize (); ++jC)  {
               SxNearN &pairJ = candidates(iTl)(jC);
               if (pairJ.done) continue;
               if (pairJ.type == pair.type)  {
                  pairJ.done = true;
                  nearestNeighbors << jC;
               }
            }
         }
      }

      // --- copy the final list of neighbors
      ssize_t nKeep = nearestNeighbors.getSize ();
      idAtom(iTl).resize (nKeep);
      ricType(iTl).resize (nKeep);
      latShift(iTl).resize (nKeep);

      if (verbose)  {
         cout << nKeep << " neighbors for "
              << elem(structure.getISpecies((int)iTl))
              << " @ " << structure(iTl) << endl;
      }
      for (ssize_t jC = 0; jC < nKeep; ++jC)  {
         SxNearN &pair = candidates(iTl)(nearestNeighbors(jC));
         idAtom(iTl)(jC) = pair.jTl;
         int type = finalType(pair.type);
         if (type == -1)
            finalType(pair.type) = type = nTypeFinal++;
         ricType(iTl)(jC) = type;
         latShift(iTl)(jC) = pair.relPos + structure(iTl)
                           - structure((int)pair.jTl);
         if (verbose)  {
            cout << "d=" << pair.relPos.norm () << " ";
            cout << elem(structure.getISpecies (int(pair.jTl))) << " @ "
                 << pair.relPos << " type=" << (type + 1) << endl;
         }
      }
   }
   nBond = nTypeFinal;

   /// --- gather final statistics
   bondTypes.resize (nBond);
   SxArray<int> nCase(nBond);
   nCase.set (0);
   SX_LOOP(iType)
      bondTypes(iType).dAvg = 0.;
   SX_LOOP(iTl) {
      if (ricType(iTl).getSize () == 0) continue;
      SX_LOOP(iN)  {
         int iType = ricType(iTl)(iN);
         nCase(iType)++;
         ssize_t jTl = idAtom(iTl)(iN);
         int is = structure.getISpecies ((int)iTl);
         int js = structure.getISpecies ((int)jTl);
         bondTypes(iType).iSpecies1 = min(is, js);
         bondTypes(iType).iSpecies2 = max(is, js);
         Coord relI = structure.constRef(jTl) + latShift(iTl)(iN)
                    - structure.constRef(iTl);
         bondTypes(iType).dAvg += relI.norm ();
      }
   }
   SX_LOOP(iType)
      bondTypes(iType).dAvg /= double(nCase(iType));

}

void SxRedundantCoords::getAngles (const SxAtomicStructure &structure)
{
   const SxArray<SxString> &elem
      = structure.atomInfo->meta.get (SxAtomicStructure::Elements);

   int nSpecies = structure.getNSpecies ();
   int nBondPairs = (nBond * (nBond + 1) / 2) * nSpecies;

   // temporary neighbor lists, sorted by bond type combinations
   SxArray<SxStack<long> > pairId_(nBondPairs);
   SxArray<SxStack<double> > cosPhi_(nBondPairs);

   // --- collect all neighbor relations
   SX_LOOP(iTl) {
      ssize_t nN = ricType(iTl).getSize ();
      if (nN <= 0) continue;
      SxArray<Coord> dir(nN);
      int iSpecies = structure.getISpecies ((int)iTl);
      SX_LOOP(iN)
         dir(iN) = (structure.constRef(idAtom(iTl)(iN)) + latShift(iTl)(iN)
                    - structure.constRef(iTl)).normalize ();

      for (int iN = 0; iN < nN; ++iN) {
         for (int jN = iN + 1; jN < ricType(iTl).getSize (); ++jN)  {
            ssize_t pairType = (getPairIdx (ricType(iTl)(iN), ricType(iTl)(jN))) * nSpecies + iSpecies;
            pairId_(pairType) << (iTl << 32) + (iN << 16) + jN;
            cosPhi_(pairType) << (dir(iN) ^ dir(jN));
         }
      }
   }

   // --- sort interatomic bond angles for each combination of bonds
   //     into sets ("types")
   //     and setup up angle list for all atoms
   angType.resize (structure.getNAtoms ());
   SX_LOOP(iTl) {
      ssize_t nBonds = ricType(iTl).getSize ();
      angType(iTl).resize (nBonds * (nBonds + 1)/2);
      for (ssize_t iB = 0; iB < nBonds; ++iB)
         angType(iTl)(iB * (iB + 1) / 2) = -1;
   }
   nAngle = 0;
   double cosSetLimit = 0.1;
   SxList<AngleType> typeList;
   SX_LOOP(iBondPair)  {
      if (pairId_(iBondPair).getSize () == 0) continue;
      SxArray<long> pairId = pairId_(iBondPair);
      SxArray<double> cosPhi = cosPhi_(iBondPair);
      int ibt = 0, jbt = int(iBondPair / nSpecies);
      while (jbt > ibt) {
         ibt++;
         jbt = int(iBondPair/nSpecies) - ibt*(ibt + 1) / 2;
      }
      if (verbose) {
         int is = int(iBondPair % nSpecies);
         cout << "bond types " << (ibt + 1) << "-" << elem(is) << "-"
              << (jbt + 1) << endl;
      }

      // --- sort by distance
      if (cosPhi.getSize () > 1) {
         SxArray<ssize_t> sortIdx = cosPhi.getSortIdx ();
         pairId.sortByIdx (sortIdx);
         cosPhi.sortByIdx (sortIdx);
      }

      // the cycle set stores sets as cyclic index loops
      // cycleSet(i) contains the index of the next member
      // the last member refers back to the first member
      // this makes it easy to join sets
      SxArray<ssize_t> cycleSet(cosPhi.getSize ());
      SX_LOOP(i) cycleSet(i) = i;

      // --- put all distances that differ little (as defined by
      //     primarySetLimit) into one set
      for (ssize_t i = 1; i < cosPhi.getSize (); i++)  {
         // cosines differ by less than limit ?
         if (cosPhi(i) - cosPhi(i-1) < cosSetLimit)  {
            // join sets
            ssize_t j = cycleSet(i);
            cycleSet(i) = cycleSet(i-1);
            cycleSet(i-1) = j;
         }
      }

      // --- acquire basic statistics on primary sets
      for (ssize_t start = 0; start < cosPhi.getSize (); start++)  {
         double dSum = 0., d2Sum = 0.;
         ssize_t nSet = 0;
         ssize_t idx = start;
         // run through cycleSet until we hit the start
         // or any part below the start (which we have already
         // encountered)
         do {
            dSum += cosPhi(idx);
            d2Sum += sqr(cosPhi(idx));
            nSet++;
            idx = cycleSet(idx);
         } while (idx > start);
         if (idx == start)  {
            // --- found a new set
            int type = nAngle;
            if (dSum/double(nSet) < -0.95)
               type = -1;
            else
               nAngle++;
            do {
               int iTl = int(pairId(idx) >> 32);
               int iN  = int(pairId(idx) & 0xffff0000) >> 16;
               int jN  = int(pairId(idx) & 0x0000ffff);
               angType(iTl)(getPairIdx(iN, jN)) = type;
               idx = cycleSet(idx);
            } while (idx != start);
            // --- output
            if (verbose)  {
               if (type >=0)
                  cout << "New set " << nAngle;
               else
                  cout << "Neglect set";
               cout << " of bond angles: phi="
                    << acos(dSum/double(nSet)) * 180. / PI;
               double deltaSum = fabs(d2Sum - sqr(dSum) / double(nSet));
               if (nSet > 1)
                  cout << " rms cos(phi)="
                       << sqrt(deltaSum / double(nSet - 1)) ;
               cout << " n=" << nSet << endl;
            }
            if (type != -1)  {
               int iTl = int(pairId(start) >> 32);
               int iN  = int(pairId(start) & 0xffff0000) >> 16;
               int jN  = int(pairId(start) & 0x0000ffff);
               typeList << AngleType ();
               typeList.last ().iSpeciesCentral = structure.getISpecies (iTl);
               int js = structure.getISpecies ((int)idAtom(iTl)(iN));
               int ks = structure.getISpecies ((int)idAtom(iTl)(jN));
               if (js > ks)  {
                  int x = js; js = ks; ks = x;
               }
               typeList.last ().jSpecies1 = js;
               typeList.last ().jSpecies2 = ks;
               typeList.last ().iBondType1 = ibt;
               typeList.last ().iBondType2 = jbt;
               typeList.last ().cosAvg = dSum / double(nSet);
            }
         }
      }
   }
   angleTypes = typeList;
}

void SxRedundantCoords::getBornVonKarmanAngles (const SxArray<int> &atoms)
{
   SxArray<int> type(nBond);
   type.set (-1);
   int nBvK = getNParam ();
   
   if (angType.getSize () != ricType.getSize ())
      angType.resize (ricType.getSize ());
   SX_LOOP(iBvK) {
      ssize_t iTl = atoms(iBvK);
      int nNeighbors = (int)idAtom(iTl).getSize ();
      int offsetBvK = nNeighbors * (nNeighbors + 1) / 2;
      if (angType(iTl).getSize () == 0)  {
         angType(iTl).resize (offsetBvK);
         angType(iTl).set (-1);
      }
      angType(iTl).resize (nNeighbors * (nNeighbors + 3)/2, true);
      for (int iN = 0; iN < nNeighbors; ++iN)  {
         int bondType = ricType(iTl)(iN);
         if (type(bondType) == -1)  {
            type(bondType) = nBvK++;
            cout << "Added Born-von-Karman bond orientation term for type=" 
                 << (bondType + 1) << " as parameter " << nBvK << endl;
         }
         angType(iTl)(offsetBvK + iN) = type(bondType);
      }
   }
   // update number of angle parameters
   nAngle = nBvK - nBond;
}


namespace {
   void collisionCrash (const SxAtomicStructure &str,
                        ssize_t iTl, ssize_t jTl, const Coord &posJ)
   {
      cout << "Atom collision" << endl;
      cout << iTl << ' ' << jTl << endl;
      cout << str.constRef(iTl) << " - " << posJ << endl;
      SX_EXIT;
   }

   inline
   Coord getDx(const SxVector<Double> &x, ssize_t iTl,  ssize_t jTl)
   {
      return  Coord::toVec3Ref(x.elements + 3 * jTl)
            - Coord::toVec3Ref(x.elements + 3 * iTl);
   }
}

SxVector<Double> SxRedundantCoords::applyH (const SxAtomicStructure &str,
                                            const SxVector<Double> &x) const
{
   SX_CHECK (str.getNAtoms () == idAtom.getSize (),
             str.getNAtoms (), idAtom.getSize ());
   SX_CHECK(x.getSize () == 3 * str.getNAtoms (),
            x.getSize (), str.getNAtoms ());
   SxVector<Double> res(x.getSize ());
   res.set (0.);
   res.reshape (max(x.nRows (),1L), max(x.nCols (),1L));
   SX_LOOP(iTl)  {
      int nN = (int)idAtom(iTl).getSize ();
      if (nN == 0) continue;
      int offsetBvK = nN * (nN + 1) / 2;
      SX_LOOP(iN)  {
         ssize_t jTl = idAtom(iTl)(iN);
         Coord posJ = str.constRef(jTl) + latShift(iTl)(iN),
               rel  = posJ - str.constRef(iTl);
         double rel2 = rel.normSqr ();
         if (rel2 < 1e-6) collisionCrash (str, iTl, jTl, posJ);
         Coord dx = getDx(x, iTl, jTl);
         double D = param(ricType(iTl)(iN));
         // get force along direction vector
         Coord dr_dx = rel * ((rel ^ dx) * D/rel2);
         // add force to result
         Coord::toVec3Ref(res.elements + 3 * iTl) += dr_dx;
         Coord::toVec3Ref(res.elements + 3 * jTl) -= dr_dx;

         // --- bond angle part
         if (angType.getSize () > 0 && angType(iTl).getSize () > 0)  {
            // the potential is E = Da * (1 - cos (phi-phi0) );
            // the Hesse at phi=phi0 is Da * (dPhi/dx_i) (dPhi/d x_j)
            // let t' = R_J - R_I, u' = R_K - R_I
            // t = t'/|t'| and u = u'/|u'|
            // cos phi = t.u
            // => (-sin phi) dphi = [d/d x_i (t.u)] d x_i
            // d/d R_j (t.u) = u.(d (t'/|t'|) / dR_J)
            //                 = u/|t'| - (u.t) t/|t'|

            double tNormInv = 1./sqrt(rel2);
            Coord t = rel * tNormInv;
            for (ssize_t kN = iN + 1; kN < idAtom(iTl).getSize (); ++kN) {
               ssize_t jParam = angType(iTl)(getPairIdx(iN,kN));
               if (jParam == -1) continue;
               ssize_t kTl = idAtom(iTl)(kN);
               Coord relK = str.constRef(kTl) + latShift(iTl)(kN)
                          - str.constRef(iTl);
               double uNormInv = 1./relK.norm ();
               Coord u = relK * uNormInv;
               double cosPhi = t ^ u, sin2Phi = (1. - cosPhi) * (1. + cosPhi);
               if (sin2Phi > 1e-12) {
                  double Da = param(nBond + jParam);
                  double sinPhiInv = -1./sqrt(sin2Phi);
                  Coord tProj = tNormInv * sinPhiInv * (u - cosPhi * t);
                  Coord uProj = uNormInv * sinPhiInv * (t - cosPhi * u);
                  Coord dxK = getDx(x, iTl, kTl);
                  double Hx = ((tProj ^ dx) + (uProj ^ dxK)) * Da;
                  Coord fT = tProj * Hx;
                  Coord fU = uProj * Hx;
                  Coord::toVec3Ref(res.elements + 3 * iTl) += fT + fU;
                  Coord::toVec3Ref(res.elements + 3 * jTl) -= fT;
                  Coord::toVec3Ref(res.elements + 3 * kTl) -= fU;
               }
            }
            if (angType(iTl).getSize () > offsetBvK)  {
               // Born-von-Karman bond-orthogonal force
               Coord dXortho = dx - rel * ((rel ^ dx) / rel2);
               Coord fOrtho = param(angType(iTl)(offsetBvK + iN)) * dXortho;
               Coord::toVec3Ref(res.elements + 3 * iTl) += fOrtho;
               Coord::toVec3Ref(res.elements + 3 * jTl) -= fOrtho;
            }
         }
      }
   }
   return res;
}

SxVector<Double>
SxRedundantCoords::getParamDeriv (const SxAtomicStructure &str,
                                  const SxVector<Double> &x) const
{
   SX_CHECK (str.getNAtoms () == idAtom.getSize (),
             str.getNAtoms (), idAtom.getSize ());
   SX_CHECK(x.getSize () == 3 * str.getNAtoms (),
            x.getSize (), str.getNAtoms ());
   ssize_t nDof = x.getSize ();
   SxVector<Double> res;
   res.reformat (nDof, nBond);
   res.set (0.);
   SX_LOOP(iTl)  {
      if (idAtom(iTl).getSize () == 0) continue;
      SX_LOOP(iN)  {
         ssize_t jTl = idAtom(iTl)(iN);
         Coord posJ = str.constRef(jTl) + latShift(iTl)(iN),
               rel  = posJ - str.constRef(iTl);
         double rel2 = rel.normSqr ();
         if (rel2 < 1e-6) collisionCrash (str, iTl, jTl, posJ);
         Coord dx = getDx(x, iTl, jTl);
         int iParam = ricType(iTl)(iN);
         // get projection onto direction vector
         Coord dr_dx = rel * ((rel ^ dx)/rel2);
         // add force to result
         Coord::toVec3Ref(res.elements + 3 * iTl + nDof * iParam) += dr_dx;
         Coord::toVec3Ref(res.elements + 3 * jTl + nDof * iParam) -= dr_dx;
      }
   }
   return res;
}

SxVector<Double>
SxRedundantCoords::getParamDerivA (const SxAtomicStructure &str,
                                   const SxVector<Double> &x) const
{
   SX_CHECK (str.getNAtoms () == idAtom.getSize (),
             str.getNAtoms (), idAtom.getSize ());
   SX_CHECK(x.getSize () == 3 * str.getNAtoms (),
            x.getSize (), str.getNAtoms ());
   SX_CHECK (nAngle > 0);
   ssize_t nDof = x.getSize ();
   SxVector<Double> res;
   res.reformat (nDof, nAngle);
   res.set (0.);
   SX_LOOP(iTl)  {
      int nN = (int)idAtom(iTl).getSize ();
      if (nN == 0) continue;
      int offsetBvK = nN * (nN + 1) / 2;
      SX_LOOP(iN)  {
         ssize_t jTl = idAtom(iTl)(iN);
         Coord posJ = str.constRef(jTl) + latShift(iTl)(iN),
               rel  = posJ - str.constRef(iTl);
         double rel2 = rel.normSqr ();
         if (rel2 < 1e-6) collisionCrash (str, iTl, jTl, posJ);
         Coord dx = getDx(x, iTl, jTl);
         // --- bond angle part
         // the potential is E = Da * (1 - cos (phi-phi0) );
         // the Hesse at phi=phi0 is Da * (dPhi/dx_i) (dPhi/d x_j)
         // let t' = R_J - R_I, u' = R_K - R_I
         // t = t'/|t'| and u = u'/|u'|
         // cos phi = t.u
         // => (-sin phi) dphi = [d/d x_i (t.u)] d x_i
         // d/d R_j (t.u) = u.(d (t'/|t'|) / dR_J)
         //                 = u/|t'| - (u.t) t/|t'|

         double tNormInv = 1./sqrt(rel2);
         Coord t = rel * tNormInv;
         for (ssize_t kN = iN + 1; kN < idAtom(iTl).getSize (); ++kN) {
            ssize_t jParam = angType(iTl)(getPairIdx(iN,kN));
            if (jParam == -1) continue;
            ssize_t kTl = idAtom(iTl)(kN);
            Coord relK = str.constRef(kTl) + latShift(iTl)(kN)
                       - str.constRef(iTl);
            double uNormInv = 1./relK.norm ();
            Coord u = relK * uNormInv;
            double cosPhi = t ^ u, sin2Phi = (1. - cosPhi) * (1. + cosPhi);
            if (sin2Phi > 1e-12) {
               double sinPhiInv = -1./sqrt(sin2Phi);
               Coord tProj = tNormInv * sinPhiInv * (u - cosPhi * t);
               Coord uProj = uNormInv * sinPhiInv * (t - cosPhi * u);
               Coord dxK = getDx(x, iTl, kTl);
               double Hx = ((tProj ^ dx) + (uProj ^ dxK));
               Coord fT = tProj * Hx;
               Coord fU = uProj * Hx;
               ssize_t offset = nDof * jParam;
               Coord::toVec3Ref(res.elements + 3 * iTl + offset) += fT + fU;
               Coord::toVec3Ref(res.elements + 3 * jTl + offset) -= fT;
               Coord::toVec3Ref(res.elements + 3 * kTl + offset) -= fU;
            }
         }
         if (angType(iTl).getSize () > offsetBvK)  {
            // Born-von-Karman bond-orthogonal force
            Coord dXortho = dx - rel * ((rel ^ dx) / rel2);
            ssize_t offset = nDof * (angType(iTl)(offsetBvK + iN) - nBond);
            Coord::toVec3Ref(res.elements + 3 * iTl + offset) += dXortho;
            Coord::toVec3Ref(res.elements + 3 * jTl + offset) -= dXortho;
         }
      }
   }
   return res;
}

SxVector<Double> SxRedundantCoords::solve (const SxAtomicStructure &str,
                                           const SxVector<Double> dF,
                                           double normWeight,
                                           double accuracy) const
{
   /*
   // --- direct solution
   ssize_t nDof = dF.getSize ();
   SxMatrix<Complex16> H(nDof, nDof);
   SxVector<Double> x(nDof);
   x.set (0.);
   SX_LOOP(i) {
      x(i) = 1.;
      H.colRef (i) <<= applyH(str, x);
      x(i) = 0.;
   }
   H += H.transpose ();
   H *= 0.5;
   SxMatrix<Complex16>::Eigensystem eig = H.eigensystem ();
   SxVector<Complex16> dots = eig.vecs.overlap (dF);
   cout << eig.vals << endl;
   double N = dots.normSqr ();
   double normWeight2 = normWeight * normWeight;

   SX_LOOP(i)  {
      //cout << eig.vals(i).re << ": " << dots(i).absSqr () / N << "->";
      dots(i) *= eig.vals(i) / (normWeight2 + eig.vals(i).absSqr ());
      //cout << dots(i).absSqr () / N << endl;
   }
   SxVector<Double> res = eig.vecs ^ dots;
   cout << "Solve: |dF|=" << dF.norm () << endl;
   cout << "Solve: |x|=" << res.norm () << endl;
   cout << "Solve: |dF - Hx|=" << (dF - applyH (str, res)).norm () << endl;
   return res;
   */

   SxVector<Double> rhs = applyH(str, dF);
   SxVector<Double> res = rhs.getCopy (), Xold;
   double R2old = 0.;
   double lambda = normWeight * normWeight; // Tikhonov regularization
   for (ssize_t it = 0; it < dF.getSize (); ++it)  {
      SxVector<Double> R = applyH(str, applyH(str, res)) - rhs;
      R += lambda * res;
      double R2 = R.normSqr ();
      if (R2 < sqr(accuracy)) { R2old = R2; break; }

      // conjugate gradient
      SxVector<Double> X;
      if (it > 0) {
         double gamma = R2/R2old;
         X = R + gamma * Xold;
      } else {
         X = R;
      }
      R2old = R2;
      Xold = X;

      // line opt
      SxVector<Double> dFx = applyH(str, applyH(str, X));
      dFx += lambda * X;
      //double lOpt = dot(R,dFx) / dFx.normSqr ();
      double lOpt = R2 / dot(dFx, X);
      res -= lOpt * X;
   }
   //cout << "R2 = " << R2old << endl;
   //cout << "Solve: |dF|=" << dF.norm () << endl;
   //cout << "Solve: |x|=" << res.norm () << endl;
   //cout << "Solve: |dF - Hx|=" << (dF - applyH (str, res)).norm () << endl;
   return res;
}

void SxRedundantCoords::paramSchlegel (const SxArray<SxString> &elem)
{
   // See B. Schlegel, Theor. Chim. Act. 66, 333 (1984). Table I
   param.resize (getNParam ());
   SxElemDB pse;
   for (ssize_t i = 0; i < nBond; i++)  {
      const BondType &bt = bondTypes(i);
      int row1 = pse.getRow (elem(bt.iSpecies1));
      if (row1 < 1 || row1 > 3) row1 = 3;
      int row2 = pse.getRow (elem(bt.iSpecies2));
      if (row2 < 1 || row2 > 3) row2 = 3;
      double A = 1.734 / 2.; // divide by 2 due to doublecounting
      double B = 0.;
      if (row1 == 1 && row2 == 1)  {
         B = -0.244;
      } else if (row1 + row2 == 3)  {
         B = 0.352; // 1--2 or 2--1
      } else if (row1 == 2 && row2 == 2)  {
         B = 1.085;
      } else if (row1 + row2 == 4)  {
         B = 0.660; // 1--3 or 3--1
      } else if (row1 + row2 == 5)  {
         B = 1.552; // 2--3 or 3--2
      } else if (row1 == 3 && row2 == 3)  {
         B = 2.068;
      } else {
         SX_EXIT;
      }
      double rB = bt.dAvg - B;
      param(i) = A / (rB * rB * rB);
   }
   for (ssize_t i = nBond; i < getNParam (); i++)  {
      const AngleType &at = angleTypes(i-nBond);
      bool withH =  pse.getAtomicNumber(elem(at.jSpecies1)) == 1
                 || pse.getAtomicNumber(elem(at.jSpecies2)) == 1;
      param(i) = withH ? 0.160 : 0.250;
   }
}

void SxRedundantCoords::paramFischer (const SxArray<SxString> &elem)
{
   // See T.H. Fischer, J. Almloef, J. Phys. Chem. 96,9768-9774 (1992), Table I
   param.resize (getNParam ());
   SxElemDB pse;
   for (ssize_t i = 0; i < nBond; i++)  {
      const BondType &bt = bondTypes(i);
      double A = 0.3601;
      double B = 1.944;
      double rcov = pse.getCovalentRadius(elem(bt.iSpecies1))
                  +  pse.getCovalentRadius(elem(bt.iSpecies1));
      param(i) = A * exp(-B * (bt.dAvg - rcov));
   }
   for (ssize_t i = nBond; i < getNParam (); i++)  {
      const AngleType &at = angleTypes(i-nBond);
      double rc1 = pse.getCovalentRadius (at.iSpeciesCentral),
             rc2 = pse.getCovalentRadius (at.jSpecies1),
             rc3 = pse.getCovalentRadius (at.jSpecies2);
      double d1 =  bondTypes(at.iBondType1).dAvg; // might be d12 or d13
      double d2 =  bondTypes(at.iBondType2).dAvg;

      double A = 0.089, B = 0.11, C = 0.44, D = -0.42;

      double rcov12 = rc1 + rc2,
             rcov13 = rc1 + rc3;
      param(i) = A
               + B/pow(rc1 * rc2, D) * exp(-C * (d1 + d2 - rcov12 - rcov13));
   }
}
