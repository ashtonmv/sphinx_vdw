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

#include <SxSpaceSeparator.h>

SxDiracVec<Double> SxSpaceSeparator::getGradRho(int dir)
{
   SX_CHECK(dir > -1 && dir < 3);

   return gradRhoComponents(dir);
}

SxDiracVec<Double> SxSpaceSeparator::getGradGradRho(int dir)
{
   SX_CHECK(dir > -1 && dir < 9);

   return gradGradRhoComponents(dir);
}

SxList<SxDiracVec<Double> > SxSpaceSeparator::computeGradRho ()
{
   // get basises
   const SxRBasis &R = rho.getBasis<SxRBasis> ();
   const SxGBasis &G = R.getGBasis ();

   SxList<SxDiracVec<Double> > result;

   //Transform to GSpace
   SxDiracVec<Complex16> rhoG = (G | rho);

   // compute grad = IG * ... and back to R Space
   result.append(R | I * G.gVec.colRef(0) * rhoG);
   result.append(R | I * G.gVec.colRef(1) * rhoG);
   result.append(R | I * G.gVec.colRef(2) * rhoG);

   return result;
}

SxList<SxDiracVec<Double> > SxSpaceSeparator::computeGradGradRho ()
{
   // get basises
   const SxRBasis &R = rho.getBasis<SxRBasis> ();
   const SxGBasis &G = R.getGBasis ();

   SxList<SxDiracVec<Double> > result;

   //Transform to GSpace
   SxDiracVec<Complex16> rhoG = (G | rho);

   // compute grad = -G^2 * ... and back to R Space
   result.append(R | -G.gVec.colRef(0) * G.gVec.colRef(0) * rhoG);
   result.append(R | -G.gVec.colRef(0) * G.gVec.colRef(1) * rhoG);
   result.append(R | -G.gVec.colRef(0) * G.gVec.colRef(2) * rhoG);
   result.append(R | -G.gVec.colRef(1) * G.gVec.colRef(0) * rhoG);
   result.append(R | -G.gVec.colRef(1) * G.gVec.colRef(1) * rhoG);
   result.append(R | -G.gVec.colRef(1) * G.gVec.colRef(2) * rhoG);
   result.append(R | -G.gVec.colRef(2) * G.gVec.colRef(0) * rhoG);
   result.append(R | -G.gVec.colRef(2) * G.gVec.colRef(1) * rhoG);
   result.append(R | -G.gVec.colRef(2) * G.gVec.colRef(2) * rhoG);

   return result;
}

SxSpaceSeparator::SxSpaceSeparator (const SxAtomicStructure &structureIn)
{
   structure = structureIn;
   PAWDist.resize(structure.getNSpecies());
   PAWDist.set(0.0);

}

SxArray<SxDiracVec<Double> > SxSpaceSeparator::voronoi (const SxRBasis &R)
{
   SX_CLOCK(Timer::Voronoi);
   SxVector3<Int> dims = R.getMesh ();
   mesh = SxMesh3D(dims);
   int nDim = mesh.product ();
   int nAtoms = structure.getNAtoms ();

   SxArray<SxDiracVec<Double> > result(nAtoms);
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      result(iAtom).resize(nDim);
      result(iAtom).set(0.0);
      result(iAtom).setBasis(R);
   }

   for (ssize_t idx = 0; idx < nDim; idx++)  {
      SxArray<int> nextAtoms = getNextAtom(idx);
      for (int iNN = 0; iNN < nextAtoms.getSize(); iNN++)  {
         result(nextAtoms(iNN))(idx) = 1.0 / double(nextAtoms.getSize ());
      }
   }

   return result;
}

SxArray<int> SxSpaceSeparator::getNextAtom (ssize_t index)
{
   // grid for neighbor search
   SxGrid grid (structure, 10);
   // get Coord
   Coord rel = mesh.getMeshVec(index,SxMesh3D::Positive)/Coord(mesh);
   Coord abs = structure.cell.relToCar (rel);

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
      double minDistSqr = neighbors.distSqr.minval(&minIdx);

      SxVector<Int> atomMarker (structure.getNAtoms ());
      atomMarker.set(0);
      for (ssize_t iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++)  {
         double delta = neighbors.distSqr(iNeighbor) - minDistSqr; 
         if (fabs(delta) < 1e-12)  {
            atomMarker(neighbors.idx(iNeighbor)) = 1;
         }
      }
      
      SxArray<int> result (atomMarker.sum ());
      ssize_t iNeighbor = 0;
      for (int iAtom = 0; iAtom < structure.getNAtoms (); iAtom++)  {
         if (atomMarker(iAtom) == 1)  {
            result(iNeighbor) = iAtom;
            iNeighbor++;
         }
      }

      return result;
}

ssize_t SxSpaceSeparator::getNextPerGrad (ssize_t index, 
      SxVector3<Double> &corr)
{
   SxVector3<Int> pos = mesh.getMeshVec(index,SxMesh3D::Positive);
   SxVector3<Double> gradRho;
   gradRho(0) = gradRhoComponents(0)(index);
   gradRho(1) = gradRhoComponents(1)(index);
   gradRho(2) = gradRhoComponents(2)(index);

   Coord gradRhoRel = structure.cell.carToRel(gradRho);

   double factor = sqrt(gradRhoRel.absSqr().maxval ());
   if (factor > 1e-6) factor = 1.0/factor;
   else factor = 0;
   
   SxVector3<Double> rGrad = factor * gradRhoRel;
   SxVector3<Int> rGrid;
   for (int i = 0; i < 3; i++)  {
      int sign = rGrad(i) < 0 ? -1 : 1;
      if (fabs(rGrad(i)) >= 0.5) rGrid(i) = sign;
   }

   corr += rGrad - rGrid;
   for (int i = 0; i < 3; i++)  {
      int sign = corr(i) < 0 ? -1 : 1;
      if (fabs(corr(i)) >= 0.5) {
         rGrid(i) += sign;
         corr(i) -= 1.0 * sign;
      }
   }

   SxVector3<Int> next = pos + rGrid;

   ssize_t result = mesh.getMeshIdx(next,SxMesh3D::Unknown);

   if (rho(result) < rho(index)) result = index;
   if (fabs(rho(result) - rho(index)) < 1e-12) result = index;

   return result;
}

SxArray<SxDiracVec<Double> > SxSpaceSeparator::bader (
      const SxDiracVec<Double> &rhoIn)
{
   SX_CLOCK(Timer::Bader);
   
   rho = 1.0 * rhoIn;
   
   // compute Hessematrix 
   gradRhoComponents = computeGradRho ();
   gradGradRhoComponents = computeGradGradRho ();
  
   const SxRBasis &R = rho.getBasis<SxRBasis> ();

   int nDim = mesh.product ();
   int nAtoms = structure.getNAtoms ();
   SxArray<SxDiracVec<Double> > result(nAtoms);
   SxDiracVec<Int> regionMap (nDim);
   regionMap.set(-1);

   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      result(iAtom).resize(nDim);
      result(iAtom).set(0.0);
      result(iAtom).setBasis(R);
   }

   for (ssize_t idx = 0; idx < nDim; idx++) {
      if (fabs(rho(idx)) > 1e-4) {
         SxVector3<Double> corr = SxVector3<Double>(0., 0., 0.);
         regionMap(idx) = findMaximum(idx, regionMap, corr, true);
         if (regionMap(idx) >= 0) result(regionMap(idx))(idx) = 1.0;
      } else {
         SxArray<int> nextAtoms = getNextAtom(idx);
         for (int iNN = 0; iNN < nextAtoms.getSize(); iNN++)  {
            result(nextAtoms(iNN)) = 1.0 / double(nextAtoms.getSize ());
         }
         regionMap(idx) = nextAtoms(0);
      }
   }

   // final refinement
   for (ssize_t idx = 0; idx < nDim; idx++) {
         SxVector3<Double> corr = SxVector3<Double>(0., 0., 0.);
         regionMap(idx) = findMaximum(idx, regionMap, corr, false);
         if (regionMap(idx) >= 0) result(regionMap(idx))(idx) = 1.0;
   }

   SxString fileName = "baderRegions.sxb";
   SxBinIO io;
   io.open      (fileName, SxBinIO::BINARY_WRITE_ONLY);
   io.writeMesh (regionMap, structure.cell, mesh);
   io.setMode   (SxBinIO::WRITE_DATA);
   io.writeMesh (regionMap, structure.cell, mesh);
   io.close();

   return result;

}

void SxSpaceSeparator::printWeights (
      const SxArray<SxDiracVec<Double> > &weights, 
      ssize_t idx)
{ 
   cout << "Weights : ";
   double sum = 0.0;
   for (ssize_t i = 0; i < weights.getSize(); i++)  {
      sum += weights(i)(idx);
      cout << weights(i)(idx) << " ";
   }
   cout << sum << endl;
}

SxArray<SxDiracVec<Double> > SxSpaceSeparator::baderTrinkle (const SxDiracVec<Double> &rhoIn)
{
   SX_CLOCK(Timer::BaderTrinkle);

   rho = 1.0 * rhoIn;

   const SxRBasis &R = rho.getBasis<SxRBasis> ();

   int nAtoms = structure.getNAtoms ();
   int nDim = mesh.product ();

   SxArray<SxDiracVec<Double> > result (structure.getNAtoms ());
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      result(iAtom).resize (nDim);
      result(iAtom).set(0.0);
      result(iAtom).setBasis(R);
   }
   
   SxVector<Int> regionMap(nDim);
   regionMap.set(-1);

   //assign PAW regions
   int nPAWSphere = 0;
   double maxRho = rho.maxval ();
   for(ssize_t idx = 0; idx < rho.getSize (); idx++)  {
      Coord pos = structure.cell.relToCar (
            mesh.getMeshVec (idx,SxMesh3D::Positive)/Coord(mesh));
      SxArray<int> nextAtoms = getNextAtom(idx);
      if (nextAtoms.getSize () == 1)  {
         int iAtom = nextAtoms(0);
         int iSpecies = structure.getISpecies(iAtom);
         double dist = structure.cell.shortestDist(pos,structure(iAtom));
         if (dist < PAWDist(iSpecies)) {
            rho(idx) += maxRho;
            regionMap(idx) = iAtom;
            result(iAtom)(idx) = 1.0;
            nPAWSphere++;
         }
      }
   }

   //sort Density
   SxArray<ssize_t> sortIdx = rho.getSortIdx();
   SxVector<Int> interior (nAtoms);
   interior.set(0);
   int nBoundary = 0;

   for(ssize_t idx = sortIdx.getSize () - 1; idx >= 0; idx--)  {
      double rhoRef = rho(sortIdx(idx));
         Coord pos = structure.cell.relToCar (
               mesh.getMeshVec (sortIdx(idx),SxMesh3D::Positive)/Coord(mesh));
         //cout << SX_SEPARATOR;
         //cout << idx << " " << sortIdx(idx) << " " << rhoRef << endl;
         //cout << "Pos Abs: ";
         //pos.print ();
      // if not assigned
      if (regionMap(sortIdx(idx)) == -1) {
         // else if big enough density
         if (rhoRef > 1e-6)  {
            SxArray<ssize_t> neighbors = findNeighbors(sortIdx(idx), 3);
            enum pointType pType = Maximum;
            double nBiggerRho = 0.0;
            for (ssize_t iNeighbor = 0; 
                  iNeighbor < neighbors.getSize (); 
                  iNeighbor++)  {
               double rhoNeighbor = rho(neighbors(iNeighbor));
               if (rhoNeighbor > rhoRef) {
           //       cout << "Bigger neighbor: rho = " << rhoNeighbor
           //          << " belonging to atom " 
           //          << regionMap(neighbors(iNeighbor))
           //        << endl;
                  nBiggerRho += 1.0;
                  if (pType == Maximum)  {
                     if (regionMap(neighbors(iNeighbor)) >= 0)  {
                        regionMap(sortIdx(idx)) =
                           regionMap(neighbors(iNeighbor));
                        pType = Interior;
                     } else pType = Boundary;
                  }
                  if (pType == Interior)  { 
                     if (regionMap(neighbors(iNeighbor)) 
                           != regionMap(sortIdx(idx)))  
                        pType = Boundary;
                  }
               }
            }
            if (pType == Maximum) {
               cout << "Maximum outside PAW Sphere at" << endl;
               pos.print ();
               cout << "Density is " << rhoRef << endl;
               SX_LOOP(iNeighbor)
                  cout << rho(neighbors(iNeighbor)) << endl;
               SX_EXIT;
            } else if (pType == Interior) {
               result(regionMap(sortIdx(idx)))(sortIdx(idx)) = 1.0;
               //cout << "So point is INTERIOR belonging to atom " 
               //   << regionMap(sortIdx(idx)) << endl;
               interior(regionMap(sortIdx(idx))) += 1;
            } else if (pType == Boundary) {
               //cout << "So point is BOUNDARY" << endl;
               nBoundary += 1;
               double max = 0.55;
               int atom = -1;
               SxDiracVec<Double> pFlux 
                  = getProbabilityFlux(sortIdx(idx),rho);
               for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
                  result(iAtom)(sortIdx(idx)) = (pFlux * result(iAtom)).sum ();
                  if (result(iAtom)(sortIdx(idx)) > max) {
                     max = result(iAtom)(sortIdx(idx));
                     atom = iAtom;
                  }
               }
               regionMap(sortIdx(idx)) = atom;
            } else {
               cout << "So point is NO TYPE ???" << endl;
               SX_EXIT;
            } 
         } else  {
            // if small density VORONOI to avoid wobling   
            SxArray<int> nextAtoms = getNextAtom(sortIdx(idx));
            //cout << "So point is VORONOI belonging to atom " 
            //     << nextAtoms << " ";
            regionMap(sortIdx(idx)) = nextAtoms(0);
            for (int iNN = 0; iNN < nextAtoms.getSize(); iNN++)  {
               result(nextAtoms(iNN))(sortIdx(idx)) = 
                  1.0 / double(nextAtoms.getSize ());
            }
            //cout << endl;
         }
      }
      //printWeights(result, sortIdx(idx));
      //cout << SX_SEPARATOR;
   }
   cout << "Sphere points " << nPAWSphere << endl;
   cout << "INTERIOR ";
   interior.print();
   cout << "BOUNDARY " << nBoundary << endl;
   SxString fileName = "baderTrinkleRegions.sxb";
   SxBinIO io;
   io.open      (fileName, SxBinIO::BINARY_WRITE_ONLY);
   io.writeMesh (regionMap, structure.cell, mesh);
   io.setMode   (SxBinIO::WRITE_DATA);
   io.writeMesh (regionMap, structure.cell, mesh);
   io.close();

   return result;
}

SxArray<ssize_t> SxSpaceSeparator::findNeighbors(ssize_t index, int shell)
{
   SxStack<ssize_t> result;

   SxVector3<Int> pos = mesh.getMeshVec(index,SxMesh3D::Positive);

   SxStack<double> dists;
   SxStack<ssize_t> indices;

   int maxShift = shell+1;

   for (int i = -maxShift; i <= maxShift; i++)  {
      for (int j = -maxShift; j <= maxShift; j++)  {
         for (int k = -maxShift; k <= maxShift; k++)  {
            Coord rel = Coord(i,j,k)/Coord(mesh);
            Coord abs = structure.cell.relToCar (rel);
            double dist = abs.norm ();
            ssize_t next = 
               mesh.getMeshIdx(pos+SxVector3<Int>(i,j,k),SxMesh3D::Unknown);
            if (next != index)  {
               dists.push(dist);
               indices.push(next);
            }
         }
      }
   }

   SxVector<Double> distVec (dists);
   SxArray<ssize_t> indexVec (indices);

   
   for (int iShell = 0; iShell < shell; iShell++)  {
      double minDist = distVec.minval ();
      double maxDist = distVec.maxval ();
      for (ssize_t iNeighbor = 0; 
            iNeighbor < indexVec.getSize(); 
            iNeighbor++)  {
         double delta = distVec(iNeighbor) - minDist;
         if (fabs(delta) < 1e-12) {
            result.push(indexVec(iNeighbor));
            distVec(iNeighbor) += maxDist;
         }
      }
   }

   return SxArray<ssize_t> (result);  
}

double SxSpaceSeparator::ramp (double u)
{
   if (u > 0) return u;
   else return 0;
}

SxDiracVec<Double> SxSpaceSeparator::getProbabilityFlux (ssize_t index,
      const SxDiracVec<Double> &rhoIn)
{
   const SxCell& cell = structure.cell;

   Coord invDim = 1.0/Coord(mesh);

   Coord iPointRel = mesh.getMeshVec (index, SxMesh3D::Positive) * invDim;
   Coord iPointAbs = cell.relToCar(iPointRel);
   int shell = 3; //Next nearest neighbors
   SxArray<ssize_t> neighbors = findNeighbors(index,shell);

   SxDiracVec<Double> result (rhoIn.getSize ());
   result.set(0.0);

   for (ssize_t iNeighbor = 0; iNeighbor < neighbors.getSize(); iNeighbor++)  {
      ssize_t j = neighbors(iNeighbor);
      Coord jPointRel = mesh.getMeshVec (j, SxMesh3D::Positive) * invDim;
      Coord jPointAbs = cell.relToCar(jPointRel); 
      double dist = cell.shortestDist(jPointAbs,iPointAbs);
      //TODO Area has to be contact area of flux
      //double area = areas(getAreaIndex(distVec));
      double area = 1.0;
      double rhoDiff = rhoIn(j) - rhoIn(index);
      result(j) = area / dist * ramp(rhoDiff);
   }

   //normalization
   double norm = result.sum ();
   if (norm > 1e-6) result /= norm;

   return result;
}

int SxSpaceSeparator::getCriticalType (ssize_t index)
{
   int result;

   // build HesseMatrix
   SxMatrix<Double> Hesse(3,3);
   Hesse(0,0) = gradGradRhoComponents(0)(index); 
   Hesse(0,1) = gradGradRhoComponents(1)(index);
   Hesse(0,2) = gradGradRhoComponents(2)(index); 
   Hesse(1,0) = gradGradRhoComponents(3)(index); 
   Hesse(1,1) = gradGradRhoComponents(4)(index); 
   Hesse(1,2) = gradGradRhoComponents(5)(index); 
   Hesse(2,0) = gradGradRhoComponents(6)(index); 
   Hesse(2,1) = gradGradRhoComponents(7)(index); 
   Hesse(2,2) = gradGradRhoComponents(8)(index); 

   SxMatrix<Double>::Eigensystem eigen = Hesse.eigensystem ();

   //omega and signatur sigma
   SxVector<Double> eigenvalues = eigen.vals.real ();
   int sigma = 0;
   int omega = 0;

   for (int i = 0; i < 3; i++) {
      if (fabs(eigenvalues(i))> 1e-12) omega++;
      if (eigenvalues(i) > 0) sigma++;
      if (eigenvalues(i) < 0) sigma--;
   }

   // local minima
   if (omega == 3 && sigma == 3) result = 3;
   // ringcritical point
   if (omega == 3 && sigma == 1) result = 2;
   // bondcritical point
   if (omega == 3 && sigma == -1) result = 1;
   // local maxima
   if (omega == 3 && sigma == -3) result = 0;

   return result;
}

int SxSpaceSeparator::findMaximum (ssize_t index, SxDiracVec<Int> &currentMap, 
                                    SxVector3<Double> &corr, 
                                    bool setTrajectory)
{
   int result = -1;

   SxVector3<Int> pos = mesh.getMeshVec(index,SxMesh3D::Positive);

   ssize_t lastIdx = index;
   SxVector3<Int> next = SxVector3<Int>(pos);
   ssize_t nextIdx = mesh.getMeshIdx(next,SxMesh3D::Positive);

   bool killTrace = false;
   SxStack<ssize_t> idxTrace;
   idxTrace.push(index);
   ssize_t counter = 0;   

   do {
      counter++;
      nextIdx = getNextPerGrad(nextIdx,corr);
      if (currentMap(nextIdx) != -1) {
         // case: next knows its maximum
         result = currentMap(nextIdx);
         killTrace = true;
      } else if (lastIdx == nextIdx) {
         int cryticalType = getCriticalType(nextIdx);
         // maximum
         if (cryticalType == 0)  {
            SxArray<int> nextAtoms = getNextAtom(nextIdx);
            mesh.getMeshVec(nextIdx,SxMesh3D::Positive).print ();
            cout << nextAtoms << endl;
            // Maximum has to be UNIQUE Atom
            SX_CHECK(nextAtoms.getSize () == 1,nextAtoms.getSize ());
            result = nextAtoms(0);
            if(setTrajectory) {
               maximaIndices.push(nextIdx);
            }
         }
         // bondcritical point
         if (cryticalType == 1)  {
            SxArray<ssize_t> neighbors = findNeighbors(nextIdx, 1);
            for (ssize_t iNN = 0; iNN < neighbors.getSize (); iNN++)  {
               result = 
                  findMaximum (neighbors(iNN), currentMap, corr, setTrajectory);
            }
            if (setTrajectory) {
               bondIndices.push(nextIdx);
            }
         }
         // ringcritical point
         if (cryticalType == 2)  {
            SxArray<ssize_t> neighbors = findNeighbors(nextIdx, 1);
            for (ssize_t iNN = 0; iNN < neighbors.getSize (); iNN++)  {
               result = 
                  findMaximum (neighbors(iNN), currentMap, corr, setTrajectory);
            }
            if (setTrajectory) {
               ringIndices.push(nextIdx);
            }
         }
         // minimum
         if (cryticalType == 3)  {
            SxArray<ssize_t> neighbors = findNeighbors(nextIdx, 1);
            for (ssize_t iNN = 0; iNN < neighbors.getSize (); iNN++)  {
               result = 
                  findMaximum (neighbors(iNN), currentMap, corr, setTrajectory);
            }
            if (setTrajectory) {
               minimaIndices.push(nextIdx);
            }
         }
         currentMap(nextIdx) = result;
         killTrace = true;
      } else {
         // trace maximum further
         if (setTrajectory) idxTrace.push(nextIdx);
         lastIdx = nextIdx;
      }
   } while (killTrace == false && counter < mesh.product ());
   
   if (counter == mesh.product ()) {
      cout << "Find Maximum has not found one!" << endl;
      cout << SxArray<ssize_t> (idxTrace) << endl;
      SX_EXIT;
   }

   while (idxTrace.isEmpty() == false) {
      ssize_t idx = idxTrace.pop ();
      currentMap(idx) = result;
   }

   return result;
}

bool SxSpaceSeparator::isBorderPoint (ssize_t index, 
      SxDiracVec<Int> &currentMap)
{
   int region = currentMap(index); 
   
   SxVector3<Int> pos = mesh.getMeshVec(index,SxMesh3D::Positive);
   SxVector3<Int> shift = SxVector3<Int> (0,0,0);
   ssize_t nextIdx;
   
   for (int i = 0; i < 3; i++)  {
      shift(i) = -1;
      nextIdx = mesh.getMeshIdx(pos+shift,SxMesh3D::Unknown);
      if (currentMap(nextIdx) != region) return true;
      
      shift(i) = 1;
      nextIdx = mesh.getMeshIdx(pos+shift,SxMesh3D::Unknown);
      if (currentMap(nextIdx) != region) return true;
   }

   return false;
}

double SxSpaceSeparator::getShareFactor(ssize_t index, 
      const SxDiracVec<Int> &atomicSpace, int iAtom)
{

   SX_CHECK(index > -1);
   SX_CHECK(iAtom > -1);
   
   double result = 0.0;
   SxVector<Int> counter (structure.getNAtoms());
   counter.set(0);
   int shell = 0;

   while (counter.sum() == 0 && shell < 3) {
      shell++;
      SxArray<ssize_t> neighbors = findNeighbors(index,shell);

      for (int iNeighbor = 0; iNeighbor < neighbors.getSize (); iNeighbor++) {
         int region = atomicSpace(neighbors(iNeighbor));
         if (region > -1) counter(region)++;
      }

      if (counter.sum() > 0) 
         result = double(counter(iAtom)) / double(counter.sum ());
   }

   return result;
            
}

SxDiracVec<Double> SxSpaceSeparator::getFilter (const SxDiracVec<Int> atomicSpace, int iAtom)
{
   const SxRBasis &R = atomicSpace.getBasis<SxRBasis> ();
   ssize_t dim = atomicSpace.getSize ();
   SxDiracVec<Double> result(dim);
   result.set(0.0);

   for (ssize_t i = 0; i < dim; i++)  {
      if (atomicSpace(i) == iAtom) result(i) = 1.0;
      if (atomicSpace(i) < 0 && iAtom > -1)  {
         result(i) = getShareFactor(i,atomicSpace,iAtom);
      }
   }
   
   result.setBasis(&R);

   return result;
}

void SxSpaceSeparator::plotDensity(const SxString &filename)
{
   RhoR density(1);
   density(0) = 1.0 * rho;
   SxString file = filename + ".sxb";
   SxRho(density).writeRho(file);
}

void SxSpaceSeparator::plotGradDensity(const SxString &filename)
{
   for (int i = 0; i < 3; i++)  {
      RhoR density(1);
      density(0) = 1.0 * gradRhoComponents(i);
      SxString file = filename + "-" + SxString(i) + ".sxb";
      SxRho(density).writeRho(file);
   }
}
