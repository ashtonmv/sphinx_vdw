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

#include <SxAOMatTB.h>

#ifndef SX_STANDALONE

SxAOMatTB::SxAOMatTB (const SxSymbolTable *tablePtr)
{
   SX_START_TIMER(InitTimer);
   if (tablePtr->containsGroup("aomatTB"))
      tablePtr = tablePtr->getGroup("aomatTB");
   else  {
      cout << "Input File does not contain aomatTB group!" << endl;
      SX_QUIT;
   }
   structure = SxAtomicStructure(tablePtr);
   speciesData.readSpecies (tablePtr);
   setKPoint (Coord(0.0,0.0,0.0));

   // --- symmetry
   try  {
      if (tablePtr->containsGroup ("symmetry"))  {
         const SxSymbolTable *sym, *symmetry;
         SxList<SymMat> symList;
         symmetry = tablePtr->getGroup ("symmetry");
         if (symmetry->containsGroup ("operator"))  {
            for ( sym  = symmetry->getGroup ("operator");
                  sym != NULL;
                  sym  = sym->nextSibling ("operator"))
            {
               symList << SymMat(sym->get("S")->toList());
            }
         }
         syms = symList;
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   cout << syms.getSize() << " Symmetries found!" << endl;

   SxGrid grid(structure, 10);
   SxNeighbors nn;
   int mode = SxNeighbors::StoreIdx 
            | SxNeighbors::StoreRel
            | SxNeighbors::IncludeZeroDistance;

   SxArray<SxStack<int> > neighborStack(structure.getNAtoms ());
   SxStack<SxMatrix<Double> > hamStack;
   SxStack<SxMatrix<Double> > overlapStack;
   SxArray<SxStack<int> > idxStack(structure.getNAtoms ());
   SxArray<SxStack<int> > symIdStack (structure.getNAtoms ());
   SxVector<Int> oriId;
   orientations.resize(structure.getNAtoms ());
   SxArray<SxStack<int> > neighborIdStack (structure.getNAtoms ());
   SxStack<Coord> genNeighborStack;
   SxArray<SxStack<Coord> > neighborCoordStack(structure.getNAtoms ());
   orbitalTypeList.resize(structure.getNSpecies());
   int idx = 0;

   //try  {
      const SxSymbolTable *species, neigborGroup;
      int is;
      int iGenNeighbor = 0;
      // --- read species data
      for (species  = tablePtr->getGroup("species"), is=0;
           species != NULL;
           species  = species->nextSibling("species"), is++)
      {
         
         // get OrbitalTypeInfo
         SxArray<int> OrbPerL (species->get("orbitalsPerL")->toIntList ());
         int lMax = (int)OrbPerL.getSize() - 1;
         orbitalTypeList(is).resize (0);
         for (int l = 0; l <= lMax; l++)  {
            for (int iOrbitalType = 0; iOrbitalType < OrbPerL(l); iOrbitalType++)  {
               orbitalTypeList(is).append(l);
               cout << l << " ";
            }
         }
         cout << endl;
      }
      for (species  = tablePtr->getGroup("species"), is=0;
            species != NULL;
            species  = species->nextSibling("species"), is++)
      {
         SxArray<int> OrbPerL (species->get("orbitalsPerL")->toIntList ());
         int lMax = (int)OrbPerL.getSize() - 1;
         
         // Is MatrixFile available ?
         SxBinIO matIO;
         int fileSpecies = -1;
         bool fileOpen = false;
         if (species->containsGroup("matFile"))  {
            const SxSymbolTable *matFile = species->getGroup("matFile");
            fileSpecies = matFile->get("iSpecies")->toInt ();
            SxString fileName = matFile->get("file")->toString ();
            try {
               matIO.open(fileName,SxBinIO::BINARY_READ_ONLY);
            } catch (SxException e)  {
               e.print ();
               SX_EXIT;
            }
            fileOpen = true;
         }

         if (species->containsGroup ("onside"))  {
            const SxSymbolTable *onsideGroup = species->getGroup ("onside");
            if (onsideGroup->contains ("orientations"))  {
               oriId = SxVector<Int> (onsideGroup->get ("orientations")
                     ->toIntList ());
               oriId -= 1; // we start from 0, but in input from 1
            }
            SxVector<Int> symId(1);
            symId.set(-1);
            if (onsideGroup->contains ("symmetries"))  {
               symId = SxVector<Int> (onsideGroup->get ("symmetries")
                     ->toIntList ());
               symId -= 1; // we start from 0, but in input from 1
            }
            SxMatrix3<Double> rotRef;
            SxMatrix<Double> overlapRef;
            SxMatrix<Double> hamRef;

            if (onsideGroup->contains ("R") && 
                onsideGroup->contains ("S") && 
                onsideGroup->contains ("H"))  {
               // read Rotation matrix
               rotRef =  SxMatrix3<Double>(onsideGroup->get("R")->toList ());
               // read reference Overlap matrix
               overlapRef = SxMatrix<Double>(onsideGroup->get("S")->toList ());
               // read reference Hamiltonian matrix
               hamRef = SxMatrix<Double>(onsideGroup->get("H")->toList ());
               int nElements = (int)overlapRef.getSize();
               int nRows = int(lround(sqrt((double)nElements)));
               if (nRows * nRows != nElements)  {
                  SX_EXIT;
               }
               overlapRef.reshape(nRows, nRows);
               hamRef.reshape(nRows, nRows);
               overlapRef = overlapRef.transpose();
               hamRef = hamRef.transpose();
               cout << SX_SEPARATOR;
               cout << "H onside: " << endl;
               for (int iRow = 0; iRow < hamRef.handle->nRows; iRow++)  {
                  for (int iCol = 0; iCol < hamRef.handle->nCols; iCol++)  {
                     sxprintf ("% .4f ", hamRef(iRow,iCol)); 
                  }
                  cout << endl;
               }
               cout << endl;
               cout << SX_SEPARATOR;
            } else if (fileOpen)  {
               Coord dist = Coord(0.,0.,0.);
               int iNeighbor = findNeighbor(matIO, fileSpecies, fileSpecies, dist);
               rotRef     = readR (matIO, iNeighbor);
               overlapRef = readS (matIO, iNeighbor);
               hamRef     = readH (matIO, iNeighbor);
               cout << SX_SEPARATOR;
               cout << iNeighbor << endl;
               cout << "H onside: " << endl;
               for (int iRow = 0; iRow < hamRef.handle->nRows; iRow++)  {
                  for (int iCol = 0; iCol < hamRef.handle->nCols; iCol++)  {
                     sxprintf ("% .4f ", hamRef(iRow,iCol)); 
                  }
                  cout << endl;
               }
               cout << endl;
               cout << SX_SEPARATOR;
            } else  {
               cout << "Either file or matElements have to be specified for onside on Species " << is << endl;
               SX_QUIT;
            }


            if (getNOrbs(is) != overlapRef.handle->nRows)  {
               cout << "Dimension Error!" << endl;
               cout << getNOrbs(is) << " orbitals expected due to orbitals per L, but Overlap has only" << endl;
               cout << overlapRef.handle->nRows << " rows in species "<< is << "." << endl;
               SX_QUIT;
            }
            if (getNOrbs(is) != overlapRef.handle->nRows)  {
               cout << "Dimension Error!" << endl;
               cout << getNOrbs(is) << " orbitals expected due to orbitals per L, but Hamiltonian has only" 
                    << endl;
               cout << overlapRef.handle->nRows << " rows in species "<< is << "." << endl;
               SX_QUIT;
            }
            // read reference neighbor position
            Coord posRef (onsideGroup->get("coords")->toList ());
            int isRef = -1;
            if (onsideGroup->contains ("atomType"))  {
               isRef = speciesData.find (onsideGroup->get("atomType")->toString ());
            }
            genNeighborStack << posRef;
            // --- find all corresponding neighbors
            SxYlm::SxClebschTable complexClebschTable = 
               SxYlm::getClebschGordan (lMax, lMax, lMax, SxYlm::ComplexYlm);
            for (int iOri = oriId.getSize () > 0 ? 0 : -1;
                  iOri < oriId.getSize ();
                  ++iOri)
            {
               SymMat O = iOri >= 0 
                  ? syms(oriId(iOri))
                  : SymMat (1,0,0, 0,1,0, 0,0,1);
               // Build Matrix for Overlap and Hamiltonian Rotation
               // 8 atom cell Si needs the -1.0 no clue why...
               SxMatrix3<Double> fullRot = 1.0 * (O ^ rotRef);
               SxArray<SxMatrix<Double> > Dmm 
                  = SxYlm::computeYlmRotMatrices (lMax, fullRot, complexClebschTable);
               int iot = 0;
               SxMatrix<Double> ham = hamRef.getCopy();
               SxMatrix<Double> overlap = overlapRef.getCopy();
               for (int iRow = 0; iRow < hamRef.handle->nRows; /*inside*/)  {
                  int rowL = orbitalTypeList(is)(iot);
                  int addRow = 2 * rowL + 1;
                  int jot = 0;
                  for (int iCol = 0; iCol < hamRef.handle->nCols; /*inside*/)  {
                     int colL = orbitalTypeList(is)(jot);
                     int addCol = 2 * colL + 1;
                     SxMatrix<Double> block, blockRot;
                     // Rotate ham
                     block = ham.getBlock(iRow,iCol,addRow,addCol);
                     blockRot = (Dmm(rowL) ^ block ^ Dmm(colL).transpose());
                     SX_CHECK ((block - (Dmm(rowL).transpose() ^ blockRot ^ Dmm(colL))).abs ().sum () < 1e-6,
                               (block - (Dmm(rowL).transpose() ^ blockRot ^ Dmm(colL))).abs ().sum ());
                     ham.setBlock(blockRot, iRow, iCol);

                     // Rotate overlap
                     block = overlap.getBlock(iRow,iCol,addRow,addCol);
                     blockRot = (Dmm(rowL) ^ block ^ Dmm(colL).transpose());
                     SX_CHECK ((block - (Dmm(rowL).transpose() ^ blockRot ^ Dmm(colL))).absSqr ().sum () < 1e-6,
                               (block - (Dmm(rowL).transpose() ^ blockRot ^ Dmm(colL))).absSqr ().sum ());
                     overlap.setBlock(blockRot, iRow, iCol);

                     iCol += addCol;
                     jot++;
                  }
                  iRow += addRow;
                  iot++;
               }
               
               hamStack << ham;
               overlapStack << overlap;


               for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {
                  bool onsideFound = false;
                  int iTlAtom = structure.getIAtom (is, ia);
                  int nSyms = (int)symId.getSize();
                  int iSym = 0;
                  while (!onsideFound && iSym < nSyms)  {
                     SymMat S = symId(iSym) >= 0 ? syms(symId(iSym)) : SymMat (1,0,0, 0,1,0, 0,0,1);
                     Coord nPos = O ^ S ^ posRef;
                     nn.compute (grid, structure, structure(is,ia) + nPos,
                           1e-1, mode);
                     if (nn.getSize () == 1)  {
                        // found a neighbor: check species
                        int jTl = nn.idx(0);
                        if (isRef < 0 || structure.getISpecies(jTl) == isRef) {
                           onsideFound = true;
                           // found a neighbor with correct species
                           neighborStack(iTlAtom) << iTlAtom; // ONSIDE !
                           neighborCoordStack(iTlAtom) << Coord (0.0, 0.0, 0.0); // ONSIDE !
                           orientations(iTlAtom) = oriId(iOri);
                           // keep track of neighbor origin
                           neighborIdStack(iTlAtom) << iGenNeighbor;
                           symIdStack(iTlAtom) << (iOri>=0 ? oriId(iOri) : -1);

                           // idx for hamiltonians and overlap
                           idxStack(iTlAtom) << idx;
                        }
                     } else if (nn.getSize () > 1)  {
                        cout << "Ambiguous neighbors for atom "
                           << speciesData.chemName(is) << (ia+1) << "@ "
                           << structure(is, ia)
                           << " near "
                           << (structure(is, ia) + nPos)
                           << endl;
                        cout << "Candidates are " << endl;
                        for (int in = 0; in < nn.getSize (); ++in)  {
                           cout << "atom " << (nn.idx(in)+1) << endl;
                        }
                        SX_EXIT;
                     }
                     iSym++;
                  }
               }
               idx++;
            }
         }

         if (species->containsGroup ("neighbor"))  {
            const SxSymbolTable *nGroup = species->getGroup ("neighbor");
            for (;
                  nGroup != NULL;
                  nGroup = nGroup->nextSibling ("neighbor"), iGenNeighbor++)
            {
               SxVector<Int> symId;
               if (nGroup->contains ("symmetries"))  {
                  symId = SxVector<Int> (nGroup->get ("symmetries")
                        ->toIntList ());
                  symId -= 1; // we start from 0, but in input from 1
               }
               
               // read reference neighbor position
               Coord posRef (nGroup->get("coords")->toList ());
               int isRef = -1;
               int neighborSpecies = (nGroup->get("species")->toInt ()) - 1;
               if (nGroup->contains ("atomType"))  {
                  isRef = speciesData.find (nGroup->get("atomType")->toString ());
               }
               genNeighborStack << posRef;

               SxMatrix3<Double> rotRef;
               SxMatrix<Double> overlapRef;
               SxMatrix<Double> hamRef;
               if (nGroup->contains ("R") && 
                   nGroup->contains ("S") && 
                   nGroup->contains ("H"))  { 
                  // read Rotation matrix
                  rotRef = SxMatrix3<Double>(nGroup->get("R")->toList ());
                  // read reference Overlap matrix
                  overlapRef = SxMatrix<Double>(nGroup->get("S")->toList ());
                  // read reference Hamiltonian matrix
                  hamRef = SxMatrix<Double>(nGroup->get("H")->toList ());
                  int nElements = (int)overlapRef.getSize();
                  int nRows = (int)sqrt((double)nElements);
                  overlapRef.reshape(nRows, nRows);
                  hamRef.reshape(nRows, nRows);
                  overlapRef = overlapRef.transpose();
                  hamRef = hamRef.transpose();
               } else if (fileOpen)  {
                  // for time beeing jSpecies = filespecies. Should be speciefied in inputFile
                  int iNeighbor = findNeighbor(matIO, fileSpecies, fileSpecies, posRef);
                  cout << "iNeighbor " << iNeighbor << endl;
                  rotRef     = readR (matIO, iNeighbor);
                  overlapRef = readS (matIO, iNeighbor);
                  hamRef     = readH (matIO, iNeighbor);
               } else  {
                  cout << "Either file or matElements have to be specified for onside on Species " << is << endl;
                  SX_QUIT;
               }


               if (getNOrbs(is) != overlapRef.handle->nRows)  {
                  cout << "Dimension Error!" << endl;
                  cout << getNOrbs(is) << " orbitals expected due to orbitals per L, but Overlap has only" << endl;
                  cout << overlapRef.handle->nRows << " rows in species "<< is << "." << endl;
                  SX_QUIT;
               }
               if (getNOrbs(is) != overlapRef.handle->nRows)  {
                  cout << "Dimension Error!" << endl;
                  cout << getNOrbs(is) << " orbitals expected due to orbitals per L, but Hamiltonian has only" << endl;
                  cout << overlapRef.handle->nRows << " rows in species "<< is << "." << endl;
                  SX_QUIT;
               }
               
               // --- find all corresponding neighbors
               SxYlm::SxClebschTable complexClebschTable = SxYlm::getClebschGordan (lMax, lMax, lMax, SxYlm::ComplexYlm);
               for (int iOri = oriId.getSize () > 0 ? 0 : -1;
                           iOri < oriId.getSize ();
                           iOri++)
               {
                  SymMat O = iOri >= 0 
                     ? syms(oriId(iOri))
                           : SymMat (1,0,0, 0,1,0, 0,0,1);
                  for (int iSym = symId.getSize () > 0 ? 0 : -1;
                        iSym < symId.getSize ();
                        ++iSym)
                  {
                     SymMat S = iSym >= 0 
                        ? syms(symId(iSym))
                        : SymMat (1,0,0, 0,1,0, 0,0,1);

                     // Build Matrix for Overlap and Hamiltonian Rotation
                     SxMatrix3<Double> fullRot = S ^ O ^ rotRef;
                     SxArray<SxMatrix<Double> > Dmm 
                        = SxYlm::computeYlmRotMatrices (lMax, fullRot, complexClebschTable);
                     int iot = 0;
                     SxMatrix<Double> ham = hamRef.getCopy();
                     SxMatrix<Double> overlap = overlapRef.getCopy();
                     for (int iRow = 0; iRow < hamRef.handle->nRows; /*inside*/)  {
                        int rowL = orbitalTypeList(is)(iot);
                        int addRow = 2 * rowL + 1;
                        int jot = 0;
                        for (int iCol = 0; iCol < hamRef.handle->nCols; /*inside*/)  {
                           int colL = orbitalTypeList(neighborSpecies)(jot);
                           int addCol = 2 * colL + 1;
                           SxMatrix<Double> block, blockRot;
                           // Rotate ham
                           block = ham.getBlock(iRow,iCol,addRow,addCol);
                           blockRot = (Dmm(rowL) ^ block ^ Dmm(colL).transpose());
                           ham.setBlock(blockRot, iRow, iCol);

                           // Rotate overlap
                           block = overlap.getBlock(iRow,iCol,addRow,addCol);
                           blockRot = (Dmm(rowL) ^ block ^ Dmm(colL).transpose());
                           overlap.setBlock(blockRot, iRow, iCol);

                           iCol += addCol;
                           jot++;
                        }
                        iRow += addRow;
                        iot++;
                     }

                     Coord pos = S ^ O ^ posRef;

                     hamStack << ham;
                     overlapStack << overlap;


                     for (int ia = 0; ia < structure.getNAtoms (is); ++ia)  {
                        int iTlAtom = structure.getIAtom (is, ia);
                        nn.compute (grid, structure, structure(is,ia) + pos,
                              1e-4, mode);
                        if (nn.getSize () == 1)  {
                           // found a neighbor: check species
                           int jTl = nn.idx(0);
                           if (isRef < 0 || structure.getISpecies(jTl) == isRef) {
                              // found a neighbor with correct species, check for orientation
                              if (oriId(iOri) == orientations(iTlAtom))  {
                                 neighborStack(iTlAtom) << jTl;
                                 neighborCoordStack(iTlAtom) << pos;
                                 // keep track of neighbor origin
                                 neighborIdStack(iTlAtom) << iGenNeighbor;
                                 symIdStack(iTlAtom) << (iSym>=0 ? symId(iSym) : -1);

                                 // idx for hamiltonians and overlap
                                 idxStack(iTlAtom) << idx;
                              }
                           }
                        } else if (nn.getSize () > 1)  {
                           cout << "Ambiguous neighbors for atom "
                              << speciesData.chemName(is) << (ia+1) << "@ "
                              << structure(is, ia)
                              << " near "
                              << (structure(is, ia) + pos)
                              << endl;
                           cout << "Candidates are " << endl;
                           for (int in = 0; in < nn.getSize (); ++in)  {
                              cout << "atom " << (nn.idx(in)+1) << endl;
                           }
                           SX_EXIT;
                        }
                     }
                     idx++;
                  }
               }
            }
         }
         matIO.close ();
      }
   /*
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   */
   genNeighbors = genNeighborStack;

   genSymId.resize (structure.getNAtoms ());
   genNeighborId.resize (structure.getNAtoms ());
   neighborIdx.resize (structure.getNAtoms ());
   hamiltonians = hamStack;
   overlaps = overlapStack;;
   neighborCoord.resize (structure.getNAtoms ());
   idxList.resize(structure.getNAtoms ());
   for (int iTlAtom = 0; iTlAtom < structure.getNAtoms (); ++iTlAtom)  {
      neighborIdx(iTlAtom) = neighborStack(iTlAtom);
      genNeighborId(iTlAtom) = neighborIdStack(iTlAtom);
      genSymId(iTlAtom) = symIdStack(iTlAtom);
      neighborCoord(iTlAtom) = neighborCoordStack(iTlAtom);
      idxList(iTlAtom) = idxStack(iTlAtom);
      cout << "atom " << (iTlAtom) << ": " 
           << neighborIdx(iTlAtom).getSize () << " neighbors ->"
           << neighborIdx(iTlAtom)
           << endl;
   }
   for (int iTlAtom = 0; iTlAtom < structure.getNAtoms (); ++iTlAtom)  {
      neighborIdx(iTlAtom) = neighborStack(iTlAtom);
      genNeighborId(iTlAtom) = neighborIdStack(iTlAtom);
      genSymId(iTlAtom) = symIdStack(iTlAtom);
      neighborCoord(iTlAtom) = neighborCoordStack(iTlAtom);
      idxList(iTlAtom) = idxStack(iTlAtom);
      cout << "atom " << (iTlAtom) << ": " 
           << idxList(iTlAtom).getSize () << " hams ->"
           << idxList(iTlAtom)
           << endl;
   }
   SX_STOP_TIMER(InitTimer);
}

int SxAOMatTB::findNeighbor (SxBinIO &io, int iSpecies, int jSpecies, Coord dist)
{
   int coordDim = 3;
   int nNeighbors = io.getDimension("nNeighbors");
   SxMatrix<Double> distFile (nNeighbors,coordDim);
   SxVector<Int> neighborSpecies (nNeighbors);
   io.read("distances", &distFile, nNeighbors, coordDim);
   io.read("neighborSpecies",&neighborSpecies, nNeighbors);
   for (int iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++)  {
      Coord d;
      d(0) = distFile.row(iNeighbor)(0); 
      d(1) = distFile.row(iNeighbor)(1); 
      d(2) = distFile.row(iNeighbor)(2); 
      if ( (neighborSpecies(iNeighbor) == jSpecies) && ((dist - d).norm() < 1e-4) ) {
         return iNeighbor;
      }
   }

   // When arriving here no item has been found... Very bad!
   cout << "Matrix not found" << endl;
   cout << "iSpecies = " << iSpecies << endl; 
   cout << "coord = " << dist << endl;
   cout << "d=" << dist.norm () << endl;
   SX_QUIT;
}

SxMatrix3<Double> SxAOMatTB::readR (SxBinIO &io, int iNeighbor)
{
   SxMatrix3<Double> result;
   SxMatrix<Double> buffer (3,3);
   int offset = 3 * iNeighbor;
   io.read("R-Matrices", &buffer, 3, 3, offset, 0);
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
         result(i,j) = buffer(i,j);

   return result;
}

SxMatrix<Double> SxAOMatTB::readS (SxBinIO &io, int iNeighbor)
{
   int nNeighbors = io.getDimension("nNeighbors");
   SxVector<Int> nNeighborOrbs (nNeighbors);
   io.read("nNeighborOrbs", &nNeighborOrbs, nNeighbors);
   int nRows = nNeighborOrbs(iNeighbor);
   int nCols = io.getDimension("nOrbs");
   SxMatrix<Double> result (nRows, nCols);
   int offset = 0;
   for(int in = 0; in < iNeighbor; in++) offset+=nNeighborOrbs(in);
   io.read("S-Matrices", &result, nRows, nCols, offset, 0);
   return result;
}

SxMatrix<Double> SxAOMatTB::readH (SxBinIO &io, int iNeighbor)
{
   int nNeighbors = io.getDimension("nNeighbors");
   SxVector<Int> nNeighborOrbs (nNeighbors);
   io.read("nNeighborOrbs", &nNeighborOrbs, nNeighbors);
   int nRows = nNeighborOrbs(iNeighbor);
   int nCols = io.getDimension("nOrbs");
   SxMatrix<Double> result (nRows, nCols);
   int offset = 0;
   for(int in = 0; in < iNeighbor; in++) offset+=nNeighborOrbs(in);
   io.read("H-Matrices", &result, nRows, nCols, offset, 0);
   
   return result;
}

void SxAOMatTB::info() const
{
   int nAtoms = (int)neighborIdx.getSize();

   // loop over all atoms
   for (int iAtom = 0; iAtom < nAtoms; iAtom++)  {
      int nNeighbors = (int)neighborIdx(iAtom).getSize();
      // loop over all neighbors
      for (int iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++)  {
         // print coord, distance, overlap, and hamiltonian
         cout << structure(iAtom) << "<=>" << (structure(iAtom)+ neighborCoord(iAtom)(iNeighbor))
              << "; d=" <<  neighborCoord(iAtom)(iNeighbor).norm () << endl << endl;
         cout << "coords = "; neighborCoord(iAtom)(iNeighbor).print();
         int idx = idxList(iAtom)(iNeighbor);
         cout << "S : " << endl; 
         for (int iRow = 0; iRow < overlaps(idx).handle->nRows; iRow++)  {
            for (int iCol = 0; iCol < overlaps(idx).handle->nCols; iCol++)  {
               sxprintf ("% .4f ", overlaps(idx)(iRow,iCol)); 
            }
            cout << endl;
         }
         cout << endl;
         cout << "H : " << endl; 
         for (int iRow = 0; iRow < hamiltonians(idx).handle->nRows; iRow++)  {
            for (int iCol = 0; iCol < hamiltonians(idx).handle->nCols; iCol++)  {
               sxprintf ("% .4f ", hamiltonians(idx)(iRow,iCol)); 
            }
            cout << endl;
         }
      }
   }
   // build and print Overlap and Ham
   int dim = getNOrbs ();
   SxVector<Complex16> psi (dim);
   SxMatrix<Complex16> H (dim,dim);
   SxMatrix<Complex16> S (dim,dim);
   for (int i = 0; i < dim; i++)  {
      psi.set(0.0);
      psi(i) = 1.0;
      H.colRef(i) << applyH(psi);
      S.colRef(i) << applyS(psi);
   }
   cout << SX_SEPARATOR;
   cout << "ik = 0" << endl;
   cout << "S: " << endl;
   S.eigenvalues ().real ().print(true);
   cout << "H: " << endl;
   H.eigenvalues ().real ().print(true);
   cout << SX_SEPARATOR;
}

int SxAOMatTB::getOffset (int iSpecies) const
{
   int result = 0;
   for (int is = 0; is < iSpecies; is++)  {
      result+= structure.getNAtoms(is) * getNOrbs(is);
   }

   return result;
}

int SxAOMatTB::getNOrbs () const
{
   int result = 0;
   for (int iSpecies = 0; iSpecies < structure.getNSpecies(); iSpecies++)  {
      result += getNOrbs(iSpecies) * structure.getNAtoms(iSpecies);
   }

   return result;
}


int SxAOMatTB::getNOrbs (int iSpecies) const
{
   SX_CHECK(orbitalTypeList.getSize() > 0);
   SX_CHECK(orbitalTypeList(iSpecies).getSize() > 0);

   int result = 0;
   for (int iot = 0; iot < orbitalTypeList(iSpecies).getSize(); iot++)  {
      result += 2 * orbitalTypeList(iSpecies)(iot) + 1;
   }

   return result;
}

SxVector<Complex16> SxAOMatTB::applyH (const SxVector<Complex16> &psi) const
{
   SX_START_TIMER(HPsiTime);
   SxVector<Complex16> result(psi.getSize());
   result.set(0.0);

   for (int iSpecies = 0; iSpecies < structure.getNSpecies(); iSpecies++) {
      for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++) {
         int iTlAtom = structure.getIAtom (iSpecies,iAtom);
         for (int iNeighbor = 0; iNeighbor < neighborIdx(iTlAtom).getSize (); iNeighbor++)  {
            double phase = kPoint ^ neighborCoord(iTlAtom)(iNeighbor);
            SxMatrix<Complex16> ham = hamiltonians(idxList(iTlAtom)(iNeighbor)) * exp(I*phase);
            int length = getNOrbs(iSpecies); // how many orbitals belong to the species of iTlAtom?
            int offset = getOffset(iSpecies) + iAtom * length; // which part of the vector belongs to this iTlAtom
            int neighborSpecies = structure.getISpecies(neighborIdx(iTlAtom)(iNeighbor));
            int neighborLength = getNOrbs(neighborSpecies);
            int neighborIAtom = neighborIdx(iTlAtom)(iNeighbor) - structure.atomInfo->offset(neighborSpecies);
            int neighborOffset = getOffset(neighborSpecies) + neighborIAtom * length;
            
            //Averaging for keeping H hermitian
            // apply H | result(iTlAtom) += H_{iTlAtom,iNeighbor} ^ c_iNeighbor
            result(SxIdx(offset,offset+length-1)) += 
               ham ^ psi(SxIdx(neighborOffset,neighborOffset+neighborLength-1));
               
            // apply H | result(iNeighbor) += H_{iTlAtom,iNeighbor}.transpose() ^ c_iTlAtom
            result(SxIdx(neighborOffset,neighborOffset+neighborLength-1)) += 
               ham.adjoint() ^ psi(SxIdx(offset,offset+length-1));
         }
      }
   }
   SX_STOP_TIMER(HPsiTime);

   return 0.5 * result;
}

SxVector<Complex16> SxAOMatTB::applyS (const SxVector<Complex16> &psi) const
{
   SX_START_TIMER(SPsiTime);
   SxVector<Complex16> result(psi.getSize());
   result.set(0.0);

   for (int iSpecies = 0; iSpecies < structure.getNSpecies(); iSpecies++) {
      for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++) {
         int iTlAtom = structure.getIAtom (iSpecies,iAtom);
         for (int iNeighbor = 0; iNeighbor < neighborIdx(iTlAtom).getSize (); iNeighbor++)  {
            double phase = kPoint ^ neighborCoord(iTlAtom)(iNeighbor);
            SxMatrix<Complex16> overlap = overlaps(idxList(iTlAtom)(iNeighbor)) * exp(I*phase);
            int length = getNOrbs(iSpecies); // how many orbitals belong to the species of iTlAtom?
            int offset = getOffset(iSpecies) + iAtom * length; // which part of the vector belongs to this iTlAtom
            int neighborSpecies = structure.getISpecies(neighborIdx(iTlAtom)(iNeighbor));
            int neighborLength = getNOrbs(neighborSpecies);
            int neighborIAtom = neighborIdx(iTlAtom)(iNeighbor) - structure.atomInfo->offset(neighborSpecies);
            int neighborOffset = getOffset(neighborSpecies) + neighborIAtom * length;
            
            // Avereraging for keeping S hermitian
            // apply S | result(iTlAtom) += H_{iTlAtom,iNeighbor} ^ c_iNeighbor
            result(SxIdx(offset,offset+length-1)) += 
               overlap ^ psi(SxIdx(neighborOffset,neighborOffset+neighborLength-1));
               
            // apply S | result(iNeighbor) += H_{iTlAtom,iNeighbor}.adjoint() ^ c_iTlAtom
               result(SxIdx(neighborOffset,neighborOffset+neighborLength-1)) += 
                  overlap.adjoint () ^ psi(SxIdx(offset,offset+length-1));
         }
      }
   }
   SX_STOP_TIMER(SPsiTime);

   return 0.5 * result;
}

void SxAOMatTB::directMin (const SxKPoints &kPoints)
{
   int dim = getNOrbs();
   SxVector<Complex16> psi (dim);

   // --- Perform diagonalisation
   SxMatrix<Double> eigVals (dim, kPoints.getNk ());
   for (int ik = 0; ik < kPoints.getNk (); ik++)  {
      // Build full H and S
      setKPoint(kPoints.getK(ik));
      SxMatrix<Complex16> H (dim,dim);
      SxMatrix<Complex16> S (dim,dim);
      for (int i = 0; i < dim; i++)  {
         psi.set(0.0);
         psi(i) = 1.0;
         H.colRef(i) << applyH(psi);
         S.colRef(i) << applyS(psi);
      }

      // Cholesky decomposition as orthogonalisation
      SxMatrix<Complex16>::Eigensystem seig = S.eigensystem ();
      if (seig.vals (0).re < 1e-12)  {
         cout << SX_SEPARATOR;
         cout << "WARNING! Remove negative eigenvalue in Overlap matrix!" << endl;
         cout << SX_SEPARATOR;
         int i = 0;
         while (seig.vals(i).re < 1e-12) i++;
         SxVector<Complex16> nonSing
            = seig.vecs (SxIdx(dim * i, dim*dim -1));
         nonSing.reshape (dim, dim-i);
         S = nonSing.adjoint () ^ S ^ nonSing;
         H = nonSing.adjoint () ^ H ^ nonSing;
      }
      SxMatrix<Complex16> LInv = S.choleskyDecomposition ().inverse ();
      // transform H in orthogonal basis
      H = LInv ^ H ^ LInv.adjoint ();
      // Solve standart eigenvalue problem
      SxVector<Double> vals = H.eigenvalues ().real () * HA2EV;
      if (vals.getSize() == dim)
         eigVals.colRef(ik) << vals;
      else  {
         for (int iRow = 0; iRow < dim; iRow++)  {
            if (iRow < vals.getSize()) eigVals.colRef(ik)(iRow) = vals(iRow);
            else eigVals.colRef(ik)(iRow) = 0;
         }
      }
   }
   eigVals = eigVals.transpose();
   SxBinIO out; 
   SxString file = "eps-aoMatTB.dat";
   out.open(file, SxBinIO::ASCII_WRITE_ONLY);
   out.writeNXYPlot(eigVals);
   out.close();
}


#else // SX_STANDALONE


int main (int argc, char** argv)
{
   initSPHInXMath ();

   cout.precision(10);

   // Command line parsing
   SxCLI cli (argc,argv);

   // Define Author
   cli.authors = "B. Lange";

   // What does the program do ?
   cli.preUsageMessage = "AAOMAT Tight-Binding";

   SxString inputFile = cli.option ("-i|--input", "file", "SPHInX input file")
                                   .toString ("input.sx");
   SxVector<Double> eAim;
   cli.option ("--eAim", "List of double Values", "Aim Energy in eV");
   cli.last ().optional = true;
   if (cli.last ().exists ())  {
      eAim = cli.last ().toDoubleList ();
   }   
   double errorAim = cli.option ("-eAim", "double Value", "Convergence Criteria")
                                   .toDouble (0.0,0.0,1.0);
   double error = cli.option ("-e", "double Value", "Convergence Criteria")
                                   .toDouble (0.0,0.0,errorAim);
   int nStates = cli.option ("-n|--nStates", "int Value", "Number of States")
                                   .toInt (1,1);
   bool keepWaves = cli.option ("--keepWaves", "bool", "Initialize next k-Point with old waves")
                                   .toBool ();
   bool keepEnergies = cli.option ("--keepEnergies", "bool", "Initialize next k-Point with old energie")
                                   .toBool ();
   bool direct = cli.option ("--direct", "bool", "Solve Hamiltonian directly")
                                   .toBool ();
   int initOnSpecies = cli.option ("--initOnSpecies", "int", "Give Species for Orbitalinitialization, other will set to zero")
                                   .toInt (-1);

   if (fabs(error) < 1e-12) error = errorAim;

   cli.finalize ();

   if (eAim.getSize () == 0)  {
      cout << "Warning: eAim not set! Now set to 1.0 eV" << endl;
      eAim.resize(1);
      eAim.set(1.0);
   }

   if (eAim.getSize () != nStates && eAim.getSize() != 1)  {
      cout << "eAim has to have one or nStates entries!" << endl;
      SX_QUIT;
   }

   //eAim to Hartree
   eAim /= HA2EV;

   // --- read input file
   SxParser parser;
   SxConstPtr<SxSymbolTable> tablePtr = parser.read (inputFile);

   SxPtr<SxAOMatTB> tbHamPtr = SxPtr<SxAOMatTB>::create(&*tablePtr);

   tbHamPtr->info();

   SxKPoints kPoints (tbHamPtr->structure.cell, &*tablePtr);

   tbHamPtr->setKPoint(kPoints.getK(0));
   
   cout << "Starting experimental minimization" << endl;
   cout << "eAim is " << eAim << endl;
   cout << "keepWaves is ";
   if (keepWaves) cout << "enabled" << endl;
   else cout << "disabled" << endl;
   cout << "keepEnergies is ";
   if (keepEnergies) cout << "enabled" << endl;
   else cout << "disabled" << endl;

   if (direct) {
      tbHamPtr->directMin (kPoints);
      tbHamPtr->setKPoint(kPoints.getK(0));
   } else {
      cout << SX_SEPARATOR;
      // Start with random Function 
      int dim = tbHamPtr->getNOrbs ();
      SxWaves init(nStates,dim, tbHamPtr);
      init.randomize();
      SxMatrix<Double> result (kPoints.getNk (),nStates);
      SxVector<Double> epsilon(nStates);
      if (eAim.getSize() == 1) epsilon.set(eAim(0));
      else for (int i = 0; i < nStates; i++)  epsilon(i) = eAim(i);
      bool opt = true;
      bool eAimFinished = false;
      for (int ik = 0; ik < kPoints.getNk (); ik++)  {
         tbHamPtr->setKPoint(kPoints.getK(ik));
         if (initOnSpecies >= 0) {
            if (initOnSpecies >= tbHamPtr->structure.getNSpecies ()) {
               cout << "InitOnSpecies higher then number of species in structure" << endl;
               SX_QUIT;
            }
            init.occupySpecies(initOnSpecies);
         }
         init.orthonormalize();

         if (!keepEnergies) {
            if (eAim.getSize() == 1) epsilon.set(eAim(0));
            else for (int i = 0; i < nStates; i++)  epsilon(i) = eAim(i);
         }
         (epsilon * HA2EV).print(true);
         SxHeS HES (tbHamPtr, init, epsilon);

         double F = 0.0, FNew = HES.getFunctional();
         opt = true;
         eAimFinished = false;

         int iSteps = 0;
         bool restart = true;
         SxWaves grad(nStates, dim, tbHamPtr), 
                 gradOld(nStates, dim, tbHamPtr), 
                 dir(nStates, dim, tbHamPtr), 
                 conjGrad(nStates, dim, tbHamPtr), 
                 conjGradOld(nStates, dim, tbHamPtr);

         while ((opt && iSteps < 1e3) || restart)  {
            iSteps++;
            F = FNew;
            grad = HES.getGradient ();

            SxVector<Double> beta (nStates);
            beta.set(0.0); 
            if (restart)  {
               // SD Direction
               conjGrad = grad;
               restart = false;
            } else  {
               // Conjugate gradients
               for (int iState = 0; iState < nStates; iState++)  {   
                  if (gradOld(iState).norm() > 1e-6)  
                     beta(iState) = dot(grad(iState),grad(iState)) 
                        / dot(gradOld(iState),gradOld(iState));
                  conjGrad(iState) = grad(iState) + beta(iState) * conjGradOld(iState);
               }
            }

            dir = conjGrad;
            // Orthogonolize directions
            for (int iState = 0; iState < nStates; iState++)  {
               for (int jState = 0; jState < nStates; jState++)  {
                  dir(iState) -= dot(dir(iState),HES.getSPsi(jState)).conj () 
                                 * HES.getPsi(jState);
               }
            }

            SxHeS HESDir = SxHeS(tbHamPtr, dir, epsilon);

            SxVector<Double> step (nStates);
            step.set(0.0);
            double c = HES.getFunctional();
            double b = 0.0, a = 0.0;
            SxVector<Int> satSteps (nStates);
            satSteps.set(0);
            for (int iState = 0; iState < nStates; iState++)  {
               // get optimal Stepwidth via parabola Fit FGuess = ax^2 + bx + c
               b = -2.0 * dot(dir(iState),grad(iState)).re;
               SxVector<Double> guess (nStates);
               guess.set(0.0);
               guess(iState) = 1.0;
               bool notSatisfied = true;
               double FGuess = 0.0;
               SxHeS HESGuess;
               while (notSatisfied  && satSteps(iState) < 20)  {
                  notSatisfied = false;
                  HESGuess = HES - guess * HESDir;
                  if (!eAimFinished) {
                     if (keepEnergies) HESGuess.orthonormalize(true);
                     else HESGuess.orthonormalize(true);
                  }
                  else HESGuess.normalize();
                  FGuess = HESGuess.getFunctional();
                  if (fabs(guess(iState)) > 1e-8) 
                     a = (FGuess - c - guess(iState) * b) / guess(iState) / guess(iState);
                  else a = 0.0;
                  if ((fabs(a) > 1e-8) && (a > 0.0)) step(iState) = -0.5 * b / a;
                  else if (a < 0)  {
                     guess(iState) *= 10.0;
                     step(iState) = 0.0;
                     notSatisfied = true;
                     restart = true;
                  }
                  else step(iState) = 0.0;
                  if (fabs(step(iState) - guess (iState)) > 1e-2)  {
                     satSteps(iState)++;
                     if (satSteps(iState) < 20) guess(iState) = step(iState);
                     notSatisfied = true;
                  }
               }

               if ((a < 0)) { restart = true; }
               //for debug
               if ((a < 0) || satSteps(iState) >=20) {
                  SxVector<Double> x (101);
                  SxVector<Double> y (101);
                  SxVector<Double> z (101);

              if (guess(iState) < 0) guess(iState) *= -1.0; 
               for (int i = 0; i < 101; i++)  {
                  x(i) = 0.1 * i * guess(iState);
                  SxVector<Double> factor (nStates);
                  factor.set(0.0);
                  factor(iState) = x(i);
                  SxHeS HESDraw = HES - factor * HESDir;
                  if (!eAimFinished)  { 
                     if (keepEnergies) HESDraw.orthonormalize(true);
                     else HESDraw.orthonormalize(false);
                  }
                  else HESDraw.normalize();
                  y(i) = HESDraw.getFunctional();
                  z(i) = a * x(i) * x(i) + b * x(i) + c;
                  if ((fabs(c - y(0)) > 1e-12) && i==0 ) {
                     cout << "Failed: y(0) not equal c" << endl; 
                     cout << "x(0) = " << x(0) << " should be 0" << endl;
                     cout << "y(0) = " << y(0) << endl;
                     cout << "c = " << c << endl;
                     HES.getPsi().checkNorm();
                     HESDraw.getPsi().checkNorm();
                     cout << "HESPsi - DrawPsi = " << (HES.getPsi(iState) - HESDraw.getPsi(iState)).norm () << endl;
                     cout << "HESHeS2Psi - DrawHeS2Psi = " << (HES.getHeS2Psi(iState) - HESDraw.getHeS2Psi(iState)).norm () << endl;
                     cout << "HESEpsilon - DrawEpsilon = " << (HES.getEpsilon()(iState) - HESDraw.getEpsilon()(iState)) << endl;
                  }
                  if (fabs(x(i) - guess(iState)) < 1e-12) {
                     cout << "CHECK GUESS CONSISTENCY" << endl;
                     if (fabs(y(i)- FGuess) > 1e-12)  {
                        cout << "Failed: y(guess) not equal FGuess" << endl; 
                        cout << "x(guess) = " << x(i) << ", y(guess) = " << y(i) << endl;
                        cout << "guess = " << guess(iState) << ", FGuess = " << FGuess << endl;
                        factor.print();
                        guess.print();
                        cout << "DrawPsi - GuessPsi = " << (HESDraw.getPsi(iState) - HESGuess.getPsi(iState)).norm () << endl;

                     }
                     if (fabs(y(i)- FGuess) > 1e-12)  {
                        cout << "Failed: y(guess) not equal z(guess)" << endl; 
                        cout << "x(guess) = " << x(i) << ", y(guess) = " << y(i) << endl;
                        cout << "x(guess) = " << x(i) << ", z(guess) = " << z(i) << endl;
                     }
                  }
               }
               SxString file = "linemin.dat";
               SxBinIO out;
               out.open(file, SxBinIO::ASCII_WRITE_ONLY);
               out.writeXYPlot(x,y);
               out.close();
               file = "lineminFit.dat";
               out.open(file, SxBinIO::ASCII_WRITE_ONLY);
               out.writeXYPlot(x,z);
               out.close();

               cout << SX_SEPARATOR;
               cout << "curent Step " << iSteps << endl;
               cout << "curent State " << iState << endl;
               cout << "kPoint = " << kPoints.getK(ik) << endl;
               cout << "Error during line minimization!" << endl;
               cout << "guess = " << guess << endl;
               cout << "step = " << step << endl;
               cout << "lineMinSteps = " << satSteps << endl;
               cout << "eAimFinished = " << eAimFinished << endl;
               cout << "a = " << a << endl;  
               cout << "b = " << b << endl;  
               cout << "c = " << c << endl;  
               cout << "beta = " << beta << endl;  
               cout << "Energies = " << HES.getEnergy() * HA2EV << endl; 
               cout << "<grad|grad> = " << grad(iState).norm () << endl; 
               cout << "<gradOld|gradOld> = " << gradOld(iState).norm () << endl; 
               cout << "<dir|dir> = " << dir(iState).norm () << endl; 
               cout << "<CG|CG> = " << conjGrad(iState).norm () << endl; 
               cout << SX_SEPARATOR;
               SX_EXIT;
            }
         }

         HES -= step * HESDir;
         if (!eAimFinished)  {
            if (keepEnergies) HES.orthonormalize (true);
            else HES.orthonormalize (false);
         }
         else HES.normalize ();
         HES.updateKappa ();
         FNew = HES.getFunctional();
         
         SxVector<Double> gNorm (nStates), ogNorm(nStates), dNorm(nStates),cgNorm(nStates);
         for (int iState = 0; iState < nStates; iState++)  {
            gNorm (iState) = grad(iState).norm();
            ogNorm (iState) = gradOld(iState).norm();
            dNorm (iState) = dir(iState).norm();
            cgNorm (iState) = conjGrad(iState).norm();
         }
         
         if (fabs(FNew - F) < errorAim && !eAimFinished)  {
            cout << SX_SEPARATOR;
            cout << "eAim reached, now optimize eigenvalues!" << endl;
            cout << SX_SEPARATOR;
            eAimFinished = true;  
            restart = true;
         }
         if (eAimFinished && (fabs(FNew - F) < error) && !restart) opt = false;

                 
         cout << SX_SEPARATOR;
         cout << "ik = " << ik << endl;
         cout << "kPoint = " << kPoints.getK(ik) << endl;
         cout << "curent Step " << iSteps << endl;
         cout << "step = " << step << endl;
         cout << "lineMinSteps = " << satSteps << endl;
         cout << "beta = " << beta << endl; 
         cout << "a = " << a << endl; 
         cout << "b = " << b << endl; 
         cout << "c = " << c << endl; 
         cout << "Energies = " << HES.getEnergy() * HA2EV << endl;
         cout << "Functional = " << HES.getFunctional () << endl;
         cout << "<grad|grad> = " << gNorm << endl; 
         cout << "<gradOld|gradOld> = " << ogNorm << endl; 
         cout << "<dir|dir> = " << dNorm << endl; 
         cout << "<CG|CG> = " << cgNorm << endl; 
         cout << SX_SEPARATOR;

         gradOld = grad;
         conjGradOld = conjGrad;

         if (eAimFinished) {
            SxMatrix<Complex16> rot = HES.setEnergyByCholesky ();
            epsilon = 1.0 * HES.getEpsilon ();
            SxWaves orig = gradOld;
            for (int i = 0; i < nStates; i++)  {
               gradOld(i).set(0.0);
               for (int j = 0; j < nStates; j++)  {
                  gradOld(i) += rot(j,i) * orig(j);
               }
            }
            FNew = HES.getFunctional ();
         }
      }
      if ( (ik == 0) && keepWaves) {
         init = HES.getPsi ();
      }

      for (int iState = 0; iState < nStates; iState++)  {
         result(ik,iState) = HES.getEnergy(iState) * HA2EV;
         cout << "ik = " << ik 
            << ", iState = " << iState 
            << ", final Residuum is " << HES.getResiduum(iState) 
            << ", eigenVal = " << result(ik,iState) 
            << ", converged after " << iSteps << " steps" << endl;
      }
   }

   SxBinIO out;
   SxString file = "eps-eAim-method.dat";
   out.open(file, SxBinIO::ASCII_WRITE_ONLY);
   out.writeNXYPlot(result);
   out.close();
   }

   printTiming ();

}

#endif /* SX_STANDALONE */
