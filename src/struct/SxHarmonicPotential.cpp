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

#include <SxHarmonicPotential.h>
#include <SxNeighbors.h>
#include <SxProjector.h>
#include <SxGaussIJ.h>
#include <SxUniqueList.h>
#include <SxSimpleParser.h>
#include <SxEquivalence.h>

double SxHarmonicPotential::betaSoft = 3.;

namespace {
   void checkSymId (const SxVector<Int> &symId, ssize_t nSym)
   {
      for (int iSym = 0; iSym < symId.getSize (); ++iSym)  {
         if (symId(iSym) >= nSym || symId(iSym) < 0) {
            cout << "Invalid symmetry id = " << (symId(iSym)+1)
                 << ". Only " << nSym << " symmetries have been defined."
                 << endl;
            SX_QUIT;
         }
      }
   }
}

SxHarmonicPotential::SxHarmonicPotential (const SxSymbolTable *table,
                                          const SxAtomicStructure &refStruct)
   : beta (1.0), energy (-1.)
{
   refStr.copy (refStruct);
   if (table->containsGroup ("harmonicPotential"))
      table = table->getGroup ("harmonicPotential");
   readSpecies (table);

   SYMBOLPARSE(table) {
      // --- symmetry
      SYMBOLGROUP("symmetry")  {
         SxStack<SymMat> symList;
         FOREACH_SYMBOLGROUP ("operator")
            symList << SymMat(SYMBOLGET("S")->toList());
         syms = symList;
      }
      beta << SYMBOLGET("beta");
      dielecConstant = SYMBOLGET("dielecConstant") || 1.;
   }

   setupQab ();
   bool refIsRelaxed = table->containsGroup ("defect");

   SxGrid grid(refStruct, 10);
   SxNeighbors nn;
   int mode = SxNeighbors::StoreIdx 
            | SxNeighbors::StoreRel
            | SxNeighbors::IncludeZeroDistance;

   ssize_t nAtoms = refStruct.getNAtoms ();
   SxArray<SxStack<int> > neighborStack(nAtoms);
   SxArray<SxStack<SxMatrix3<Double> > > hesseStack (nAtoms);
   SxArray<SxStack<Coord> > dist0Stack (nAtoms);
   SxArray<SxStack<int> > symIdStack (nAtoms);
   SxArray<SxStack<int> > neighborIdStack (nAtoms);
   SxStack<Coord> genNeighborStack;
   SxStack<SxMatrix<Double> > matSymStack, matSymQStack;
   genQ.resize (nAtoms);
   genQ.set (-1);
   genQSymId.resize (nAtoms);
   genQSymId.set (-1);

   SxPtr<SxSymGroup> fullSym;

   double dist = 0.5;
   SYMBOLPARSE(table) {
      int iGenNeighbor = 0;
      refIsRelaxed << SYMBOLGET("relaxedRef");

      // --- read in equivalence tables
      SxEquivalence equiv;
      SYMBOLGROUP ("equivalence") {
         equiv.read (SYMBOLGROUP_TABLE);
         equiv.mapSyms (syms);
         SX_LOOP(ia) {
            if (equiv.mappingRot(ia) < 0) {
               cout << "Equivalence symmetries do not map" << endl;
               SX_QUIT;
            }
         }
         SX_LOOP(iEq)  {
            for (int iSym = 0; iSym < equiv.matSym(iEq).getSize (); ++iSym)
               if (equiv.matSym(iEq)(iSym) < 0) {
                  cout << "Equivalence symmetries do not map" << endl;
                  SX_QUIT;
               }
         }
      }
      if (equiv.getSize () > 0)  {
         if (HAVE_SYMBOL("equivalenceIds"))  {
            equiv.equivId = SYMBOLGET("equivalenceIds");
            equiv.mappingRot = SYMBOLGET("mappingSym");
         } else {
            dist = SYMBOLGET ("maxDist") || 0.5;
            cout << "Mapping supercell to equivalence cell" << endl;
            equiv = equiv.mapToComplex (refStruct, dist);
         }
      }

      SxArray<SxMatrix3<Double> > bulkH;
      if (HAVE_SYMBOL("bulkHesse"))  {
         SxList<double> vals = SYMBOLGET("bulkHesse")->toList ();
         ssize_t nH = vals.getSize () / 9;
         if (nH * 9 != vals.getSize ())  {
            cout << "Error: bulkHesse has wrong format" << endl;
            cout << "Should be nHesse x 3 x 3" << endl;
            SX_QUIT;
         }
         SxList<double>::Iterator it = vals.begin ();
         bulkH.resize (nH);
         SX_LOOP3(iH,i,j) bulkH(iH)(j,i) = *it++;
      }

      // --- read species data
      int is = 0;
      FOREACH_SYMBOLGROUP("species")  {

         valenceCharge(is) = SYMBOLGET("valenceCharge") || 0.;

         // --- equivalence-class dependent charge tensors/charges
         SYMBOLGROUP("equivalenceClass") {
            int eqId = SYMBOLGET("id");
            if (eqId < 1 || eqId > equiv.getSize ())  {
               cout << "Invalid equivalence id = " << eqId << endl;
               SX_QUIT;
            }
            eqId--;
            SxMatrix3<Double> Qeff;
            if (HAVE_SYMBOL("chargeTensor"))  {
               Qeff = SxMatrix3<Double> (SYMBOLGET("chargeTensor")->toList ());
            } else if (HAVE_SYMBOL("charge")) {
               double Q = SYMBOLGET("charge");
               Qeff = SxMatrix3<Double> (Q, 0., 0., 0., Q, 0., 0., 0., Q);
            }

            // --- assign charge tensor to all atoms of this equivalence class
            bool anyAtom = false;
            for (int ia = 0; ia < refStr.getNAtoms (is); ++ia)  {
               int iTlAtom = refStr.atomInfo->offset(is) + ia;
               if (equiv.equivId (iTlAtom) == eqId)  {
                  const SxMatrix3<Double> &S = syms(equiv.mappingRot(iTlAtom));
                  Qab(iTlAtom) = S ^ Qeff ^ S.transpose ();
                  genQ(iTlAtom) = (int)matSymQStack.getSize ();
                  genQSymId(iTlAtom) = equiv.mappingRot(iTlAtom);
                  cout << (iTlAtom + 1) << ' ';
                  anyAtom = true;
               }
            }
            if (anyAtom)
               cout << "=> tensor " <<  matSymStack.getSize () << endl;

            // --- tensor symmetries
            SxList<int> QsymId;
            SX_LOOP(iSym) QsymId << equiv.matSym(eqId)(iSym);
            cout << "Charge tensor " << matSymStack.getSize () << ": ";
            matSymQStack << getMatrixConstraint (QsymId, false);
         }

         // --- interatomic force constants
         FOREACH_SYMBOLGROUP ("neighbor")  {
            dist = SYMBOLGET ("maxDist") || 0.5;

            SxVector<Int> symId;
            if (HAVE_SYMBOL("symmetries"))  {
               symId = SxVector<Int> (SYMBOLGET ("symmetries")->toIntList ());
               symId -= 1; // we start from 0, but in input from 1
               checkSymId (symId, syms.getSize ());
            }
            // read reference Hesse matrix
            SxMatrix3<Double> hesseRef;
            if (HAVE_SYMBOL("Hesse"))
               hesseRef = SxMatrix3<Double>(SYMBOLGET("Hesse")->toList ());
            else if (iGenNeighbor < bulkH.getSize ())
               hesseRef = bulkH(iGenNeighbor);
            else
               hesseRef.set (0.);
            // read reference neighbor position
            Coord posRef (SYMBOLGET("coords")->toList ());
            int isRef = find (SYMBOLGET("atomType") || "");
            genNeighborStack << posRef;

            if (HAVE_SYMBOL("matSyms"))  {
               SxList<int> matSymId = SYMBOLGET("matSyms")->toIntList ();
               for (int iSymM = 0; iSymM < matSymId.getSize (); ++iSymM)  {
                  matSymId(iSymM) -= 1;
                  SymMat S = syms(matSymId(iSymM));
                  if (((S ^ posRef) - posRef).normSqr () > 1e-10)  {
                     cout << "Matrix symmetry does not preserve position"
                          << endl;
                     cout << "Neighbor is " << posRef << endl;
                     cout << "Offending symmetry is " << (matSymId(iSymM)+1)
                          << ": " << S << endl;
                     cout << "This is inconsistent input." << endl;
                     SX_QUIT;
                  }
               }
               matSymStack << getMatrixConstraint (matSymId);
            } else {
               matSymStack << getMatrixConstraint (SxList<int> ());
            }

            if (HAVE_SYMBOL("equivalenceId"))  {
               // --- add neighbors for all atoms belonging to iEq class
               int iEq = SYMBOLGET("equivalenceId")->toInt () - 1;

               for (int ia = 0; ia < refStr.getNAtoms (is); ++ia)  {
                  int iTl = refStruct.getIAtom (is, ia);
                  // check that this atom has the right equivalence class
                  if (equiv.equivId (iTl) != iEq) continue;

                  // rotates that maps representative atom to iTl
                  const SymMat &rotI = syms(equiv.mappingRot (iTl));

                  // --- loop over local symmetries
                  for (int iSym = symId.getSize () > 0 ? 0 : -1;
                       iSym < symId.getSize ();
                       ++iSym)
                  {
                     SymMat S = iSym >= 0
                              ? (rotI ^ syms(symId(iSym)))
                              : rotI;
                     Coord jCoord = refStruct(iTl) + (S ^ posRef);
                     nn.compute (grid, refStruct, jCoord, dist, mode);
                     if (nn.getSize () == 1)  {
                        // found a neighbor: check species
                        int jTl = nn.idx(0);
                        if (isRef < 0 || refStruct.getISpecies(jTl) == isRef) {
                           Coord d0 = refStruct(jTl)-refStruct(iTl);
                           // find symmetry id for S
                           int jSym;
                           for (jSym = 0; jSym < syms.getSize (); ++jSym)
                              if (S == syms(jSym)) break;

                           neighborStack(iTl) << jTl;
                           hesseStack(iTl) << (S ^ hesseRef ^ S.transpose ());
                           dist0Stack(iTl) << d0;
                           neighborIdStack(iTl) << iGenNeighbor;
                           symIdStack(iTl) << jSym;
                        }
                     } else if (nn.getSize () > 1)  {
                        cout << "Ambiguous neighbors for atom "
                             << chemName(is) << (ia+1) << "@ "
                             << refStruct(is, ia) << " near " << jCoord << endl;
                        cout << "Candidates are " << endl;
                        for (int in = 0; in < nn.getSize (); ++in)  {
                           cout << "atom " << (nn.idx(in)+1) << endl;
                        }
                        SX_QUIT;
                     }
                  }
               }
            } else {
               // --- find all corresponding neighbors
               for (int iSym = symId.getSize () > 0 ? 0 : -1;
                    iSym < symId.getSize ();
                    ++iSym)
               {
                  SymMat S = iSym >= 0
                           ? syms(symId(iSym))
                           : SymMat (1,0,0, 0,1,0, 0,0,1);
                  SxMatrix3<Double> hesseMat = S ^ hesseRef ^ S.transpose ();
                  Coord pos = S ^ posRef;

                  for (int ia = 0; ia < refStruct.getNAtoms (is); ++ia)  {
                     int iTlAtom = refStruct.getIAtom (is, ia);
                     nn.compute (grid, refStruct, refStruct(is,ia) + pos,
                                 dist, mode);
                     if (nn.getSize () == 1)  {
                        // found a neighbor: check species
                        int jTl = nn.idx(0);
                        if (isRef < 0 || refStruct.getISpecies(jTl) == isRef) {
                           // found a neighbor with correct species
                           neighborStack(iTlAtom) << jTl;
                           hesseStack(iTlAtom) << hesseMat;
                           if (refIsRelaxed)  {
                              dist0Stack(iTlAtom) << (refStruct(jTl) - refStruct(is,ia));
                           } else {
                              // pos with additional lattice shift
                              Coord posR = refStruct(jTl) - nn.relPositions(0)
                                         - refStruct(is,ia);
                              dist0Stack(iTlAtom) << posR;
                           }
                           // keep track of neighbor origin
                           neighborIdStack(iTlAtom) << iGenNeighbor;
                           symIdStack(iTlAtom) << (iSym>=0 ? symId(iSym) : -1);
                        }
                     } else if (nn.getSize () > 1)  {
                        cout << "Ambiguous neighbors for atom "
                             << chemName(is) << (ia+1) << "@ "
                             << refStruct(is, ia)
                             << " near "
                             << (refStruct(is, ia) + pos)
                             << endl;
                        cout << "Candidates are " << endl;
                        for (int in = 0; in < nn.getSize (); ++in)  {
                           cout << "atom " << (nn.idx(in)+1) << endl;
                        }
                        SX_EXIT;
                     }
                  }
               }
            }
            iGenNeighbor++;
         }
         is++;
      }
      if (bulkH.getSize () > 0 && bulkH.getSize () != iGenNeighbor)  {
         cout << "Error: bulkHesse size does not fit number of neighbor types"
              << endl;
         cout << "bulkHesse contains " << bulkH.getSize ()
              << " Hesse matrices" << endl;
         cout << "Number of neighbor types is " << iGenNeighbor << endl;
         SX_QUIT;
      }
   }
   genNeighbors = genNeighborStack;
   matSyms = matSymStack;

   genSymId.resize (nAtoms);
   genNeighborId.resize (nAtoms);
   neighborIdx.resize (nAtoms);
   hesse.resize (nAtoms);
   dist0.resize (nAtoms);
   for (int iTlAtom = 0; iTlAtom < nAtoms; ++iTlAtom)  {
      neighborIdx(iTlAtom) = neighborStack(iTlAtom);
      hesse(iTlAtom) = hesseStack(iTlAtom);
      dist0(iTlAtom) = dist0Stack(iTlAtom);
      genNeighborId(iTlAtom) = neighborIdStack(iTlAtom);
      genSymId(iTlAtom) = symIdStack(iTlAtom);
      cout << "atom " << (iTlAtom+1) << ": " 
           << neighborIdx(iTlAtom).getSize () << " neighbors ->"
           << neighborIdx(iTlAtom)
           << endl;
   }
   SxArray<SxMatrix3<Double> > QabBulk;

   SxArray<int> genId;
   if (table->containsGroup ("defect"))  {
      QabBulk = Qab;
      Coord center;
      SxVector<Int> symId;
      SxList<Coord> dNeighbor;
      SxList<SxMatrix3<Double> > bornQ;
      SxList<SxArray<SxMatrix3<Double> > > defectHesse;
      SYMBOLPARSE (table->getGroup ("defect")) {
         center = Coord (SYMBOLGET("center")->toList ());
         if (HAVE_SYMBOL("symmetries"))  {
            symId = SxVector<Int> (SYMBOLGET("symmetries")->toIntList ());
            symId -= 1; // we start from 0, but in input from 1
            checkSymId (symId, syms.getSize ());
         }
         
         FOREACH_SYMBOLGROUP("neighbor") {
            dNeighbor << Coord (SYMBOLGET("coords")->toList ());
            SxMatrix3<Double> Qeff (SYMBOLGROUP_TABLE->get("chargeTensor", true)
                                    ->toList ());
            bornQ << Qeff;
            SxMatrix<Double> hesses;
            SxArray<SxMatrix3<Double> > hesseList;
            if (HAVE_SYMBOL("defectHesse"))  {
               hesses = SYMBOLGET("defectHesse")->toList ();
               int nHesses = (int)hesses.getSize () / 9;
               if (hesses.getSize () != nHesses * 9) { SX_EXIT; }
               hesseList.resize (nHesses);
               for (int iH = 0; iH < nHesses; ++iH)
                  for (int i = 0; i < 9; ++i)
                     hesseList(iH)(i / 3, i % 3) = hesses (i + 9 * iH);
            }
            defectHesse << hesseList;
         }
      }
      nGenSets.resize (int(dNeighbor.getSize ()));

      genId.resize (nAtoms);
      genId.set (-1);

      int nDefAtoms = 0;
      // --- find and update all defect atoms
      for (int in = 0; in < dNeighbor.getSize (); ++in)  {
         cout << "Defect atom type " << (in+1) << endl;
         SxList<int> idn, idSym, symN;
         // --- find all symmetry-equivalent defect atoms
         Coord posRef0;
         for (int iiSym = -1; iiSym < symId.getSize (); ++iiSym)
         {
            // the actual symmetry number (or -1 for identity)
            int iSym = iiSym >= 0 ? symId(iiSym) : -1;
            SymMat S = iSym >= 0 
                     ? syms(iSym)
                     : SymMat (1,0,0, 0,1,0, 0,0,1);
            Coord pos = S ^ dNeighbor(in);

            nn.compute (grid, refStruct, center + pos, dist, mode);

            if (nn.getSize () == 1)  {
               int jTlAtom = nn.idx(0);
               if (genId(jTlAtom) == -1)  {
                  idn << jTlAtom;
                  idSym << iSym;
                  genId(jTlAtom) = in;
                  if (iSym == -1) 
                     posRef0 = pos + nn.relPositions (0);
                  cout << "New neighbor: " << (jTlAtom+1) << endl;
                  // set charge tensor
                  Qab(jTlAtom) = S ^ bornQ(in) ^ S.transpose ();
                  genQ(jTlAtom) = (int)matSymQStack.getSize ();
                  genQSymId(jTlAtom) = iSym;
               } else {
                  // --- this leaves the representative defect atom in place
                  if (jTlAtom == idn(0))  {
                     symN << iSym;
                  }
                  // --- check charge tensor
                  SxMatrix3<Double> dQ =  Qab(jTlAtom)
                                       - (S ^ bornQ(in) ^ S.transpose ());
                  if (dQ.absSqr ().sum () > 1e-6)  {
                     cout << "Charge tensor inconsistency: " << endl;
                     cout << "First match: " << Qab(jTlAtom) << endl;
                     cout << "This match: " << (S ^ bornQ(in) ^ S.transpose ())
                          << endl;
                     cout << "Difference: " << dQ << endl;
                     SX_QUIT;
                  }
               }
            } else {
               cout << "Failed to resolve defect neighbor @ "
                    << center + pos << endl;
               cout << "Candidates are";
               if (nn.getSize () > 0) cout << " missing.";
               cout << endl;
               for (int jn = 0; jn < nn.getSize (); ++jn)  {
                  cout << "atom " << (nn.idx(jn)+1) << endl;
               }
               SX_EXIT;
            }
         }
         nDefAtoms += int(idn.getSize ());
         int iTlAtom = idn(0);
         cout << "Charge tensor " << matSymQStack.getSize () << ": ";
         matSymQStack << getMatrixConstraint (symN, false);

         cout << genNeighborId(iTlAtom) << endl;
         cout << genSymId(iTlAtom) << endl;

         // --- shift interactions with non-defect atoms to defect atoms
         for (int ia = 0; ia < neighborIdx.getSize (); ++ia)  {
            if (genId(ia) >= 0)  continue;
            for (int jn = 0; jn < neighborIdx(ia).getSize (); ++jn) {
               if (genId(neighborIdx(ia)(jn)) >= 0)  {
                  shiftHesseToNeighbor (ia, jn);
                  jn--;
               }
            }
         }

         int nNewN = int(genNeighborId(iTlAtom).getSize ());
         cout << "=> " << nNewN << " neighbors" << endl;

         // --- symmetry-reduce neighbors of the representative defect atom
         cout << posRef0 << endl;
         SxAtomicStructure relN (nNewN);
         for (int idN = 0; idN < nNewN; ++idN)  {
            // try to estimate position wrt defect center
            // assumption: idN is within Wigner-Seitz cell of iTlAtom
            Coord pos = refStruct(neighborIdx(iTlAtom)(idN))
                      - refStruct(iTlAtom);
            refStruct.cell.map (&pos, SxCell::WignerSeitz);
            relN.ref (idN) = posRef0 + pos;
         }
         SxArray<int> representative(nNewN), idSymN(nNewN);
         representative.set (-1);
         SxList<int> genNId;
         // loop over neighbors of the defect atom
         for (int idN = 0; idN < nNewN; ++idN)  {
            if (representative(idN) >= 0) continue;
            representative(idN) = int(genNId.getSize ());
            genNId << idN;
            idSymN(idN) = -1;
            cout << genNId.getSize () << " -> " << relN(idN) << endl;
            SxList<int> repMatSyms;
            // loop over the symmetries that leave defect atom in place
            for (int iiSym = 0; iiSym < symN.getSize (); ++iiSym)  {
               SymMat S = (symN(iiSym) == -1) ? SymMat(1.,0.,0.,0.,1.,0.,0.,0.,1.) 
                                              : syms(symN(iiSym));
               Coord rotPos = S ^ relN(idN);
               for (int kn = 0; kn < nNewN; ++kn)  {
                  if (representative(kn) >= 0) continue;
                  if ((rotPos - relN(kn)).normSqr () < dist*dist)  {
                     representative(kn) = int(genNId.getSize ()) - 1;
                     idSymN(kn) = symN(iiSym);
                  }
               }
               if ((rotPos - relN(idN)).normSqr () < dist*dist)  {
                  repMatSyms << symN(iiSym);
               }
            }
            matSymStack << getMatrixConstraint (repMatSyms);
         }
         int nNewGen = int(genNId.getSize ());

         // --- new set of generating neighbors
         cout << "=> " << nNewGen << " independent interactions" << endl;
         // --- setup defect Hesse matrix
         if (defectHesse(in).getSize () > 0)  {
            if (defectHesse(in).getSize () != nNewGen)  {
               cout << "Inconsistency in size of Hesse!"
                    << endl;
               cout << "Should be 9 x " << nNewGen << endl;
               cout << "but is 9 x " << defectHesse(in).getSize ()
                    << endl;
               SX_QUIT;
            }
         } else {
            defectHesse(in).resize (nNewGen);
            for (int iGen = 0; iGen < nNewGen; ++iGen)
               defectHesse(in)(iGen).set (0.);
         }

         nGenSets(in) = nNewGen;
         int nGen = int(genNeighbors.getSize ());
         genNeighbors.resize (nGen + nNewGen, true);
         for (int iGen = 0; iGen < nNewGen; ++iGen)  {
            Coord genRef = genNeighbors(genNeighborId(iTlAtom)(genNId(iGen)));
            
            SxMatrix3<Double> gS(1.,0.,0.,0.,1.,0.,0.,0.,1.);
            int iSym = genSymId(iTlAtom)(genNId(iGen)),
                nSym = int(syms.getSize ());
            if (iSym >= 0)  {
               if (iSym < nSym)  {
                  gS = syms(iSym);
               } else {
                  // special code: this was shifted from a non-defect atom
                  gS = syms(iSym - nSym).transpose ();
                  genRef = -genRef;
               }
            } else if (iSym == -2) {
               genRef = -genRef;
            }
            genNeighbors(nGen + iGen) = gS ^ genRef;
         }

         for (int kn = 0; kn < nNewN; ++kn)  {
            // --- reset generating info for representative defect atom
            genNeighborId(iTlAtom)(kn) = nGen + representative(kn);
            genSymId(iTlAtom)(kn) = idSymN(kn);

            SxMatrix3<Double> H = defectHesse(in)(representative(kn));
            if (genSymId(iTlAtom)(kn) >= 0)  {
               SymMat S = syms(genSymId(iTlAtom)(kn));
               H = S ^ H ^ S.transpose ();
            }
            hesse(iTlAtom)(kn) = H;
         }
         SxArray2<int> symMulTable = SxSymGroup::getSymMulTable(syms);

         // --- update other defect atoms
         for (int i = 1; i < idn.getSize (); ++i)  {
            int jTlAtom = idn(i);
            int iSym = idSym(i);
            SX_CHECK (iSym >= 0);

            if (neighborIdx(jTlAtom).getSize () != nNewN)  {
               cout << "Inconsistency for atom " << (jTlAtom+1) << endl;
               cout << "Should have " << nNewN << " neighbors, but has "
                    << neighborIdx(jTlAtom).getSize () << " neighbors." << endl;
               SX_EXIT;
            }

            for (int jn = 0; jn < nNewN; ++jn)  {
               // --- generating info
               genNeighborId(jTlAtom)(jn) = genNeighborId(iTlAtom)(jn);
               int repSym = genSymId(iTlAtom)(jn);
               if (repSym >= 0)
                  genSymId(jTlAtom)(jn) = symMulTable(iSym, repSym);
               else
                  genSymId(jTlAtom)(jn) = iSym;

               // --- rotated Hesse matrix
               SymMat S = syms(iSym);
               hesse(jTlAtom)(jn) = S^hesse(iTlAtom)(jn)^S.transpose ();

               // --- find right neighbor
               Coord pos = refStruct.getAtom(jTlAtom) 
                         + (syms(genSymId(jTlAtom)(jn)) 
                            ^ genNeighbors(genNeighborId(jTlAtom)(jn)));

               nn.compute (grid, refStruct, pos, dist, mode);
               if (nn.getSize () != 1)  { 
                  cout << nn.relPositions << endl;
                  cout << "Expected to find " 
                       << " atom near " << pos << endl;
                  cout << "as a neighbor of "
                       << chemName(refStruct.getISpecies (jTlAtom))
                       << " @ " << refStruct.getAtom(jTlAtom)
                       << " at a distance of "
                       << (syms(genSymId(jTlAtom)(jn))
                          ^ genNeighbors(genNeighborId(jTlAtom)(jn)))
                       << endl;
                  cout << "obtained from "
                       << genNeighborId(jTlAtom)(jn) << ": "
                       << genNeighbors(genNeighborId(jTlAtom)(jn))
                       << " via symmetry " << genSymId(jTlAtom)(jn) + 1<< endl;
                  cout << "Extended search distMax=" << dist * 3. << endl;
                  nn.compute (grid, refStruct, pos, dist * 3., 
                              mode | SxNeighbors::StoreAbs);
                  if (nn.getSize () > 0)  {
                     SX_LOOP(iac) { 
                        cout << chemName(nn.absPositions.getISpecies ((int)iac))
                             << " @ " << nn.absPositions(iac)
                             << " d=" << nn.relPositions(iac).norm () << endl;
                     }
                  }
                  cout << genSymId(iTlAtom)(jn) << endl;
                  Coord pos1 = refStruct.getAtom(iTlAtom),
                        pos2 = refStruct.getAtom (neighborIdx(iTlAtom)(jn));

                  cout << "--- Data for the reference defect atom / neighbor" << endl;
                  if  (genSymId(iTlAtom)(jn) >= 0)
                     pos1 += (syms(genSymId(iTlAtom)(jn))
                                ^  genNeighbors(genNeighborId(iTlAtom)(jn)));
                  else
                     pos1 += genNeighbors(genNeighborId(iTlAtom)(jn));
                  cout << pos1 << " " << pos2 << ": "
                       << refStruct.cell.getMapped (pos2-pos1, SxCell::Origin)
                       << " d="
                       << refStruct.cell.getMapped (pos2-pos1, SxCell::Origin).norm ()
                       << endl;
                  cout << "Stopping due to defect neighbor mismatch." << endl;
                  SX_EXIT; 
               }
               neighborIdx(jTlAtom)(jn) = nn.idx(0);

               // pos with additional lattice shift
               if (refIsRelaxed)  {
                  dist0(jTlAtom)(jn) = refStruct(nn.idx(0)) - refStruct(jTlAtom);
               } else {
                  dist0(jTlAtom)(jn) = refStruct(nn.idx(0)) - nn.relPositions(0)
                                     - refStruct(jTlAtom);
               }
            }
         }

         // --- no duplicate interactions among equivalent defect neighbors
         for (int i = 0; i < idn.getSize (); ++i)  {
            int ia = idn(i);
            for (int jn = 0; jn < neighborIdx(ia).getSize (); ++jn) {
               int ja = neighborIdx(ia)(jn);
               if (genId(ja) == in && ja > ia)  {
                  shiftHesseToNeighbor (ia, jn);
                  jn--;
               }
            }
         }
      }
      // update list of matrix symmetry constraints 
      matSyms = matSymStack;
      matSymQ = matSymQStack;

      cout << SX_SEPARATOR;
      cout << "| Defect @ " << center << endl;
      cout << "| Number of primitive neighbors: " << dNeighbor.getSize () 
           << endl;
      cout << "| Number of all neighbors: " << nDefAtoms << endl;
      cout << SX_SEPARATOR;
   }

   setupDeltaQ (refStruct, QabBulk, genId);

   if (!noCoulomb)  {
      double gCut = -log(1e-16)/(betaSoft*betaSoft/4.);
      SxMesh3D mesh = SxGBasis::getCommensurableMesh (gCut, refStruct.cell);
      cout << "Preparing G basis ... ";
      gBasis = SxPtr<SxGBasis>::create (mesh, refStruct, gCut);
      cout << gBasis->ng << " G vectors" << endl;
   }
}

SxMatrix<Double> SxHarmonicPotential::getParamSpace () const
{
   int nFreeParam = 0;
   for (int i = 0; i < matSyms.getSize (); ++i)  {
      nFreeParam += (int)matSyms(i).nCols ();
   }
   SxMatrix<Double> res(9 * matSyms.getSize (), nFreeParam);
   res.set (0.);
   for (int i = 0, iFree = 0; i < matSyms.getSize (); ++i)  {
      for (int j = 0; j < matSyms(i).nCols (); ++j, ++iFree)  {
         res.colRef (int(iFree))(SxIdx(9 * i, 9 * i + 8))
            <<= matSyms(i).colRef (j);
      }
   }
   return res;
}

void
SxHarmonicPotential::setupDeltaQ (const SxAtomicStructure &refStruct,
                                  const SxArray<SxMatrix3<Double> > &QabBulk,
                                  SxArray<int> genId)
{
   noCoulomb = true;
   int nAtoms = refStruct.getNAtoms ();

   // --- find anisotropic atoms
   qIso.resize (nAtoms);
   mapAniso.resize (nAtoms);
   int nAniso = 0;
   for (int ia = 0; ia < nAtoms; ++ia)  {
      const SxMatrix3<Double> Q = Qab(ia);
      mapAniso(ia) = 0;
      qIso(ia) = Q.trace () / 3.;
      if (fabs(qIso(ia)) > 1e-10) noCoulomb=false;
      for (int i = 0; i < 3; ++i)  {
         for (int j = 0; j < 3; ++j)  {
            if (i == j)  {
               if (fabs(Q(i,i) - qIso(ia)) > 1e-10) mapAniso(ia) = 1;
            } else {
               if (fabs(Q(i,j))            > 1e-10) mapAniso(ia) = 1;
            }
         }
      }
      /*
      if (QabBulk.getSize () == Qab.getSize () && !mapAniso(ia))  {
         // --- mark atom as anisotropic if it is so in the bulk
         double QBulk = QabBulk(ia).trace () / 3.;
         SX_LOOP2(i(3),j(3))  {
            if (i==j) {
               if (fabs(QabBulk(ia)(i,i) - QBulk) > 1e-10) mapAniso(ia) = 1;
            } else {
               if (fabs(QabBulk(ia)(i,j))         > 1e-10) mapAniso(ia) = 1;
            }
         }
      }
      */
      if (mapAniso(ia)) nAniso++;
   }
   if (nAniso > 0) noCoulomb = false;

   // --- find neighbors
   SxArray<SxUniqueList<int> > frameList(nAtoms);
   SxArray<SxStack<Coord> >    posList(nAtoms);
   {
      SxGrid grid(refStruct, 10);
      SxNeighbors nn;
      double rad = 5.;
      int mode = SxNeighbors::StoreIdx | SxNeighbors::StoreRel;

      for (int ia = 0; ia < nAtoms; ++ia)  {
         if (mapAniso(ia) == 1)  {
            nn.compute (grid, refStruct, refStruct.getAtom(ia), rad, mode);
            for (int in = 0; in < nn.getSize (); ++in)  {
               int ja = nn.idx(in);
               frameList(ia) << ja;
               frameList(ja) << ia;
               if (frameList(ia).getSize () > (ssize_t)posList(ia).getSize ())
                  posList(ia) << nn.relPositions(in);
               if (frameList(ja).getSize () > (ssize_t)posList(ja).getSize ())
                  posList(ja) << -nn.relPositions(in);
               if (mapAniso(ja) == 0)  {
                  mapAniso(ja) = 2;
                  nAniso++;
               }
            }
         }
      }
   }
   cout << "nAniso = " << nAniso << endl;

   // --- get maps: allAtoms <-> aniso
   idAniso.resize (nAniso);
   frameAtoms.resize (nAniso);
   frameAtomPos.resize (nAniso);
   framePositions.resize (nAniso);
   for (int ia = 0, iAniso = 0; ia < nAtoms; ++ia)  {
      if (mapAniso(ia))  {
         mapAniso(ia) = iAniso;
         idAniso(iAniso) = ia;
         frameAtoms(iAniso) = frameList(ia);
         frameAtomPos(iAniso) = posList(ia);

         // --- position of atom with respect to reference frame
         ssize_t nFrame = frameAtoms(iAniso).getSize ();
         framePositions(iAniso) = refStruct(ia);
         for (int iF = 0; iF < nFrame; ++iF)
            framePositions(iAniso) -= refStruct(frameAtoms(iAniso)(iF)) 
                                    / double(nFrame);

         iAniso++;
      } else {
         mapAniso(ia) = -1;
      }
   }

   // no further setup needed if there are no anisotropic atoms
   if (nAniso == 0) return;

   // --- setup traceless anisotropy
   SxArray<SxMatrix3<Double> > dQ(nAniso);
   SxMatrix3<Double> sumQ(0,0,0,0,0,0,0,0,0);
   for (int iAniso = 0; iAniso < nAniso; ++iAniso)  {
      int ia = idAniso(iAniso);
      dQ(iAniso) = Qab(ia);
      for (int i = 0; i < 3; ++i) dQ(iAniso)(i,i) -= qIso(ia);
      sumQ += dQ(iAniso);
   }
   cout << "dQ sum " << sumQ << endl;

   // --- setup weights
   SxMatrix<Double> wik(nAniso, nAniso);
   wik.set (0.);
   for (int iAniso = 0; iAniso < nAniso; ++iAniso)  {
      int ia = idAniso(iAniso);
      ssize_t nFrame = frameList(ia).getSize ();
      for (int jn = 0; jn < nFrame; ++jn)
         wik(iAniso, mapAniso(frameList(ia)(jn))) = 1./double(nFrame);
   }

   //cout << "wik " << wik << endl;
   //cout << "wik ev " << wik.eigenvalues () << endl;
   //cout << wik.eigenvectors ().colRef(nAniso - 1) << endl;
   SxMatrix<Double> changeQ = -wik.transpose ();
   for (int i = 0; i < nAniso; ++i) changeQ(i,i) += 1.;

   //cout << changeQ << endl;
   //cout << "changeQ: " << changeQ.eigenvalues () << endl;

   // --- calculate frame charge tensors
   QFrame.resize (nAniso);
   for (int i = 0; i < 3; ++i)  {
      for (int j = 0; j < 3; ++j)  {
         SxVector<Double> tensorQ(nAniso);
         for (int iAniso = 0; iAniso < nAniso; ++iAniso)
            tensorQ(iAniso) = dQ(iAniso)(i,j);
         tensorQ = changeQ.solve (tensorQ);
         for (int iAniso = 0; iAniso < nAniso; ++iAniso)
            QFrame(iAniso)(i,j) = tensorQ(iAniso);
      }
   }
   for (int iAniso = 0; iAniso < nAniso; ++iAniso)  {
      SxMatrix3<Double> Q = QFrame(iAniso);
      for (int jAniso = 0; jAniso < nAniso; ++jAniso)  {
         Q -= wik(jAniso, iAniso) * QFrame(jAniso);
      }
      cout << "ia = " << (idAniso(iAniso)+1) << endl;
      cout << "Q real  = " << Qab(idAniso(iAniso)) << endl;
      cout << "Q aniso = " << dQ(iAniso) << endl;
      cout << "Q frame = " << QFrame(iAniso) << endl;
      cout << "dQ test = " << (dQ(iAniso) - Q) << endl;
   }
   if (QabBulk.getSize () == Qab.getSize ())  {
      SxArray<SxMatrix3<Double> > QFrameBulk(nAniso);
      SX_LOOP2(i,j)  {
         SxVector<Double> tensorQ(nAniso);
         SX_LOOP(iAniso) {
            int ia = idAniso(iAniso);
            tensorQ(iAniso) = QabBulk(ia)(i,j);
            if (i==j) tensorQ(iAniso) -= QabBulk(ia).trace ()/3.;
         }
         tensorQ = changeQ.solve (tensorQ);
         SX_LOOP(iAniso)
            QFrameBulk(iAniso)(i,j) = tensorQ(iAniso);
      }
      SX_LOOP(iAniso)  {
         int ia = idAniso(iAniso);
         if (genId(ia) < 0 &&
             (QFrame(iAniso) - QFrameBulk(iAniso)).absSqr ().sum () > 1e-10)
         {
            cout << "Modified anisotropy for bulk atom " << (ia+1)
                 << endl;
         }
      }

   }
}

SxMatrix<Double> 
SxHarmonicPotential::getMatrixConstraint (const SxList<int> &matSymId,
                                          bool symmetric) const
{
   int nSymM = int(matSymId.getSize ());
   SxMatrix<Double> U;
   if (symmetric)  {
      U.reformat (9,6);
      U.set (0.);
      // free parameters: diagonal elements
      for (int i = 0; i < 3; ++i) U(i+3*i,i) = 1.;
      // free parameters: off-diagonal elements
      U(1+2*3,3) = U(2+1*3,3) = sqrt(0.5);
      U(0+2*3,4) = U(2+0*3,4) = sqrt(0.5);
      U(0+1*3,5) = U(1+0*3,5) = sqrt(0.5);
   } else {
      U.reformat (9,9);
      SX_LOOP2(i,j) U(i,j) = (i==j) ? 1. : 0.;
   }
   SxMatrix<Double> M(9,9);
   for (int iSymM = 0; iSymM < nSymM; ++iSymM)  {
      M.set (0.);
      int iSym = matSymId(iSymM);
      SymMat S = syms(iSym);
      for (int j1 = 0; j1 < 3; ++j1)  {
         for (int k1 = 0; k1 < 3; ++k1)  {
            // delta term
            M(j1 + 3 * k1, j1 + 3 * k1) += 2.;

            for (int j2 = 0; j2 < 3; ++j2)  {
               for (int k2 = 0; k2 < 3; ++k2)  {
                  M(j1 + 3 * k1, j2 + 3 * k2) -= S(j1, j2) * S(k1, k2)
                                               + S(j2, j1) * S(k2, k1);
               }
            }
         }
      }
      SxMatrix<Double>::Eigensystem eig = SxMatrix<Double>(U.adjoint () ^ M ^ U).eigensystem ();
      eig.vecs = U ^ eig.vecs;
      int nDof = 0;
      for (int i = 0; i < eig.vecs.nCols (); ++i)
         if (fabs (eig.vals(i).re) < 1e-9)
            U.colRef (nDof++) <<= eig.vecs.colRef (i);
      if (nDof == 0) {
         cout << "Symmetry catastrophy: no degrees of freedom left..." << endl;
         SX_EXIT;
      }
      U.reshape (U.getSize ());
      U.resize (9 * nDof, true);
      U.reshape (9, nDof);
   }
   if (nSymM > 1)
      cout << "Computing symmetry constraints from " << nSymM << " symmetries: "
           << U.nCols () << " free parameters" << endl;
   //cout << M.eigenvalues ().real () << endl;
   return U;
}

void SxHarmonicPotential::shiftHesseToNeighbor (int ia, int in)
{
   int nSym = int(syms.getSize ());
   int ja = neighborIdx(ia)(in);
   int nOld = int(neighborIdx(ja).getSize ());

   // avoid duplicating this interaction at ja
   bool duplicate = false;
   for (int jn = 0; jn < nOld; ++jn)  {
      if (neighborIdx(ja)(jn) == ia)  {
         duplicate=true;
         int iGen = genNeighborId(ia)(in),
             jGen = genNeighborId(ja)(jn);
         if (iGen != jGen)  {
            cout << "WARNING: duplicate interactions arise from different "
                 << "origins: " << (iGen+1) << "<=>" << (jGen+1) << endl;
         }
         break;
      }
   }
   if (!duplicate)  {
      int jn;
      for (jn = 0; jn < nOld; ++jn)  {
         if (genNeighborId(ja)(jn) < genNeighborId(ia)(in)) continue;
         if (genNeighborId(ja)(jn) > genNeighborId(ia)(in)) break;
         if (genSymId(ja)(jn) % nSym < genSymId(ia)(in)) continue;
         break;
      }
      // insert new neighbor at jn
      cout << "Relocating interaction from " << (ia+1) << ":"
           << (in+1) << " to " << (ja+1) << ":" << (jn+1) << endl;

      genNeighborId(ja).resize (nOld + 1, true);
      genSymId(ja).resize (nOld + 1, true);

      neighborIdx(ja).resize (nOld + 1, true);
      hesse(ja).resize (nOld + 1, true);
      dist0(ja).resize (nOld + 1, true);

      // --- shift elements >= jn upwards
      for (int kn = nOld; kn > jn; kn--)  {
         genNeighborId(ja)(kn) = genNeighborId(ja)(kn-1);
         genSymId(ja)(kn) = genSymId(ja)(kn-1);
         neighborIdx(ja)(kn) = neighborIdx(ja)(kn-1);
         hesse(ja)(kn) = hesse(ja)(kn-1);
         dist0(ja)(kn) = dist0(ja)(kn-1);
      }

      // --- insert ia as neighbor jn to atom ja
      genNeighborId(ja)(jn) = genNeighborId(ia)(in);
      if (genSymId(ia)(in) >= 0)
         genSymId(ja)(jn) = genSymId(ia)(in) + nSym; // signalize inverse
      else
         genSymId(ja)(jn) = -2;
      neighborIdx(ja)(jn) = ia;
      hesse(ja)(jn) = hesse(ia)(in).transpose ();
      dist0(ja)(jn) = -dist0(ia)(in);
   } else {
      cout << "Removing duplicate interaction " << (ia+1) << ":"
           << (in+1) << " to " << (ja+1) << endl;
   }

   // --- remove neighbor in from atom ia
   for (int kn = in + 1; kn < neighborIdx(ia).getSize (); ++kn)  {
      genNeighborId(ia)(kn-1) =  genNeighborId(ia)(kn);
      genSymId(ia)(kn-1) = genSymId(ia)(kn);
      neighborIdx(ia)(kn-1) = neighborIdx(ia)(kn);
      hesse(ia)(kn-1) = hesse(ia)(kn);
      dist0(ia)(kn-1) = dist0(ia)(kn);
   }
   int nNew = int(neighborIdx(ia).getSize ()) - 1;
   genNeighborId(ia).resize (nNew, true);
   genSymId(ia).resize (nNew, true);
   neighborIdx(ia).resize (nNew, true);
   hesse(ia).resize (nNew, true);
   dist0(ia).resize (nNew, true);
}

void SxHarmonicPotential::setExtraForce (const SxAtomicStructure &f,
                                       const SxAtomicStructure &x)
{
   if (f.getSize () == 0)  {
      extraForce = f;
      return;
   }
   extraForce.copy (f);
   if (x.getNAtoms () > 0)  {
      SX_CHECK (refStr.getSize () > 0);
      SxGrid grid(refStr, 10);
      SxConstPtr<SxAtomInfo> map = refStr.match (grid, x);
      extraForce.replaceInfo (map);
   }
}
                  

SxHarmonicPotential::~SxHarmonicPotential ()
{
   // empty
}

void SxHarmonicPotential::init (const SxAtomicStructure &,
                                const SxSymbolTable *)
{
   SX_EXIT;
}


void SxHarmonicPotential::setupQab ()
{
   SX_CHECK (Qab.getSize () == 0);
   SxMatrix3<Double> one (1.,0.,0., 0.,1.,0., 0., 0., 1.);
   /*
   epsInv = (1./dielecConstant) * one;
   dielecConstant = 1.;
   */

   Qab.resize (refStr.getNAtoms ());
   for (int is = 0, iTlAtom = 0; is < getNSpecies (); ++is)  {
      for (int ia = 0; ia < refStr.getNAtoms (is); ++ia, ++iTlAtom)  {
         Qab(iTlAtom) = valenceCharge(is) * one;
      }
      if (fabs(valenceCharge(is)) > 1e-12) valenceCharge(is) = 1.;
   }
}

SxAtomicStructure
SxHarmonicPotential::getMu (const SxAtomicStructure &tau) const
{
   int nAniso = (int)idAniso.getSize ();
   if (!nAniso) return SxAtomicStructure ();
   SxAtomicStructure mu(nAniso);
   for (int iAniso = 0; iAniso < nAniso; ++iAniso)  {
      Coord d = tau(idAniso(iAniso));
      ssize_t nFrame = frameAtoms(iAniso).getSize ();
      for (int iF = 0; iF < nFrame; ++iF)
         d -= tau(frameAtoms(iAniso)(iF)) / double(nFrame);
      d -= framePositions(iAniso);
      mu.ref(iAniso) = QFrame(iAniso) ^ d;
   }
   cout << "mu=" << mu << endl;
   return mu;
}

SxAtomicStructure
SxHarmonicPotential::getElForces (const SxAtomicStructure &tau)
{
   SxAtomicStructure forces = tau.getNewStr ();
   if (noCoulomb)  {
      energy = 0.;
      forces.set (Coord(0.,0.,0.));
      return forces;
   }
   int nAniso = (int)idAniso.getSize ();
   if (nAniso >= 0)  {
      SxGBasis &G = *gBasis;
      G.changeTau (tau);
      
      // --- calculate dipoles
      SxAtomicStructure mu = getMu (tau);

      // --- set up rho
      PsiG rho(G), V(G), muG, gradV;
      rho.set (0.);
      muG = SxDiracMat<Complex16> (G.ng, 3);
      muG.setBasis (G);
      muG.set (0.);
      for (int is = 0, iTl = 0; is < tau.getNSpecies (); ++is)  {
         for (int ia = 0; ia < tau.getNAtoms (is); ++ia, ++iTl)  {
            // no charge tensor
            int iAniso = mapAniso(iTl);
            if (fabs(qIso(iTl)) < 1e-12 && iAniso < 0) continue;

            PsiG T = G.getPhaseFactors(is,ia);

            // isotropic charge
            rho.plus_assign_ax (qIso(iTl), T);

            // anisotropy
            if (iAniso >= 0)  {
               for (int i = 0; i < 3; ++i)
                  muG.colRef (i).plus_assign_ax (mu(iAniso)(i), T);
            }
         }
      }
      if (fabs(rho(0).re) > 1e-8)  {
         cout << "Warning: net charge=" << rho(0).re << endl;
      }
      //rho(0) = 0.;

      // correction for using betaSoft instead of beta come below
      SxDiracVec<Double> shape = exp(-G.g2*(betaSoft*betaSoft/4.))
                               / sqrt(tau.cell.volume);
      rho *= shape;
      for (int i = 0; i < 3; ++i)
         muG.colRef(i) *= shape;

      SxIdx noZero(1,G.ng-1);
      V = rho.getCopy ();
      for (int i = 0; i < 3; ++i)
         V.plus_assign_ax (-I, muG.colRef (i) * G.gVec.colRef(i));
      V *= (FOUR_PI / dielecConstant);
      V(noZero) /= G.g2(noZero);
      V(0) = -rho(0).re * FOUR_PI * betaSoft*betaSoft/4. / dielecConstant;

      gradV = SxDiracMat<Complex16> (G.ng, 3);
      gradV.setBasis (G);
      for (int i = 0; i < 3; ++i)  {
         gradV.colRef(i) <<= (V * G.gVec.colRef(i));
         //gradV(0,i) = -I * muG(0,i)/3.
         //             * FOUR_PI/dielecConstant; // G=0 spherical average
      }
      gradV *= I;

      double monopoleSelf = 1./(betaSoft*sqrt(TWO_PI)) / dielecConstant;
      double dipoleSelf = 1./(3.*betaSoft*betaSoft*betaSoft*sqrt(TWO_PI)*dielecConstant);
      energy = 0.5 * dot (rho, V).re            // charge energy
             - qIso.normSqr () * monopoleSelf;  // charge self-interaction
      if (nAniso > 0)  {
         energy += 0.5 * dot (muG, gradV).re        // dipole energy
                 - mu.absSqr ().sum () * dipoleSelf;// dipole self-interaction
      }

      // this is the G^2 term of rho with the 1/G^2 of V
      energy += 0.5 * rho(0).re * V(0).re;

      sxprintf ("energySoft=%.16f\n",energy);
      
      // --- Gaussian width correction
      SxGrid grid (tau, 3);
      SxNeighbors nn;
      int neighborParams = SxNeighbors::StoreRel
                         | SxNeighbors::StoreIdx
                        // | SxNeighbors::IncludeZeroDistance
                         ;
      double rcut = 13. * max(beta, betaSoft);
      SxYlm::SxClebschTable cg = SxYlm::getClebschGordan (2,3,2,SxYlm::RealYlm);
      double beta2 = beta*beta, betap2 = betaSoft * betaSoft;
      double energy2 = 0.;


      // --- FORCES + Gaussian width correction

      // force = - Q(R) [d phi_R(G)/dR] shape(G) V(G) 
      //         - mu(R,i) [d phi_R(G)/ d R] shape(G) gradV(G,i)
      //         - [d mu(R',i) / d R] phi_R(G) shape(G) gradV(G,i)
      //         + [d mu(R',i) / d R] mu(R', i) * dipoleSelf;

      forces.set (Coord(0.,0.,0.));

      SxGaussIJ gaussIJ;
      for (int is = 0, iTl = 0; is < tau.getNSpecies (); ++is)  {
         for (int ia = 0; ia < tau.getNAtoms (is); ++ia, ++iTl)  {
            int iAniso = mapAniso(iTl);
            // no charge tensor
            if (fabs(qIso(iTl)) < 1e-12 && iAniso < 0) continue;

            PsiG Ts = G.getPhaseFactors(is,ia) * shape;
            Coord dEdMu(0., 0., 0.);

            // - Q(R) [d phi_R(G)/dR] shape(G) V(G) 
            for (int iDir = 0; iDir < 3; ++iDir)  {
               forces.ref(is,ia)(iDir) 
                  += qIso(iTl) * dot(Ts,  V * G.gVec.colRef(iDir)).im;
            }

            // --- now correct for using betaSoft instead of beta...
            nn.compute (grid, tau, tau.getAtom (iTl), rcut, neighborParams);
            for (int ja = 0; ja < nn.relPositions.getNAtoms(); ++ja)  {
               int jTl = nn.idx (ja);
               //int js = tau.getISpecies (jTl);
               Coord r = nn.relPositions.getAtom(ja);
               gaussIJ.setDelta (r, 2. * beta2, 2. * betap2, 3);
               SxMatrix<Double> hIJ = gaussIJ.compute (1, 1, cg);
               // note: need a prefactor of (2L + 1)!!/sqrt(2L+1) due to
               // normalization convention of the gaussIJ kernel for expressions
               // involving G^L
               hIJ *= 1. / FOUR_PI / dielecConstant;

               Coord f;

               // monopole-monopole correction
               energy2 += 0.5 * qIso(iTl) * qIso(jTl) * hIJ(0,0);
               for (int i = 1; i < 4; i++)
                  f(i % 3) = 0.5 * sqrt(3.) * qIso(iTl) * qIso(jTl) * hIJ(i,0);

               if (iAniso >= 0)  {
                  // dipole-monopole (note: lm ordering is y z x)
                  for (int i = 1; i < 4; i++)  {
                     energy2 -= sqrt(3.) * qIso(jTl) * mu(iAniso)(i % 3)
                              * hIJ(i,0);
                     dEdMu(i % 3) -= sqrt(3.) * qIso(jTl) * hIJ(i,0);
                     for (int j = 1; j < 4; j++)  {
                        f(j % 3) += 3. * qIso(jTl) * mu(iAniso)(i % 3)
                                  * hIJ(i,j);
                     }
                  } 
                  if (mapAniso(jTl) >= 0)  {
                     // --- dipole-dipole energy
                     int jAniso = mapAniso(jTl);
                     for (int i = 1; i < 4; i++)  {
                        for (int j = 1; j < 4; j++)  {
                           double dE = hIJ(i,j) * mu(jAniso)(j % 3);
                           energy2 += 1.5 * dE * mu(iAniso)(i % 3);
                           dEdMu(i % 3) += 3. * dE;
                        }
                     }
                     // --- dipole-dipole forces
                     gaussIJ.compute (1, 1, cg, 1); // prepare forces
                     for (int i = 1; i < 4; i++)  {
                        for (int j = 1; j < 4; j++)  {
                           f += 1.5 / FOUR_PI / dielecConstant
                                * mu(jAniso)(j % 3) * mu(iAniso)(i % 3) 
                                * gaussIJ.getForce (i, j, cg);
                        }
                     }
                  }
               }

               forces.ref(iTl) += f;
               forces.ref(jTl) -= f;
            }

            if (iAniso >= 0)  {
               // - mu(R,i) [d phi_R(G)/ d R] shape(G) gradV(G,i)
               PsiG muGradV = mu(iAniso)(0) * gradV.colRef(0)
                            + mu(iAniso)(1) * gradV.colRef(1)
                            + mu(iAniso)(2) * gradV.colRef(2);
               for (int iDir = 0; iDir < 3; ++iDir)  {
                  forces.ref(is,ia)(iDir) 
                     += dot(Ts,  muGradV * G.gVec.colRef(iDir)).im;
               }
               // - [d mu(R',i) / d R] phi_R(G) shape(G) gradV(G,i)
               for (int iDir = 0; iDir < 3; ++iDir)  {
                  dEdMu(iDir) += dot(Ts, gradV.colRef(iDir)).re;
               }
               // + [d mu(R',i) / d R] mu(R', i) * dipoleSelf;
               dEdMu -= 2. * mu(iAniso) * dipoleSelf;

               // d mu(R',i) / d R(,j) = Q(R')_ij (delta(R,R') - w_Frame(R',R))
               Coord dEdFrame = QFrame(iAniso).transpose () ^ dEdMu;

               forces.ref(is,ia) -= dEdFrame;
               ssize_t nFrame = frameAtoms(iAniso).getSize ();
               for (int iF = 0; iF < nFrame; ++iF)  {
                  int jAtom = frameAtoms(iAniso)(iF);
                  forces.ref (jAtom) += dEdFrame / double(nFrame);
               }

            }
         }
      }

      energy += energy2;
      sxprintf ("energy=%.16f\n", energy);

      return forces;
   }
   SX_EXIT;

   forces.set (Coord(0.,0.,0.));
   energy = 0.;
   
   // --- electrostatic part
   SxGBasis &G = *gBasis;
   G.changeTau (tau);
   PsiG rho(G), V(G);
   rho.set (0.);
   for (int is = 0; is < getNSpecies (); ++is)  {
      rho.plus_assign_ax (valenceCharge(is), G.structureFactors(is));
   }
   if (fabs(rho(0).re) > 1e-8)  {
      cout << "Warning: net charge=" << rho(0).re << endl;
   }
   rho(0) = 0.;
   SxDiracVec<Double> shape = exp(-G.g2*(beta*beta/4.)) / sqrt(tau.cell.volume);
   rho *= shape;

   SxIdx noZero(1,G.ng-1);
   V(noZero) = FOUR_PI/dielecConstant
             * rho(noZero)/G.g2(noZero);
   V(0) = 0.;
   VALIDATE_VECTOR(V);

   energy += 0.5 * dot (rho, V).re;

   // --- external potential
   if (vExt.getSize () > 0)  {
      SX_CHECK (vExt.getBasisPtr () == &G);
      energy += dot(rho, vExt).re;
      V += vExt;
   }

   // --- electrostatic and external forces
   SxDiracMat<Complex16> VG(G.ng, 3);
   for (int idir = 0; idir < 3; ++idir)  {
      VG.colRef(idir) <<= G.gVec.colRef(idir) * V;
   }
   for (int is = 0; is < tau.getNSpecies (); ++is)  {
      for (int ia = 0; ia < tau.getNAtoms (is); ++ia)  {
         PsiG charge = shape * G.getPhaseFactors(is,ia);
         for (int idir = 0; idir < 3; ++idir)  {
            forces.ref(is,ia)(idir) -= valenceCharge(is) 
                                     * dot(VG.colRef(idir), charge).im;
         }
      }
   }
   return forces;
}


SxAtomicStructure
SxHarmonicPotential::getForces (const SxAtomicStructure &tau,
                                const SxSymbolTable *)
{
   SX_CHECK (tau.getSize () == neighborIdx.getSize (),
             tau.getSize (), neighborIdx.getSize ());
   // electrostatic part
   SxAtomicStructure forces = getElForces (tau);
   cout << "El. Force drift: " << forces.sum () << endl;

   // --- calculate local interactions
   for (int ia = 0; ia < tau.getNAtoms (); ++ia)  {
      for (int in = 0; in < neighborIdx(ia).getSize (); ++in)  {
         int idN = neighborIdx(ia)(in);
         Coord x = tau(ia) - tau(idN) + dist0(ia)(in);
         Coord df = hesse(ia)(in) ^ x;
         /*
         if (x.normSqr () > 1e-10)  {
            cout << (ia) << "<->" << idN << endl;
            cout << hesse(ia)(in) << endl;
            cout << df << endl;
         }
         */

         // energy calculations
         energy += 0.5 * (x ^ df);

         // forces
         forces.ref(ia)  -= df;
         forces.ref(idN) += df;
      }
   }

   if (extraForce.getSize () > 0)  {
      cout << "Adding external forces" << endl;
      forces += extraForce;
      energy -= dot(extraForce.coordRef (), (tau - refStr).coordRef ());
   }

   cout << "Force drift: " << forces.sum () << endl;
   sxprintf ("energy = %.16f\n", energy);
   return forces;
}

PrecEnergy SxHarmonicPotential::getEnergy () const
{
   return energy;
}


/// Set extra potential
void SxHarmonicPotential::setExtraPot (const SxMeshR &vR, bool screened)
{
   if (Qab.getSize () > 0)  {
      cout << "External potential not implemented for anisotropic electrostatics" << endl;
      SX_EXIT;
   }
   vExt = *gBasis | vR;
   if (!screened)  {
      vExt /= dielecConstant;
   }
}


SxMatrix<Double> 
SxHarmonicPotential::paramGrad (const SxAtomicStructure &tau)
{
   SX_CHECK (tau.getSize () == neighborIdx.getSize (),
             tau.getSize (), neighborIdx.getSize ());
   // --- gradient from local interactions
   int nParam = int(genNeighbors.getSize ()) * 9;
   SxMatrix<Double> res (3 * tau.getNAtoms (), nParam);
   res.set (0.);

   //cout << SX_SEPARATOR;
   for (int ia = 0; ia < tau.getNAtoms (); ++ia)  {
      for (int in = 0; in < neighborIdx(ia).getSize (); ++in)  {
         int idN = neighborIdx(ia)(in);
         Coord x = tau(ia) - tau(idN) + dist0(ia)(in);
         if (x.normSqr () < 1e-12) continue;

         int iGen = genNeighborId(ia)(in);

         //cout << "iGen=" << iGen << " " << (ia+1) << "/"
         //     << (idN+1) << ": ";

         //Coord df = hesse(ia)(in) ^ x;
         // df(a)/d H(b,c) = delta(a,b) x(c)
         // df(a)/ d G(d,e) = sum_b,c delta(a,b) x(c) S^t(d,b) S(c,e)
         //                 = sum_c S(a,d) S^t(e,c) x(c)

         int iSym = genSymId(ia)(in);
         if (iSym >= 0)  {
            SymMat S = syms(iSym);
            x  = S.transpose () ^ x;
            //cout << x << endl;
            //double x0 = x(0), x1 = x(1), x2 = x(2);
            for (int d = 0; d < 3; ++d)  {
               for (int e = 0; e < 3; ++e)  {
                  int iParam = 9 * iGen + 3 * e + d;
                  // => forces.ref(ia)  -= df;
                  for (int a = 0; a < 3; ++a)
                     res(3 * ia + a, iParam)  -= S(a,d) * x(e);
                  // => forces.ref(idN) += df;
                  for (int a = 0; a < 3; ++a)
                     res(3 * idN + a, iParam) += S(a,d) * x(e);
               }
               /*
               double S0 = S(0,d), S1 = S(1,d), S2 = S(2,d);
               double S0x0 = S0 * x0;
               double S0x1 = S0 * x1;
               double S0x2 = S0 * x2;
               double S1x0 = S1 * x0;
               double S1x1 = S1 * x1;
               double S1x2 = S1 * x2;
               double S2x0 = S2 * x0;
               double S2x1 = S2 * x1;
               double S2x2 = S2 * x2;
               int off0 = 3 * tau.getNAtoms () * (9 * iGen + d),
                   off1 = 3 * tau.getNAtoms () * (9 * iGen + d + 3),
                   off2 = 3 * tau.getNAtoms () * (9 * iGen + d + 6);
               res(3 * ia + 0 + off0) -= S0x0;
               res(3 * ia + 1 + off0) -= S1x0;
               res(3 * ia + 2 + off0) -= S2x0;
               res(3 * ia + 0 + off1) -= S0x1;
               res(3 * ia + 1 + off1) -= S1x1;
               res(3 * ia + 2 + off1) -= S2x1;
               res(3 * ia + 0 + off2) -= S0x2;
               res(3 * ia + 1 + off2) -= S1x2;
               res(3 * ia + 2 + off2) -= S2x2;
               res(3 * idN + 0 + off0) += S0x0;
               res(3 * idN + 1 + off0) += S1x0;
               res(3 * idN + 2 + off0) += S2x0;
               res(3 * idN + 0 + off1) += S0x1;
               res(3 * idN + 1 + off1) += S1x1;
               res(3 * idN + 2 + off1) += S2x1;
               res(3 * idN + 0 + off2) += S0x2;
               res(3 * idN + 1 + off2) += S1x2;
               res(3 * idN + 2 + off2) += S2x2;
               */
            }
         } else {
            //cout << x << endl;
            // S = E
            for (int e = 0; e < 3; ++e)  {
               int iParam = 9 * iGen + 3 * e /* + a */;
               for (int a = 0; a < 3; ++a)  {
                  // => forces.ref(ia)  -= df;
                  res(3 * ia + a, iParam + a)  -= x(e);
                  // => forces.ref(idN) += df;
                  res(3 * idN + a, iParam + a) += x(e);
               }
            }
         }
      }
   }
   return res;

}


SxMatrix<Double>
SxHarmonicPotential::paramQ (int iTlAtom, const Coord &field) const
{
   SX_CHECK (iTlAtom >= 0 && iTlAtom < genQ.getSize (), iTlAtom,
             genQ.getSize ());
   int g = genQ(iTlAtom);
   int iSym = genQSymId(iTlAtom);
   SxMatrix<Double> res(9,3);
   if (iSym == -1)  {
      res.set (0.);
      // S is identity operation => b=c, a=d
      for (int a = 0; a < 3; a++)
         for (int b = 0; b < 3; b++)
            res(3 * b + a, a) = field(b);
   } else {
      SymMat S = syms(iSym);
      Coord rotE = field ^ S; // \sum b S_bc E_b
      SX_LOOP3(a,c,d)
         res(3 * c +  d, ssize_t(a)) = S(a,d) * rotE(c);
   }
   return matSymQ(g).transpose () ^ res;
}

void SxHarmonicPotential::setHesse (const SxVector<Double> &param)
{
   SX_CHECK (refStr.getSize () == neighborIdx.getSize (),
             refStr.getSize (), neighborIdx.getSize ());
   SX_CHECK (param.getSize () == 9 * genNeighbors.getSize (),
             param.getSize (), genNeighbors.getSize ());

   SxArray<SxMatrix3<Double> > hRef(genNeighbors.getSize ());
   SX_LOOP3(iGen, i(3), j(3)) hRef(iGen)(i,j) = param(9 * iGen + 3 * j + i);

   for (int ia = 0; ia < refStr.getNAtoms (); ++ia)  {
      for (int in = 0; in < genNeighborId(ia).getSize (); ++in)  {
         int iGen = genNeighborId(ia)(in);
         int iSym = genSymId(ia)(in);
         if (iSym >= 0)  {
            const SymMat &S = syms(iSym);
            hesse(ia)(in) = S ^ hRef(iGen) ^ S.transpose ();
         } else {
            hesse(ia)(in) = hRef(iGen);
         }
      }
   }
   if (extraForce.getSize () == 0)  {
      SxAtomicStructure f = getForces (refStr);
      if (f.coords.normSqr () > 1e-12 * f.getNAtoms ())
         extraForce = -f;
   }
}

inline ssize_t idx21(ssize_t a, ssize_t b)
{
   SX_CHECK (a >= 0 && a < 3, a);
   SX_CHECK (b >= 0 && b < 3, b);
   return (a==b ? a : (2+a+b));
}

enum HesseTimer { GaussIJ, HesseGIJ, HesseNN, HesseRS, WabTime, TiTime, WjPrepare,
                  HesseEl, HesseShort,  HesseTotal};
SX_REGISTER_TIMERS (HesseTimer)
{
   regTimer (GaussIJ, "gaussIJ");
   regTimer (HesseGIJ, "Hesse gaussIJ");
   regTimer (HesseNN, "Hesse NN");
   regTimer (HesseRS, "Hesse RS");
   regTimer (WabTime, "Hesse Wab(k)");
   regTimer (WjPrepare,  "Wj");
   regTimer (HesseEl,  "Hesse el.");
   regTimer (HesseShort,  "Hesse short-range");
   regTimer (HesseTotal,  "Hesse total");
}

SxMatrix<Complex16> 
SxHarmonicPotential::getHesseEl (const Coord &kVec,
                                 const SxAtomicStructure &tau)
{
   if ((tau.cell - refStr.cell).absSqr ().sum () > 1e-10)  {
      cout << endl << SX_SEPARATOR;
      cout << "Invalid cell: " << tau.cell << endl;
      cout << "Should be: " << refStr.cell << endl;
      SX_EXIT;
   }
   int nDof = 3 * tau.getNAtoms ();
   SxMatrix<Complex16> res(nDof, nDof);
   res.set (0.);
   if (noCoulomb) return res;

   // --- electrostatic part
   SX_CLOCK (HesseEl);
   Coord kVecM = tau.cell.getReciprocalCell ()
                 .getMapped (kVec, SxCell::WignerSeitz);
   bool gammaPoint = (kVec.normSqr () < 1e-20);
   if (kVecM.normSqr () < 0.999 * kVec.normSqr () && !gammaPoint)  {
      cout << "WARNING: k-point " << kVec 
           << "outside first Brillouin zone" << endl;
   }
   if (gammaPoint)  {
      cout << "Warning: computing Hesse matrix for Gamma point without LO-TO splitting"
           << endl;
   }

   int nAniso = (int)idAniso.getSize ();

   SxGBasis &G = *gBasis;
   G.changeTau (tau);
   
   // --- calculate dipoles
   SxAtomicStructure mu;
   if (nAniso > 0) mu = getMu (tau);

   // --- set up rho
   PsiG rho(G), V(G), muG, gradV;
   rho.set (0.);
   muG = SxDiracMat<Complex16> (G.ng, 3);
   muG.setBasis (G);
   muG.set (0.);
   for (int is = 0, iTl = 0, iAniso = 0; is < tau.getNSpecies (); ++is)  {
      for (int ia = 0; ia < tau.getNAtoms (is); ++ia, ++iTl)  {
         // no charge tensor
         if (fabs(qIso(iTl)) < 1e-12 && mapAniso(iTl) == -1) continue;

         PsiG T = G.getPhaseFactors(is,ia);

         // isotropic charge
         rho.plus_assign_ax (qIso(iTl), T);

         // anisotropy
         if (mapAniso(iTl) == iAniso)  {
            for (int i = 0; i < 3; ++i)
               muG.colRef (i).plus_assign_ax (mu(iAniso)(i), T);
            iAniso++;
         }
      }
   }
   bool chargedSystem = fabs(rho(0).re) > 1e-8;
   bool chargeCorrection = true;
   if (chargedSystem)  {
      cout << "Warning: net charge=" << rho(0).re << endl;
      if (!gammaPoint && chargeCorrection)  {
         cout << "Using charge corrections." << endl;
      }
   }
   //rho(0) = 0.;

   // correction for using betaSoft instead of beta come below
   // note: we actually include the Gaussian factor of the charge here, too.
   // so our shape here is the square of the single Gaussian shape
   SxDiracVec<Double> shape = exp(-G.g2*(betaSoft*betaSoft/2.))
                            / tau.cell.volume;
   rho *= shape;
   for (int i = 0; i < 3; ++i)
      muG.colRef(i) *= shape;

   SxIdx noZero(1,G.ng-1);
   V = rho.getCopy ();
   for (int i = 0; i < 3; ++i)
      V.plus_assign_ax (-I, muG.colRef (i) * G.gVec.colRef(i));
   V *= (FOUR_PI / dielecConstant);
   V(noZero) /= G.g2(noZero);
   V(0) = -rho(0).re * FOUR_PI * betaSoft*betaSoft/2. / dielecConstant;

   gradV = SxDiracMat<Complex16> (G.ng, 3);
   gradV.setBasis (G);
   for (int i = 0; i < 3; ++i)  {
      gradV.colRef(i) <<= (V * G.gVec.colRef(i));
      //gradV(0,i) = -I * muG(0,i)/3.
      //             * FOUR_PI/dielecConstant; // G=0 spherical average
   }
   gradV *= I;

   double dipoleSelf = 1./(3.*betaSoft*betaSoft*betaSoft*sqrt(TWO_PI)*dielecConstant);

   SxDiracVec<Complex16> Gk;
   Gk.reformat (G.ng, 3);
   Gk.setBasis (G);
   SX_LOOP(i)
      Gk.colRef(i) <<= G.gVec.colRef(i) + kVec(i);
   
   // --- Gaussian width correction
   SxGrid grid (tau, 3);
   SxNeighbors nn;
   int neighborParams = SxNeighbors::StoreRel
                      | SxNeighbors::StoreIdx
                     // | SxNeighbors::IncludeZeroDistance
                      ;
   double rcut = 13. * max(beta, betaSoft);
   SxYlm::SxClebschTable cg = SxYlm::getClebschGordan (3,4,3,SxYlm::RealYlm);
   double beta2 = beta*beta, betap2 = betaSoft * betaSoft;

   // --- G+k - dependent vectors: G+k, |G+k|^2, exp(-|G+k|)^2 etc.
   SxDiracVec<Double> gk2 = kVec.normSqr () + G.g2;
   SX_LOOP(iDir) gk2 += 2. * G.gVec.colRef(iDir) * kVec(iDir);

   SxDiracVec<Double> shapeK2 = exp(-gk2*(betaSoft*betaSoft/2.))
                              *(FOUR_PI / (tau.cell.volume * dielecConstant));
   shapeK2(noZero) /= gk2(noZero);
   if (!gammaPoint) shapeK2(0) /= gk2(0);
   else             shapeK2(0) = 0.;

   SxDiracVec<Double> shapeK2ab(G, 6);
   for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
         shapeK2ab.colRef (idx21(a,b)) <<= shapeK2
                                         * Gk.colRef(a) * Gk.colRef(b);

   // --- we block on i and j for performance in the algebra (larger matrices)
   int iBlockSize = 64 * 2;
   int jBlockSize = 64;
   for (int iBlock = 0; iBlock * iBlockSize < tau.getNAtoms (); ++iBlock)  {
      int iFrom = iBlock * iBlockSize;
      int iNext = min(iFrom + iBlockSize, tau.getNAtoms ());

      // --- collect phase factors for the current i-block
      PsiGI allTi;
      allTi.reformat (G.ng, iNext - iFrom);
      for (int iTl = iFrom; iTl < iNext; ++iTl)  {
         int ia, is = tau.getISpecies (iTl, &ia);
         allTi.colRef(iTl - iFrom) <<= G.getPhaseFactors (is,ia);
      }

      for (int jBlock = 0; jBlock * jBlockSize < tau.getNAtoms (); ++jBlock)  {
         int jFrom = jBlock * jBlockSize;
         int jNext = min(jFrom + jBlockSize, tau.getNAtoms ());

         // --- collect Wj = Tj * Gk(a) Gk(b) * shape2 for current j-block
         SX_START_TIMER (WjPrepare);
         PsiGI allWj;
         allWj.reformat (G.ng, 6 * (jNext - jFrom));
         for (int jTl = jFrom; jTl < jNext; ++jTl)  {
            int ja, js = tau.getISpecies (jTl, &ja);
            PsiG Tj = G.getPhaseFactors(js,ja);
            for (int ab = 0; ab < 6; ++ab)
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
               for (int ig = 0; ig < G.ng; ++ig)
                  allWj(ig, 6 * (jTl - jFrom) + ab)
                     = Tj(ig) * shapeK2ab(ig, ab);
         }
         SX_STOP_TIMER (WjPrepare);

         // --- calculate W^{ij}_ab(k) for current i/j blocks
         SX_START_TIMER (WabTime);
         PsiG allWabJ = allTi.overlap (allWj);
         SX_STOP_TIMER (WabTime);

// <<<<<<<< shift indent to the left, now come the actual i/j loops

   SxGaussIJ gaussIJ;
   for (int js = 0, jTl = 0; js < tau.getNSpecies (); ++js)  {
      for (int ja = 0; ja < tau.getNAtoms (js); ++ja, ++jTl)  {

         // perform loop body only for current j-block
         if (jTl < jFrom || jTl >= jNext) continue;

         // no charge tensor
         if (fabs(qIso(jTl)) < 1e-12 && mapAniso(jTl) < 0) continue;
         
         if (tau.getNAtoms () > 100) (cout << ".").flush ();

         Coord muJ(0.,0.,0.);
         int jAniso = mapAniso(jTl);
         if (jAniso >= 0) muJ = mu(jAniso);

         // --- auxiliary data for the reciprocal part
         PsiG Tj = G.getPhaseFactors(js,ja);

         PsiG Wj = allWj(SxIdx (G.ng * 6 * (jTl - jFrom), 
                                G.ng * (6 * (jTl - jFrom + 1)) - 1) ); 
         Wj.reshape (G.ng, 6);
         int nr = int(allWabJ.nRows ());
         PsiG allWab = allWabJ(SxIdx (nr * 6 * (jTl - jFrom),
                                      nr * (6 * (jTl - jFrom + 1)) - 1) );
         allWab.reshape (nr, 6);

         // --- reciprocal space part
         for (int is = 0, iTl = 0; is < tau.getNSpecies (); ++is)  {
            for (int ia = 0; ia < tau.getNAtoms (is); ++ia, ++iTl)  {
               // perform loop body only for current j-block
               if (iTl < iFrom || iTl >= iNext) continue;

               // no charge tensor
               if (fabs(qIso(iTl)) < 1e-12 && mapAniso(iTl) < 0) continue;

               Coord muI(0.,0.,0.);
               int iAniso = mapAniso(iTl);
               if (iAniso >= 0) muI = mu(iAniso);

               //PsiG Ti = G.getPhaseFactors(is,ia);
               PsiG Ti = allTi.colRef (iTl - iFrom);

               //PsiG Wab = Ti.adjoint () ^ Wj;
               PsiG Wab = allWab.row (iTl - iFrom);
               // d^2 E/d R^k d mu^n
               PsiG d2EdRdMu(6); d2EdRdMu.set (0.);

               // Q/Q term of H-1 (Eq. 64 in implementation notes)
               SX_LOOP2(a(3), b(3))
                  res(3*iTl+a, 3*jTl+b) += qIso(iTl) * qIso(jTl)
                                           * Wab(idx21(a,b));

               // Q+mu/mu term of H-1 (Eq. 64 in implementation notes)
               if (muJ.normSqr () > 1e-12)  {
                  // muJ part: G(j) * mu(J)
                  PsiG TiMuJ(G.ng, 0.);
                  SX_LOOP(a) TiMuJ += muJ(a) * Gk.colRef(a);
                  if (muI.normSqr () > 1e-12)  {
                     // Q_i + I \sum_b G(b) mu_i(b)
                     PsiG QMuI(G.ng);
                     QMuI.set (qIso(iTl));
                     SX_LOOP(b) QMuI.plus_assign_ax (I*muI(b), Gk.colRef(b));
                     TiMuJ *= QMuI;
                  } else {
                     // Q_i
                     TiMuJ *= qIso(iTl);
                  }
                  TiMuJ *= Ti;
                  PsiG Wab2 = (TiMuJ.adjoint () ^ Wj);
                  SX_LOOP2(a(3),b(3))
                     res(3*iTl+a, 3*jTl+b) += I * Wab2(idx21(a,b));
               }
               // --- further mu/X terms in H-1 and H-2
               if (muI.normSqr () > 1e-12) {
                  PsiG TmuI(G.ng, 0.);
                  SX_LOOP(b) TmuI.plus_assign_ax (-I*muI(b), Gk.colRef(b));
                  TmuI *= Ti;
                  PsiG Wab2 = TmuI.adjoint () ^ Wj;
                  // mu/Q term of H-1 (Eq. 64 in implementation notes)
                  SX_LOOP2(a(3),b(3))
                     res(3*iTl+a, 3*jTl+b) -= qIso(jTl) * Wab2(idx21(a,b));
                  
                  // --- mu_k/\tilde Q^n term of H-2 (Eq. 65)
                  if (jAniso >= 0) d2EdRdMu += Wab2;
               }

               if (jAniso >= 0)  {
                  // Q^k term in Eq. (65) in implementation notes
                  d2EdRdMu += qIso(iTl) * Wab;
                  
                  addH23 (iTl, jTl, &res, d2EdRdMu, kVec);
               }

               // --- H-4 Eq. (66)
               if (iAniso >= 0 && jAniso >= 0)  {
                  SxMatrix3<Complex16> d2EdMu2;
                  SX_LOOP2(a2,b2) d2EdMu2(a2,b2) = Wab(idx21(a2,b2));
                  // dipole self energy, Eq. 70 on the diagonal 
                  if (iTl == jTl)  {
                     for (int aa = 0; aa < 3; ++aa)
                        d2EdMu2(aa,aa) -= 2. * dipoleSelf;
                  }
                  addH4 (iTl, jTl, &res, d2EdMu2, kVec);
               }
            }
         }
         if (iBlock > 0) continue;

         // --- diagonal (delta_kl in Eq. 64)
         
         // QV term of diagonal
         SX_LOOP2(a,b)  {
            PsiG Vab = G.gVec.colRef(a) * G.gVec.colRef (b) * V;
            res(3*jTl + a, 3*jTl + b) -= qIso(jTl) * dot (Tj, Vab);
         }
         if (chargedSystem && !gammaPoint && chargeCorrection)  {
            // charge correction
            Coord kNorm = kVec / kVec.norm ();
            double prefac = FOUR_PI * qIso(jTl) * rho(0).re / dielecConstant;
            SX_LOOP2(a,b)  {
               res(3*jTl + a, 3*jTl + b) -= prefac * kNorm(a) * kNorm(b);
            }
         }

         // mu gradV term of diagonal
         if (muJ.normSqr () > 1e-12)  {
            PsiG muE(G);
            muE.set (0.);
            SX_LOOP(iDir) muE += muJ(iDir) * gradV.colRef(iDir);
            SX_LOOP2(a,b)  {
               PsiG muEab = G.gVec.colRef(a) * G.gVec.colRef (b) * muE;
               res(3*jTl + a, 3*jTl + b) += dot (Tj, muEab);
            }
         }

         // --- wkl term in H-2 / H-3
         if (jAniso >= 0) {
            PsiG diagMu(6);
            for (int aa = 0; aa < 3; ++aa)
               diagMu(aa) = dot(Tj, G.gVec.colRef(aa) * gradV.colRef(aa));
            diagMu(3) = dot(Tj, G.gVec.colRef(0) * gradV.colRef(1));
            diagMu(4) = dot(Tj, G.gVec.colRef(0) * gradV.colRef(2));
            diagMu(5) = dot(Tj, G.gVec.colRef(1) * gradV.colRef(2));
            diagMu *= I;
            if (chargedSystem && !gammaPoint && chargeCorrection)  {
               // charge correction
               Coord kNorm = kVec / kVec.norm ();
               double prefac = FOUR_PI * rho(0).re / dielecConstant;
               for (int aa = 0; aa < 3; ++aa)
                  diagMu(aa) -= prefac * sqr(kNorm(aa));
               diagMu(3) -= prefac * kNorm(0) * kNorm(1); 
               diagMu(4) -= prefac * kNorm(0) * kNorm(2); 
               diagMu(5) -= prefac * kNorm(1) * kNorm(2); 
            }
               
            addH23 (jTl, jTl, &res, diagMu, kVec);

         }
         
         // --- real-space part (Gaussian width correction)
         // --- now correct for using betaSoft instead of beta...
         //SX_START_TIMER (HesseNN);
         nn.compute (grid, tau, tau.getAtom (jTl), rcut, neighborParams);
         //SX_STOP_TIMER (HesseNN);
         SX_CLOCK (HesseRS);
         for (int ia = 0; ia < nn.relPositions.getNAtoms(); ++ia)  {
            int iTl = nn.idx (ia);
            int iAniso = mapAniso(iTl);
            // no charge tensor?
            if (fabs(qIso(iTl)) < 1e-12 && iAniso < 0) continue;

            Coord muI(0.,0.,0.);
            if (iAniso >= 0) muI = mu(iAniso);

            Coord r = nn.relPositions.getAtom(ia);
            SxMatrix<Double> hIJ;
            //SX_START_TIMER(GaussIJ);
            if (iAniso >= 0 || jAniso >= 0)  {
               gaussIJ.setDelta (r, 2. * beta2, 2. * betap2, 4);
               hIJ = gaussIJ.compute (1, 1, cg, 2);
            } else {
               gaussIJ.setDelta (r, 2. * beta2, 2. * betap2, 2);
               hIJ = gaussIJ.compute (0, 0, cg, 2);
            }
            //SX_STOP_TIMER(GaussIJ);
            // note: need a prefactor of (2L + 1)!!/sqrt(2L+1) due to
            // normalization convention of the gaussIJ kernel for expressions
            // involving G^L

            // minus sign: r = r(i) - r(j) = -r_ij
            SxComplex16 phase = exp (-I * (kVec ^ r));

            if (jAniso >= 0 && iAniso >= 0)  {
               SxMatrix3<Complex16> d2EdMu2;
               for (int a = 1; a < 4; a++)
                  for (int b = 1; b < 4; b++)
                     d2EdMu2(a % 3, b % 3) = 3. * hIJ(a,b) * phase;
               d2EdMu2 /= FOUR_PI * dielecConstant;
               addH4 (iTl, jTl, &res, d2EdMu2, kVec); // Eq. (69)
            }

            if (jAniso >= 0)  {
               SxMatrix3<Double> d2EdRdMu;
               for (int b2 = 1; b2 < 4; ++b2)  {
                  // --- Q^k term of Eq. (68)
                  {
                     Coord f = gaussIJ.getForce (0, b2, cg);
                     for (int a = 0; a < 3; ++a)
                        d2EdRdMu(a,b2 % 3) = sqrt(3.) * qIso(iTl) * f(a);
                  }

                  // --- mu^k term of Eq. (68)
                  for (int a2 = 1; a2 < 4; ++a2) {
                     Coord f = gaussIJ.getForce (a2, b2, cg);
                     f *= 3. * muI(a2 % 3);
                     for (int a = 0; a < 3; ++a)
                        d2EdRdMu(a, b2 % 3) += f(a);
                  }
               }
               // pack indices for addH23
               PsiG d2EdRdMuPacked(6);
               SX_LOOP2(a,b)
                  d2EdRdMuPacked(idx21(a,b)) = d2EdRdMu(a,b);
               d2EdRdMuPacked /= FOUR_PI * dielecConstant;
               addH23 (iTl, jTl, &res, phase * d2EdRdMuPacked, kVec);
               addH23 (jTl, jTl, &res, -d2EdRdMuPacked, kVec);
            }

            // --- 2nd derivative from moving multipoles (Eq. 67)
            SxMatrix3<Double> Hij;

            // Q-Q term in Eq. (67) in implementation notes
            Hij = qIso(iTl) * qIso(jTl) * gaussIJ.getHesse (0, 0, cg);

            // Q-mu term in Eq. (67) in implementation notes
            for (int a2 = 1; a2 < 4; ++a2)  {
               double Qmu = qIso(iTl) * muJ(a2 % 3) - qIso(jTl) * muI(a2 % 3);
               if (fabs(Qmu) > 1e-12)
                  Hij -= sqrt(3.) * Qmu * gaussIJ.getHesse (0, a2, cg);
            }
            // mu-mu term in Eq. (67) in implementation notes
            if (muI.normSqr () > 1e-12 && muJ.normSqr () > 1e-12)  {
               for (int a2 = 1; a2 < 4; ++a2)  {
                  for (int b2 = 1; b2 < 4; ++b2)  {
                     double muMu = muI(a2 % 3) * muJ(b2 % 3);
                     Hij += 3. * muMu * gaussIJ.getHesse (a2, b2, cg);
                  }
               }
            }
            Hij /= FOUR_PI * dielecConstant;
            SX_LOOP2(a,b)  {
               res(3*iTl+a, 3*jTl+b) -= Hij(a,b) * phase;
               res(3*iTl+a, 3*iTl+b) += Hij(a,b);
            }
         }
      }
   }
// >>>>>> end of indent shift 
      } // jBlock
   } // iBlock
   if (tau.getNAtoms () > 100) (cout << endl).flush ();

   return res;
}

void SxHarmonicPotential::addH4 (int iTl, int jTl,
                                 SxMatrix<Complex16> *resPtr,
                                 const SxMatrix3<Complex16> &d2EdMu2,
                                 const Coord &kVec)
{
   // i=m, j=n in writeup
   SX_CHECK (resPtr);
   SxMatrix<Complex16> &res = *resPtr;
   int iAniso = mapAniso(iTl);
   int jAniso = mapAniso(jTl);
   SxMatrix3<Complex16> QQW 
      = SxMatrix3<Complex16>(QFrame(iAniso)).transpose ()
      ^ d2EdMu2 
      ^ SxMatrix3<Complex16>(QFrame(jAniso));

   ssize_t nFrameI = frameAtoms(iAniso).getSize ();
   ssize_t nFrameJ = frameAtoms(jAniso).getSize ();
   // delta(i,k) delta(j,l) term in w_ik w_jl
   SX_LOOP2(a,b)
      res(3*iTl+a, 3*jTl+b) += QQW(a,b); 
   SX_LOOP(iF) {
      int kTl = frameAtoms(iAniso)(iF);
      double kdI = kVec ^ frameAtomPos(iAniso)(iF);
      // put the complex conjugate in exponent
      SxComplex16 wkI = -exp(-I*kdI) / double(nFrameI);

      // delta(j,l) term 
      SX_LOOP2(a,b)
         res(3*kTl+a, 3*jTl+b) += wkI * QQW(a,b);
      // others
      SX_LOOP(jF) {
         int lTl = frameAtoms(jAniso)(jF);
         double kdJ = kVec ^ frameAtomPos(jAniso)(jF);
         SxComplex16 wkJ = -exp(I*kdJ) / double(nFrameJ);
         if (iF == 0)  {
            // delta(i,k) term 
            SX_LOOP2(a,b)
               res(3*iTl+a, 3*lTl+b) += wkJ * QQW(a,b);
         }
         // other terms of Eq. (66)
         SX_LOOP2(a,b)
            res(3*kTl+a, 3*lTl+b) += wkI * wkJ * QQW(a,b);
      }
   }
}

void SxHarmonicPotential::addH23 (int iTl, int jTl,
                                  SxMatrix<Complex16> *resPtr,
                                  const PsiG &d2EdRdMu,
                                  const Coord &kVec)
{
   SX_CHECK (resPtr);
   SxMatrix<Complex16> &res = *resPtr;
   int jAniso = mapAniso(jTl);
   SxMatrix3<Complex16> QW; QW.set (0.);
   SX_LOOP3(a,b,b2)
      QW(a,b) += QFrame(jAniso)(b2,b)
               * d2EdRdMu(idx21(a,b2));
   // delta(n,l) term in w(n,l)
   SX_LOOP2(a,b)  {
      res(3*iTl+a, 3*jTl+b) += QW(a,b); 
      res(3*jTl+b, 3*iTl+a) += QW(a,b).conj (); 
   }
   // frame atoms
   ssize_t nFrame = frameAtoms(jAniso).getSize ();
   for (int iF = 0; iF < nFrame; ++iF)  {
      int kTl = frameAtoms(jAniso)(iF);
      double kd = kVec ^ frameAtomPos(jAniso)(iF);
      SxComplex16 wk = -exp(I*kd) / double(nFrame);
      SX_LOOP2(a,b)  {
         // H-2
         res(3*iTl+a, 3*kTl+b) += wk * QW(a,b);
         // H-3
         res(3*kTl+b, 3*iTl+a) += (wk * QW(a,b)).conj ();
      }
   }
}

SxMatrix<Complex16> SxHarmonicPotential::getHesse (const Coord &kVec,
                                                   const SxAtomicStructure &tau)
{
   SX_CLOCK (HesseTotal);
   SxMatrix<Complex16> res = getHesseEl (kVec, tau);
   
   SX_CLOCK (HesseShort);
   // --- calculate local interactions
   for (int ia = 0; ia < tau.getNAtoms (); ++ia)  {
      for (int in = 0; in < neighborIdx(ia).getSize (); ++in)  {
         int idN = neighborIdx(ia)(in);

         // --- k-dependent phase
         int iGen = genNeighborId(ia)(in);
         int iSym = genSymId(ia)(in);
         Coord dist;
         if (iSym < 0)  {
            dist = genNeighbors(iGen);
         } else {
            dist = syms(iSym) ^ genNeighbors(iGen);
         }
         ///* --- needs to be tested: phase correction for displaced structures
         {
            // calculate how dist differs from actual coordinate difference
            // in tau
            Coord d = tau(idN) - tau(ia);
            d -= dist;
            tau.cell.map (&d, SxCell::Origin);
            dist += d;
         }
         SxComplex16 phase = exp(I * (dist ^ kVec));

         // --- update Hessian
         for (int b = 0; b < 3; ++b)  {
            for (int a = 0; a < 3; ++a)  {
               double Hab = hesse(ia)(in)(a,b);
               res(3 * ia  + a, 3 * ia  + b)  += Hab;
               res(3 * idN + a, 3 * ia  + b)  -= Hab * phase.conj ();
               res(3 * ia  + a, 3 * idN + b)  -= Hab * phase;
               res(3 * idN + a, 3 * idN + b)  += Hab;
            }
         }
      }
   }

   if (tau.getNAtoms () < 100) cout << "H=" << res << endl;
   return res;
}

