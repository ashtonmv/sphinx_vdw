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

#include <SxEquivalence.h>
#include <SxNeighbors.h>
#include <SxSimpleParser.h>

void SxEquivalence::setup (const SxAtomicStructure &primStrIn)
{
   // deep copy of atoms
   primStr.copy (primStrIn);
   // deep copy of symmetries
   primStr.cell.symGroupPtr
      = SxPtr<SxSymGroup>::create (*primStrIn.cell.symGroupPtr);
   SX_CHECK (primStr.cell.symGroupPtr);
   const SxSymGroup &sym = *primStr.cell.symGroupPtr;
   SxGrid primGrid (primStr, 10);
   // check that symmetry group is primitive
   SX_CHECK (sym.getNPrimitive () == sym.getSize (),
             sym.getNPrimitive (), sym.getSize ());

   // --- Generate equivalence classes of primitive structure
   equivId.resize (primStr.getNAtoms ());
   mappingRot.resize (primStr.getNAtoms ());
   equivId.set (-1);
   int nEquiv = 0;
   SxList<SxList<int> > matSymList;
   SX_LOOP (ia)  {
      if (equivId(ia) >= 0) continue;
      matSymList.append (SxList<int> ());

      // --- loop over symmetries to find equivalent atoms
      for (int iSym = 0; iSym < sym.getSize (); ++iSym)  {
         int ja = primStr.find (sym(iSym) ^ primStr(ia), primGrid);
         if (ja < 0) { cout << "Symmetry problem" << endl; SX_EXIT; }
         if (equivId(ja) == -1)  {
            mappingRot(ja) = iSym;
            equivId(ja) = nEquiv;
         }
         // atom maps onto itself -> matrix symmetry
         if (ia == ja) matSymList.last () << iSym;
      }

      nEquiv++;
   }

   // --- setup matSym from temporary
   matSym.resize (nEquiv);
   for (int iEq = 0; iEq < nEquiv; ++iEq)  {
      matSym(iEq) = matSymList.first ();
      matSymList.removeFirst ();
   }
}

// --- print out an index list (+1 for each int)
static void fprintList(FILE *fp, const SxArray<int> &list)
{
   sxfprintf (fp, "[%d", list(0) + 1);
   for (int i = 1; i < list.getSize (); i++)
      sxfprintf (fp, ", %d", list(i) + 1);
   sxfprintf (fp, "];");
}

void SxEquivalence::fprintsx (FILE *fp) const
{
   sxfprintf (fp, "equivalence {\n");
   sxfprintf (fp, "\n   // --- equivalence reference structure\n");
   primStr.fprint (fp, SxAtomicStructure::PrintSymmetries);
   sxfprintf (fp, "\n   // --- equivalence classes\n");
   SX_LOOP (iEq)  {
      sxfprintf (fp, "   equivalenceClass { /* id = %d */ ", int(iEq + 1));
      if (matSym(iEq).getSize () > 1)  {
         sxfprintf (fp, "matSym = ");
         fprintList (fp, matSym(iEq));
      }
      sxfprintf (fp, " }\n");
      sxfprintf (fp, "   // => member atoms:");
      SX_LOOP(ia) if (equivId(ia) == iEq) sxfprintf(fp, " %d", int(ia+1));
      sxfprintf (fp, "\n");
   }
   sxfprintf (fp, "   // ---\n");

   // --- print out mappings to equivalence classes (and their representative)
   sxfprintf (fp, "\n   // equivalence class ids for all atoms\n");
   sxfprintf (fp, "   equivalenceIds = ");
   fprintList (fp, equivId);
   sxfprintf (fp, "\n");
   sxfprintf (fp, "   // symmetry id that maps from representative atom, for all atoms\n");
   sxfprintf (fp, "   mappingSym = ");
   fprintList (fp, mappingRot);
   sxfprintf (fp, "\n}\n\n");
}

void SxEquivalence::read (const SxSymbolTable *table)
{
   SYMBOLPARSE(table)  {
      SYMBOLGROUP("equivalence")  {
         primStr = SxAtomicStructure (SYMBOLGROUP_TABLE);
         int nEquiv = SYMBOL_COUNT("equivalenceClass");
         matSym.resize (nEquiv);
         int iEq = 0;
         FOREACH_SYMBOLGROUP("equivalenceClass")
            matSym(iEq++) << SYMBOLGET("matSym");
         equivId = SYMBOLGET("equivalenceIds");
         mappingRot = SYMBOLGET ("mappingSym");
      }
   }
   SX_LOOP (ia) equivId(ia) -= 1;
   SX_LOOP (ia) mappingRot(ia) -= 1;
   SX_LOOP (iEq) {
      if (matSym(iEq).getSize () > 0)
         SX_LOOP(iSym) matSym(iEq)(iSym) -= 1;
   }
}

SxEquivalence
SxEquivalence::mapToComplex (const SxAtomicStructure &supercell, double dist)
const
{
   // make sure that equivId refers to primStr
   // if not, you must use reset to get back to the original primStr
   SX_CHECK (equivId.getSize () == primStr.getNAtoms (),
             equivId.getSize (), primStr.getNAtoms ());
   SX_CHECK (equivId.getSize () == mappingRot.getSize ());
   SxEquivalence res;
   res.primStr = primStr;
   res.matSym  = matSym;
   res.equivId   .resize (supercell.getNAtoms ());
   res.mappingRot.resize (supercell.getNAtoms ());

   SxGrid primGrid (primStr, 10);
   SxNeighbors nn;
   int mode = SxNeighbors::StoreIdx
            | SxNeighbors::StoreAbs
            | SxNeighbors::IncludeZeroDistance;
   SX_LOOP(ia)  {
      int iaPrim = primStr.find (supercell(ia), primGrid);
      nn.compute (primGrid, primStr, supercell(ia), dist, mode);
      int is = supercell.getISpecies ((int)ia);
      if (nn.getSize () == 1) {
         iaPrim = nn.idx(0);
      } else if (nn.absPositions.getNAtoms (is) == 1)  {
         iaPrim = nn.idx (nn.absPositions.getIAtom(is, 0));
      } else {
         cout << "Warn: failed to map equivalence for atom " << (ia + 1)
              << " @ " << supercell(ia)
              << " to primitive cell" << endl;
         if (nn.getSize () > 1)  {
            cout << "Candidates are " << endl;
            SX_LOOP(in)
               cout << "primitive atom @ " << primStr(nn.idx(in)) << endl;
            SX_QUIT;
         }
      }

      if (iaPrim < 0) {
         // atom does not map
         res.equivId(ia) = -1;
         res.mappingRot(ia) = -1;
      } else {
         res.equivId(ia) = equivId(iaPrim);
         res.mappingRot(ia) = mappingRot(iaPrim);
      }
   }
   return res;
}

void SxEquivalence::mapSyms (const SxArray<SymMat> &newSym)
{
   SxArray<int> map(primStr.cell.symGroupPtr->getSize ());
   SX_LOOP (iSym)  {
      const SymMat &S = primStr.cell.symGroupPtr->getRot ((int)iSym);
      map(iSym) = -1;
      SX_LOOP (jSym)  {
         if ((newSym(jSym) - S).absSqr ().sum () < primStr.cell.epsSym)  {
            map(iSym) = (int)jSym;
            break;
         }
      }
   }
   SX_LOOP(ia)  {
      mappingRot(ia) = map(mappingRot(ia));
   }
   SX_LOOP(iEq)  {
      if (matSym(iEq).getSize () > 0)
         SX_LOOP(iSym) matSym(iEq)(iSym) = map (matSym(iEq)(iSym) );
   }
}
