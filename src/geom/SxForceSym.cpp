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

#include <SxForceSym.h>
#include <SxGrid.h>

void SxForceSym::setup (const SxAtomicStructure &str)
{
   SX_CHECK(str.cell.symGroupPtr);
   SX_CHECK(str.atomInfo);
   const SxSymGroup &S = *str.cell.symGroupPtr;

   // get and check number of symmetries
   int nSym = S.getSize ();
   if (nSym <= 1) return;

   permutations.resize (nSym);
      
   // --- set up auxiliary grid for structure matching
   const int atomsPerCell = 3; // rough guess from very few test runs
   SxGrid grid (str, atomsPerCell);

   // compute permutations
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      permutations(iSym) = str.match (grid, S(iSym) ^ str);
   }

   // done (unbelievable simple, no?)
   
}

SxAtomicStructure SxForceSym::operator* (const SxAtomicStructure &forces) const
{
   // --- checks
   int nSym = getNSymmetries ();
   if (nSym == 0 && forces.cell.symGroupPtr->getSize () <= 1)
      return SxAtomicStructure(forces, SxAtomicStructure::Copy);
   SX_CHECK(forces.cell.symGroupPtr);
   const SxSymGroup &S = *forces.cell.symGroupPtr;
   
   SX_CHECK (nSym > 0); // must be initialized
   SX_CHECK (S.getSize () == nSym, S.getSize (), nSym);
   SX_CHECK (forces.atomInfo);
   SX_CHECK (permutations.getSize() == nSym,permutations.getSize(),nSym);


   SxAtomicStructure fSym(forces.cell, forces.getNAtoms (), forces.atomInfo);
   fSym.set (Coord(0.,0.,0.));

   // sum over symmetries
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      SX_CHECK (permutations(iSym).getPtr () != NULL);
      SX_CHECK (permutations(iSym)->isChild (forces.atomInfo));
      fSym += (S.getRot(iSym) ^ forces)  // rotate 
              .replaceInfo (permutations(iSym)); // and permute
      
      // remark: the actual permutation is done in operator+=
   }
   
   // average over symmetries
   fSym /= nSym; 

   return fSym;
}


bool SxForceSym::checkStr (const SxAtomicStructure &str) const
{
   // --- checks
   int nSym = getNSymmetries ();
   SX_CHECK (nSym > 0); // must be initialized
   SX_CHECK(str.cell.symGroupPtr);
   const SxSymGroup &syms = *str.cell.symGroupPtr;
   SX_CHECK (syms.getSize () == nSym, syms.getSize (), nSym);
   SX_CHECK (str.atomInfo);

   int nAtoms = str.getNAtoms ();
   Coord d;

   SxVector<Int>::Iterator mapIt;
   SxSymOp S;
   for (int iSym = 0; iSym < nSym; ++iSym)  {
      // more checks
      SX_CHECK (nAtoms == permutations(iSym)->parentMap.getSize (),
                nAtoms, permutations(iSym)->parentMap.getSize ());
      SX_CHECK (permutations(iSym)->isChild(str.atomInfo));

      mapIt = permutations(iSym)->parentMap.begin ();
      // check atom by atom with permutation
      S = syms(iSym);
      for (int ia = 0; ia < nAtoms; ++ia)  {
         d= str.constRef(*mapIt++) - (S ^ str.constRef(ia));
         str.cell.map (&d, SxCell::Origin);
         if (   fabs(d(0)) > str.epsEqual
             || fabs(d(1)) > str.epsEqual
             || fabs(d(2)) > str.epsEqual)
            return false;
      }
   }
   return true;
}
