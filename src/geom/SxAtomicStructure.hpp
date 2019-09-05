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

/// Rotate atoms around origin
inline SxAtomicStructure operator^ (const SxMatrix<TPrecTauR> &rot, 
                                    const SxAtomicStructure &str)
{
   SX_CHECK (!str.isCreationMode ());
   SX_CHECK (rot.nRows () == 3 && rot.nCols () == 3,
             rot.nRows (), rot.nCols ());
   SxAtomicStructure res (str, SxAtomicStructure::Reference);
   res.coords = rot ^ str.coords;
   // semi-bug in vector class: for 1 atom coords looks like a vector
   if (str.getNAtoms () == 1) 
      res.coords.reshape (3,1 /* nAtoms */);
   return res;
}

/// Symmetry-transform atoms
inline SxAtomicStructure operator^ (const SxSymOp &op, 
                                    const SxAtomicStructure &str)
{
   SX_CHECK (!str.isCreationMode ());
   SxAtomicStructure res (str.getNewStr ());
   for (int ia = 0; ia < str.getNAtoms (); ++ia)
      res.ref(ia) = op ^ str.constRef(ia);
   return res;
}

/// Rotate atoms around origin
inline SxAtomicStructure operator^ (const SymMat &S, 
                                    const SxAtomicStructure &str)
{
   SX_CHECK (!str.isCreationMode ());
   SxAtomicStructure res (str.getNewStr ());
   for (int ia = 0; ia < str.getNAtoms (); ++ia)
      res.ref(ia) = S ^ str.constRef(ia);
   return res;
}

SxAtomicStructure & 
SxAtomicStructure::operator^= (const SxMatrix<TPrecTauR> &rot)
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (rot.nRows () == 3 && rot.nCols () == 3,
             rot.nRows (), rot.nCols ());
   coords = rot ^ coords;
   // semi-bug in vector class: for 1 atom coords looks like a vector
   if (nTlAtoms == 1) 
      coords.reshape (3,1 /* nAtoms */);
   return *this;
}



inline SxAtomicStructure SxAtomicStructure::operator- () const
{
   SX_CHECK(!isCreationMode ());
   SX_CHECK (getNAtoms () > 0);
   SxAtomicStructure res (*this, SxAtomicStructure::Reference);
   res.set(-coords,SxAtomicStructure::Reference);
   return res;
}

inline SxAtomicStructure SxAtomicStructure::operator* (double s) const
{
   SX_CHECK(!isCreationMode ());
   SX_CHECK (getNAtoms () > 0);
   SxAtomicStructure res (*this, SxAtomicStructure::Reference);
   res.set(s * coords,SxAtomicStructure::Reference);
   return res;
}



SxVector3<TPrecTauR>& SxAtomicStructure::ref (ssize_t iTlAtom)
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (coords.getSize () == 3 * nTlAtoms,
             coords.getSize (), 3 * nTlAtoms);
   SX_CHECK (iTlAtom >= 0 && iTlAtom < nTlAtoms, iTlAtom, nTlAtoms);
   return Coord::toVec3Ref(coords.elements + 3 * iTlAtom);
}

const SxVector3<TPrecTauR>& SxAtomicStructure::constRef (ssize_t iTlAtom) const
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (coords.getSize () == 3 * nTlAtoms,
             coords.getSize (), 3 * nTlAtoms);
   SX_CHECK (iTlAtom >= 0 && iTlAtom < nTlAtoms, iTlAtom, nTlAtoms);
   return Coord::toVec3Ref(coords.elements + 3 * iTlAtom);
}

SxVector3<TPrecTauR>& SxAtomicStructure::ref (ssize_t iSpecies, ssize_t iAtom)
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK(atomInfo);
   SX_CHECK (iSpecies < getNSpecies () && iSpecies >= 0,
             iSpecies, getNSpecies ());
   SX_CHECK (iAtom < atomInfo->nAtoms(iSpecies) && iAtom >= 0,
             iAtom, atomInfo->nAtoms(iSpecies));
   return ref(iAtom + atomInfo->offset(iSpecies));

}

const SxVector3<TPrecTauR> &
SxAtomicStructure::constRef (ssize_t iSpecies, ssize_t iAtom) const
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK(atomInfo);
   SX_CHECK (iSpecies < getNSpecies () && iSpecies >= 0,
             iSpecies, getNSpecies ());
   SX_CHECK (iAtom < atomInfo->nAtoms(iSpecies) && iAtom >= 0,
             iAtom, atomInfo->nAtoms(iSpecies));
   return constRef(iAtom + atomInfo->offset(iSpecies));
}

void SxAtomicStructure::setAtom (int iSpecies, int iAtom, 
                                 const Coord &pos)
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (atomInfo);
   SX_CHECK (iSpecies < getNSpecies () && iSpecies >= 0,
             iSpecies, getNSpecies ());
   SX_CHECK (iAtom < atomInfo->nAtoms(iSpecies) && iAtom >= 0,
             iAtom, atomInfo->nAtoms(iSpecies));
   ref(iSpecies,iAtom) = pos;
}

// ---- comparison operator
bool SxAtomicStructure::operator== (const SxAtomicStructure &x) const
{
   return isEqual(x, NULL);
}

bool SxAtomicStructure::operator!= (const SxAtomicStructure &x) const
{
   return !(*this == x);
} 

int SxAtomicStructure::getIAtom (int iSpecies, int iAtom) const 
{
   SX_CHECK (atomInfo);
   SX_CHECK (iSpecies < atomInfo->nSpecies, iSpecies, atomInfo->nSpecies);
   SX_CHECK (iAtom < atomInfo->nAtoms(iSpecies),
             iAtom, atomInfo->nAtoms(iSpecies));
   return atomInfo->offset(iSpecies) + iAtom;
}

SxAtomicStructure SxAtomicStructure::getNewStr () const
{
   SxAtomicStructure res(nTlAtoms, atomInfo);
   res.cell = cell; // may not be initialized
   res.epsEqual = epsEqual;
   return res;
}

