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

#include <SxAtomicStructure.h>
#include <SxSortedList.h>
#include <SxRotation.h>

#include <SxGrid.h>
#include <SxNeighbors.h>
#include <SxStack.h>

SxAtomicStructure::SxAtomicStructure ()
   : nTlAtoms(0),
     creationStatus(Unknown),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   /* empty */
}

SxAtomicStructure::SxAtomicStructure (const SxCell& inCell)
   : nTlAtoms(0),
     cell(inCell),
     creationStatus(Unknown),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   /* empty */
}

SxAtomicStructure::SxAtomicStructure (const SxAtomicStructure& in,
                                      enum DataHandling handling)
   : nTlAtoms(in.nTlAtoms),
     atomInfo(in.atomInfo),
     creationStatus(in.creationStatus),
     epsEqual(in.epsEqual)
{
   SX_CHECK (in.creationStatus != ListConstruction);
   // don't initialize since cell need not be set up...
   cell = in.cell;
   if (nTlAtoms > 0)  {
      if (handling == Reference)  {
         coords   = in.coords;
      }  else  {
         coords.copy (in.coords);
      }
   }
}

SxAtomicStructure::SxAtomicStructure (const SxSymbolTable* table,
                                      const SxString &coordName)
   : nTlAtoms(0),
     creationStatus(ListConstruction),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   SX_CHECK (table);
   const SxSymbolTable *str = NULL, *species, *atom;
   const SxString strName = "structure";
   const SxString specName = "species";
   const SxString atomName = "atom";
   const SxString relName = "relative";
   Coord coord;
   SxList<SxString> labels;
   SxList<SxString> elements;
   bool haveElements = false;
   try  {
      // get structure group
      str = (table->name == strName) ? table : table->getGroup (strName);

      // --- get cell
      if (str->contains("cell")) cell = SxCell(str);

      // --- get epsEqual
      if (str->contains("epsEqual"))  {
         epsEqual = str->get ("epsEqual")->toReal ();
      }

      // --- get atomic position
      for (species = str->getGroup (specName);
           species != NULL;
           species = species->nextSibling (specName))
      {
         newSpecies ();
         if (species->contains ("element"))
            elements << species->get("element")->toString ();
         else
            elements << SxString ();
         if (elements.last ().getSize () > 0) haveElements = true;
         SxList<SxSymbolTable *>::ConstIterator itAtom;
         for (itAtom  = species->children.begin ();
              itAtom != species->children.end ();
            ++itAtom)
         {
            if ( (*itAtom)->name != atomName ) continue;

            atom = *itAtom;
            // get coords
            if (atom->contains(coordName, true))
               coord = Coord(atom->get(coordName)->toList ());
            else
               coord.set (0.);

            // relative coords ?
            bool hasRelName = atom->contains (relName, true);
            if (hasRelName && isPeriodic ())
               cell.changeToCar (&coord);
            else if (hasRelName && !isPeriodic ())  {
               cout << " relative coordinates can't be used without a cell !\n";
               SX_QUIT;
            }

            // append atom
            addAtom (coord);

            if (atom->contains ("label", true))  {
               if (labels.getSize () == 0) labels.resize (nTlAtoms - 1);
               SxString label = atom->get("label")->toString ();
               if (label.contains ('|'))  {
                  cout << "Illegal label: '" << label << "'." << endl;
                  cout << "Atom labels must not contain bars |." << endl;
                  SX_QUIT; 
               }
               labels << label;
            } else if (labels.getSize () > 0)  {
               labels << SxString ();
            }
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   endCreation ();

   // --- set labels
   if (labels.getSize () > 0)  {
      SX_CHECK (labels.getSize () == nTlAtoms,
                labels.getSize (), nTlAtoms);
      atomInfo->meta.attach (Labels, SxArray<SxString>(labels));
   }
   if (haveElements)
      atomInfo->meta.attach (Elements, SxArray<SxString> (elements));

   cout << SX_SEPARATOR;
   if (!str->containsGroup("symmetry"))  {
      if (isPeriodic () && coordName == "coords")  {
         // --- perform automatic symmetry search
         cout << "| Automatic search for symmetry elements\n";
         cell.symGroupPtr = SxPtr<SxSymGroup>::create(getSymmetries ());
      }  else  {
         // use the identity
         cell.symGroupPtr = SxPtr<SxSymGroup>::create(SxList<SymMat> ()
                            << SymMat (1.,0.,0.,0.,1.,0.,0.,0.,1.));
      }
   }

   cout << "| Symmetry elements:\n";
   for (int iSym = 0; iSym < cell.symGroupPtr->getNSymmorphic (); ++iSym)  {
      cout << "| S=" << cell.symGroupPtr->getSymmorphic (iSym) << endl;
   }
   cout << SX_SEPARATOR;
}

SxAtomicStructure::SxAtomicStructure (const SxCell &inCell,
                                      const SxList<SxList<Coord> > &atoms)
   : cell(inCell),
     creationStatus(Unknown),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   set (atoms);
}

SxAtomicStructure::SxAtomicStructure (const SxBinIO &io)
   : epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   read (io);
}

SxAtomicStructure::SxAtomicStructure (const SxList<Coord> &coord)
   : nTlAtoms(int(coord.getSize ())),
     creationStatus(Finished),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   if (nTlAtoms == 0) return;

   // allocate memory
   coords.reformat(3, nTlAtoms);
   // copy coordinates
   PrecTauR *pos = coords.elements;
   SxList<Coord>::ConstIterator coordIt = coord.begin ();
   for (int ia = 0; ia < nTlAtoms; ++ia, pos+=3)
      Coord::toVec3Ref(pos) = *coordIt++;
}

SxAtomicStructure::SxAtomicStructure (const SxArray<Coord> &coord)
   : nTlAtoms(int(coord.getSize ())),
     creationStatus(Finished),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   if (nTlAtoms == 0) return;

   // allocate memory
   coords.reformat(3, nTlAtoms);
   // copy coordinates
   PrecTauR *pos = coords.elements;
   for (int ia = 0; ia < nTlAtoms; ++ia, pos+=3)
      Coord::toVec3Ref(pos) = coord(ia);
}

SxAtomicStructure::SxAtomicStructure (const SxCell &cellIn,
                                      int nAtoms,
                                      const SxAtomInfo::ConstPtr &atomInfoIn)
   : nTlAtoms(nAtoms),
     cell(cellIn),
     atomInfo(atomInfoIn),
     creationStatus(Unknown),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   SX_CHECK (!atomInfo || atomInfo->nAtoms.sum () == nAtoms);
   if (nAtoms > 0)  {
      coords.reformat (3, nAtoms);
      creationStatus = Finished;
   }
}

SxAtomicStructure::SxAtomicStructure (int nAtoms,
                                      const SxAtomInfo::ConstPtr &atomInfoIn)
   : nTlAtoms(nAtoms),
     atomInfo(atomInfoIn),
     creationStatus(Unknown),
     epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   SX_CHECK (!atomInfo || atomInfo->nAtoms.sum () == nAtoms);
   if (nAtoms > 0)  {
      coords.reformat (3, nAtoms);
      creationStatus = Finished;
   }
}

SxAtomicStructure::SxAtomicStructure (const SxVector<TPrecTauR> &coord,
                                      enum DataHandling handling)
 : nTlAtoms((int)coord.getSize () / 3),
   creationStatus(Finished),
   epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   SX_CHECK(handling == Reference || handling == Copy);
   if (handling == Reference)
   {
      coords = coord;
      // no reshape because vector class doesn't allow to change
      // shape only for referencing vector (we don't want to change the
      // shape of coord, since it is const)
      // if (nTlAtoms > 0) coords.reshape (3, nTlAtoms);
   } else {
      coords.copy(coord);
      if (nTlAtoms > 0) coords.reshape (3, nTlAtoms);
   }
   SX_CHECK (coord.getSize () % 3 == 0);
   if (nTlAtoms == 0) creationStatus = Unknown;
}

SxAtomicStructure::SxAtomicStructure (const SxAtomInfo::ConstPtr &atomInfoIn,
                                      const SxAtomicStructure &from)
 : atomInfo(atomInfoIn),
   creationStatus(Finished),
   epsEqual(SX_EPS_STRUCT_DEFAULT)
{
   cell = from.cell; // cell need not be set up
   epsEqual = from.epsEqual;
   SX_CHECK (atomInfoIn);
   SX_CHECK (from.atomInfo);
   SX_CHECK (atomInfoIn->getParent () == from.atomInfo);
   nTlAtoms = (int)atomInfo->parentMap.getSize ();
   SX_CHECK (atomInfo->nAtoms.sum () == nTlAtoms,
             atomInfo->nAtoms.sum (), nTlAtoms);
   if (nTlAtoms == 0) return;

   coords.reformat (3, nTlAtoms);
   SxVector<Int>::Iterator it = atomInfo->parentMap.begin ();
   for (int ia = 0; ia < nTlAtoms; ++ia,++it)
      ref(ia) = from.constRef(*it);
}

SxAtomicStructure &SxAtomicStructure::operator= (const SxAtomicStructure &in)
{
   SX_CHECK (!in.isCreationMode ());
   SX_CHECK (!isCreationMode ());
   if ( &in == this )  return *this;
   cell     = in.cell;
   coords   = in.coords;
   nTlAtoms = in.nTlAtoms;
   epsEqual = in.epsEqual;
   atomInfo = in.atomInfo;
   creationStatus = nTlAtoms == 0 ? Unknown : Finished;
   return *this;
}

void SxAtomicStructure::operator<<= (const SxAtomicStructure &in)
{
   SX_CHECK (!in.isCreationMode ());
   SX_CHECK (!isCreationMode ());
   if ( &in == this || nTlAtoms == 0)  return;
   if (   (in.atomInfo == atomInfo) // standard case: same atomInfo
       || (!atomInfo || !in.atomInfo))  // or atomInfo is lacking
   {
      // --- in and this are compatible
      SX_CHECK (nTlAtoms == in.nTlAtoms, nTlAtoms, in.nTlAtoms);
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords <<= in.coords;
   } else if (in.atomInfo->isChild (atomInfo))  {
      // --- in is derived from this
      SxVector<Int> idx = in.atomInfo->getIdxMap (atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      int na = in.nTlAtoms;
      for (int ia = 0; ia < na; ++ia, ++it)
         ref(*it) = in.constRef (ia);
   } else if (atomInfo->isChild (in.atomInfo))  {
      // --- this is derived from in
      SxVector<Int> idx = atomInfo->getIdxMap (in.atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia, ++it)
         ref(ia) = in.constRef (*it);
   } else if ( *in.atomInfo == *atomInfo)  {
      // --- in and this are size compatible, but not related
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords <<= in.coords;
   } else {
      // mismatch between this and in
      SX_EXIT;
   }

}

SxAtomicStructure
SxAtomicStructure::operator+ (const SxAtomicStructure &in) const
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (in.creationStatus == Finished);

   if (   (in.atomInfo == atomInfo) // standard case: same atomInfo
       || (!atomInfo || !in.atomInfo))  // or atomInfo is lacking
   {
      // --- in and this are compatible
      SX_CHECK (nTlAtoms == in.nTlAtoms, nTlAtoms, in.nTlAtoms);
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      SxAtomicStructure res ((atomInfo ? *this : in), Reference);
      res.coords = coords + in.coords;
      return res;
   } else if (in.atomInfo->isChild (atomInfo))  {
      // --- in is derived from this
      SX_EXIT;
      // not implemented, which atomInfo should be used?
   } else if (atomInfo->isChild (in.atomInfo))  {
      // --- this is derived from in
      SX_EXIT;
      // not implemented, which atomInfo should be used?
   } else if ( *in.atomInfo == *atomInfo)  {
      // size compatible, but not related
      SX_EXIT;
      // not implemented, I don't see any use for this
   } else {
      // mismatch between this and in
      SX_EXIT;
   }
   SX_EXIT;
   return SxAtomicStructure ();

}

SxAtomicStructure
SxAtomicStructure::operator+ (const Coord &trans) const
{
   SX_CHECK (creationStatus == Finished);

   SxAtomicStructure res(getNewStr ());
   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      res.ref(iAtom) = constRef(iAtom) + trans;
   return res;
}

SxAtomicStructure
SxAtomicStructure::operator- (const SxAtomicStructure &in) const
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (in.creationStatus == Finished);

   if (   (in.atomInfo == atomInfo) // standard case: same atomInfo
       || (!atomInfo || !in.atomInfo))  // or atomInfo is lacking
   {
      // --- in and this are compatible
      SX_CHECK (nTlAtoms == in.nTlAtoms, nTlAtoms, in.nTlAtoms);
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      SxAtomicStructure res ((atomInfo ? *this : in), Reference);
      res.coords = coords - in.coords;
      return res;
   } else if (in.atomInfo->isChild (atomInfo))  {
      // --- in is derived from this
      SX_EXIT;
      // not implemented, which atomInfo should be used?
   } else if (atomInfo->isChild (in.atomInfo))  {
      // --- this is derived from in
      SX_EXIT;
      // not implemented, which atomInfo should be used?
   } else if ( *in.atomInfo == *atomInfo)  {
      // size compatible, but not related
      SX_EXIT;
      // not implemented, I don't see any use for this
   } else {
      // mismatch between this and in
      SX_EXIT;
   }

   SX_EXIT;
   return SxAtomicStructure ();
}

SxAtomicStructure
SxAtomicStructure::operator- (const Coord &trans) const
{
   SX_CHECK (creationStatus == Finished);

   SxAtomicStructure res (getNewStr ());
   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      res.ref(iAtom) = constRef(iAtom) - trans;
   return res;
}

SxAtomicStructure &
SxAtomicStructure::operator+= (const SxAtomicStructure &in)
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (in.creationStatus == Finished);


   if (   (in.atomInfo == atomInfo) // standard case: same atomInfo
       || (!atomInfo || !in.atomInfo))  // or atomInfo is lacking
   {
      // --- in and this are compatible
      SX_CHECK (nTlAtoms == in.nTlAtoms, nTlAtoms, in.nTlAtoms);
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords += in.coords;
   } else if (in.atomInfo->isChild (atomInfo))  {
      // --- in is derived from this
      SxVector<Int> idx = in.atomInfo->getIdxMap (atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      int na = in.nTlAtoms;
      for (int ia = 0; ia < na; ++ia, ++it)
         ref(*it) += in.constRef (ia);
   } else if (atomInfo->isChild (in.atomInfo))  {
      // --- this is derived from in
      SxVector<Int> idx = atomInfo->getIdxMap (in.atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia, ++it)
         ref(ia) += in.constRef (*it);
   } else if ( *in.atomInfo == *atomInfo)  {
      // --- in and this are size compatible, but not related
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords += in.coords;
   } else {
      // mismatch between this and in
      SX_EXIT;
   }

   return *this;
}

SxAtomicStructure &
SxAtomicStructure::operator+= (const Coord &trans)
{
   SX_CHECK (creationStatus == Finished);

   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      ref(iAtom) += trans;
   return *this;
}

SxAtomicStructure &
SxAtomicStructure::operator-= (const SxAtomicStructure &in)
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (in.creationStatus == Finished);


   if (   (in.atomInfo == atomInfo) // standard case: same atomInfo
       || (!atomInfo || !in.atomInfo))  // or atomInfo is lacking
   {
      // --- in and this are compatible
      SX_CHECK (nTlAtoms == in.nTlAtoms, nTlAtoms, in.nTlAtoms);
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords -= in.coords;
   } else if (in.atomInfo->isChild (atomInfo))  {
      // --- in is derived from this
      SxVector<Int> idx = in.atomInfo->getIdxMap (atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      int na = in.nTlAtoms;
      for (int ia = 0; ia < na; ++ia, ++it)
         ref(*it) -= in.constRef (ia);
   } else if (atomInfo->isChild (in.atomInfo))  {
      // --- this is derived from in
      SxVector<Int> idx = atomInfo->getIdxMap (in.atomInfo);
      SxVector<Int>::Iterator it = idx.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia, ++it)
         ref(ia) -= in.constRef (*it);
   } else if ( *in.atomInfo == *atomInfo)  {
      // --- in and this are size compatible, but not related
      SX_CHECK (in.coords.getSize() == coords.getSize(),
                in.coords.getSize(),   coords.getSize());
      coords -= in.coords;
   } else {
      // mismatch between this and in
      SX_EXIT;
   }

   return *this;
}

SxAtomicStructure &
SxAtomicStructure::operator-= (const Coord &trans)
{
   SX_CHECK (creationStatus == Finished);

   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      ref(iAtom) -= trans;
   return *this;
}

// ---- mapping

SxAtomicStructure
SxAtomicStructure::operator% (const SxCell &mapCell) const
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (mapCell.volume > 0.);

   SxAtomicStructure res (getNewStr ());
   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      res.ref(iAtom) = mapCell.getMapped(constRef(iAtom));
   return res;
}

SxAtomicStructure & SxAtomicStructure::map(const SxCell &mapCell,
                                           enum SxCell::Mapping mode)
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (mapCell.volume > 0.);

   for (int iAtom = 0; iAtom < nTlAtoms; ++iAtom)
      mapCell.map(&ref(iAtom), mode);
   return *this;
}

// --------------- setting
void SxAtomicStructure::newSpecies ()
{
   if (creationStatus != ListConstruction && nTlAtoms == 0) startCreation ();
   SX_CHECK(creationStatus == ListConstruction);
   posList.append(SxList<Coord> ());
}

SxAtomicStructure& SxAtomicStructure::set (const Coord &fillVal)
{
   SX_CHECK (creationStatus == Finished);
   for (int i=0; i < nTlAtoms; i++)  coords.colRef(i) <<= fillVal;
   return *this;
}

void SxAtomicStructure::importStack(const SxStack<Coord> &stack, int is)
{
   // This is a hack to get access to SxStack's protected member
   // function exportStack
   // I don't think it is worth specializing SxStack<Coord> just to
   // become its friend
   class StackAccess : public SxStack<Coord> {
      public:
         void exportStack (Coord *target, size_t n) const {
            SxStack<Coord>::exportStack (target, n);
         }
   };
   const StackAccess &stackAlias = static_cast<const StackAccess&>(stack);
   if (is == -1)  {
      SX_CHECK(nTlAtoms == (int)stack.getSize (),
               nTlAtoms, stack.getSize ());
      stackAlias.exportStack((Coord*)coords.elements, nTlAtoms);
   } else {
      SX_CHECK (atomInfo);
      SX_CHECK (getNAtoms(is) == (int)stack.getSize (),
                getNAtoms(is), stack.getSize ());
      stackAlias.exportStack((Coord*)coords.elements + atomInfo->offset(is),
                              getNAtoms(is));
   }
}

// ------ accessing atoms

SxAtomicStructure SxAtomicStructure::repeat (int na1, int na2, int na3, bool symmetrize_) const
{
   SX_CHECK (na1 > 0, na1);
   SX_CHECK (na2 > 0, na1);
   SX_CHECK (na3 > 0, na1);

   return repeat (SxVector3<Int> (na1, na2, na3), symmetrize_);
}

SxAtomicStructure SxAtomicStructure::repeat (const SxVector3<Int> &repVec, bool symmetrize_) const
{
   SX_CHECK (repVec(0) > 0, repVec(0));
   SX_CHECK (repVec(1) > 0, repVec(1));
   SX_CHECK (repVec(2) > 0, repVec(2));
   SX_CHECK (cell.volume > 0.);

   SxAtomicStructure repStruct;
   int is, ia;
   SxVector3<Int> a;

   repStruct.cell = CellMat (repVec(0)*cell.col(0),
                             repVec(1)*cell.col(1),
                             repVec(2)*cell.col(2)).transpose();
   repStruct.cell.epsSym = cell.epsSym;

   for (is = 0; is < getNSpecies (); is++)  {
      repStruct.newSpecies ();
      for (ia = 0; ia < getNAtoms(is); ia++)  {
         for (a(0) = 0; a(0) < repVec(0); a(0)++)  {
            for (a(1) = 0; a(1) < repVec(1); a(1)++)  {
               for (a(2) = 0; a(2) < repVec(2); a(2)++)  {
                  repStruct.addAtom( (*this)(is,ia) + (cell^a) );
               }
            }
         }
      }
   }
   repStruct.endCreation ();
   if (symmetrize_)
      repStruct.cell.symGroupPtr = SxPtr<SxSymGroup>::create (repStruct.getSymmetries ());
   if (atomInfo->meta.contains (Elements))
      repStruct.atomInfo->meta.attach (Elements, atomInfo->meta
                                       .get<SxArray<SxString> > (Elements));
   return repStruct;
}

SxAtomicStructure
SxAtomicStructure::repeat (const SxMatrix3<Int> &repMatrix, bool symmetrize_) const
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (cell.volume > 0.);

   // --- repeat atoms
   SxAtomicStructure res;
   res.cell = cell ^ repMatrix;

   SxCell genCell = cell;
   SxVector3<Int> mesh = genCell.computeGen (res.cell), n;
   SX_CHECK (abs(repMatrix.determinant ()) == mesh.product (),
             repMatrix.determinant (), mesh.product ());

   for (n(0) = 0; n(0) < mesh(0); ++n(0))
      for (n(1) = 0; n(1) < mesh(1); ++n(1))
         for (n(2) = 0; n(2) < mesh(2); ++n(2))
            res.append ((*this) + genCell.relToCar (n));
   res.endCreation ();
   if (symmetrize_)
      res.cell.symGroupPtr = SxPtr<SxSymGroup>::create (res.getSymmetries ());
   if (atomInfo->meta.contains (Elements))
      res.atomInfo->meta.attach (Elements, atomInfo->meta
                                           .get<SxArray<SxString> > (Elements));
   return res;
}

void SxAtomicStructure::copy (const SxAtomicStructure &in)
{
   SX_CHECK (!in.isCreationMode ());
   SX_CHECK (this != &in);
   cell   = in.cell;
   nTlAtoms = in.nTlAtoms;
   atomInfo = in.atomInfo;
   SX_CHECK(!atomInfo || (atomInfo->nAtoms.sum () == nTlAtoms),
            atomInfo->nAtoms.sum (), nTlAtoms);
   // make sure any previous data is released
   coords = SxMatrix<TPrecTauR> ();
   // do the coordinate copy
   coords.copy (in.coords);
   creationStatus = Finished;
}

void SxAtomicStructure::set (const SxList<SxList<Coord> > &atoms)
{
   SX_CHECK (!isCreationMode ());
   if (atoms.getSize () == 0)  {
      creationStatus = Unknown;
      nTlAtoms = 0;
      coords = SxMatrix<TPrecTauR> ();
      atomInfo = SxConstPtr<SxAtomInfo> ();
      return;
   }

   SxPtr<SxAtomInfo> aInfo = SxAtomInfo::create ();
   aInfo->nSpecies = int(atoms.getSize ());
   aInfo->nAtoms.resize (aInfo->nSpecies);

   nSpecies = aInfo->nSpecies;
   nTlAtoms = 0;
   SxList<SxList<Coord> >::ConstIterator specIt = atoms.begin ();
   for (int is = 0; is < aInfo->nSpecies; ++is, ++specIt)
      nTlAtoms += aInfo->nAtoms(is) = int((*specIt).getSize ());
   aInfo->setupOffset ();
   if (atomInfo.getPtr() != aInfo.getPtr()) atomInfo = aInfo;

   // allocate memory
   coords.reformat(3, nTlAtoms);

   // copy coordinates
   PrecTauR *pos = coords.elements;
   SxList<Coord>::ConstIterator coordIt;
   specIt = atoms.begin ();
   for (int is = 0; is < aInfo->nSpecies; ++is, ++specIt)  {
      coordIt = (*specIt).begin ();
      for (int ia = 0; ia < aInfo->nAtoms(is); ++ia, pos+=3)  {
         Coord::toVec3Ref(pos) = *coordIt++;
      }
   }
   creationStatus = Finished;

}

void SxAtomicStructure::set (const SxVector<Double> &pos,
                             enum DataHandling handling)
{
   SX_CHECK (&pos != &coords);
   SX_CHECK(!isCreationMode ());
   if (nTlAtoms == 0 && pos.getSize () == 0) return;
   SX_CHECK (handling == Reference || handling == Copy);
   SX_CHECK(pos.nRows () == 3 || handling == Copy, pos.nRows ());
   nTlAtoms = (int)pos.nCols ();

   if (handling == Reference)  {
      coords = pos;
   } else {// (handling = Copy)
      SX_CHECK (pos.getSize () % 3 == 0, pos.getSize ());
      nTlAtoms = (int)pos.getSize () / 3;
      coords = SxMatrix<Double> (3, nTlAtoms);
      coords <<= pos;
   }
   if (atomInfo && (atomInfo->nAtoms.sum () != nTlAtoms))
      atomInfo = SxConstPtr<SxAtomInfo> ();
   creationStatus = Finished;
}

void SxAtomicStructure::startCreation ()
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (posList.getSize () == 0, posList.getSize ());
   if (nTlAtoms > 0)  {
      if (atomInfo)  {
         for (int is = 0; is < atomInfo->nSpecies; ++is)  {
            posList.append(SxList<Coord> ());
            for (int ia = 0; ia < atomInfo->nAtoms(is); ++ia)  {
               posList.last ().append (getAtom(is,ia));
            }
         }
      } else {
         posList.append (SxList<Coord> ());
         SxList<Coord> &list = posList.last ();
         for (int i = 0; i < nTlAtoms; i++)
            list.append (getAtom(i));
      }
      coords = SxMatrix<TPrecTauR> ();
   }

   creationStatus = ListConstruction;
}

void SxAtomicStructure::endCreation ()
{
   if (creationStatus == Unknown)  {
      SX_CHECK (posList.getSize () == 0);
      return;
   }
   SX_CHECK (creationStatus == ListConstruction);
   creationStatus = Finished;
   set (posList);
   posList.resize (0);
}

const SxArray<SxString>& SxAtomicStructure::getLabels () const
{
   SX_CHECK (atomInfo);
   if (!atomInfo->meta.contains (Labels)) {
      cout << "No labels in atomic structure!" << endl;
      SX_QUIT;
   }
   return atomInfo->meta.get (Labels);
}

void SxAtomicStructure::resize (int nAtoms)
{
   SX_CHECK (!isCreationMode ());
   nTlAtoms = nAtoms;
   if (nAtoms > 0)
      coords.reformat (3, nAtoms);
   else
      coords = SxMatrix<TPrecTauR> ();
   // forget atom info if it doesn't match nAtoms
   if (atomInfo && atomInfo->nAtoms.sum () != nAtoms)
      atomInfo = SxConstPtr<SxAtomInfo> ();
   creationStatus = Finished;
}

void SxAtomicStructure::addAtom (const Coord &pos)
{
   if (creationStatus != ListConstruction && nTlAtoms == 0)  {
      startCreation ();
   }
   SX_CHECK(creationStatus == ListConstruction);
   if (posList.getSize () == 0) newSpecies ();
   posList.last ().append (pos);
   nTlAtoms++;
}

void SxAtomicStructure::addAtom (int iSpecies, const Coord &pos)
{
   if (creationStatus != ListConstruction && nTlAtoms == 0)  {
      startCreation ();
      newSpecies ();
   }
   SX_CHECK(creationStatus == ListConstruction);
   SX_CHECK (iSpecies >= 0 && iSpecies < posList.getSize (),
             iSpecies, posList.getSize ());
   posList(iSpecies).append (pos);

   nTlAtoms++;
}

void SxAtomicStructure::removeAtom (int iSpecies, int iAtom)
{
   SX_CHECK (creationStatus == ListConstruction);
   SX_CHECK (iSpecies >= 0 && iSpecies < posList.getSize (),
             iSpecies, posList.getSize ());
   SX_CHECK (iAtom < posList(iSpecies).getSize () && iAtom >= 0,
             posList(iSpecies).getSize (), iAtom);

   posList(iSpecies).remove (iAtom);
   nTlAtoms--; // one atom less

}

void SxAtomicStructure::append (const SxAtomicStructure &in)
{
   if (creationStatus != ListConstruction && nTlAtoms == 0) startCreation ();
   SX_CHECK (creationStatus == ListConstruction);
   int is,ia;
   for (is = 0; is < in.getNSpecies (); ++is)  {
      if (is == posList.getSize ()) newSpecies ();
      for (ia = 0; ia < in.getNAtoms(is); ++ia)
         addAtom (is, in(is,ia));
   }
}

// -------------- comparison
bool SxAtomicStructure::isEqual (const SxAtomicStructure &x,
                                 SxVector<Int>* mappingPtr,
                                 const bool debug) const
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (atomInfo); // todo: remove this constraint
   SX_CHECK (!atomInfo || getNSpecies () > 0);
   //SX_CHECK (atomInfo == x.atomInfo); // is this constraint necessary?
   // Same number of atoms
   if (x.getNAtoms () != nTlAtoms) { cout << "nAtom diff." << endl; return false; }
   double eps = epsEqual;
   if (epsEqual < x.epsEqual)  {
      cout << "Warning from str.isEqual: epsEquals differ!" << endl;
      eps = x.epsEqual;
   }//CHECK to be safe

   // --- a bunch of counters and temporary stuff that we need
   SxArray<bool> found (nTlAtoms); // atoms found in x
   found.set (false);
   int iSpecies, i, iAtom, jAtom, jStart, dim;
   int nSpecies = getNSpecies ();
   Coord position, dist;
   // do not map distance into cell if cell is uninitialized
   bool periodicCoords = cell.volume > 0.;

   // check that *mappingPtr has proper size
   SX_CHECK (!mappingPtr || mappingPtr->getSize () == nTlAtoms,
             mappingPtr->getSize (), nTlAtoms);

   // --- loop over species
   iAtom = 0;
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      jStart = iAtom;

      // --- loop over own atoms in species
      for (i = 0; i < getNAtoms(iSpecies); i++, iAtom++)  {
         position = getAtom (iAtom);

         // --- loop over atoms of same species in x
         for (jAtom = jStart; jAtom < jStart + getNAtoms(iSpecies); jAtom++)  {

            // atoms in x can be found only once (speeds up the thing)
            if (found(jAtom)) continue;

            // get distance vector
            dist = position - x.getAtom(jAtom);

            // map distance vector to be between -1/2 and +1/2 relative coords
            if (periodicCoords) cell.map (&dist, SxCell::Origin);

            // --- check that coordinates are zero
            for (dim = 0; dim < 3; dim++)
               if (fabs(dist(dim)) >= eps) break;

            if (dim == 3)  {
               // all coordinates differ by less than epsEqual
               found(jAtom) = true;
               if (mappingPtr) (*mappingPtr)(iAtom) = jAtom;
               break; // break j loop
            }
         } // for j
         if (jAtom == jStart + getNAtoms(iSpecies))  {
            // ith atom is not found
            if (debug) cout << iAtom << "not found! @" << position << endl;
            return false;
         }
      } // for i
   } // for iSpecies
   return true;
}

SxAtomInfo::ConstPtr
SxAtomicStructure::match (const SxGrid            &grid,
                          const SxAtomicStructure &toBeMatched) const
{
   SxVector<Int> mapping(toBeMatched.nTlAtoms);
   SxVector<Int>::Iterator mapIt = mapping.begin ();
   int id;

   // --- do the matching
   for (int ia = 0; ia < toBeMatched.nTlAtoms; ++ia)  {
      id = find(toBeMatched.constRef(ia), grid);

      // breaking criterium 1: atom not found
      if (id == -1) { return SxAtomInfo::ConstPtr (); }

      // breaking criterium 2: species mismatch
      if (toBeMatched.atomInfo)  {
         if (toBeMatched.getISpecies (ia) != getISpecies(id))
         { return SxAtomInfo::ConstPtr (); }
      }
      *mapIt++ = id;
   }

   // --- create atomInfo
   SxAtomInfo::Ptr newInfo = SxAtomInfo::derive(atomInfo);
   if (toBeMatched.atomInfo)  {
      // copy species info from toBeMatched
      newInfo->set (toBeMatched.atomInfo->nAtoms);
      /* do not transfer chemNames by default
      if (toBeMatched.atomInfo->meta.contains (Elements))
         newInfo->meta.attach (Elements, toBeMatched.atomInfo->meta
                                         .get<SxArray<SxString > > (Elements));
      */
   } else {
      // UNTESTED! I don't see any reasonable application, anyway.

      //newInfo->nAtoms.resize (toBeMatched.nTlAtoms);
      //Does not work if nTlAtoms < nSpecies in this structure. Should be
      newInfo->nAtoms.resize (getNSpecies ());
      //since this is the true reference in the end.

      newInfo->nAtoms.set (0);
      mapIt = mapping.begin ();
      // --- generate species info from total atom id's
      int lastSpecies = -1, iSpecies;
      for (int ia = 0; ia < toBeMatched.nTlAtoms; ++ia, ++mapIt)  {
         iSpecies = getISpecies(/* id = */ *mapIt);
         if (iSpecies < lastSpecies)  {
            // this must never happen: species misordering
            // you need to reorder toBeMatched
            cout << SX_SEPARATOR;
            cout << "| Warning:" << endl;
            cout << "| Species misordering in SxAtomicStructure::match()" << endl;
            cout << "| Ok, if you know what you are doing." << endl;
            cout << SX_SEPARATOR;
            //SX_EXIT;
         }
        
         newInfo->nAtoms(iSpecies)++;
         lastSpecies = iSpecies;
      }
      newInfo->set(newInfo->nAtoms);
   }
   newInfo->parentMap = mapping;

   return newInfo;
}

int SxAtomicStructure::find (const Coord &x, const SxGrid &grid) const
{
   SX_CHECK (grid.nAtPerGridCell.sum () == nTlAtoms,
             grid.nAtPerGridCell.sum (), nTlAtoms);

   Coord xRel = grid.gridCell.carToRel (x), delta;
   Coord dx = grid.gridCell.getBoundingBox (sqrt(3.) * epsEqual);
   SxVector3<Int> from = floor(xRel - dx),
                  to   = floor(xRel + dx),
                  meshVec;
   for (int dim=0; dim < 3; ++dim)
      if (fabs(xRel(dim)-round(xRel(dim))) < 1e-4)  {
         // atom at boundary: include neighbor grid cells
         from(dim)--; to(dim)++;
      }
   int iGrid;
   SxVector<Int>::Iterator idxIt, idxEnd;
   for (meshVec(0) = from(0); meshVec(0) <= to(0); meshVec(0)++) {
      for (meshVec(1) = from(1); meshVec(1) <= to(1); meshVec(1)++) {
         for (meshVec(2) = from(2); meshVec(2) <= to(2); meshVec(2)++) {
            // get grid cell
            iGrid = (int)grid.gridMesh.getMeshIdx(meshVec, SxMesh3D::Unknown);
            // --- loop over atoms in grid cell
            idxIt  = grid.getIdxIt (iGrid);
            idxEnd = grid.computeItEnd (iGrid);
            for ( ; idxIt != idxEnd; ++idxIt)  {
               // compute absolute position
               delta = x - constRef(*idxIt);
               grid.regCell.map (&delta, SxCell::Origin);
               // --- check distance
               if (    fabs(delta(0)) < epsEqual
                    && fabs(delta(1)) < epsEqual
                    && fabs(delta(2)) < epsEqual)
               {
                  // found atom
                  return *idxIt;
               }
            }
         }
      }
   }
   // atom not found
   return -1;
}

int SxAtomicStructure::getISpecies (int iTlAtom, int *iAtomPtr) const
{
   SX_CHECK(iTlAtom < nTlAtoms && iTlAtom >= 0, iTlAtom, nTlAtoms);
   SX_CHECK (atomInfo);
   int iSpecies = 0;
   while (iTlAtom >= atomInfo->nAtoms(iSpecies))
      iTlAtom -= atomInfo->nAtoms(iSpecies++);
   SX_CHECK (iSpecies < atomInfo->nSpecies, iSpecies, atomInfo->nSpecies);
   
   if (iAtomPtr) *iAtomPtr = iTlAtom;

   return iSpecies;
}

SxVector3<Double> SxAtomicStructure::sum () const
{
   SxVector3<Double> res(0.);
   for (int iAtom = 0; iAtom < nTlAtoms; iAtom++)
      res += constRef(iAtom);
   return res;
}

SxVector<TPrecTauR> SxAtomicStructure::absSqr () const
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (nTlAtoms > 0, nTlAtoms);
   SxVector<TPrecTauR> r2(nTlAtoms);
   SxVector<TPrecTauR>::Iterator r2It = r2.begin ();
   for (int i = 0; i < nTlAtoms; i++)
      *r2It++ = constRef(i).normSqr();

   VALIDATE_VECTOR(r2);
   return r2;
}


SxList<SymMat> SxAtomicStructure::getSymmetries () const
{
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (cell.volume >= 0.); // initialized cell

   SxList<SymMat> symmetries;
   SxList<SymMat>::Iterator sym;
   const SxAtomicStructure &original = *this;

   // get lattice symmetries
   SxList<SymMat> latSym = cell.getLatticeSymmetries ();

   // --- try them

   const int gridThreshold = 100,
             atomsPerGridCell = 3;

   const SxArray<SxString> *labels
      = atomInfo->meta.getPtr<SxArray<SxString> > (Labels);

   if (nTlAtoms < gridThreshold && !labels)  {
      for (sym = latSym.begin (); sym != latSym.end (); sym++)
         if (original == (*sym ^ original) )
            symmetries << *sym;
   } else {
      // use grids for matching
      int nGridCells =  nTlAtoms / atomsPerGridCell;
      SxGrid grid(*this, SxGrid::suggestMesh(cell,nGridCells) );
      SxConstPtr<SxAtomInfo> rotInfo;

      bool warn = true;
      for (sym = latSym.begin (); sym != latSym.end (); sym++)
         if ( rotInfo = match(grid, *sym ^ original) )  {
            bool isSym = true;
            if (labels)  {
               // --- check that all labels are the same
               for (int ia = 0; ia < nTlAtoms; ++ia)  {
                  int ja = rotInfo->parentMap(ia);
                  if ((*labels)(ia) != (*labels)(ja))  {
                     if (warn)
                        cout << "| Warning: Labels break symmetry!" << endl;
                     warn = false;
                     isSym = false;
                     break;
                  }
               }
            }
            if (isSym) symmetries << *sym;
         }
   }

   // --- check that symmetries are a group
   //cout << "| Validating symmetry group ..." << endl;
   int nSym = int(symmetries.getSize ());
   SxArray<bool> isFound (nSym);
   SymMat result;
   int i,j,k;
   for (i = 0; i < nSym; i++)  {
      isFound.set (false);
      for (j = 0; j < nSym; j++)  {
         result = symmetries(i) ^ symmetries(j);
         // find result
         for (k = 0; k < nSym; k++)  {
            if (isFound(k)) continue; // each result can be obtained only once
            if ((symmetries(k) - result).absSqr ().sum () > cell.epsSym)
               continue;
            isFound(k) = true;
            break;
         }
         if (k == nSym)  {
            cout << "Symmetry inconsistency error: symmetries are no group\n";
            cout << "Offending multiplication:" << endl;
            cout << "cartesian coordinates:" << endl;
            cout << symmetries(i) << " ^ ";
            cout << symmetries(j) << " = ";
            cout << result << endl;
            cout << "relative coordinates:" << endl;
            cout << cell.carToRel(symmetries(i)) << " ^ ";
            cout << cell.carToRel(symmetries(j)) << " = ";
            cout << cell.carToRel(result) << endl;
            cout << "This is a bug in the symmetry finding routine. Report it!";
            SX_EXIT;
         }
      } // for j
   } // for i

   return symmetries;
}

void SxAtomicStructure::print (const SxSpeciesData &speciesData) const
{
   SX_CHECK (atomInfo);
   cout << endl;
   cout << SX_SEPARATOR;
   cout << "| STRUCTURE" << endl;
   cout << SX_SEPARATOR;

   if (cell.volume > 0.)  {

      const CellMat B (cell.getReciprocalCell ());

      sxprintf ("|    a1:  %10.6f %10.6f %10.6f  Bohr\n",
                         cell(0,0), cell(1,0), cell(2,0));
      sxprintf ("|    a2:  %10.6f %10.6f %10.6f  Bohr\n",
                         cell(0,1), cell(1,1), cell(2,1));
      sxprintf ("|    a3:  %10.6f %10.6f %10.6f  Bohr\n",
                         cell(0,2), cell(1,2), cell(2,2));
      cout << "|    Omega: " << cell.volume << " Bohr^3" << endl;
      cout << SX_SEPARATOR;
      cout << "| Reciprocal lattice vectors:" << endl;
      sxprintf ("|    b1:  %10.6f %10.6f %10.6f    1/Bohr\n",
                         B(0,0), B(1,0), B(2,0));
      sxprintf ("|    b2:  %10.6f %10.6f %10.6f    1/Bohr\n",
                         B(0,1), B(1,1), B(2,1));
      sxprintf ("|    b3:  %10.6f %10.6f %10.6f    1/Bohr\n",
                         B(0,2), B(1,2), B(2,2));
      cout << SX_SEPARATOR;
      cout << "|  Definition of super cell:" << endl;
   }
   int iSpecies, iAtom;
   Coord tau;
   for (iSpecies=0; iSpecies < getNSpecies (); iSpecies++)  {
      sxprintf ("|    %s, zv=%.2f\n", speciesData.elementName(iSpecies).ascii(),
                                speciesData.valenceCharge(iSpecies));
      for (iAtom=0; iAtom < getNAtoms(iSpecies); iAtom++)  {
         tau = getAtom (iSpecies, iAtom);
         sxprintf ("|    %d:    %10.6f %10.6f %10.6f\n", iAtom + 1,
                 tau(0), tau(1), tau(2));
      }
      cout << "|" << endl;
   }
   // ---  symmetry operations
   if (cell.symGroupPtr)  {
      SxArray<SymMat> symmetries = cell.symGroupPtr->getSymmorphic ();
      cout << SX_SEPARATOR;
      cout << "| Symmetries:\n";
      int iOp, nOps = int(symmetries.getSize ());
      for (iOp=0; iOp < nOps; ++iOp)  {
         sxprintf ("|    %2d:  %s\n", iOp+1,
                 SxRotation::getType(symmetries(iOp)).identifier.ascii());
      }
   }

   cout << SX_SEPARATOR;
   cout << endl;
   cout.flush ();

}

// This is the true fprint routine
void SxAtomicStructure::fprint (FILE *output,
                                int printOptions,
                                const SxAtomicStructure &forces) const
{
   SX_CHECK (output);
   SX_CHECK (atomInfo);
   bool printForces = forces.getNAtoms () > 0;
   SX_CHECK (!printForces || (getNAtoms () == forces.getNAtoms ()),
             getNAtoms (), forces.getNAtoms ());
   SX_CHECK (!printForces || atomInfo == forces.atomInfo);
   const SxArray<SxVector3<Int> > *sticky
      = atomInfo->meta.getPtr<SxArray<SxVector3<Int> > > (StickyFilter);
   const SxArray<SxString>* elements
      = atomInfo->meta.getPtr<SxArray<SxString> >(Elements);
   bool printSticky = false;
   if (sticky)   {
      // --- check if any atom is actually fixed
      for (int i = 0; i < sticky->getSize (); ++i)
         if (!((*sticky)(i) == SxVector3<Int> (1,1,1)) )
            printSticky = true;
   }

   if (!(printOptions & CustomGroupName)) sxfprintf(output, "structure  {\n");
   sxfprintf(output, "   cell = [[% 12.8f, % 12.8f, % 12.8f],\n"
                   "           [% 12.8f, % 12.8f, % 12.8f],\n"
                   "           [% 12.8f, % 12.8f, % 12.8f]];\n",
                   cell(0,0), cell(1,0), cell(2,0),
                   cell(0,1), cell(1,1), cell(2,1),
                   cell(0,2), cell(1,2), cell(2,2));
   if (cell.epsSym != EPS_SYM_DEFAULT)
      sxfprintf(output, "   epsSym = %.4g;\n", cell.epsSym);
   if (printOptions & PrintSymmetries)  {
      SX_CHECK (cell.symGroupPtr);
      cell.symGroupPtr->fprintsx (output, true);
   }
   int is, ia, iAtom = 0;
   Coord coord;
   SxString stickyString("");
   const SxArray<SxString> *labels
      = atomInfo->meta.getPtr<SxArray<SxString> > (Labels);
   if (!printSticky) sxfprintf (output,"   movable;\n");
   for (is = 0; is < getNSpecies (); is++)  {
      sxfprintf(output, "   species  {\n");
      if (elements && is < elements->getSize() && (*elements)(is).getSize() > 0)
         sxfprintf(output, "      element=\"%s\";\n", (*elements)(is).ascii ());
      // --- atoms
      for (ia = 0; ia < getNAtoms(is); ia++, ++iAtom)  {
         coord = getAtom (is, ia);
         sxfprintf(output,"      atom {coords = [%12.8f, %12.8f, %12.8f];",
                 coord(0), coord(1), coord(2));
         // --- label
         if (labels && (*labels)(iAtom).getSize () > 0)
            sxfprintf (output," label = \"%s\";", (*labels)(iAtom).ascii ());
         // --- movable
         if (printSticky)  {
            if ((*sticky)(iAtom) == SxVector3<Int> (1, 1, 1))  {
               sxfprintf(output," movable;");
            } else  {
               for (int d = 0; d < 3; ++d)
                  if ((*sticky)(iAtom)(d))
                     sxfprintf(output, " movable%c;", 'X'+d);
            }
         }
         if (printForces)  {
            coord = forces.getAtom (is, ia);
            sxfprintf(output,"\n            force  = [% .10f,% .10f,% .10f];",
                    coord(0), coord(1), coord(2));
         }
         sxfprintf(output, " }\n");
      }
      sxfprintf(output, "   }\n");
   }
   if (!(printOptions & CustomGroupName)) sxfprintf(output, "}\n");
}

void SxAtomicStructure::fprintRel (FILE *output) const
{
   SX_CHECK (output);
   SX_CHECK (atomInfo);
   sxfprintf(output, "structure  {\n");
   sxfprintf(output, "   cell = [[%12.8f, %12.8f, %12.8f],\n"
                   "           [%12.8f, %12.8f, %12.8f],\n"
                   "           [%12.8f, %12.8f, %12.8f]];\n",
                   cell(0,0), cell(1,0), cell(2,0),
                   cell(0,1), cell(1,1), cell(2,1),
                   cell(0,2), cell(1,2), cell(2,2));
   Coord coord;
   sxfprintf (output,"   movable;\n");
   const SxArray<SxString> *labels
      = atomInfo->meta.getPtr<SxArray<SxString> > (Labels);
   const SxArray<SxString>* elements 
      = atomInfo->meta.getPtr<SxArray<SxString> >(Elements);
   for (int is = 0, iTlAtom = 0; is < getNSpecies (); is++)  {
      sxfprintf(output, "   species  {\n");
      if (elements && (*elements)(is).getSize () > 0)
         sxfprintf(output, "      element=\"%s\";\n", (*elements)(is).ascii());
      for (int ia = 0; ia < getNAtoms(is); ia++, iTlAtom++)  {
         coord = cell.carToRel(getAtom (is, ia));
         sxfprintf(output,
                 "      atom {coords = [% .12f, % .12f, % .12f]; relative;",
                 coord(0), coord(1), coord(2));
         if (labels && (*labels)(iTlAtom).getSize () > 0)
            sxfprintf (output," label = \"%s\";", (*labels)(iTlAtom).ascii ());
         sxfprintf (output, " }\n");
      }
      sxfprintf(output, "   }\n");
   }
   sxfprintf(output, "}\n");
}

void SxAtomicStructure::write (const SxString &fileName) const
{
   try {
      SxBinIO io(fileName, SxBinIO::BINARY_WRITE_ONLY);
      write (io);
      io.setMode (SxBinIO::WRITE_DATA);
      write (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxAtomicStructure::write (SxBinIO &io) const
{
   SX_CHECK (   io.mode == SxBinIO::BINARY_WRITE_ONLY
             || io.mode == SxBinIO::BINARY_WRITE_LARGE
             || io.mode == SxBinIO::BINARY_WRITE_PARALLEL, io.mode);
   SX_CHECK (creationStatus == Finished);
   SX_CHECK (atomInfo);

   // --- put chemical symbols into one string
   SxString elemString;
   if (atomInfo->meta.contains (Elements))  {
      const SxArray<SxString> &chemName = atomInfo->meta.get (Elements);
      if (chemName.getSize () == getNSpecies ())  {
         for (int iSpecies = 0; iSpecies < getNSpecies (); iSpecies++)
            elemString+= chemName(iSpecies) + SxString(',');
         // remove final ','
         elemString.resize (elemString.getSize () - 1, true);
      }
   }

   SxString labels;
   if (atomInfo->meta.contains (Labels))  {
      SxArray<SxString> &theLabels = atomInfo->meta.get (Labels);
      labels = SxString::join (theLabels, '|');
   }

   try  {
      // --- create dimensions
      io.addDimension ("nSpecies", getNSpecies ());
      io.addDimension ("nAllAtoms", nTlAtoms);
      io.addDimension ("xyz", 3);

      // TODO  Verify data distrition for parallel IO.
      if ((io.ncMode == SxBinIO::WRITE_HEADER) || (SxLoopMPI::me() == 0))
      {
         // --- write data
         io.write ("nAtoms", atomInfo->nAtoms, "nSpecies");
         io.write ("tau", SxMatrix<TPrecTauR> (coords.transpose ()),
               "nAllAtoms", "xyz");
         // write species names if available
         if (elemString.getSize () > 0) io.write ("chemNames", elemString);
         // write labels if available
         if (labels.getSize () > 0) io.write ("atomLabels", labels);
         // write cell if set up
         if (cell.volume > 0.) cell.write (io);
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxAtomicStructure::read (const SxBinIO &io)
{
   SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);
   SxPtr<SxAtomInfo> aInfo = SxAtomInfo::create ();

   try  {
      // read cell
      if (io.contains("cell")) cell.read (io);

      // --- get dimensions
      aInfo->nSpecies = io.getDimension ("nSpecies");
      nTlAtoms = io.getDimension ("nAllAtoms");

      // --- read data
      aInfo->nAtoms.resize (aInfo->nSpecies);
      io.read ("nAtoms", &aInfo->nAtoms, aInfo->nSpecies);
      aInfo->setupOffset ();

      if (io.contains ("atomLabels"))  {
         SxString labels;
         io.read ("atomLabels", &labels);
         SxArray<SxString> labelList = labels.tokenize ('|', true);
         aInfo->meta.attach (Labels, labelList);
         if (labelList.getSize () != nTlAtoms)  {
            cout << labels << " => " << labelList << endl;
            cout << "Inconsistent labels in sxb file!" << endl;
            SX_EXIT;
         }
      }
      if (io.contains ("chemNames"))  {
         SxString elemString;
         io.read ("chemNames", &elemString);
         SxArray<SxString> chemNames = elemString.tokenize (',');
         chemNames.resize (aInfo->nSpecies, true);
         aInfo->meta.attach (Elements, chemNames);
      }

      coords.reformat (nTlAtoms, 3);
      io.read ("tau", &coords, nTlAtoms, 3);
      coords = coords.transpose ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // --- clear unused vars
   creationStatus = Finished;
   posList.resize (0);
   atomInfo = aInfo;

   // compute symmetries if not yet there
   if (cell.volume > 0. && !cell.symGroupPtr)
      cell.symGroupPtr = SxPtr<SxSymGroup>::create(getSymmetries ());

}

// --- sorting

/** \brief Sortable vector (coordinate-wise: z, y, x)
  */
class SxSortableVector3
{
   public:
   SxVector3<Double> vec;
   double epsEqual;
   /// Constructor from vector
   SxSortableVector3 (SxVector3<Double> v,
                      const double epsE = SX_EPS_STRUCT_DEFAULT)
      : vec (v), epsEqual (epsE) {}
   /// Empty constructor for lists etc.
   SxSortableVector3 () {}
   /// Smaller operator
   bool operator< (SxSortableVector3 v) const
   {
      if (this->vec(2) < v.vec(2) - epsEqual) return true;
      if (this->vec(2) > v.vec(2) + epsEqual) return false;
      if (this->vec(1) < v.vec(1) - epsEqual) return true;
      if (this->vec(1) > v.vec(1) + epsEqual) return false;
      if (this->vec(0) < v.vec(0) - epsEqual) return true;
      return false;
   }
   /// Larger operator
   bool operator> (const SxSortableVector3 &v) const
   {
      return (v < *this);
   }
};

void SxAtomicStructure::sort (bool derive)
{
   SX_CHECK (creationStatus == Finished); // not that logical ...
   SxSortedList<SxSortableVector3> newAtoms;
   SxList<SxSortableVector3>::Iterator it;
   SxList<int> idMap;
   SxAtomInfo::Ptr newInfo;
   if (derive)  {
      newInfo = SxAtomInfo::derive (atomInfo);
      newInfo->nAtoms.copy (atomInfo->nAtoms);
      newInfo->setupOffset ();
      newInfo->meta = atomInfo->meta;
      newInfo->parentMap.resize (nTlAtoms);
   }
   for (int is = 0; is < getNSpecies (); is++)  {
      int nAtom = getNAtoms(is);
      // sort by z,y,x
      newAtoms.resize (0);
      for (int ia = 0; ia < nAtom; ia++)  {
         ssize_t pos = newAtoms.append (SxSortableVector3((*this)(is,ia), epsEqual));
         if (derive) { idMap.insert (pos, ia + atomInfo->offset(is)); }
      }
      // copy sorted list
      it = newAtoms.begin ();
      for (int ia = 0; ia < nAtom; ia++,++it)
         setAtom(is, ia, (*it).vec);
      if (derive)  {
         for (int ia = 0; ia < nAtom; ++ia)
            newInfo->parentMap(ia + atomInfo->offset(is)) = idMap(ia);
         idMap.resize (0);
      }
   }
   if (derive) atomInfo = newInfo;
}

SxAtomicStructure
SxAtomicStructure::operator* (const SxVector<TPrecTauR> &v) const
{
   SX_CHECK (creationStatus == Finished);
   if (v.getSize () == nTlAtoms)  {
      // --- atom-specific factors
      SxAtomicStructure res (getNewStr ());
      SxVector<TPrecTauR>::Iterator vIt = v.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia)
         res.ref(ia) = constRef(ia) * *vIt++;
      return res;
   } else {
      // --- species-specific factors
      SX_CHECK (atomInfo);
      SX_CHECK (v.getSize() == getNSpecies(), v.getSize(), getNSpecies());

      SxAtomicStructure res (*this, SxAtomicStructure::Copy);
      int iSpecies, nSpecies = getNSpecies ();
      SxIdx range;
      for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
         // get atom index range
         range = getRange (iSpecies);
         // switch to single coordinate ranges
         range.start *=3;
         range.end = range.end * 3 + 2;
         // scale these coordinates
         res.coords(range) *= v(iSpecies);
      }

      return res;
   }
}

SxAtomicStructure
SxAtomicStructure::operator/ (const SxVector<TPrecTauR> &v) const
{
   SX_CHECK (creationStatus == Finished);
   if (v.getSize () == nTlAtoms)  {
      // --- atom-specific factors
      SxAtomicStructure res (getNewStr ());
      SxVector<TPrecTauR>::Iterator vIt = v.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia)
         res.ref(ia) = constRef(ia) / *vIt++;
      return res;
   } else {
      // --- species-specific factors
      SX_CHECK (atomInfo);
      SX_CHECK (v.getSize() == getNSpecies(), v.getSize(), getNSpecies());

      SxAtomicStructure res (*this, SxAtomicStructure::Copy);
      int iSpecies, nSpecies = getNSpecies ();
      SxIdx range;
      for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
         // get atom index range
         range = getRange (iSpecies);
         // switch to single coordinate ranges
         range.start *=3;
         range.end = range.end * 3 + 2;
         // scale these coordinates
         res.coords(range) /= v(iSpecies);
      }

      return res;
   }
}

SxAtomicStructure &
SxAtomicStructure::operator*= (const SxVector<TPrecTauR> &v)
{
   SX_CHECK (creationStatus == Finished);
   if (v.getSize () == nTlAtoms)  {
      // --- atom-specific factors
      SxVector<TPrecTauR>::Iterator vIt = v.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia)
         ref(ia) *=  *vIt++;
   } else {
      // --- species-specific factors
      SX_CHECK (atomInfo);
      SX_CHECK (v.getSize() == getNSpecies(), v.getSize(), getNSpecies());

      int iSpecies, nSpecies = getNSpecies ();
      SxIdx range;
      for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
         // get atom index range
         range = getRange (iSpecies);
         // switch to single coordinate ranges
         range.start *=3;
         range.end = range.end * 3 + 2;
         // scale these coordinates
         coords(range) *= v(iSpecies);
      }
   }
   return *this;
}

SxAtomicStructure &
SxAtomicStructure::operator/= (const SxVector<TPrecTauR> &v)
{
   SX_CHECK (creationStatus == Finished);
   if (v.getSize () == nTlAtoms)  {
      // --- atom-specific factors
      SxVector<TPrecTauR>::Iterator vIt = v.begin ();
      for (int ia = 0; ia < nTlAtoms; ++ia)
         ref(ia) /=  *vIt++;
   } else {
      // --- species-specific factors
      SX_CHECK (atomInfo);
      SX_CHECK (v.getSize() == getNSpecies(), v.getSize(), getNSpecies());

      int iSpecies, nSpecies = getNSpecies ();
      SxIdx range;
      for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
         // get atom index range
         range = getRange (iSpecies);
         // switch to single coordinate ranges
         range.start *=3;
         range.end = range.end * 3 + 2;
         // scale these coordinates
         coords(range) /= v(iSpecies);
      }
   }
   return *this;
}

SxAtomicStructure operator- (const Coord &trans, const SxAtomicStructure &str)
{
   SX_CHECK (! str.isCreationMode ());

   SxAtomicStructure res (str.nTlAtoms, str.atomInfo);
   res.cell = str.cell;
   res.epsEqual = str.epsEqual;
   for (int iAtom = 0; iAtom < str.getNAtoms (); ++iAtom)
      res.ref(iAtom) = trans - str.constRef(iAtom);
   return res;
}

std::ostream &operator<< (std::ostream &s, const SxAtomicStructure &str)
{
   SX_CHECK (!str.isCreationMode ());
   if (str.nTlAtoms == 0) return s;
   if (str.atomInfo)  {
      for (int is = 0; is < str.getNSpecies (); ++is)  {
         s << "species " << (is+1) << std::endl;
         for (int ia = 0; ia < str.getNAtoms (is); ++ia)  {
            s << "atom " << (ia+1) << ": " << str(is,ia) << std::endl;
         }
      }
   } else {
      for (int ia = 0; ia < str.nTlAtoms; ++ia)
         s << "atom " << (ia+1) << ": " << str(ia) << std::endl;
   }
   return s;
}

SxVector<TPrecTauR> SxAtomicStructure::operator^(const Coord &x) const
{
   if (nTlAtoms == 0) return SxVector<TPrecTauR> ();
   SxVector<TPrecTauR> res(nTlAtoms);
   SxVector<TPrecTauR>::Iterator it = res.begin ();
   for (int ia = 0; ia < nTlAtoms; ++ia)
      *it++ = x ^ constRef(ia);
   VALIDATE_VECTOR (res);
   return res;
}

SxVector<TPrecTauR>
SxAtomicStructure::operator^(const SxAtomicStructure &x) const
{
   if (nTlAtoms == 0 || x.nTlAtoms) return SxVector<TPrecTauR> ();
   SX_CHECK (this != &x); // use absSqr () instead
   SxVector<TPrecTauR> res(nTlAtoms, x.nTlAtoms);
   SxVector<TPrecTauR>::Iterator it = res.begin ();

   for (int ja = 0; ja < x.nTlAtoms; ++ja)
      for (int ia = 0; ia < nTlAtoms; ++ia)
         *it++ = constRef(ia) ^ x.constRef(ja);

   VALIDATE_VECTOR (res);
   return res;
}

SxCell SxAtomicStructure::getPrimitiveCell () const
{
   SxList<Coord> translations;
   // always add the current basis vectors
   // this garantuees a fall-back to the current cell
   translations << cell.basis(0) << cell.basis(1) << cell.basis(2);

   int is = -1, ia, nAtoms;
   if (atomInfo) {
      // pick species with fewest atoms
      nAtoms = atomInfo->nAtoms.minval (&is);
      SX_CHECK (is >= 0 && is < getNSpecies () );
   } else {
      nAtoms = nTlAtoms;
   }

   cout << SX_SEPARATOR
        << "| Primitive cell search:" << endl
        << SX_SEPARATOR;
   cout.flush ();
   // try out all possible translation between 1st and n-th atom
   Coord t;
   SxCell primCell = cell;
   SxGrid grid(*this, 3);
   int nRealChecks = 0, nPrm = 0, nRed = 0, nFail = 0, nComm = 0;
   SxList<Coord> failed;
   for (ia = 1; ia < nAtoms; ia++)  {

      if (translations.getSize () > 20)  {
         // construct preliminary primitive cell
         int nDim = primCell.setFromVectors (translations);
         if (nDim < 3)  {
            cout << translations;
            cout << "Internal error: primitive cell collapsed" << endl;
            SX_EXIT;
         }
         // new translation list
         translations.resize (0);
         for (int d = 0; d < 3; ++d)
            translations << cell.basis(d) << primCell.basis(d);
         // cout << "ia=" << (ia+1) << "; primCell=" << primCell << endl;
         // cout << nRealChecks << "/" << nPrm << '/' << nRed << '/'
         //     << nFail << endl;

         // clean up failed list, retain only unique failures
         SxList<Coord>::Iterator iF, jF;
         int i;
         bool known;
         for (i = 0, iF = failed.begin (); iF != failed.end (); ++iF, ++i)  {
            known = false;
            for (jF = failed.begin (); !known && jF != iF; ++jF)
               known = (primCell.getMapped(*iF -*jF, SxCell::Origin)
                       .normSqr () < epsEqual);
            if (known)
               failed.remove (i--); // TODO: failed.remove (*iF);
         }
      }

      // get translation
      if (is == -1)
         t = getAtom(ia) - getAtom(0);
      else
         t = getAtom(is,ia) - getAtom(is,0);
      cell.map (&t, SxCell::Origin);

      // check for zero translation
      if (t.normSqr () < epsEqual)  {
         Coord a, b;
         if (is == -1)  {
            a = getAtom(0); b = getAtom(ia);
         } else {
            a = getAtom(is,0); b = getAtom(is,ia);
         }
         cout << "Two atoms found on the same place when trying to "
                 "find primitive cell." << endl;
         cout << "Lattice: " << cell << endl;
         cout << "Atom A(is=" << is << " ia=" << 0  << "): " << a << endl;
         cout << "Atom B(is=" << is << " ia=" << ia << "): " << b << endl;
         SX_QUIT;
      }

      if (primCell.getMapped (t, SxCell::Origin).normSqr () < epsEqual)  {
         // this is part of the current guess for primitive lattice
         nPrm++;
         continue;
      }

      // --- Do simple commensurability test
      // Ap ... the basis of the primitive cell
      // As ... the basis of the supercell
      // S  ... supercell representation in relative coordinates (integer)
      //    As = Ap ^ S    <=>    Ap = As ^ (S^-1)
      // det(S) is number of primitive cells in supercell
      // NB: det(S) ^ (S^-1) is integer (follows from determinantal
      //     representation of inverse)
      // obviously nAtoms = N * det(S)
      // where N is the number of atoms in the primitive cell
      // Any lattice vector of the primitive cell is given by
      //    t = Ap ^ n
      // Then
      //    nAtoms * t = As ^ (N * det(S) * S^-1)
      // is a lattice vector of the supercell.
      if (cell.getMapped(nAtoms * t, SxCell::Origin).normSqr () > epsEqual)  {
         nComm++;
         continue;
      }

      {
         // Test if t can be constructed from known translations
         // by reducing it to zero when subtracting known translations
         // this helps a lot once the primitive vectors are found
         Coord u = t;
         SxList<Coord>::Iterator it  = translations.begin (),
                                 end = translations.end ();
         for (; it != end; ++it)  {
            if ((*it).normSqr () < 1.999 * (u ^ (*it)))  {
               // the 1.999 instead of 2 avoids numerical noise problems
               // |v|^2 < 2 u^v <=> |u-v|^2 < |u|^2
               // replace u by u-v
               u -= *it;
               if (u.normSqr () < epsEqual) break;
               // --- start over again
               it = translations.begin ();
            }
         }
         if (u.normSqr () < epsEqual)  {
            // remove numerical noise by enforcing commensurability (see above)
            t -= cell.getMapped(nAtoms * t, SxCell::Origin) / nAtoms;
            translations << t << (-t);
            nRed++;
            continue;
         }
      }

      // Check for known failures
      {
         SxList<Coord>::Iterator it, end = failed.end ();
         bool unknown = true;
         for (it = failed.begin (); unknown && it != end; ++it)
            unknown = (primCell.getMapped(t - *it, SxCell::Origin)
                      .normSqr () >= epsEqual);
         if (!unknown) {
            nFail++;
            continue;
         }
      }

      // do real check
      nRealChecks++;
      if (match(grid, (*this) + t))  {
         // yippieh, this is a translational symmetry
         // remove numerical noise by enforcing commensurability (see above)
         t -= cell.getMapped(nAtoms * t, SxCell::Origin) / nAtoms;
         translations << t << (-t);
      } else {
         failed << t;
      }
   }
   primCell.setFromVectors (translations);
   cout << "| Deduced success (prim-test): " << nPrm << endl
        << "| Deduced success (red-test) : " << nRed << endl
        << "| Deduced failure (comm-test): " << nComm << endl
        << "| Deduced failure (fail-list): " << nFail << endl
        << "| Checked                    : " << nRealChecks << endl
        << "| Total                      : " << (nAtoms - 1) << endl;
   cout << "| Supercell contains " << round(cell.volume / primCell.volume)
        << " primitive cells." << endl;
   (cout << SX_SEPARATOR).flush ();

   // --- improve accuracy by resetting primCell from cell
   {
      SxMatrix3<Double> rel; // supercell in relative coordinates of primitive cell
      for (int d = 0; d < 3; ++d)  {

         // without checks:
         // rel.setCol(d, round(primCell.carToRel(cell.basis(d))));

         // --- with consistency check
         rel.setCol(d, primCell.carToRel(cell.basis(d)));
         for (int i = 0; i < 3; ++i)  {
            if (fabs(rel(i,d) - round(rel(i,d))) > 1e-2 /* generous */)  {
               cout << "Internal accuracy error in "
               << __FILE__ << ':' << __LINE__ << endl
               << "Please inform the developers; attach the input file used!"
               << endl;
               SX_EXIT;
            }
            rel(i,d) = round(rel(i,d)); // throw away numerical noise
         }

      }
      // --- now try to align the primitive cell vectors to input cell
      {
         rel = rel.transpose ();
         for (int i = 0; i < 3; i++)  {
            for (int j = (i + 1) % 3; j != i; j = (j + 1) % 3)  {
               if (fabs(rel(i,j)) > 1e-10)  {
                  double others = 0.;
                  for (int k = 0; k < 3; ++k)  {
                     if (k != i) others += fabs(rel(k,j));
                     if (k != j) others += fabs(rel(i,k));
                  }
                  if (others < 1e-10)  {
                     // swap rows i and j
                     Coord oldI = rel.col(i);
                     rel.setCol (i, rel.col(j));
                     rel.setCol (j, oldI);
                  }
               }
            }
         }
         // --- optimize orientation
         for (int i = 0; i < 3; i++)
            if (rel(i,i) < -1e-10) rel.setCol(i,-rel.col(i));
         rel = rel.transpose ();
      }
      // set primitive cell from supercell, including the symmetryGroupPtr
      primCell = cell ^ rel.inverse ();
      primCell.symGroupPtr = cell.symGroupPtr;

      if (primCell.determinant () < -1e-10)  {
         // --- try to get right-handed orientation by switching sign
         //     of non-aligned basis vector. If all vectors are aligned,
         //     the left-handed orientation of the input cell is kept.
         for (int i = 0; i < 3; ++i)  {
            if (fabs(rel(i,(i+1) % 3)) + fabs(rel(i,(i+2) % 3)) > 1e-10)  {
               for (int j = 0; j < 3; ++j) primCell(j,i) = -primCell(j,i);
               primCell.setup ();
               break;
            }
         }
      }
   }

   return primCell;
}

SxAtomicStructure SxAtomicStructure::getPrimStr (const SxCell &primCell) const
{
   // check that primCell is really a primitive cell
   for (int d = 0; d < 3; ++d)  {
      Coord a = primCell.carToRel(cell.basis(d));
      if ((round(a)-a).normSqr () > epsEqual)  {
         cout << "cell is no supercell of primCell" << endl;
         cout << "cell basis " << (d+1) << " is " << a
              << " in rel. coords of primCell " << endl;
         SX_EXIT;
      }
   }

   SxAtomicStructure str(*this, Copy);
   str.cell = primCell;
   str.map ();

   int nPrim       = int(lround(cell.volume/primCell.volume)),
       nAtomsPrim  = getNAtoms () / nPrim;
   SX_CHECK (getNAtoms () % nAtomsPrim == 0, getNAtoms (), nAtomsPrim);

   SxPtr<SxAtomInfo> aInfo;
   if (atomInfo)  {
      aInfo = SxPtr<SxAtomInfo>::create (getNSpecies ());
      for (int is = 0; is < getNSpecies (); ++is)
         aInfo->nAtoms(is) = atomInfo->nAtoms(is) / nPrim;
      // aInfo->set(atomInfo->nAtoms / nPrim);
      aInfo->setupOffset ();
      SX_CHECK (aInfo->nAtoms.sum () == nAtomsPrim,
                aInfo->nAtoms.sum (), nAtomsPrim);
      if (atomInfo->meta.contains (Elements))  {
         aInfo->meta.attach (Elements, atomInfo->meta
                                       .get<SxArray<SxString> > (Elements));
      }
   }

   SxAtomicStructure res(primCell, nAtomsPrim, aInfo);
   res.epsEqual = epsEqual;

   SxGrid grid(str, 10);
   SxNeighbors nn;
   SxArray<bool> done(getNAtoms ());
   done.set(false);
   double maxDist = 3. * epsEqual; // i.e. small
   SxVector<Int>::Iterator it, end;
   int iAtomPrim = 0;
   for (int ia = 0; ia < getNAtoms (); ++ia)  {
      if (done(ia)) continue;

      // find equivalent atoms
      nn.compute (grid, str, getAtom(ia),  maxDist,
                  SxNeighbors::StoreIdx | SxNeighbors::IncludeZeroDistance);
      SX_CHECK (nn.getSize () == nPrim, nn.getSize (), nPrim);
      end = nn.idx.end ();
      for (it = nn.idx.begin (); it != end; ++it) done(*it) = true;

      // set primitive structure from any representative
      res.ref(iAtomPrim) = constRef(nn.idx(0));
      primCell.map(&res.ref(iAtomPrim));
      ++iAtomPrim;
   }
   SX_CHECK (iAtomPrim == nAtomsPrim, iAtomPrim, nAtomsPrim);
   return res;

}

SxAtomicStructure SxAtomicStructure::splitSpecies (const SxCell &redCell) const
{
   // check that primCell is really a primitive cell
   for (int d = 0; d < 3; ++d)  {
      Coord a = redCell.carToRel(cell.basis(d));
      if ((round(a)-a).normSqr () > epsEqual)  {
         cout << "cell is no supercell of primCell" << endl;
         cout << "cell basis " << (d+1) << " is " << a
              << " in rel. coords of primCell " << endl;
         SX_EXIT;
      }
   }

   int redFactor = int(lround(cell.volume/redCell.volume));
   int nSpecNew = getNAtoms () / redFactor;

   SxVector<Int> nAtomsNew(nSpecNew);
   nAtomsNew.set (0);
   SxPtr<SxAtomInfo> newInfo = SxAtomInfo::derive (atomInfo);
   newInfo->parentMap.resize (getNAtoms ());
   
   SxArray<Coord> refAtoms(nSpecNew);
   SxAtomicStructure res(cell, getNAtoms ());
   res.epsEqual = epsEqual;
   int nRef = 0;
   for (int ia = 0; ia < getNAtoms (); ++ia)  {
      for (int ja = 0; ja <= nRef; ++ja)  {
         if (ja == nRef)  {
            // --- new reference atom
            refAtoms(nRef++) = constRef(ia);
            int newIdx = ja * redFactor + nAtomsNew(ja)++;
            // set map
            newInfo->parentMap(newIdx) = ia;
            // set coordinate
            res.ref(newIdx) = constRef(ia);
            break;
         }
         Coord diff = constRef(ia) - refAtoms(ja);
         redCell.map (&diff, SxCell::Origin);
         if (diff.norm () < epsEqual) {
            SX_CHECK (nAtomsNew(ja) < redFactor, nAtomsNew(ja), redFactor);
            int newIdx = ja * redFactor + nAtomsNew(ja)++;
            // set map
            newInfo->parentMap(newIdx) = ia;
            // set coordinate
            res.ref(newIdx) = constRef(ia);
            break;
         }
      }
   }
   SX_CHECK ((nAtomsNew + -redFactor).sqr ().sum () == 0);	
   newInfo->set (nAtomsNew);
   res.atomInfo = newInfo;
   return res;
}


SxSymGroup SxAtomicStructure::getSymGroup (const SymSearchMode mode, const
      SxCell &inCell) const
{
   SX_CHECK (!isCreationMode ());
   SX_CHECK (atomInfo); // TODO: not really needed

   // --- 1st step: get the primitive cell
   SxCell primCell = (inCell.volume > 0. ? inCell : getPrimitiveCell ());

   // --- 2nd step: find primitive structure symmetries
   SxAtomicStructure primStr = getPrimStr (primCell), rotStr;
   SX_CHECK (mode == FullSymmetries || mode == SupercellCompatible);
   SxArray<SymMat> rots = (mode == FullSymmetries)
                        ? primCell.getLatticeSymmetries ()
                        : cell.getLatticeSymmetries ();
   SxGrid grid(primStr, 3);
   SxList<SxSymOp> syms;
   Coord t;

   // cout << "eps=" << epsEqual << endl;

   // work on smallest species
   int is;
   int nAtoms = primStr.atomInfo->nAtoms.minval(&is);

   for (int iSym = 0; iSym < rots.getSize (); ++iSym)  {
      // --- capture rare case where the primitive cell has reduced
      //     symmetries compared to full cell's lattice
      if (!primCell.isLatticeSymmetry (rots(iSym))) continue;

      // rotate
      rotStr = rots(iSym) ^ primStr;
      // try possible translations
      for (int ia = 0; ia < nAtoms; ++ia)  {
         // idea: make rotStr(is,0) equivalent to (is,ia)
         t = primStr.constRef(is, ia) - rotStr.constRef(is, 0);
         t = primCell.getMapped(t,SxCell::Origin);
         // important to identify symmorphic symmetries as such:
         if (t.normSqr () < 3. * sqr(epsEqual)) t.set (0.);
         // check for symmetry matching
         if (primStr.match(grid, rotStr + t))
            syms << SxSymOp(rots(iSym), t);
         /*
         else
            cout << rots(iSym) << " " << t << " " << (ia+1) << endl
                 << (rotStr + t) << endl;
         */
      }
   }
   return SxSymGroup (primCell, syms, cell);
}

SxAtomicStructure SxAtomicStructure::getInequivalentAtoms (SxArray<int> *equivalentListPtr) const
{
   SxSymGroup symGroup = getSymGroup (SxAtomicStructure::SupercellCompatible);
   symGroup.setToPrimitive ();
   SxCell primCell = symGroup.getPrimCell ();
   primCell.symGroupPtr = SxPtr<SxSymGroup>::create (symGroup);
   return getInequivalentAtoms (primCell, equivalentListPtr);
}

SxAtomicStructure SxAtomicStructure::getInequivalentAtoms (
      const SxCell &primCell,
      SxArray<int> *equivalentListPtr) const
{
   SX_CHECK (primCell.symGroupPtr);
   SxAtomicStructure primStr = getPrimStr (primCell);
   SxGrid grid (primStr, 3);

   int nAtom = primStr.getNAtoms ();

   SxArray<bool> done(nAtom);
   done.set (false);
   SxArray<int> primEquivList (nAtom);
   int idx, nSym = primCell.symGroupPtr->getSize ();
   SxStack<int>   idxList;
   SxStack<Coord> ineqAtoms;
   int nInEqAtoms = 0;
   SxSymOp symOp;
   for (int ia = 0; ia < nAtom; ++ia)  {
      if (done(ia)) continue;
      for (int iSym = 0; iSym < nSym; ++iSym)  {
         symOp = primCell.symGroupPtr->operator() (iSym);
         idx = primStr.find (symOp ^ primStr.constRef (ia), grid);
         if (idx < 0 || idx < ia) {
            primCell.symGroupPtr->print ();
            cout << "idx= " << (idx+1) << " < ia= " << (ia+1) << endl;
            cout << "sym " << (iSym+1) << ": " << symOp << endl;
            cout << primStr.getAtom (ia) << " -> "
                 << (symOp ^ primStr.getAtom (ia))
                 << endl;
            // oops, this is not a symmetry of this structure
            cout << "Internal error: symmetry mismatch" << endl;
            cout << "Please inform the developers." << endl;
            SX_EXIT;
         }
         if (idx == ia && !done(ia))  {
            nInEqAtoms++;
            idxList << ia;
            ineqAtoms << primStr.constRef (ia);
         }
         primEquivList(idx) = ia;
         done(idx) = true;
      }
   }

   // --- map primitive equivalent list to full cell
   if (equivalentListPtr != NULL)  {
      SxArray<int> &equivList = *equivalentListPtr;
      for (int ia = 0; ia < getNAtoms(); ++ia)  {
         idx = primStr.find (constRef(ia),grid);
         Coord idxPos = primStr.constRef(primEquivList(idx));
         SxGrid grid2 (*this, 3);
         int iAtom = find(idxPos,grid2);
         equivList(ia) = iAtom;
      }
      // --- sort list
      for (int ia = 0; ia < getNAtoms(); ++ia)  {
         int equivAtom = equivList(ia);
         if (equivAtom > ia) {
            int newRef = ia;
            for (int ja = 0; ja < getNAtoms(); ++ja)  {
               if (equivList(ja) == equivAtom)
                  equivList(ja) = newRef;
            }
         }
      }
   }


   SxAtomicStructure res (primCell, nInEqAtoms);
   res.epsEqual = epsEqual;
   // copy coordinate data
   res.importStack (ineqAtoms);

   if (primStr.atomInfo)  {
      // --- create species and mapping info
      SxPtr<SxAtomInfo> info = SxAtomInfo::derive (primStr.atomInfo);
      res.atomInfo = info;

      if (primStr.atomInfo->meta.contains (Elements))
         info->meta.attach (Elements, primStr.atomInfo->meta
                                      .get <SxArray<SxString> > (Elements));

      // mapping
      info->parentMap = idxList;

      // species
      info->nAtoms.set (0);
      for (int ia = 0; ia < nInEqAtoms; ++ia)
         info->nAtoms (getISpecies(ia))++;
      info->setupOffset ();
   }
   return res;
}

SxAtomicStructure SxAtomicStructure::cut (const Coord lower, const Coord upper)  const
{
   SxAtomicStructure cutStructure;
   cutStructure.cell = 1.0 * cell;
   for (int iSpecies = 0; iSpecies < getNSpecies (); iSpecies++)   {
      cutStructure.newSpecies ();
      for (int iAtom = 0; iAtom < getNAtoms(iSpecies); iAtom++)   {
         bool add = true;
         Coord diff = getAtom(iSpecies, iAtom) - lower;
         for (int iDir = 0; iDir < 3; iDir++) if (diff(iDir) < 0) add = false;
         diff = upper - getAtom(iSpecies, iAtom);
         for (int iDir = 0; iDir < 3; iDir++) if (diff(iDir) < 0) add = false;
         if (add) cutStructure.addAtom( (*this)(iSpecies, iAtom));
      }
   }
   cutStructure.endCreation ();

   return cutStructure;

}
