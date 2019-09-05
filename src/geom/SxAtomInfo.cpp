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

#include<SxAtomInfo.h>

void SxAtomInfo::resize (int newNSpecies, bool keep)
{
   SX_CHECK (newNSpecies >= 0, newNSpecies);
   nSpecies = newNSpecies;
   nAtoms.resize (nSpecies, keep, 0);
   offset.resize (nSpecies, keep, 0);
   if (!keep)  { 
      nAtoms.set (0);
      offset.set(0); 
   }
}

void SxAtomInfo::set(const SxVector<Int> &nAtomsIn)
{
   SX_CHECK(nAtomsIn.getSize () > 0, nAtomsIn.getSize ());
   
   nSpecies = (int)nAtomsIn.getSize ();
   
   if (&nAtomsIn != &nAtoms)  {
      // --- copy nAtoms
      nAtoms.resize (nSpecies);
      nAtoms.copy (nAtomsIn);
   }

   setupOffset ();
}

void SxAtomInfo::setupOffset ()
{
   // setup offset
   offset.resize (nSpecies);

   int sum = 0;
   for (int is = 0; is < nSpecies; ++is)  {
      offset(is) = sum;
      sum += nAtoms(is);
   }
   
}

bool SxAtomInfo::operator== (const SxAtomInfo &in) const
{
   if (&in == this) return true;
   if (nAtoms.getSize () != in.nAtoms.getSize ()) return false;
   if (nAtoms.getSize () == 0) return true;

   SxVector<Int>::Iterator it1 = nAtoms.begin (),
                           it2 = in.nAtoms.begin ();
   for ( ; it1 != nAtoms.end (); ++it1, ++it2)
      if (*it1 != *it2) return false;

   return true;
      
}

SxAtomInfo::Ptr SxAtomInfo::derive(const SxAtomInfo::ConstPtr &from)
{
   SX_CHECK(from);
   Ptr res = Ptr::create (from->nSpecies);
   res->parent = from;
   return res;
}

SxAtomInfo::Ptr 
SxAtomInfo::getCopy(const SxAtomInfo::ConstPtr &newParent) const
{
   SX_CHECK (!newParent || isChild(newParent));
   Ptr res = Ptr::create (nSpecies);
   res->nAtoms.copy (nAtoms);
   res->offset.copy (offset);
   res->meta = meta;
   res->parent = newParent;
   if (newParent == parent)  {
      res->parentMap.copy (parentMap);
   } else if (newParent)  {
      res->parentMap = getIdxMap (newParent);
   } else {
      // no parent -> no map
   }
   return res;
}

bool SxAtomInfo::isChild(const SxAtomInfo::ConstPtr &in) const
{
   ConstPtr ptr = parent;
   while (ptr)  {
      SX_CHECK (ptr.getPtr () != this);
      if (ptr == in) return true;
      ptr = ptr->parent;
   }
   return false;
}

SxVector<Int> SxAtomInfo::getIdxMap(const SxAtomInfo::ConstPtr &in) const
{
   if (parent == in) return parentMap;
   
   SX_CHECK(isChild(in));
   SX_CHECK(parentMap.getSize () == nAtoms.sum ());
   
   SxVector<Int> res;
   res.copy (parentMap);
   SxVector<Int>::Iterator it;
   ConstPtr ptr = parent;
   while (ptr)  {
      SX_CHECK (ptr.getPtr () != this);
      if (ptr == in) return res;
      for (it = res.begin (); it != res.end (); ++it)
         *it = ptr->parentMap (*it);
      ptr = ptr->parent;
   }
   SX_EXIT;
   return res;
}

