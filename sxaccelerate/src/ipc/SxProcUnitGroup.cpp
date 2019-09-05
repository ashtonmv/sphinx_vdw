// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#include <SxProcUnitGroup.h>
#include <SxProcUnit.h>
#include <SxString.h>

SxProcUnitGroup::SxProcUnitGroup ()
   : groupId (SxProcUnitGroup::getNextGroupId())
{
   // empty
}

SxProcUnitGroup::SxProcUnitGroup (const SxPtr<SxProcUnitGroup> &parent_) 
   : parent (parent_),
     groupId (SxProcUnitGroup::getNextGroupId())
{
   // empty
}

SxProcUnitGroup::~SxProcUnitGroup ()
{
   if (groupId == 0)  return;

// // --- remove nested groups
// SxList<SxPtr<SxProcUnitGroup> >::Iterator gIt;
// for (gIt = groups.begin(); gIt != groups.end(); ++gIt)  {
//    removeGroup (*gIt);
// }
// SX_CHECK (groups.getSize() == 0);

   // --- reparent process units of dying group
   SxList<SxPtr<SxProcUnit> >::Iterator tIt;
   for (tIt = procUnits.begin (); tIt != procUnits.end(); ++tIt)  {
      (*tIt)->reparent (parent);
   }

   // --- deregister from parent
// if (parent)  {
//    parent->groups.removeElement (getThis());
// }

// if (groupId == 0 && procUnit.getSize() == 0)  {
//    sxprintf ("   shutting down procUnit mgr...\n");
//    
// }
}


ssize_t SxProcUnitGroup::getNextGroupId ()
{
   static ssize_t groupId = -1;
   groupId++;
   return groupId;
}

// --------------------------------------------------------------------


void SxProcUnitGroup::removeGroup (const SxPtr<SxProcUnitGroup> &grp)
{
   SX_CHECK (grp);
   SX_CHECK (groups.contains (grp), groupId, grp->groupId);

   // --- remove all nested groups
   SxList<SxPtr<SxProcUnitGroup> >::Iterator grpIt;
   for (grpIt = grp->groups.begin(); grpIt != grp->groups.end(); ++grpIt)
      grp->removeGroup (*grpIt);
   SX_CHECK (grp->groups.getSize() == 0);

   // --- reparent process units of dying group
   SxList<SxPtr<SxProcUnit> >::Iterator tIt;
   for (tIt = grp->procUnits.begin (); tIt != grp->procUnits.end(); ++tIt)  {
      grp->procUnits.removeElement (*tIt);
      (*tIt)->reparent (getThis());
   }
   SX_CHECK (grp->procUnits.getSize() == 0);

   grp->parent = SxPtr<SxProcUnitGroup> ();
   groups.removeElement (grp);
}


void SxProcUnitGroup::registerProcUnit (const SxPtr<SxProcUnit> &u)
{
   SX_CHECK (!procUnits.contains (u));
   SX_CHECK (!getTopLevelGroup()->getProcUnits(true).contains(u));
   procUnits << u;
}


void SxProcUnitGroup::deregisterProcUnit (const SxPtr<SxProcUnit> &u)
{
   SX_CHECK (procUnits.contains (u));
   procUnits.removeElement (u);
}


SxList<SxPtr<SxProcUnit> > SxProcUnitGroup::getProcUnits (bool recursive) const
{
   SxList<SxPtr<SxProcUnit> > res = procUnits;
   if (recursive)  {
      SxList<SxPtr<SxProcUnitGroup> >::ConstIterator grp;
      for (grp = groups.begin(); grp != groups.end(); ++grp)  {
         res << (*grp)->getProcUnits (recursive);
      }
   }
   return res;
}



void SxProcUnitGroup::print (int indent) const
{
   SxString indentStr;
   for (int i=0; i < indent; ++i) indentStr += " ";
   cout << indentStr << "ProcUnitGroup: " << groupId << endl;
   cout << indentStr << "   proc units:  " << procUnits.getSize() << endl;
   SxList<SxPtr<SxProcUnitGroup> >::ConstIterator it;
   for (it = groups.begin(); it != groups.end(); ++it)
      (*it)->print (indent+3);
}
