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

#include <SxProcUnit.h>
#include <SxProcUnitGroup.h>
#include <SxThreadPool.h>

SxProcUnit::SxProcUnit ()
   : unitId (0)
{
   // empty
}

SxProcUnit::SxProcUnit (const SxPtr<SxProcUnitGroup> &group_)
   : unitId (0)
{
   if (group_)  {
      group = group_;
      group->registerProcUnit (getThis());
   }
}

void SxProcUnit::setThreadPool (const SxPtr<SxThreadPool> &pool_)
{
   SX_CHECK (pool_.getPtr ());
   pool = pool_;
}

SxProcUnit::~SxProcUnit ()
{
   if (group)  group->deregisterProcUnit (getThis());
}

//void SxProcUnit::setId (ssize_t id_)
//{
//   unitId = id_;
//}


//ssize_t SxProcUnit::getId () const
//{
//   return unitId;
//}


void SxProcUnit::reparent (const SxPtr<SxProcUnitGroup> &newGroup)
{
   group = newGroup;
   if (group)  group->registerProcUnit (getThis());
}
