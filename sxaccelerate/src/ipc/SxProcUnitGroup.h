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

#ifndef _SX_PROC_UNIT_GROUP_H_
#define _SX_PROC_UNIT_GROUP_H_

#include <SxIPC.h>
#include <SxProcUnit.h>
#include <SxPtr.h>
#include <SxList.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_IPC SxProcUnitGroup : public SxThis<SxProcUnitGroup>
{
   public:
      SxProcUnitGroup ();
      SxProcUnitGroup (const SxPtr<SxProcUnitGroup> &parent_);
      virtual ~SxProcUnitGroup ();


      SxList<SxPtr<SxProcUnit> > getProcUnits (bool recursive=false) const;

      virtual void start ()=0;
      virtual void barrier ()=0;

      virtual SxPtr<SxProcUnitGroup> addGroup (int maxItems=0)=0;
      void removeGroup (const SxPtr<SxProcUnitGroup> &);

      void print (int indent=0) const;

   protected:
      friend class SxProcUnit;

      SxPtr<SxProcUnitGroup>          parent;
      SxList<SxPtr<SxProcUnitGroup> > groups;
      SxList<SxPtr<SxProcUnit> >      procUnits;

      ssize_t groupId;

      SxPtr<SxProcUnitGroup> getGroup ();

      virtual void registerProcUnit   (const SxPtr<SxProcUnit> &);
      virtual void deregisterProcUnit (const SxPtr<SxProcUnit> &);

      virtual SxPtr<SxProcUnitGroup> getTopLevelGroup ()=0;

      static ssize_t getNextGroupId ();
};

#endif /* _SX_PROC_UNIT_GROUP_H_ */
