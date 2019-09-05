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

#ifndef _SX_PROC_UNIT_H_
#define _SX_PROC_UNIT_H_

#include <SxIPC.h>
#include <SxPtr.h>
#include <SxSignals.h>

class SxProcUnitGroup;
class SxThreadPool;

/** \brief ...

    \b SxProcUnit = S/PHI/nX Process Unit

    ....

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_IPC SxProcUnit : public SxThis<SxProcUnit>
{
   public:

      SxProcUnit ();
      SxProcUnit (const SxPtr<SxProcUnitGroup> &);
      virtual ~SxProcUnit ();

      void setThreadPool (const SxPtr<SxThreadPool> &pool_);

      virtual void start ()=0;
      virtual void wait ()=0;

      //void setId (ssize_t id);
      //ssize_t getId () const;

   signals:

      SxSignal<SxProcUnit *, const char *> SX_SIGNAL (sigStarted);
      SxSignal<SxProcUnit *, const char *> SX_SIGNAL (sigFinished);

   protected:

      friend class SxTask;
      friend class SxProcUnitGroup;

      ssize_t                 unitId;
      SxPtr<SxProcUnitGroup>  group;
      SxPtr<SxThreadPool>     pool; // pool container

      void reparent (const SxPtr<SxProcUnitGroup> &);

};

#endif /* _SX_PROC_UNIT_H_ */
