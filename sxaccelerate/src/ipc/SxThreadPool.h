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

#ifndef _SX_THREAD_POOL_H_
#define _SX_THREAD_POOL_H_

#include <SxArray.h>
#include <SxPtr.h>
#include <SxThreadCond.h>
#include <SxThread.h>
#include <SxSystemThread.h>
#include <SxIPC.h>

/** \brief Thread Pool

    \b SxThreadPool = S/PHI/nX Thread Pool

    ....

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_IPC SxThreadPool : public SxThis<SxThreadPool>
{
   public:
      SxThreadPool (int nThreads=0);
      virtual ~SxThreadPool ();

      /** \return the number of resources (System threads) */
      ssize_t getSize () const;

      // --- task manager
      void submit (const SxPtr<SxThread> &); // #exception
      void operator<< (const SxPtr<SxThread> &); // #exception
      bool contains (const SxPtr<SxThread> &thread_) const;
      bool isRunning (const SxPtr<SxThread> &thread_) const;
      void wait (const SxPtr<SxThread> &thread_);
      void waitAll ();
      void stop(); // stop all pool threads
      ssize_t getNTasks () const;

      static int getNProcs ();

   protected:
      bool containsWaiting (const SxPtr<SxThread> &thread_) const;
      bool containsRunning (const SxPtr<SxThread> &thread_) const;
      void resizePool (ssize_t nThreads);
      void shutdown ();

      class PoolThread : public SxSystemThread
      {
         public:
            SxThreadPool  *pool;
            int            id;

            PoolThread (SxThreadPool *pool, int id);
            virtual ~PoolThread ();
            virtual void main ();
      };

      SxArray<SxPtr<PoolThread> >             pool;
      SxList<SxPtr<SxThread> >                tasks;
      SxList<SxList<SxPtr<SxThread> >::Node*> waitingTasks;
      SxList<SxList<SxPtr<SxThread> >::Node*> runningTasks;
      mutable SxMutex                         mutex;
      mutable SxThreadCond                    condition;
      mutable SxThreadCond                    waitCondition;
      bool                                    terminateCmd;
};

#endif /* _SX_THREAD_POOL_H_ */
