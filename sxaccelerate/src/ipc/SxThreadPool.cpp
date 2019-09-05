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

#include <SxThreadPool.h>
#include <SxThread.h>
#include <SxConfig.h>
#include <SxException.h>
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#ifdef MACOSX
#  include <sys/sysctl.h>
#endif


SxThreadPool::PoolThread::PoolThread (SxThreadPool *pool_, int id_)
   : SxSystemThread (),
     pool(pool_),
     id(id_)
{
   SX_TRACE ();
   SX_DBG_MSG (SxString("PoolThread #" + SxString(id)).ascii ());
}

SxThreadPool::PoolThread::~PoolThread ()
{
   SX_TRACE ();
   SX_DBG_MSG (SxString("PoolThread #" + SxString(id)).ascii ());
}

void SxThreadPool::PoolThread::main ()
{
   SX_TRACE ();
   SX_DBG_MSG (SxString("PoolThread #" + SxString(id)).ascii ());
   bool workingLoop = true;

   pool->mutex.lock ();

   while (workingLoop)  {
      while (!pool->terminateCmd && pool->waitingTasks.getSize () < 1)  {
         SX_DBG_MSG ("PoolThread #" << id << " waiting for new data");
         pool->condition.wait (&pool->mutex);
      }
      SX_DBG_MSG ("PoolThread #" << id << " new data present");

      if (pool->terminateCmd)  {
         SX_DBG_MSG ("PoolThread #" << id << " terminate command");
         workingLoop = false;

      }  else  {
         SX_DBG_MSG ("PoolThread #" << id << " start working");
         // --- move from waiting tasks to running tasks
         SxList<SxPtr<SxThread> >::Node* task = pool->waitingTasks.first ();
         pool->waitingTasks.removeFirst ();
         pool->runningTasks << task;

         // --- run the thread
         SxPtr<SxThread> thread = task->elem;
         SX_CHECK (thread.getPtr ());
         pool->mutex.unlock ();
            thread->mainWrap ();
         pool->mutex.lock ();
         SX_DBG_MSG ("PoolThread #" << id << " work finished");

         // --- remove the task
         pool->runningTasks.removeElement (task);
         pool->tasks.removeItem (task);

         pool->waitCondition.wakeAll ();
      }
   }

   SX_DBG_MSG ("PoolThread #" << id << " terminated");
   pool->mutex.unlock ();

}

// default for SxThreadMgr was : int offset=0, int multiplier=1
SxThreadPool::SxThreadPool (int nThreads_)
   : SxThis<SxThreadPool>(),
     terminateCmd (false)
{
   SX_TRACE ();

   if (nThreads_ < 1)  {
      nThreads_ = SxThreadPool::getNProcs ();
   }

   SX_DBG_MSG ("creating thread pool nThreads=" << nThreads_);

   resizePool (nThreads_);
}

SxThreadPool::~SxThreadPool ()
{
   SX_TRACE ();

   shutdown ();

   ssize_t i;
   for (i=0; i < pool.getSize(); ++i)  {
      SX_CHECK (pool(i).getPtr ());
      pool(i)->wait ();
   }
   SX_DBG_MSG ("thread pool terminated");
}

// --- Expandable thread pools
//bool isExpandable () const { return expandable; }
//void resize (ssize_t nThreads);
//void SxThreadPool::resize (ssize_t nThreads_)
//{
//   SX_CHECK (expandable);
//   resizePool (nThreads_);
//}

ssize_t SxThreadPool::getSize () const
{
   SX_TRACE ();
   return pool.getSize();
}

void SxThreadPool::submit (const SxPtr<SxThread> &thread_)
{
   SX_TRACE ();

   if (getSize () < 1)  {
      SX_THROW ("Can't submit tasks to empty threadpool.");
   }
   if (!thread_.getPtr ())  {
      SX_THROW ("Can't submit empty process unit.");
   }

   SX_MUTEX (mutex)  {
      tasks << thread_;
      waitingTasks << tasks.lastElement;
      condition.wakeOne ();
   }
   // connect the threadpool object to thread
   thread_->setThreadPool (getThis());
}

void SxThreadPool::operator<< (const SxPtr<SxThread> &thread_)
{
   SX_TRACE ();
   submit (thread_);
}

bool SxThreadPool::contains (const SxPtr<SxThread> &thread_) const
{
   SX_TRACE ();
   bool result = false;

   SX_MUTEX (mutex)  {
      result = containsWaiting (thread_) || containsRunning (thread_);
   }
   return result;
}

bool SxThreadPool::isRunning (const SxPtr<SxThread> &thread_) const
{
   SX_TRACE ();
   bool result = false;

   SX_MUTEX (mutex)  {
      result = containsWaiting (thread_) || containsRunning (thread_);
   }

   return result;
}

void SxThreadPool::wait (const SxPtr<SxThread> &thread_)
{
   SX_TRACE ();

   SX_MUTEX (mutex)  {
      while (containsWaiting (thread_) || containsRunning (thread_))  {
         SX_DBG_MSG ("Task is waiting for result");
         waitCondition.wait (&mutex);
      }
   }
}

void SxThreadPool::waitAll ()
{
   SX_TRACE ();

   SX_MUTEX (mutex) {
      while ((waitingTasks.getSize () + runningTasks.getSize ()) > 0) {
         waitCondition.wait (&mutex);
      }
   }
}

ssize_t SxThreadPool::getNTasks () const
{
   SX_TRACE ();

   ssize_t result = 0;

   SX_MUTEX (mutex)  {
      result = waitingTasks.getSize () + runningTasks.getSize ();
   }

   return result;
}

int SxThreadPool::getNProcs ()
{
   SX_TRACE ();

   int nProcs = 0;

   // --- try SX_THREADS first
   char *str = getenv ("SX_THREADS");
   if (str) {
      nProcs = atoi (str);
      if (nProcs >= 1)  {
         SX_DBG_MSG ("SX_THREADS=" << nProcs);
         return nProcs;
      }
   }

   // --- use OS calls to find the number of cores
#  ifdef WIN32
      SYSTEM_INFO info;
      GetSystemInfo (&info);
      return info.dwNumberOfProcessors;
#  else
#    ifdef _SC_NPROCESSORS_ONLN
      // or _SC_NPROCESSORS_CONF
      nProcs = static_cast<int>(sysconf (_SC_NPROCESSORS_ONLN));
      SX_DBG_MSG ("_SC_NPROCESSORS_ONLN=" << nProcs);
      return nProcs;
#    elif defined(CTL_HW) && defined(HW_AVAILCPU) && defined(HW_NCPU)
      size_t len = sizeof(nProcs); 
      int name[2];
      name[0] = CTL_HW;
      name[1] = HW_AVAILCPU;
      sysctl (name, 2, &nProcs, &len, NULL, 0);

      if (nProcs < 1)  {
         name[1] = HW_NCPU;
         sysctl (name, 2, &nProcs, &len, NULL, 0);

         if (nProcs < 1)  {
            nProcs = 1;
         }
      }
      SX_DBG_MSG ("sysctl(CPU)=" << nProcs);
      return nProcs;
#    else
      return 1;
#    endif
#  endif
}

// lock mutex before calling
bool SxThreadPool::containsWaiting (const SxPtr<SxThread> &thread_) const
{
   SX_TRACE ();

   SxList<SxList<SxPtr<SxThread> >::Node*>::ConstIterator it;
   for (it = waitingTasks.begin(); it != waitingTasks.end(); ++it)  {
      SxPtr<SxThread> thread = (*it)->elem;
      if (thread == thread_)  {
         return true;
      }
   }

   return false;
}

// lock mutex before calling
bool SxThreadPool::containsRunning (const SxPtr<SxThread> &thread_) const
{
   SX_TRACE ();

   SxList<SxList<SxPtr<SxThread> >::Node*>::ConstIterator it;
   for (it = runningTasks.begin(); it != runningTasks.end(); ++it)  {
      SxPtr<SxThread> thread = (*it)->elem;
      if (thread == thread_)  {
         return true;
      }
   }

   return false;
}

void SxThreadPool::resizePool (ssize_t nThreads_)
{
   SX_TRACE ();

   SX_DBG_MSG ("resizing thread pool from " << pool.getSize ()
                << " to " << nThreads_ << " threads");

   ssize_t nOld = pool.getSize ();
   SX_CHECK (nThreads_ >= nOld, nThreads_, nOld);

   if (nThreads_ == nOld)  {
      return;
   }

   pool.resize (nThreads_, true);

   for (ssize_t i = nOld; i < nThreads_; ++i)  {
      pool(i) = SxPtr<PoolThread>::create (this, static_cast<int>(i));
      pool(i)->start ();
   }
}

void SxThreadPool::shutdown ()
{
   SX_TRACE ();

   SX_MUTEX (mutex)  {
      terminateCmd = true;
      condition.wakeAll ();
   }
}

void SxThreadPool::stop ()
{
   SX_TRACE ();

   SX_MUTEX(mutex) {
      terminateCmd = true;
      condition.wakeAll();
   }
   ssize_t i;
   for (i = 0; i < pool.getSize(); ++i) {
      SX_CHECK(pool(i).getPtr());
      pool(i)->wait();
   }
}

#if 0
// --- handle priority for tasks

SxArray<SxList<SxTask> >  tasks;  // :priority

void SxThreadPool::submit (const SxTask &task)
{
      for (int i=0; i < 100; i++)  {
         SX_CHECK (tasks(i).getSize() == 0, tasks(i).getSize());
      }
      tasks(task.getPriority()) << task;
      nTasks++;
}

bool SxThreadPool::processNextTask (int poolId, bool useCond)
{
  // SX_MUTEX (qMutex) {
      if (nTasks > 0)  {

         SxList<SxTask>::Iterator it;
         SxArray<int> rank(tasks.getSize());
         for (int p=0; p < SxTask::MaxPriority; ++p)  {
            rank(p) = -1;

            if (tasks(p).getSize() > 0)
               if (tasks(p).first().spoolingExpired())  {
                  tasks(p).removeFirst ();
                  nTasks--;
               }

            if (tasks(p).getSize() > 0)
               rank(p) = tasks(p).first().updatePriority();
         }

         // --- determine task with highest rank
         int maxVal    = -1;
         int maxTaskId = -1;
         for (int p=0; p < SxTask::MaxPriority; ++p)  {
            if (rank(p) > maxVal)  {
               maxVal    = rank(p);
               maxTaskId = p;
            }
         }

         // --- spool job
         if (maxTaskId != -1)  {
            cout << "SPOOLING... " << poolId << endl;
            SxPtr<SxThread> unit = tasks(maxTaskId).first().getTask();
            tasks(maxTaskId).removeFirst ();
            nTasks--;

            pool(poolId)->thread = unit;
            if (useCond)
               SX_THREAD_SET_COND(pool(poolId)->runCond);
            cout << "SPOOLING... " << poolId << " DONE\n";
            return true;
         } else {
            cout << "NOT SPOOLING...\n";
         }
         cout.flush ();
      }
   //}
      return false;
}

#endif /* if 0 */

