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

#ifndef _SX_MUTEX_H_
#define _SX_MUTEX_H_

#include <SxIPC.h>
#include <SxConfig.h>

#ifdef WIN32
#  include <windows.h>
#else
#  include <pthread.h>
#endif /* WIN32 */

/** \brief Mutex - Mutual exclusion lock

    Exclusive access by a thread to a variable or set of variables.

    \b SxMutex = S/PHI/nX Mutex

    SxMutex uses critical section on windows. It Seems that pthread
    mutex in pthread part of SxMutex can do both critical section and mutex.
    Critical sections are faster (25x from one source) but they can not be
    shared by multiple processes as locking mechanism.

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_IPC SxMutex
{
   public:
      SxMutex ();
     ~SxMutex ();

      void lock ();
     
      /** Tries to lock the mutex.
          \return true if the lock was obtained by this call or false if the
                  mutex is already held by someone else */
      bool trylock ();
     
      void unlock ();

   protected:
      friend class SxThreadCond;
#     ifdef WIN32
         //HANDLE mutex; slower, but inter-process
         CRITICAL_SECTION mutex;
#     else
         pthread_mutex_t mutex;
#     endif
};

class SX_EXPORT_IPC SxMutexBlock
{
   public:
      SxMutex *mutex;
      SxMutexBlock (SxMutex *m) : mutex(m) { mutex->lock();   }
     ~SxMutexBlock ()                      { mutex->unlock(); }
};

#define SX_SERIALIZE      SX_BLOCK(SxMutex)
#define SX_MUTEX(mutex)   SX_BLOCK_ARG(SxMutexBlock,&mutex)

#endif /* _SX_MUTEX_H_ */
