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

#ifndef _SX_THREAD_COND_H_
#define _SX_THREAD_COND_H_

#include <SxIPC.h>
#include <SxMutex.h>

/** \brief Condition Variable

    Synchronization of threads based upon data value.

    \b SxThreadCond = S/PHI/nX Thread Condition

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_IPC SxThreadCond
{
   public:
      SxThreadCond ();
     ~SxThreadCond ();

      /** Calling thread blocks and waits on the condition.
      
\code
   SxMutex mutex;
   SxThreadCond condition;
   bool hasNewData = false;
   
   SxSystemThread A  {
      // --- wait up to 5 seconds for new data
      bool timeout = false;
      double t = SxTime::getRealTime ();
      mutex.lock ();
         while (!timeout && !hasNewData)  {
            timeout = condition.wait (&mutex, t + 5.);
         }   
      mutex.unlock ();
   }
   
   SxSystemThread B  {
      // --- signal new data to thread A after 1 second
      SxTime::msleep (1000);
      mutex.lock ();
         hasNewData = true;
         condition.wakeAll ();
      mutex.unlock ();
   }
\endcode
           \return true if time limit exceeded */
      bool wait (SxMutex *mutex, double absoluteTimeSec=-1.);

      /** Unblocks one thread waiting on the condition (if there is any). */
      void wakeOne ();
      
      /** Unblocks all threads waiting on the condition. */
      void wakeAll ();
      

   protected:
#     ifdef WIN32
         CONDITION_VARIABLE cond;
#     else
         pthread_cond_t cond;
         pthread_mutex_t mutex;
         bool wait (double sec);
#     endif /* WIN32 */
};

#if 0
/** 
  author Sixten Boeck, boeck@mpie.de */
class SxThreadCondWrapper
{
   public:
      SxThreadCond *cond;
      bool emit;

      SxThreadCondWrapper (SxThreadCond &cond_, bool emit_) 
         : cond(&cond_), emit(emit_)
      {
         cond->lock ();
         if (!emit) cond->wait ();
         cond->unlock ();
      }

      ~SxThreadCondWrapper ()
      {
         cond->lock ();
         if (emit)  cond->release ();
         cond->unlock ();
      }
};

#define SX_UNIQ_THR_COND_NAME(x)   sxCondWrapper ## x
#define SX_UNIQ_THR_COND(x)        SxThreadCondWrapper SX_UNIQ_THR_COND_NAME(x)

#define SX_THREAD_SET_COND(cond)   { SX_UNIQ_THR_COND(__LINE__) (cond, true); }
#define SX_THREAD_WAIT_COND(cond)  { SX_UNIQ_THR_COND(__LINE__) (cond, false); }
#endif 

#endif /* _SX_THREAD_COND_H_ */
