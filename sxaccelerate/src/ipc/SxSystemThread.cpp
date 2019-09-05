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

#include <SxSystemThread.h>
#include <SxString.h>
#include <SxError.h>
#include <SxTime.h>
#include <unistd.h>
#ifdef WIN32
#  include <process.h>
#endif

SxSystemThread::SxSystemThread ()
  : running (false),
    hasLastException (false),
    hasThreadId (false)
{
   SX_TRACE ();
   // empty
}

SxSystemThread::~SxSystemThread ()
{
   SX_TRACE ();
   if (isRunning ())  SX_EXIT;
}

void SxSystemThread::start (Priority priority)
{
   SX_TRACE ();
   SX_CHECK (!hasThreadId);

   try {
#  ifdef WIN32
      // --- prepare user data
      userData = (void **)HeapAlloc(GetProcessHeap(), HEAP_ZERO_MEMORY, 
                          sizeof(void *));
      *userData = this;


      // --- create thread
      // Note: Under MingW we have to use _beginthreadex() instead of 
      //       CreateThread() to avoid memleaks!
      threadHandle = (HANDLE)_beginthreadex (
                        NULL, 0, &loopCallback, userData,
                        0, &threadId.id); 
#  else
      int newPriority=0;

      // --- setup thread priority
#     ifndef SX_ANDROID
      if (priority == SameAsParent)  {
         pthread_attr_setinheritsched (&threadAttribs, PTHREAD_INHERIT_SCHED);
      } else  {
         // --- get default schedule policy
         int schedPolicy;
         int error = pthread_attr_getschedpolicy (&threadAttribs, &schedPolicy);
         SX_CHECK (!error);
         schedPolicy = SCHED_RR;
         pthread_attr_setschedpolicy (&threadAttribs, SCHED_RR);
   
         // --- get allowed priority range
         int minPriority = sched_get_priority_min (schedPolicy);
         int maxPriority = sched_get_priority_max (schedPolicy);
         SX_CHECK (minPriority > 0);
         SX_CHECK (maxPriority > 0);

         if      (priority == Idle)     newPriority = minPriority;
         else if (priority == Highest)  newPriority = maxPriority;
         else {
            newPriority = ((maxPriority - minPriority) / maxPriority) * priority
                        +   minPriority;
         
         }

         sched_param schedParam;
         schedParam.sched_priority = newPriority;

         pthread_attr_setinheritsched (&threadAttribs, PTHREAD_EXPLICIT_SCHED);
         pthread_attr_setschedparam (&threadAttribs, &schedParam);
      }
#     endif /* SX_ANDROID */
   

      int err = pthread_create (&threadId.id, NULL,
                                (void *(*)(void *))loopCallback,
                                this);
      if (err)  {
         // TODO: error
         return;
      }
#  endif
   } catch (...) {
      SX_EXIT; // main routine of derived thread class needs try-catch block
   }
   hasThreadId = true;
}

void SxSystemThread::launcher ()
{
   SX_TRACE ();

   SX_MUTEX (mutex) {
      running = true;
      hasLastException = false;
   }

   try {
      main ();
   }  catch (SxException e) {
      SX_MUTEX (mutex) {
         hasLastException = true;
         lastException = e;
      }
   }  catch (...) {
#  ifdef WIN32
      DWORD err = GetLastError ();
      SxString msg = "System thread interrupted unexpectedly: " 
         + sxstrerror (err) 
         + "(" + sxprintf ("0x%x", err) + ")";
#  else
      SxString msg = "System thread interrupted unexpectedly";
#  endif
      SX_MUTEX (mutex) {
         hasLastException = true;
         lastException = SxException(msg.elements, __FILE__, __LINE__);
      }
   }
   SX_MUTEX (mutex) {
      running = false;
   }
}

bool SxSystemThread::hasException () const
{
   SX_TRACE ();
   bool res = false;
   SX_MUTEX (mutex) {
      res = hasLastException;
   }
   return res;
}

SxException SxSystemThread::getException () const
{
   SX_TRACE ();
   SxException e;
   SX_MUTEX (mutex) {
      e = lastException;
   }
   return e;
}

#ifdef WIN32
unsigned __stdcall SxSystemThread::loopCallback (void *userData)
#else
void SxSystemThread::loopCallback (void *userData)
#endif
{
   SX_TRACE ();
#  ifdef WIN32
      SxSystemThread *obj = *static_cast<SxSystemThread **>(userData);
      SX_CHECK (obj);
      obj->launcher ();
      _endthreadex (0);
      return 0;
#  else
      SxSystemThread *obj = static_cast<SxSystemThread *>(userData);
      SX_CHECK (obj);
      obj->launcher ();
      pthread_exit(0);
#  endif
}


void SxSystemThread::wait ()
{
   SX_TRACE ();
   if (hasThreadId)  {
#     ifdef WIN32
         if (threadHandle != NULL)  {
            DWORD result = WaitForSingleObject (threadHandle, INFINITE);
            if (result != WAIT_OBJECT_0)   {
               cout << "Wait for single object unexpectedly quit: " << result << endl;
               if (result == WAIT_FAILED)  {
                  DWORD err = GetLastError ();
                  cout << "Wait for tid " + SxString(threadId.id)
                     + " failed: " + sxstrerror (err);
               }
            }
            SX_MUTEX (waitMutex)  {
               if (userData) HeapFree (GetProcessHeap (), 0, userData);
               userData = NULL;
               if (threadHandle) CloseHandle (threadHandle);
               threadHandle = NULL;
            } 
         }
#     else
         pthread_join (threadId.id, NULL);
#     endif
      SX_MUTEX (waitMutex) {
         hasThreadId = false;
         threadId = SxThreadId ();
      }
   }
}


// void SxSystemThread::terminate ()
// {
//    running = false;
// #  ifdef WIN32
//       _endthreadex (0);
// #  else
//       pthread_exit (0);
// #  endif
// }


bool SxSystemThread::isRunning () const
{
   SX_TRACE ();
   bool res = false;
   SX_MUTEX (mutex) {
      res = running;
   }
   return res;
}

void SxSystemThread::main ()
{
   SX_TRACE ();
   SX_CHECK (mainLambda);
   mainLambda ();
}
