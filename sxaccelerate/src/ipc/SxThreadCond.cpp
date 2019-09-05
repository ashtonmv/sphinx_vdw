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

#include <SxThreadCond.h>
#include <SxTime.h>
#include <errno.h>


SxThreadCond::SxThreadCond () 
{
#  ifdef WIN32
      InitializeConditionVariable (&cond);
#  else
      int err = pthread_mutex_init (&mutex, NULL);
      SX_CHECK (err == 0);
      err = pthread_cond_init (&cond, NULL);
      SX_CHECK (err == 0);
#  endif /* WIN32 */
}

SxThreadCond::~SxThreadCond ()
{
#  ifdef WIN32
      // empty
#  else
      int err = pthread_cond_destroy (&cond);
      SX_CHECK (err == 0);
      err = pthread_mutex_destroy (&mutex);
      SX_CHECK (err == 0);
#  endif /* WIN32 */
}


bool SxThreadCond::wait (SxMutex *mutex_, double sec_)
{
   SX_CHECK (mutex_);
#  ifdef WIN32
      DWORD dwmsec = 0;
      if (sec_ < 0.)  {
         dwmsec = INFINITE;
      }  else  {
         // --- abstime to time-out in milliseconds
         double t = sec_ - SxTime::getRealTime ();
         if (t > 0.)  {
            dwmsec = static_cast<DWORD>(t * 1e3);
         }
      }
      if (!SleepConditionVariableCS (&cond, &mutex_->mutex, dwmsec))  {
         DWORD err = GetLastError ();
         if (err == ERROR_TIMEOUT)  {
            return true;
         }
         if (err)  {
            std::cout << "ERROR: " << sxstrerror (err) << std::endl;
            SX_EXIT;
         }
      }
      // --- no timeout
      return false;
#  else
      int err = pthread_mutex_lock (&mutex);
      SX_CHECK (err == 0);
   
      mutex_->unlock ();
      bool result = wait (sec_);
      mutex_->lock ();
   
      return result;
#  endif /* WIN32 */
}

#ifndef WIN32
bool SxThreadCond::wait (double sec_)
{
   bool timeout = false; // true if time limit exceeded
   
      int err = 0;
      
      if (sec_ < 0.)  {
         err = pthread_cond_wait (&cond, &mutex);
      }  else  {
         // --- time interval
//         time_t sec = (int)(sec_);
//         unsigned long msec = (unsigned long)((sec_ - (double)sec) * 1e3);
//
//         struct timeval tp;
//         SX_CHECK_ERR (gettimeofday (&tp, NULL));
//
//         struct timespec deadline;
//         deadline.tv_nsec = tp.tv_usec * 1000L + msec * 1000000L;
//         deadline.tv_sec  = tp.tv_sec + sec + (deadline.tv_nsec / 1000000000L);
//         deadline.tv_nsec %= 1000000000L;

         // --- absolute time
         struct timespec deadline;
         deadline.tv_sec  = (time_t)(sec_);
         deadline.tv_nsec = (long)((sec_ - (double)deadline.tv_sec) * 1e9);
            
         err = pthread_cond_timedwait (&cond, &mutex, &deadline);
      }

      if (err == ETIMEDOUT)  {
         timeout = true;
      }  else if (err == EINVAL)  {
         SX_CHECK (err == 0);
      }
      // --- do not EXIT on undocumented errors (err > 0)
      
      err = pthread_mutex_unlock (&mutex);
      SX_CHECK (err == 0);

   return timeout;
}
#endif /* UNIX */

void SxThreadCond::wakeOne ()
{
#  ifdef WIN32
      WakeConditionVariable (&cond);
#  else
      int err = pthread_mutex_lock (&mutex);
      SX_CHECK (err == 0);
      err = pthread_cond_signal (&cond);
      SX_CHECK (err == 0);
      err = pthread_mutex_unlock (&mutex);
      SX_CHECK (err == 0);
#  endif /* WIN32 */
}


void SxThreadCond::wakeAll ()
{
#  ifdef WIN32
      WakeAllConditionVariable (&cond);
#  else
      int err = pthread_mutex_lock (&mutex);
      SX_CHECK (err == 0);
      err = pthread_cond_broadcast (&cond);
      SX_CHECK (err == 0);
      err = pthread_mutex_unlock (&mutex);
      SX_CHECK (err == 0);
#  endif /* WIN32 */      
}
