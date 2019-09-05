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

#include <SxMutex.h>
#include <SxError.h>
#include <errno.h>


SxMutex::SxMutex ()
{
#  ifdef WIN32
      InitializeCriticalSection (&mutex);
      //mutex = CreateMutex (NULL, FALSE, NULL);
      //if (!mutex)  {
      //   SX_CHECK_WIN_ERR (GetLastError ());
      //}
#  else
      int err = pthread_mutex_init (&mutex, NULL);
      SX_CHECK (err == 0);
#  endif
}


SxMutex::~SxMutex ()
{
#  ifdef WIN32
      DeleteCriticalSection (&mutex);
      //if (!CloseHandle (mutex))  {
      //   SX_CHECK_WIN_ERR (GetLastError ());
      //}
#  else
      int err = pthread_mutex_destroy (&mutex);
      SX_CHECK (err == 0);
#  endif
}


void SxMutex::lock ()
{
#  ifdef WIN32
      EnterCriticalSection (&mutex);
      //DWORD err = WaitForSingleObject (mutex, INFINITE);
      //if (err == WAIT_FAILED)  {
      //   SX_CHECK_WIN_ERR (GetLastError ());
      //}  else if (err == WAIT_ABANDONED)  {
      //   cout << "ERROR: WaitForSingleObject: mutex was not released by the "
      //           "thread that owned the mutex before that thread terminated"
      //        << endl;
      //   SX_EXIT;
      //}
#  else
      int err = pthread_mutex_lock (&mutex);
      SX_CHECK (err == 0, err);
#  endif
}


bool SxMutex::trylock ()
{
#  ifdef WIN32
      if (TryEnterCriticalSection (&mutex))  {
         return true;
      }
      return false;
      //DWORD err = WaitForSingleObject (mutex, 0);
      //if (err == WAIT_FAILED || err == WAIT_TIMEOUT)  {
      //   // --- It was not possible to lock the mutex
      //   return false;
      //}  else if (err == WAIT_ABANDONED)  {
      //   cout << "ERROR: WaitForSingleObject: mutex was not released by the "
      //           "thread that owned the mutex before that thread terminated"
      //        << endl;
      //   SX_EXIT;
      //}
      //// --- The mutex is locked by this call
      //return true;
#  else
      int err = pthread_mutex_trylock (&mutex);
      if (err == 0)  {
         // --- The mutex is locked by this call
         return true;
      }  else if (err == EBUSY)  {
         // --- The mutex was already locked by someone else.
         return false;
      }  else  {
         SX_CHECK (err == 0);
      }

      // --- It was not possible to lock the mutex.
      return false;
#  endif
}


void SxMutex::unlock ()
{
#  ifdef WIN32
      LeaveCriticalSection (&mutex);
      //if (!ReleaseMutex (mutex))  {
      //   SX_CHECK_WIN_ERR (GetLastError ());
      //}
#  else
      int err = pthread_mutex_unlock (&mutex);
      SX_CHECK (err == 0);
#  endif
}
