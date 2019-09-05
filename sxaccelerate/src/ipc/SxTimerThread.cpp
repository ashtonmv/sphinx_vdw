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

#include <SxTimerThread.h>
#include <SxTime.h>

SxTimerThread::SxTimerThread ()
   : intervalSec(-1.),
     timeoutSec(-1),
     stopCmd(false)
{
   // empty
}

SxTimerThread::~SxTimerThread ()
{
   // empty
}

void SxTimerThread::setInterval (double intervalSec_)
{
   intervalSec = intervalSec_;
}

void SxTimerThread::setTimeout (double timeoutSec_)
{
   timeoutSec = timeoutSec_;
}

void SxTimerThread::start (Priority priority)
{
   SX_MUTEX (mutex)  {
      stopCmd = false;
   }
   SxSystemThread::start (priority);
}

void SxTimerThread::stop ()
{
   SX_MUTEX (mutex)  {
      stopCmd = true;
      condition.wakeAll ();
   }
}

void SxTimerThread::main ()
{
   double timeStart = SxTime::getRealTime ();
   double t1 = timeStart;
   double timeWait = intervalSec;
   if (intervalSec < 0. || (timeoutSec > 0. && intervalSec > timeoutSec))  {
      timeWait = timeoutSec;
   }
   
   mutex.lock ();
      while (!stopCmd)  {
         // --- do not wait longer than timeout
         if (timeoutSec > 0. && t1 + timeWait > timeStart + timeoutSec)  {
            timeWait = timeStart + timeoutSec - t1;
         }
         // --- wait on condition (stop) for max time timeWait seconds
         //     timeWait is either interval or remaining time till the timeout
         condition.wait (&mutex, t1 + timeWait);
         if (stopCmd)  {
            break;
         }
         
         double t2 = SxTime::getRealTime ();
         // --- interval
         if (intervalSec > 0. && t2 + 1e-3 > t1 + intervalSec)  {
            mutex.unlock ();
               interval ();
               signalInterval.send (t2);
            mutex.lock ();
            t1 = t2;
         }
         // --- timeout
         if (timeoutSec > 0. && t2 + 1e-3 > timeStart + timeoutSec)  {
            mutex.unlock ();
               timeout ();
               signalTimeout.send (t2);
            mutex.lock ();
            break;
         }
      }
   mutex.unlock ();
}

