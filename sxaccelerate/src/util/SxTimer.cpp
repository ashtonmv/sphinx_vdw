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

#include <SxTimer.h>
#include <stdio.h>
#include <SxError.h>
#include <SxTime.h>
#include <iostream>

#ifndef WIN32
#  include <sys/time.h>
#  include <sys/resource.h>
#else
#  include <windows.h>
#endif
#include <stdio.h>

SX_EXPORT_UTIL SxMap<SxString, int> registeredTimerEnums;

SxTimer::SxTimer ()
{
   // empty
}

SxTimer::SxTimer (int nTimers_, bool realTime_)
{
   char *str = getenv ("SX_WALL_TIMER");
   if (str)  {
      SxUtil::getGlobalObj (); // note: make sure it is initialized
      SxString mode(str);
      if (mode.getSize () == 0 || mode == "on")  {
         //if (!realTime_)  {
         //   std::cout << "Overriding timing mode: using wall-clock time."
         //             << endl;
         //}
         realTime_=true;
      } else if (mode == "off")  {
         //if (realTime_)  {
         //   std::cout << "Overriding timing mode: using total CPU time."
         //             << endl;
         //}
         realTime_=false;
      } else {
         std::cout << "Unknown timing mode '" << mode 
                   << "' requested in environment variable SX_WALL_TIMER.\n"
                   << "Must be 'on' or 'off'." << endl;
         SX_EXIT;
      }
   } else if (SxUtil::getGlobalObj ().nProcs > 1) {
      //if (!realTime_)  {
      //   std::cout << "Overriding timing mode for multi-threaded run: "
      //                "using wall-clock time" << endl;
      //}
      realTime_ = true;
   }
   init (nTimers_, realTime_);
}


SxTimer::~SxTimer ()
{
   // empty
}


void SxTimer::init (int nTimers_, bool realTime_)
{
   SX_CHECK (accTime.getSize() == 0);

   nTimers  = nTimers_;
   realTime = realTime_;
   accTime.resize (nTimers);
   procTime.resize (nTimers);
   nCalls.resize (nTimers);
   isRunning.resize (nTimers);
   timerName.resize (nTimers);
   for (int i=0; i < nTimers; i++)   reset (i);
}

void SxTimer::resize (int newSize)
{
   SX_CHECK (newSize > 0);
   if (accTime.getSize () == 0)  {
      init (newSize);
      return;
   }
   SX_CHECK (newSize > nTimers, newSize, nTimers);

   // --- resize and keep previous data
   accTime.resize   (newSize, true);
   procTime.resize  (newSize, true);
   nCalls.resize    (newSize, true);
   isRunning.resize (newSize, true);
   timerName.resize (newSize, true);

   // reset new timers
   int i = nTimers;
   nTimers = newSize;
   for (; i < newSize; i++) reset (i);
}

void SxTimer::start (int iTimer)
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);
// SX_CHECK      (!isRunning(iTimer));

   procTime(iTimer)  = getClock ();
   isRunning(iTimer) = true;
   nCalls(iTimer)++;
}

void SxTimer::stop (int iTimer)
{
   SX_CHECK (isRunning(iTimer));
   refresh ();
   isRunning(iTimer) = false;
}

void SxTimer::refresh ()
{
   double now = getClock ();

   double timeDiff;
   for (int i=0; i < nTimers; i++)  {
      timeDiff = now - procTime(i);
      procTime(i) = now;
      if (isRunning(i) && timeDiff > 0.) accTime(i) += timeDiff;
   }
}

void SxTimer::reset (int iTimer)
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   accTime(iTimer)  = 0.;
   procTime(iTimer)  = getClock ();
   isRunning(iTimer) = false;
   nCalls(iTimer)    = 0;
}

bool SxTimer::isIdle (int iTimer) const
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   return !isRunning(iTimer);
}


double SxTimer::getTime (int iTimer) const
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   return accTime(iTimer);
}


SxString SxTimer::toMinutes (double t)
{
   SX_CHECK (t >= 0, t);

   int m = (int)(t / 60);
   int s = (int)(t - 60 * m);
   SxString minutes (m);
   SxString seconds = s < 10 ? SxString("0") + s : s;

   return minutes + ":" + seconds;
}


void SxTimer::setName (int iTimer, const SxString &name)
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);
   SX_CHECK (name.getSize() <= 20, name.getSize());  // see output format

   timerName(iTimer) = name;
   // update getMaxTimerId for global timer
   if (&getGlobalTimer () == this && getMaxTimerId () <= iTimer)
      getMaxTimerId () = iTimer;
}


double SxTimer::getAvgTime (int iTimer) const
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   return accTime(iTimer) / (double)( !nCalls(iTimer) ? 1 : nCalls(iTimer) ) ;
}


int SxTimer::getNCalls (int iTimer) const
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   return nCalls(iTimer);
}

SxString SxTimer::getName (int iTimer) const
{
   SX_CHECK (iTimer >= 0 && iTimer < nTimers, iTimer, nTimers);

   if (timerName(iTimer) == "")  return "(unknown)";
   else                          return timerName(iTimer);
}



double SxTimer::getClock () const
{
   if (!realTime)  {
      return SxTime::getProcessTime ();
   }  else  {
      return SxTime::getRealTime ();
   }
}


void SxTimer::print (int iTlTimer, bool final) const
{
   cout << endl;
   cout << SX_SEPARATOR;
   cout << "| TIMING RESULTS (" << (realTime ? "wall-clock" 
                                             : "user"      ) << " time)\n";
   cout << SX_SEPARATOR;
   cout << "|  Timer name                   "
        << "   Accumulated Time    Average Time   #Calls\n";
   cout << "|\n";
   SxString str;
   if (final) str="Final "; else str="";
   if (iTlTimer > -1 && getTime(iTlTimer) > 0.)  {
      for (int i=0; i < nTimers; i++)  {
         if (getNCalls(i) > 0)  {
            sxprintf ("|  %-26s  %12.3f %5.1f%%  %14.6f     %d\n",
                    (str+getName(i)).ascii(),
                    getTime(i),
                    100. * getTime(i) / getTime(iTlTimer),
                    getAvgTime(i),
                    getNCalls(i));
         }
      }
   }  else  {
      for (int i=0; i < nTimers; i++)  {
         if (getNCalls(i) > 0)  {
            sxprintf ("|  %-26s  %12.3f         %14.6f      %d\n",
                    (str+getName(i)).ascii(),
                    getTime(i),
                    getAvgTime(i),
                    getNCalls(i));
         }
      }
   }
   cout << SX_SEPARATOR;
#  ifndef NDEBUG
   cout << "| DEBUG MODE: Timings are not meaningful!\n";
   cout << SX_SEPARATOR;
# endif
   cout << endl;
}


void SxTimer::printUsage () const
{
#  ifndef WIN32
      struct rusage u;
      getrusage (RUSAGE_SELF, &u);

      cout << SX_SEPARATOR;
      cout << "| PROGRAM STATISTICS\n";
      cout << SX_SEPARATOR;
      double uTime = double(u.ru_utime.tv_sec)
                   + 0.000001*double(u.ru_utime.tv_usec);
      if (uTime > 0. )
         cout << "| User time:                           "
              << uTime << " sec\n";
      if (u.ru_maxrss > 0)
         cout << "| Max. resident set size:              "
              << u.ru_maxrss << " bytes\n";
      if (u.ru_ixrss > 0)
         cout << "| Integral shared memory size:         "
              << u.ru_ixrss << " bytes\n";
      if (u.ru_idrss > 0)
         cout << "| Integral unshared memory size:       "
              << u.ru_idrss << " bytes\n";
      if (u.ru_isrss > 0)
         cout << "| Integral unshared stack size:        "
              << u.ru_isrss << " bytes\n";
      if (u.ru_minflt > 0)
         cout << "| Page reclaims:                       "
              << u.ru_minflt << endl;
      // Page faults: trials of writing accesses on RO areas,
      //              usually blocked by (I/O) threads
      if (u.ru_majflt > 0)
         cout << "| Page faults:                         "
              << u.ru_majflt << endl;
      if (u.ru_nswap > 0)
         cout << "| Number of swaps                      "
              << u.ru_nswap << endl;
      if (u.ru_inblock > 0)
         cout << "| Number of block input operations:    "
              << u.ru_inblock << endl;
      if (u.ru_oublock > 0)
         cout << "| Number of block output operations:   "
              << u.ru_oublock << endl;
      if (u.ru_msgsnd > 0)
         cout << "| Messages sent:                       "
              << u.ru_msgsnd << endl;
      if (u.ru_msgrcv > 0)
         cout << "| Messages received:                   "
              << u.ru_msgrcv << endl;
      if (u.ru_nsignals > 0)
         cout << "| Signals received:                    "
              << u.ru_nsignals << endl;
      if (u.ru_nvcsw > 0)
         cout << "| Voluntary context switches:          "
              << u.ru_nvcsw << endl;
      if (u.ru_nivcsw > 0)
         cout << "| Involuntary context switches:        "
              << u.ru_nivcsw << endl;
      cout << SX_SEPARATOR;
      cout << endl;
#  endif
}

/// Start the global timer
namespace Timer  {
   static class StartTotal {
      public:
         StartTotal ()  {
            SX_START_TIMER (Total);
         }
   } startTotalNow;
}

void SX_EXPORT_UTIL printTiming (bool restart, bool final)
{
   SX_STOP_TIMER (Timer::Total);
   SxTimer &sxGlobalTimer = SxTimer::getGlobalTimer ();
   sxGlobalTimer.print (SxTimer::getTimerId(Timer::Total),final);
   if (restart)  { SX_START_TIMER (Timer::Total); }
   else          { sxGlobalTimer.printUsage (); }
}
