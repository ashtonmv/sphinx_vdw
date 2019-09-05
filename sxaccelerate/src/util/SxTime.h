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

#ifndef _SX_TIME_H_
#define _SX_TIME_H_

#include <SxUtil.h>
#include <stdlib.h>
#include <SxString.h>
#include <SxTypeDefs.h>

#ifndef WIN32
#  include <sys/resource.h>
#  include <sys/time.h>
#endif

#ifdef WIN32
#  include <time.h>
#  include <windows.h> /* Sleep() */
#endif

#ifdef WIN32
#  ifndef MSVC
#     ifndef ctime_r
#        include <string.h>
         // see MINGW pthread.h:
#        define ctime_r(_clock, _buf) \
                 (  strcpy( (_buf), ::ctime( (_clock) ) ), \
                    (_buf) )
#     endif
#  endif
#endif

#ifdef MACOSX
#  include <mach/mach_time.h>
#endif

/** A time utility
    \see sxtime.cpp
    \ingroup group_os
    \author Sixten Boeck, boeck@mpie.de
    \author Vaclav Bubnik, bubnik@mpie.de */
class SxTime
{
   public:

   /** Get wall clock time in seconds.
       Human perception of the passage of time.
       \return the number of seconds since the Epoch. */
   static uint64_t getEpoch ();

   /** Get wall clock time in fractional seconds.
       \return the number of seconds since the Epoch, with high precision. */
   static double getRealTime ();

   /** Get process time in seconds.
       CPU time spent executing user instructions.
       \return the number of seconds spent in process. */
   static double getProcessTime ();

   /** Get wall clock time stamp.
       \return OS and CPU dependent value usually with no time units.

       \par Example:
       \code
         uint64_t start, end, diff;

         start = SxTime::getRealTimeTicks ();
         // work ...
         end = SxTime::getRealTimeTicks ();
         diff = SxTime::ticksToNanoseconds (start, end);

         printf ("work time %d ns\n", (int)diff);
       \endcode */
   static uint64_t getRealTimeTicks ();

   /** Get the difference between two time stamps in nanoseconds.
       \return the number of nanoseconds. */
   static uint64_t ticksToNanoseconds (uint64_t start, uint64_t end);

   /** Causes the calling Thread to sleep.
       The original usleep() [micro sleep] is obsolete. */
   static void msleep (unsigned int milliseconds);

   /** Causes the calling Thread to sleep. */
   static void sleep (double seconds);

   /** Platform independent version of ctime.
       example: "Thu Jan  1 01:00:00 1970" or "Thu Jan 01 01:00:00 1970" */
   static SxString ctime (long secondsSince1970);

   /** Platform independent version of strftime
     */
   static SxString strftime (const SxString &fmt);
};


inline uint64_t SxTime::getEpoch ()
{
#  ifdef WIN32
      FILETIME tp;
      GetSystemTimeAsFileTime (&tp);
      uint64_t stamp = ((uint64_t)tp.dwHighDateTime << 32) + tp.dwLowDateTime;
      const uint64_t stampSeconds = (uint64_t)stamp * 1e-7; // stamp(1) == 100ns
      const uint64_t stamp1970 = 11644473600; // year 1970-1601 seconds: time()
      return stampSeconds - stamp1970; // n seconds since January 1 1970
#  else
      struct timeval tp;
      gettimeofday (&tp, NULL);        // n seconds since January 1 1970
      return (uint64_t)tp.tv_sec;
#  endif
}

inline double SxTime::getRealTime ()
{
#  ifdef WIN32
      FILETIME tp;
      GetSystemTimeAsFileTime (&tp);
      uint64_t stamp = ((uint64_t)tp.dwHighDateTime << 32) + tp.dwLowDateTime;
      const double stampSeconds = (double)stamp * 1e-7; // stamp(1) == 100ns
      const double stamp1970 = 11644473600.0; // year 1970-1601 seconds: time()
      return stampSeconds - stamp1970; // n seconds since January 1 1970
#  else
      struct timeval tp;
      gettimeofday (&tp, NULL);        // n seconds since January 1 1970
      return ((double)tp.tv_sec) + (1e-6 * (double)tp.tv_usec);
#  endif
}


inline double SxTime::getProcessTime ()
{
#ifdef WIN32
   FILETIME create, exit, kernel, user;
   GetProcessTimes (GetCurrentProcess(), &create, &exit, &kernel, &user);
   uint64_t stamp = ((uint64_t)user.dwHighDateTime << 32) + user.dwLowDateTime;
   double stampSeconds = (double)stamp * 1e-7; // stamp(1) == 100ns
   return stampSeconds;

   //return (((double)clock ()) / CLOCKS_PER_SEC);
#else
   struct rusage usage;
   getrusage (RUSAGE_SELF, &usage);
   return (double)usage.ru_utime.tv_sec + 1e-6 * (double)usage.ru_utime.tv_usec;
#endif
}


inline uint64_t SxTime::getRealTimeTicks ()
{
#ifdef MACOSX
   return mach_absolute_time ();
//#elif defined(__linux__)
//   // requires librt, LDFLAGS += -lrt
//   struct timespec tp;
//   clock_gettime (CLOCK_REALTIME, &tp); // returns success:0, error:-1
//   return tp.tv_sec * 1e9 + tp.tv_nsec;
#elif defined(WIN32)
   FILETIME tp;
   GetSystemTimeAsFileTime (&tp);
   uint64_t stamp = ((uint64_t)tp.dwHighDateTime << 32) + tp.dwLowDateTime;
   return stamp * 100;
#else
   struct timeval tp;
   gettimeofday (&tp, NULL);
   return (uint64_t)((double)tp.tv_sec * 1e9 + (double)tp.tv_usec * 1e3);
#endif
}


inline uint64_t SxTime::ticksToNanoseconds (uint64_t start, uint64_t end)
{
#ifdef MACOSX
   //
   // http://developer.apple.com/mac/library/qa/qa2004/qa1398.html
   //
   mach_timebase_info_data_t tinfo = {0, 0};
   mach_timebase_info (&tinfo);
   double diffns = (end - start) * (double)tinfo.numer / tinfo.denom;
   return static_cast<uint64_t> (diffns);
#else
   return end - start;
#endif
}


inline void SxTime::msleep (unsigned int milliseconds)
{
#  ifdef WIN32
      // --- Sleep(): Causes the calling Thread to sleep.
      //              The parameter is the number of milliseconds to sleep.
      //              sleep() and usleep() are not available
      Sleep (milliseconds);
#  else
      struct timespec req;// = {0};
      time_t sec = (time_t)((long)milliseconds / 1000L);
      long msec = (long)milliseconds - ((long)sec*1000L);
      req.tv_sec = sec; 
      req.tv_nsec = msec * 1000000L;
      // --- nanosleep(): Suspends the execution of the calling thread.
      //                  struct timespec {
      //                     time_t tv_sec;   seconds
      //                     long   tv_nsec;  nanoseconds
      //                  }
      //                  nanosleep (2):
      //                  usleep() [micro sleep] is obsolete
      nanosleep (&req, &req);
#  endif
}


inline void SxTime::sleep (double seconds)
{
   SX_CHECK (seconds >= 0., seconds);
   SxTime::msleep ((unsigned int)(seconds * 1e3));
}

// long type should be enough until 2038
inline SxString SxTime::ctime (long secondsSince1970)
{
   SX_CHECK (secondsSince1970 >= 0, secondsSince1970);
   time_t time = secondsSince1970;
#  ifdef MSVC
      const int DATE_SIZE=26;  // see _wctime_s man page (msdn)
      wchar_t expBuf[DATE_SIZE];
      char res[DATE_SIZE];

      _wctime_s (expBuf, DATE_SIZE, &time);
      for (int c=0; c < DATE_SIZE; ++c)  res[c] = expBuf[c];
      if (res[24]=='\n') res[24]='\0';
      return SxString(res);
#  else
      char buf[26];
      ctime_r(&time, buf);
      if (buf[24]=='\n') buf[24]='\0';
      return SxString(buf);
#  endif
}

inline SxString SxTime::strftime (const SxString &fmt)
{
   struct tm tmTime;
#  ifdef WIN32
      __int64 lTime = 0;
      _time64 (&lTime);
      //_gmtime64_s (&time, &lTime);   // without timezone correction
      _localtime64_s (&tmTime, &lTime);  // uses local timezone
#  else
      time_t lTime = time (NULL);
      //gmtime_r (&lTime, &tmTime);    // without timezone correction
      localtime_r (&lTime, &tmTime);   // uses local timezone
#  endif
   // YYYY-MM-DD HH:MM:SS

   char buffer[1024];
   ::strftime (buffer, 1024, fmt.ascii(), (const struct tm *)&tmTime);
   return SxString (buffer);
}

#endif /* _SX_TIME_H_ */
