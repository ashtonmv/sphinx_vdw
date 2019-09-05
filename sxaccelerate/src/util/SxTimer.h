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

#ifndef _SX_TIMER_H_
#define _SX_TIMER_H_

#include <SxString.h>
#include <SxConfig.h>
#include <SxUtil.h>
#include <SxArray.h>
#include <SxMap.h>
#include <typeinfo>

/** The following map is needed to circumvent problems if
    SxTimer::getStartId<T> is instantiated more than once, e.g. in MSVS DLL's.
    Only the first instance gets the start index from the global timer object,
    all others acquire it from the map.
*/
extern SX_EXPORT_UTIL SxMap<SxString, int> registeredTimerEnums;

/** A stop watch class.

    \ingroup  group_os
    \author   Sixten Boeck
  */
class SX_EXPORT_UTIL SxTimer
{
   public:

      int nTimers;
      SxArray<double>   accTime;
      SxArray<double>   procTime;
      SxArray<int>      nCalls;
      SxArray<bool>     isRunning;
      SxArray<SxString> timerName;
      bool              realTime;

      SxTimer ();
      /** Instantiate a set of new timers. One can choose whether the
          user time or the realtime should be taken. */
      explicit SxTimer (int nTimers, bool realTime=false);
      ~SxTimer ();

      void init    (int nTimersIn, bool realTime=false);
      void start   (int iTimer=0);
      void stop    (int iTimer=0);
      void refresh ();
      void reset   (int iTimer=0);

      void setName (int, const SxString &);

      void resize (int newSize);

      int getSize () const { return nTimers; }

      double getTime (int iTimer=0) const;
      /** return time string in format mm:ss */
      static SxString toMinutes (double);
      double getAvgTime (int iTimer=0) const;
      int    getNCalls (int iTimer=0) const;
      bool   isIdle (int iTimer=0) const;
      SxString getName (int iTimer=0) const;

      void print (int iTlTimer=-1, bool final=false) const;
      void printUsage () const;

      /** This is the (platform dependent) function call to read out the
          current time. Denpding on the system different routines will be
          used in order to get the user time:
          - Unix: User time taken by this process or real time
          - Windows: Approximation of the processor time used this program

          Note that in the windows version an overflow of clock() occurs
          every 72 minutes. However, as long as the at least one timer
          is being updated within this period this class overcomes this
          problem internally. */
      double getClock () const;

      ///@name Global timer
      ///@{
      inline static SxTimer& getGlobalTimer ()
      {
         static SxTimer globalTimer (64); // start with 64 timer slots
         return globalTimer;
      }

   protected:
      /// Get the number of registered timers
      inline static int &getMaxTimerId ()  {
         static int maxTimerId = -1;
         return maxTimerId;
      }

   public:
      /// Get the offset in the timer list for particular enum type
      template<class T>
      static int getStartId ()  {
         static int startId 
            = registeredTimerEnums.containsKey (typeid(T).name ())
            // 2nd instance: get startId from map
            ? registeredTimerEnums(typeid(T).name ())
            // 1st instance: get startId from getMaxTimerId and set it in map
            : ( registeredTimerEnums(typeid(T).name ()) = getMaxTimerId () + 1);
         return startId;
      }

      /// Get timer id for an enum
      template <class T> inline static int getTimerId (const T &);
      ///@}

};


/**
  \brief This is a timer enum registration class for the global timer

  \example
  \code
  namespace Timer {
     enum MyTimerEnum { Clock1, Clock2, Fast, Slow};
     }
REGISTER_TIMERS(MyTimerEnum)
{
   using namespace Timer;
   regTimer (Clock1, "Clock 1");
   regTimer (Clock2, "Clock 2");
   regTimer (Fast,   "Clock for the fast part");
   regTimer (Slow,   "Clock for the slow part");
}
   \endcode

   It is recommended to put timer enums into namespace Timer for readability,
   but this is not technically required.

   \note The enum must not have assigned values.
  */
template <class T>
class SxTimerNames  {
   public:
      /** \brief Constructor
          \note must be template-specialized for each enum using REGISTER_TIMERS
          */
      SxTimerNames ();
   protected:
      /// Register name associated with enum id
      static void regTimer (const T& id, const SxString &name);
      friend class SxTimer;
      /** \brief Update and check id (for debug mode)
          @param id Enum id to be checked (or updated)
          @return false if this higher than all ids ever seen by checkId,
                  true otherwise
        */
      static inline bool checkId (const T& id)  {
         static int maxId = 0;
         bool ok = (int(id) <= maxId);
         if (!ok) maxId = id;
         return ok;
      }
};

#define SX_REGISTER_TIMERS(T) template <> inline SxTimerNames<T>::SxTimerNames ()

template<class T>
void SxTimerNames<T>::regTimer (const T &id, const SxString &name)
{
   SX_CHECK (int(id) >= 0, int(id));

   // --- do not use getTimerId here to avoid infinite loop
   int uniqueId = SxTimer::getStartId<T>() + int(id);

   checkId (id); // update maxId

   SxTimer &globalTimer = SxTimer::getGlobalTimer ();
   // resize globalTimer if necessary
   // if resizing, increase size but at least factor 2 (resizing is slow)
   if (globalTimer.getSize () <= uniqueId)  {
      int minSize = uniqueId + 1, twiceCurrent = globalTimer.getSize () * 2;
      globalTimer.resize ((minSize > twiceCurrent) ? minSize : twiceCurrent);
   }
   globalTimer.setName (uniqueId, name);
}

template<class T>
int SxTimer::getTimerId (const T &id)
{
   // this initializes the names
   static SxTimerNames<T> xx;
   // this makes sure we don't go below startId
   SX_CHECK ((int) id >= 0, id);
   // if code stops here, an unregistered id was used
   SX_CHECK (SxTimerNames<T>::checkId (id), id);
   return getStartId<T> () + int(id);
}

/** \brief Time taking class */
class SxTakeTime {
   int id;
   public:
      SxTakeTime (int id_) : id(id_) {
         // start timer if not yet in use
         if (SxTimer::getGlobalTimer ().isIdle (id))
            SxTimer::getGlobalTimer ().start (id);
         else
            id = -1;
      }
      ~SxTakeTime () {
         if (id >= 0) SxTimer::getGlobalTimer ().stop (id);
      }
};

namespace Timer  { enum TotalTime { Total }; }
SX_REGISTER_TIMERS (Timer::TotalTime) { regTimer (Timer::Total, "Total Time"); }

#   define SX_START_TIMER(id)  SxTimer::getGlobalTimer ()\
                            .start(SxTimer::getTimerId (id))
#   define SX_STOP_TIMER(id)   SxTimer::getGlobalTimer ()\
                            .stop (SxTimer::getTimerId (id))
// Do NOT ask about the complicated definition for CLOCK.
// It has to do with the right order of expansion and concatenation.
#   define UNIQUE_TIMER_NAME(x) timer_line ## x
#   define UNIQUE_TIMER(x)  SxTakeTime UNIQUE_TIMER_NAME(x)
#   define SX_CLOCK(id)        UNIQUE_TIMER(__LINE__) (SxTimer::getTimerId(id));
#   define GETTIME(id)      (SxTimer::getGlobalTimer ()\
                            .getTime(SxTimer::getTimerId (id)))

/** Print out the timing
    @param restart If true, this is an intermediate timing, so keep the total
           time timer running. Otherwise, also print out OS-provided program
           statistics.
    @param final If true, precedes all timer names by Final
  */
void SX_EXPORT_UTIL printTiming (bool restart = false, bool final=false);

#endif // _SX_TIMER_H_
