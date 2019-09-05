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

#ifndef _SX_TIMER_THREAD_H_
#define _SX_TIMER_THREAD_H_

#include <SxSystemThread.h>
#include <SxThreadCond.h>

/** \brief Timer

\code
class Timer : public SxTimerThread
{
   public:
      virtual void interval () { std::cout << "interval" << std::endl; }
      virtual void timeout () { std::cout << "timeout" << std::endl; }
};

int main ()
{
   Timer timer;
   timer.setInterval (1.);
   timer.setTimeout (2.);

   timer.start ();
   timer.wait ();
   
   // prints:     time
   //  interval    1 s
   //  interval    2 s
   //  timeout     2 s
      
   return 0;
}
\endcode

    \author Vaclav Bubnik, bubnik@gemmantics.com */
class SX_EXPORT_IPC SxTimerThread : public SxSystemThread
{
   public:
      SxTimerThread ();
      virtual ~SxTimerThread ();
      
      void setInterval (double);
      void setTimeout (double);
      
      virtual void start (Priority priority=SameAsParent);
      void stop ();
      
      virtual void main ();


      virtual void interval () { /* empty */ };
      virtual void timeout () { /* empty */ };
 
   signals:
      SxSignal<double, const char *> SX_SIGNAL(signalInterval);
      SxSignal<double, const char *> SX_SIGNAL(signalTimeout);

   protected:
      double intervalSec; // delay in seconds between each call (if set)
      double timeoutSec;  // how many seconds to run the timer (if set)
      bool stopCmd;
      SxMutex mutex;
      SxThreadCond condition;
};

#endif /* _SX_TIMER_THREAD_H_ */
