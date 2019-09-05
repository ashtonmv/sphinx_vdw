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

#include <SxThread.h>
#include <SxError.h>
#include <SxTime.h>
#include <SxThreadPool.h>

SxThread::SxThread (const SxString &tag_,
                    int priority_,
                    double maxRunTime_,
                    double maxSpoolTime_)
   : SxThis<SxThread>(),
     tag(tag_),
     priority(priority_),
     maxRunTime(maxRunTime_),
     maxSpoolTime(maxSpoolTime_),
     timer(2, true),
     prioPoints(0)
{
   SX_TRACE ();
   SX_CHECK (priority >=0 && priority < MaxPriority, priority);
   timer.start (Spooling);
}

SxThread::SxThread (const SxString            &tag_,
                    const SxPtr<SxThreadPool> &pool_)
   : SxThis<SxThread>(),
     tag(tag_),
     priority(0),
     maxRunTime(-1.),
     maxSpoolTime(-1.),
     timer(2, true),
     prioPoints(0)
{
   SX_TRACE ();
   SX_CHECK (pool_.getPtr ());
   pool = pool_;
   timer.start (Spooling);
}

SxThread::~SxThread ()
{
   SX_TRACE ();
}

void SxThread::setThreadPool (const SxPtr<SxThreadPool> &pool_)
{
   SX_CHECK (pool_.getPtr ());
   pool = pool_;
}

int SxThread::updatePriority ()
{
   ++prioPoints;
   return prioPoints - priority;
}

bool SxThread::spoolingExpired () const
{
   if (maxSpoolTime < 0)  return false;
   return timer.getTime(Spooling) > maxSpoolTime;
}

bool SxThread::runningExpired () const
{
   if (maxRunTime < 0)  return false;
   return timer.getTime(Running) > maxRunTime;
}

void SxThread::start ()
{
   SX_TRACE ();
   SX_CHECK (pool.getPtr ());
   SX_CHECK (!pool->contains ( getThis () ));
   timer.start (Running);
   pool->submit (getThis ());
}


void SxThread::wait ()
{
   SX_TRACE ();
   SX_CHECK (pool.getPtr ());
   pool->wait (getThis ());
}


void SxThread::terminate ()
{
   SX_TRACE ();
   SX_EXIT;
}


bool SxThread::isRunning () const
{
   SX_TRACE ();
   SX_CHECK (pool.getPtr ());
   return pool->isRunning (getThis ());
}


void SxThread::sleep (unsigned int seconds)
{
   SX_TRACE ();
   millisleep (1000U * seconds);
}


void SxThread::millisleep (unsigned int mseconds)
{
   SX_TRACE ();
   SxTime::msleep (mseconds);
}

void SxThread::main ()
{
   SX_TRACE ();
   SX_EXIT;  // this function is not pure abstract to allow static linkage
}

void SxThread::handleException (const SxException &e)
{
   SX_TRACE ();
   cerr << "Exception Stack for SxThread with tag=" << tag
        << endl << e.toString (SxException::DebugStack);

   SX_EXIT;  // this function is not pure abstract to allow static linkage
}

void SxThread::mainWrap ()
{
   SX_TRACE ();

   sigStarted.send (this);

   try {
      main ();
   } catch (SxException e) {
      handleException (e);
   } catch (...) {
      SX_THROW ("caught unknown exception in SxThread with tag=" + tag);
   }

   sigFinished.send (this);
}
