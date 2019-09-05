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


#include <SxTimeOut.h>


SxTimeOut::SxTimeOut (SxThread *task_, double sec) 
   : SxThread (),
     task (task_),
     timer (1, true),
     timeout (sec),
     exceeded (false)

{
   // empty
}


SxTimeOut::~SxTimeOut ()
{
   // empty
}

void SxTimeOut::main ()
{
   timer.reset ();
   timer.start ();
   task->start ();
   unsigned int interval = 50;  // msec;
   exceeded = false;
   do {
      SxThread::millisleep (interval);
      timer.refresh ();
      if ( timer.getTime () >= timeout)  exceeded = true;
   } while (!exceeded && task->isRunning());
   timer.stop ();
   task->terminate ();
}

bool SxTimeOut::hasExceeded () const
{
   return exceeded;
}
