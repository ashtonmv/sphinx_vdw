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

#ifndef _SX_CRASH_LISTENER_H_
#define _SX_CRASH_LISTENER_H_

#include <SxCrashHandler.h>
#include <SxUtil.h>

/** \brief Crash Listener

    \b SxCrash = SPHInX Crash and Exit Listener

    Classes derived from SxCrashListener can react on an sudden program
    termination caused by an unexpected SX_EXIT or SX_CHECK*. The child class
    simply overloads the purely abstract function handleCrash.
    For example, a class which creates temporary files should removed those
    in case of a crash.

    The usage is analog to that of SxEmergency.

    \par Example:
\code
   #include <SxCrashListener.h>

   class SxMyClass : public SxCrashListener
   {
      public:
         SxMyClass () : SxCrashListener () { }
         void handleCrash () { cout << "SxMyClass::handleCrash" << endl; }
   };
\endcode

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_UTIL SxCrashListener
{
   public:
      SxCrashListener () { SxCrashHandler::getGlobalInstance().commit (this);     }
      virtual ~SxCrashListener () { SxCrashHandler::getGlobalInstance().deregister (this); }

      virtual void handleCrash ()=0;
};

#endif /* _SX_CRASH_LISTENER_H_ */
