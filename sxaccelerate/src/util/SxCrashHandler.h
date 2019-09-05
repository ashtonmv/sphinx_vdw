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

#ifndef _SX_CRASH_HANDLER_H_
#define _SX_CRASH_HANDLER_H_

#include <SxList.h>
#include <SxUtil.h>

class SxCrashListener;

/** \brief Crash and exit handler

    \b SxCrashHandler = SPHInX Crash and Exit Handler

    This handler manages a set of registered receivers/callbacks which 
    are called in case of
    - SX_EXIT statement or
    - SX_CHECK* statement
    was executed.

    An useful application of this handler class is for example a cleanup
    of temporary files after a crash.

    Note, that only SPHInX-assert crashes can be trapped here. Here
    segfault or bus error is not handled here.
    In this case SxSignalHandler might be used.

    The usage of SxCrashHandler is analog to SxSignalHandler. Please
    also read the documentation of SxCrashListener.

    \sa SxCrashListener
    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_UTIL SxCrashHandler
{
   public:
      SxCrashHandler ();
     ~SxCrashHandler ();

      void commit (SxCrashListener *);
      void deregister (SxCrashListener *);
      /** \brief Identifier if this is the first signal handler

          In threadded applications every thread catches a crash and
          hence, the routine ::trap would be called multiple times. In order
          to avoid this it is necessary to call ::trap only for the first
          signal handler.
          The variable ::firstHandler is initialized from ::initSxUtil. */
     bool firstHandler;

     static void trapCrash ();

     static SxCrashHandler &getGlobalInstance ();

   protected:
     SxList<SxCrashListener *> listeners;
};

#endif /* _SX_CRASH_HANDLER_H_ */
