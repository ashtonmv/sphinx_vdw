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
#ifndef _SX_SIGNAL_HANDLER_H_
#define _SX_SIGNAL_HANDLER_H_

#include <SxList.h>
#include <SxIPC.h>

class SxEmergency;


/**
  This object manager has been designed for stop a program gracefully.
  @author  Sixten Boeck
  @ingroup group_os
  */
class SX_EXPORT_IPC SxSignalHandler
{
   public:
      SxSignalHandler ();

      static SxSignalHandler &getGlobalInstance ();

      void commit (SxEmergency *obj);
      void deregister (SxEmergency *obj);
      
      static void trap (int);

      /** \brief Identifier if this is the first signal handler

          In threadded applications every thread catches a signal and
          hence, the routine ::trap would be called multiple times. In order
          to avoid this it is necessary to call ::trap only for the first
          signal handler.
          The variable ::firstHandler is initialized from ::initSPHInXCOM. */
      static bool firstHandler;

   protected:
      SxList<SxEmergency *> objList;
};

#endif /* _SX_SIGNAL_HANDLER_H_ */
