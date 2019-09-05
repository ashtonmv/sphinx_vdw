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
#ifndef _SX_EMERGENCY_H_
#define _SX_EMERGENCY_H_

#include <SxSignalHandler.h>

/**
  @ingroup Communication
  */
class SX_EXPORT_IPC SxEmergency
{
   public:
      SxEmergency () { SxSignalHandler::getGlobalInstance().commit (this); }
      virtual ~SxEmergency () { };

      virtual void dumpData (int)=0;
      
      /** \brief Allows to ignore SIGINT and to receive multiple SIGINTs.

\code
bool ServerEmergency::interrupt ()
{
   static int n = 0; // counter to ignore SIGINT once
   if (n > 0)  {
      return true;
   }
   n++;
   
   cout << "connect to all clients and cancel commands" << endl;
   
   return false;
}
\endcode
*/
      virtual bool interrupt () { return true; }
};

#endif /* _NG_EMERGENCY_H_ */
