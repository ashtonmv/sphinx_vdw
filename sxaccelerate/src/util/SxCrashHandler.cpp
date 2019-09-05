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

#include <SxCrashHandler.h>
#include <SxCrashListener.h>
#include <SxString.h>
#include <stdio.h>
#include <stdlib.h>

SxCrashHandler &SxCrashHandler::getGlobalInstance ()
{
   static SxCrashHandler handler;
   return handler;
}

SxCrashHandler::SxCrashHandler ()
{
   firstHandler = true;
}

SxCrashHandler::~SxCrashHandler ()
{
   // empty
}


void SxCrashHandler::commit (SxCrashListener *obj)
{
   listeners << obj;
}


void SxCrashHandler::deregister (SxCrashListener *obj)
{
   listeners.removeElement (obj);
}

void SxCrashHandler::trapCrash ()
{
   if (getGlobalInstance().firstHandler)  {
      // --- in case of many threads call this function only once
      getGlobalInstance().firstHandler = false;

      fflush (NULL); // flush all open streams

      if (getGlobalInstance().listeners.getSize() > 0)  {
         cout << endl << endl;
         cout << SX_SEPARATOR;

         sxprintf ("| CRASH HANDLER: Commencing termination of program.\n");
         sxprintf ("|                Shutting down %d object(s).\n",
                 (int)getGlobalInstance().listeners.getSize());

         SxList<SxCrashListener *>::Iterator it;
         for (it  = getGlobalInstance().listeners.begin(); 
              it != getGlobalInstance().listeners.end(); 
              it++)  
         {
            (*it)->handleCrash ();
         }
         sxprintf ("| CRASH HANDLER: Shutdown sequence finished.\n\n");
         cout << SX_SEPARATOR;
      }
   }
}
