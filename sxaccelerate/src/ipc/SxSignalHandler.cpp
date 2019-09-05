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
#include <SxSignalHandler.h>
#include <SxEmergency.h>
#include <SxString.h>
#include <SxConfig.h>
#include <signal.h>

bool SxSignalHandler::firstHandler = true;


SxSignalHandler &SxSignalHandler::getGlobalInstance ()
{
   static SxSignalHandler globalInstance;
   return globalInstance;
}


SxSignalHandler::SxSignalHandler ()
{
   if (signal (SIGINT, SxSignalHandler::trap) == SIG_ERR)  {  // ^C
      sxprintf ("Can't intialize signal handler (SIGINT).\n");
   }
   if (signal (SIGTERM, SxSignalHandler::trap) == SIG_ERR)  {  // kill
      sxprintf ("Can't intialize signal handler (SIGTERM).\n");
   }
#ifndef WIN32
   if (signal (SIGQUIT, SxSignalHandler::trap) == SIG_ERR)  {  // quit
      sxprintf ("Can't intialize signal handler (SIGQUIT).\n");
   }
   if (signal (SIGUSR1, SxSignalHandler::trap) == SIG_ERR)  {  // usr1
      sxprintf ("Can't intialize signal handler (SIGUSR1).\n");
   }
   if (signal (SIGUSR2, SxSignalHandler::trap) == SIG_ERR)  {  // usr2
      sxprintf ("Can't intialize signal handler (SIGUSR2).\n");
   }
#endif
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxSignalHandler::commit (SxEmergency *obj)
{
   objList << obj;
}


void SxSignalHandler::deregister (SxEmergency *obj)
{
   objList.removeElement (obj);
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void SxSignalHandler::trap (int signalId)
{
   if (signalId == SIGINT)  {
      // --- allow to ignore SIGINT
      SxSignalHandler *signalHandler = &getGlobalInstance();
      if (signalHandler->objList.getSize () > 0)  {
         bool doNotInterrupt = false;
         SxList<SxEmergency *>::Iterator it;
         for (it  = signalHandler->objList.begin(); 
              it != signalHandler->objList.end(); 
              it++)  
         {
            if (!(*it)->interrupt ())  {
               doNotInterrupt = true;
            }
         }
         if (doNotInterrupt)  {
            return;
         }
      }   
   }
   
   if (firstHandler)  {
      // --- in case of many threads call this function only once
      firstHandler = false;

      fflush (NULL); // flush all open streams

      SxSignalHandler *signalHandler = &getGlobalInstance();

      if (signalHandler->objList.getSize () > 0)  {
         SxList<SxEmergency *>::Iterator it;
         for (it  = signalHandler->objList.begin(); 
              it != signalHandler->objList.end(); 
              it++)  
         {
            //sxprintf ("SxSignalHandler::trap disabled\n");
            (*it)->dumpData (signalId);
         }
      }
      //cout << "\nProcess terminated due to external request.\n";

      fflush (NULL); // flush all open streams
      _Exit(1);
   }
}

