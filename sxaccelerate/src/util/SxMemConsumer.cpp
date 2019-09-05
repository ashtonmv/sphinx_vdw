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

#include <SxMemConsumer.h>
#include <SxMemMonitor.h>
#include <SxError.h>


SxMemConsumer::SxMemConsumer ()
#ifdef USE_MEMTRACKING
   : className(NULL),
     memCounterId(0)
#endif /* USE_MEMTRACKING */
{
   // empty
}

#ifdef USE_MEMTRACKING   
SxMemConsumer::SxMemConsumer (const char *className_)
   : className(className_),
     memCounterId(0)
#else   
SxMemConsumer::SxMemConsumer (const char *)
#endif /* USE_MEMTRACKING */
{
   // empty
}

SxMemConsumer::~SxMemConsumer ()
{
   // empty
}

#ifdef USE_MEMTRACKING
void SxMemConsumer::trackMemory (size_t nBytes) const
{
// sxprintf ("update: %d bytes\n", nBytes);
   if (memCounterId)  {
      SX_CHECK (memCounterId != -1);
//    sxprintf ("updating %d = %d bytes\n", memCounterId, nBytes);
      sxMemMonitor.update (memCounterId, nBytes);
   }
}




void SxMemConsumer::registerMemObserver (const char *name, 
                                         const char *file, int line)
                                         
{
   memCounterId = sxMemMonitor.registerMemObserver (name, file, line);
// sxprintf ("register: [%s]%s:%d -> %d\n", name, file, line, memCounterId);
}



void SxMemConsumer::setMemCounterId (size_t id)
{
   SX_CHECK (memCounterId == 0, memCounterId);
   memCounterId = id;
}

#endif /* USE_MEMTRACKING */
