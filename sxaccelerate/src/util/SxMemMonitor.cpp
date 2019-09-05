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

#include <SxMemMonitor.h>
#include <stdio.h>


SxMemMonitor sxMemMonitor;


SxMemMonitor::SxMemMonitor ()
   : nCounters(100)
{
   maxMemCounter.resize (nCounters);
   classNames.resize (nCounters);
   fileLines.resize  (nCounters);
   for (int i=0; i < nCounters; ++i)  {
      maxMemCounter(i) = 0;
      classNames(i) = "";
      fileLines(i)  = "";
   }
}


SxMemMonitor::~SxMemMonitor ()
{
   // empty
}


int SxMemMonitor::getFreeId ()
{
   return ++lastId;
}


int SxMemMonitor::registerMemObserver (const char *className, 
                                       const char *file, int line)
{
   int id;
   SxString fileLine = SxString() + file + ":" + line;
   if ( !consumerMap.containsKey (fileLine) )  {
      id = getFreeId ();
      consumerMap(fileLine) = id;
   }  else  {
      id = consumerMap(fileLine);
   }

   SX_CHECK (id >= 0 && id < nCounters, id, nCounters); 
   classNames(id) = className;
   fileLines(id)  = fileLine;

   return id;
}


void SxMemMonitor::update (int id, size_t nBytes) const
{
   if (nBytes > maxMemCounter(id))  maxMemCounter(id) = nBytes;
}


SxString SxMemMonitor::getMemoryString (ssize_t nBytes) const
{
   SxString res;
   if (nBytes <    1024)  {
      res = SxString (nBytes)                  + " By";
   }  else if (nBytes < 1048576)  {
      res = SxString (static_cast<double>(nBytes)/1024., "%.1f")    + " KB";
   }  else  {
      res = SxString (static_cast<double>(nBytes)/1048576., "%.1f") + " MB";
   }
   return res;
}


void SxMemMonitor::print ()
{
   sxMemMonitor.printStatistics ();
}

void SxMemMonitor::printStatistics ()
{
#ifdef USE_MEMTRACKING
   sxprintf ("\n+---------------------------------------"
           "--------------------------------------\n");
   sxprintf ("| MEMORY CONSUMPTION STATISTICS\n");
   sxprintf ("+-------------------------------------+"
           "---------------------------+-----------\n");
   sxprintf ("| %-35s | %-25s | %s\n", 
           "Consumer", "File:Line", "Memory");
   sxprintf ("+-------------------------------------+"
           "---------------------------+-----------\n");
   SxString name;
   // --- prepare sorting by consumed memory, reverse order
   SxArray<ssize_t> idxArray = maxMemCounter.getSortIdx();

   // --- print out statistics
   int i, j;
   for (j=0; j < nCounters; ++j)  {
      i = j; //idxArray(nCounters-j-1);
      if (classNames(i).getSize() > 0)  {

         sxprintf ("| %-35s ", classNames(i).ascii());
         sxprintf ("| %-25s ", fileLines(i).ascii());

         // --- use proper units to save space in print line
         sxprintf ("| %s\n", getMemoryString (maxMemCounter(i)).ascii());
      }
   }
   sxprintf ("+-------------------------------------+"
           "---------------------------+-----------\n");
#endif /* USE_MEMTRACKING */
}
