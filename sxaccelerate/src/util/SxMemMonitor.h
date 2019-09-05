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

#ifndef _SX_MEM_MONITOR_H_
#define _SX_MEM_MONITOR_H_

#include <SxString.h>
#include <SxArray.h>
#include <SxMap.h>
#include <SxUtil.h>
#include <stdio.h>


/** \brief monitors memory consumption of all memory consumer classes

    \b SxMemMonitor = SFHIngX Memory Monitor

    This class manages all memory consumers (SxMemConsumer) and print
    the memory statistics.

    For more information see also in the description of SxMemConsumer.

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_UTIL SxMemMonitor
{
   public:

      /** \brief Default constructor

          The default constructor allocates the program global memory 
          counters. */
      SxMemMonitor ();

     ~SxMemMonitor ();

      /** \brief returns the next free memory observer id

          All observers are organized in arrays. This function returns
          the number of the next free slot and increments the internal
          counter lastId. */
      int getFreeId ();

      /** \brief Returns the index in the maxMemCounter array 

          A memory consumer observer is identified by the observer
          object's address and the line number where the observer is
          being registered. This function tries to find such an
          observer-line pair in the local list. If it succeeded it
          returns the index in the maxMemCounter array.
          Otherwise a new observer-line pair is being created and
          assigned with a new empty free id. */
      int registerMemObserver (const char *className,
                               const char *file, int line);

      /** \brief update memory information about a consumer

          This function checks whether the provided memory count of an
           memory observer exceeds the current one. If so the higher value
           is taken. 
           \param id       Memory observer id
           \param nBytes   currently occupied space of the observer
       */
      void update (int id, size_t nBytes) const;


      static void print ();

      /** \brief Print the memory statistics - sorted by memory usage. 
       
          Don't use this function directly. Instead call SxMemMonitor::print.
       */
      void printStatistics ();

   protected:

      int lastId;
      const int nCounters;

      SxMap<SxString, int>     consumerMap;    // {file:line} |-> idx
      mutable SxArray<size_t>  maxMemCounter;  // mutable, see
                                               // SxMemConsumer::trackMemory
      SxArray<SxString>        classNames;
      SxArray<SxString>        fileLines;


      SxString getMemoryString (ssize_t nBytes) const;

};

extern SxMemMonitor sxMemMonitor;

#endif /* _SX_MEM_MONITOR_H_ */
