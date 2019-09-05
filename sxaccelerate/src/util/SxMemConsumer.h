
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

#ifndef _SX_MEM_CONSUMER_H_
#define _SX_MEM_CONSUMER_H_

#include <SxConfig.h>
#include <SxUtil.h>
#include <stdio.h>


template<class T> 
size_t getNBytes (const T &) 
{ 
   return sizeof(T); 
}
//template<> size_t getNBytes<SxComplex8> (const SxComplex8 &) 
//{ return 2 * sizeof(float); }

/** \brief keeps track of memory consumption

    \b SxMemConsumer = SFHIngX Memory Consumption Monitor

    This class keeps track of memory allocations of all classes derived
    from this base class.

    \par Example: Tracking the memory consumption of a member variable

    Assuming a class A which should keep track of its memory consumption:

\code
   class A : public SxMemConsumer
   {
      public:
         SxList<int> dataA;

         A () : SxMemConsumer ("A")
         {
            TRACK_MEMORY (dataA);
            dataA.resize (2)
         }
    };
\endcode

    In the above example the class A has been derived from SxMemConsumer.
    That creates a new memory oberver in the memory monitor SxMemMonitor.
    In order to generate an useful memory statistics printout in the constructor
    the class name "A" has been provided.
    The variable which should be analyzed during the program run is dataA.
    It can be switched in using the macro TRACK_MEMORY. Once it has been
    called all memory allocations/deletions/reallocations are analyzed and
    the largest memory consumption is being stored in the memory monitor.

    \par Example: Tracking the memory consumption of a local variable

    The macro TRACK_MEMORY works also for local and temporary variables
    as shown in the following example:

\code
   class B : public SxMemConsumer
   {
      public:
         A           a;
         SxList<int> dataB;

         B () : SxMemConsumer ("B")
         {
            TRACK_MEMORY (a);
         }

         void foo ()
         {
            SxList<int> tmpList;
            TRACK_MEMORY (tmpList);
         }
         
   };
\endcode

   Here, in foo(), the local variable tmpList is being analyzed during the
   entire program run.

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_UTIL SxMemConsumer
{
   public:
      // constructor used by SxArray, SxVector, SxList
      SxMemConsumer (); 
      // constructor used by actual monitored classes
      SxMemConsumer (const char *className_);
     ~SxMemConsumer ();

#     ifdef USE_MEMTRACKING

      void trackMemory (size_t) const;

      /** \brief Track memory allocations of a variable

          Do not use this function directly. Instead use the macro
          TRACK_MEMORY.

          \par Example
\code
   class A : public SxMemConsumer
   {
      SxList<int>  data;
      A() : SxMemConsumer ("A")
      {
         TRACK_MEMORY (data);
      }
   };
\endcode
      */
      void registerMemObserver (const char *name, const char *file, int line);


   protected:

      const char *className;
      mutable int memCounterId;  // trackMemory might be called from const()
                                 // functions of derived classes -> mutable!

      void setMemCounterId (size_t);
#     endif /* USE_MEMTRACKING */
};

#ifdef USE_MEMTRACKING
/** \brief Keep track of the memory consumption of a variable

    This macro can be used in all classes derived from SxMemConsumer. 
    Please find a detailed documentation about its usage with some examples
    in the documentation of SxMemConsumer */
#  define TRACK_MEMORY(var)         var.registerMemObserver(#var,            \
                                                        __FILE__, __LINE__); \
                                    var.trackMemory (getNBytes(var));
/** \brief  Update a single memory observer

    Call this macro after an element of a recursively defined datatype such
    as SxList<SxArray>, has been resized. */
#  define UPDATE_MEMORY(var)    var.trackMemory (getNBytes(var));    
/** \brief Update memory information of an memory observer

    This macro is used for the service classes which are able to track
    memory information (SxList, SxArray, ...).  */
#  define TRACK_MALLOC(var,nElem)   trackMemory (getNBytes(var) * nElem);
#else
#  define TRACK_MEMORY(var)         ((void)0)
#  define UPDATE_MEMORY(var)    ((void)0)
#  define TRACK_MALLOC(T, nElem)    ((void)0)
#endif /* USE_MEMTRACKING */

#endif /* _SX_MEM_CONSUMER_H_ */
