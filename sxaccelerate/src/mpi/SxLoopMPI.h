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
#ifndef _SX_LOOPMPI_
#define _SX_LOOPMPI_

#include <SxConfig.h>
#include <SxAutoLoop.h>

#include <SxParallelHierarchy.h>
#include <SxTaskGroup.h>
#include <SxStack.h>
#include <SxNArray.h>

/// Whether only task group master or all ranks participate in data transfer
enum SxMPITaskSelection { TaskGroupMaster, TaskGroupAll };
/// Which MPI task group level participates in data transfer
enum SxMPILevelSelection {
   TopLevel, CurrentLevel, SpecifiedName, LevelSiblings
};

// --- debug macros ---
/*

#define SX_CHECKPOINT(id) \
{ \
   SxLoopMPI::barrier(); \
   sxprintf("\nSXCHECKPOINT : %s:%d %s\n\n", \
         __FUNCTION__, __LINE__, SxString(id).ascii()); \
   fflush(stdout); \
   SxLoopMPI::barrier(); \
}

#define SX_PRINTVAR(var) \
{ \
   std::cout << __FUNCTION__ << ":" << __LINE__ << ": " \
             << #var << " = " << var << std::endl << std::flush; \
}

*/

#define SX_CHECKPOINT(id) {}
#define SX_PRINTVAR(var) {}
// ------



#ifdef USE_MPI
#include <mpi.h>
#define USE_LOOPMPI
#  define SX_MPI_MASTER_ONLY  if (SxLoopMPI::me () == 0)
#  define SX_MPI_SLAVE_ONLY   if (SxLoopMPI::me () != 0)

#  define SX_NO_MPI \
   if (SxLoopMPI::nr () > 1)  { \
      sxprintf ("This part of the code does not work MPI-parallel.\n");\
      SX_EXIT; \
   }

// Build a unique name for SxMPILevel instantiation, compare <SxTimer.h>.
#   define MPI_LVL_INSTANCE_NAME(x) mpi_lvl_line_ ## x
#   define MPI_LVL_INSTANCE(x)      MPI_LVL_INSTANCE_NAME(x)
#   define SX_MPI_LEVEL(id)         SxMPILevel MPI_LVL_INSTANCE(__LINE__) (id)
#   define SX_MPI_SOURCE SxLoopMPI::OverloadSource  MPI_LVL_INSTANCE(__LINE__)
#   define SX_MPI_TARGET SxLoopMPI::OverloadTarget MPI_LVL_INSTANCE(__LINE__)

#else
#  define SX_MPI_MASTER_ONLY
#  define SX_MPI_SLAVE_ONLY   if (false)

#   define SX_NO_MPI (void)(0)
#   define SX_MPI_LEVEL(id) (void)(0)
#   define SX_MPI_SOURCE SxLoopMPI_nothing
#   define SX_MPI_TARGET SxLoopMPI_nothing
   inline void SxLoopMPI_nothing (const char *, SxMPITaskSelection) {}
   inline void SxLoopMPI_nothing (SxMPILevelSelection,
                                  SxMPITaskSelection = TaskGroupAll) {}
#endif


class SX_EXPORT_MPI SxMPILevel
{
   public:
      SxMPILevel(const SxString &);
      ~SxMPILevel();
};

class SxAutoLevel;

class SX_EXPORT_MPI SxLoopMPI
{
   private:
      /** \brief Semaphore to trigger the use of the parallel 3D FFT,
       *         required to avoid "double" parallelization in certain regions.
       **/
      static bool doParallelFft;

      /** \brief Indicates whether MPI is initialized */
      static bool initialized;

      /// Shutdown class to call finalize automatically
      static class ShutDown {
         public:
            ~ShutDown () {
               SxLoopMPI::finalize ();
            }
      } shutDown;

      /** \brief Stack to store the current level of the hierarchical parallelization  */
      static SxStack<SxTaskGroup*> mpiLevel;

      // SxMPILevel needs to access mpiLevel stack
      friend class SxMPILevel;
      // SxAutoLevel needs to access mpiLevel stack
      friend class SxAutoLevel;

   public:
      /// Get the current MPI level
      static SxMpiTaskGroup* getLevel () {
         if (mpiLevel.isEmpty ()) return NULL;
         return dynamic_cast<SxMpiTaskGroup*>(mpiLevel.top ());
      }


   public:

      /** \brief Initialize the MPI environment */
      static void init (int argc, char **argv);

      /** \brief Terminate the MPI environment and the program */
      static void finalize ();


      // --- wrappers for standard MPI operations ---

      /** \brief Return the MPI rank of the current process */
      static int me ();

      /** \brief Return the number of MPI ranks */
      static int nr ();

      /** \brief Determine if the current process is responsible
       *         for a piece of work specified by an index i */
      static bool myWork (int i);

      /** Never ever use SxAutoLoop on MPI-loops!
          The autoloop limits cannot be detected if they are in
          a no-work MPI task.
        */
      static bool myWork (const SxAutoLoop &) { SX_EXIT; }

      /** \brief Determine which rank has to do a piece of work
       *         specified by an index i */
      static int whoseWork (int i);

      /** \brief Return the MPI rank of the current process */
      static int invMap ();

      /** \brief Logical OR */
      static bool lor (bool in);

      /** \brief Calculate the sum over in from all ranks and return the
       *         result (on all ranks) */
      static int sum (int in);

      /** \brief Calculate the sum over in from all ranks and return the
       *         result (on all ranks) */
      static double sum (double in);

      /** \brief This is the placeholder for arbitrary summable types
        If some type should support MPI summation, please add a
        customized wrapper for the C-pointer implementation sum
        as a template specialization of SxLoopMPI::sum.
        \example
        An implementation for SxArray<SxComplex16>
\code
template<>
inline void SxLoopMPI::sum (SxArray<SxComplex16> &inout) 
{
#ifdef USE_MPI
   // note: SxComplex16 is memory-compatible to double[2]
   sum ((double*)inout.elements, (double*)inout.elements, inout.getSize () * 2);
#endif
}
\endcode
        */
      template<class T>
      static void sum (T &inout);

      /** \brief Distribute the double "in" from rank "source" to all ranks */
      static double bcast (double in, int source);

      /** \brief Distribute the integer "in" from rank "source" to all ranks */
      static int bcast (int in, int source);

      /** \brief This is the placeholder for arbitrary type broadcasting
        If some type should support MPI broadcasting, please add a
        customized wrapper for the C-pointer implementation bcast
        as a template specialization of SxLoopMPI::bcast.
        \example
        An implementation for SxArray<SxComplex16>
\code
template<>
inline void SxLoopMPI::bcast (SxArray<SxComplex16> &inout, int source) 
{
#ifdef USE_MPI
   // note: SxComplex16 is memory-compatible to double[2]
   bcast ((double*)inout.elements, inout.getSize () * 2, source);
#endif
}
\endcode
        */
      template<class T>
      static void bcast (T &inout, int source);

      /** \brief Wait until all MPI ranks have arrived, then continue */
      static void barrier ();


      /** \brief Method to switch parallel FFTs ON  */
      static void setParallelFft (bool mode);

      /** \brief Method to get the current status for parallel FFTs */
      static bool parallelFft ();

      /** \brief Hierarchical parallelization / put level onto stack */
      static void pushMpiLevel(SxTaskGroup *);

      /** \brief Hierarchical parallelization / take level from stack */
      static SxTaskGroup* popMpiLevel();


#ifdef USE_MPI
      /* Native MPI datatypes are required by NetCDF4 parallel IO in <SxBinIo.cpp> */
      /** \brief Returns the global MPI communicator */
      static MPI_Comm MpiCommWorld();
      /** \brief Returns MPI_INFO_NULL */
      static MPI_Info MpiInfoNull();
#endif

   // ----- MPI data source and target selector classes
   private:
      /// Source and target levels for sum and broadcast
      class TaskSelection {
         public:
            /// Name of the task group level (NULL if not specified by name)
            const char* level;
            /// Whether only master or all ranks participate
            SxMPITaskSelection   who;
            /// Whether group is specified by name or by something else
            SxMPILevelSelection  levelType;
            /// Constructor
            TaskSelection ()
               : level(NULL), who(TaskGroupAll), levelType(TopLevel)
            { /* empty */ }
      };

      /// The MPI data source
      static TaskSelection MPIsource;
      /// The MPI data target
      static TaskSelection MPItarget;

   public:
      /// Overload temporarily the MPI data source
      class OverloadSource {
         private:
            TaskSelection oldValue;
         public:
            /// Constructor: specify name of task group level
            OverloadSource (const char *levelName, SxMPITaskSelection who)
               : oldValue (MPIsource)
            {
               MPIsource.level     = levelName;
               MPIsource.levelType = SpecifiedName;
               MPIsource.who       = who;
            }

            /// Constructor: specify type of task group level
            OverloadSource (SxMPILevelSelection levelIn, 
                            SxMPITaskSelection who = TaskGroupAll)
               : oldValue (MPIsource)
            {
               SX_CHECK (levelIn != SpecifiedName);
               MPIsource.levelType = levelIn;
               MPIsource.who       = who;
            }

            /// Restore previous MPI data source
            ~OverloadSource ()
            {
               MPIsource = oldValue;
            }
      };

      /// Overload temporarily the MPI data target
      class OverloadTarget {
         private:
            TaskSelection oldValue;
         public:
            /// Constructor: specify name of task group level
            OverloadTarget (const char *levelName, SxMPITaskSelection who)
               : oldValue (MPItarget)
            {
               MPItarget.level     = levelName;
               MPItarget.levelType = SpecifiedName;
               MPItarget.who       = who;
            }

            /// Constructor: specify type of task group level
            OverloadTarget (SxMPILevelSelection levelIn, 
                            SxMPITaskSelection who)
               : oldValue (MPItarget)
            {
               SX_CHECK (levelIn != SpecifiedName);
               MPItarget.levelType = levelIn;
               MPItarget.who       = who;
            }

            /// Restore previous MPI data target
            ~OverloadTarget ()
            {
               MPItarget = oldValue;
            }
      };

   protected:
      /// Try to find the task group specified in lvl
      static SxMpiTaskGroup* getLevel(const TaskSelection &lvl);

      /// Find the target once the level for MPIsource is known
      static SxMpiTaskGroup* getTarget(SxTaskGroup *source);

      /** \brief The  main implementation of the sum routine
          This routine takes care of MPIsource and MPItarget for
          multilevel summations.
        */
      template<class T>
      static void sum(T* in, T* out, ssize_t nElem);
   public:
      /** \brief The  main implementation of the bcast routine
          This routine takes care of MPIsource and MPItarget for
          broadcasts.

          @Example:
          Within a certain level "k-points", broadcast value X to all
          MPI tasks belonging to the same (logical) sibling within that group
          (i.e., all MPI tasks working on the same k-point). The source
          of the broadcast is the task group's master (i.e. task 0 within
          the sibling group).
          
          \code
   SX_MPI_LEVEL("k-points");
   SX_MPI_SOURCE (CurrentLevel, TaskGroupMaster);
   SX_MPI_TARGET (CurrentLevel, TaskGroupAll);
   bcast (X, 0);
          \endcode
            
          @Example:
          Within a certain level "k-points", broadcast value X to all
          MPI tasks belonging to all siblings within that group
          (i.e., all MPI tasks working on k-points). The source
          of the broadcast is one siblings task group's master (i.e.
          task 0 within the sibling group).
          \code
   SX_MPI_LEVEL("k-points");
   SX_MPI_SOURCE (CurrentLevel, TaskGroupMaster);
   SX_MPI_TARGET (LevelSiblings, TaskGroupAll);
   // X has been computed for index ik by one sibling group.
   bcast (X, SxLoopMPI::whoseWork(ik));
   res(ik) = X;
          \endcode
          
          @Example:
          Assume a hierarchy with a high-level parallelization over
          k-points and a lower-level parallelization over states.
          Within the level "states", broadcast value X to all
          MPI tasks working on states for the same k-point. The source
          of the broadcast is one siblings task group's master (i.e.
          task 0 within the sibling group).
          \code
   SX_MPI_LEVEL("states");
   SX_MPI_SOURCE (CurrentLevel, TaskGroupMaster);
   SX_MPI_TARGET (LevelSiblings, TaskGroupAll);
   // X has been computed for index iState by one sibling group.
   bcast (X, SxLoopMPI::whoseWork(iState));
   res(iState) = X;
          \endcode
          
          @Example:
          Assume a hierarchy with a high-level parallelization over
          k-points and a lower-level parallelization over states.
          Within the level "states", broadcast value X to all
          MPI task groups working on states for the same k-point. The source
          of the broadcast is one siblings task group's master (i.e.
          task 0 within the sibling group). The target is each siblings
          master (who may further manipulate the data on its own
          before any further synchronization occurs).
          \code
   SX_MPI_LEVEL("states");
   SX_MPI_SOURCE (CurrentLevel, TaskGroupMaster);
   SX_MPI_TARGET (LevelSiblings, TaskGroupMaster);
   // X has been computed for index iState by one sibling group.
   bcast (X, SxLoopMPI::whoseWork(iState));
   res(iState) = X;
          \endcode
          
          @Example:
          Broadcast a value X from one MPI task to all others (global)
          \code
   SX_MPI_SOURCE (TopLevel);
   SX_MPI_TARGET (TopLevel);
   bcast (X, sourceId);
          \endcode
          
        */
      template<class T>
      static void bcast(T* inout, ssize_t nElem, int sourceId);
};

/** \brief This class introduces additional MPI levels
           according to the specifications given to it in the
           append routine.
    \example
    \code
SxAutoLevel ()
   .append ("A", nA) // level A loops over nA items
   .append ("B", nB, 0.2) // level B loops over nB items, but 20% of the time
                          // cannot be parallelized (e.g. pre-loop setup)
   .append ("C", nC);
  */
class SX_EXPORT_MPI SxAutoLevel {
   private:
      /// \brief Auxiliary class to store the level description from append
      class Workload {
         public:
            SxString name;
            int work;
            double serial;
            Workload (int workIn, double serialIn)
               : work(workIn), serial(serialIn) {}
            Workload () {}
      };

      /// List of desired levels
      SxList<Workload> levels;

      /** \brief Estimate optimal efficiency of a given level (and below) for
                 a given number of MPI tasks. This routine calls itself
                 recursively to optimize the efficiency of lower levels.
          @param nMember number of MPI tasks available
          @param iLevel  index n the list of levels
          @param optSiblings work space: optimal number of MPI tasks
                 (iLevel:nMembers-1)
          @param optEff work space: best efficiencies (iLevel:nMembers-1)
        */
      double getEfficiency (int nMember, int iLevel,
                            SxArray2<int> &optSiblings,
                            SxArray2<double> &optEff);
      bool finished;
      SxString context;
   public:
      SxAutoLevel () : finished(false), context(SxString("auto")) { }
      SxAutoLevel (const SxString &context_) : finished(false), context(context_) { }

      /** \brief Register an MPI level for automatic addition
          @param name           name of the level
          @param work           loop size over which this level parallelizes
          @param serialFraction relative amount of work related to this
                                loop that can be parallelized by outer 
                                loops only (if known).

          The current strategy attempts to reduce idle time for the
          last round of indices, when the later MPI tasks have one less
          index to do than the first ones.
          Large values of serialFraction will move parallelization to
          outer loops even if the idle time is larger there.
        */
      inline SxAutoLevel& append (const char *name,
                                  int work,
                                  double serialFraction = 0.);

      /** Do the actual optimization. Done automatically in destructor. */
      void finalize ();
      ~SxAutoLevel () 
      {
#ifdef USE_LOOPMPI
         if (!finished) finalize ();
#endif
      }
};

#ifdef USE_LOOPMPI
SxAutoLevel& 
SxAutoLevel::append (const char *name,
                         int work,
                         double serialFraction)
{
   SX_CHECK (!finished);
   levels << Workload(work, serialFraction);
   levels.last ().name = name;
   return *this;
}
#else
SxAutoLevel& SxAutoLevel::append (const char *, int , double)
{
   return *this;
}
#endif

#ifndef USE_MPI
// --- empty inline implementations for MPI-specific operations when
//     MPI is not used
inline int    SxLoopMPI::sum (int in) { return in; }
inline double SxLoopMPI::sum (double in) { return in; }
inline double SxLoopMPI::bcast (double in, int) { return in; }
inline int    SxLoopMPI::bcast (int in, int) { return in; }
inline void SxLoopMPI::barrier () { }
template<class T> inline void SxLoopMPI::bcast(T*, ssize_t, int)
{
   // empty
}
/*
// --- the following template functions should never be instantiated
//     without MPI
template<class T>
inline void SxLoopMPI::sum (T* in, T* out, ssize_t nElem)
{
   if (in == out) return;
   // copy data
   memcpy (out, in, nElem * sizeof(T));
}
*/
#endif

/// broadcast contents of a bool array (must be appropriately sized)
template<>
inline void SxLoopMPI::bcast<SxArray<bool> > (SxArray<bool> &array, int source)
{
   SX_CHECK (array.getSize () > 0); // assuming that the array is properly sized!
   bcast (reinterpret_cast<char*>(array.elements),
          array.getSize () * sizeof(bool), source);
}

#endif // _SX_LOOPMPI_
