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

#ifndef _SX_THREAD_H_
#define _SX_THREAD_H_

#include <SxString.h>
#include <SxConfig.h>
#include <SxTimer.h>
#include <SxPtr.h>
#include <SxIPC.h>
#include <SxSigSlots.h>

class SxThreadPool;

/** \brief Interface for pre-spawned threads in a pool.

    Each SxThread is part of a SxThreadPool which is managing the actual
    operating system threads (SxSystemThread).

    The pool will pre-spawn such threads so that SxThread does not have the
    cost of starting and stopping a operating system thread. Creating and
    disposing SxThread is light-weight. Execution speed depends on size of the
    pool and its load.

    A new thread object is created by deriving from SxThread and overloading
    of the SxThread::main function. The thread is started using SxThread::start
    and can be joined with other threads with SxThread::wait.

    \par Example:
\code
   class MyThread : public SxThread
   {
      public:
         int id;

         MyThread (int id_) : SxThread (), id(id_) { }

         virtual void main ()  {
            for (int i=0; i < 10; ++i) {
               printf ("This is thread id %d, i=%d\n", id, i); fflush (stdout);
               SxThread::sleep(1);
            }
         }
   };

   int main ()
   {
      SxThreadPool pool;
      MyThread t1(1), t2(2);
      t1.setThreadPool (pool);
      t2.setThreadPool (pool);
      t1.start ();
      t2.start ();
      t1.wait ();
      t2.wait ();
      printf ("Back in main thread\n");
      return 0;
   }
\endcode

    \ingroup group_os
    \author  Sixten Boeck, boeck@sfhingx.de */
class SX_EXPORT_IPC SxThread : public SxThis<SxThread>
{
   public:

      enum Priority { MaxPriority = 100 }; //  0 - lowest priority
                                           // 99 - highest priority

      /** \brief Create a new thread

          Specify a tag to identify this particular task. Its useful to
          relate exceptions. It may also be used to pass flags.

          In order to start the thread the SxThread::start function has to
          be called explicitly. */
      SxThread (const SxString &tag=0,
                int priority_=0,
                double maxRunTime_=-1.,
                double maxSpoolTime_=-1.);

      SxThread (const SxString &tag, const SxPtr<SxThreadPool> &pool_);

      /** \brief Destructor.

          The destrcutor does not terminate a thread! Waiting for the
          end of a thread is done using the barrier function SxThread::wait. */
      virtual ~SxThread ();

      void setThreadPool (const SxPtr<SxThreadPool> &pool_);

      /** \brief Start a new thread

          This function sets up the thread environment and execute the
          SxThread::main function of the derived class.
          \param priority  The priority used for scheduling the new thread */
      void start ();

      /** \brief Block a thread and wait its termination.

           This function creates a barrier and waits until the thread has
           (been) terminated.

           \par Example:
\code
   MyThread t1, t2;

   // --- multi-threadding part starts here
   t1.start (); t2.start ();
   t1.end();    t2.end();
   // --- end of multi-threadding part
\endcode
       */
      void wait ();


      /** \brief Terminates a thread

          This function terminates the thread by calling the corresponding
          kill function of the underlying thread library. */
      void terminate ();

      /** \brief Slot for computation routine of a thread

          When the thread is executed from SxThread::start this virtual
          function is called. When defining a new thread type by deriving
          from SxThread, override this function. It should contain the
          computational part which runs threadded.

          The tag member is available an can be used for identification and/or
          inspected for flags. */
      virtual void main ();

      /** \brief Called when a exception occurs in main()

           Override this funtion and implement exception handling specific to
           the main routine. Note that tag is available to identify the
           causing thread. */
      virtual void handleException (const SxException &e);

      /** send signals before and after calling main()
        */
      virtual void mainWrap ();

      /** \brief Sleep for some seconds

        Calling this function makes the current process/thread sleep until
        a certain time has elapsed. This function can be applied to threads
        or to entire processes.

        \par Example 1
        Making a process sleep for 5 seconds
\code
   int main ()
   {
      for (int i=0; i < 100; ++i)  {
         printf ("i=%d\n", i);
         SxThread::sleep(5);
      }
   }
\endcode

       \par Example 2
       Making a single thread sleep for 1 second for observing changes of
       a variable (which might be controlled from the slave thread.
\code
   class MyMasterThread : public SxThread
   {
      public:
         MyMasterThread () : SxThread (), f(false)  {  }

         void setFlag (bool f)  {  flag = f; }

         virtual void main ()
         {
             for (;;)  {
                if (flag)  doSomeThing ();
                SxThread::sleep(1);
             }
         }
      protected:
         bool flag;
   };

   class MySlaveThread : public SxThread
   {
      public:
         MySlaveThread (MyMasterThread *m) : SxThread (), master(m);

         virtual void main ()
         {
            // do some heavy calculation and communicate with the master
            for (i=0; i<10000; ++i)  {
               // eat CPU time
               if (...)  master->setFlag (true);  // call master process
            }
         }

      protected:
         MyMasterThread *master;
   };

   int main ()
   {
      SxThreadPool pool;
      MyMasterThread tMaster;
      MySlaveThread  tSlave(&tMaster);
      tMaster.setThreadPool (pool.getThis());
      tSlave.setThreadPool (pool.getThis());
      tMaster.start ();
      tSlave.start ();
      tMaster.wait ();
      tSlave.wait ();
      return 0;
   }
\endcode */
      static void sleep (unsigned int seconds);
      static void millisleep (unsigned int mseconds);

      bool isRunning () const;

      int getPriority () const { return priority; }
      bool spoolingExpired () const;
      bool runningExpired () const;

      int updatePriority ();

   signals:

      SxSignal<SxThread *, const char *> SX_SIGNAL (sigStarted);
      SxSignal<SxThread *, const char *> SX_SIGNAL (sigFinished);

   protected:
      const SxString tag;  // to identify and/or pass flags
      SxPtr<SxThreadPool> pool; // pool container

      // Task variables
      enum Timer { Spooling=0, Running=1 };

      int priority;
      double maxRunTime;
      double maxSpoolTime;

      SxTimer timer;

      int prioPoints;
};

#endif /* _SX_THREAD_H_ */
